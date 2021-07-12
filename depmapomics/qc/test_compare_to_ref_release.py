import numpy as np
import pandas as pd
import pytest
from depmapomics.qc.config import (CORRELATION_THRESHOLDS, LEGACY_PATCH_FLAGS,
                                      FILE_ATTRIBUTES, FILES_RELEASED_BEFORE,
                                      REFERENCE_RELEASE,
                                      SKIP_ARXSPAN_COMPARISON, VIRTUAL_RELEASE)
from taigapy import TaigaClient

tc = TaigaClient()

FILE_ATTRIBUTES_PAIRED = [x for x in FILE_ATTRIBUTES if x['file'] in FILES_RELEASED_BEFORE]


def tsv2csv(df):
    df.to_csv('/tmp/data.tsv', index=False)
    df = pd.read_csv('/tmp/data.tsv', sep='\t')
    return df

####### FIXTURES ####
def get_both_releases_from_taiga(file):
    data1 = tc.get(name=REFERENCE_RELEASE['name'], file=file, version=REFERENCE_RELEASE['version'])
    data2 = tc.get(name=VIRTUAL_RELEASE['name'], file=file, version=VIRTUAL_RELEASE['version'])
    if LEGACY_PATCH_FLAGS['tsv2csv']:
        # some older taiga data formats (probably 20q1 and older) are tsv and deprecated
        data1 = tsv2csv(data1)

    if LEGACY_PATCH_FLAGS['rename_column']:
        # in 21q1 CCLE_mutations file this column was renamed
        data1.rename(columns={'Tumor_Allele': 'Tumor_Seq_Allele1'}, inplace=True)
    return data1, data2

def merge_dataframes(file, merge_columns):
    # TODO: figure out how to call the data fixture instead
    data1, data2 = get_both_releases_from_taiga(file)
    data_merged = pd.merge(data1, data2, on=merge_columns, indicator=True, how='outer')
    return data_merged

@pytest.fixture(scope='module')
def data(request):
    return get_both_releases_from_taiga(request.param)

@pytest.fixture(scope='module')
def dataframes_merged(request):
    merged_df =  merge_dataframes(request.param[0], request.param[1])
    return merged_df


##### TESTS ############
# PARAMS_compare_column_names = [x['file'] for x in FILE_ATTRIBUTES if not x['ismatrix']]
PARAMS_compare_column_names = [x['file'] for x in FILE_ATTRIBUTES_PAIRED]
@pytest.mark.parametrize('data', PARAMS_compare_column_names, indirect=['data'])
@pytest.mark.compare
def test_compare_column_names(data):
    data1, data2 = data
    assert set(data1.columns) == set(data2.columns)


PARAMS_matrix_correlations = [(x['file'],
    CORRELATION_THRESHOLDS['CCLE_gene_cn'] if (x['file'] == 'CCLE_gene_cn') else CORRELATION_THRESHOLDS['all_expressions'])
    for x in FILE_ATTRIBUTES_PAIRED if x['ismatrix'] & (x['omicssource']=='RNA')]
@pytest.mark.parametrize('method', ['spearman',
                         pytest.param('pearson',
                         marks=pytest.mark.skip(reason='Pearson can be sensitive to outliers'))])
@pytest.mark.parametrize('axisname', ['persample',
                        pytest.param('pergene',
                        marks=pytest.mark.skip(reason='If persample fails usually this fails too'))])
@pytest.mark.parametrize('data, threshold', PARAMS_matrix_correlations, indirect=['data'])
@pytest.mark.compare
def test_matrix_correlations(data, threshold, axisname, method):
    axis = 0 if axisname == 'pergene' else 1 if axisname == 'persample' else None
    data1, data2 = data
    corrs = data1.corrwith(data2, axis=axis, drop=True, method=method)
    # TODO: tests for NAs instead of dropping them. It seems like NA can happen if there are all zeros in one of the vectors, so let's dropna for now
    corrs.dropna(inplace=True)
    assert (corrs >= threshold).all(), 'the cases which did not pass the tests are:\n{}'.format(corrs[corrs<threshold].sort_values())


PARAMS_fraction_of_unequl_columns_from_merged_file = [((x['file'], x['merge_cols']), x['expected_changed_cols'])
                                                      for x in FILE_ATTRIBUTES_PAIRED if 'merge_cols' in x]
@pytest.mark.parametrize('dataframes_merged, expected_changed_cols', PARAMS_fraction_of_unequl_columns_from_merged_file,
                         indirect=['dataframes_merged'], ids=[x[0][0] for x in PARAMS_fraction_of_unequl_columns_from_merged_file])
@pytest.mark.compare
def test_fraction_of_unequal_columns_from_merged_file(dataframes_merged, expected_changed_cols):
    dataframes_merged.drop([x+'_x' for x in expected_changed_cols], inplace=True, axis=1)
    dataframes_merged.drop([x+'_y' for x in expected_changed_cols], inplace=True, axis=1)
    cols = list(set([x[:-2] for x in dataframes_merged.columns if x.endswith('_x') | x.endswith('_y')]))
    dataframe_merge_both = dataframes_merged[dataframes_merged['_merge'] == 'both']
    dataframe_merge_both.set_index('DepMap_ID', inplace=True)
    unequal_values = pd.DataFrame(index=dataframe_merge_both.index, columns=cols)
    cols_dtype = dataframe_merge_both[[col+'_x' for col in cols]].dtypes
    equal_nonNA = lambda a, b: (a == b) | ((a != a) & (b != b)) # this is a regular equality tests with the exception that NA==NA
    almost_equal_nonNA = lambda a, b: np.isclose(a, b) | ((a != a) & (b != b)) # this is a regular equality tests with the exception that NA==NA
    for col in cols:
        if cols_dtype[col+'_x'] == np.dtype('float64'): # otherwise very close values will look different
            unequal_values[col] = ~almost_equal_nonNA(dataframe_merge_both[col+'_x'], dataframe_merge_both[col+'_y'])
        else:
            unequal_values[col] = ~equal_nonNA(dataframe_merge_both[col+'_x'], dataframe_merge_both[col+'_y'])

    unequal_columns = unequal_values.agg(['mean', 'sum']).T
    unequal_columns.sort_values('mean', ascending=False, inplace=True)
    unequal_columns['sum'] = unequal_columns['sum'].astype(int)
    unequal_columns = unequal_columns[unequal_columns['sum'] > 0]
    unequal_columns.rename(columns={'mean': 'freq', 'sum': 'count'}, inplace=True)

    unequal_values_sum = unequal_values.groupby('DepMap_ID').sum()
    unequal_values_sum = unequal_values_sum[(unequal_values_sum > 0).any(axis=1)]
    unequal_values_sum = unequal_values_sum.loc[:, (unequal_values_sum > 0).all()]

    assert unequal_columns.empty, 'fraction of unequal values in each column that are expected to be equal:\n{}\
        \n\ncell lines affected by these changes:\n{}'.format(unequal_columns, unequal_values_sum)


PARAMS_compare_nan_fractions = [x['file'] for x in FILE_ATTRIBUTES_PAIRED if not x['ismatrix']]
@pytest.mark.parametrize('arxspans', ['allcells', 'sharedcells'])
@pytest.mark.parametrize('data', PARAMS_compare_nan_fractions, indirect=True)
@pytest.mark.compare
def test_compare_nan_fractions(data, arxspans, atol=1e-2):
    data1, data2 = data
    if arxspans == 'sharedcells':
        # subset data1 and data2 by shared arxspan IDs
        arxspans_ids = set(data1['DepMap_ID']) & set(data2['DepMap_ID'])
        data1 = data1[data1['DepMap_ID'].isin(arxspans_ids)]
        data2 = data2[data2['DepMap_ID'].isin(arxspans_ids)]
    elif arxspans != 'allcells':
        raise Exception('unknown value for arxspans')

    nan_fractions = pd.concat([data1.isnull().mean(), data2.isnull().mean()], axis=1, join='inner', keys=['old', 'new'])
    nan_fractions.sort_values('new', ascending=False, inplace=True)
    nan_fractions = nan_fractions[~np.isclose(nan_fractions['old'], nan_fractions['new'], atol=atol)]
    assert nan_fractions.empty, 'the following NA fractions are different according to tolerance {:.1e}:\n{}'.format(atol, nan_fractions)


PARAMS_compare_column_dtypes = [x['file'] for x in FILE_ATTRIBUTES_PAIRED]
@pytest.mark.parametrize('method', ['pd_dtypes', 'map_type'])
@pytest.mark.parametrize('data', PARAMS_compare_column_dtypes, indirect=['data'])
@pytest.mark.compare
def test_compare_column_dtypes(data, method):
    data1, data2 = data
    if method == 'pd_dtypes': # per column dtype
        get_dtypes = lambda df: df.dtypes
    elif method == 'map_type': # per element type
        get_dtypes = lambda df: df.apply(lambda x: x.dropna().map(type).unique()).T[0]
    else:
        raise Exception('bad values for dtype method')

    dtypes_compare = pd.concat([get_dtypes(data1), get_dtypes(data2)], axis=1, keys=['reference', 'virtual'])
    dtypes_compare.dropna(inplace=True)
    dtypes_compare_nonmatching = dtypes_compare.apply(lambda x: x[0]!=x[1], axis=1)
    dtypes_compare = dtypes_compare[dtypes_compare_nonmatching]
    assert dtypes_compare.empty, 'the following columns have changed types between the releases:\n{}'.format(dtypes_compare)


PARAMS_compare_cell_lines = [(x['file'], 'index' if x['ismatrix'] else 'DepMap_ID') for x in FILE_ATTRIBUTES_PAIRED]
@pytest.mark.skipif(SKIP_ARXSPAN_COMPARISON, reason="skipping, since it is normal to have arxspan differences between the releases")
@pytest.mark.parametrize('data, arxspan_col', PARAMS_compare_cell_lines, indirect=['data'])
@pytest.mark.compare
def test_compare_cell_lines_released(data, arxspan_col):
    data1, data2 = data
    if arxspan_col == 'index':
        arxspans1 = set(data1.index)
        arxspans2 = set(data2.index)
    elif arxspan_col == 'DepMap_ID':
        arxspans1 = set(data1['DepMap_ID'])
        arxspans2 = set(data2['DepMap_ID'])
    else:
        raise Exception('unknown value provided for arxspan_col')
    assert arxspans1 == arxspans2#, 'lines added:\n{}\nlines removed:\n {}'.format(', '.join(arxspans2-arxspans1), ', '.join(arxspans1-arxspans2))


@pytest.mark.skipif([1 for x in FILE_ATTRIBUTES_PAIRED if x['file']=='CCLE_segment_cn'] == [], reason='skipped by user')
@pytest.mark.parametrize('data', ['CCLE_segment_cn'], indirect=['data'])
@pytest.mark.compare
def test_source_changes(data):
    data1, data2 = data
    source1 = data1.groupby('DepMap_ID')['Source'].apply(lambda x: x.iloc[0])
    source2 = data2.groupby('DepMap_ID')['Source'].apply(lambda x: x.iloc[0])
    source_changes = pd.concat([source1, source2], axis=1, keys=['old', 'new'])
    source_changes_matrix = source_changes.groupby(['old', 'new']).size().unstack(fill_value=0)
    source_changes = source_changes[source_changes['old'] != source_changes['new']]

    assert source_changes.empty, 'the following cell lines have had a source change in CCLE_segment_cn:\n{}\n\nSource change matrix (counting cell lines):\n {}'.format(source_changes, source_changes_matrix)



# TODO: implement the assert almost equal version of the tests
# from pandas.testing import assert_frame_equal
# PARAMS_matrix_correlations = [('CCLE_expression_full', 0.01), ('CCLE_gene_cn', 0.01)]
# # @pytest.mark.parametrize('axisname', ['pergene', 'persample'])
# @pytest.mark.parametrize('data, rtol', PARAMS_matrix_correlations, indirect=['data'])
# def test_matrix_equality(data, rtol):
#     # axis = 0 if axisname == 'pergene' else 1 if axisname == 'persample' else None
#     data1, data2 = data
#     row = set(data1.index) & set(data2.index)
#     col = set(data1.columns) & set(data2.columns)
#     data1_ = data1.loc[row, col].T
#     data2_ = data2.loc[row, col].T
#     assert_frame_equal(data1_, data2_, rtol=rtol)
