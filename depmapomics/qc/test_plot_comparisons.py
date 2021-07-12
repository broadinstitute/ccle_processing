import matplotlib.pyplot as plt
import pandas as pd
import pytest
import seaborn as sns
from depmapomics.qc.config import PLOTS_OUTPUT_FILENAME_PREFIX
from depmapomics.qc.test_compare_to_ref_release import (
    FILE_ATTRIBUTES_PAIRED, REFERENCE_RELEASE, VIRTUAL_RELEASE, data,
    get_both_releases_from_taiga)

NEW_TO_OLD_CORRELATION_THRESHOLD = 0.95
SHARED_DATA_CORRELATION_THRESHOLD = 0.95
MIN_SAMPLESIZE_FOR_CORR = 10

def get_data_stack(file, number_of_points=1000000, random_state=0):
    data1, data2 = get_both_releases_from_taiga(file)

    row = set(data1.index) & set(data2.index)
    col = set(data1.columns) & set(data2.columns)

    data1_stack = data1.loc[row, col].stack()
    data2_stack = data2.loc[row, col].stack()

    if number_of_points > 0:
        data1_stack = data1_stack.sample(number_of_points, random_state=random_state)
    data2_stack = data2_stack.loc[data1_stack.index]

    data_stack = pd.concat([data1_stack, data2_stack], axis=1)

    cols = ['{:s}.{:d}'.format(REFERENCE_RELEASE['name'], REFERENCE_RELEASE['version']),
            '{:s}.{:d}'.format(VIRTUAL_RELEASE['name'], VIRTUAL_RELEASE['version'])]
    data_stack.columns = cols
    data_stack.reset_index(inplace=True)
    data_stack.rename(columns={'level_0': 'DepMap_ID', 'level_1': 'gene'}, inplace=True)
    return data_stack, cols


@pytest.fixture(scope='module')
def data_stack(request):
    data_stack, cols = get_data_stack(request.param)
    return data_stack, cols


@pytest.fixture(scope='function')
def CCLE_gene_cn_with_source_change():
    CCLE_gene_cn_12, cols = get_data_stack('CCLE_gene_cn')

    names = [REFERENCE_RELEASE['name'], VIRTUAL_RELEASE['name']]

    CCLE_segment_cn_1, CCLE_segment_cn_2 = get_both_releases_from_taiga('CCLE_segment_cn')

    sources = pd.merge(CCLE_segment_cn_1[['DepMap_ID', 'Source']].drop_duplicates(),
         CCLE_segment_cn_2[['DepMap_ID', 'Source']].drop_duplicates(),
         on='DepMap_ID', suffixes=['_'+names[0], '_'+names[1]])
    del CCLE_segment_cn_1, CCLE_segment_cn_2
    sources['source_change'] = sources.apply(lambda x: '{:s} -> {:s}'.format(x['Source_'+names[0]], x['Source_'+names[1]]), axis=1)
    sources['source_has_changed'] = (sources['Source_'+names[0]] != sources['Source_'+names[1]])

    CCLE_gene_cn_12 = pd.merge(CCLE_gene_cn_12, sources, on='DepMap_ID')
    return CCLE_gene_cn_12, cols


@pytest.mark.skipif([1 for x in FILE_ATTRIBUTES_PAIRED if x['file']=='CCLE_gene_cn'] == [], reason='skipped by user')
@pytest.mark.plot
def test_plot_gene_cn_comparison(CCLE_gene_cn_with_source_change):

    CCLE_gene_cn_12, cols = CCLE_gene_cn_with_source_change

    plt.figure(figsize=(20,10))
    sns.scatterplot(data=CCLE_gene_cn_12, x=cols[0], y=cols[1],
                    hue='source_change', style='source_has_changed', alpha=0.5, cmap='Tab20')

    output_img_file = PLOTS_OUTPUT_FILENAME_PREFIX + '{}-vs-{}.png'.format(cols[1], cols[0])
    print('saved to {}'.format(output_img_file))
    plt.savefig(output_img_file, bbox_inches='tight')
    plt.close()
    corr = CCLE_gene_cn_12[cols].corr().iloc[0, 1]
    assert corr > SHARED_DATA_CORRELATION_THRESHOLD



PARAMS_plot_per_gene_means = [(x['file'], x) for x in FILE_ATTRIBUTES_PAIRED if x['ismatrix']]
@pytest.mark.parametrize('data, file_attr', PARAMS_plot_per_gene_means, indirect=['data'])
@pytest.mark.plot
def test_plot_per_gene_means(data, file_attr):
    data1, data2 = data

    new_lines = set(data2.index) - set(data1.index)
    old_lines = set(data2.index) & set(data1.index)

    stats_old = data2.loc[old_lines].mean()
    stats_new =data2.loc[new_lines].mean()

    data_compare_stats = pd.concat([stats_old, stats_new], keys=['old', 'new'], axis=1)
    corr = data_compare_stats.corr().iloc[0, 1]

    if file_attr['omicssource'] == 'RNA':
        data_compare_stats['ERCC'] = data_compare_stats.index.map(lambda x: x.startswith('ERCC-'))
        hue = 'ERCC'
    else:
        hue = None
    print(data_compare_stats.head())
    sns.scatterplot(data=data_compare_stats, x='old', y='new', hue=hue)

    minmax = (data_compare_stats.min().min(), data_compare_stats.max().max())
    plt.xlabel('per gene values in the new release\n averaged across cell lines shared between the new and old release\n(n={:d})'.format(len(old_lines)))
    plt.ylabel('per gene values in the new release\n averaged across cell lines only in the new release\n(n={:d})'.format(len(new_lines)))
    plt.plot(minmax, minmax, 'r--')
    plt.title('{}: corr = {:.2f}'.format(file_attr['file'], corr), fontsize=12)
    output_img_file = PLOTS_OUTPUT_FILENAME_PREFIX + 'newlines_per_gene_means_{}.png'.format(file_attr['file'])
    print('saved to {}'.format(output_img_file))
    plt.savefig(output_img_file, bbox_inches='tight')
    plt.close()
    assert (corr > NEW_TO_OLD_CORRELATION_THRESHOLD) | (len(stats_new) >= MIN_SAMPLESIZE_FOR_CORR)



PARAMS_plot_matrix_comparison = [(x['file'], x['file']) for x in FILE_ATTRIBUTES_PAIRED if x['ismatrix']]
@pytest.mark.parametrize('data_stack, file', PARAMS_plot_matrix_comparison, indirect=['data_stack'])
@pytest.mark.plot
def test_plot_matrix_comparison(data_stack, file):
    data_stack_df, cols = data_stack
    corr = data_stack_df.corr().iloc[0, 1]
    sns.scatterplot(data=data_stack_df, x=cols[0], y=cols[1])
    minmax = (data_stack_df[cols].min().min(), data_stack_df[cols].max().max())
    plt.plot(minmax, minmax, 'r--')
    plt.xlabel(cols[0])
    plt.ylabel(cols[1])
    plt.title('{}: corr = {:.2f}'.format(file, corr), fontsize=12)
    output_img_file = PLOTS_OUTPUT_FILENAME_PREFIX + 'sharedlines_per_gene_values_{}.png'.format(file)
    print('saved to {}'.format(output_img_file))
    plt.savefig(output_img_file, bbox_inches='tight')
    plt.close()
    assert corr > SHARED_DATA_CORRELATION_THRESHOLD
