import re

import numpy as np
import pytest
from depmapomics.qc.config import (FILE_ATTRIBUTES, VIRTUAL_RELEASE)
from taigapy import TaigaClient

tc = TaigaClient()

@pytest.fixture(scope='module')
def data(request):
    return tc.get(name=VIRTUAL_RELEASE['name'], file=request.param, version=VIRTUAL_RELEASE['version'])


PARAMS_wrong_columns = [x['file'] for x in FILE_ATTRIBUTES]
@pytest.mark.parametrize('data', PARAMS_wrong_columns, indirect=True)
@pytest.mark.format
def test_wrong_columns(data):
    wrong_columns = {'Unnamed: 0'}
    assert wrong_columns & set(data.columns) == set()


PARAMS_test_symbol_and_entrezid_in_column = [x['file'] for x in FILE_ATTRIBUTES if x['ismatrix'] and x['gene_id']=='entrez']
@pytest.mark.parametrize('data', PARAMS_test_symbol_and_entrezid_in_column, indirect=True)
@pytest.mark.format
def test_symbol_and_entrezid_in_column(data):
    matches = data.columns.map(lambda x: re.match(r'([a-zA-Z0-9-_/.]+)\s\((\d+)\)$', x))
    assert  matches.notnull().all(), \
        'some columns do not follow the symbol (entrez id) format. The first few are: \n{}'.format(data.columns[matches.isnull()][:20])


PARAMS_test_symbol_and_enstid_in_column = [x['file'] for x in FILE_ATTRIBUTES if x['ismatrix'] and x['gene_id']=='enst']
@pytest.mark.parametrize('data', PARAMS_test_symbol_and_enstid_in_column, indirect=True)
@pytest.mark.format
def test_symbol_and_enstid_in_column(data):
    p1 = r'(?:[a-zA-Z0-9-_/.]+)' # only gene symbol
    p2 = r'ENST(?:\d{11})' # only ensembl transcript id
    p3 = p1 + r'\s\(' + p2 + r'\)' # gene symbol (ensembl id)
    p4 = r'ERCC-(?:\d{5})' # ERCC
    p  = '|'.join([p2, p3, p4])
    matches = data.columns.map(lambda x: re.fullmatch(p, x))
    assert  matches.notnull().all(), \
        'some columns do not follow the symbol (ensembl transcript id) format. The first few are: \n{}'.format(data.columns[matches.isnull()][:20])


PARAMS_test_symbol_and_ensgid_in_column = [x['file'] for x in FILE_ATTRIBUTES if x['ismatrix'] and x['gene_id']=='ensg']
@pytest.mark.parametrize('data', PARAMS_test_symbol_and_ensgid_in_column, indirect=True)
@pytest.mark.format
def test_symbol_and_ensgid_in_column(data):
    p1 = r'(?:[a-zA-Z0-9-_/.]+)' # only gene symbol
    p2 = r'ENSG(?:\d{11})' # only ensembl gene id
    p3 = p1 + r'\s\(' + p2 + r'\)' # gene symbol (ensembl id)
    p4 = r'ERCC-(?:\d{5})' # ERCC
    p  = '|'.join([p2, p3, p4])
    matches = data.columns.map(lambda x: re.fullmatch(p, x))
    assert  matches.notnull().all(), \
        'some columns do not follow the symbol (ensembl gene id) format. The first few are: \n{}'.format(data.columns[matches.isnull()][:20])


PARAMS_test_arxspan_ids = [x['file'] for x in FILE_ATTRIBUTES if not x['ismatrix']]
@pytest.mark.parametrize('data', PARAMS_test_arxspan_ids, indirect=True)
@pytest.mark.format
def test_arxspan_ids(data):
    assert 'DepMap_ID' in data.columns
    column = data['DepMap_ID']
    matches = column.map(lambda x: re.match(r'ACH-[\d]{6}$', x))
    assert  matches.notnull().all(), \
        'some rows do not follow the ACH-xxxxxx format. The first few are: \n{}'.format(column.index[matches.isnull()][:20])


PARAMS_test_null_values = [pytest.param(x['file'], marks=pytest.mark.xfail(strict=True)) if x['hasNA'] else x['file']  for x in FILE_ATTRIBUTES if x['ismatrix']]
@pytest.mark.parametrize('data', PARAMS_test_null_values, indirect=True)
@pytest.mark.format
def test_null_values(data):
    assert data.notnull().values.all(), 'null values identified in the matrix'


PARAMS_test_matrix_datatypes = [x['file'] for x in FILE_ATTRIBUTES if x['ismatrix']]
@pytest.mark.parametrize('data', PARAMS_test_matrix_datatypes, indirect=True)
@pytest.mark.format
def test_matrix_datatypes(data):
    datatypes = set(data.dtypes)
    assert len(datatypes) == 1
    assert list(datatypes)[0] == np.dtype('float64')


@pytest.mark.parametrize('data', ['CCLE_expression_full'], indirect=True)
@pytest.mark.format
def test_expression_logtransform(data):
    assert data.min().min() == 0
    assert not data.sum(axis=1).map(lambda x: np.isclose(x, 1e6, rtol=1e-4)).all(), 'expression data is not log-transformed'


@pytest.mark.parametrize('data', ['CCLE_segment_cn', 'CCLE_mutations'], indirect=True)
@pytest.mark.format
def test_chromosome_names(data):
    matches = data['Chromosome'].drop_duplicates().map(lambda x: re.match(r'^\d+|X|Y|M$', x))
    assert matches.notnull().all()


@pytest.mark.format
@pytest.mark.skip(reason='not implemented')
def test_cnv_logtransform():
    # TODO: implement some direct way to test for CNV log transformation
    assert True
