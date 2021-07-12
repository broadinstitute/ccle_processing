# there are some issues in the older versions of omics data on virtual that this flag deals with
# some columns got renamed in the mutation file and some data was uploaded as tsv instead of csv
LEGACY_PATCH_FLAGS = {'rename_column': False, 'tsv2csv': False}

# release ids on taiga
TENTATIVE_VIRTUAL = {'name': 'tentative-virtual-d84e', 'version': 15}
VIRTUAL_RELEASE = TENTATIVE_VIRTUAL # new release
# VIRTUAL_RELEASE = {'name': 'internal-21q2-9d16', 'version': 6} # new release
REFERENCE_RELEASE = {'name': 'internal-21q2-9d16', 'version': 17} # old release used as ground truth
# REFERENCE_RELEASE = {'name': 'internal-21q1-4fc4', 'version': 39} # old release used as ground truth
# VIRTUAL_RELEASE = {'name': 'internal-21q1-4fc4', 'version': 39} # old release used as ground truth
# REFERENCE_RELEASE = TENTATIVE_VIRTUAL # old release used as ground truth
# REFERENCE_RELEASE = {'name': 'internal-20q4-2540', 'version': 47}; LEGACY_PATCH_FLAGS = {'rename_column': True, 'tsv2csv': True} # old release used as ground truth
# REFERENCE_RELEASE = {'name': 'internal-20q1-f1a0', 'version': 15}; LEGACY_PATCH_FLAGS = {'rename_column': True, 'tsv2csv': True} # old release used as ground truth
# REFERENCE_RELEASE = {'name': 'internal-20q2-7f46', 'version': 18} # old release used as ground truth

# these are the columns that if merged with an older release (assuming that old data was not altered),
# should uniquely identify each row of the file to find equal values in each column
FUSIONS_MERGE_COLS = ['DepMap_ID', 'LeftGene', 'RightGene', 'LeftBreakpoint', 'RightBreakpoint']
SEGMENT_CN_MERGE_COLS = ['DepMap_ID', 'Chromosome', 'Start', 'End']
MUTATIONS_MERGE_COLS = ['DepMap_ID', 'Chromosome', 'Start_position', 'End_position', 'Tumor_Seq_Allele1']

# if there are new files that were added in the previous release, add them here
FILES_RELEASED_BEFORE = ['CCLE_expression', 'CCLE_expression_proteincoding_genes_expected_count',
                         'CCLE_RNAseq_transcripts', 'CCLE_expression_transcripts_expected_count',
                         'CCLE_expression_full', 'CCLE_RNAseq_reads', 'CCLE_fusions', 'CCLE_fusions_unfiltered',
                         'CCLE_gene_cn', 'CCLE_segment_cn', 'CCLE_mutations']

# correlation thresholds above which we consider two releases as 'similar'
CORRELATION_THRESHOLDS = {'CCLE_gene_cn': 0.99, 'all_expressions': 0.99999}

SKIP_ARXSPAN_COMPARISON = True # set to False if you want to test whether some arxspans were added/removed

PLOTS_OUTPUT_FILENAME_PREFIX = '/tmp/plots_' # location/prefix for saving output plots

# all the file attributes
FILE_ATTRIBUTES = [
    {'file': 'CCLE_expression', 'ismatrix': True, 'hasNA': False, 'gene_id': 'entrez', 'omicssource':'RNA'},
    {'file': 'CCLE_expression_proteincoding_genes_expected_count', 'ismatrix': True, 'hasNA': False, 'gene_id': 'entrez', 'omicssource':'RNA'},
    {'file': 'CCLE_RNAseq_transcripts', 'ismatrix': True, 'hasNA': False, 'gene_id': 'enst', 'omicssource':'RNA'},
    {'file': 'CCLE_expression_transcripts_expected_count', 'ismatrix': True, 'hasNA': False, 'gene_id': 'enst', 'omicssource':'RNA'},
    {'file': 'CCLE_expression_full', 'ismatrix': True, 'hasNA': False, 'gene_id': 'ensg', 'omicssource':'RNA'},
    {'file': 'CCLE_RNAseq_reads', 'ismatrix': True, 'hasNA': False, 'gene_id': 'ensg', 'omicssource':'RNA'},
    {'file': 'CCLE_fusions', 'ismatrix': False, 'omicssource':'RNA', 'merge_cols': FUSIONS_MERGE_COLS, 'expected_changed_cols':['CCLE_count']},
    {'file': 'CCLE_fusions_unfiltered', 'ismatrix': False, 'omicssource':'RNA', 'merge_cols': FUSIONS_MERGE_COLS, 'expected_changed_cols':['CCLE_count']},
    {'file': 'CCLE_ssGSEA', 'ismatrix': True, 'hasNA': False,'omicssource':'RNA', 'gene_id': None},
    {'file': 'CCLE_gene_cn', 'ismatrix': True, 'hasNA': True, 'gene_id': 'entrez', 'omicssource':'DNA'},
    {'file': 'CCLE_segment_cn', 'ismatrix': False, 'omicssource':'DNA', 'merge_cols': SEGMENT_CN_MERGE_COLS, 'expected_changed_cols':[]},
    {'file': 'CCLE_mutations', 'ismatrix': False, 'omicssource':'DNA', 'merge_cols': MUTATIONS_MERGE_COLS, 'expected_changed_cols':[]},
    {'file': 'CCLE_mutations_bool_damaging', 'ismatrix': True, 'hasNA': False, 'gene_id': 'entrez', 'omicssource':'DNA'},
    {'file': 'CCLE_mutations_bool_hotspot', 'ismatrix': True, 'hasNA': False, 'gene_id': 'entrez', 'omicssource':'DNA'},
    {'file': 'CCLE_mutations_bool_otherconserving', 'ismatrix': True, 'hasNA': False, 'gene_id': 'entrez', 'omicssource':'DNA'},
    {'file': 'CCLE_mutations_bool_nonconserving', 'ismatrix': True, 'hasNA': False, 'gene_id': 'entrez', 'omicssource':'DNA'}
]

# comment/uncomment to use all/subset of files for testing
FILE_ATTRIBUTES = [x for x in FILE_ATTRIBUTES if x['omicssource'] in ['RNA']]

# the following information is used to create a tentative virtual
MUTATIONS_TAIGA_ID = 'mutations-latest-ed72'
FUSIONS_TAIGA_ID = 'fusions-95c9'
EXPRESSION_TAIGA_ID = 'expression-d035'
CN_TAIGA_ID = 'cn-achilles-version-06ca'

TAIGA_IDS_LATEST = {
    MUTATIONS_TAIGA_ID:[
        ('CCLE_mutations', 'all_somatic_mutations_all_21Q2_depmapversion'),
        ('CCLE_mutations_bool_damaging', 'all_somatic_mutations_boolmatrix_fordepmap_damaging'),
        ('CCLE_mutations_bool_nonconserving', 'all_somatic_mutations_boolmatrix_fordepmap_othernoncons'),
        ('CCLE_mutations_bool_otherconserving', 'all_somatic_mutations_boolmatrix_fordepmap_othercons'),
        ('CCLE_mutations_bool_hotspot', 'all_somatic_mutations_boolmatrix_fordepmap_hotspot')
    ],
    FUSIONS_TAIGA_ID:[
        ('CCLE_fusions_unfiltered', 'fusions_latest'),
        ('CCLE_fusions', 'filteredfusions_latest')
    ],
    EXPRESSION_TAIGA_ID:[
        ('CCLE_expression_full', 'genes_tpm_logp1'),
        ('CCLE_RNAseq_transcripts', 'transcripts_tpm_logp1'),
        ('CCLE_RNAseq_reads', 'genes_expected_count'),
        ('CCLE_expression', 'proteincoding_genes_tpm_logp1'),
        ('CCLE_expression_proteincoding_genes_expected_count', 'proteincoding_genes_expected_count'),
        ('CCLE_expression_transcripts_expected_count', 'transcripts_expected_count'),
        ('CCLE_ssGSEA', 'gene_sets_all')
    ],
    CN_TAIGA_ID:[
        ('CCLE_gene_cn', 'all_21Q2_gene_cn'),
        ('CCLE_segment_cn', 'all_21Q2_segment')
    ]
}
