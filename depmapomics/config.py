import numpy as np

########################## GENERIC PARAMS

CACHE_PATH = '~/.depmapomics/'
TMP_PATH = '/tmp/'
ENSEMBL_SERVER_V = "http://nov2020.archive.ensembl.org/biomart"

SHEETCREDS = '../.credentials.json'
MY_ID = '~/.client_secret.json'
MYSTORAGE_ID = "~/.storage.json"

SHEETNAME = 'ccle sample tracker'

GCS_PAYER_PROJECT = 'broad-firecloud-ccle'

TAIGA_ETERNAL = 'depmap-a0ab'

REFSHEET_URL = "https://docs.google.com/spreadsheets/d/1Pgb5fIClGnErEqzxpU7qqX6ULpGTDjvzWwDN8XUJKIY"

PRIVACY_RELEASE = "https://docs.google.com/spreadsheets/d/115TUgA1t_mD32SnWAGpW9OKmJ2W5WYAOs3SuSdedpX4"
DEPMAP_PV = "https://docs.google.com/spreadsheets/d/1uqCOos-T9EMQU7y2ZUw4Nm84opU5fIT1y7jet1vnScE"

POTENTIAL_LIST = "https://docs.google.com/spreadsheets/d/1YuKEgZ1pFKRYzydvncQt9Y_BKToPlHP-oDB-0CAv3gE"

DEPMAP_TAIGA = "arxspan-cell-line-export-f808"

SAMPLEID="DepMap_ID"

SAMPLESETNAME = "21Q3"

RELEASE = SAMPLESETNAME

isCCLE = True

doCleanup = True

LINES_TO_RELEASE = [
    "ACH-000023",
    "ACH-000145",
    "ACH-000205",
    "ACH-000345",
    "ACH-000359",
    "ACH-000399",
    "ACH-000409",
    "ACH-000515",
    "ACH-000664",
    "ACH-000744",
    "ACH-000992",
    "ACH-001032",
    "ACH-001373",
    "ACH-001393",
    "ACH-001410",
    "ACH-001558",
    "ACH-001679",
    "ACH-001696",
    "ACH-002022",
    "ACH-002060",
    "ACH-002659",
    "ACH-002687",
    "ACH-002707",
]

PREV_VIRTUAL = {}

#20Q3
#PREV_VIRTUAL={}
#PREV_VIRTUAL['public'] = 'public-20q3-3d35'
#PREV_VIRTUAL['dmc'] = 'dmc-20q3-deprecated-never-released--5f55'
#PREV_VIRTUAL['internal'] = 'internal-20q3-00d0'

#20Q4
#PREV_VIRTUAL={}
#PREV_VIRTUAL['public'] = 'public-20q4-a4b3'
#PREV_VIRTUAL['dmc'] = 'dmc-20q4-fcf4'
#PREV_VIRTUAL['ibm'] = 'ibm-20q4-269f'
#PREV_VIRTUAL['internal'] = 'internal-20q4-2540'

#21Q1
PREV_VIRTUAL['public'] = 'public-21q1-4b39'
PREV_VIRTUAL['ibm'] = 'ibm-21q1-abd9'
PREV_VIRTUAL['dmc'] = 'dmc-21q1-0e11'
PREV_VIRTUAL['internal'] = 'internal-21q1-4fc4'

CHANGES = """

"""

RUN_NOTEBOOKS = ['WGS_CCLE.ipynb', 'RNA_CCLE.ipynb']

UPLOAD_NOTEBOOK = ['DepMap_Upload.ipynb']

NOTEBOOKS = RUN_NOTEBOOKS+UPLOAD_NOTEBOOK

VIRTUAL_FOLDER = "8d9c4c0691154a1f86b1b6e67c3fb683"

RELEASE = SAMPLESETNAME.lower()

############## TERRA

HG38BAMCOL = ['internal_bam_filepath',
              "internal_bai_filepath"]

############## LOADING

TORAISE=["ACH-001195"],

TO_UPDATE = {'primary_disease': ['Primary Disease'],
            'sex': ['CCLF Gender'],
            'primary_site': ['Sample Collection Site'],
            'subtype': ['lineage_subtype'],
            'subsubtype': ['lineage_sub_subtype'],
            'origin': ['lineage'],
            #'source': ['Flagship'],
            'parent_cell_line': ["Parental ID"],
            'comments': ['Comments'],
            'mediatype': ['Culture Medium', 'Culture Type'],
            'stripped_cell_line_name': ['Stripped Cell Line Name'],
            "cellosaurus_id": ["RRID"],
            "age": ["CCLF Age"]}

MAXAGE = '2021-01-01'

EXTRACT_TO_CHANGE = {'from_arxspan_id': 'participant'}

REPLACE = {'T': 'Tumor', 'N': 'Normal', 'm': 'Unknown', 'L': 'Unknown'}

values = ['legacy_bam_filepath', 'legacy_bai_filepath'],

filetypes = ['bam', 'bai']

MATCH=['ACH-', 'CDS-']

## old GP storage buckets
#rnaworkspace2 = "broad-firecloud-ccle/CCLE_DepMap_RNAseq"
#rnaworkspace4 = "broad-genomics-delivery/Cancer_Cell_Line_Factory_CCLF_RNAseq"
#rnaworkspace5 = "nci-mimoun-bi-org/CCLF_RNA_2_0"
#rnaworkspace3 = "broad-genomics-delivery/CCLE_DepMap_RNAseq"
#rnaworkspace1 = "broad-genomics-delivery/Getz_IBM_CellLines_RNASeqData"

## curent GP buckets
rnaworkspace6 = "terra-broad-cancer-prod/CCLE_DepMap_RNAseq"
rnaworkspace7 = "terra-broad-cancer-prod/Getz_IBM_CellLines_RNASeqData"

## and their correesponding sample source
#rnasource1 = "ibm"
#rnasource2 = "ccle"
#rnasource3 = "ccle"
#rnasource4 = "cclf"
#rnasource5 = "cclf"

rnasource6 = "ccle"
rnasource7 = "ibm"

## our working workspace (reference)
RNAWORKSPACE = "broad-firecloud-ccle/DepMap_hg38_RNAseq"

## curent WGS GP buckets
wgsworkspace1 = "terra-broad-cancer-prod/DepMap_WGS"
wgsworkspace2 = "terra-broad-cancer-prod/Getz_IBM_CellLines_WGS"

## and their corresponding sample source
wgssource1 = "ccle"
wgssource2 = "ibm"

WGSWORKSPACE = "broad-firecloud-ccle/DepMap_WGS_CN"
WESCNWORKSPACE = "broad-firecloud-ccle/DepMap_WES_CN_hg38"
WESMUTWORKSPACE = "broad-firecloud-ccle/DepMap_Mutation_Calling_CGA_pipeline"

############## DNAseq

WGS_methods = [
  "gatk/PreProcessingForVariantDiscovery_GATK4/8",
  "GP-TAG/Manta_SomaticSV/9",
  "gkugener/ArrayOfFilesToTxt/1",
  "vdauwera/BamToUnmappedRGBams/4",
  "gatk/CNV_Somatic_Pair_Workflow/9",
  "gkugener/Aggregate_CN_seg_files/2"
]

CNWES_methods = [
  "gatk/PreProcessingForVariantDiscovery_GATK4/8",
  "GP-TAG/Manta_SomaticSV/9",
  "gkugener/ArrayOfFilesToTxt/1",
  "vdauwera/BamToUnmappedRGBams/4",
  "gatk/CNV_Somatic_Pair_Workflow/9",
  "gkugener/Aggregate_CN_seg_files/2"
]

PROCQC = ["allelic_counts_tumor", "delta_MAD_tumor", "denoised_MAD_tumor",
          "scaled_delta_MAD_tumor", "denoised_copy_ratios_lim_4_plot_tumor",
          "denoised_copy_ratios_plot_tumor", "modeled_segments_plot_tumor",
          "gatk_cnv_all_plots", "lego_plotter_pngs", "copy_number_qc_report",
          "ffpe_OBF_figures", "mut_legos_html", "oxoG_OBF_figures",
          "tumor_bam_base_distribution_by_cycle_metrics",
          "tumor_bam_converted_oxog_metrics"]

BAMQC = ["duplication_metrics", "bqsr_report",
         "tumor_bam_alignment_summary_metrics",
         "tumor_bam_bait_bias_summary_metrics",
         "tumor_bam_gc_bias_summary_metrics",
         "tumor_bam_hybrid_selection_metrics",
         "tumor_bam_insert_size_histogram",
         "tumor_bam_insert_size_metrics",
         "tumor_bam_pre_adapter_summary_metrics",
         "tumor_bam_quality_by_cycle_metrics",
         "tumor_bam_quality_distribution_metrics",
         "tumor_bam_quality_yield_metrics"]

KNOWN_DROP = ['CDS-R22IHj', 'CDS-xMnTwN']

############## CN

COLRENAMING = {'CONTIG': 'Chromosome',
               'START': 'Start',
               'END': 'End',
               'end': 'End',
               'seqnames': 'Chromosome',
               'start': 'Start',
               'Sample': SAMPLEID,
               'NUM_POINTS_COPY_RATIO': 'Num_Probes',
               'MEAN_LOG2_COPY_RATIO': 'Segment_Mean',
               'CALL': 'Status'}

SOURCE_RENAME = {'CCLF': 'Broad', 'CHORDOMA': 'Chordoma',
                'SANGER': 'Sanger', 'IBM': 'Broad', "CCLE2": "Broad",
                np.nan: 'Broad', 'DEPMAP': 'Broad',
                'IBM WES': "Broad WES", 'Broad CCLF': "Broad WES"}

wrongwes_arxspan = {'ACH-001189', 'ACH-002303', 'ACH-002315', 'ACH-002341', 'ACH-001011',
                    'ACH-001108', 'ACH-001187', 'ACH-002875', 'ACH-002874',
                    "ACH-001955",  # chordoma lines
                    "ACH-001956",
                    "ACH-001957"}

toreprocess = ["CDS-C2RlCj", "CDS-8GqFo5"]

############## Mutations

MUTATION_GROUPS = {
  "other conserving": ["5'Flank", "Intron", "IGR", "3'UTR", "5'UTR"],
  "other non-conserving": ["In_Frame_Del", "In_Frame_Ins", "Stop_Codon_Del",
                            "Stop_Codon_Ins", "Missense_Mutation", "Nonstop_Mutation"],
  'silent': ['Silent'],
  "damaging": ['De_novo_Start_OutOfFrame', 'Frame_Shift_Del', 'Frame_Shift_Ins',
                'Splice_Site', 'Start_Codon_Del', 'Start_Codon_Ins', 'Start_Codon_SNP', 'Nonsense_Mutation']
}

MUTCOL_DEPMAP = ['Hugo_Symbol', 'Entrez_Gene_Id', 'NCBI_Build', 'Chromosome',
                 'Start_position', 'End_position', 'Strand', 'Variant_Classification',
                 'Variant_Type', 'Reference_Allele', 'Tumor_Allele', 'dbSNP_RS',
                 'dbSNP_Val_Status', 'Genome_Change', 'Annotation_Transcript',
                 SAMPLEID, 'cDNA_Change', 'Codon_Change', 'Protein_Change', 'isDeleterious',
                 'isTCGAhotspot', 'TCGAhsCnt', 'isCOSMIChotspot', 'COSMIChsCnt',
                 'ExAC_AF', "Variant_annotation", 'CGA_WES_AC', 'HC_AC',
                 'RD_AC', 'RNAseq_AC', 'SangerWES_AC', 'WGS_AC']


############## FUSION

FUSION_COLNAME = ['FusionName', 'JunctionReadCount',
                  'SpanningFragCount', 'SpliceType', 'LeftGene', 'LeftBreakpoint',
                  'RightGene', 'RightBreakpoint', 'LargeAnchorSupport', 'FFPM',
                  'LeftBreakDinuc', 'LeftBreakEntropy', 'RightBreakDinuc',
                  'RightBreakEntropy', 'annots']

FUSION_RED_HERRING = ['GTEx_recurrent', 'DGD_PARALOGS', 'HGNC_GENEFAM',
                      'Greger_Normal', 'Babiceanu_Normal', 'ConjoinG', 'NEIGHBORS']

############## EXPRESSION

## genomic annotations (v35)
GENCODE = 'ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_35/gencode.v35.annotation.gff3.gz'

RNAGSPATH38="gs://cclebams/rnasq_hg38/"

STARBAMCOLTERRA = [ "internal_bam_filepath", 'internal_bai_filepath']

RSEM_TRANSCRIPTS = ['rsem_transcripts_expected_count',
                    'rsem_transcripts_tpm']

RSEMFILENAME_GENE=["genes_tpm",
                   "genes_expected_count"]

RSEMFILENAME_TRANSCRIPTS=["transcripts_tpm", "transcripts_expected_count"]

RSEMFILENAME = RSEMFILENAME_GENE+RSEMFILENAME_TRANSCRIPTS

SSGSEAFILEPATH = "data/genesets/msigdb.v7.2.symbols.gmt"

PATHTOGENEPY = "../"

RNASEQC_THRESHOLDS_LOWQUAL = {'minmapping': 0.85, 'minendmapping': 0.75, 'minefficiency': 0.75,
                              'maxendmismatch': 0.02, 'maxmismatch': 0.02, 'minhighqual': 0.8,
                              'minexon': 0.7, "maxambiguous": 0.05, "maxsplits": 0.1,
                              "maxalt": 0.2, "maxchim": 0.05, "minreads": 20000000,
                              "minlength": 80, "maxgenes": 35000, "mingenes": 12000}


RNASEQC_THRESHOLDS_FAILED = {'minmapping': 0.7, 'minendmapping': 0.66, 'minefficiency': 0.6,
                             'maxendmismatch': 0.02, 'maxmismatch': 0.02, 'minhighqual': 0.7,
                             'minexon': 0.66, "maxambiguous": 0.1, "maxsplits": 0.1,
                             "maxalt": 0.5, "maxchim": 0.2, "minreads": 20000000,
                             "minlength": 80, "maxgenes": 35000, "mingenes": 10000}

PREVIOUS_QC_FAIL = ['CDS-12DTEw', 'CDS-9hv1zM', 'CDS-A6GSeQ', 'CDS-aWlMRt', 'CDS-B1ywOH', 'CDS-BixxtG',
                    'CDS-DRM3l2', 'CDS-jOlYT4', 'CDS-KMhiT9', 'CDS-M6mnMA', 'CDS-pYwECX', 'CDS-v6E624',
                    'CDS-vxTqNJ', 'CDS-YxtmkI',"CDS-fk564T","CDS-kU30H5","CDS-G0F5f5","CDS-ABH0uZ"]

###################### README

#20Q3
#prev_virtual['public'] = 'public-20q3-3d35'
#prev_virtual['dmc'] = 'dmc-20q3-deprecated-never-released--5f55'
#prev_virtual['internal'] = 'internal-20q3-00d0'

#20Q4
# prev_virtual['public'] = 'public-20q4-a4b3'
# prev_virtual['dmc'] = 'dmc-20q4-fcf4'
# prev_virtual['ibm'] = 'ibm-20q4-269f'
# prev_virtual['internal'] = 'internal-20q4-2540'

#21Q1
PREV_VIRTUAL = {
  'public': 'public-21q1-4b39',
  'ibm': 'ibm-21q1-abd9',
  'dmc': 'dmc-21q1-0e11',
  'internal': 'internal-21q1-4fc4',
}

TAIGA_MUTATION = "mutations-latest-ed72"

TAIGA_CN="cn-latest-d8d4"

TAIGA_EXPRESSION="expression-d035"

TAIGA_FUSION = "fusions-95c9"

datasets = ['internal', 'ibm', 'dmc', 'public']

README_folder = "../depmap-release-readmes/"

README_currentfolder = README_folder+"release-"+RELEASE+"/"

README_changes = """


"""

FUSIONreadme = """
# Fusions

PORTAL TEAM SHOULD NOT USE THIS: There are lines here that should not make it even to internal.

/!\ This is the most up to date version of the CCLE CN data.

## Annotations

Description: Gene fusions derived from RNAseq data.

Rows: cell lines, IDs contained in the column DepMap_ID

Unfiltered data contains all output fusions, while the filtered data uses the filters suggested by the star fusion docs. These filters are:
- FFPM > 0.1 -  a cutoff of 0.1 means&nbsp;at least 1 fusion-supporting RNAseq fragment per 10M total reads
- Remove known false positives, such as GTEx recurrent fusions and certain paralogs
- Genes that are next to each other
- Fusions with mitochondrial breakpoints
- Removing fusion involving mitochondrial chromosomes or HLA genes
- Removed common false positive fusions (red herring annotations as described in the STAR-Fusion docs)
- Recurrent fusions observed in CCLE across cell lines (in more than 10% of our samples)
- Removed fusions where SpliceType='INCL_NON_REF_SPLICE' and LargeAnchorSupport='NO_LDAS' and FFPM < 0.1
- FFPM < 0.05
"""


RNAseqreadme = """
# RNAseq

PORTAL TEAM SHOULD NOT USE THIS: There are lines here that should not make it even to internal.

/!\ This is the most up to date version of the CCLE RNA data.

## Annotations:

transcriptions (Transcripts rpkm):

genes (gene rpkm):
__Rows__:
__Columns__:
Counts (gene counts):
__Rows__:
__Columns__:
Gene level CN data:
__Rows__:
__Columns__:
 DepMap cell line IDs
 gene names in the format HGNC\_symbol (Entrez\_ID)
DepMap\_ID, Chromosome, Start, End, Num\_Probes, Segment\_Mean
"""


CNreadme = """
# Copy Number

PORTAL TEAM SHOULD NOT USE THIS: There are lines here that should not make it even to internal.

/!\ This is the most up to date version of the CCLE CN data.

# Notations:

all: everything

allWES: all data comes from the WExomeS samples we posses

allWGS: all data comes from the WGenomeS samples we posses

withreplicates: if we have two different sequencing from a sample, we kept both, see the depmap sample tracker for annotations [https://docs.google.com/spreadsheets/d/1XkZypRuOEXzNLxVk9EOHeWRE98Z8_DBvL4PovyM01FE](https://docs.google.com/spreadsheets/d/1XkZypRuOEXzNLxVk9EOHeWRE98Z8_DBvL4PovyM01FE). this dataset is more geared toward QC or in-depth analysis of a particular cell line.

merged: everything from both WGS and WES

latest: only the latest sequencing versions of the samples were kept


Gene level CN data:

__Rows__: cell line IDs

__Columns__: gene names in the format HGNC\_symbol (Entrez\_ID)

Segment level data:

__Columns__: DepMap\_ID, Chromosome, Start, End, Segment\_Mean, Num\_Probes, Calls"""


Achillesreadme = """
# Copy Number

Combined segment and gene-level CN calls from Broad WES, Sanger WES, and Broad SNP. Relative CN, log2(x+1) transformed.

PORTAL TEAM SHOULD NOT USE THIS: There are lines here that should not make it even to internal. Must use subsetted dataset instead. These data will not make it on the portal starting 19Q1. With the DMC portal, there is new cell line release prioritization as to which lines can be included, so a new taiga dataset will be created containing CN for the portal.

These data are generated for Achilles to pull from to run CERES.


Gene level CN data:

__Rows__: DepMap cell line IDs

__Columns__: gene names in the format HGNC\_symbol (Entrez\_ID)

Segment level data:

__Columns__: DepMap\_ID, Chromosome, Start, End, Num\_Probes, Segment\_Mean"""


Mutationsreadme = """
# Mutations

PORTAL TEAM SHOULD NOT USE THIS: There are lines here that should not make it even to internal.

/!\ This is the most up to date version of the CCLE Mutatios data.
The data is most likely of a better quality that what is on other folder. It is however in beta version as not all changes have either been confirmed or accepted by the DepMap Ops and the DepMap Portal Team.

# Notations:

all: every cell lines we have

WES: all data comes from the WExomeS samples we posses

WGS: all data comes from the WGenomeS samples we posses

withreplicates: if we have two different sequencing from a sample, we kept both, see the depmap sample tracker for annotations [https://docs.google.com/spreadsheets/d/1XkZypRuOEXzNLxVk9EOHeWRE98Z8_DBvL4PovyM01FE](https://docs.google.com/spreadsheets/d/1XkZypRuOEXzNLxVk9EOHeWRE98Z8_DBvL4PovyM01FE). this dataset is more geared toward QC or in-depth analysis of a particular cell line.

merged: everything from both WGS and WES

latest: only the latest sequencing versions of the samples were kept

genes (gene rpkm):
__Rows__:
__Columns__:
Counts (gene counts):
__Rows__:
__Columns__:
Gene level CN data:
__Rows__:
__Columns__:
 DepMap cell line IDs
 gene names in the format HGNC\_symbol (Entrez\_ID)
DepMap\_ID, Chromosome, Start, End, Num\_Probes, Segment\_Mean
 """
