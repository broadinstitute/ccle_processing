from gsheets import Sheets
from depmapomics.tests.archive import get_cell_line_set as gcl
import pickle
import pandas as pd

lines_to_release = 'https://docs.google.com/spreadsheets/d/1-MrXF2O-255vUPb5vbPO3cSJ1Cfp8yibTL88KtSwfxM/edit#gid=0'


taiga_dict = {'21q1': {'public': ['public-21q1-4b39', 33],
                    'dmc': ['dmc-21q1-0e11', 30],
                    'ibm': ['ibm-21q1-abd9', 26],
                    'internal': ['internal-21q1-4fc4', 39]},
              '21q2': {'public': ['public-21q2-110d', 9],
                    'dmc': ['dmc-21q2-27e1', 10],
                    'ibm': ['ibm-21q2-9ed1', 12],
                    'internal': ['internal-21q2-9d16', 14]}
                    }
quarters = ['21q1', '21q2']

recompute = False
if recompute:
    taiga_dict_expanded = gcl.propagate_taiga_dict_with_filenames(taiga_dict)
    arxspan_dict = gcl.get_all_arxspans(taiga_dict_expanded, verbose=True)
    pickle.dump(arxspan_dict, open('arxspan_dict.pkl', 'bw'))
else:
    arxspan_dict = pickle.load(open('arxspan_dict.pkl', 'br'))

lines_to_remove = {'ACH-002010', 'ACH-000314'}
# lines_to_remove = set()
lines_to_release_21q1 = gcl.get_expected_lines(lines_to_release)
text = gcl.pretty_print_diff(arxspan_dict, lines_to_release_21q1, lines_to_remove=lines_to_remove, quarters = quarters, savefile=True)

gcl.check_if_fusion_has_expression_released(arxspan_dict, quarter=quarters[1])
gcl.check_acciddental_release(arxspan_dict, quarter = quarters[1])
gcl.plot_diff_heatmap(arxspan_dict, lines_to_release_21q1, lines_to_remove=lines_to_remove, quarters = quarters)


# lines_to_drop_dict = {
#     'ACH-002217': "no bam found: Sanger",
#     'ACH-002335': "no bam found: Chordoma",
#     'ACH-002378': "no bam found: Sanger",
#     'ACH-001956': "failed qc wes (chordoma)",
#     'ACH-001955': "failed qc wes (chordoma)",
#     'ACH-001957': "failed qc wes (chordoma)",
#     'ACH-001011': 'wrong new (should only have HC)',
#     'ACH-001108': 'wrong new (should only have HC)',
#     'ACH-001187': 'wrong new (should only have HC)',
#     'ACH-001189': 'wrong (from 20Q4)',
#     'ACH-002303': 'wrong (from 20Q4)',
#     'ACH-002315': 'wrong (from 20Q4)',
#     'ACH-001675': 'wrong (from 20Q4)',
#     'ACH-003000': 'wrong engineered',
#     'ACH-002875': 'wrong engineered',
#     'ACH-002874': 'wrong engineered'
# }
