from gsheets import Sheets
from depmapomics.qc.archive import get_cell_line_set as gcl
import pickle
import pandas as pd

lines_to_release_21q1_sheet = 'https://docs.google.com/spreadsheets/d/1YuKEgZ1pFKRYzydvncQt9Y_BKToPlHP-oDB-0CAv3gE/edit#gid=0'

taiga_dict = {'21q1': {'public': ['public-21q1-4b39', 17],
                    'dmc': ['dmc-21q1-0e11', 18],
                    'ibm': ['ibm-21q1-abd9', 20],
                    'internal': ['internal-21q1-4fc4', 21]},
            #  '20q3': {'public': ['public-20q3-3d35', 33],
            #         'dmc': ['dmc-20q3-deprecated-never-released--5f55', 18],
            #         'internal': ['internal-20q3-00d0', 9]},
             '20q4': {'public': ['public-20q4-a4b3', 46],
             'dmc': ['dmc-20q4-fcf4', 38],
             'internal': ['internal-20q4-2540', 45]}
                    }


def main(recompute=False, lines_to_remove={'ACH-001189', 'ACH-002303', 'ACH-002315', 'ACH-002341', 'ACH-002359'}):

  if recompute:
      taiga_dict_expanded = gcl.propagate_taiga_dict_with_filenames(taiga_dict)
      arxspan_dict = gcl.get_all_arxspans(taiga_dict_expanded, verbose=True)
      # arxspan_dict['20q3']['ibm'] = arxspan_dict['20q3']['dmc']
      arxspan_dict['20q4']['ibm'] = arxspan_dict['20q4']['dmc']
      pickle.dump(arxspan_dict, open('arxspan_dict.pkl', 'bw'))
  else:
      arxspan_dict = pickle.load(open('arxspan_dict.pkl', 'br'))

  
  # lines_to_remove = set()
  lines_to_release_21q1 = gcl.get_expected_lines(lines_to_release_21q1_sheet)
  text = gcl.pretty_print_diff(arxspan_dict, lines_to_release_21q1, lines_to_remove=lines_to_remove, quarters = ['20q4', '21q1'], savefile=True)

  gcl.check_if_fusion_has_expression_released(arxspan_dict, quarter='21q1')
  gcl.check_acciddental_release(arxspan_dict, quarter = '21q1')
  gcl.plot_diff_heatmap(arxspan_dict, lines_to_release_21q1, lines_to_remove=lines_to_remove, quarters = ['20q4', '21q1'])


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



