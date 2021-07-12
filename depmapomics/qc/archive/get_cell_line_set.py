from gsheets import Sheets
import pandas as pd
from taigapy import TaigaClient
import seaborn as sns
import matplotlib.pyplot as plt

tc = TaigaClient()

DEFAULT_FILENAMES=['CCLE_expression', 'CCLE_expression_full', 'CCLE_fusions',
                   'CCLE_fusions_unfiltered', 'CCLE_gene_cn', 'CCLE_mutations',
                   'CCLE_RNAseq_reads', 'CCLE_RNAseq_transcripts', 'CCLE_segment_cn']

def tcget(name=None, version=None, file=None):
    # get arxspan ids from the Taiga file
    # this function avoids overloading the memory because the entire variable is local to it
    tc_df = tc.get(name=name, version=version, file=file)

    if 'DepMap_ID' in tc_df.columns:
        arxspan_ids = set(tc_df['DepMap_ID'])
    else:
        arxspan_ids = set(tc_df.index)
    assert all([x.startswith('ACH-') for x in arxspan_ids])
    return arxspan_ids

def applyfunc_to_json(json_dict, func,  verbose=False):
    arxspans = {}
    if type(json_dict)==dict:
        output_dict = {}
        for k, v in json_dict.items():
            if verbose:
                print(k)
            output_dict[k] = applyfunc_to_json(v, func, verbose=verbose)
        return output_dict
    else:
        return func(json_dict)

def get_expected_lines(sheets_url, merge_ibm_dmc=True):
    sheets_obj = Sheets.from_files('~/.client_secret.json', '~/.storage.json')
    sheets = sheets_obj.get(sheets_url).sheets
    release = sheets[0].to_frame(header=0, index_col=None)
    release.columns = release.columns.str.lower()
    if merge_ibm_dmc:
        cols = release.columns
        ibm = pd.DataFrame(set(release[['dmc', 'ibm']].stack().values), columns=['ibm'])
        # ibm = pd.concat([release['dmc'], release['ibm']]
        #         ).dropna().drop_duplicates().to_frame('ibm')
        release = pd.concat([release[['internal', 'dmc', 'public']], ibm], axis=1)
        lines_to_release = release[cols]
    return release

def get_release_diffs(arxspan_dict, lines_to_release, lines_to_remove=set(), quarters=None):
    release_diff = {}
    release_diff_reverse = {}
    for portal, portal_dict in arxspan_dict[quarters[0]].items():
        release_diff[portal] = {}
        release_diff_reverse[portal] = {}
        for file in portal_dict.keys():
            expected_lines = set(lines_to_release[portal].dropna()) | arxspan_dict[quarters[0]][portal][file]
            expected_lines = expected_lines-lines_to_remove

            release_diff[portal][file] = arxspan_dict[quarters[1]][portal][file] - expected_lines
            release_diff_reverse[portal][file] = expected_lines - arxspan_dict[quarters[1]][portal][file]
    return release_diff, release_diff_reverse

def pool_arxspans_per_portal(arxspan_dict):
    arxspan_dict_per_portal = {}
    for quarter, qurter_dict in arxspan_dict.items():
        portal_arxspans_dict = {}
        for portal, portal_dict in qurter_dict.items():
            portal_arxspans = set()
            for file in portal_dict.keys():
                portal_arxspans = portal_arxspans | arxspan_dict[quarter][portal][file]
            portal_arxspans_dict[portal] = portal_arxspans
        arxspan_dict_per_portal[quarter] = portal_arxspans_dict
    return arxspan_dict_per_portal


def get_release_diff_pooled(arxspan_dict, lines_to_release, lines_to_remove=set(), quarters=None):
    arxspans = pool_arxspans_per_portal(arxspan_dict)
    arxspan_diff = {}
    arxspan_revdiff = {}
    for portal in arxspans[quarters[0]].keys():
        arxspans_expected = (arxspans[quarters[0]][portal] | set(lines_to_release[portal].dropna()))
        arxspans_expected = arxspans_expected-lines_to_remove
        arxspan_diff[portal] = arxspans[quarters[1]][portal] - arxspans_expected
        arxspan_revdiff[portal] = arxspans_expected - arxspans[quarters[1]][portal]
    return arxspan_diff, arxspan_revdiff

def get_all_arxspans(taiga_dict_expanded, verbose=False):
    arxspan_dict = applyfunc_to_json(taiga_dict_expanded, lambda x: tcget(name=x[0], version=x[1], file=x[2]), verbose=True)
    return arxspan_dict

def propagate_taiga_dict_with_filenames(taiga_dict, files=DEFAULT_FILENAMES):
    return applyfunc_to_json(taiga_dict, lambda x:{file:x+[file] for file in files})

def pretty_print(release_diffs):
    def _convert_to_text(x):
        if x == set():
            return 'None'
        else:
            return ', '.join(x)
    text = ''
    release_diffs_text = applyfunc_to_json(release_diffs, lambda x: _convert_to_text(x))
    for portal, portal_dict in release_diffs_text.items():
        text += portal + ':\n'
        if type(portal_dict) == dict:
            for file in portal_dict.keys():
                text+= '\t{}: {}\n'.format(file, release_diffs_text[portal][file])
            text+= '\n'
        else:
            text+= '\t{}\n'.format(release_diffs_text[portal])
    return(text)

def convert_diff_to_boolmatrix(release_diffs):
    diff_df = pd.DataFrame(release_diffs)
    diff_df = diff_df.stack()
    arxspans = list(set(diff_df.map(list).sum()))
    diff_df = diff_df.to_frame().apply(lambda arx_set: [x in arx_set[0] for x in arxspans], result_type='expand', axis=1)
    diff_df.columns = arxspans
    diff_df = diff_df.T.astype(int)
    return diff_df

def get_release_diff_boolmatrices(arxspan_dict, lines_to_release, lines_to_remove=set(), quarters = None):
    release_diffs, release_diffs_reverse = get_release_diffs(arxspan_dict, lines_to_release, lines_to_remove=lines_to_remove, quarters = quarters)
    release_diffs_mat = convert_diff_to_boolmatrix(release_diffs)
    release_diffs_reverse_mat = convert_diff_to_boolmatrix(release_diffs_reverse)
    return release_diffs_mat, release_diffs_reverse_mat

def check_if_fusion_has_expression_released(arxspan_dict, quarter=None):
    for portal in ['public', 'dmc', 'ibm', 'internal']:
        assert (arxspan_dict[quarter][portal]['CCLE_fusions_unfiltered'] - arxspan_dict[quarter][portal]['CCLE_expression']) == set()
        assert (arxspan_dict[quarter][portal]['CCLE_fusions'] - arxspan_dict[quarter][portal]['CCLE_expression']) == set()
    print('all fusion files have accompanying CCLE_expression')


def check_acciddental_release(arxspan_dict, quarter = None,
    portal_comparisons = [['dmc', 'ibm'], ['public', 'dmc'],
                        ['internal', 'ibm'], ['internal', 'dmc']]):
    for portal_pairs in portal_comparisons:
        nodiff=True
        for file in DEFAULT_FILENAMES:
            diff = arxspan_dict[quarter][portal_pairs[0]][file] - arxspan_dict[quarter][portal_pairs[1]][file]
            if diff!=set():
                nodiff=False
                print('lines in {} but not {} for {}: {}'.format(*portal_pairs, file, diff))
        if nodiff:
            print('all lines in {} are also in {}'.format(*portal_pairs))
        print('\n')


def plot_diff_heatmap(arxspan_dict, lines_to_release, lines_to_remove=set(), quarters = None, width=15, height_scale=2):
    release_diffs_mat, release_diffs_reverse_mat = get_release_diff_boolmatrices(arxspan_dict, lines_to_release, lines_to_remove=lines_to_remove, quarters = quarters)

    plt.figure(figsize=(width, int(release_diffs_mat.shape[0]/height_scale+1)))
    sns.heatmap(release_diffs_mat, cbar=False)
    plt.title('{} vs {}: files added (white: true)'.format(quarters[1], quarters[0]));
    plt.savefig('{}-vs-{}-files-added.png'.format(quarters[1], quarters[0]), bbox_inches='tight')

    plt.figure(figsize=(width, int(release_diffs_reverse_mat.shape[0]/height_scale+1)))
    sns.heatmap(release_diffs_reverse_mat, cbar=False)
    plt.title('{} vs {}: files dropped (white: true)'.format(quarters[1], quarters[0]));
    plt.savefig('{}-vs-{}-files-dropped.png'.format(quarters[1], quarters[0]), bbox_inches='tight')

def pretty_print_diff(arxspan_dict, lines_to_release, lines_to_remove=set(), quarters = None, savefile=False):
    release_diffs, release_diffs_reverse = get_release_diffs(arxspan_dict, lines_to_release, lines_to_remove=lines_to_remove, quarters = quarters)
    release_diffs_pooled, release_diffs_reverse_pooled = get_release_diff_pooled(arxspan_dict, lines_to_release, lines_to_remove=lines_to_remove, quarters = quarters)
    text= 'lines added ({} compared to {})\n'.format(*quarters[::-1])
    text+= pretty_print(release_diffs)
    text += '\n\tlines pooled across portals\n'
    text += pretty_print(release_diffs_pooled)

    text += '\n'+'_'*20 + '\n'

    text += 'lines dropped ({} compared to {})\n'.format(*quarters[::-1])
    text += pretty_print(release_diffs_reverse)
    text += '\n\tlines pooled across portals\n'
    text += pretty_print(release_diffs_reverse_pooled)
    if savefile:
        with open('diff_{}_vs_{}.txt'.format(*quarters[::-1]), 'w') as f:
            f.write(text)

    return text
