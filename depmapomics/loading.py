# -*- coding: utf-8 -*-
# Jérémie Kalfon
# for BroadInsitute
# in 2019


####
#
# HELPER FUNC  ######################################
#
#
from gsheets import Sheets
from genepy.google.google_sheet import dfToSheet
import pandas as pd
import numpy as np
import dalmatian as dm
from taigapy import TaigaClient
tc = TaigaClient()
from depmapomics import tracker
from depmapomics.config import *
from depmapomics import terra as myterra
from genepy import terra
from genepy import sequencing as seq
from genepy.utils import helper as h
from genepy.google import gcp

#####################
# Const Variables
#####################

# chromosome list
CHROMLIST = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8',
             'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15',
             'chr16', 'chr17', 'chr18', 'chr19', 'chr20',
             'chr21', 'chr22', 'chrX']

# default values in the GP workspaces and our sample tracker (to change if you use another workspace/
# sample tracker)
extract_defaults = {
    'name': 'sample_alias',
    'bai': 'crai_or_bai_path',
    'bam': 'cram_or_bam_path',
    'ref_bam': 'legacy_bam_filepath',
    'ref_type': 'datatype',
    'ref_bai': "legacy_bai_filepath",
    'version': 'version',
    'primary_disease': 'primary_disease',
    'ref_arxspan_id': 'arxspan_id',
    'ref_name': 'stripped_cell_line_name',
    'source': 'source',
    'size': 'size',
    "legacy_size":"legacy_size",
    'from_arxspan_id': 'individual_alias',
    'ref_id': 'sample_id',
    'PDO_id':'PDO',
    "update_time":"update_time",
    'from_patient_id': 'individual_alias',
    'patient_id': 'participant_id',
    'ref_date': 'date_sequenced',
    'hs_hs_library_size': 'hs_hs_library_size',
    'hs_het_snp_sensitivity': 'hs_het_snp_sensitivity',
    'hs_mean_bait_coverage': 'hs_mean_bait_coverage',
    'hs_mean_target_coverage': 'hs_mean_target_coverage',
    'hs_on_target_bases': 'hs_on_target_bases',
    'total_reads': 'total_reads',
    'release_date': 'sequencing_date',
    'hash': 'crc32c_hash',
    'legacy_hash': 'legacy_crc32c_hash',
    'mean_depth': 'mean_depth'
}

# minimum bam file size in bytes for each sequencing type
MINSIZES = {
    'rna': 2000000000,
    'wes': 3000000000,
    'wgs': 50000000000,
}

# known cell lines that are from the same patient
samepatient = [["ACH-000635", "ACH-000717", "ACH-000864", "ACH-001042", "ACH-001547"],
               ["ACH-002291", "ACH-001672"],
               ["ACH-001706", "ACH-001707"]]

# known duplicate arxspan-ids
dup = {"ACH-001620": "ACH-001605",
       "ACH-001621": "ACH-001606"}

# rename ccle_name TODO: ask becky what to do
rename = {"PEDS117": "CCLFPEDS0009T"}


#####################
# Loading Functions
#####################

def GetNewCellLinesFromWorkspaces(wmfroms, sources, stype, maxage, refurl="",
                                  addonly=[], match='ACH', extract={},
                                  extract_defaults=extract_defaults, wto=None, refsamples=None,
                                  participantslicepos=10, accept_unknowntypes=False,
                                  rename=dict(), recomputehash=False, my_id = MY_ID,
                                  mystorage_id = MYSTORAGE_ID,):
  """
  As GP almost always upload their data to a data workspace. we have to merge it to our processing workspace

  Will merge samples from a set of data workspaces to a processing workspace on Terra. Will only
  get a subset of the metadata and rename it.
  Will find out the duplicates based on the file size.
  Can also upload the bam files to a google storage bucket

  Args:
  -----
    wto: str the workspace where you want to create the tsvs
    wfroms: list[str] the workspaces where the Samples to add are stored
    sources: list[str] the corresponding source names
    stype: str sequencing type
    maxage: str earliest date of the bam file upload to be considered new
    refurl: str(url) the reference url for the cell line tracker spreadsheet (only if no refsamples)
    match: list[str]|str the possible values that a sample id need to contain to be considered valid
    refsamples: pdDataFrame with columns matching values is in "extract" for the right keys (see "extract_default")
    participantslicepos: int the length of the sample id string
    accept_unknowntypes: bool whether or not the sample type column for that sample can be different from "Tumor"
    rename: dict(str:str) mapping a wrong arxpand_id to a good arxspan id for known cases of misslabelling
    recomputehash: bool whether or not to recompute the hash of the bam file when loading it
    addonly: list of sample id that you only want to add
    extract: if you want to specify what values should refer to which column names
      dict{
      'name':
      'bai':
      'bam':
      'source':
      'from_arxspan_id':
      ...} (see extract_defaults)
    extract_defaults: the full default dict to specificy what values should refer to which column names

  Returns:
  -------
    samples: a dataframe with the samples that were resolved by the tool (we still need to add some more annotations)
    pairs: the corresponding pair from matching known normals with known tumors
    wrongssamples: a dataframe containing samples that passed most QCs but couldn't be resolved

  Raise:
  -----
    Exception: when no new samples in this matrix
  """
  extract.update(extract_defaults)
  if type(match) is str and match:
    match = [match]
  if refurl:
    print('refsamples is overrided by a refurl')
    sheets = Sheets.from_files(my_id, mystorage_id)
    refsamples = sheets.get(refurl).sheets[0].to_frame(index_col=0)
  if refsamples is None:
    if wto is None:
      raise ValueError('missing refsamples or refworkspace (wto)')
    wto = dm.WorkspaceManager(wto)
    print('we do not have refsamples data. Using the wto workspace sample data instead')
    refsamples = wto.get_samples()
    # TODO: update directly the df if data is not already in here)
    refsamples[extract['ref_arxspan_id']] = [a.split('_')[0] for a in refsamples[extract['ref_arxspan_id']] if type(a) is str]
    if extract['hash'] not in refsamples.columns:
      refsamples[extract['hash']] = [gcp.extractHash(val) for val in gcp.lsFiles(
          [i for i in refsamples[extract["ref_bams"]] if type(i) is str and str(i) != 'NA'], "-L", 200)]
    if extract['size'] not in refsamples.columns:
      refsamples['size'] = [gcp.extractSize(i)[1] for i in gcp.lsFiles(refsamples[extract['bam']].tolist(), '-al', 200)]
    if extract['release_date'] not in refsamples.columns:
      refsamples[extract["ref_bams"]] = seq.getBamDate(refsamples[extract["ref_bams"]])
  refsamples[extract['release_date']] = list(h.datetoint(refsamples[extract["release_date"]].values, '/'))
  if stype not in set(refsamples[extract['ref_type']]):
    h.ask("we have never seen this type: " + stype + ", in the reference, continue?")
  # do NOT make refids a set; we use the num of occurences as way to determine what number to add to the sample id
  # filter refids to only include those that include the strings in the 'match' argument
  refsamples = refsamples[refsamples.index.str.contains('|'.join(match))]
  for match_substring in match:
    refsamples.index = [match_substring + i.split(match_substring)[-1] if match_substring in i else i for i in refsamples.index]
  refsamples.index = [i[:participantslicepos] for i in refsamples.index]
  print("Getting sample infos...")
  if type(sources) is str:
    sources = [sources]
  if type(wmfroms) is str:
    wmfroms = [wmfroms]
  sampless = pd.DataFrame()
  wrongsampless = pd.DataFrame()
  for source, wmfrom in zip(sources, wmfroms):
    broken_bams = []
    wmfrom = dm.WorkspaceManager(wmfrom)
    samples = wmfrom.get_samples().replace(np.nan, '', regex=True).reset_index()
    # keep samples that contain the match requirement (e.g. ACH for DepMap IDs)

    print("\nThe shape of the sample tsv from " + str(wmfrom) + ": " + str(samples.shape))

    # remove true duplicates from consideration
    print("Identifying any true duplicates by checking file hashes (this runs for each data source)...")
    print("This step can take a while as we need to use gsutil to check the size of each potential duplicate...")
    dups_to_remove = []
    # check for broken bam files; if broken, then remove from consideration
    # need to check for broken filepaths before checking if the sample is in Terra so that we don't
    # add a broken file path for a new participant
    foundfiles = gcp.lsFiles(samples[extract['bam']])
    broken_bams = set(samples[extract['bam']]) - set(foundfiles)
    print('These ' + str(len(broken_bams)) + ' bam file path do not exist: ' + str(broken_bams))

    wrongsamples = samples[(~samples[extract['bam']].isin(broken_bams)) & (~samples[extract['from_arxspan_id']].str.contains('|'.join(match)))]
    wrongsamples = extractFromWorkspace(wrongsamples, stype, recomputehash, extract)
    if wrongsamples is not None:
      wrongsamples = mapSamples(wrongsamples, source, extract)
      wrongsampless = pd.concat([wrongsampless, wrongsamples], sort=False)
    samples = samples[(~samples[extract['bam']].isin(broken_bams)) & (samples[extract['from_arxspan_id']].str.contains('|'.join(match)))]
    # getting correct arxspan id
    if samples is None:
      continue
    samples = extractFromWorkspace(samples, stype, recomputehash, extract)
    if samples is None:
      continue
    samples = mapSamples(samples, source, extract)
    samples = resolveFromWorkspace(samples, refsamples[refsamples[extract['ref_type']] == stype], match, participantslicepos,
                                   accept_unknowntypes, addonly, extract)
    if samples is None:
      continue
    sampless = pd.concat([sampless, samples], sort=False)

  if len(sampless) == 0:
    print("no new data available")
    return sampless, pd.DataFrame()

  sampless = assessAllSamples(sampless, refsamples, stype, rename, extract)
  # creating pairs
  pairs = myterra.setupPairsFromSamples(sampless, refsamples[refsamples[extract['ref_type']] == stype], extract)
  # I am trying to remove duplicates from samples without arxspan ids to then look more into them
  # and see if I have to get data for them or if I should just throw them out
  toremov = set()
  for k, val in wrongsampless.iterrows():
    withsamesize = wrongsampless[wrongsampless[extract["legacy_size"]] == val[extract["legacy_size"]]]
    if (val[extract["legacy_size"]] in sampless[extract["legacy_size"]].tolist()) or (val[extract["legacy_size"]] in refsamples[extract["size"]]):
      toremov.add(k)
    if len(withsamesize) > 1:
      for l, _ in withsamesize.iloc[1:].iterrows():
        toremov.add(l)
    #elif len(refsamples[refsamples[extract['size']] == withsamesize[extract["size"]][0]]):
      #toremov.add(k)
  for i in toremov:
    wrongsampless = wrongsampless.drop(i)
  for i, v in wrongsampless.iterrows():
    if not gcp.exists(v[extract['ref_bam']]):
      print(v.ccle_name)
      wrongsampless = wrongsampless.drop(i)
  a = len(sampless)
  sampless = deleteClosest(sampless,refsamples, extract['legacy_size'], extract['legacy_size'], extract['ref_arxspan_id'])
  sampless = deleteClosest(sampless,refsamples, extract['legacy_size'], extract['legacy_size'], extract['ref_arxspan_id'])
  print('removed: '+str(a-len(sampless))+" samples from size alone (too similar to a replicate)")
  wrongsampless = wrongsampless[~wrongsampless[extract['legacy_size']].isin(set(refsamples[extract['legacy_size']]))]
  wrongsampless = wrongsampless[~wrongsampless[extract['legacy_size']].isin(set(refsamples[extract['legacy_size']]))]
  wrongsampless = deleteClosest(
      wrongsampless, refsamples, extract['legacy_size'], extract['legacy_size'], extract['ref_arxspan_id'])
  wrongsampless = deleteClosest(
      wrongsampless, refsamples, extract['legacy_size'], extract['legacy_size'], extract['ref_arxspan_id'])
  #removing duplicate PDOs
  a = len(sampless)
  wrongsampless = wrongsampless[~wrongsampless[extract['PDO_id']].isin(set(refsamples[extract['PDO_id']]))]
  sampless = sampless[~sampless[extract['PDO_id']].isin(
      set(refsamples[extract['PDO_id']]))]
  print('removed: '+str(a-len(sampless)) +
        " samples with duplicat PDO ids ")
  # removing anything too old
  a = len(sampless)
  wrongsampless = wrongsampless[wrongsampless[extract['update_time']] > maxage]
  sampless = sampless[sampless[extract['update_time']]>maxage]
  print('removed: '+str(a-len(sampless))+" samples that have not changed since last time (likely\
     duplicate having been removed)")
  return sampless, pairs, wrongsampless


def deleteClosest(sampless, refsamples, size='legacy_size', ref_size='legacy_size', arxspid='arxspan_id'):
  """
  for a list of samples and a tracker, will find the index of the sample with the closest size

  if this sample is the same cell line, it will judge it to be a duplicate and remove it

  Args:
  -----
    sampless: pd dataframes of samples with at least arxspan ids and sizes
    refsamples: pd dataframe representing a sample tracker
    size: str colname of size in the sample list
    ref_size: str colname of size in the sample tracker
    arxspid: str colnme of sample ids
    extract: if you want to specify what values should refer to which column names
      dict{
      'name':
      'bai':
      'bam':
      'source':
      'from_arxspan_id':
      ...} (see extract_defaults)

  Returns:
  --------
    samples: pd dataframe the filtered sample list
  """
  sizes = refsamples[ref_size].tolist()
  for k,v in sampless.iterrows():
    if type(v[size]) is int:
      val = refsamples.iloc[sizes.index(h.closest(sizes, v[size]))]
      if val[arxspid] == v[arxspid]:
        sampless = sampless.drop(v.name)
  return sampless


def extractFromWorkspace(samples, stype, recomputeTime=True, recomputesize=True, 
                         recomputedate=True, recompute_hash=True, extract={}):
  """
  Extract more information from a list of samples found on GP workspaces

  Args:
  -----
    samples: pd dataframes of samples with at least arxspan ids and sizes
    stype: str sequencing type
    recomputeTime: bool whether to recompute the date of upload of the bam file
    recomputesize: bool whether to recompute the of the bam file
    recomputehash: bool whether to recompute the of the bam file
    extract: if you want to specify what values should refer to which column names
      dict{
      'name':
      'bai':
      'bam':
      'source':
      'from_arxspan_id':
      ...} (see extract_defaults)

  Returns:
  --------
    samples: pd dataframe the filtered sample list
  """
  extract.update(extract_defaults)
  if extract['legacy_hash'] not in samples.columns or recompute_hash:
    samples[extract['hash']] = [gcp.extractHash(val) for val in gcp.lsFiles(samples[extract["bam"]].tolist(), "-L", 200)]
  lis = gcp.lsFiles(samples[extract['bam']].tolist(), '-al', 200)
  if extract['legacy_size'] not in samples.columns or recomputesize:
    samples[extract['legacy_size']] = [gcp.extractSize(i)[1] for i in lis]
  if extract['update_time'] not in samples.columns or recomputeTime:
    samples[extract['update_time']] = [gcp.extractTime(i) for i in lis]
  todrop = []
  for k, val in samples.iterrows():
    if val[extract['legacy_size']] < MINSIZES[stype]:
      todrop.append(k)
      print("too small size, removing sample: " + str(val[extract["from_arxspan_id"]]))
  samples = samples.drop(index=todrop)
  # getting the date released
  if len(samples) == 0:
    return None
  if extract['release_date'] not in samples.columns or recomputedate:
    samples[extract["release_date"]] = seq.getBamDate(samples[extract["bam"]])
  samples[extract['release_date']] = list(h.datetoint(samples[extract['release_date']].values))
  return samples


def mapSamples(samples, source, extract={}):
  """
  Convert samples from a list of GP workspaces to something being able to be merged with the sample tracker

  Args:
  -----
    samples: pd dataframes of samples with at least arxspan ids and sizes
    extract: if you want to specify what values should refer to which column names
      dict{
      'name':
      'bai':
      'bam':
      'source':
      'from_arxspan_id':
      ...} (see extract_defaults)
    source:

  Returns:
  --------
    samples: pd dataframe the filtered sample list
  """
  # creating unique ids
  samples[extract['ref_id']] = ['CDS-' + h.randomString(stringLength=6, stype='all', withdigits=True) for _ in range(len(samples))]
  samples[extract['patient_id']] = ['PT-' + h.randomString(stringLength=8, stype='all', withdigits=True) for _ in range(len(samples))]
  samples.reset_index(drop=True, inplace=True)
  samples[extract['source']] = source

  # renamings
  samples = samples.rename(columns={extract['bam']: extract['ref_bam'],
                                    extract['bai']: extract['ref_bai'],
                                    extract['name']: extract['ref_name'],
                                    extract['from_arxspan_id']: extract['ref_arxspan_id']
                                    }).set_index(extract["ref_id"], drop=True)
  # subsetting
  print(samples.columns)
  samples = samples[[extract['ref_bam'], extract['ref_bai'], extract['ref_name'], extract["ref_arxspan_id"], extract["release_date"], extract["patient_id"], extract['legacy_size'], extract['PDO_id'], extract['update_time']]]
  return samples


def resolveFromWorkspace(samples, refsamples, match, participantslicepos=10, accept_unknowntypes=True, addonly=[], extract={}):
  """
  Filters our list by trying to find duplicate in our dataset and remove any sample that isn't tumor

  Args:
  -----
    match: list[str]|str the possible values that a sample id need to contain to be considered valid
    participantslicepos: int the length of the sample id string
    accept_unknowntypes: bool whether or not the sample type column for that sample can be different from "Tumor"
    refsamples: pd dataframe representing a sample tracker
    samples: pd dataframes of samples with at least arxspan ids and sizes
    extract: if you want to specify what values should refer to which column names
      dict{
      'name':
      'bai':
      'bam':
      'source':
      'from_arxspan_id':
      ...} (see extract_defaults)

  Returns:
  --------
    samples: pd dataframe the filtered sample list
  """
  extract.update(extract_defaults)
  prevlen = len(samples)
  for match_substring in match:
    samples[extract['ref_arxspan_id']] = [(match_substring + i.split(match_substring)[-1]) if match_substring in i else i for i in samples[extract['ref_arxspan_id']]]
  samples[extract['ref_arxspan_id']] = [i[:participantslicepos] for i in samples[extract['ref_arxspan_id']]]
  print('we found and removed ' + str(prevlen - len(samples)) + ' samples which did not match our id names: ' + str(match))

  tolookfor = [val[extract['ref_bam']] for _, val in samples.iterrows() if val[extract['ref_arxspan_id']] in set(refsamples[extract['ref_arxspan_id']])]
  print("found " + str(len(tolookfor)) + ' likely replicate')
  sample_hash = {gcp.extractSize(val)[1]: gcp.extractSize(val)[0] for val in gcp.lsFiles(tolookfor, "-la")}
  dups_to_remove = [sample_hash[a] for a in set(sample_hash.keys()) & set(refsamples[extract['legacy_size']])]
  dups_to_remove.extend([sample_hash[a] for a in set(sample_hash.keys()) & set(refsamples[extract['legacy_size']])])
  # remove the duplicates from consideration
  print("Len of samples before removal: " + str(len(samples)))
  print("Dups from this workspace has len " + str(len(dups_to_remove)) + ":\n " + str(dups_to_remove))
  # remove the samples with broken bam filepaths from consideration
  samples = samples[~samples[extract['ref_bam']].isin(dups_to_remove)]

  print("Len of samples after removal: " + str(len(samples)))
  if len(samples) == 0:
    return None

  # if only add some samples
  if len(addonly) > 0:
    samples = samples[samples[extract['ref_arxspan_id']].isin(addonly)]

  # unknown types
  if 'sample_type' in samples.columns:
    if not accept_unknowntypes:
      samples = samples[samples['sample_type'].isin(['Tumor'])]
  return samples


def assessAllSamples(sampless, refsamples, stype, rename={}, extract={}):
  """
  Will look for matching lines and duplicates in our sample tracker and compute version and patient information

  Args:
  -----
    refsamples: pd dataframe representing a sample tracker
    stype: str sequencing type
    rename: dict(str:str) mapping a wrong arxpand_id to a good arxspan id for known cases of misslabelling
    samples: pd dataframes of samples with at least arxspan ids and sizes
    extract: if you want to specify what values should refer to which column names
      dict{
      'name':
      'bai':
      'bam':
      'source':
      'from_arxspan_id':
      ...} (see extract_defaults)
  Returns:
  --------
    samples: pd daataframe the filtered sample list
  """
  extract.update(extract_defaults)
  rename.update(dup)
  prevlen = len(sampless)
  sampless[extract['ref_type']] = stype

  # checking no duplicate in the buckets
  for k, val in sampless.iterrows():
    withsamesize = sampless[sampless[extract["legacy_size"]] == val[extract["legacy_size"]]]
    if len(withsamesize) > 1:
      if len(withsamesize[withsamesize[extract["ref_name"]] == val[extract["ref_name"]]]) < 2:
        raise ValueError('we have duplicate samples with different names!')
      else:
        for l, v in withsamesize.iloc[1:].iterrows():
          sampless = sampless.drop(l)
  print("we had " + str(prevlen - len(sampless)) + " duplicates in the release buckets")
  # check: currently, below lines prevent forcekeep from working on true duplicates
  # (aka same size  file). Need to think about how to bring back the forcekeep functionality
  sampless[extract['ref_arxspan_id']] = [rename[name] if name in rename else name for name in sampless[extract['ref_arxspan_id']]]
  names = []
  # need to keep track of whether we're adding more than one new entry for a given sample id
  subrefsamples = refsamples[refsamples[extract['ref_type']] == stype]
  for k, val in sampless.iterrows():
    val = val[extract['ref_arxspan_id']]
    names.append(val)
    sampless.loc[k, extract['version']] = len(subrefsamples[subrefsamples[extract['ref_arxspan_id']] == val]) + names.count(val)
  sampless[extract['version']] = sampless[extract['version']].astype(int)

  sampless[extract['patient_id']] = [val[extract['patient_id']] if
                                     len(refsamples[refsamples[extract['ref_arxspan_id']] == val[extract['ref_arxspan_id']]]) < 1 else
                                     refsamples[refsamples[extract['ref_arxspan_id']] == val[extract['ref_arxspan_id']]][extract['patient_id']].values[0]
                                     for i, val in sampless.iterrows()]

  return sampless


def completeFromMasterSheet(samples, notfound, toupdate=TO_UPDATE, 
                            my_id=MY_ID,
                            pv_index=SAMPLEID,
                            master_index="arxspan_id",
                            pv_tokeep = ['Culture Type', 'Culture Medium'],
                            mystorage_id=MYSTORAGE_ID,
                            masterfilename="ACH",
                            nanslist=['None', 'nan', 'Unknown', None, np.nan],
                            depmap_pv=DEPMAP_PV,
                            depmap_taiga=DEPMAP_TAIGA):
  """complete the missing sample information from a given DepMap Ops MasterSheet

  Args:
      samples ([type]): [description]
      notfound ([type]): [description]
      toupdate ([type], optional): [description]. Defaults to TO_UPDATE.
      my_id ([type], optional): [description]. Defaults to MY_ID.
      pv_index ([type], optional): [description]. Defaults to SAMPLEID.
      master_index (str, optional): [description]. Defaults to "arxspan_id".
      pv_tokeep (list, optional): [description]. Defaults to ['Culture Type', 'Culture Medium'].
      mystorage_id ([type], optional): [description]. Defaults to MYSTORAGE_ID.
      masterfilename (str, optional): [description]. Defaults to "ACH".
      nanslist (list, optional): [description]. Defaults to ['None', 'nan', 'Unknown', None].
      depmap_pv ([type], optional): [description]. Defaults to DEPMAP_PV.
      depmap_taiga ([type], optional): [description]. Defaults to DEPMAP_TAIGA.

  Returns:
      [type]: [description]
  """
  di = {k: [] for k, _ in toupdate.items()}

  sheets = Sheets.from_files(my_id, mystorage_id)

  depmap_pv = sheets.get(depmap_pv).sheets[0].to_frame(header=2)
  depmap_pv = depmap_pv.drop(depmap_pv.iloc[:1].index)
  depmap_pv = depmap_pv.drop(depmap_pv.iloc[:1].index).set_index(
      pv_index, drop=True)[pv_tokeep]
  depmap_master = tc.get(name=depmap_taiga, file=masterfilename).set_index(
      master_index, drop=True)
  depmap_master = depmap_master.join(depmap_pv)
  unk = []
  for k, v in samples.loc[notfound].iterrows():
    a = depmap_master[depmap_master.index == v.arxspan_id]
    if len(a) == 0:
      for k, v in toupdate.items():
        di[k].append('')
      # for these samples I will need to check and manually add the data in the list
      print("no data found for "+str(v.arxspan_id))
      unk.append(k)
      continue
    for k, j in toupdate.items():
      l = []
      if len(j) > 0:
        ll = depmap_master.loc[v.arxspan_id, j] if len(
          a) == 1 else depmap_master.loc[v.arxspan_id, j].iloc[0]
        for val in ll.tolist():
          if val not in nanslist:
            l.append(val)
        di[k].append(", ".join(l))
  for k, v in di.items():
    samples.loc[notfound, k] = v
  return samples, unk

def loadWES(samplesetname, 
            workspaces=[
            "terra-broad-cancer-prod/CCLE_DepMap_WES",
            "terra-broad-cancer-prod/Getz_IBM_CellLines_Exomes"],
            sources=["ccle","ibm"],
            maxage=MAXAGE,
            baits = 'ice',
            stype = "wes",
            **kwargs):
  """
  function to load WES data from GP workspaces

  @see load()
  """
  return load(samplesetname=samplesetname, workspaces=workspaces,
              sources=sources, maxage=maxage, baits=baits, stype=stype, **kwargs)

def loadWGS(samplesetname, 
            workspaces=[
                wgsworkspace1,
                wgsworkspace2],
            sources=[wgssource1,wgssource2],
            maxage=MAXAGE,
            baits = 'genome',
            stype = "wgs",
            **kwargs):
  """
  function to load WGS data from GP workspaces

  @see load()
  """
  return load(samplesetname=samplesetname, workspaces=workspaces,
              sources=sources, maxage=maxage, baits=baits, stype=stype, **kwargs)

def loadRNA(samplesetname=SAMPLESETNAME,
            workspaces=[
                rnaworkspace6,
                rnaworkspace7],
            sources=[rnasource6, rnasource7],
            maxage=MAXAGE,
            baits='polyA',
            stype="rna", **kwargs):
  """
  function to load RNA data from GP workspaces 

  @see load()
  """
  return load(samplesetname=samplesetname, workspaces=workspaces,
              sources=sources, maxage=maxage, baits=baits, stype=stype, **kwargs)

def load(samplesetname, workspaces,
         sources,
         maxage,
         baits,
         stype,
        toupdate=TO_UPDATE,
        pv_index=SAMPLEID,
        master_index="arxspan_id",
        my_id=MY_ID,
        creds=SHEETCREDS,
        mystorage_id=MYSTORAGE_ID,
        refsheet_url = REFSHEET_URL,
        depmappvlink = DEPMAP_PV,
        extract_to_change = EXTRACT_TO_CHANGE,
        # version 102
        match = MATCH,
        pv_tokeep=['Culture Type', 'Culture Medium'],
        masterfilename="ACH",
        nanslist=['None', 'nan', 'Unknown', None, np.nan],
        depmap_taiga=DEPMAP_TAIGA,
        toraise=TORAISE,
        participantslicepos=10, accept_unknowntypes=True,
        recomputehash=True):
  """function to load and extract data from the GP workspaces

  Args:
      samplesetname ([type]): [description]
      workspaces ([type]): [description]
      sources ([type]): [description]
      maxage ([type]): [description]
      baits ([type]): [description]
      stype ([type]): [description]
      toupdate ([type], optional): [description]. Defaults to TO_UPDATE.
      pv_index ([type], optional): [description]. Defaults to SAMPLEID.
      master_index (str, optional): [description]. Defaults to "arxspan_id".
      my_id ([type], optional): [description]. Defaults to MY_ID.
      creds ([type], optional): [description]. Defaults to SHEETCREDS.
      mystorage_id ([type], optional): [description]. Defaults to MYSTORAGE_ID.
      refsheet_url ([type], optional): [description]. Defaults to REFSHEET_URL.
      depmappvlink ([type], optional): [description]. Defaults to DEPMAP_PV.
      extract_to_change ([type], optional): [description]. Defaults to EXTRACT_TO_CHANGE.
      match ([type], optional): [description]. Defaults to MATCH.
      pv_tokeep (list, optional): [description]. Defaults to ['Culture Type', 'Culture Medium'].
      masterfilename (str, optional): [description]. Defaults to "ACH".
      nanslist (list, optional): [description]. Defaults to ['None', 'nan', 'Unknown'].
      depmap_taiga ([type], optional): [description]. Defaults to DEPMAP_TAIGA.
      toraise ([type], optional): [description]. Defaults to TORAISE.
      participantslicepos (int, optional): [description]. Defaults to 10.
      accept_unknowntypes (bool, optional): [description]. Defaults to True.
      recomputehash (bool, optional): [description]. Defaults to True.

  Raises:
      ValueError: [description]

  Returns:
      [type]: [description]
  """
  release = samplesetname
  sheets = Sheets.from_files(my_id, mystorage_id)
  ccle_refsamples = sheets.get(refsheet_url).sheets[0].to_frame(index_col=0)

  ## Adding new data

  # we will be missing "primary disease","sm_id", "cellosaurus_id", "gender, "age", "primary_site", "primary_disease", "subtype", "subsubtype", "origin", "comments"
  #when SMid: match== 
  samples, _, noarxspan = GetNewCellLinesFromWorkspaces(stype=stype, 
                                                        maxage=maxage, refurl=refsheet_url, 
                                                        wmfroms=workspaces,
                                                        sources=sources, match=match, 
                                                        participantslicepos=participantslicepos, 
                                                        accept_unknowntypes=accept_unknowntypes, 
                                                        extract=extract_to_change, 
                                                        recomputehash=recomputehash)

  ### finding back arxspan
  noarxspan = tracker.retrieveFromCellLineName(noarxspan, ccle_refsamples, 
  datatype=stype, depmappvlink=depmappvlink, extract=extract_to_change)

  #assess any potential issues
  samples = pd.concat([samples, noarxspan[noarxspan.arxspan_id!='0']], sort=False)
  noarxspan = noarxspan[noarxspan.arxspan_id=='0']

  extract_defaults.update(extract_to_change)
  samples = assessAllSamples(
      samples, ccle_refsamples, stype=stype, rename={}, extract=extract_defaults)

  # converting date to str
  samples[extract_defaults['release_date']] = [h.inttodate(i) for i in
      samples[extract_defaults['release_date']]]

  if len(noarxspan) > 0:
    print("we found "+str(len(noarxspan))+" samples without arxspan_ids!!")
    noarxspan = noarxspan.sort_values(by = 'stripped_cell_line_name')
    dfToSheet(samples, "depmap ALL samples not found", creds)
    noarxspan.to_csv('temp/noarxspan_'+stype+'_' + release + '.csv')
    if h.askif("Please review the samples (on 'depmap samples not found') and write yes once \
      finished, else write no to quit and they will not be added"):
      updated_samples = sheets.get(
        "https://docs.google.com/spreadsheets/d/1yC3brpov3JELvzNoQe3eh0W196tfXzvpa0jUezMAxIg").sheets[
          0].to_frame().set_index('sample_id')
      samples = pd.concat([samples, updated_samples], sort=False)
  samples, notfound = tracker.updateFromTracker(samples, ccle_refsamples)
  
  for val in toraise:
    if val in samples['arxspan_id'].tolist():
      raise ValueError('some samples were amongst the known wrong samples')
  
  samples['baits'] = baits
  if len(samples.loc[notfound])>0:
    print("we found some samples where we could not get annotations. \
      trying to infer it from depmap master sheet and arxspan export")
    samples, unk = completeFromMasterSheet(samples, notfound, 
      toupdate=toupdate, 
      my_id=my_id,
      pv_index=pv_index,
      master_index=master_index,
      mystorage_id=mystorage_id,
      pv_tokeep=pv_tokeep,
      masterfilename=masterfilename,
      nanslist=nanslist,
      depmap_pv=depmappvlink,
      depmap_taiga=depmap_taiga)
    if len(unk) > 0:
      print("some samples could still not be inferred")
      dfToSheet(samples.loc[notfound], "depmap samples not found", creds)
      samples.loc[notfound].to_csv('temp/notfound_'+stype+'_'+release+'.csv')
      if h.askif("Please review the samples (on 'depmap samples not found') and write yes once \
        finished, else write no to quit and they will not be added"):
        updated_samples = sheets.get(
          "https://docs.google.com/spreadsheets/d/1yC3brpov3JELvzNoQe3eh0W196tfXzvpa0jUezMAxIg").sheets[
          0].to_frame().set_index('sample_id')
        samples.loc[updated_samples.index, updated_samples.columns] = updated_samples.values
  
  dfToSheet(samples, 'depmap ALL samples found', secret=creds)
  samples.to_csv('temp/new_'+stype+'_'+release+'.csv')
  return samples


def updateWES(samples, samplesetname, bucket="gs://cclebams/wes/",
                name_col="index", values=['legacy_bam_filepath', 'legacy_bai_filepath'], 
                filetypes=['bam', 'bai'],
                my_id=MY_ID,
                mystorage_id=MYSTORAGE_ID,
                refworkspace=WESMUTWORKSPACE,
                cnworkspace=WESCNWORKSPACE,
                stype= "wes",
                baits= 'ICE',
                extract = {},
                creds = SHEETCREDS,
                sampletrackername =SHEETNAME,
                refsheet_url = REFSHEET_URL,):
  """[summary]

  Args:
      samples ([type]): [description]
      samplesetname ([type]): [description]
      bucket (str, optional): [description]. Defaults to "gs://cclebams/wes/".
      name_col (str, optional): [description]. Defaults to "index".
      values (list, optional): [description]. Defaults to ['legacy_bam_filepath', 'legacy_bai_filepath'].
      filetypes (list, optional): [description]. Defaults to ['bam', 'bai'].
      my_id ([type], optional): [description]. Defaults to MY_ID.
      mystorage_id ([type], optional): [description]. Defaults to MYSTORAGE_ID.
      refworkspace ([type], optional): [description]. Defaults to WESMUTWORKSPACE.
      cnworkspace ([type], optional): [description]. Defaults to WESCNWORKSPACE.
      stype (str, optional): [description]. Defaults to "wes".
      baits (str, optional): [description]. Defaults to 'ICE'.
      extract (dict, optional): [description]. Defaults to {}.
      creds ([type], optional): [description]. Defaults to SHEETCREDS.
      sampletrackername ([type], optional): [description]. Defaults to SHEETNAME.
      refsheet_url ([type], optional): [description]. Defaults to REFSHEET_URL.
  """

  # uploading to our bucket (now a new function)
  terra.changeToBucket(samples, bucket, name_col=name_col,
                        values=values, filetypes=filetypes, catchdup=True, test=False)
  
  extract.update(extract_defaults)
  sheets = Sheets.from_files(my_id, mystorage_id)
  ccle_refsamples = sheets.get(refsheet_url).sheets[0].to_frame(index_col=0)

  names=[]
  subccle_refsamples = ccle_refsamples[ccle_refsamples['datatype'] == stype]
  for k, val in samples.iterrows():
    val = val["arxspan_id"]
    names.append(val)
    samples.loc[k, 'version'] = len(subccle_refsamples[subccle_refsamples['arxspan_id'] == val]) + names.count(val)
  samples['version'] = samples['version'].astype(int)

  ccle_refsamples = ccle_refsamples.append(samples, sort=False)

  dfToSheet(ccle_refsamples,sampletrackername, secret=creds)

  pairs = myterra.setupPairsFromSamples(samples, subccle_refsamples, extract)

  #uploading new samples to mut
  refwm = dm.WorkspaceManager(refworkspace)
  refwm = refwm.disable_hound()
  refwm.upload_samples(samples)
  refwm.upload_entities('pairs', pairs)
  refwm.update_pair_set(pair_set_id=samplesetname, pair_ids=pairs.index)
  sam = refwm.get_samples()

  pair = refwm.get_pairs()
  refwm.update_pair_set(pair_set_id='all', pair_ids=pair.index)

  refwm.update_pair_set(pair_set_id='all_'+baits, pair_ids=pair[pair["case_sample"].isin(
      [i for i in sam[(sam['baits'] == baits) | (sam['baits'].isna())].index.tolist() if i != 'nan'])].index)

  #creating a sample set
  refwm.update_sample_set(sample_set_id=samplesetname, sample_ids=samples.index)
  refwm.update_sample_set(sample_set_id='all', sample_ids=[
                          i for i in sam.index.tolist() if i != 'nan'])

  refwm.update_sample_set(sample_set_id='all_'+baits, sample_ids=[i for i in sam[(
      sam['baits'] == baits) | (sam['baits'].isna())].index.tolist() if i != 'nan'])

  #and CN
  cnwm = dm.WorkspaceManager(cnworkspace)
  cnwm = cnwm.disable_hound()
  cnwm.upload_samples(samples)
  cnwm.upload_entities('pairs', pairs)
  cnwm.update_pair_set(pair_set_id=samplesetname, pair_ids=pairs.index)
  sam = cnwm.get_samples()

  pair = cnwm.get_pairs()
  cnwm.update_pair_set(pair_set_id='all', pair_ids=pair.index)
  cnwm.update_pair_set(pair_set_id='all_'+baits, pair_ids=pair[pair["case_sample"].isin(
      [i for i in sam[(sam['baits'] == baits) | (sam['baits'].isna())].index.tolist() if i != 'nan'])].index)
  #creating a sample set
  cnwm.update_sample_set(sample_set_id=samplesetname, sample_ids=samples.index)
  cnwm.update_sample_set(sample_set_id='all', sample_ids=[
                        i for i in sam.index.tolist() if i != 'nan'])
  cnwm.update_sample_set(sample_set_id='all_'+baits, sample_ids=[i for i in sam[(
      sam['baits'] == baits) | (sam['baits'].isna())].index.tolist() if i != 'nan'])


def update(samples, stype, bucket, refworkspace, samplesetname=SAMPLESETNAME,
          name_col="index", values=['legacy_bam_filepath', 'legacy_bai_filepath'],
          filetypes=['bam', 'bai'],
          my_id=MY_ID,
          mystorage_id=MYSTORAGE_ID,
          creds=SHEETCREDS,
          sampletrackername=SHEETNAME,
          refsheet_url=REFSHEET_URL,):
  """update the samples on a depmapomics terra processing workspace

  Args:
      samples ([type]): [description]
      samplesetname ([type]): [description]
      stype ([type]): [description]
      bucket ([type]): [description]
      refworkspace ([type]): [description]
      name_col (str, optional): [description]. Defaults to "index".
      values (list, optional): [description]. Defaults to ['legacy_bam_filepath', 'legacy_bai_filepath'].
      filetypes (list, optional): [description]. Defaults to ['bam', 'bai'].
      my_id ([type], optional): [description]. Defaults to MY_ID.
      mystorage_id ([type], optional): [description]. Defaults to MYSTORAGE_ID.
      creds ([type], optional): [description]. Defaults to SHEETCREDS.
      sampletrackername ([type], optional): [description]. Defaults to SHEETNAME.
      refsheet_url ([type], optional): [description]. Defaults to REFSHEET_URL.
  """

  # uploading to our bucket (now a new function)
  terra.changeToBucket(samples, bucket, name_col=name_col,
                       values=values, filetypes=filetypes, catchdup=True, dryrun=False)

  samplesetname
  sheets = Sheets.from_files(my_id, mystorage_id)
  ccle_refsamples = sheets.get(refsheet_url).sheets[0].to_frame(index_col=0)

  names = []
  subccle_refsamples = ccle_refsamples[ccle_refsamples['datatype'] == stype]
  for k, val in samples.iterrows():
    val = val["arxspan_id"]
    names.append(val)
    samples.loc[k, 'version'] = len(
        subccle_refsamples[subccle_refsamples['arxspan_id'] == val]) + names.count(val)
  samples['version'] = samples['version'].astype(int)

  ccle_refsamples = ccle_refsamples.append(samples, sort=False)
  dfToSheet(ccle_refsamples, sampletrackername, secret=creds)

  #uploading new samples to mut
  refwm = dm.WorkspaceManager(refworkspace).disable_hound()
  refwm.upload_samples(samples)
  sam = refwm.get_samples()

  #creating a sample set
  refwm.update_sample_set(sample_set_id=samplesetname,
                          sample_ids=samples.index)
  refwm.update_sample_set(sample_set_id='all', sample_ids=[
                          i for i in sam.index.tolist() if i != 'nan'])
