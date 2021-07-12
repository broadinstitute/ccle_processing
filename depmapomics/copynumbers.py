#cn.py

import numpy as np
from gsheets import Sheets
from depmapomics import tracker as track
from depmapomics import utils
from depmapomics.qc import cn
from depmapomics.config import *
from IPython.display import Image, display
import dalmatian as dm
import pandas as pd
from taigapy import TaigaClient
tc = TaigaClient()
import os
from genepy import mutations as mut
from genepy.utils import helper as h
from genepy import terra
from genepy.google import gcp
from genepy import rna
import matplotlib.pyplot as plt

def renameColumns(df):
  """
  rename some of the main columns names from RSEM, GATK.. to more readable column names

  Args:
  -----
    df: the df to rename

  Returns:
  ------
    df the renamed df
  """
  return df.rename(columns=COLRENAMING)

def loadFromGATKAggregation(refworkspace,  sortby=[SAMPLEID, 'Chromosome', "Start", "End"], 
                            save_output='', doCleanup=True,
                            todrop=[], showPlots=False, colname="combined_seg_file",
                            plotColname="modeled_segments_plot_tumor", tempFolder="temp/",
                            toremove=["readgroup_ubams", ],
                            sampleset="all", colRenaming=COLRENAMING):
  """
   """
  wm = dm.WorkspaceManager(refworkspace)
  if save_output:
    terra.saveConfigs(refworkspace, os.path.join(save_output, 'terra/'))
  
  if doCleanup:
    print('cleaning up')
    for val in toremove:
        wm.disable_hound().delete_entity_attributes('samples', toremove)
    a = wm.get_samples()
    e = []
    for i in a[toremove].values.tolist():
        if i is not np.nan:
            e.extend(i)
    gcp.rmFiles(e)
  
  segments = pd.read_csv(wm.get_entities(
        'sample_set').loc[sampleset, colname], sep='\t').rename(columns=colRenaming)
 
  segments = segments[~segments[SAMPLEID].isin(todrop)].reset_index(drop=True)
  if "chr" in segments['Chromosome'][0]:
     segments['Chromosome'] = [i[3:] for i in segments['Chromosome']]
  segments.Segment_Mean = 2**segments.Segment_Mean
  segments.Start = segments.Start.astype(int)
  segments.End = segments.End.astype(int)
  segments.loc[segments[segments.Chromosome.isin(
    ['X', 'Y'])].index, 'Segment_Mean'] = segments[segments.Chromosome.isin(['X', 'Y'])]['Segment_Mean']/2
  segments = segments.sort_values(by=sortby)
  
  print("loading "+ str(len(set(segments[SAMPLEID])))+ " rows")
  if showPlots:
    # plotting results of CN calls for this new sample set
    for i, (k, val) in enumerate(wm.get_samples().loc[wm.get_sample_sets().loc[
        sampleset].samples].iterrows()):
      plot = val[plotColname]
      os.system('gsutil cp '+plot+' '+tempFolder)
      print(k)
      print(val['arxspan_id'], val['sex'])
      if i > 30:
        continue
      display(Image(os.path.join(tempFolder,plot.split('/')[-1])))
  return segments

  
def updateTracker(tracker, selected, samplesetname, samplesinset, lowqual, newgs='',
                  sheetcreds = SHEETCREDS,
                  sheetname=SHEETNAME, procqc=[], bamqc=[], refworkspace=None,
                  onlycol=['internal_bam_filepath', 'internal_bai_filepath'],
                  ):
  
  # updating locations of bam files and extracting infos
  if newgs and refworkspace is not None:

    res, _=terra.changeGSlocation(refworkspace, newgs=newgs, onlycol=onlycol,
                                  entity='sample', keeppath=False, dry_run=False, 
                                  onlysamples=samplesinset)
    tracker.loc[res.index.tolist()][['legacy_size', 'legacy_crc32c_hash']
                                      ] = tracker.loc[
                                        res.index.tolist()][
                                          ['size', 'crc32c_hash']].values
    tracker.loc[res.index.tolist()][HG38BAMCOL]=res[onlycol[:2]].values
    tracker.loc[res.index.tolist(), 'size']=[gcp.extractSize(
      i)[1] for i in gcp.lsFiles(res[onlycol[0]].tolist(), '-l')]
    tracker.loc[res.index.tolist(), 'crc32c_hash']=[gcp.extractHash(
      i) for i in gcp.lsFiles(res[onlycol[0]].tolist(), '-L')]
    tracker.loc[res.index.tolist(), 'md5_hash']=gcp.catFiles(
      dm.WorkspaceManager(refworkspace).get_samples().loc[
        samplesinset, 'analysis_ready_bam_md5'].tolist(), cut=32)

  # computing QC
  print('looking for QC..')
  dataProc={}
  dataBam={}
  if procqc and refworkspace is not None:
    dataProc=getQC(workspace=refworkspace, only=samplesinset, qcname=procqc)
  if bamqc and refworkspace is not None:
    dataBam=getQC(workspace=refworkspace, only=samplesinset, qcname=bamqc)
  for k,v in dataProc.items():
    if k =='nan':
      continue
    a = tracker.loc[k,'processing_qc']
    a = '' if a is np.nan else a
    tracker.loc[k,'processing_qc'] = str(v) + ',' + a
  for k,v in dataBam.items():
    if k =='nan':
      continue
    a = tracker.loc[k,'bam_qc']
    a = '' if a is np.nan else a
    tracker.loc[k,'bam_qc'] = str(v) + ',' + a
  
  tracker.loc[tracker[tracker.datatype.isin(['wes',"wgs"])].index, samplesetname]=0
  len(selected)
  tracker.loc[selected, samplesetname]=1
  tracker.loc[samplesinset, ['low_quality', 'blacklist', 'prioritized']]=0
  tracker.loc[lowqual,'low_quality']=1
  tracker.loc[lowqual,'blacklist']=1
  dfToSheet(tracker, sheetname, secret=sheetcreds)
  print("updated the sheet, please reactivate protections")


def managingDuplicates(samples, failed, datatype, tracker):
  """removes duplicates and solves failed data

  Args:
      failed (list): [description]
      datatype (str): [description]
      tracker (dataframe[datatype, prioritized, arxspan_id, index]): [description]
      samples (list): [description]

  Returns:
      [type]: [description]
      [type]: [description]
  """
  # selecting the right arxspan id (latest version)
  renaming = track.removeOlderVersions(names=samples, 
    refsamples=tracker[tracker.datatype == datatype], priority="prioritized")

  # reparing QC when we have a better duplicate
  ref=pd.DataFrame(
    tracker[tracker.datatype == datatype]['arxspan_id'])
  replace=0
  for val in failed:
    if val in list(renaming.keys()):
      a=ref[ref.arxspan_id == ref.loc[val].arxspan_id].index
      for v in a:
        if v not in failed:
          renaming[v]=renaming.pop(val)
          replace += 1
          break
  print('could replace:')
  print(replace)
  return renaming


def postProcess(refworkspace, sampleset='all', save_output="", doCleanup=True,  sortby=[
        SAMPLEID, 'Chromosome', "Start", "End"], todrop=[], priority=[],
        genechangethresh=0.025, segmentsthresh=2000, ensemblserver=ENSEMBL_SERVER_V,
        source_rename={}, useCache=False):
  """post process an aggregated CN segment file, the CCLE way

  (usually a CN segment file from the Aggregate_CN terra workflow)

  Args:
      refworkspace ([type]): [description]
      sampleset (str, optional): [description]. Defaults to 'all'.
      save_output (str, optional): [description]. Defaults to "".
      doCleanup (bool, optional): [description]. Defaults to True.
      sortby (list, optional): [description]. Defaults to [ SAMPLEID, 'Chromosome', "Start", "End"].
      todrop (list, optional): [description]. Defaults to [].
      priority (list, optional): [description]. Defaults to [].
      genechangethresh (float, optional): [description]. Defaults to 0.025.
      segmentsthresh (int, optional): [description]. Defaults to 2000.
      ensemblserver ([type], optional): [description]. Defaults to ENSEMBL_SERVER_V.
      source_rename (dict, optional): [description]. Defaults to {}.
      useCache (bool, optional): [description]. Defaults to False.

  Returns:
      [type]: [description]
  """
  h.createFoldersFor(save_output)
  print('loading CN from Terra')
  segments = loadFromGATKAggregation(
      refworkspace, sampleset=sampleset, sortby=sortby, todrop=todrop, doCleanup=doCleanup)
  print('making gene level copy number')

  mybiomart = utils.generateGeneNames(
      ensemble_server=ensemblserver, useCache=useCache,
      attributes=['start_position', 'end_position', "chromosome_name"])
  mybiomart = mybiomart.rename(columns={'start_position': 'start',
                                        'end_position': 'end',
                                        'chromosome_name': 'Chromosome',
                                        })
  mybiomart['Chromosome'] = mybiomart['Chromosome'].astype(str)
  mybiomart = mybiomart.sort_values(by=['Chromosome', 'start', 'end'])
  mybiomart = mybiomart[mybiomart['Chromosome'].isin(
      set(segments['Chromosome']))]
  mybiomart['gene_name'] = [i['hgnc_symbol'] +
                            ' (' + str(i['entrezgene_id']).split('.')[0] +
                            ')' for _, i in mybiomart.iterrows()]
  mybiomart = mybiomart.drop_duplicates('gene_name', keep="first")
  genecn = mut.toGeneMatrix(mut.manageGapsInSegments(segments), mybiomart)

  # validation step
  print('summary of the gene cn data:')
  print(genecn.values.min(), genecn.values.mean(), genecn.values.max())
  mut.checkGeneChangeAccrossAll(genecn, thresh=genechangethresh)
  failed = mut.checkAmountOfSegments(segments, thresh=segmentsthresh)

  print("failed our QC")
  print(failed)
  if source_rename:
    segments = segments.replace({'Source': source_rename})
  if save_output:
    h.listToFile(failed, save_output+'failed.txt')
  # subsetting
  segments = segments[~segments[SAMPLEID].isin(
      (set(failed) | set(todrop))-set(priority))].reset_index(drop=True)
  genecn = genecn[~genecn.index.isin(
      (set(failed) | set(todrop))-set(priority))]

  #saving
  print('saving files')
  segments.to_csv(save_output+ 'segments_all.csv', index=False)
  genecn.to_csv(save_output+ 'genecn_all.csv')
  print("done")

  return segments, genecn, failed


def CCLEPostProcessing(wesrefworkspace=WESCNWORKSPACE, wgsrefworkspace=WGSWORKSPACE,
                       samplesetname=SAMPLESETNAME, AllSamplesetName='all',
                       my_id=MY_ID, mystorage_id=MYSTORAGE_ID,
                       sheetcreds=SHEETCREDS, sheetname=SHEETNAME,
                       refsheet_url=REFSHEET_URL, todrop=KNOWN_DROP,
                       prevgenecn='ccle',
                       taiga_dataset=TAIGA_CN, dataset_description=CNreadme,
                       subsetsegs=[SAMPLEID, 'Chromosome',
                                   'Start', 'End', 'Segment_Mean',
                                   'Num_Probes', 'Status', 'Source'],
                       bamqc=BAMQC,
                       procqc=PROCQC,
                       source_rename=SOURCE_RENAME,
                      **kwargs):
  """the full CCLE Copy Number post processing pipeline (used only by CCLE)

  see postprocessing() to reproduce our analysis

  Args:
      wesrefworkspace ([type]): [description]
      wgsrefworkspace ([type]): [description]
      samplesetname ([type]): [description]
      AllSamplesetName (str, optional): [description]. Defaults to ''.
      my_id ([type], optional): [description]. Defaults to MY_ID.
      mystorage_id ([type], optional): [description]. Defaults to MYSTORAGE_ID.
      sheetcreds ([type], optional): [description]. Defaults to SHEETCREDS.
      sheetname ([type], optional): [description]. Defaults to SHEETNAME.
      refsheet_url ([type], optional): [description]. Defaults to REFSHEET_URL.
      prevgenecn ([type], optional): [description]. Defaults to tc.get(name=TAIGA_ETERNAL, file='CCLE_gene_cn').
      taiga_dataset (str, optional): [description]. Defaults to TAIGA_CN.
      dataset_description ([type], optional): [description]. Defaults to CNreadme.
      subsetsegs (list, optional): [description]. Defaults to [SAMPLEID, 'Chromosome', 'Start', 'End', 'Segment_Mean', 'Num_Probes', 'Status', 'Source'].
      bamqc ([type], optional): [description]. Defaults to BAMQC.
      procqc ([type], optional): [description]. Defaults to PROCQC.
      source_rename ([type], optional): [description]. Defaults to SOURCE_RENAME.
  """
  if prevgenecn is 'ccle':
    prevgenecn = tc.get(name=TAIGA_ETERNAL, file='CCLE_gene_cn')

  sheets = Sheets.from_files(my_id, mystorage_id)
  tracker = sheets.get(refsheet_url).sheets[0].to_frame(index_col=0)

  wesrefwm = dm.WorkspaceManager(wesrefworkspace)
  wgsrefwm = dm.WorkspaceManager(wgsrefworkspace)  

  # doing wes
  print('doing wes')
  folder=os.path.join("temp", samplesetname, "wes_")
  priority=tracker[(tracker.datatype=='wes')&(tracker.prioritized == 1)].index.tolist()
  todropwes=todrop+tracker[(tracker.datatype=='wes')&(tracker.blacklist == 1)].index.tolist()
  wessegments, genecn, wesfailed = postProcess(wesrefworkspace, AllSamplesetName if AllSamplesetName else samplesetname,
              todrop=todropwes,
              save_output=folder,
              priority=priority, **kwargs)
  
  wesrenaming = managingDuplicates(set(wessegments[SAMPLEID]), (set(
      wesfailed) - set(priority)) | set(todropwes), "wes", tracker)
  h.dictToFile(wesrenaming, folder+"sample_renaming.json")
  print('renaming')
  wespriosegments = wessegments[wessegments[SAMPLEID].isin(set(wesrenaming.keys()))].replace(
    {SAMPLEID: wesrenaming}).reset_index(drop = True)
  wespriogenecn = genecn[genecn.index.isin(set(wesrenaming.keys()))].rename(index=wesrenaming)

  # annotating source
  for v in set(wespriosegments.DepMap_ID):
    wespriosegments.loc[wespriosegments[wespriosegments.DepMap_ID == v].index,
                 'Source'] = tracker[tracker.index == v].source.values[0]
    wespriosegments.Source = wespriosegments.Source.replace(source_rename)
    wespriosegments.Source += ' WES'
  
  #saving prio
  wespriosegments.to_csv(folder+"segments_all_latest.csv", index=False)
  wespriogenecn.to_csv(folder+"genecn_all_latest.csv")
  
  # doing wgs
  print('doing wgs')
  folder=os.path.join("temp", samplesetname, "wgs_")
  priority=tracker[(tracker.datatype=='wgs')&(tracker.prioritized == 1)].index.tolist()
  todropwgs=todrop+tracker[(tracker.datatype=='wgs')&(tracker.blacklist == 1)].index.tolist()
  wgssegments, genecn, wgsfailed = postProcess(wgsrefworkspace, AllSamplesetName if AllSamplesetName else samplesetname,
              todrop=todropwgs,
              save_output=folder,
              priority=priority, **kwargs)

  wgsrenaming = managingDuplicates(set(wgssegments[SAMPLEID]), (set(
      wgsfailed) - set(priority)) | set(todropwgs), "wgs", tracker)  

  h.dictToFile(wgsrenaming, folder+"sample_renaming.json")

  print('renaming')
  wgspriosegments = wessegments[wessegments[SAMPLEID].isin(set(wgsrenaming.keys()))].replace(
    {SAMPLEID: wgsrenaming}).reset_index(drop = True)
  wgspriogenecn = genecn[genecn.index.isin(set(wgsrenaming.keys()))].rename(index=wgsrenaming)
  
  # annotating source
  for v in set(wgspriosegments.DepMap_ID):
    wgspriosegments.loc[wgspriosegments[wgspriosegments.DepMap_ID == v].index,
                 'Source'] = tracker[tracker.index == v].source.values[0]
    wgspriosegments.Source = wgspriosegments.Source.replace(source_rename)
    wgspriosegments.Source += ' WGS'

  #saving prio
  wgspriosegments.to_csv(folder+ "segments_all_latest.csv", index=False)
  wgspriogenecn.to_csv(folder+ "genecn_all_latest.csv")
   
  print('comparing to previous version')
  #h.compareDfs(priosegments, tc.get(name='depmap-a0ab', file='CCLE_segment_cn'))
  h.compareDfs(wespriogenecn, prevgenecn)
  
  #adding to the sample tracker the sequencing that were selected and the ones that failed QC
  selected = {j:i for i,j in wesrenaming.items()}
  selected.update({j:i for i,j in wgsrenaming.items()})
  try:
    wgssamplesinset=[i['entityName'] for i in wesrefwm.get_entities(
            'sample_set').loc[samplesetname].samples]
    updateTracker(tracker, selected, wgssamplesinset, samplesetname,
                  list(wgsfailed), sheetcreds=sheetcreds, sheetname=sheetname, bamqc=bamqc, procqc=procqc)
  except:
    print('no wgs for this sampleset')
  
  try:
    wessamplesinset=[i['entityName'] for i in wgsrefwm.get_entities(
        'sample_set').loc[samplesetname].samples]
    updateTracker(tracker, selected, wessamplesinset, samplesetname,
                  list(wesfailed), sheetcreds=sheetcreds, sheetname=sheetname, bamqc=bamqc, procqc=procqc)
  except:
    print('no wes for this sampleset')
  
  #merging WES/WGS
  folder=os.path.join("data", samplesetname, '')
  mergedsegments =  wgspriosegments.append(wespriosegments[~wespriosegments[SAMPLEID].isin(
    set(wgspriosegments[SAMPLEID]))])[subsetsegs]
  mergedgenecn =  wgspriogenecn.append(wespriogenecn[~wespriogenecn.index.isin(
    set(wgspriogenecn.index))])

  mergedgenecn.to_csv(folder+ "merged_genecn_all.csv")
  mergedsegments.to_csv(folder+ "merged_segments_all_.csv",index=False)

  #uploading to taiga
  print('uploading to taiga')
  tc.update_dataset(changes_description="new "+samplesetname+" release! (removed misslabellings, see changelog)",
                    dataset_permaname=taiga_dataset,
                    upload_files=[
                      {
                        "path": "temp/"+samplesetname+"/wes_segments_all_latest.csv",
                        "format": "TableCSV",
                        "encoding": "utf-8"
                      },
                      {
                        "path": "temp/"+samplesetname+"/wes_genecn_all_latest_.csv",
                        "format": "NumericMatrixCSV",
                        "encoding": "utf-8"
                      },
                      {
                        "path": "temp/"+samplesetname+"/wes_segments_all.csv",
                        "format": "TableCSV",
                        "encoding": "utf-8"
                      },
                      {
                        "path": "temp/"+samplesetname+"/wes_genecn_all.csv",
                        "format": "NumericMatrixCSV",
                        "encoding": "utf-8"
                      },
                      {
                        "path": "temp/"+samplesetname+"/merged_genecn_all.csv",
                        "format": "NumericMatrixCSV",
                        "encoding": "utf-8"
                      },
                      {
                        "path": "temp/"+samplesetname+"/merged_segments_all.csv",
                        "format": "NumericMatrixCSV",
                        "encoding": "utf-8"
                      },
                      {
                        "path": "temp/"+samplesetname+"/wgs_segments_all.csv",
                        "format": "TableCSV",
                        "encoding": "utf-8"
                      },
                      {
                        "path": "temp/"+samplesetname+"/wgs_genecn_all.csv",
                        "format": "NumericMatrixCSV",
                        "encoding": "utf-8"
                      },
                      {
                        "path": "temp/"+samplesetname+"/wgs_segments_all_latest.csv",
                        "format": "TableCSV",
                        "encoding": "utf-8"
                      },
                      {
                        "path": "temp/"+samplesetname+"/wgs_genecn_all_latest.csv",
                        "format": "NumericMatrixCSV",
                        "encoding": "utf-8"
                      },
                    ],
                    dataset_description=dataset_description)
  print("done")
  return wespriosegments, wgspriosegments


def ProcessForAchilles(wespriosegs, wgspriosegs, samplesetname=SAMPLESETNAME, bad=["ACH-001011",
                        "ACH-001108",
                        "ACH-001187",
                        "ACH-002291"  # added for some reason?
                        # much more than that..
                        "ACH-002010",
                        "ACH-000314"], taiga_legacy_loc='depmap-wes-cn-data--08f3',
                       taiga_legacy_filename='legacy_segments',
                       taiga_dataset="cn-wes-achilles-4dcd",
                       dataset_description=Achillesreadme,
                       cytobandloc='data/hg38_cytoband.gz', 
                       gene_mapping=pd.read_csv('data/genemapping_19Q1.csv'),
                       prevsegments=tc.get(name=TAIGA_ETERNAL, file='CCLE_segment_cn'),
                       prevgenecn=(
                           2**tc.get(name=TAIGA_ETERNAL, file='CCLE_gene_cn'))-1,
                       gene_expected_count=tc.get(name=TAIGA_ETERNAL, 
                        file='CCLE_expression_proteincoding_genes_expected_count')):
  # load legacy_segments
  legacy_segments=tc.get(
    name=taiga_legacy_loc, file=taiga_legacy_filename).drop(columns='Unnamed: 0')
  legacy_segments['Status']='U'
  legacy_segments.loc[legacy_segments[legacy_segments.Chromosome.str.contains("chr")].index, "Chromosome"]=[
  i[3:] for i in legacy_segments[legacy_segments.Chromosome.str.contains("chr")].Chromosome]

  onlyinleg = set(legacy_segments[SAMPLEID]) - \
      (set(wespriosegs[SAMPLEID]) | (set(wgspriosegs[SAMPLEID])))
  #samegenes = set(prevgenecn.columns) & set(priogenecn.columns)

  onlyinleg=onlyinleg - set(bad)
  print('found samples that are only in the legacy datasets')
  print(onlyinleg)
  # merging 
  print('merging wes/wgs/legacy')
  mergedsegments=wespriosegs[~wespriosegs[SAMPLEID].isin(list(onlyinleg))].append(
    legacy_segments[legacy_segments[SAMPLEID].isin(list(onlyinleg))]).reset_index(drop=True)
  
  mergedsegments=wgspriosegs.append(
    mergedsegments[~mergedsegments[SAMPLEID].isin(set(wgspriosegs[SAMPLEID]))])

  mergedsegments=mergedsegments[[SAMPLEID, 'Chromosome', 'Start', 'End', 
  'Segment_Mean', 'Num_Probes', 'Status', 'Source']].sort_values(by=
    [SAMPLEID, 'Chromosome', 'Start', 'End']).reset_index(drop=True)
  
  #setting amplification status to U for X chromosome as it is artificially 
  #amplified in female samples:
  mergedsegments.loc[mergedsegments[mergedsegments.Chromosome ==
                                    "X"].index, 'Status']='U'
  #making the gene cn matrix
  print( 'making the gene cn')
  cyto=pd.read_csv(cytobandloc, sep='\t',
                    names=['chrom', 'start', 'end', 'loc', 'stains']).iloc[:-1]
  cyto['chrom']=[i[3:] for i in cyto['chrom']]
  gene_mapping['Chromosome'] = gene_mapping['Chromosome'].astype(str)
  gene_mapping = gene_mapping.sort_values(by=['Chromosome', 'start', 'end'])
  gene_mapping = gene_mapping[gene_mapping['Chromosome'].isin(set(mergedsegments['Chromosome']))]
  gene_mapping['gene_name'] = [i['symbol'] +
                            ' (' + str(i['ensembl_id']).split('.')[0] +
                            ')' for _, i in gene_mapping.iterrows()]
  mergedsegments=mut.manageGapsInSegments(mergedsegments, cyto=cyto)
  mergedgenecn=mut.toGeneMatrix(
      mergedsegments, gene_mapping, ).apply(lambda x: np.log2(1+x))

  # some QC
  print('copy number change with previous release')
  cn.plotCNchanges(mergedgenecn, prevgenecn.apply(
      lambda x: np.log2(1+x)), mergedsegments, prevsegments)
  
  if(mergedgenecn.values.max() > 100):
    print("\n\n\nTOO HIGH, not LOG2 transformed!")
  if(len(mergedgenecn.index) > len(set(mergedgenecn.index))):
    print("Duplicate CL, not reprioritized well!")
  

  # computing relationship with RNAseq
  print('correlation with RNAseq:')
  _, ax=plt.subplots()
  rna.rnaseqcorrelation(mergedgenecn.fillna(
      0), gene_expected_count.fillna(0), ax, name="current")
  rna.rnaseqcorrelation(prevgenecn[prevgenecn.index.isin(
      mergedgenecn.index.tolist())], gene_expected_count.fillna(0), ax, name="prev")

  h.compareDfs(mergedgenecn, prevgenecn)
  #h.compareDfs(mergedsegments, tc.get(name='depmap-a0ab', file='CCLE_segment_cn'))

  # saving 
  print('saving')
  mergedgenecn.to_csv('temp/'+samplesetname+'/achilles_gene_cn.csv')
  mergedsegments.to_csv('temp/'+samplesetname+'/achilles_segment.csv', index=False)

  #saving to taiga
  print('uploading to taiga')

  tc.update_dataset(changes_description="updated to new " + samplesetname + 
    " release! (updated from relabelling see google drive file for more info)",
                    dataset_permaname=taiga_dataset,
                    upload_files=[
                      {
                          "path": 'temp/'+samplesetname+'/achilles_segment.csv',
                          "format": "TableCSV",
                          "encoding": "utf-8"
                      },
                      {
                          "path": 'temp/'+samplesetname+'/achilles_gene_cn.csv',
                          "format": "NumericMatrixCSV",
                          "encoding": "utf-8"
                      },
                    ],
                    dataset_description=dataset_description)
  print("done")
