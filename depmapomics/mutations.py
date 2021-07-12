from genepy import terra
from genepy.utils import helper as h
from genepy import mutations as mut
import os
import dalmatian as dm
import pandas as pd
from genepy.google.gcp import cpFiles
import numpy as np
from collections import Counter

from depmapomics.config import *

import dalmatian as dm
import pandas as pd
from taigapy import TaigaClient
tc = TaigaClient()
from gsheets import Sheets


def download_maf_from_workspace(refwm, sample_set_ids=['all_ice', 'all_agilent'],
                                output_maf='/tmp/mutation_filtered_terra_merged.txt'):
  """
  TODO: javad to document
  """
  sample_sets = refwm.get_sample_sets()
  dfs = []
  for sample_set_id in sample_sets.index.intersection(sample_set_ids):
    cpFiles(sample_sets.loc[sample_set_id, 'filtered_CGA_MAF_aggregated'],
                '/tmp/tmp.txt', payer_project_id='broad-firecloud-ccle', verbose=False)
    df = pd.read_csv('/tmp/tmp.txt', sep='\t', low_memory=False)
    dfs.append(df)
  dfs_concat = pd.concat(dfs)
  dfs_concat.to_csv(output_maf, index=False, sep='\t')
  return dfs_concat


def annotateLikelyImmortalized(maf, sample_col=SAMPLEID,
                               genome_change_col="Genome_Change", TCGAlocs=['TCGAhsCnt', 'COSMIChsCnt'],
                               max_recurrence=0.05, min_tcga_true_cancer=5):
  maf['is_likely_immortalization'] = False
  leng = len(set(maf[sample_col]))
  tocheck = []
  for k, v in Counter(maf[genome_change_col].tolist()).items():
    if v > max_recurrence*leng:
      tocheck.append(k)
  for val in list(set(tocheck)-set([np.nan])):
    if np.nan_to_num(maf[maf[genome_change_col] == val][TCGAlocs], 0).max() < min_tcga_true_cancer:
      maf.loc[maf[maf[genome_change_col]
                  == val].index, 'is_likely_immortalization'] = True
  return maf


def addAnnotation(maf, NCBI_Build='37', Strand="+"):
  """
  adds NCBI_Build and Strand annotation on the whole maf file
  """
  maf['NCBI_Build'] = NCBI_Build
  maf['Strand'] = Strand
  return maf


def add_variant_annotation_column(maf):
    rename = {}
    for k, v in MUTATION_GROUPS.items():
      for e in v:
        rename[e] = k
    maf['Variant_annotation'] = [rename[i] for i in maf['Variant_Classification'].tolist()]
    return maf


def managingDuplicates(samples, failed, datatype, tracker):
  # selecting the right arxspan id (latest version)
  renaming = tracker.removeOlderVersions(names=samples,
                                         refsamples=tracker[tracker.datatype == datatype], 
                                         priority="prioritized")

  # reparing QC when we have a better duplicate
  ref = pd.DataFrame(
      tracker[tracker.datatype == datatype]['arxspan_id'])
  replace = 0
  for val in failed:
    if val in list(renaming.keys()):
      a = ref[ref.arxspan_id == ref.loc[val, 'arxspan_id']].index
      for v in a:
        if v not in failed:
          renaming[v] = renaming.pop(val)
          replace += 1
          break
  print('could replace:')
  print(replace)
  return renaming


def postProcess(refworkspace, sampleset='all', mutCol="mut_AC", save_output="", doCleanup=False,
        rename_cols={"i_ExAC_AF": "ExAC_AF", "Tumor_Sample_Barcode": SAMPLEID, 
        "Tumor_Seq_Allele2": "Tumor_Allele"},):
  """post process an aggregated MAF file the CCLE way

  (usually a MAF file from the Aggregate_MAF Terra worklflow)

  Args:
      refworkspace ([type]): [description]
      sampleset (str, optional): [description]. Defaults to 'all'.
      mutCol (str, optional): [description]. Defaults to "mut_AC".
      save_output (str, optional): [description]. Defaults to "".
      doCleanup (bool, optional): [description]. Defaults to False.
      rename_cols (dict, optional): [description]. Defaults to {"i_ExAC_AF": "ExAC_AF", "Tumor_Sample_Barcode": SAMPLEID, "Tumor_Seq_Allele2": "Tumor_Allele"}.

  Returns:
      [type]: [description]
  """
  h.createFoldersFor(save_output)
  print('loading from Terra')
  if save_output:
    terra.saveConfigs(refworkspace, save_output + 'config/')
  refwm = dm.WorkspaceManager(refworkspace)
  mutations = pd.read_csv(refwm.get_sample_sets().loc[sampleset, 'filtered_CGA_MAF_aggregated'], sep='\t') 
  mutations = mutations.rename(columns=rename_cols).drop(
      columns=['Center', 'Tumor_Seq_Allele1'])

  mutations[mutCol] = [str(i[0]) + ':' + str(i[1]) for i in np.nan_to_num(mutations[
      ['t_alt_count', 't_ref_count']].values, 0).astype(int)]
  mutations = mut.filterCoverage(mutations, loc=[mutCol], sep=':',cov=2)
  mutations = mut.filterAllelicFraction(mutations, loc=[mutCol], sep=':',frac=0.1)
  mutations = addAnnotation(mutations, NCBI_Build='37', Strand="+")
  mutations = annotateLikelyImmortalized(mutations,
                                          TCGAlocs=['TCGAhsCnt', 'COSMIChsCnt'], 
                                          max_recurrence=0.05, min_tcga_true_cancer=5)

  mutations.to_csv(save_output + 'somatic_mutations_all.csv', index=None)
  print('done')
  return mutations


def CCLEPostProcessing(wesrefworkspace=WESMUTWORKSPACE, wgsrefworkspace=WGSWORKSPACE, 
                      samplesetname=SAMPLESETNAME, todrop=KNOWN_DROP,
                       AllSamplesetName='all', doCleanup=False,
                       my_id=MY_ID, mystorage_id=MYSTORAGE_ID,
                       refsheet_url=REFSHEET_URL,
                       taiga_description=Mutationsreadme, taiga_dataset=TAIGA_MUTATION,
                       mutation_groups=MUTATION_GROUPS,
                       prev='ccle',
                       minfreqtocall=0.05,
                       **kwargs):
  """the full CCLE mutations post processing pipeline (used only by CCLE)

  see postprocess() to reproduce our analysis

  Args:
      wesrefworkspace ([type]): [description]
      wgsrefworkspace ([type]): [description]
      samplesetname ([type]): [description]
      AllSamplesetName (str, optional): [description]. Defaults to 'all'.
      doCleanup (bool, optional): [description]. Defaults to False.
      my_id ([type], optional): [description]. Defaults to MY_ID.
      mystorage_id ([type], optional): [description]. Defaults to MYSTORAGE_ID.
      refsheet_url ([type], optional): [description]. Defaults to REFSHEET_URL.
      taiga_description ([type], optional): [description]. Defaults to Mutationsreadme.
      taiga_dataset (str, optional): [description]. Defaults to TAIGA_MUTATION.
      mutation_groups ([type], optional): [description]. Defaults to MUTATION_GROUPS.
      prev ([type], optional): [description]. Defaults to tc.get(name=TAIGA_ETERNAL, file='CCLE_mutations').
      minfreqtocall (float, optional): [description]. Defaults to 0.05.
  """  
  if prev == 'ccle':
    prev = tc.get(name=TAIGA_ETERNAL, file='CCLE_mutations')
  if doCleanup:
    #TODO:
    val = ""
    #gcp.rmFiles('gs://fc-secure-012d088c-f039-4d36-bde5-ee9b1b76b912/$val/**/call-tumorMM_Task/*.cleaned.bam')
  # sometimes it does not work so better check again

  # doing wes
  print('doing wes')
  folder=os.path.join("temp", samplesetname, "wes_")
  
  wesmutations = postProcess(wesrefworkspace, AllSamplesetName if AllSamplesetName else samplesetname,
                             save_output=folder, doCleanup=True, mutCol="CGA_WES_AC", **kwargs)

  # renaming
  print('renaming')
  #wesrenaming = track.removeOlderVersions(names=set(
  #    wesmutations[SAMPLEID]), refsamples=wesrefwm.get_samples(),
  #    arxspan_id="arxspan_id", version="version", priority=priority)
  
  wesrenaming = h.fileToDict(folder+"sample_renaming.json")
  
  wesmutations = wesmutations[wesmutations[SAMPLEID].isin(wesrenaming.keys())].replace({
      SAMPLEID: wesrenaming})
  wesmutations.to_csv(folder + 'somatic_mutations_latest.csv', index=None)

  # doing wgs
  print('doing wgs')
  folder=os.path.join("temp", samplesetname, "wgs_")
  
  wgsmutations = postProcess(wgsrefworkspace, "allcurrent",  # AllSamplesetName if AllSamplesetName else samplesetname,
                         save_output=folder, doCleanup=True, mutCol="CGA_WES_AC", **kwargs)

  # renaming
  print('renaming')
  #wgsrenaming = track.removeOlderVersions(names=set(
  #    wesmutations[SAMPLEID]), refsamples=wgsrefwm.get_samples(),
  #    arxspan_id="arxspan_id", version="version", priority=priority)

  wgsrenaming = h.fileToDict(folder+"sample_renaming.json")

  wgsmutations = wgsmutations[wgsmutations[SAMPLEID].isin(wgsrenaming.keys())].replace({
      SAMPLEID: wgsrenaming})
  wgsmutations.to_csv(folder + 'somatic_mutations_latest.csv', index=None)

  # merge
  print('merging')
  folder = os.path.join("temp", samplesetname, "merged_")
  toadd = set(wgsmutations[SAMPLEID]) - set(wesmutations[SAMPLEID])
  priomutations = wesmutations.append(
      wgsmutations[wgsmutations[SAMPLEID].isin(toadd)]).reset_index(drop=True)
  #normals = set(ccle_refsamples[ccle_refsamples.primary_disease=="normal"].arxspan_id)
  #mutations = mutations[~mutations[SAMPLEID].isin(normals)]
  priomutations.to_csv(folder+'somatic_mutations.csv', index=False)

  #making binary mutation matrices
  print("creating mutation matrices")
  # binary mutations matrices
  mut.mafToMat(priomutations[(priomutations.isDeleterious)]).astype(
      int).T.to_csv(folder+'somatic_mutations_boolmatrix_deleterious.csv')
  mut.mafToMat(priomutations[~(priomutations.isDeleterious | priomutations.isCOSMIChotspot | 
                               priomutations.isTCGAhotspot | 
                               priomutations['Variant_Classification'] == 'Silent')]).astype(int).T.to_csv(
      folder+'somatic_mutations_boolmatrix_other.csv')
  mut.mafToMat(priomutations[(priomutations.isCOSMIChotspot | priomutations.isTCGAhotspot)]).astype(
    int).T.to_csv(folder+'somatic_mutations_boolmatrix_hotspot.csv')


  # genotyped mutations matrices
  mut.mafToMat(priomutations[(priomutations.isDeleterious)], mode="genotype",
              minfreqtocall=minfreqtocall).T.to_csv(folder+'somatic_mutations_matrix_deleterious.csv')
  mut.mafToMat(priomutations[~(priomutations.isDeleterious | priomutations.isCOSMIChotspot | 
                               priomutations.isTCGAhotspot |
                               priomutations['Variant_Classification'] == 'Silent')], 
                mode="genotype", minfreqtocall=minfreqtocall).T.to_csv(
                  folder+'somatic_mutations_matrix_other.csv')
  mut.mafToMat(priomutations[(priomutations.isCOSMIChotspot | priomutations.isTCGAhotspot)],
              mode="genotype", minfreqtocall=minfreqtocall).T.to_csv(
                folder+'somatic_mutations_matrix_hotspot.csv')

  # adding lgacy datasetss
  print('add legacy datasets')
  legacy_hybridcapture = tc.get(name='mutations-da6a', file='legacy_hybridcapture')
  legacy_raindance = tc.get(name='mutations-da6a', file='legacy_raindance')
  legacy_rna = tc.get(name='mutations-da6a', file='legacy_rna')
  legacy_wes_sanger = tc.get(name='mutations-da6a', file='legacy_wes_sanger')
  legacy_wgs_exoniconly = tc.get(name='mutations-da6a', file='legacy_wgs_exoniconly')

  merged = mut.mergeAnnotations(
      priomutations, legacy_hybridcapture, "HC_AC", useSecondForConflict=True, dry_run=False)
  merged = mut.mergeAnnotations(
    merged, legacy_raindance, "RD_AC", useSecondForConflict=True, dry_run=False)
  merged = mut.mergeAnnotations(
    merged, legacy_wgs_exoniconly, "WGS_AC", useSecondForConflict=False, dry_run=False)
  merged = mut.mergeAnnotations(
    merged, legacy_wes_sanger, "SangerWES_AC", useSecondForConflict=False, dry_run=False)
  merged = mut.mergeAnnotations(
    merged, legacy_rna, "RNAseq_AC", useSecondForConflict=False, dry_run=False)

  merged = merged[merged['tumor_f'] > 0.05]
  merged = annotateLikelyImmortalized(merged, TCGAlocs=[
                                                'TCGAhsCnt', 'COSMIChsCnt'], max_recurrence=0.05, 
                                                min_tcga_true_cancer=5)
  print("changing variant annotations")
  rename = {}
  for k,v in mutation_groups.items():
    for e in v:
      rename[e] = k
  merged['Variant_annotation'] = [rename[i] for i in merged['Variant_Classification'].tolist()]

  print('compare to previous release')
  a = set(merged[SAMPLEID])
  # tc.get(name='internal-20q2-7f46', version=18, file='CCLE_mutations')
  b = set(prev[SAMPLEID])
  print("new lines:")
  print(a-b)
  print('lost lines:')
  print(b-a)

  # making a depmap version
  # removing immortalized ffor now 
  merged = merged[~merged['is_likely_immortalization']]
  #reverting to previous versions
  merged = merged[MUTCOL_DEPMAP].rename(columns={
      "Tumor_Allele": "Tumor_Seq_Allele1"})
  merged.to_csv(folder+'somatic_mutations_withlegacy.csv', index=False)

  # making binary matrices
  merged = merged[merged['Entrez_Gene_Id'] != 0]
  merged['mutname'] = merged['Hugo_Symbol'] + " (" + merged["Entrez_Gene_Id"].astype(str) + ")"
  mut.mafToMat(merged[(merged.Variant_annotation=="damaging")], mode='bool', 
    mutNameCol="mutname").astype(int).T.to_csv(folder+'somatic_mutations_boolmatrix_fordepmap_damaging.csv')
  mut.mafToMat(merged[(merged.Variant_annotation=="other conserving")], mode='bool', 
    mutNameCol="mutname").astype(int).T.to_csv(folder+'somatic_mutations_boolmatrix_fordepmap_othercons.csv')
  mut.mafToMat(merged[(merged.Variant_annotation=="other non-conserving")], mode='bool', 
    mutNameCol="mutname").astype(int).T.to_csv(folder+'somatic_mutations_boolmatrix_fordepmap_othernoncons.csv')
  mut.mafToMat(merged[(merged.isCOSMIChotspot | merged.isTCGAhotspot)], mode='bool', 
    mutNameCol="mutname").astype(int).T.to_csv(folder+'somatic_mutations_boolmatrix_fordepmap_hotspot.csv')

  # uploading to taiga
  tc.update_dataset(changes_description="new "+samplesetname+" release!",
                    dataset_permaname=taiga_dataset,
                    upload_files=[
                      # for depmap
                        {
                            "path": folder+"somatic_mutations_boolmatrix_fordepmap_hotspot.csv",
                            "format": "NumericMatrixCSV",
                            "encoding": "utf-8"
                        },
                        {
                            "path": folder+"somatic_mutations_boolmatrix_fordepmap_othernoncons.csv",
                            "format": "NumericMatrixCSV",
                            "encoding": "utf-8"
                        },
                        {
                            "path": folder+"somatic_mutations_boolmatrix_fordepmap_damaging.csv",
                            "format": "NumericMatrixCSV",
                            "encoding": "utf-8"
                        },
                      # genotyped
                        {
                            "path": folder+"somatic_mutations_matrix_hotspot.csv",
                            "format": "NumericMatrixCSV",
                            "encoding": "utf-8"
                        },
                        {
                            "path": folder+"somatic_mutations_matrix_other.csv",
                            "format": "NumericMatrixCSV",
                            "encoding": "utf-8"
                        },
                        {
                            "path": folder+"somatic_mutations_matrix_deleterious.csv",
                            "format": "NumericMatrixCSV",
                            "encoding": "utf-8"
                        },
                      # new
                        {
                            "path": folder+"somatic_mutations_boolmatrix_fordepmap_hotspot.csv",
                            "format": "NumericMatrixCSV",
                            "encoding": "utf-8"
                        },
                        {
                            "path": folder+"somatic_mutations_boolmatrix_fordepmap_othernoncons.csv",
                            "format": "NumericMatrixCSV",
                            "encoding": "utf-8"
                        },
                        {
                            "path": folder+"somatic_mutations_boolmatrix_fordepmap_damaging.csv",
                            "format": "NumericMatrixCSV",
                            "encoding": "utf-8"
                        },
                        {
                            "path": folder+"somatic_mutations_withlegacy.csv",
                            "format": "TableCSV",
                            "encoding": "utf-8"
                        },
                        {
                            "path": folder+"somatic_mutations.csv",
                            "format": "TableCSV",
                            "encoding": "utf-8"
                        },
                        {
                            "path": 'temp/'+samplesetname+"/wes_somatic_mutations_all.csv",
                            "format": "TableCSV",
                            "encoding": "utf-8"
                        },
                        {
                            "path": 'temp/'+samplesetname+"/wgs_somatic_mutations_all.csv",
                            "format": "TableCSV",
                            "encoding": "utf-8"
                        },
                    ],
                    upload_async=False,
                    dataset_description=taiga_description)

def analyzeUnfiltered(workspace=WGSWORKSPACE, allsampleset='all', folder="temp/",
                      subsetcol=[SAMPLEID, 'Hugo_Symbol', 'Entrez_Gene_Id',
                                 'Chromosome', 'Start_position', 'End_position',
                                 'Variant_Classification', 'Variant_Type', 'Reference_Allele',
                                 'Tumor_Allele', 'dbSNP_RS', 'dbSNP_Val_Status', 'Genome_Change',
                                 'Annotation_Transcript', 'cDNA_Change', 'Codon_Change',
                                 'HGVS_protein_change',  'Protein_Change',
                                 't_alt_count', 't_ref_count', 'tumor_f', 'CGA_WES_AC'],
                      taiga_dataset=TAIGA_MUTATION,):

  print("retrieving unfiltered")
  ####### WES
  res = dm.WorkspaceManager(workspace).get_sample_sets()
  unfiltered = pd.read_csv(res.loc[allsampleset, 'unfiltered_CGA_MAF_aggregated'], sep='\t',
  encoding='L6',na_values=["__UNKNOWN__",'.'], engine='c', dtype=str)
  unfiltered['somatic'] = unfiltered['somatic'].replace('nan','False')
  unfiltered['HGNC_Status'] = unfiltered['HGNC_Status'].replace('nan','Unapproved')
  unfiltered['judgement'] = unfiltered['judgement'].replace('nan','REMOVE')
  unfiltered = unfiltered.rename(columns={"i_ExAC_AF":"ExAC_AF",
  "Tumor_Sample_Barcode":SAMPLEID,
  "Tumor_Seq_Allele2":"Tumor_Allele"}).drop(columns=['Tumor_Seq_Allele1'])
  unfiltered['CGA_WES_AC'] = [str(i[0]) + ':' + str(i[1]) for i in np.nan_to_num(
    unfiltered[['t_alt_count','t_ref_count']].values.astype(float),0).astype(int)]
  toremove = []
  subunfilt = unfiltered.iloc[:10000]
  for i, val in enumerate(unfiltered.columns):
    h.showcount(i,len(unfiltered.columns))
    if len(set(subunfilt[val])-set(['nan']))==1:
      if len(set(unfiltered[val])-set(['nan']))==1:
        toremove.append(val)
  unfiltered = unfiltered.drop(columns=set(toremove))
  toint =  ["Start_position", "End_position"]
  for val in toint:
    unfiltered[val]  = unfiltered[val].astype(int)
  unfiltered.to_csv(folder+'mutation_somatic_unfiltered_withreplicates.csv.gz', index=False)
  unfiltered = unfiltered[subsetcol]
  unfiltered.to_csv(folder+'mutation_somatic_unfiltered_withreplicates_subseted.csv.gz', index=False)
  os.system('gunzip '+folder+'mutation_somatic_unfiltered_withreplicates.csv.gz')
  del unfiltered
  tc.update_dataset(changes_description="adding unfiltered mutations",
                  dataset_permaname=taiga_dataset,
                  upload_files=[  
                    {
                        "path": folder+"mutation_somatic_unfiltered_withreplicates_subseted.csv",
                        "format": "TableCSV",
                        "encoding": "utf-8"
                    },
                    {
                        "path": folder+"mutation_somatic_unfiltered_withreplicates.csv",
                        "format": "TableCSV",
                        "encoding": "utf-8"
                    },
                  ],
                  add_all_existing_files=True,
                  upload_async=False)
