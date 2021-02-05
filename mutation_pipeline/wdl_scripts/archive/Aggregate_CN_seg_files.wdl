## DepMap workflow CN for HG38

This workspace contains the workflows used to realign WES samples to hg38 and call copy number for CCLE to be used in DepMap. This workspace contains annotations for all of the Broad WES samples as well as the Sanger WES samples. It is updated quarterly with the new samples that need to be processed for releases.

##### Copy number {#copy-number}

To generate the copy number dataset:



*   **BamToUnmappedRGBams_MC** vdauwera/BamToUnmappedRGBamsSnapshot ID: 3
*   **Generate_uBAM_File_List** gkugener/ArrayOfFilesToTxtSnapshot ID: 1
*   **Realign_WES_GATK4** gatk/PreProcessingForVariantDiscovery_GATK4Snapshot ID: 7
*   **CNV_sample_XX** gatk/CNV_Somatic_Pair_WorkflowSnapshot ID: 9
*   **Aggregate_CN_seg_files** gkugener/Aggregate_CN_seg_filesSnapshot ID: 2

This output file for download will be saved under the sample set under the combined_seg_file attribute.

There are several other tasks in this workspace. In brief:



*   **CNV_Somatic_Panel_Workflow_Agilent_XX** gatk/CNV_Somatic_Panel_WorkflowSnapshot ID: 11. This task was used in this workspace to generate the Sanger PON. In the Sanger dataset, there is a set of 40 normal cell lines samples (cell lines derived from matched normal tissue). We can use these to generate a PON to normalize to rather than using the Agilent PON we use for the other CCLE cell lines. This leads to less noisy results. HOWEVER, results using the PON from this workflow should not use the X chromosome, as the sanger normals are not exclusively female or male (it is likely a mix).
*   **SANGER_PON_CNV_sample_XX** gatk/CNV_Somatic_Pair_WorkflowSnapshot ID: 9. Same as the CNV_sample_XX_gatk, except that is uses the Sanger based PON. Should be used only for the Sanger cell lines.
*   **Sanger_PON_Aggregate_CN_seg_files** gkugener/Aggregate_CN_seg_filesSnapshot ID: 2. Aggregates the segment files for the samples that were run using the Sanger PON based CNV workflow.