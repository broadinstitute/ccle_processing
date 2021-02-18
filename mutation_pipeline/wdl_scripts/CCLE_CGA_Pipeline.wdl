# if there are any renamings in the directory structures, these paths need to be changed
import "https://raw.githubusercontent.com/BroadInstitute/ccle_processing/mutation_pipeline/CGA_WES_CCLE.wdl" as CGA_pipeline
import "https://raw.githubusercontent.com/BroadInstitute/ccle_processing/mutation_pipeline/common_variant_filter.wdl" as CommonVariantFilter
import "https://raw.githubusercontent.com/BroadInstitute/ccle_processing/mutation_pipeline/filterMAF_on_CGA.wdl" as filterMAF_on_CGA


workflow CCLE_CGA_Pipeline {
    call CGA_pipeline.CGA_Production_Analysis_Workflow {}
    call filter

    output {
        #### CGA_pipeline ######


        ####### QC Tasks Outputs #######
        # Copy Number QC Report files
        File tumor_bam_lane_list=CGA_pipeline.CopyNumberReportQC_Task.tumorBamLaneList
        File normal_bam_lane_list=CopyNumberReportQC_Task.normalBamLaneList
        File tumor_bam_read_coverage_lane=CopyNumberReportQC_Task.tumorRCL
        File normal_bam_read_coverage_lane=CopyNumberReportQC_Task.normalRCL
        File copy_number_qc_report=CopyNumberReportQC_Task.CopyNumQCReport
        File copy_number_qc_report_png=CopyNumberReportQC_Task.CopyNumQCReportPNG
        File copy_number_qc_mix_ups=CopyNumberReportQC_Task.CopyNumQCMixUps
        # Picard Multiple Metrics Task - NORMAL BAM
        File? normal_bam_alignment_summary_metrics=normalMM_Task.alignment_summary_metrics
        File? normal_bam_bait_bias_detail_metrics=normalMM_Task.bait_bias_detail_metrics
        File? normal_bam_bait_bias_summary_metrics=normalMM_Task.bait_bias_summary_metrics
        File? normal_bam_base_distribution_by_cycle=normalMM_Task.base_distribution_by_cycle
        File? normal_bam_base_distribution_by_cycle_metrics=normalMM_Task.base_distribution_by_cycle_metrics
        File? normal_bam_gc_bias_detail_metrics=normalMM_Task.gc_bias_detail_metrics
        File? normal_bam_gc_bias=normalMM_Task.gc_bias
        File? normal_bam_gc_bias_summary_metrics=normalMM_Task.gc_bias_summary_metrics
        File? normal_bam_insert_size_histogram=normalMM_Task.insert_size_histogram
        File? normal_bam_insert_size_metrics=normalMM_Task.insert_size_metrics
        File? normal_bam_pre_adapter_detail_metrics=normalMM_Task.pre_adapter_detail_metrics
        File? normal_bam_pre_adapter_summary_metrics=normalMM_Task.pre_adapter_summary_metrics
        File? normal_bam_quality_by_cycle=normalMM_Task.quality_by_cycle
        File? normal_bam_quality_by_cycle_metrics=normalMM_Task.quality_by_cycle_metrics
        File? normal_bam_quality_distribution=normalMM_Task.quality_distribution
        File? normal_bam_quality_distribution_metrics=normalMM_Task.quality_distribution_metrics
        File? normal_bam_quality_yield_metrics=normalMM_Task.quality_yield_metrics
        File? normal_bam_converted_oxog_metrics=normalMM_Task.converted_oxog_metrics
        File? normal_bam_hybrid_selection_metrics=normalMM_Task.hsMetrics
        # Picard Multiple Metrics Task - TUMOR BAM
        File? tumor_bam_alignment_summary_metrics=tumorMM_Task.alignment_summary_metrics
        File? tumor_bam_bait_bias_detail_metrics=tumorMM_Task.bait_bias_detail_metrics
        File? tumor_bam_bait_bias_summary_metrics=tumorMM_Task.bait_bias_summary_metrics
        File? tumor_bam_base_distribution_by_cycle=tumorMM_Task.base_distribution_by_cycle
        File? tumor_bam_base_distribution_by_cycle_metrics=tumorMM_Task.base_distribution_by_cycle_metrics
        File? tumor_bam_gc_bias_detail_metrics=tumorMM_Task.gc_bias_detail_metrics
        File? tumor_bam_gc_bias=tumorMM_Task.gc_bias
        File? tumor_bam_gc_bias_summary_metrics=tumorMM_Task.gc_bias_summary_metrics
        File? tumor_bam_insert_size_histogram=tumorMM_Task.insert_size_histogram
        File? tumor_bam_insert_size_metrics=tumorMM_Task.insert_size_metrics
        File? tumor_bam_pre_adapter_detail_metrics=tumorMM_Task.pre_adapter_detail_metrics
        File? tumor_bam_pre_adapter_summary_metrics=tumorMM_Task.pre_adapter_summary_metrics
        File? tumor_bam_quality_by_cycle=tumorMM_Task.quality_by_cycle
        File? tumor_bam_quality_by_cycle_metrics=tumorMM_Task.quality_by_cycle_metrics
        File? tumor_bam_quality_distribution=tumorMM_Task.quality_distribution
        File? tumor_bam_quality_distribution_metrics=tumorMM_Task.quality_distribution_metrics
        File? tumor_bam_quality_yield_metrics=tumorMM_Task.quality_yield_metrics
        File? tumor_bam_converted_oxog_metrics=tumorMM_Task.converted_oxog_metrics
        File? tumor_bam_hybrid_selection_metrics=tumorMM_Task.hsMetrics
        # Cross-Sample Contamination Task
        File contamination_data=ContEST_Task.contamDataFile
        File contestAFFile=ContEST_Task.contestAFFile
        File contest_base_report=ContEST_Task.contestBaseReport
        File contest_validation=ContEST_Task.validationOutput
        Float fracContam=ContEST_Task.fracContam
        # Cross Check Lane Fingerprints Task
        File cross_check_fingprt_metrics=CrossCheckLaneFingerprints_Task.crossCheckMetrics
        File cross_check_fingprt_report=CrossCheckLaneFingerprints_Task.crossCheckReport
        Float cross_check_fingprt_min_lod_value=CrossCheckLaneFingerprints_Task.crossCheckMinLODValue
        String cross_check_fingprt_min_lod_lanes=CrossCheckLaneFingerprints_Task.crossCheckMinLODLanes
        ####### Mutation Calling Tasks Outputs #######
        # MutectFC_Task
        File mutect_force_call_cs=MutectFC_Task.mutectfc_cs
        File mutect_force_call_pw=MutectFC_Task.mutectfc_pw
        File mutect_force_call_cw=MutectFC_Task.mutectfc_cw
        # Strelka
        File strelka_passed_indels=Strelka.Strelka_passed_indels
        File strelka_passed_snvs=Strelka.Strelka_passed_snvs
        File strelka_all_indels=Strelka.Strelka_all_indels
        File strelka_all_snvs=Strelka.Strelka_all_snvs
        # Gather MuTect1 power and coverage wiggle files
        File MuTect1_merged_power_wig=GatherWIGFiles_Task.MuTect1_merged_power_wig
        File MuTect1_merged_coverage_wig=GatherWIGFiles_Task.MuTect1_merged_coverage_wig
        # Gathered MuTect1 and MuTect2 calls stats
        File MUTECT1_CS_SNV=Gather_Task.MUTECT1_CS_SNV
        File MUTECT2_VCF_ALL=Gather_Task.MUTECT2_VCF_ALL
        File MUTECT2_VCF_INDELS=Gather_Task.MUTECT2_VCF_INDELS
        # deTiN (Tumor in Normal)
        Float? TiN=DeTiN_Task.TiN
        Int? deTiN_number_added_SSNVs=DeTiN_Task.number_added_SSNV
        String? TiN_CI=DeTiN_Task.TiN_CI
        File? deTiN_call_stats=DeTiN_Task.deTiN_call_stats
        File? deTiN_indels=DeTiN_Task.deTiN_indels
        File? deTiN_SSNVs_plot=DeTiN_Task.deTiN_SSNVs_plot
        File? deTiN_aSCNA_model=DeTiN_Task.aSCNA_model
        File? deTiN_aSCNA_kmeans_RSS_plot=DeTiN_Task.deTiN_aSCNA_kmeans_RSS_plot
        File? deTiN_aSCNA_scatter_plot=DeTiN_Task.deTiN_aSCNA_scatter_plot
        File? deTiN_TiN_modes_plot=DeTiN_Task.deTiN_TiN_modes_plot
        File? deTiN_segments=DeTiN_Task.deTiN_segments
        # Variant Effector Predictor Task
        File MUTECT1_VEP_annotated_vcf=VEP_Task.MUTECT1_VEP_annotated_vcf
        File MUTECT2_VEP_annotated_vcf=VEP_Task.MUTECT2_VEP_annotated_vcf
        File STRELKA_VEP_annotated_vcf=VEP_Task.STRELKA_VEP_annotated_vcf
        # Oncotator Output
        File mutect1_snv_mutect2_indel_strelka_indel_annotated_maf=Oncotate_Task.WXS_Mutation_M1_SNV_M2_INDEL_Strelka_INDEL_annotated_maf
        ####### Filtering Tasks Outputs #######
        # Orientation Bias Filter - OxoG
        Float oxoG_OBF_q_val=oxoGOBF.q_val
        File oxoG_OBF_figures=oxoGOBF.OBF_figures
        File oxoG_OBF_passed_mutations=oxoGOBF.WXS_Mutation_OBF_filtered_maf
        File oxoG_OBF_passed_and_rejected_mutations=oxoGOBF.WXS_Mutation_OBF_unfiltered_maf
        Int oxoG_OBF_number_mutations_passed=oxoGOBF.num_passed_mutations
        Int oxoG_OBF_number_mutations_rejected=oxoGOBF.num_rejected_mutations
        # Orientation Bias Filter - FFPE
        Float ffpe_OBF_q_val=ffpeOBF.q_val
        File ffpe_OBF_figures=ffpeOBF.OBF_figures
        File ffpe_OBF_passed_mutations=ffpeOBF.WXS_Mutation_OBF_filtered_maf
        File ffpe_OBF_passed_and_rejected_mutations=ffpeOBF.WXS_Mutation_OBF_unfiltered_maf
        Int ffpe_OBF_number_mutations_passed=ffpeOBF.num_passed_mutations
        Int ffpe_OBF_number_mutations_rejected=ffpeOBF.num_rejected_mutations
        # MAFPoNFilter
        Array[File] filter_passed_mutations=MAFPonFilter.passMaf
        Array[File] filter_passed_and_rejected_mutations=MAFPonFilter.allMaf
        Array[Int] number_mutations_passed=MAFPonFilter.num_passed_mutations
        Array[Int] number_mutations_rejected=MAFPonFilter.num_rejected_mutations
        # Blat Re-Aligner
        File blat_passed_mutations=blat.blat_results
        #File blat_debug_results=blat.debug_results
        File blat_all_maf=blat.blat_all_maf
        Int blat_number_mutations_passed=blat.num_passed_mutations
        Int blat_number_mutations_rejected=blat.num_rejected_mutations
        # Merge MAF File Task
        File filters_passed_merged_intersection_maf=merge_mafs_task.merged_intersection_maf
        File filters_passed_merged_union_maf=merge_mafs_task.merged_union_maf
        # Mutation Validator
        #File mutation_validator_pileup_preprocessing=mutation_validator.pileup_preprocessing_txt
        File mutation_validator_validated_maf=mutation_validator.validated_maf
        ####### Copy Number - GATK CNV & Allelic CapSeg #######
        File gatk_cnv_coverage_file=gatk_acnv_only.gatk_cnv_coverage_file
        File gatk_cnv_seg_file=gatk_acnv_only.gatk_cnv_seg_file
        File gatk_cnv_tn_coverage=gatk_acnv_only.gatk_cnv_tn_coverage
        File gatk_cnv_pre_tn_coverage=gatk_acnv_only.gatk_cnv_pre_tn_coverage
        File gatk_het_ad_normal=gatk_acnv_only.gatk_het_ad_normal
        File gatk_het_ad_tumor=gatk_acnv_only.gatk_het_ad_tumor
        Array[File] gatk_cnv_all_plots=gatk_acnv_only.gatk_cnv_all_plots
        File alleliccapseg_plot=gatk_acnv_only.alleliccapseg_plot
        File alleliccapseg_tsv=gatk_acnv_only.alleliccapseg_tsv
        Float alleliccapseg_skew=gatk_acnv_only.alleliccapseg_skew
        ####### Absolute #######
        File? absolute_highres_plot=absolute.absolute_highres_plot
        File? absolute_rdata=absolute.absolute_rdata
        ####### Lego Plot ######
        Array[File]? lego_plotter_ais=lego_plotter_task.ais
        Array[File]? lego_plotter_pngs=lego_plotter_task.pngs
        Array[File]? lego_plotter_figs=lego_plotter_task.figs
        Array[File]? lego_plotter_pss=lego_plotter_task.pss
        File? mut_legos_html=lego_plotter_task.mut_legos_html


        ##### CommonVariantFilter
        File commonfilter_annotated_maf = commonfilterTask.annotatedMAF
        File commonfilter_passed_maf = commonfilterTask.passedMAF
        File commonfilter_rejected_maf = commonfilterTask.rejectedMAF
        Int commonfilter_considered_count = commonfilterTask.consideredCount
        Int commonfilter_pass_count = commonfilterTask.passCount
        Int commonfilter_reject_count = commonfilterTask.rejectCount
    }
}
