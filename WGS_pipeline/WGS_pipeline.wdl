import "https://raw.githubusercontent.com/broadinstitute/ccle_processing/master/WGS_pipeline/BamToUnmappedRGBams.wdl" as BamToUnmappedRGBams 
import "https://raw.githubusercontent.com/broadinstitute/ccle_processing/master/WGS_pipeline/ArrayOfFilesToTxt.wdl" as ArrayOfFilesToTxt 
import "https://raw.githubusercontent.com/broadinstitute/ccle_processing/master/WGS_pipeline/PreProcessingForVariantDiscovery_GATK4.wdl" as PreProcessingForVariantDiscovery_GATK4 
import "https://raw.githubusercontent.com/broadinstitute/ccle_processing/master/WGS_pipeline/CNV_Somatic_Workflow_on_Sample.wdl" as CNV_Somatic_Workflow_on_Sample 
import "https://raw.githubusercontent.com/broadinstitute/ccle_processing/master/WGS_pipeline/cnn-variant-filter.wdl" as cnn_variant_filter 
import "https://raw.githubusercontent.com/broadinstitute/ccle_processing/master/WGS_pipeline/Manta_SomaticSV.wdl" as Manta_SomaticSV 
import "https://raw.githubusercontent.com/broadinstitute/ccle_processing/master/WGSmut_pipeline/CGA_WES_CCLE_Characterization_Pipeline_v0.1_Jul2019_copy.wdl" as CGA_WES_CCLE_Characterization_Pipeline_v0 
import "https://raw.githubusercontent.com/broadinstitute/ccle_processing/master/WGSmut_pipeline/common_variant_filter.wdl" as common_variant_filter 
import "https://raw.githubusercontent.com/broadinstitute/ccle_processing/master/WGSmut_pipeline/filterMAF_on_CGA_pipeline.wdl" as filterMAF_on_CGA_pipeline 


workflow WGS_pipeline {

	#BamToUnmappedRGBams
	#"broadinstitute/genomes-in-the-cloud:2.3.1-1504795437"
	File ref_fasta
	File ref_fasta_index

	File input_bam
	File input_bam_index
	
	String picard_path
	String picard_docker

	Int preemptible_tries

	#ArrayOfFilesToTxt
	String sample_name #this.name

	#PreProcessingForVariantDiscovery_GATK4
	#"broadinstitute/gatk:4.beta.3"
	#"broadinstitute/genomes-in-the-cloud:2.3.0-1501082129"
	String ref_name

	File ref_dict

	File dbSNP_vcf
	File dbSNP_vcf_index

	Array[File] known_indels_sites_VCFs
	Array[File] known_indels_sites_indices

	String gotc_docker
	String gatk_docker
	String python_docker

	String gotc_path

	#CNV_Somatic_Workflow_on_Sample
	#us.gcr.io/broad-gatk/gatk:4.1.5.0
	File common_sites
	File intervals

	#cnn_variant_filter
	#"broadinstitute/gatk:4.1.4.0"
	# Manta

	#CGA_Production_Analysis_Workflow
	#gs://ccle_default_params/references/gatk-package-4.0.5.1-local.jar
	#gs://getzlab-workflows-reference_files-oa/gatk-protected.jar
	File ref_fasta_H19
	File ref_fasta_H19_index
	File ref_fasta_H19_dict

	call BamToUnmappedRGBams.BamToUnmappedRGBamsWf as BamToUnmappedRGBamsWf {
		input:
			ref_fasta=ref_fasta,
			ref_fasta_index=ref_fasta_index,
			input_bam=input_bam,
			picard_path=picard_path,
			preemptible_tries=preemptible_tries 
	}

	call ArrayOfFilesToTxt.CreateTxt as CreateTxt {
		input:
			array_of_files=BamToUnmappedRGBamsWf.sortsam_out,
			list_name=sample_name
	}

	call PreProcessingForVariantDiscovery_GATK4.PreProcessingForVariantDiscovery_GATK4 as PreProcessingForVariantDiscovery_GATK4 {
		input:
			sample_name=sample_name,
			ref_name=ref_name,
			flowcell_unmapped_bams_list=CreateTxt.file_list_name,
			ref_dict=ref_dict,
			ref_fasta=ref_fasta,
			ref_fasta_index=ref_fasta_index,
			dbSNP_vcf=dbSNP_vcf,
			dbSNP_vcf_index=dbSNP_vcf_index,
			known_indels_sites_VCFs=known_indels_sites_VCFs,
			known_indels_sites_indices=known_indels_sites_indices,
			gotc_docker=gotc_docker,
			python_docker=python_docker,
			gotc_path=gotc_path,
			picard_path=picard_path,
			preemptible_tries=preemptible_tries
	}

	call CNV_Somatic_Workflow_on_Sample.CNVSomaticPairWorkflow as CNVSomaticPairWorkflow {
		input:
			common_sites=common_sites,
			intervals=intervals,
			ref_fasta=ref_fasta,
			ref_fasta_dict=ref_dict,
			ref_fasta_fai=ref_fasta_index,
			tumor_bam=PreProcessingForVariantDiscovery_GATK4.analysis_ready_bam,
			tumor_bam_idx=PreProcessingForVariantDiscovery_GATK4.analysis_ready_bam_index
	}

	call cnn_variant_filter.Cram2FilteredVcf as Cram2FilteredVcf {
		input:
			input_file=PreProcessingForVariantDiscovery_GATK4.analysis_ready_bam,
			input_file_index=PreProcessingForVariantDiscovery_GATK4.analysis_ready_bam_index,
			reference_fasta=ref_fasta,
			reference_dict=ref_dict,
			reference_fasta_index=ref_fasta_index,
			output_prefix=sample_name
	}

	call Manta_SomaticSV.MantaSomaticSV as MantaSomaticSV {
		input:
			sample_name=sample_name,
			tumor_bam=PreProcessingForVariantDiscovery_GATK4.analysis_ready_bam,
			tumor_bam_index=PreProcessingForVariantDiscovery_GATK4.analysis_ready_bam_index,
			ref_fasta=ref_fasta,
			ref_fasta_index=ref_fasta_index,
			is_cram=false
	}

	call CGA_WES_CCLE_Characterization_Pipeline_v0.CGA_Production_Analysis_Workflow as CGA_Production_Analysis_Workflow {
		input:
			tumorBam=input_bam,
			tumorBamIdx=input_bam_index,
			pairName=sample_name,
			refFasta=ref_fasta_H19,
			refFastaIdx=ref_fasta_H19_index,
			refFastaDict=ref_fasta_H19_dict,
	}

	call common_variant_filter.CommonVariantFilter as CommonVariantFilter {
		input:
			sampleId=sample_name,
			maf=CGA_Production_Analysis_Workflow.mutation_validator_validated_maf
	}

	call filterMAF_on_CGA_pipeline.filterMaf as filterMaf {
		input:
			sample_id=sample_name,
			inMAFfn=CommonVariantFilter.commonfilter_passed_maf
	}

	output {
		#BamToUnmappedRGBams
		#createTxT
		#PreProcessingForVariantDiscovery_GATK4
		File duplication_metrics = PreProcessingForVariantDiscovery_GATK4.duplication_metrics
		File bqsr_report = PreProcessingForVariantDiscovery_GATK4.bqsr_report
		File analysis_ready_bam = PreProcessingForVariantDiscovery_GATK4.analysis_ready_bam
		File analysis_ready_bam_index = PreProcessingForVariantDiscovery_GATK4.analysis_ready_bam_index
		File analysis_ready_bam_md5 = PreProcessingForVariantDiscovery_GATK4.analysis_ready_bam_md5
		#CNVSomaticPairWorkflow
		File preprocessed_intervals = CNVSomaticPairWorkflow.preprocessed_intervals
		File read_counts_entity_id_tumor = CNVSomaticPairWorkflow.read_counts_entity_id_tumor
		File read_counts_tumor = CNVSomaticPairWorkflow.read_counts_tumor
		File allelic_counts_entity_id_tumor = CNVSomaticPairWorkflow.allelic_counts_entity_id_tumor
		File allelic_counts_tumor = CNVSomaticPairWorkflow.allelic_counts_tumor
		File denoised_copy_ratios_tumor = CNVSomaticPairWorkflow.denoised_copy_ratios_tumor
		File standardized_copy_ratios_tumor = CNVSomaticPairWorkflow.standardized_copy_ratios_tumor
		File het_allelic_counts_tumor = CNVSomaticPairWorkflow.het_allelic_counts_tumor
		File normal_het_allelic_counts_tumor = CNVSomaticPairWorkflow.normal_het_allelic_counts_tumor
		File copy_ratio_only_segments_tumor = CNVSomaticPairWorkflow.copy_ratio_only_segments_tumor
		File copy_ratio_legacy_segments_tumor = CNVSomaticPairWorkflow.copy_ratio_legacy_segments_tumor
		File allele_fraction_legacy_segments_tumor = CNVSomaticPairWorkflow.allele_fraction_legacy_segments_tumor
		File modeled_segments_begin_tumor = CNVSomaticPairWorkflow.modeled_segments_begin_tumor
		File copy_ratio_parameters_begin_tumor = CNVSomaticPairWorkflow.copy_ratio_parameters_begin_tumor
		File allele_fraction_parameters_begin_tumor = CNVSomaticPairWorkflow.allele_fraction_parameters_begin_tumor
		File modeled_segments_tumor = CNVSomaticPairWorkflow.modeled_segments_tumor
		File copy_ratio_parameters_tumor = CNVSomaticPairWorkflow.copy_ratio_parameters_tumor
		File allele_fraction_parameters_tumor = CNVSomaticPairWorkflow.allele_fraction_parameters_tumor
		File called_copy_ratio_segments_tumor = CNVSomaticPairWorkflow.called_copy_ratio_segments_tumor
		File called_copy_ratio_legacy_segments_tumor = CNVSomaticPairWorkflow.called_copy_ratio_legacy_segments_tumor
		File denoised_copy_ratios_plot_tumor = CNVSomaticPairWorkflow.denoised_copy_ratios_plot_tumor
		File denoised_copy_ratios_lim_4_plot_tumor = CNVSomaticPairWorkflow.denoised_copy_ratios_lim_4_plot_tumor
		File standardized_MAD_tumor = CNVSomaticPairWorkflow.standardized_MAD_tumor
		Float standardized_MAD_value_tumor = CNVSomaticPairWorkflow.standardized_MAD_value_tumor
		File denoised_MAD_tumor = CNVSomaticPairWorkflow.denoised_MAD_tumor
		Float denoised_MAD_value_tumor = CNVSomaticPairWorkflow.denoised_MAD_value_tumor
		File delta_MAD_tumor = CNVSomaticPairWorkflow.delta_MAD_tumor
		Float delta_MAD_value_tumor = CNVSomaticPairWorkflow.delta_MAD_value_tumor
		File scaled_delta_MAD_tumor = CNVSomaticPairWorkflow.scaled_delta_MAD_tumor
		Float scaled_delta_MAD_value_tumor = CNVSomaticPairWorkflow.scaled_delta_MAD_value_tumor
		File modeled_segments_plot_tumor = CNVSomaticPairWorkflow.modeled_segments_plot_tumor
		File oncotated_called_file_tumor = CNVSomaticPairWorkflow.oncotated_called_file_tumor
		File oncotated_called_gene_list_file_tumor = CNVSomaticPairWorkflow.oncotated_called_gene_list_file_tumor
		#cnn-variant-filter
		File cnn_filtered_vcf = Cram2FilteredVcf.cnn_filtered_vcf
		File cnn_filtered_vcf_index= Cram2FilteredVcf.cnn_filtered_vcf_index
		#MantaSomaticSV
		File candidate_indel_vcf= MantaSomaticSV.candidate_indel_vcf
		File candidate_indel_vcf_index= MantaSomaticSV.candidate_indel_vcf_index
		File candidate_sv_vcf= MantaSomaticSV.candidate_sv_vcf
		File candidate_sv_vcf_index= MantaSomaticSV.candidate_sv_vcf_index
		File germline_sv_vcf= MantaSomaticSV.germline_sv_vcf
		File germline_sv_vcf_index= MantaSomaticSV.germline_sv_vcf_index
		File somatic_sv_vcf= MantaSomaticSV.somatic_sv_vcf
		File somatic_sv_vcf_index= MantaSomaticSV.somatic_sv_vcf_index
		#CGA_Production_Analysis_Workflow
		####### QC Tasks Outputs #######
		# Copy Number QC Report files
		File tumor_bam_lane_list = CGA_Production_Analysis_Workflow.tumor_bam_lane_list
		File tumor_bam_read_coverage_lane = CGA_Production_Analysis_Workflow.tumor_bam_read_coverage_lane
		File copy_number_qc_report = CGA_Production_Analysis_Workflow.copy_number_qc_report
		File copy_number_qc_report_png = CGA_Production_Analysis_Workflow.copy_number_qc_report_png
		File copy_number_qc_mix_ups = CGA_Production_Analysis_Workflow.copy_number_qc_mix_ups
		# Picard Multiple Metrics Task - TUMOR BAM
		File? tumor_bam_bam_validation = CGA_Production_Analysis_Workflow.tumor_bam_bam_validation
		File? tumor_bam_alignment_summary_metrics = CGA_Production_Analysis_Workflow.tumor_bam_alignment_summary_metrics
		File? tumor_bam_bait_bias_detail_metrics = CGA_Production_Analysis_Workflow.tumor_bam_bait_bias_detail_metrics
		File? tumor_bam_bait_bias_summary_metrics = CGA_Production_Analysis_Workflow.tumor_bam_bait_bias_summary_metrics
		File? tumor_bam_base_distribution_by_cycle = CGA_Production_Analysis_Workflow.tumor_bam_base_distribution_by_cycle
		File? tumor_bam_base_distribution_by_cycle_metrics = CGA_Production_Analysis_Workflow.tumor_bam_base_distribution_by_cycle_metrics
		File? tumor_bam_gc_bias_detail_metrics = CGA_Production_Analysis_Workflow.tumor_bam_gc_bias_detail_metrics
		File? tumor_bam_gc_bias = CGA_Production_Analysis_Workflow.tumor_bam_gc_bias
		File? tumor_bam_gc_bias_summary_metrics = CGA_Production_Analysis_Workflow.tumor_bam_gc_bias_summary_metrics
		File? tumor_bam_insert_size_histogram = CGA_Production_Analysis_Workflow.tumor_bam_insert_size_histogram
		File? tumor_bam_insert_size_metrics = CGA_Production_Analysis_Workflow.tumor_bam_insert_size_metrics
		File? tumor_bam_pre_adapter_detail_metrics = CGA_Production_Analysis_Workflow.tumor_bam_pre_adapter_detail_metrics
		File? tumor_bam_pre_adapter_summary_metrics = CGA_Production_Analysis_Workflow.tumor_bam_pre_adapter_summary_metrics
		File? tumor_bam_quality_by_cycle = CGA_Production_Analysis_Workflow.tumor_bam_quality_by_cycle
		File? tumor_bam_quality_by_cycle_metrics = CGA_Production_Analysis_Workflow.tumor_bam_quality_by_cycle_metrics
		File? tumor_bam_quality_distribution = CGA_Production_Analysis_Workflow.tumor_bam_quality_distribution
		File? tumor_bam_quality_distribution_metrics = CGA_Production_Analysis_Workflow.tumor_bam_quality_distribution_metrics
		File? tumor_bam_quality_yield_metrics = CGA_Production_Analysis_Workflow.tumor_bam_quality_yield_metrics
		File? tumor_bam_converted_oxog_metrics = CGA_Production_Analysis_Workflow.tumor_bam_converted_oxog_metrics
		File? tumor_bam_hybrid_selection_metrics = CGA_Production_Analysis_Workflow.tumor_bam_hybrid_selection_metrics
		# Cross-Sample Contamination Task
		File contamination_data = CGA_Production_Analysis_Workflow.contamination_data
		File contestAFFile = CGA_Production_Analysis_Workflow.contestAFFile
		File contest_base_report = CGA_Production_Analysis_Workflow.contest_base_report
		File contest_validation = CGA_Production_Analysis_Workflow.contest_validation
		Float fracContam = CGA_Production_Analysis_Workflow.fracContam
		# Cross Check Lane Fingerprints Task
		File cross_check_fingprt_metrics = CGA_Production_Analysis_Workflow.cross_check_fingprt_metrics
		File cross_check_fingprt_report = CGA_Production_Analysis_Workflow.cross_check_fingprt_report
		Float cross_check_fingprt_min_lod_value = CGA_Production_Analysis_Workflow.cross_check_fingprt_min_lod_value
		String cross_check_fingprt_min_lod_lanes = CGA_Production_Analysis_Workflow.cross_check_fingprt_min_lod_lanes
		####### Mutation Calling Tasks Outputs #######
		# MutectFC_Task
		File mutect_force_call_cs = CGA_Production_Analysis_Workflow.mutect_force_call_cs
		File mutect_force_call_pw = CGA_Production_Analysis_Workflow.mutect_force_call_pw
		File mutect_force_call_cw = CGA_Production_Analysis_Workflow.mutect_force_call_cw
		# Strelka
		File strelka_passed_indels = CGA_Production_Analysis_Workflow.strelka_passed_indels
		File strelka_passed_snvs = CGA_Production_Analysis_Workflow.strelka_passed_snvs
		File strelka_all_indels = CGA_Production_Analysis_Workflow.strelka_all_indels
		File strelka_all_snvs = CGA_Production_Analysis_Workflow.strelka_all_snvs
		# Gather MuTect1 power and coverage wiggle files
		File MuTect1_merged_power_wig = CGA_Production_Analysis_Workflow.MuTect1_merged_power_wig
		File MuTect1_merged_coverage_wig = CGA_Production_Analysis_Workflow.MuTect1_merged_coverage_wig
		# Gathered MuTect1 and MuTect2 calls stats
		File MUTECT1_CS_SNV = CGA_Production_Analysis_Workflow.MUTECT1_CS_SNV
		File MUTECT2_VCF_ALL = CGA_Production_Analysis_Workflow.MUTECT2_VCF_ALL
		File MUTECT2_VCF_INDELS = CGA_Production_Analysis_Workflow.MUTECT2_VCF_INDELS        
		# Variant Effector Predictor Task
		File MUTECT1_VEP_annotated_vcf = CGA_Production_Analysis_Workflow.MUTECT1_VEP_annotated_vcf
		File MUTECT2_VEP_annotated_vcf = CGA_Production_Analysis_Workflow.MUTECT2_VEP_annotated_vcf
		File STRELKA_VEP_annotated_vcf = CGA_Production_Analysis_Workflow.STRELKA_VEP_annotated_vcf
		# Oncotator Output
		File mutect1_snv_mutect2_indel_strelka_indel_annotated_maf = CGA_Production_Analysis_Workflow.mutect1_snv_mutect2_indel_strelka_indel_annotated_maf       
		####### Filtering Tasks Outputs #######
		# Orientation Bias Filter - OxoG
		Float oxoG_OBF_q_val = CGA_Production_Analysis_Workflow.oxoG_OBF_q_val
		File oxoG_OBF_figures = CGA_Production_Analysis_Workflow.oxoG_OBF_figures
		File oxoG_OBF_passed_mutations = CGA_Production_Analysis_Workflow.oxoG_OBF_passed_mutations
		File oxoG_OBF_passed_and_rejected_mutations = CGA_Production_Analysis_Workflow.oxoG_OBF_passed_and_rejected_mutations
		Int oxoG_OBF_number_mutations_passed = CGA_Production_Analysis_Workflow.oxoG_OBF_number_mutations_passed
		Int oxoG_OBF_number_mutations_rejected = CGA_Production_Analysis_Workflow.oxoG_OBF_number_mutations_rejected
		# Orientation Bias Filter - FFPE
		Float ffpe_OBF_q_val = CGA_Production_Analysis_Workflow.ffpe_OBF_q_val
		File ffpe_OBF_figures = CGA_Production_Analysis_Workflow.ffpe_OBF_figures
		File ffpe_OBF_passed_mutations = CGA_Production_Analysis_Workflow.ffpe_OBF_passed_mutations
		File ffpe_OBF_passed_and_rejected_mutations = CGA_Production_Analysis_Workflow.ffpe_OBF_passed_and_rejected_mutations
		Int ffpe_OBF_number_mutations_passed = CGA_Production_Analysis_Workflow.ffpe_OBF_number_mutations_passed
		Int ffpe_OBF_number_mutations_rejected = CGA_Production_Analysis_Workflow.ffpe_OBF_number_mutations_rejected
		# MAFPoNFilter
		Array[File] filter_passed_mutations = CGA_Production_Analysis_Workflow.filter_passed_mutations
		Array[File] filter_passed_and_rejected_mutations = CGA_Production_Analysis_Workflow.filter_passed_and_rejected_mutations
		Array[Int] number_mutations_passed = CGA_Production_Analysis_Workflow.number_mutations_passed
		Array[Int] number_mutations_rejected = CGA_Production_Analysis_Workflow.number_mutations_rejected
		# Blat Re-Aligner
		File blat_passed_mutations = CGA_Production_Analysis_Workflow.blat_passed_mutations
		#File blat_debug_results=blat.debug_results
		File blat_all_maf = CGA_Production_Analysis_Workflow.blat_all_maf
		Int blat_number_mutations_passed = CGA_Production_Analysis_Workflow.blat_number_mutations_passed
		Int blat_number_mutations_rejected = CGA_Production_Analysis_Workflow.blat_number_mutations_rejected
		# Merge MAF File Task
		File filters_passed_merged_intersection_maf = CGA_Production_Analysis_Workflow.filters_passed_merged_intersection_maf
		File filters_passed_merged_union_maf = CGA_Production_Analysis_Workflow.filters_passed_merged_union_maf
		# Mutation Validator
		#File mutation_validator_pileup_preprocessing=mutation_validator.pileup_preprocessing_txt
		File mutation_validator_validated_maf = CGA_Production_Analysis_Workflow.mutation_validator_validated_maf
		####### Copy Number - GATK CNV & Allelic CapSeg #######
		File gatk_cnv_tn_coverage = CGA_Production_Analysis_Workflow.gatk_cnv_tn_coverage
		File gatk_cnv_pre_tn_coverage = CGA_Production_Analysis_Workflow.gatk_cnv_pre_tn_coverage
		File gatk_het_ad_tumor = CGA_Production_Analysis_Workflow.gatk_het_ad_tumor 
		File alleliccapseg_plot = CGA_Production_Analysis_Workflow.alleliccapseg_plot
		File alleliccapseg_tsv = CGA_Production_Analysis_Workflow.alleliccapseg_tsv
		Float alleliccapseg_skew = CGA_Production_Analysis_Workflow.alleliccapseg_skew
		####### Absolute #######
		File? absolute_highres_plot = CGA_Production_Analysis_Workflow.absolute_highres_plot
		File? absolute_rdata = CGA_Production_Analysis_Workflow.absolute_rdata
		####### Lego Plot ######
		Array[File]? lego_plotter_ais = CGA_Production_Analysis_Workflow.lego_plotter_ais
		Array[File]? lego_plotter_pngs = CGA_Production_Analysis_Workflow.lego_plotter_pngs
		Array[File]? lego_plotter_figs = CGA_Production_Analysis_Workflow.lego_plotter_figs
		Array[File]? lego_plotter_pss = CGA_Production_Analysis_Workflow.lego_plotter_pss
		File? mut_legos_html = CGA_Production_Analysis_Workflow.mut_legos_html
		# CommonVariantFilter
		File annotatedMAF=CommonVariantFilter.commonfilter_annotated_maf
		File passedMAF=CommonVariantFilter.commonfilter_passed_maf
		File rejectedMAF=CommonVariantFilter.commonfilter_rejected_maf
		String consideredCount=CommonVariantFilter.commonfilter_considered_count
		String passCount=CommonVariantFilter.commonfilter_pass_count
		String rejectCount=CommonVariantFilter.commonfilter_reject_count
		#filterMaf
		File outMAFannotatedFN=filterMaf.outMAFannotatedFN
		File outMAFfn=filterMaf.outMAFfn
	}
}