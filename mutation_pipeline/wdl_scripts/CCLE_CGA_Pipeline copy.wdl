# if there are any renamings in the directory structures, these paths need to be changed
import "https://raw.githubusercontent.com/BroadInstitute/ccle_processing/mutation_pipeline/CGA_WES_CCLE.wdl" as CGA_pipeline
import "https://raw.githubusercontent.com/BroadInstitute/ccle_processing/mutation_pipeline/common_variant_filter.wdl" as CommonVariantFilter
import "https://raw.githubusercontent.com/BroadInstitute/ccle_processing/mutation_pipeline/filterMAF_on_CGA.wdl" as filterMAF_on_CGA


workflow CCLE_CGA_Pipeline {
    ## GETZ LAB CGA WES CHARACTERIZATION WORKFLOW
    ## Copyright (c) 2017-2018, Broad Institute, Inc. and The General Hospital Corporation. All rights reserved.
    ## Copyright (c) 2017-2018, Contributors and authors of the pipeline and WDL code: Liudmila Elagina,
    ## Ignaty Leshchiner, Chip Stewart,  Chet Birger,  Ruslana Frazer, Eddie Salinas, Gad Getz.
    ## All rights reserved.
    ##
    ## LICENSING :
    ## This script is released under the CGA WES Characterization License (see
    ## https://docs.google.com/document/d/1D3UYYIBk-iIicrYhzqOVt4B51Apd_EZr9WCJuuaCJis).
    ## Note that the programs it calls may be subject to different licenses.  Users are
    ## responsible for checking that they are authorized to run all programs before running
    ## this script.  Please see the CGA WES Characterization License for a list of the
    ## programs called by this script and the locations of their respective licenses.



    # WORKFLOW INPUT PARAMS
    # Configuration file with optional parameters (json format)
    File cga_pipeline_config

    # Pair Input
    # sample tumor BAM file (see https://samtools.github.io/hts-specs/SAMv1.pdf)
    File tumorBam
    # sample normal BAM file (see https://samtools.github.io/hts-specs/SAMv1.pdf)
    File normalBam
    # sample normal BAI file (BAM indexed) (see samtools index command http://www.htslib.org/doc/samtools.html)
    File tumorBamIdx
    # sample normal BAI file (BAM indexed) (see samtools index command http://www.htslib.org/doc/samtools.html)
    File normalBamIdx
    # a string for the name of the pair under analysis used for naming output files
    String pairName
    # a string for the name of the tumor sample under analysis used for naming output files
    String caseName
    # a string for the name of the normal sample under analysis used for naming output files
    String ctrlName
    # Reference Files
    # list of read groups to exclude from the analysis in MuTect1 and MuTect_FC tasks
    File readGroupBlackList
    # the FASTA file for the appropriate genome build (Reference sequence file)
    File refFasta
    # the FASTA file index for the reference genome (see http://www.htslib.org/doc/faidx.html)
    File refFastaIdx
    # the FASTA file dictionary for the reference genome (see https://broadinstitute.github.io/picard/command-line-overview.html#CreateSequenceDictionary)
    File refFastaDict
    # an interval list file that contains the locations of the targets
    File targetIntervals
    # an interval list file that contains the locations of the baits used
    File baitIntervals
    # VCF format dbSNP file, used to exclude regions around known polymorphisms from analysis by some PROGRAMs;
    # PROGRAMs whose CLP doesn't allow for this argument will quietly ignore it
    File DB_SNP_VCF
    # index file of VCF file of DB SNP variants
    File DB_SNP_VCF_IDX
    # catalogue of somatic mutations in VCF format
    File cosmicVCF
    # panel of normals
    File MuTectNormalPanel
    # TSV file of chromsomal annotation ; chr, start, end, band, stain
    File cytoBandFile
    # GATK Jar file
    File GATK4_JAR

    # Loading optional parameters from configuration file
    Map[String, String] runtime_params = read_json(cga_pipeline_config)

    # List of PONs for MAF filtering in task MafPonFilter
    File PONs_list
    Array[Object] PONs_data = read_objects(PONs_list)

    # COMPUTE FILE SIZE
    Int tumorBam_size   = ceil(size(tumorBam,   "GB") + size(tumorBamIdx,    "GB"))
    Int normalBam_size  = ceil(size(normalBam,  "GB") + size(normalBamIdx,   "GB"))
    Int db_snp_vcf_size = ceil(size(DB_SNP_VCF, "GB") + size(DB_SNP_VCF_IDX, "GB"))
    Int refFasta_size   = ceil(size(refFasta,   "GB") + size(refFastaDict,   "GB") + size(refFastaIdx, "GB"))
    Int gatk4_jar_size  = ceil(size(GATK4_JAR,  "GB"))

    # Does the sample already have picard metrics computed
    Boolean hasPicardMetrics_tumor           = false
    Boolean hasPicardMetrics_normal          = false
    # Should we compute picard metrics anyway, even if they exist
    Boolean forceComputePicardMetrics_tumor  = true
    Boolean forceComputePicardMetrics_normal = true
    # Avoids running DeTiN task when no matched normal is provided
    Boolean runDeTiN = false

    ## Common variant filter
    meta {
        author: "Brendan Reardon"
        email: "breardon@broadinstitute.org"
        laboratory: "Van Allen Lab"
        institution: "Dana-Farber Cancer Institute, Broad Institute of MIT & Harvard"
        github: "https://github.com/vanallenlab/common_variant_filter"
        license: "MIT License"
    }

    String CommonVariantFilter.sampleId

    Int CommonVariantFilter.min_exac_ac = 10
    Int CommonVariantFilter.min_filter_depth = 0
    Boolean? CommonVariantFilter.disable_whitelist = false
    Boolean? CommonVariantFilter.filter_noncoding = false

    Int CommonVariantFilter.RAM = 8
    Int CommonVariantFilter.SSD = 60
    Int CommonVariantFilter.preemptible = 3

    String CommonVariantFilter.docker_tag = "1.0.3"

    ##############################

    call CGA_pipeline.CopyNumberReportQC_Task as CopyNumberReportQC_Task {
        input:
            tumorBam=tumorBam,
            tumorBamIdx=tumorBamIdx,
            normalBam=normalBam,
            normalBamIdx=normalBamIdx,
            readGroupBlackList=readGroupBlackList,
            tumorBam_size=tumorBam_size,
            normalBam_size=normalBam_size,
            diskGB_buffer=runtime_params["CopyNumberReportQC_Task.diskGB_buffer"],
            diskGB_boot=runtime_params["CopyNumberReportQC_Task.diskGB_boot"],
            preemptible=runtime_params["CopyNumberReportQC_Task.preemptible"],
            memoryGB=runtime_params["CopyNumberReportQC_Task.memoryGB"],
            cpu=runtime_params["CopyNumberReportQC_Task.cpu"]
    }

    # ContEst is a method for estimating the amount of cross-sample contamination in next generation sequencing data.
    # Using a Bayesian framework, contamination levels are estimated from array based genotypes and sequencing reads.
    call CGA_pipeline.ContEST_Task as ContEST_Task {
        input:
            tumorBam=tumorBam,
            tumorBamIdx=tumorBamIdx,
            normalBam=normalBam,
            normalBamIdx=normalBamIdx,
            pairName=pairName,
            refFasta=refFasta,
            refFastaIdx=refFastaIdx,
            refFastaDict=refFastaDict,
            refFasta_size=refFasta_size,
            targetIntervals=targetIntervals,
            tumorBam_size=tumorBam_size,
            normalBam_size=normalBam_size,
            diskGB_buffer=runtime_params["ContEST_Task.diskGB_buffer"],
            diskGB_boot=runtime_params["ContEST_Task.diskGB_boot"],
            preemptible=runtime_params["ContEST_Task.preemptible"],
            memoryGB=runtime_params["ContEST_Task.memoryGB"],
            cpu=runtime_params["ContEST_Task.cpu"]
    }

    # Program to check that all read groups within the set of BAM files appear to come from the same individual.
    call CGA_pipeline.CrossCheckLaneFingerprints_Task as CrossCheckLaneFingerprints_Task {
        input:
            tumorBam=tumorBam,
            normalBam=normalBam,
            tumorBamIdx=tumorBamIdx,
            normalBamIdx=normalBamIdx,
            pairName=pairName,
            GATK4_JAR=GATK4_JAR,
            gatk4_jar_size=gatk4_jar_size,
            tumorBam_size=tumorBam_size,
            normalBam_size=normalBam_size,
            validationStringencyLevel=runtime_params["CrossCheckLaneFingerprints_Task.validationStringencyLevel"],
            diskGB_buffer=runtime_params["CrossCheckLaneFingerprints_Task.diskGB_buffer"],
            diskGB_boot=runtime_params["CrossCheckLaneFingerprints_Task.diskGB_boot"],
            preemptible=runtime_params["CrossCheckLaneFingerprints_Task.preemptible"],
            memoryGB=runtime_params["CrossCheckLaneFingerprints_Task.memoryGB"],
            cpu=runtime_params["CrossCheckLaneFingerprints_Task.cpu"]
    }

    #Picard tasks (tumor and normal)
    ###################################
    # The task runs 3 tools:
    # ValidateSamFile, CollectMultipleMetrics and CollectWgsMetrics
    # ValidateSamFile makes sure the the given file is constructed correctly.
    # CollectMultipleMetrics collects multiple classes of metrics. This 'meta-metrics' tool runs one or more of the metrics collection modules at the same time
    # to cut down on the time spent reading in data from input files.
    # Available modules include CollectAlignmentSummaryMetrics, CollectInsertSizeMetrics, QualityScoreDistribution, MeanQualityByCycle,
    # CollectBaseDistributionByCycle, CollectGcBiasMetrics, RnaSeqMetrics, CollectSequencingArtifactMetrics, and CollectQualityYieldMetrics.
    # CollectWgsMetrics adds coverage statistics for WGS files, on top of CollectMultipleMetrics.

    # tumor
    if (forceComputePicardMetrics_tumor || !hasPicardMetrics_tumor) {
        call CGA_pipeline.PicardMultipleMetrics_Task as tumorMM_Task {
            input:
                bam=tumorBam,
                bamIndex=tumorBamIdx,
                sampleName=caseName,
                refFasta=refFasta,
                DB_SNP_VCF=DB_SNP_VCF,
                DB_SNP_VCF_IDX=DB_SNP_VCF_IDX,
                GATK4_JAR=GATK4_JAR,
                targetIntervals=targetIntervals,
                baitIntervals=baitIntervals,
                refFasta_size=refFasta_size,
                gatk4_jar_size=gatk4_jar_size,
                db_snp_vcf_size=db_snp_vcf_size,
                bam_size=tumorBam_size,
                diskGB_buffer=runtime_params["PicardMultipleMetrics_Task.diskGB_buffer"],
                diskGB_boot=runtime_params["PicardMultipleMetrics_Task.diskGB_boot"],
                preemptible=runtime_params["PicardMultipleMetrics_Task.preemptible"],
                memoryGB=runtime_params["PicardMultipleMetrics_Task.memoryGB"],
                cpu=runtime_params["PicardMultipleMetrics_Task.cpu"]
        }
    }

    #####################################
    #normal
    if (forceComputePicardMetrics_normal || !hasPicardMetrics_normal) {
        call CGA_pipeline.PicardMultipleMetrics_Task as normalMM_Task {
            input:
                bam=normalBam,
                bamIndex=normalBamIdx,
                sampleName=ctrlName,
                refFasta=refFasta,
                DB_SNP_VCF=DB_SNP_VCF,
                DB_SNP_VCF_IDX=DB_SNP_VCF_IDX,
                targetIntervals=targetIntervals,
                baitIntervals=baitIntervals,
                GATK4_JAR=GATK4_JAR,
                refFasta_size=refFasta_size,
                db_snp_vcf_size=db_snp_vcf_size,
                gatk4_jar_size=gatk4_jar_size,
                bam_size=normalBam_size,
                diskGB_buffer=runtime_params["PicardMultipleMetrics_Task.diskGB_buffer"],
                diskGB_boot=runtime_params["PicardMultipleMetrics_Task.diskGB_boot"],
                memoryGB=runtime_params["PicardMultipleMetrics_Task.memoryGB"],
                preemptible=runtime_params["PicardMultipleMetrics_Task.preemptible"],
                cpu=runtime_params["PicardMultipleMetrics_Task.cpu"]
        }
    }

    #####################################

    # PREPARE FOR SCATTER
    call CGA_pipeline.CallSomaticMutations_Prepare_Task as CallSomaticMutations_Prepare_Task {
        input:
            refFasta=refFasta,
            refFastaIdx=refFastaIdx,
            refFastaDict=refFastaDict,
            targetIntervals=targetIntervals,
            nWay=runtime_params["CallSomaticMutations_Prepare_Task.nWay"],
            diskGB_boot=runtime_params["CallSomaticMutations_Prepare_Task.diskGB_boot"],
            preemptible=runtime_params["CallSomaticMutations_Prepare_Task.preemptible"]
    }

    #SCATTER AND ANALYZE
    scatter (idx in CallSomaticMutations_Prepare_Task.scatterIndices) {
            # Identification of somatic point mutations in next generation sequencing data of cancer genomes.
            call CGA_pipeline.Mutect1_Task as Mutect1_Task {
                input:
                    tumorBam=tumorBam,
                    tumorBamIdx=tumorBamIdx,
                    normalBam=normalBam,
                    normalBamIdx=normalBamIdx,
                    pairName=pairName,
                    caseName=caseName,
                    ctrlName=ctrlName,
                    fracContam= if runtime_params["Mutect1_Task.contest_value"] != "" then runtime_params["Mutect1_Task.contest_value"] else ContEST_Task.fracContam,
                    mutectIntervals=CallSomaticMutations_Prepare_Task.interval_files[idx],
                    refFasta=refFasta,
                    refFastaIdx=refFastaIdx,
                    refFastaDict=refFastaDict,
                    DB_SNP_VCF=DB_SNP_VCF,
                    DB_SNP_VCF_IDX=DB_SNP_VCF_IDX,
                    cosmicVCF=cosmicVCF,
                    readGroupBlackList=readGroupBlackList,
                    MuTectNormalPanel=MuTectNormalPanel,
                    refFasta_size=refFasta_size,
                    db_snp_vcf_size=db_snp_vcf_size,
                    tumorBam_size=tumorBam_size,
                    normalBam_size=normalBam_size,
                    downsampleToCoverage=runtime_params["Mutect1_Task.downsampleToCoverage"],
                    diskGB_buffer=runtime_params["Mutect1_Task.diskGB_buffer"],
                    diskGB_boot=runtime_params["Mutect1_Task.diskGB_boot"],
                    preemptible=runtime_params["Mutect1_Task.preemptible"],
                    memoryGB=runtime_params["Mutect1_Task.memoryGB"],
                    cpu=runtime_params["Mutect1_Task.cpu"]
            }

            call CGA_pipeline.Mutect2_Task as Mutect2_Task {
                input:
                    tumorBam=tumorBam,
                    tumorBamIdx=tumorBamIdx,
                    normalBam=normalBam,
                    normalBamIdx=normalBamIdx,
                    pairName=pairName,
                    caseName=caseName,
                    ctrlName=ctrlName,
                    fracContam= if runtime_params["Mutect2_Task.contest_value"] != "" then runtime_params["Mutect2_Task.contest_value"] else ContEST_Task.fracContam,
                    mutectIntervals=CallSomaticMutations_Prepare_Task.interval_files[idx],
                    refFasta=refFasta,
                    refFastaIdx=refFastaIdx,
                    refFastaDict=refFastaDict,
                    readGroupBlackList=readGroupBlackList,
                    MuTectNormalPanel=MuTectNormalPanel,
                    GATK4_JAR=GATK4_JAR,
                    refFasta_size=refFasta_size,
                    tumorBam_size=tumorBam_size,
                    normalBam_size=normalBam_size,
                    gatk4_jar_size=gatk4_jar_size,
                    diskGB_buffer=runtime_params["Mutect2_Task.diskGB_buffer"],
                    diskGB_boot=runtime_params["Mutect2_Task.diskGB_boot"],
                    preemptible=runtime_params["Mutect2_Task.preemptible"],
                    memoryGB=runtime_params["Mutect2_Task.memoryGB"],
                    cpu=runtime_params["Mutect2_Task.cpu"]
            }
    }

    # MuTect is run in force-call mode to search for somatic variants at a set of specified loci for clinical relevance
    call CGA_pipeline.MutectFC_Task as MutectFC_Task {
        input:
            tumorBam=tumorBam,
            tumorBamIdx=tumorBamIdx,
            normalBam=normalBam,
            normalBamIdx=normalBamIdx,
            pairName=pairName,
            caseName=caseName,
            ctrlName=ctrlName,
            fracContam= if runtime_params["MutectFC_Task.contest_value"] != "" then runtime_params["MutectFC_Task.contest_value"] else ContEST_Task.fracContam,
            refFasta=refFasta,
            refFastaIdx=refFastaIdx,
            refFastaDict=refFastaDict,
            DB_SNP_VCF=DB_SNP_VCF,
            DB_SNP_VCF_IDX=DB_SNP_VCF_IDX,
            cosmicVCF=cosmicVCF,
            readGroupBlackList=readGroupBlackList,
            MuTectNormalPanel=MuTectNormalPanel,
            refFasta_size=refFasta_size,
            db_snp_vcf_size=db_snp_vcf_size,
            tumorBam_size=tumorBam_size,
            normalBam_size=normalBam_size,
            downsampleToCoverage=runtime_params["MutectFC_Task.downsampleToCoverage"],
            diskGB_buffer=runtime_params["MutectFC_Task.diskGB_buffer"],
            diskGB_boot=runtime_params["MutectFC_Task.diskGB_boot"],
            preemptible=runtime_params["MutectFC_Task.preemptible"],
            memoryGB=runtime_params["MutectFC_Task.memoryGB"],
            cpu=runtime_params["MutectFC_Task.cpu"]
    }

    # Strelka is an analysis package designed to detect somatic SNVs and small indels from the aligned sequencing reads
    # of matched tumor-normal samples.
    call CGA_pipeline.Strelka as Strelka {
        input:
            tumorBam=tumorBam,
            tumorBamIdx=tumorBamIdx,
            normalBam=normalBam,
            normalBamIdx=normalBamIdx,
            pairName=pairName,
            refFasta=refFasta,
            refFastaIdx=refFastaIdx,
            refFastaDict=refFastaDict,
            refFasta_size=refFasta_size,
            tumorBam_size=tumorBam_size,
            normalBam_size=normalBam_size,
            diskGB_buffer=runtime_params["Strelka.diskGB_buffer"],
            diskGB_boot=runtime_params["Strelka.diskGB_boot"],
            preemptible=runtime_params["Strelka.preemptible"],
            memoryGB=runtime_params["Strelka.memoryGB"],
            cpu=runtime_params["Strelka.cpu"]
    }

    # Gather outputs of MuTect1 and MuTect2
    call CGA_pipeline.Gather_Task as Gather_Task {
        input :
            mutect1_cs=Mutect1_Task.mutect1_cs,
            mutect2_cs=Mutect2_Task.mutect2_cs,
            pairName=pairName,
            diskGB_buffer=runtime_params["Gather_Task.diskGB_buffer"],
            diskGB_boot=runtime_params["Gather_Task.diskGB_boot"],
            preemptible=runtime_params["Gather_Task.preemptible"],
            memoryGB=runtime_params["Gather_Task.memoryGB"],
            cpu=runtime_params["Gather_Task.cpu"]
    }

    # Gather power and coverage wiggle files from MuTect1 and zip output
    call CGA_pipeline.GatherWIGFiles_Task as GatherWIGFiles_Task {
        input:
            pairName=pairName,
            mutect1_pw=Mutect1_Task.mutect1_pw,
            mutect1_cw=Mutect1_Task.mutect1_cw,
            diskGB_buffer=runtime_params["GatherWIGFiles_Task.diskGB_buffer"],
            diskGB_boot=runtime_params["GatherWIGFiles_Task.diskGB_boot"],
            memoryGB=runtime_params["GatherWIGFiles_Task.memoryGB"],
            preemptible=runtime_params["GatherWIGFiles_Task.preemptible"],
            cpu=runtime_params["GatherWIGFiles_Task.cpu"]
    }

    # GATK ACNV, an allelic copy-number variation method built on the Genome Analysis Toolkit.
    # ACNV is a tool for detecting somatic copy-number activity from whole exome and whole genome sequencing data by segmenting
    # the genome into regions of constant copy number and estimating copy ratio and minor-allele fraction in those regions.
    call CGA_pipeline.gatk_acnv_only as gatk_acnv_only {
        input:
            tumorBam=tumorBam,
            tumorBamIdx=tumorBamIdx,
            normalBam=normalBam,
            normalBamIdx=normalBamIdx,
            pairName=pairName,
            refFasta=refFasta,
            refFastaIdx=refFastaIdx,
            refFastaDict=refFastaDict,
            refFasta_size=refFasta_size,
            tumorBam_size=tumorBam_size,
            normalBam_size=normalBam_size,
            mutect1_call_stats=Gather_Task.MUTECT1_CS_SNV,
            diskGB_buffer=runtime_params["gatk_acnv_only.diskGB_buffer"],
            diskGB_boot=runtime_params["gatk_acnv_only.diskGB_boot"],
            preemptible=runtime_params["gatk_acnv_only.preemptible"],
            memoryGB=runtime_params["gatk_acnv_only.memoryGB"],
            cpu=runtime_params["gatk_acnv_only.cpu"]
    }

    # DeTiN estimates tumor in normal (TiN) based on tumor and matched normal sequencing data.
    # The estimate is based on both candidate SSNVs and aSCNAs.
    # DeTiN then applies the joint TiN estimate to reclassify SSNVs and InDels as somatic or germline.
    if (runDeTiN) {
        call CGA_pipeline.DeTiN_Task as DeTiN_Task {
            input :
                MUTECT1_CS=Gather_Task.MUTECT1_CS_SNV,
                MUTECT2_INDELS=Gather_Task.MUTECT2_VCF_INDELS,
                seg_file=gatk_acnv_only.alleliccapseg_tsv,
                tumor_hets=gatk_acnv_only.gatk_het_ad_tumor,
                normal_hets=gatk_acnv_only.gatk_het_ad_normal,
                pairName=pairName,
                release_version=runtime_params["DeTiN_Task.release_version"],
                Mutation_prior=runtime_params["DeTiN_Task.Mutation_prior"],
                TiN_prior=runtime_params["DeTiN_Task.TiN_prior"],
                diskGB_buffer=runtime_params["DeTiN_Task.diskGB_buffer"],
                diskGB_boot=runtime_params["DeTiN_Task.diskGB_boot"],
                preemptible=runtime_params["DeTiN_Task.preemptible"],
                memoryGB=runtime_params["DeTiN_Task.memoryGB"],
                cpu=runtime_params["DeTiN_Task.cpu"]
        }
    }

    # VEP determines the effect of your variants (SNPs, insertions, deletions, CNVs or structural variants) on genes, transcripts, and protein sequence, as well as regulatory regions. Simply input the coordinates of your variants and the nucleotide changes to find out the:
    call CGA_pipeline.VEP_Task as VEP_Task {
        input:
            MUTECT1_CS=Gather_Task.MUTECT1_CS_SNV,
            MUTECT2_VCF=Gather_Task.MUTECT2_VCF_INDELS,
            STRELKA_VCF=Strelka.Strelka_passed_indels,
            pairName=pairName,
            caseName=caseName,
            ctrlName=ctrlName,
            diskGB_buffer=runtime_params["VEP_Task.diskGB_buffer"],
            diskGB_boot=runtime_params["VEP_Task.diskGB_boot"],
            preemptible=runtime_params["VEP_Task.preemptible"],
            memoryGB=runtime_params["VEP_Task.memoryGB"],
            cpu=runtime_params["VEP_Task.cpu"]
    }

    # Oncotator is a tool for annotating human genomic point mutations and indels with data relevant to cancer researchers.
    call CGA_pipeline.Oncotate_Task as Oncotate_Task {
        input :
            MUTECT1_CS=VEP_Task.MUTECT1_VEP_annotated_filtered_vcf,
            MUTECT2_INDELS=VEP_Task.MUTECT2_VEP_annotated_filtered_vcf,
            STRELKA_INDELS=VEP_Task.STRELKA_VEP_annotated_filtered_vcf,
            pairName=pairName,
            caseName=caseName,
            ctrlName=ctrlName,
            diskGB_buffer=runtime_params["Oncotate_Task.diskGB_buffer"],
            diskGB_boot=runtime_params["Oncotate_Task.diskGB_boot"],
            preemptible=runtime_params["Oncotate_Task.preemptible"],
            memoryGB=runtime_params["Oncotate_Task.memoryGB"],
            cpu=runtime_params["Oncotate_Task.cpu"]
    }

    # Detects and screens out OxoG artifacts from a set of SNV calls.
    # Oxidation of guanine to 8-oxoguanine is one of the most common pre-adapter artifacts associated with genomic library preparation,
    # arising from a combination of heat, shearing, and metal contaminates in a sample).
    # The 8-oxoguanine base can pair with either cytosine or adenine, ultimately leading to G→T transversion mutations during PCR amplification.
    # CC. -> CA.
    # .GG -> .TG <= DNA F1R2 (Context - ".GG", REF Allele - "G", ALT Allele - "T")
    call CGA_pipeline.OrientationBias_filter_Task as oxoGOBF {
        input:
            stub="oxog",
            tumorBam=select_first([tumorMM_Task.Bam,tumorBam]),
            tumorBamIdx=select_first([tumorMM_Task.Bai,tumorBamIdx]),
            pairName=pairName,
            detailMetrics=tumorMM_Task.pre_adapter_detail_metrics,
            MAF=Oncotate_Task.WXS_Mutation_M1_SNV_M2_INDEL_Strelka_INDEL_annotated_maf,
            GATK4_JAR=GATK4_JAR,
            refFasta=refFasta,
            refFasta_size=refFasta_size,
            tumorBam_size=tumorBam_size,
            gatk4_jar_size=gatk4_jar_size,
            diskGB_buffer=runtime_params["OrientationBias_filter_Task.diskGB_buffer"],
            diskGB_boot=runtime_params["OrientationBias_filter_Task.diskGB_boot"],
            preemptible=runtime_params["OrientationBias_filter_Task.preemptible"],
            memoryGB=runtime_params["OrientationBias_filter_Task.memoryGB"],
            cpu=runtime_params["OrientationBias_filter_Task.cpu"]
    }

    # Detects and screens out FFPE artifacts from a set of SNV calls.
    # FFPE introduces multiple types of DNA damage including deamination, which converts cytosine to uracil and leads to downstream mispairing
    # in PCR: C>T/G>A. Because deamination occurs prior to ligation of palindromic Illumina adapters, likely deamination artifacts will have
    # a read orientation bias. The FFPE Filter Task uses this read orientation to identify artifacts and calculate a Phred scaled Q-score for FFPE artifacts.
    # .CG -> .TG <= DNA F1R2 (Context - ".CG", REF Allele - "C", ALT Allele - "T")
    # CG. -> CA.
    call CGA_pipeline.OrientationBias_filter_Task as ffpeOBF {
        input:
            stub="ffpe",
            tumorBam=select_first([tumorMM_Task.Bam,tumorBam]),
            tumorBamIdx=select_first([tumorMM_Task.Bai,tumorBamIdx]),
            pairName=pairName,
            detailMetrics=tumorMM_Task.pre_adapter_detail_metrics,
            MAF=Oncotate_Task.WXS_Mutation_M1_SNV_M2_INDEL_Strelka_INDEL_annotated_maf,
            GATK4_JAR=GATK4_JAR,
            refFasta=refFasta,
            refFasta_size=refFasta_size,
            tumorBam_size=tumorBam_size,
            gatk4_jar_size=gatk4_jar_size,
            diskGB_buffer=runtime_params["OrientationBias_filter_Task.diskGB_buffer"],
            diskGB_boot=runtime_params["OrientationBias_filter_Task.diskGB_boot"],
            preemptible=runtime_params["OrientationBias_filter_Task.preemptible"],
            memoryGB=runtime_params["OrientationBias_filter_Task.memoryGB"],
            cpu=runtime_params["OrientationBias_filter_Task.cpu"]
    }

    # MAFPoNFilter uses a likelihood model to compare somatic mutations against a Panel of Normals (PoN)
    # in order to screen out somatic mutations. The PoN represents sequencing conditions in the case sample,
    # including germline variants and technical artifacts. Refer to the Panel of Normals section for more information.
    scatter (pon_object in PONs_data) {
        call MAFPonFilter{
            input:
                MAFFile=Oncotate_Task.WXS_Mutation_M1_SNV_M2_INDEL_Strelka_INDEL_annotated_maf,
                pairName=pairName,
                cytoBandFile=cytoBandFile,
                PONFile=pon_object.pon_url,
                stub=pon_object.pon_name,
                TOTNStr=runtime_params["MAFPonFilter.TOTNStr"],
                NMIN=runtime_params["MAFPonFilter.NMIN"],
                THRESH=pon_object.pon_threshold,
                CODING_ONLY=runtime_params["MAFPonFilter.CODING_ONLY"],
                MIN_ALT_COUNT=runtime_params["MAFPonFilter.MIN_ALT_COUNT"],
                public_release=runtime_params["MAFPonFilter.public_release"],
                diskGB_buffer=runtime_params["MAFPonFilter.diskGB_buffer"],
                diskGB_boot=runtime_params["MAFPonFilter.diskGB_boot"],
                preemptible=runtime_params["MAFPonFilter.preemptible"],
                memoryGB=runtime_params["MAFPonFilter.memoryGB"],
                cpu=runtime_params["MAFPonFilter.cpu"]
        }
    }

    call CGA_pipeline.blat as blat {
        input:
            tumorBam=tumorBam,
            tumorBamIdx=tumorBamIdx,
            MAF=Oncotate_Task.WXS_Mutation_M1_SNV_M2_INDEL_Strelka_INDEL_annotated_maf,
            pairName=pairName,
            tumorBam_size=tumorBam_size,
            diskGB_buffer=runtime_params["blat.diskGB_buffer"],
            diskGB_boot=runtime_params["blat.diskGB_boot"],
            preemptible=runtime_params["blat.preemptible"],
            memoryGB=runtime_params["blat.memoryGB"],
            cpu=runtime_params["blat.cpu"]
    }

    call CGA_pipeline.merge_mafs_task as merge_mafs_task {
        input:
            oxoGOBF_maf=oxoGOBF.WXS_Mutation_OBF_unfiltered_maf_with_annotations,
            ffpeOBF_maf=ffpeOBF.WXS_Mutation_OBF_unfiltered_maf_with_annotations,
            pon_filtered_mafs=MAFPonFilter.allMaf_with_annotations,
            blat_maf=blat.blat_all_maf_with_annotations,
            pairName=pairName,
            diskGB_buffer=runtime_params["merge_mafs_task.diskGB_buffer"],
            diskGB_boot=runtime_params["merge_mafs_task.diskGB_boot"],
            preemptible=runtime_params["merge_mafs_task.preemptible"],
            memoryGB=runtime_params["merge_mafs_task.memoryGB"],
            cpu=runtime_params["merge_mafs_task.cpu"]
    }

    call CGA_pipeline.mutation_validator as mutation_validator {
        input:
            pairName=pairName,
            MAF=merge_mafs_task.merged_intersection_maf,
            tumorBam=tumorBam,
            normalBam=normalBam,
            tumorBamIdx=tumorBamIdx,
            normalBamIdx=normalBamIdx,
            tumorBam_size=tumorBam_size,
            normalBam_size=normalBam_size,
            maf_type=runtime_params["mutation_validator.maf_type"],
            diskGB_buffer=runtime_params["mutation_validator.diskGB_buffer"],
            diskGB_boot=runtime_params["mutation_validator.diskGB_boot"],
            preemptible=runtime_params["mutation_validator.preemptible"],
            memoryGB=runtime_params["mutation_validator.memoryGB"],
            cpu=runtime_params["mutation_validator.cpu"]
    }

    if (merge_mafs_task.found_snv) {
    # Estimate purity/ploidy, and from that compute absolute copy-number and mutation multiplicities.
        call CGA_pipeline.absolute as absolute {
            input:
                maf=mutation_validator.validated_maf,
                seg_file=gatk_acnv_only.alleliccapseg_tsv,
                skew=gatk_acnv_only.alleliccapseg_skew,
                pairName=pairName,
                diskGB_buffer=runtime_params["absolute.diskGB_buffer"],
                diskGB_boot=runtime_params["absolute.diskGB_boot"],
                preemptible=runtime_params["absolute.preemptible"],
                memoryGB=runtime_params["absolute.memoryGB"],
                cpu=runtime_params["absolute.cpu"]

        }


        call CGA_pipeline.lego_plotter_task as lego_plotter_task {
            input:
                maf=mutation_validator.validated_maf,
                pairName=pairName,
                covString=runtime_params["lego_plotter_task.covString"],
                diskGB_buffer=runtime_params["lego_plotter_task.diskGB_buffer"],
                diskGB_boot=runtime_params["lego_plotter_task.diskGB_boot"],
                preemptible=runtime_params["lego_plotter_task.preemptible"],
                memoryGB=runtime_params["lego_plotter_task.memoryGB"],
                cpu=runtime_params["lego_plotter_task.cpu"]
        }

        ############ Extra filters beyond CGA ####

        #### CommonVariantFilter #####
        call CommonVariantFilter.commonfilterTask as commonfilterTask {
            input: sampleId=CommonVariantFilter.sampleId,
                maf=mutation_validator.validated_maf, #maf,
                min_exac_ac=CommonVariantFilter.min_exac_ac,
                min_filter_depth=CommonVariantFilter.min_filter_depth,
                disable_whitelist=CommonVariantFilter.disable_whitelist,
                filter_noncoding=CommonVariantFilter.filter_noncoding,
                RAM=CommonVariantFilter.RAM,
                SSD=CommonVariantFilter.SSD,
                preemptible=CommonVariantFilter.preemptible,
                docker_tag=CommonVariantFilter.docker_tag
        }

    }

    output {
        #### CGA_pipeline ######


        ####### QC Tasks Outputs #######
        # Copy Number QC Report files
        File tumor_bam_lane_list=CopyNumberReportQC_Task.tumorBamLaneList
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
