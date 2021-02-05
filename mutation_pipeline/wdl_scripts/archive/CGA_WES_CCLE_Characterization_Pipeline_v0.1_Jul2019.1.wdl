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

workflow CGA_Production_Analysis_Workflow {

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
   
##############################

    call CopyNumberReportQC_Task {
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
    call ContEST_Task {
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
    call CrossCheckLaneFingerprints_Task {
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
        call PicardMultipleMetrics_Task as tumorMM_Task {
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
                validationStringencyLevel=runtime_params["PicardMultipleMetrics_Task.validationStringencyLevel"],
                run_clean_sam=runtime_params["PicardMultipleMetrics_Task.run_clean_sam"],
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
        call PicardMultipleMetrics_Task as normalMM_Task {
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
                validationStringencyLevel=runtime_params["PicardMultipleMetrics_Task.validationStringencyLevel"],
                run_clean_sam=runtime_params["PicardMultipleMetrics_Task.run_clean_sam"],
                diskGB_buffer=runtime_params["PicardMultipleMetrics_Task.diskGB_buffer"],
                diskGB_boot=runtime_params["PicardMultipleMetrics_Task.diskGB_boot"],
                memoryGB=runtime_params["PicardMultipleMetrics_Task.memoryGB"], 
                preemptible=runtime_params["PicardMultipleMetrics_Task.preemptible"],
                cpu=runtime_params["PicardMultipleMetrics_Task.cpu"]
        }
    }

    #####################################

    # PREPARE FOR SCATTER
    call CallSomaticMutations_Prepare_Task {
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
            call Mutect1_Task {
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

            call Mutect2_Task {
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
    call MutectFC_Task {
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
    call Strelka {
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
    call Gather_Task {
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
    call GatherWIGFiles_Task {
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
    call gatk_acnv_only {
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
        call DeTiN_Task {
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
    call VEP_Task {
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
    call Oncotate_Task {
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
    call OrientationBias_filter_Task as oxoGOBF {
        input:
            stub="oxog",
            tumorBam=tumorBam,
            tumorBamIdx=tumorBamIdx,
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
    call OrientationBias_filter_Task as ffpeOBF {
        input:
            stub="ffpe",
            tumorBam=tumorBam,
            tumorBamIdx=tumorBamIdx,
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

    call blat {
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

    call merge_mafs_task {
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

    call mutation_validator {
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
        call absolute {
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


        call lego_plotter_task {
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
    }

    output {
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
        File? normal_bam_bam_validation=normalMM_Task.bam_validation
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
        File? tumor_bam_bam_validation=tumorMM_Task.bam_validation
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
    }
}


# TASKS DEFINITION

task CopyNumberReportQC_Task {

    # TASK INPUT PARAMS
    File tumorBam
    File tumorBamIdx
    File normalBam
    File normalBamIdx
    File regionFile
    File readGroupBlackList
    File captureNormalsDBRCLZip

    # FILE SIZE
    Int tumorBam_size
    Int normalBam_size
    
    # RUNTIME INPUT PARAMS
    String preemptible
    String diskGB_boot
    String diskGB_buffer
    String memoryGB
    String cpu

    # DEFAULT VALUES
    String default_cpu = "1"
    String default_memoryGB = "7"
    String default_preemptible = "1"
    String default_diskGB_boot = "15"
    String default_diskGB_buffer = "20"

    # COMPUTE MEMORY SIZE
    Int machine_memoryGB = if memoryGB != "" then memoryGB else default_memoryGB
    Int command_memoryGB = machine_memoryGB - 1
    
    # COMPUTE DISK SIZE
    Int machine_diskGB_buffer = if diskGB_buffer != "" then diskGB_buffer else default_diskGB_buffer
    Int diskGB = ceil(tumorBam_size + normalBam_size + size(regionFile, "G") + size(readGroupBlackList, "G") + size(captureNormalsDBRCLZip, "G")*2
                    + machine_diskGB_buffer)

    parameter_meta {
        tumorBam : "sample tumor  BAM file"
        tumorBamIdx : "sample tumor BAI file (BAM indexed)"
        normalBam : "sample normal BAM file"
        normalBamIdx : "sample normal BAI file (BAM indexed)"
        regionFile : ""
        readGroupBlackList : ""
        captureNormalsDBRCLZip : ""
    }

    command <<<

        # e - exit when a command fails
        # u - exit when script tries to use undeclared variables
        # x - trace what gets executed
        # o - exit status of the last command that threw a non-zero exit code is returned
        set -euxo pipefail

        # Copy Number QC for Capture
        #Make lane lists for tumor and normal
        #MakeLaneList_12        
        java "-Xmx${command_memoryGB}g" -jar /usr/local/bin/MakeLaneList.jar ${tumorBam}  case.lanelist.txt ;
        java "-Xmx${command_memoryGB}g" -jar /usr/local/bin/MakeLaneList.jar ${normalBam} control.lanelist.txt ;

        # mMke region coverages per lane for tumor
        #RegionCovPerLane_18
        for CHROM in `seq 1 24` ; do
            echo "Working on $CHROM for ${tumorBam}"
            SCATTER_DIR="scatter.$CHROM" ;
            mkdir -v $SCATTER_DIR
            OUTPATH=`echo -ne "$SCATTER_DIR/chr$CHROM.case.rcl"` ;
            echo "OUTPATH is $OUTPATH"
            echo "chrom is $CHROM"
            java -jar /usr/local/bin/RegionCovPerLane.jar ${tumorBam} ${regionFile} $OUTPATH $CHROM
        done ;
        wait ;
        /usr/local/bin/rcl.gather.sh  case.rcl ;
        rm -rf scatter.* ;

        #make region coverages per lane for control
        #RegionCovPerLane_18
        for CHROM in `seq 1 24` ; do
            echo "Working on $CHROM for ${normalBam}"
            SCATTER_DIR="scatter.$CHROM" ;
            mkdir -v $SCATTER_DIR
            OUTPATH=`echo -ne "$SCATTER_DIR/chr$CHROM.control.rcl"` ;
            echo "OUTPATH is $OUTPATH"
            echo "chrom is $CHROM"
            java -jar /usr/local/bin/RegionCovPerLane.jar ${normalBam} ${regionFile} $OUTPATH $CHROM
        done ;
        wait ;
        /usr/local/bin/rcl.gather.sh  control.rcl ;
        rm -rf scatter.* ;

        #run the mat-lab based report
        #CopyNumberQCReport_27
        cp -vfr /CopyNumberQCReport_27/unzip/* .

        #Make a file of files paths
        #This command aims to list the zip files contents, filtering out all but the file paths and then write the paths to a file
python <<CODE
from zipfile import ZipFile

with open('capture_normals_db_wdl', 'w') as writer:
    with ZipFile("${captureNormalsDBRCLZip}", 'r') as f:
        names = f.namelist()
        for name in names:
            writer.write(name + '\n')
CODE

        unzip ${captureNormalsDBRCLZip}
        ./run_fh_CopyNumberQCReport.sh /opt/MATLAB/MATLAB_Compiler_Runtime/v710/ PairCopyNumQCReport \
        case.rcl control.rcl case.lanelist.txt control.lanelist.txt \
        ${readGroupBlackList} ${regionFile} capture_normals_db_wdl NA NA
        #zip up the output
        zip CopyNumQCout.zip report.html num.mixups.txt PairCopyNumQCReport* *.png

    >>>

    runtime {
        docker         : "gcr.io/broad-getzlab-workflows/cga_production_pipeline:v0.2.ccle"
        bootDiskSizeGb : if diskGB_boot != "" then diskGB_boot else default_diskGB_boot
        preemptible    : if preemptible != "" then preemptible else default_preemptible
        cpu            : if cpu != "" then cpu else default_cpu
        disks          : "local-disk ${diskGB} HDD"
        memory         : machine_memoryGB + "GB"
    }

    output {
        #lane lists
        File tumorBamLaneList="case.lanelist.txt"
        File normalBamLaneList="control.lanelist.txt"
        File tumorRCL="case.rcl"
        File normalRCL="control.rcl"
        #CopyNumQC
        File CopyNumQCOutZip="CopyNumQCout.zip"
        File CopyNumQCReport="report.html"
        File CopyNumQCReportPNG="PairCopyNumQCReport_CopyNumberQC.png"
        File CopyNumQCMixUps="num.mixups.txt"
    }
}


task PicardMultipleMetrics_Task {

    # TASK INPUT PARAMS
    File bam
    File bamIndex
    String sampleName
    File targetIntervals
    File baitIntervals
    File refFasta
    File DB_SNP_VCF
    File DB_SNP_VCF_IDX
    File GATK4_JAR

    String validationStringencyLevel
    String run_clean_sam

    # FILE SIZE
    Int bam_size
    Int refFasta_size
    Int gatk4_jar_size
    Int db_snp_vcf_size

    # RUNTIME INPUT PARAMS
    String preemptible
    String diskGB_boot
    String diskGB_buffer
    String memoryGB
    String cpu

    # DEFAULT VALUES
    String default_cpu = "1"
    String default_memoryGB = "10"
    String default_preemptible = "1"
    String default_diskGB_boot = "15"
    String default_diskGB_buffer = "20"
    String default_stringencyLevel = "LENIENT"
    String default_run_clean_sam = false

    # COMPUTE MEMORY SIZE
    Int machine_memoryGB = if memoryGB != "" then memoryGB else default_memoryGB
    Int command_memoryGB = machine_memoryGB - 1
    
    # COMPUTE DISK SIZE
    Int machine_diskGB_buffer = if diskGB_buffer != "" then diskGB_buffer else default_diskGB_buffer
    Int diskGB = ceil(bam_size + refFasta_size + gatk4_jar_size + db_snp_vcf_size + machine_diskGB_buffer)

    String stringencyLevel = if validationStringencyLevel != "" then validationStringencyLevel else default_stringencyLevel
    String clean_sam_flag = if run_clean_sam != "" then run_clean_sam else default_run_clean_sam

    parameter_meta {
        bam : "sample (normal or tumor) BAM file"
        bamIndex : "sample (normal or tumor) BAI file (BAM indexed)"
        sampleName : "sample (normal or tumor) name, prefix for output"
        refFasta : "FASTA file for the appropriate genome build (Reference sequence file)"
        DB_SNP_VCF : "VCF format dbSNP file, used to exclude regions around known polymorphisms from analysis by some PROGRAMs"
        DB_SNP_VCF_IDX : "dbSNP indexed file"
    }

    command <<<

        set -euxo pipefail

        /usr/local/jre1.8.0_73/bin/java "-Xmx${command_memoryGB}g" -jar ${GATK4_JAR} ValidateSamFile \
        --INPUT ${bam} \
        --OUTPUT "${sampleName}.bam_validation" \
        --MODE VERBOSE \
        --IGNORE_WARNINGS true \
        --REFERENCE_SEQUENCE ${refFasta} \
        --VALIDATION_STRINGENCY ${stringencyLevel}

        if [ "${clean_sam_flag}" = true ] ;
        then
            # Run bam through CleanSam to set MAPQ of unmapped reads to zero
            /usr/local/jre1.8.0_73/bin/java "-Xmx${command_memoryGB}g" -jar ${GATK4_JAR} CleanSam \
            --INPUT ${bam} \
            --OUTPUT ${sampleName}.unmapped_reads_cleaned.bam
        fi

        /usr/local/jre1.8.0_73/bin/java "-Xmx${command_memoryGB}g" -jar ${GATK4_JAR} CollectMultipleMetrics \
        --INPUT ${bam} \
        --OUTPUT ${sampleName}.multiple_metrics \
        --REFERENCE_SEQUENCE ${refFasta} \
        --DB_SNP ${DB_SNP_VCF} \
        --PROGRAM CollectAlignmentSummaryMetrics \
        --PROGRAM CollectInsertSizeMetrics \
        --PROGRAM QualityScoreDistribution \
        --PROGRAM MeanQualityByCycle \
        --PROGRAM CollectBaseDistributionByCycle \
        --PROGRAM CollectSequencingArtifactMetrics \
        --PROGRAM CollectQualityYieldMetrics \
        --PROGRAM CollectGcBiasMetrics
        
        #Extract OxoG metrics from generalized artifacts metrics. 
        # This tool extracts 8-oxoguanine (OxoG) artifact metrics from the output of CollectSequencingArtifactsMetrics 
        # (a tool that provides detailed information on a variety of
        # artifacts found in sequencing libraries) and converts them to the CollectOxoGMetrics tool's output format. This
        # conveniently eliminates the need to run CollectOxoGMetrics if we already ran CollectSequencingArtifactsMetrics in our
        # pipeline.
        /usr/local/jre1.8.0_73/bin/java -jar ${GATK4_JAR} ConvertSequencingArtifactToOxoG \
        --INPUT_BASE "${sampleName}.multiple_metrics" \
        --OUTPUT_BASE "${sampleName}.multiple_metrics.converted" \
        --REFERENCE_SEQUENCE ${refFasta} \
        --VALIDATION_STRINGENCY ${stringencyLevel}

        #zip up reports for QC Nozzle report
        zip picard_multiple_metrics.zip ${sampleName}.multiple_metrics.*

        # Collect WES HS metrics
        /usr/local/jre1.8.0_73/bin/java "-Xmx${command_memoryGB}g" -jar ${GATK4_JAR} CollectHsMetrics \
        --INPUT ${bam} \
        --BAIT_INTERVALS ${targetIntervals} \
        --TARGET_INTERVALS ${baitIntervals} \
        --OUTPUT "${sampleName}.HSMetrics.txt" \
        --VALIDATION_STRINGENCY ${stringencyLevel}

    >>>

    runtime {
        docker: "gcr.io/broad-getzlab-workflows/cga_production_pipeline:v0.2.ccle"
        bootDiskSizeGb: if diskGB_boot != "" then diskGB_boot else default_diskGB_boot
        preemptible: if preemptible != "" then preemptible else default_preemptible
        cpu: if cpu != "" then cpu else default_cpu
        disks: "local-disk ${diskGB} HDD"
        memory: machine_memoryGB + "GB"
    }

    output {
        File bam_validation="${sampleName}.bam_validation"
        File metricsReportsZip="picard_multiple_metrics.zip"
        File alignment_summary_metrics="${sampleName}.multiple_metrics.alignment_summary_metrics"
        File bait_bias_detail_metrics="${sampleName}.multiple_metrics.bait_bias_detail_metrics"
        File bait_bias_summary_metrics="${sampleName}.multiple_metrics.bait_bias_summary_metrics"
        File base_distribution_by_cycle="${sampleName}.multiple_metrics.base_distribution_by_cycle.pdf"
        File base_distribution_by_cycle_metrics="${sampleName}.multiple_metrics.base_distribution_by_cycle_metrics"
        File gc_bias_detail_metrics="${sampleName}.multiple_metrics.gc_bias.detail_metrics"
        File gc_bias="${sampleName}.multiple_metrics.gc_bias.pdf"
        File gc_bias_summary_metrics="${sampleName}.multiple_metrics.gc_bias.summary_metrics"
        File insert_size_histogram="${sampleName}.multiple_metrics.insert_size_histogram.pdf"
        File insert_size_metrics="${sampleName}.multiple_metrics.insert_size_metrics"
        #File oxog_metrics="${sampleName}.multiple_metrics.oxog_metrics"
        File pre_adapter_detail_metrics="${sampleName}.multiple_metrics.pre_adapter_detail_metrics"
        File pre_adapter_summary_metrics="${sampleName}.multiple_metrics.pre_adapter_summary_metrics"
        File quality_by_cycle="${sampleName}.multiple_metrics.quality_by_cycle.pdf"
        File quality_by_cycle_metrics="${sampleName}.multiple_metrics.quality_by_cycle_metrics"
        File quality_distribution="${sampleName}.multiple_metrics.quality_distribution.pdf"
        File quality_distribution_metrics="${sampleName}.multiple_metrics.quality_distribution_metrics"
        File quality_yield_metrics="${sampleName}.multiple_metrics.quality_yield_metrics"
        File converted_oxog_metrics="${sampleName}.multiple_metrics.converted.oxog_metrics"
        File hsMetrics="${sampleName}.HSMetrics.txt"
    }
}


task ContEST_Task {
        
    # TASK INPUT PARAMS
    File tumorBam
    File tumorBamIdx
    File normalBam
    File normalBamIdx
    File refFasta
    File refFastaIdx
    File refFastaDict
    File targetIntervals
    
    File SNP6Bed
    File HapMapVCF
    String pairName

    # FILE SIZE
    Int tumorBam_size
    Int normalBam_size
    Int refFasta_size

    # RUNTIME INPUT PARAMS
    String preemptible
    String diskGB_boot
    String diskGB_buffer
    String memoryGB
    String cpu

    # DEFAULT VALUES
    String default_cpu = "1"
    String default_memoryGB = "10"
    String default_preemptible = "1"
    String default_diskGB_boot = "15"
    String default_diskGB_buffer = "20"

    # COMPUTE MEMORY SIZE
    Int machine_memoryGB = if memoryGB != "" then memoryGB else default_memoryGB
    Int command_memoryGB = machine_memoryGB - 1
    
    # COMPUTE DISK SIZE
    Int machine_diskGB_buffer = if diskGB_buffer != "" then diskGB_buffer else default_diskGB_buffer
    Int diskGB = ceil(tumorBam_size + normalBam_size + refFasta_size
                + size(targetIntervals, "G") + size(SNP6Bed, "G") + size(HapMapVCF, "G")
                + machine_diskGB_buffer)

    parameter_meta {
        tumorBam : "sample tumor BAM file"
        tumorBamIdx : "sample tumor BAI file (indexed BAM file)"
        normalBam : "sample normal BAM file"
        normalBamIdx : "sample normal BAI file (indexed BAM file)"
        pairName : "sample name, prefix for output"
        refFasta : "the FASTA file for the appropriate genome build (Reference sequence file)"
        refFastaIdx : "FASTA file index for the reference genome"
        refFastaDict : "FASTA file dictionary for the reference genome"
        exomeIntervals : ""
        SNP6Bed : ""
        HapMapVCF : "the population allele frequencies for each SNP in HapMap"
    }

    command <<<

        set -euxo pipefail
         
        java "-Xmx${command_memoryGB}g" -Djava.io.tmpdir=/tmp -jar /usr/local/bin/GenomeAnalysisTK.jar \
        -T ContEst \
        -I:eval ${tumorBam} \
        -I:genotype ${normalBam} \
        -L ${targetIntervals} \
        -L ${SNP6Bed} \
        -isr INTERSECTION \
        -R ${refFasta} \
        -l INFO \
        -pf ${HapMapVCF} \
        -o contamination.af.txt \
        --trim_fraction 0.03 \
        --beta_threshold 0.05 \
        -br contamination.base_report.txt \
        -mbc 100 \
        --min_genotype_depth 30 \
        --min_genotype_ratio 0.8

        python /usr/local/bin/extract_contamination.py contamination.af.txt fraction_contamination.txt \
        contamination_validation.array_free.txt ${pairName}

        #Contamination validation/consensus
        python /usr/local/populateConsensusContamination_v26/contaminationConsensus.py \
        --pass_snp_qc false \
        --output contest_validation.output.tsv \
        --annotation contamination_percentage_consensus_capture \
        --array contamination_validation.array_free.txt \
        --noarray contamination_validation.array_free.txt

    >>>

    runtime {
        docker: "gcr.io/broad-getzlab-workflows/cga_production_pipeline:v0.2.ccle"
        bootDiskSizeGb: if diskGB_boot != "" then diskGB_boot else default_diskGB_boot
        preemptible: if preemptible != "" then preemptible else default_preemptible
        cpu: if cpu != "" then cpu else default_cpu
        disks: "local-disk ${diskGB} HDD"
        memory: machine_memoryGB + "GB"
    }

    output {
        File contamDataFile="contamination_validation.array_free.txt"
        File contestAFFile="contamination.af.txt"
        File contestBaseReport="contamination.base_report.txt"
        File validationOutput="contest_validation.output.tsv"
        Float fracContam=read_float("fraction_contamination.txt")
    }
}


task CrossCheckLaneFingerprints_Task {

    # TASK INPUT PARAMS
    File tumorBam
    File normalBam
    File tumorBamIdx
    File normalBamIdx
    String pairName
    File HaplotypeDBForCrossCheck
    File GATK4_JAR

    String validationStringencyLevel

    # FILE SIZE
    Int tumorBam_size
    Int normalBam_size
    Int gatk4_jar_size

    # RUNTIME INPUT PARAMS
    String preemptible
    String diskGB_boot
    String diskGB_buffer
    String memoryGB
    String cpu

    # DEFAULT VALUES
    String default_cpu = "1"
    String default_memoryGB = "3"
    String default_preemptible = "1"
    String default_diskGB_boot = "15"
    String default_diskGB_buffer = "20"
    String default_stringencyLevel = "LENIENT"

    # COMPUTE MEMORY SIZE
    Int machine_memoryGB = if memoryGB != "" then memoryGB else default_memoryGB
    Int command_memoryGB = machine_memoryGB - 1
    
    # COMPUTE DISK SIZE
    Int machine_diskGB_buffer = if diskGB_buffer != "" then diskGB_buffer else default_diskGB_buffer
    Int diskGB = ceil(tumorBam_size + normalBam_size + gatk4_jar_size + size(HaplotypeDBForCrossCheck, "G") 
                    + machine_diskGB_buffer)

    String stringencyLevel = if validationStringencyLevel != "" then validationStringencyLevel else default_stringencyLevel

    parameter_meta {
        tumorBam : "sample tumor BAM file"
        tumorBamIdx : "sample tumor BAI file (indexed BAM file)"
        normalBam : "sample normal BAM file"
        normalBamIdx : "sample normal BAI file (indexed BAM file)"
        pairName : "a string for the name of the pair under analysis used for naming output files"
        HaplotypeDBForCrossCheck : ""
        validationStringencyLevel : ""
    }

    command <<<

        set -euxo pipefail

        #drop from haplotypeDB seq entries which aren't in BAM if there are any found
        PREPPED_HAPLOTYPE_DB=PreppedHaplotypeDB.txt
        /usr/local/bin/filter_not_in_bam_dict.pl ${normalBam} ${HaplotypeDBForCrossCheck} $PREPPED_HAPLOTYPE_DB

        #CrosscheckLaneFingerprints[version=9]
        mkdir -v tmp
        /usr/local/jre1.8.0_73/bin/java "-Xmx${command_memoryGB}g" -jar ${GATK4_JAR} CrosscheckFingerprints \
        -I ${tumorBam} \
        -I ${normalBam} \
        -H $PREPPED_HAPLOTYPE_DB \
        --TMP_DIR `pwd`/tmp \
        --QUIET false \
        --EXIT_CODE_WHEN_MISMATCH 0 \
        --OUTPUT crosscheck.metrics \
        --VALIDATION_STRINGENCY ${stringencyLevel} 

        #Produce crosscheck.stats.txt file for making the html report
        grep -v "#" crosscheck.metrics | sed 1d > crosscheck.metrics.stripped
        
        python /usr/local/bin/crosscheck_report.py crosscheck.metrics.stripped

    >>>

    runtime {
        docker: "gcr.io/broad-getzlab-workflows/cga_production_pipeline:v0.2.ccle"
        bootDiskSizeGb: if diskGB_boot != "" then diskGB_boot else default_diskGB_boot
        preemptible: if preemptible != "" then preemptible else default_preemptible
        cpu: if cpu != "" then cpu else default_cpu
        disks: "local-disk ${diskGB} HDD"
        memory: machine_memoryGB + "GB"
    }

    output {
        File crossCheckMetrics="crosscheck.metrics"
        File crossCheckReport="report.html"
        Float crossCheckMinLODValue=read_float("crosscheck_min_lod_value.txt")
        String crossCheckMinLODLanes=read_string("crosscheck_min_lod_lanes.txt")
    }
}


task CallSomaticMutations_Prepare_Task {

    # TASK INPUT PARAMS
    File targetIntervals
    File refFasta
    File refFastaIdx
    File refFastaDict

    String nWay

    # RUNTIME INPUT PARAMS
    String preemptible
    String diskGB_boot

    # DEFAULT VALUES
    String default_preemptible = "1"
    String default_diskGB_boot = "15"
    String default_nWay = "10"

    String input_nWay = if nWay != "" then nWay else default_nWay

    parameter_meta {
        nWay : "Number of ways to scatter (MuTect1 and MuTect2)"
        mutectIntervals : "a list of genomic intervals over which MuTect1 will operate"
        refFasta : "FASTA file for the appropriate genome build (Reference sequence file)"
        refFastaIdx : "FASTA file index for the reference genome"
        refFastaDict : "FASTA file dictionary for the reference genome"
    }

    command <<<

        set -euxo pipefail

        # Calculate disk size for all shards of the mutect run
        SIZE_FILE=split_base_sizes_disk.dat

        # Create list of indices for the scatter job
        seq 0 $((${input_nWay}-1)) > indices.dat

        # Run the prepare task that splits the .interval_list file into subfiles
        java -jar /usr/local/bin/GatkScatterGatherPrepare.jar . ${input_nWay} \
        --intervals ${targetIntervals} --reference_sequence ${refFasta}

    >>>

    runtime {
        docker: "gcr.io/broad-getzlab-workflows/cga_production_pipeline:v0.2.ccle"
        bootDiskSizeGb: if diskGB_boot != "" then diskGB_boot else default_diskGB_boot
        preemptible: if preemptible != "" then preemptible else default_preemptible
        memory: "1 GB"
    }

    output {
        Array[File] interval_files=glob("gatk-scatter.*")
        Array[Int] scatterIndices=read_lines("indices.dat")
    }
}


task Mutect1_Task {

    # TASK INPUT PARAMS
    File tumorBam
    File normalBam
    File tumorBamIdx
    File normalBamIdx
    String pairName
    String caseName
    String ctrlName
    File mutectIntervals
    File DB_SNP_VCF
    File DB_SNP_VCF_IDX
    File cosmicVCF
    File readGroupBlackList
    File MuTectNormalPanel
    File refFasta
    File refFastaIdx
    File refFastaDict
    String fracContam
    
    String downsampleToCoverage

    # FILE SIZE
    Int tumorBam_size
    Int normalBam_size
    Int refFasta_size
    Int db_snp_vcf_size

    # RUNTIME INPUT PARAMS
    String preemptible
    String diskGB_boot
    String diskGB_buffer
    String memoryGB
    String cpu

    # DEFAULT VALUES
    String default_cpu = "1"
    String default_memoryGB = "15"
    String default_preemptible = "1"
    String default_diskGB_boot = "15"
    String default_diskGB_buffer = "20"
    String default_downsampleToCoverage = "99999"

    # COMPUTE MEMORY SIZE
    Int machine_memoryGB = if memoryGB != "" then memoryGB else default_memoryGB
    Int command_memoryGB = machine_memoryGB - 1
   
    # COMPUTE DISK SIZE
    Int machine_diskGB_buffer = if diskGB_buffer != "" then diskGB_buffer else default_diskGB_buffer
    Int diskGB = ceil(tumorBam_size + normalBam_size + refFasta_size + db_snp_vcf_size
                + size(mutectIntervals, "G") + size(cosmicVCF, "G") + size(readGroupBlackList, "G")
                + size(MuTectNormalPanel, "G") + select_first([diskGB_buffer, default_diskGB_buffer]))

    String downsample = if downsampleToCoverage != "" then downsampleToCoverage else default_downsampleToCoverage

    parameter_meta {
        tumorBam : "sample tumor BAM file"
        tumorBamIdx : "sample tumor BAI file (indexed BAM file)"
        normalBam : "sample normal BAM file"
        normalBamIdx : "sample normal BAI file (indexed BAM file)"
        pairName : "a string for the name of the pair under analysis used for naming output files"
        caseName : "tumor sample name, prefix for output"
        ctrlName : "normal sample name, prefix for output"
        mutectIntervals : "a list of genomic intervals over which MuTect1 will operate"
        DB_SNP_VCF : "VCF format dbSNP file, used to exclude regions around known polymorphisms from analysis by some PROGRAMs"
        cosmicVCF : " catalogue of somatic mutations in VCF format"
        readGroupBlackList : ""
        MuTectNormalPanel : ""
        refFasta : "FASTA file for the appropriate genome build (Reference sequence file)"
        refFastaIdx : "FASTA file index for the reference genome"
        refFastaDict : "FASTA file dictionary for the reference genome"
        fracContam : "contamination fraction (output from ContEST)"
        downsampleToCoverage : ""
    }

    command <<<

        set -euxo pipefail

        #variable for normal panel
        NORMAL_PANEL_FLAG_AND_VAL=""
        if [ -s "${MuTectNormalPanel}" ] ; then
            NORMAL_PANEL_FLAG_AND_VAL="--normal_panel ${MuTectNormalPanel}" ;
        fi ;

        #mutect 1
        java "-Xmx${command_memoryGB}g" -jar /usr/local/bin/muTect-1.1.6.jar --analysis_type MuTect \
        -L ${mutectIntervals} \
        --normal_sample_name ${ctrlName} \
        -I:normal  ${normalBam} \
        --tumor_sample_name ${caseName} \
        -I:tumor ${tumorBam} \
        --reference_sequence ${refFasta} \
        --fraction_contamination ${fracContam} \
        --dbsnp ${DB_SNP_VCF} \
        --cosmic ${cosmicVCF} \
        --read_group_black_list ${readGroupBlackList} \
        --out ${pairName}.MuTect1.call_stats.txt \
        --coverage_file ${pairName}.MuTect1.coverage.wig.txt \
        --power_file ${pairName}.MuTect1.power.wig.txt \
        --downsample_to_coverage ${downsample} \
        $NORMAL_PANEL_FLAG_AND_VAL

    >>>

    runtime {
        docker: "gcr.io/broad-getzlab-workflows/cga_production_pipeline:v0.2.ccle"
        bootDiskSizeGb: if diskGB_boot != "" then diskGB_boot else default_diskGB_boot
        preemptible: if preemptible != "" then preemptible else default_preemptible
        cpu: if cpu != "" then cpu else default_cpu
        disks: "local-disk ${diskGB} HDD"
        memory: machine_memoryGB + "GB"
    }

    output {
        File mutect1_cs="${pairName}.MuTect1.call_stats.txt"
        File mutect1_pw="${pairName}.MuTect1.power.wig.txt"
        File mutect1_cw="${pairName}.MuTect1.coverage.wig.txt"
    }
}


task Mutect2_Task {

    # TASK INPUT PARAMS
    File tumorBam
    File normalBam
    File tumorBamIdx
    File normalBamIdx
    String pairName
    String caseName
    String ctrlName
    File mutectIntervals
    File readGroupBlackList
    File MuTectNormalPanel
    File refFasta
    File refFastaIdx
    File refFastaDict
    String fracContam
    File GATK4_JAR

    # FILE SIZE
    Int tumorBam_size
    Int normalBam_size
    Int refFasta_size
    Int gatk4_jar_size

    # RUNTIME INPUT PARAMS
    String preemptible
    String diskGB_boot
    String diskGB_buffer
    String memoryGB
    String cpu

    # DEFAULT VALUES
    String default_cpu = "1"
    String default_memoryGB = "15"
    String default_preemptible = "1"
    String default_diskGB_boot = "15"
    String default_diskGB_buffer = "20"

    # COMPUTE MEMORY SIZE
    Int machine_memoryGB = if memoryGB != "" then memoryGB else default_memoryGB
    Int command_memoryGB = machine_memoryGB - 1
   
    # COMPUTE DISK SIZE
    Int machine_diskGB_buffer = if diskGB_buffer != "" then diskGB_buffer else default_diskGB_buffer
    Int diskGB = ceil(tumorBam_size + normalBam_size + refFasta_size + gatk4_jar_size
                + size(mutectIntervals, "G") + size(readGroupBlackList, "G")
                + size(MuTectNormalPanel, "G") + machine_diskGB_buffer)

    parameter_meta {
        tumorBam : "sample tumor BAM file"
        tumorBamIdx : "sample tumor BAI file (indexed BAM file)"
        normalBam : "sample normal BAM file"
        normalBamIdx : "sample normal BAI file (indexed BAM file)"
        pairName : "a string for the name of the pair under analysis used for naming output files"
        caseName : "tumor sample name, prefix for output"
        ctrlName : "normal sample name, prefix for output"
        mutectIntervals : "a list of genomic intervals over which MuTect2 will operate" 
        DB_SNP_VCF : ""
        readGroupBlackList : ""
        MuTectNormalPanel : ""
        refFasta : "FASTA file for reference genome"
        refFastaIdx : "FASTA file index for the reference genome"
        refFastaDict : "FASTA file dictionary for the reference genome"
        fracContam : "contamination fraction (output from ContEST)"        
        GATK4_JAR : ""
    }

    command <<<

        set -euxo pipefail

        #variable for normal panel
        NORMAL_PANEL_FLAG_AND_VAL=""
        if [ -s "${MuTectNormalPanel}" ] ; then
            BZ="${MuTectNormalPanel}.gz"
            #bgzip the file and index it
            bgzip ${MuTectNormalPanel}
            tabix $BZ 
            NORMAL_PANEL_FLAG_AND_VAL="--normal_panel $BZ" ;
        fi ;

        #MuTect2 wants names that match those in the BAMs so grab them from the BAMs
        /usr/local/jre1.8.0_73/bin/java "-Xmx${command_memoryGB}g" -jar ${GATK4_JAR} GetSampleName \
        -I ${tumorBam} -O tumorName.txt
        /usr/local/jre1.8.0_73/bin/java "-Xmx${command_memoryGB}g" -jar ${GATK4_JAR} GetSampleName \
        -I ${normalBam} -O normalName.txt
        NORMAL_NAME=`cat normalName.txt`
        TUMOR_NAME=`cat tumorName.txt`

        #mutect 2 ----- gatk4
        /usr/local/jre1.8.0_73/bin/java "-Xmx${command_memoryGB}g" -jar ${GATK4_JAR} Mutect2 \
        --input ${normalBam} \
        --input ${tumorBam} \
        --tumor-sample "$TUMOR_NAME" \
        --normal-sample "$NORMAL_NAME" \
        --reference ${refFasta} \
        --panel-of-normals $BZ \
        --contamination-fraction-to-filter ${fracContam} \
        --intervals ${mutectIntervals} \
        --output ${pairName}.MuTect2.call_stats.unfiltered.unaf.txt

        #filter the variants
        /usr/local/jre1.8.0_73/bin/java -jar -Xmx4g ${GATK4_JAR} FilterMutectCalls \
        -O ${pairName}.MuTect2.call_stats.filtered.unaf.txt -V ${pairName}.MuTect2.call_stats.unfiltered.unaf.txt

        python /usr/local/bin/process_af.py "${pairName}.MuTect2.call_stats.filtered.unaf.txt"  \
        "${pairName}.MuTect2.call_stats.txt" "$TUMOR_NAME" "$NORMAL_NAME" "${caseName}" "${ctrlName}"

    >>>

    runtime {
        docker: "gcr.io/broad-getzlab-workflows/cga_production_pipeline:v0.2.ccle"
        bootDiskSizeGb: if diskGB_boot != "" then diskGB_boot else default_diskGB_boot
        preemptible: if preemptible != "" then preemptible else default_preemptible
        cpu: if cpu != "" then cpu else default_cpu
        disks: "local-disk ${diskGB} HDD"
        memory: machine_memoryGB + "GB"
    }

    output {
        File mutect2_cs="${pairName}.MuTect2.call_stats.txt"
    }
}


task MutectFC_Task {

    # TASK INPUT PARAMS
    File tumorBam
    File normalBam
    File tumorBamIdx
    File normalBamIdx
    String pairName
    String caseName
    String ctrlName
    File mutectIntervals
    File DB_SNP_VCF
    File DB_SNP_VCF_IDX
    File cosmicVCF
    File readGroupBlackList
    File MuTectNormalPanel
    File refFasta
    File refFastaIdx
    File refFastaDict
    String fracContam
    
    String downsampleToCoverage

    # FILE SIZE
    Int tumorBam_size
    Int normalBam_size
    Int refFasta_size
    Int db_snp_vcf_size

    # RUNTIME INPUT PARAMS
    String preemptible
    String diskGB_boot
    String diskGB_buffer
    String memoryGB
    String cpu

    # DEFAULT VALUES
    String default_cpu = "1"
    String default_memoryGB = "15"
    String default_preemptible = "1"
    String default_diskGB_boot = "15"
    String default_diskGB_buffer = "20"
    String default_downsampleToCoverage = "1000"

    # COMPUTE MEMORY SIZE
    Int machine_memoryGB = if memoryGB != "" then memoryGB else default_memoryGB
    Int command_memoryGB = machine_memoryGB - 1
   
    # COMPUTE DISK SIZE
    Int machine_diskGB_buffer = if diskGB_buffer != "" then diskGB_buffer else default_diskGB_buffer
    Int diskGB = ceil(tumorBam_size + normalBam_size + refFasta_size + db_snp_vcf_size
                + size(mutectIntervals, "G") + size(cosmicVCF, "G")
                + size(readGroupBlackList, "G") + size(MuTectNormalPanel, "G") + machine_diskGB_buffer)

    String downsample = if downsampleToCoverage != "" then downsampleToCoverage else default_downsampleToCoverage

    parameter_meta {
        tumorBam : "sample tumor BAM file"
        tumorBamIdx : "sample tumor BAI file (indexed BAM file)"
        normalBam : "sample normal BAM file"
        normalBamIdx : "sample normal BAI file (indexed BAM file)"
        pairName : "a string for the name of the pair under analysis used for naming output files"
        caseName : "tumor sample name, prefix for output"
        ctrlName : "normal sample name, prefix for output"
        mutectIntervals : "interval list of targets"
        DB_SNP_VCF : "VCF format dbSNP file, used to exclude regions around known polymorphisms from analysis by some PROGRAMs"
        cosmicVCF : " catalogue of somatic mutations in VCF format"
        readGroupBlackList : ""
        MuTectNormalPanel : ""
        refFasta : "FASTA file for reference genome"
        refFastaIdx : "FASTA file index for the reference genome"
        refFastaDict : "FASTA file dictionary for the reference genome"
        fracContam : "contamination fraction (output from ContEST)"
        downsampleToCoverage : ""
    }

    command <<<
        
        set -euxo pipefail

        #variable for normal panel
        NORMAL_PANEL_FLAG_AND_VAL=""
        if [ -s "${MuTectNormalPanel}" ] ; then
            NORMAL_PANEL_FLAG_AND_VAL="--normal_panel ${MuTectNormalPanel}" ;
        fi ;

        java "-Xmx${command_memoryGB}g" -jar /usr/local/bin/muTect-1.1.6.jar --analysis_type MuTect \
        -L ${mutectIntervals} \
        --normal_sample_name ${ctrlName} \
        -I:normal  ${normalBam} \
        --tumor_sample_name ${caseName} \
        -I:tumor ${tumorBam}  \
        --reference_sequence ${refFasta} \
        --fraction_contamination ${fracContam} \
        --dbsnp ${DB_SNP_VCF} \
        --cosmic ${cosmicVCF} \
        --force_output \
        --read_group_black_list ${readGroupBlackList} \
        --out ${pairName}.MuTect1.call_stats.txt \
        --coverage_file ${pairName}.MuTect1.coverage.wig.txt \
        --power_file ${pairName}.MuTect1.power.wig.txt \
        --downsample_to_coverage ${downsample} \
        $NORMAL_PANEL_FLAG_AND_VAL

    >>>

    runtime {
        docker: "gcr.io/broad-getzlab-workflows/cga_production_pipeline:v0.2.ccle"
        bootDiskSizeGb: if diskGB_boot != "" then diskGB_boot else default_diskGB_boot
        preemptible: if preemptible != "" then preemptible else default_preemptible
        cpu: if cpu != "" then cpu else default_cpu
        disks: "local-disk ${diskGB} HDD"
        memory: machine_memoryGB + "GB"
    }

    output {
        File mutectfc_cs="${pairName}.MuTect1.call_stats.txt"
        File mutectfc_pw="${pairName}.MuTect1.power.wig.txt"
        File mutectfc_cw="${pairName}.MuTect1.coverage.wig.txt"
    }
}


task Strelka {

    # TASK INPUT PARAMS
    File tumorBam
    File tumorBamIdx
    File normalBam
    File normalBamIdx
    File refFasta
    File refFastaDict
    File refFastaIdx
    File config
    String pairName

    # FILE SIZE
    Int tumorBam_size
    Int normalBam_size
    Int refFasta_size

    # RUNTIME INPUT PARAMS
    String preemptible
    String diskGB_boot
    String diskGB_buffer
    String memoryGB
    String cpu

    # DEFAULT VALUES
    String default_cpu = "1"
    String default_memoryGB = "25"
    String default_preemptible = "1"
    String default_diskGB_boot = "15"
    String default_diskGB_buffer = "20"

    # COMPUTE MEMORY SIZE
    Int machine_memoryGB = if memoryGB != "" then memoryGB else default_memoryGB
   
    # COMPUTE DISK SIZE
    Int machine_diskGB_buffer = if diskGB_buffer != "" then diskGB_buffer else default_diskGB_buffer
    Int diskGB = ceil(tumorBam_size + normalBam_size + refFasta_size + machine_diskGB_buffer)

    parameter_meta {
        tumorBam : "sample tumor BAM file"
        tumorBamIdx : "sample tumor BAI file (indexed BAM file)"
        normalBam : "sample normal BAM file"
        normalBamIdx : "sample normal BAI file (indexed BAM file)"
        pairName : "a string for the name of the pair under analysis used for naming output files"
        refFasta : "FASTA file for reference genome"
        refFastaIdx : "FASTA file index for the reference genome"
        refFastaDict : "FASTA file dictionary for the reference genome"
        config : "Strelka configuration file"
    }

    command <<<

        set -euxo pipefail

python_cmd="
import subprocess,os
def run(cmd):
    subprocess.check_call(cmd,shell=True)

run('ln -sT `pwd` /opt/execution')
run('ln -sT `pwd`/../inputs /opt/inputs')

##########################

run('ln -vs ${refFasta} ')
run('ln -vs ${refFastaDict} ')
run('ln -vs ${refFastaIdx} ')
REFERENCE = os.path.abspath(os.path.basename('${refFasta}'))

run('ln -vs ${tumorBam} tumor.bam')
run('ln -vs ${tumorBamIdx} tumor.bai')
run('ln -vs ${normalBam} normal.bam')
run('ln -vs ${normalBamIdx} normal.bai')
TBAM = os.path.abspath('tumor.bam')
NBAM = os.path.abspath('normal.bam')

run('ln -vs ${config} ')
CONFIG = os.path.basename('${config}')

# somthing to chew on in to kick off stdout...
run('ls -latrh ')
run('samtools idxstats ' + TBAM)
run('samtools idxstats ' + NBAM)

LIBDIR = '/opt/src/RunStrelka'

run('bash /opt/src/RunStrelka/runStrelka_nolink.sh  ' + LIBDIR + '  ' + TBAM + '  ' + NBAM  + ' ' + REFERENCE + ' ${pairName}  ' + CONFIG )

run('ls -latrh ')

run('tar cvfz ${pairName}.strelka.tar.gz ${pairName}.all.somatic.indels.vcf ${pairName}.passed.somatic.indels.vcf ${pairName}.all.somatic.snvs.vcf ${pairName}.passed.somatic.snvs.vcf')

"
        echo "$python_cmd"
        python -c "$python_cmd"

    >>>

    runtime {
        docker: "gcr.io/broad-getzlab-workflows/cga_production_pipeline:v0.2.ccle"
        bootDiskSizeGb: if diskGB_boot != "" then diskGB_boot else default_diskGB_boot
        preemptible: if preemptible != "" then preemptible else default_preemptible
        cpu: if cpu != "" then cpu else default_cpu
        disks: "local-disk ${diskGB} HDD"
        memory: machine_memoryGB + "GB"
    }

    output {
        File Strelka_targz="${pairName}.strelka.tar.gz"
        File Strelka_passed_indels="${pairName}.passed.somatic.indels.vcf"
        File Strelka_passed_snvs="${pairName}.passed.somatic.snvs.vcf"
        File Strelka_all_indels="${pairName}.all.somatic.indels.vcf"
        File Strelka_all_snvs="${pairName}.all.somatic.snvs.vcf"
    }
}


task Gather_Task {

    # TASK INPUT PARAMS
    Array[File] mutect1_cs
    Array[File] mutect2_cs
    String pairName

    # RUNTIME INPUT PARAMS
    String preemptible
    String diskGB_boot
    String diskGB_buffer
    String memoryGB
    String cpu

    # DEFAULT VALUES
    String default_cpu = "1"
    String default_memoryGB = "15"
    String default_preemptible = "1"
    String default_diskGB_boot = "15"
    String default_diskGB_buffer = "20" # need to increase buffer for storing wiggle files
    
    # COMPUTE MEMORY SIZE
    Int machine_memoryGB = if memoryGB != "" then memoryGB else default_memoryGB
   
    # COMPUTE DISK SIZE
    Int machine_diskGB_buffer = if diskGB_buffer != "" then diskGB_buffer else default_diskGB_buffer
    Int diskGB = ceil(machine_diskGB_buffer)

    parameter_meta {
        mutect1_cs : "MuTect1 call stats"
        mutect2_cs : "MuTect2 call stats (vcf)"
        pairName: "a string for the name of the pair under analysis used for naming output files"
    }

    command <<<

        set -euxo pipefail

python <<CODE

mutect1_cs_file_array = '${sep="," mutect1_cs}'.split(",")
mutect2_cs_file_array = '${sep="," mutect2_cs}'.split(",")
mutect1_cs_list = open('mutect1_cs_list.txt', 'w')
mutect2_cs_list = open('mutect2_cs_list.txt', 'w')

for i in range(len(mutect1_cs_file_array)):
    mutect1_cs_list.write(mutect1_cs_file_array[i] + '\n')
for i in range(len(mutect2_cs_file_array)):
    mutect2_cs_list.write(mutect2_cs_file_array[i] + '\n')
mutect1_cs_list.close()
mutect2_cs_list.close()
CODE

        ######## Files #################

        # MuTect1 Files
        MUTECT1_CS="${pairName}.MuTect1.call_stats.txt"
        # MuTect2 Files
        MUTECT2_CS="${pairName}.MuTect2.call_stats.vcf"
        MUTECT2_INDELS="${pairName}.MuTect2.call_stats.indels.vcf"
        MUTECT2_OTHER="${pairName}.MuTect2.call_stats.other.vcf"

        ######## Merge MuTect1 and MuTect2 Call Stats #################

        # Gather MuTect1 call stats
        python /usr/local/bin/merge_callstats.py "mutect1_cs_list.txt" $MUTECT1_CS
        # Gather MuTect2 call stats
        python /usr/local/bin/merge_callstats.py "mutect2_cs_list.txt" $MUTECT2_CS
        # Filter out indels and SNPs for MuTect2 results
        python /usr/local/bin/filter_indels.py $MUTECT2_CS $MUTECT2_INDELS $MUTECT2_OTHER

    >>>

    runtime {
        docker: "gcr.io/broad-getzlab-workflows/cga_production_pipeline:v0.2.ccle"
        bootDiskSizeGb: if diskGB_boot != "" then diskGB_boot else default_diskGB_boot
        preemptible: if preemptible != "" then preemptible else default_preemptible
        cpu: if cpu != "" then cpu else default_cpu
        disks: "local-disk ${diskGB} HDD"
        memory: machine_memoryGB + "GB"
    }

    output {
        File MUTECT1_CS_SNV="${pairName}.MuTect1.call_stats.txt"
        File MUTECT2_VCF_ALL="${pairName}.MuTect2.call_stats.vcf"
        File MUTECT2_VCF_INDELS="${pairName}.MuTect2.call_stats.indels.vcf"
    }
}


task DeTiN_Task {

    # TASK INPUT PARAMS
    File MUTECT1_CS
    File MUTECT2_INDELS
    File tumor_hets
    File seg_file
    File normal_hets
    String pairName
    File exac_pickle

    String TiN_prior
    String Mutation_prior
    String release_version

    # RUNTIME INPUT PARAMS
    String preemptible
    String diskGB_boot
    String diskGB_buffer
    String memoryGB
    String cpu

    # DEFAULT VALUES
    String default_cpu = "1"
    String default_memoryGB = "15"
    String default_preemptible = "1"
    String default_diskGB_boot = "15"
    String default_diskGB_buffer = "20"
    String default_TiN_prior = "0.2"
    String default_Mutation_prior = "0.05"
    String default_release_version = ""

    # COMPUTE MEMORY SIZE
    Int machine_memoryGB = if memoryGB != "" then memoryGB else default_memoryGB
   
    # COMPUTE DISK SIZE
    Int machine_diskGB_buffer = if diskGB_buffer != "" then diskGB_buffer else default_diskGB_buffer
    Int diskGB = ceil(size(seg_file, "G") + size(tumor_hets, "G") + size(normal_hets, "G")
                    + size(exac_pickle, "G") + size(MUTECT1_CS, "G") + size(MUTECT2_INDELS, "G") + machine_diskGB_buffer)

    String Prior_TiN = if TiN_prior != "" then TiN_prior else default_TiN_prior
    String Prior_Mutation = if Mutation_prior != "" then Mutation_prior else default_Mutation_prior
    String deTiN_version = if release_version != "" then release_version else default_release_version

    parameter_meta {
        pairName: "a string for the name of the pair under analysis used for naming output files"
        seg_file : "filename pointing to the input - either a HAPSEG file or a segmentation file"
        tumor_hets : "a tab delimited file with read counts from the tumor sample."
        normal_hets : "a tab delimited file with read counts from the normal sample."
        exac_pickle : ""
        TiN_prior : ""
        Mutation_prior : ""
    }

    command <<<

        set -euxo pipefail

        ########### If specific version of deTiN provided pull it from GitHub ###########################
        if [ "${deTiN_version}" != "" ] ; then
            git clone https://github.com/broadinstitute/deTiN.git
            git -C deTiN/ checkout tags/${deTiN_version}
        fi ;

        ########### Run deTiN on MuTect 1 & 2 ###########################
        # generate null outputs
        touch ${pairName}.TiN_estimate.txt
        touch ${pairName}.TiN_estimate_CI.txt
        touch ${pairName}.deTiN_aSCNAs.txt
        touch ${pairName}.deTiN_SSNVs.txt
        touch ${pairName}.deTiN_indels.txt
        touch ${pairName}_SSNVs_plot.png
        touch ${pairName}_KmeansEval_plot.png
        touch ${pairName}_TiN_models_plot.png
        touch ${pairName}_KmeansEval_scatter_plot.png

        python /root/deTiN/deTiN/deTiN.py --mutation_data_path ${MUTECT1_CS} --cn_data_path ${seg_file} \
        --tumor_het_data ${tumor_hets} --normal_het_data ${normal_hets} --exac_data_path ${exac_pickle} \
        --output_name ${pairName} --TiN_prior ${Prior_TiN} --mutation_prior ${Prior_Mutation} \
        --indel_data_path ${MUTECT2_INDELS} --indel_data_type "MuTect2"
        
    >>>

    runtime {
        docker: "gcr.io/broad-getzlab-workflows/cga_production_pipeline:v0.2.ccle"
        bootDiskSizeGb: if diskGB_boot != "" then diskGB_boot else default_diskGB_boot
        preemptible: if preemptible != "" then preemptible else default_preemptible
        cpu: if cpu != "" then cpu else default_cpu
        disks: "local-disk ${diskGB} HDD"
        memory: machine_memoryGB + "GB"
    }

    output {
        Float TiN=read_float("${pairName}.TiN_estimate.txt")
        Int number_added_SSNV=read_int("${pairName}.number_of_SSNVs_added.txt")
        String TiN_CI=read_string("${pairName}.TiN_estimate_CI.txt")
        File deTiN_call_stats="${pairName}.deTiN_SSNVs.txt"
        File deTiN_indels="${pairName}.deTiN_indels.txt"
        File deTiN_SSNVs_plot="${pairName}_SSNVs_plot.png"
        File deTiN_aSCNA_kmeans_RSS_plot="${pairName}_KmeansEval_plot.png"
        File deTiN_aSCNA_scatter_plot="${pairName}_KmeansEval_scatter_plot.png"
        File deTiN_TiN_modes_plot="${pairName}_TiN_models_plot.png"        
        File aSCNA_model="${pairName}_TiN_hets_aSCNA_model.png"
        File deTiN_segments="${pairName}.deTiN_aSCNAs.txt"
    }
}

task VEP_Task {

    # TASK INPUT PARAMS
    File MUTECT1_CS
    File MUTECT2_VCF
    File STRELKA_VCF
    String pairName
    String caseName
    String ctrlName
    File VEP_File    
    File GNOMAD_FILE
    File GNOMAD_FILE_IDX
    File oncoDBTarBall_JustRef

    # RUNTIME INPUT PARAMS
    String preemptible
    String diskGB_boot
    String diskGB_buffer
    String memoryGB
    String cpu

    # DEFAULT VALUES
    String default_cpu = "1"
    String default_memoryGB = "10"
    String default_preemptible = "1"
    String default_diskGB_boot = "15"
    String default_diskGB_buffer = "50"

    # COMPUTE MEMORY SIZE
    Int machine_memoryGB = if memoryGB != "" then memoryGB else default_memoryGB
    Int command_memoryGB = machine_memoryGB - 1

    # COMPUTE DISK SIZE
    Int machine_diskGB_buffer = if diskGB_buffer != "" then diskGB_buffer else default_diskGB_buffer
    Int diskGB = ceil(size(MUTECT1_CS, "G") + size(MUTECT2_VCF, "G") + size(STRELKA_VCF, "G")
                + size(GNOMAD_FILE, "G")*5 + size(GNOMAD_FILE_IDX, "G") + size(VEP_File, "G") * 5 
                + machine_diskGB_buffer)

    parameter_meta {
        MUTECT1_CS : ""
        MUTECT2_VCF : ""
        STRELKA_VCF : ""
        pairName: "a string for the name of the pair under analysis used for naming output files"
        caseName : "tumor sample name, prefix for output"
        ctrlName : "normal sample name, prefix for output"
        VEP_File : ""
        GNOMAD_FILE : ""      
    }

    command <<<

        set -x

        ############## Pre-process MuTect1 call stats #################################
        # Running Oncotator to Convert MuTect1 callstats to VCF 

        MUTECT1_CS_PASSED="${pairName}.MuTect1.call_stats.passed.txt"
        MUTECT1_CS_REJECTED="${pairName}.MuTect1.call_stats.rejected.txt"
        MUTECT1_CS_MAFLITE="${pairName}.MuTect1.call_stats.maflite"
        MUTECT1_VCF="${pairName}.MuTect1.call_stats.maflite.annotated.vcf"

        # Filter MuTect1 mutation calls that passed filter
        python /usr/local/bin/filter_passed_mutations.py ${MUTECT1_CS} $MUTECT1_CS_PASSED $MUTECT1_CS_REJECTED "KEEP"
        # Convert MuTect1 call stats to MafLite
        python /usr/local/bin/callstats_to_maflite.py $MUTECT1_CS_PASSED $MUTECT1_CS_MAFLITE

        ########################### Unzip Oncotator database ###############################

        #find TARBALL type 
        # TODO Find better way to get extension
        TYPE=`echo 'if("${oncoDBTarBall_JustRef}"=~m/z$/) { print "GZ" ; } else { print "TAR" ; } '| perl` ;

        #obtain the name of the directory for oncodb and unpack based on the TYPE
        if [ "$TYPE" == "GZ" ] ; then
        # TODO Find better way to get name of the file
            ONCO_DB_DIR_NAME=`gunzip -c ${oncoDBTarBall_JustRef} |tar -tf /dev/stdin|head -1` ;
            tar -xzf ${oncoDBTarBall_JustRef}
        else
            ONCO_DB_DIR_NAME=`tar -tf ${oncoDBTarBall_JustRef} |head -1` ;
            tar -xf ${oncoDBTarBall_JustRef} ;
        fi ;

        ################## Annotate MuTect1, MuTect2, Strelka call stats ###########################

        # Annotate MuTect1 call stats (MAFLITE to VCF)
        # --infer-onps \        
        /root/oncotator_venv/bin/oncotator \
        -i MAFLITE \
        -o VCF \
        --db-dir `pwd`/$ONCO_DB_DIR_NAME \
        --infer_genotypes yes \
        -a normal_barcode:${ctrlName} \
        -a tumor_barcode:${caseName} \
        $MUTECT1_CS_PASSED $MUTECT1_VCF hg19

        ############### Run Variant Effector Predictor (VEP)  #########################

        #make a link from the home directory to the current directory to avoid running out of disk space on the boot disk
        mkdir -v vep_data_dir
        #delete the existing directory first to make a successful link
        rm -rf ~/.vep
        ln -vs `pwd`/vep_data_dir ~/.vep

        #In either case unpack the data into the home directory where VEP expects to find it
        IS_ZIP=`echo ${VEP_File}|grep -Pic '\.zip$'` ;
        if [ "$IS_ZIP" -eq "1" ] ;
        then
            #it's a zip file
            unzip -d ~ ${VEP_File} ;
        else
            #tar ball
            tar -C ~ -xvzf ${VEP_File}
        fi ;

        # VEP for MuTect1
        MUTECT1_VEP="${pairName}.mutect1.vep_annotated.vcf"
        MUTECT1_VEP_filtered="${pairName}.mutect1.vep_annotated.filtered.vcf"
        /ensembl-tools-release-83/ensembl-tools-release-83/scripts/variant_effect_predictor/variant_effect_predictor.pl \
        --custom ${GNOMAD_FILE},gnomADg,vcf,exact,0,GT,AC,AF,AN,AF_AFR,AF_AMR,AF_ASJ,AF_EAS,AF_FIN,AF_NFE,AF_OTH \
        --input_file $MUTECT1_VCF \
        --output_file $MUTECT1_VEP \
        --vcf \
        --symbol \
        --cache \
        --offline \
        --failed 1

        python /usr/local/bin/vep_filter_germline.py $MUTECT1_VEP $MUTECT1_VEP_filtered

        # VEP for MuTect2
        MUTECT2_VEP="${pairName}.mutect2.vep_annotated.vcf"
        MUTECT2_VEP_filtered="${pairName}.mutect2.vep_annotated.filtered.vcf"
        /ensembl-tools-release-83/ensembl-tools-release-83/scripts/variant_effect_predictor/variant_effect_predictor.pl \
        --custom ${GNOMAD_FILE},gnomADg,vcf,exact,0,GT,AC,AF,AN,AF_AFR,AF_AMR,AF_ASJ,AF_EAS,AF_FIN,AF_NFE,AF_OTH \
        --input_file ${MUTECT2_VCF} \
        --output_file $MUTECT2_VEP \
        --vcf \
        --symbol \
        --cache \
        --offline \
        --failed 1

        python /usr/local/bin/vep_filter_germline.py $MUTECT2_VEP $MUTECT2_VEP_filtered

        # VEP for Strelka
        STRELKA_VEP="${pairName}.strelka.vep_annotated.vcf"
        STRELKA_VEP_filtered="${pairName}.strelka.vep_annotated.filtered.vcf"
        /ensembl-tools-release-83/ensembl-tools-release-83/scripts/variant_effect_predictor/variant_effect_predictor.pl \
        --custom ${GNOMAD_FILE},gnomADg,vcf,exact,0,GT,AC,AF,AN,AF_AFR,AF_AMR,AF_ASJ,AF_EAS,AF_FIN,AF_NFE,AF_OTH \
        --input_file ${STRELKA_VCF} \
        --output_file $STRELKA_VEP \
        --vcf \
        --symbol \
        --cache \
        --offline \
        --failed 1

        python /usr/local/bin/vep_filter_germline.py $STRELKA_VEP $STRELKA_VEP_filtered

    >>>

    runtime {
        docker: "gcr.io/broad-getzlab-workflows/cga_production_pipeline:v0.2.ccle"
        bootDiskSizeGb: if diskGB_boot != "" then diskGB_boot else default_diskGB_boot
        preemptible: if preemptible != "" then preemptible else default_preemptible
        cpu: if cpu != "" then cpu else default_cpu
        disks: "local-disk ${diskGB} HDD"
        memory: machine_memoryGB + "GB"
    }

    output {
        File MUTECT1_VEP_annotated_vcf="${pairName}.mutect1.vep_annotated.vcf"
        File MUTECT2_VEP_annotated_vcf="${pairName}.mutect2.vep_annotated.vcf"
        File STRELKA_VEP_annotated_vcf="${pairName}.strelka.vep_annotated.vcf"  
        File MUTECT1_VEP_annotated_filtered_vcf="${pairName}.mutect1.vep_annotated.filtered.vcf"
        File MUTECT2_VEP_annotated_filtered_vcf="${pairName}.mutect2.vep_annotated.filtered.vcf"
        File STRELKA_VEP_annotated_filtered_vcf="${pairName}.strelka.vep_annotated.filtered.vcf"       
    }
}


task GatherWIGFiles_Task {

    # TASK INPUT PARAMS
    Array[File] mutect1_pw
    Array[File] mutect1_cw
    String pairName
    
    # RUNTIME INPUT PARAMS
    String preemptible
    String diskGB_boot
    String diskGB_buffer
    String memoryGB
    String cpu

    # DEFAULT VALUES
    String default_cpu = "1"
    String default_memoryGB = "15"
    String default_preemptible = "1"
    String default_diskGB_boot = "15"
    String default_diskGB_buffer = "20"

    # COMPUTE MEMORY SIZE
    Int machine_memoryGB = if memoryGB != "" then memoryGB else default_memoryGB
   
    # COMPUTE DISK SIZE
    Int machine_diskGB_buffer = if diskGB_buffer != "" then diskGB_buffer else default_diskGB_buffer
    Int diskGB = ceil(machine_diskGB_buffer)

    parameter_meta {
        mutect1_pw : "List of MuTect1 power wiggle files"
        mutect1_cw : "List of MuTect1 coverage wiggle files"
        pairName: "a string for the name of the pair under analysis used for naming output files"
    }

    command <<<

        set -euxo pipefail

python <<CODE

mutect1_pw_file_array = '${sep="," mutect1_pw}'.split(",")
mutect1_cw_file_array = '${sep="," mutect1_cw}'.split(",")

with open('mutect1_pw_list.txt', 'w') as mutect1_pw_list:
    for i in range(len(mutect1_pw_file_array)):
        mutect1_pw_list.write(mutect1_pw_file_array[i] + '\n')
with open('mutect1_cw_list.txt', 'w') as mutect1_cw_list:
    for i in range(len(mutect1_cw_file_array)):
        mutect1_cw_list.write(mutect1_cw_file_array[i] + '\n')

CODE

        ######## Files #################

        # MuTect1 Files
        MUTECT1_CW="${pairName}.MuTect1.merged.coverage.wig.txt"
        # MuTect1 Files
        MUTECT1_PW="${pairName}.MuTect1.merged.power.wig.txt"
        
        ######## Merge MuTect1 power and coverage wiggle files #################

        # Gather MuTect1 power wig files
        python /usr/local/bin/merge_wig_files.py "mutect1_pw_list.txt" $MUTECT1_PW
        # Gather MuTect1 coverage wig files
        python /usr/local/bin/merge_wig_files.py "mutect1_cw_list.txt" $MUTECT1_CW

        ########### Zip merged files ###############

        zip "${pairName}.MuTect1.merged.power.wig.zip" "$MUTECT1_PW"
        zip "${pairName}.MuTect1.merged.coverage.wig.zip" "$MUTECT1_CW"

    >>>

    runtime {
        docker         : "gcr.io/broad-getzlab-workflows/cga_production_pipeline:v0.2.ccle"
        bootDiskSizeGb : if diskGB_boot != "" then diskGB_boot else default_diskGB_boot
        preemptible    : if preemptible != "" then preemptible else default_preemptible
        cpu            : if cpu != "" then cpu else default_cpu
        disks          : "local-disk ${diskGB} HDD"
        memory         : machine_memoryGB + "GB"
    }

    output {
        File MuTect1_merged_power_wig="${pairName}.MuTect1.merged.power.wig.zip"
        File MuTect1_merged_coverage_wig="${pairName}.MuTect1.merged.coverage.wig.zip"
    }
}


task Oncotate_Task {

    # TASK INPUT PARAMS
    File MUTECT1_CS
    File MUTECT2_INDELS
    File STRELKA_INDELS
    String pairName
    String caseName
    String ctrlName
    File oncoDBTarBall

    # RUNTIME INPUT PARAMS
    String preemptible
    String diskGB_boot
    String diskGB_buffer
    String memoryGB
    String cpu

    # DEFAULT VALUES
    String default_cpu = "1"
    String default_memoryGB = "15"
    String default_preemptible = "1"
    String default_diskGB_boot = "15"
    String default_diskGB_buffer = "20"

    # COMPUTE MEMORY SIZE
    Int machine_memoryGB = if memoryGB != "" then memoryGB else default_memoryGB
   
    # COMPUTE DISK SIZE
    Int machine_diskGB_buffer = if diskGB_buffer != "" then diskGB_buffer else default_diskGB_buffer
    Int diskGB = ceil(size(oncoDBTarBall, "G") * 5 + machine_diskGB_buffer)

    parameter_meta {
        MUTECT1_CS : ""
        MUTECT2_INDELS : ""
        STRELKA_INDELS : ""
        pairName: "a string for the name of the pair under analysis used for naming output files"
        caseName : "tumor sample name, prefix for output"
        ctrlName : "normal sample name, prefix for output"
        oncoDBTarBall : ""
    }

    command <<<

        set -x

        ######## Files #################

        # MuTect1 files
        MUTECT1_CS_PASSED="${pairName}.MuTect1.call_stats.passed.txt"
        MUTECT1_CS_REJECTED="${pairName}.MuTect1.call_stats.rejected.txt"
        MUTECT1_CS_MAFLITE="${pairName}.MuTect1.call_stats.maflite"
        MUTECT1_CS_ANNOTATED_MAF="${pairName}.MuTect1.call_stats.annotated.maf"

        # MuTect2 files
        MUTECT2_INDELS_PASSED="${pairName}.MuTect2.call_stats.passed.vcf"
        MUTECT2_INDELS_REJECTED="${pairName}.MuTect2.call_stats.rejected.vcf"
        MUTECT2_INDELS_ANNOTATED_MAF="${pairName}.MuTect2.call_stats.indels.annotated.maf"
        
        # Strelka files
        STRELKA_REFORMATTED_VCF="${pairName}.Strelka.call_stats.re_formatted.vcf"
        STRELKA_ANNOTATED_MAF="${pairName}.Strelka.call_stats.annotated.maf"

        # Merged files
        MUTECT2_STRELKA_MERGED_MAF="${pairName}.MuTect2_INDEL.Strelka_INDEL.merged.maf" 
        MUTECT1_MUTECT2_STRELKA_MERGED_MAF="${pairName}.MuTect1_SNV.MuTect2_INDEL.Strelka_INDEL.annotated.maf"
        MUTECT1_MUTECT2_STRELKA_ANNOTATED_VCF="${pairName}.MuTect1_SNV.MuTect2_INDEL.Strelka_INDEL.annotated.vcf"

        ######## Processing call stats before annotation #################

        # Filter MuTect1 mutation calls that passed filter
        python /usr/local/bin/filter_passed_mutations.py ${MUTECT1_CS} $MUTECT1_CS_PASSED $MUTECT1_CS_REJECTED "PASS"
        # Convert MuTect1 call stats to MafLite
        #

        # Filter MuTect2 mutation calls that passed filter
        python /usr/local/bin/filter_passed_mutations.py ${MUTECT2_INDELS} $MUTECT2_INDELS_PASSED $MUTECT2_INDELS_REJECTED "PASS"

        # Edit Strelka VCF (adding AD and AF, replacing TUMOR/NORMAL with caseName/ctrlName)
        python /usr/local/bin/strelka_allelic_count.py ${STRELKA_INDELS} $STRELKA_REFORMATTED_VCF "${caseName}" "${ctrlName}"

        ########################### Unzip Oncotator database ###############################

        #find TARBALL type 
        # TODO Find better way to get extension
        TYPE=`echo 'if("${oncoDBTarBall}"=~m/z$/) { print "GZ" ; } else { print "TAR" ; } '| perl` ;

        #obtain the name of the directory for oncodb and unpack based on the TYPE
        if [ "$TYPE" == "GZ" ] ; then
        # TODO Find better way to get name of the file
            ONCO_DB_DIR_NAME=`gunzip -c ${oncoDBTarBall} |tar -tf /dev/stdin|head -1` ;
            tar -xzf ${oncoDBTarBall}
        else
            ONCO_DB_DIR_NAME=`tar -tf ${oncoDBTarBall} |head -1` ;
            tar -xf ${oncoDBTarBall} ;
        fi ;

        ################## Annotate MuTect1, MuTect2, Strelka call stats ###########################

        # Annotate MuTect1 call stats (VCF to TCGAMAF)
        # --infer-onps \        
        /root/oncotator_venv/bin/oncotator -i VCF --db-dir `pwd`/$ONCO_DB_DIR_NAME \
        --skip-no-alt \
        --longer-other-tx \
        -a normal_barcode:${ctrlName} \
        -a tumor_barcode:${caseName} \
        $MUTECT1_CS_PASSED $MUTECT1_CS_ANNOTATED_MAF hg19
        

        # Annotate MuTect2 call stats (VCF to TCGAMAF)
        /root/oncotator_venv/bin/oncotator -i VCF --db-dir `pwd`/$ONCO_DB_DIR_NAME \
        --longer-other-tx \
        --skip-no-alt \
        -a normal_barcode:${ctrlName} \
        -a tumor_barcode:${caseName} \
        $MUTECT2_INDELS_PASSED $MUTECT2_INDELS_ANNOTATED_MAF hg19

        # Annotate Strelka call stats (VCF to TCGAMAF)
        /root/oncotator_venv/bin/oncotator -i VCF --db-dir `pwd`/$ONCO_DB_DIR_NAME \
        --longer-other-tx \
        --skip-no-alt \
        -a normal_barcode:${ctrlName} \
        -a tumor_barcode:${caseName} \
        $STRELKA_REFORMATTED_VCF $STRELKA_ANNOTATED_MAF hg19

        ####################### Merge MuTect1, MuTect2, Strelka annotated TCGAMAFs into one file ####################

        # Merge Strelka and MuTect2 indels (only Strelka calls selected, MuTect2 and Strelka common calls are annotated)
        python /usr/local/bin/merge_strelka_mutect2.py $STRELKA_ANNOTATED_MAF $MUTECT2_INDELS_ANNOTATED_MAF $MUTECT2_STRELKA_MERGED_MAF
        # Merge MuTect1 SNVs and Strelka & MuTect2 indels into one MAF file
        python /usr/local/bin/maf_merge.py $MUTECT1_CS_ANNOTATED_MAF $MUTECT2_STRELKA_MERGED_MAF $MUTECT1_MUTECT2_STRELKA_MERGED_MAF

    >>>

    runtime {
        docker: "gcr.io/broad-getzlab-workflows/cga_production_pipeline:v0.2.ccle"
        bootDiskSizeGb: if diskGB_boot != "" then diskGB_boot else default_diskGB_boot
        preemptible: if preemptible != "" then preemptible else default_preemptible
        cpu: if cpu != "" then cpu else default_cpu
        disks: "local-disk ${diskGB} HDD"
        memory: machine_memoryGB + "GB"
    }

    output {
        File WXS_Mutation_M1_SNV_M2_INDEL_Strelka_INDEL_annotated_maf="${pairName}.MuTect1_SNV.MuTect2_INDEL.Strelka_INDEL.annotated.maf"               
    }
}


task OrientationBias_filter_Task {

    # TASK INPUT PARAMS
    File tumorBam
    File tumorBamIdx
    File? detailMetrics
    File MAF
    String pairName
    String stub
    File refFasta
    File GATK4_JAR

    # FILE SIZE
    Int tumorBam_size
    Int refFasta_size
    Int gatk4_jar_size

    # RUNTIME INPUT PARAMS
    String preemptible
    String diskGB_boot
    String diskGB_buffer
    String memoryGB
    String cpu

    # DEFAULT VALUES
    String default_cpu = "1"
    String default_memoryGB = "7"
    String default_preemptible = "1"
    String default_diskGB_boot = "15"
    String default_diskGB_buffer = "20"

    Boolean found_detailMetrics = defined(detailMetrics)

    # COMPUTE MEMORY SIZE
    Int machine_memoryGB = if memoryGB != "" then memoryGB else default_memoryGB
    Int command_memoryGB = machine_memoryGB - 1

    # COMPUTE DISK SIZE
    Int machine_diskGB_buffer = if diskGB_buffer != "" then diskGB_buffer else default_diskGB_buffer
    Int diskGB = ceil(tumorBam_size + refFasta_size + gatk4_jar_size + size(MAF, "G") 
                      + machine_diskGB_buffer)

    parameter_meta {        
        tumorBam: "sample tumor BAM file"
        tumorBamIdx : "sample tumor BAI file (indexed BAM)"
        MAF : "filename pointing to a mutation annotation format (MAF) file (data for somatic point mutations)"
        detailMetrics : "output from the Picard Multiple Metrics CollectSequencingArtifactMetrics run with the settings here; this allows for passthrough instead of recomputing"       
        refFasta : "Reference that was used to align BAM"
        stub : "string used to indicate in the output the effect name"
    }

    command <<<

        set -euxo pipefail

        SNV_MAF="${pairName}.snv.maf"
        INDEL_MAF="${pairName}.indel.maf"
        python /usr/local/bin/split_maf_indel_snp.py -i ${MAF} -o $SNV_MAF -f Variant_Type -v "SNP|DNP|TNP|MNP"
        python /usr/local/bin/split_maf_indel_snp.py -i ${MAF} -o $INDEL_MAF -f Variant_Type -v "INS|DEL"

        ################################

        python /usr/local/bin/get_context_ref_alt_alleles.py -s ${stub} -i ${pairName}

        CONTEXT=$( cat "${pairName}.context.${stub}.txt" )
        ARTIFACT_ALLELE=$( cat "${pairName}.artifact_allele.${stub}.txt" )

        if [ -s "${detailMetrics}" ] ;
        then
            DETAIL_METRICS_FILE=${detailMetrics}
        else
            /usr/local/jre1.8.0_73/bin/java "-Xmx${command_memoryGB}g" -jar ${GATK4_JAR} CollectSequencingArtifactMetrics \
            --INPUT ${tumorBam} \
            --OUTPUT ${pairName} \
            --REFERENCE_SEQUENCE ${refFasta}
            DETAIL_METRICS_FILE="${pairName}.pre_adapter_detail_metrics"  
        fi ;

        # Now parsing metrics for Q value
        python /usr/local/orientationBiasFilter/annotate_orientationBiasQ.py -i ${pairName} -m $DETAIL_METRICS_FILE -c $CONTEXT -a $ARTIFACT_ALLELE

        Q=$( cat "${pairName}.orientation_BiasQ.txt" )

        # Appending Q value to MAF
        python /usr/local/orientationBiasFilter/AppendAnnotation2MAF.py -i ${pairName} -m $SNV_MAF -f ${stub}_Q -v $Q

        OUTPUT_INTERVAL_FILE="${pairName}.intervals"
        OUTPUT_INFO_FILE="${pairName}.orientation_info.txt" 

        python /usr/local/orientationBiasFilter/write_interval_file.py "${pairName}.${stub}_Q.maf.annotated" $OUTPUT_INTERVAL_FILE       

        java "-Xmx${command_memoryGB}g" -jar /usr/local/orientationBiasFilter/GenomeAnalysisTK.jar \
        --analysis_type OxoGMetrics \
        -R ${refFasta} \
        -I ${tumorBam} \
        -L $OUTPUT_INTERVAL_FILE \
        -o $OUTPUT_INFO_FILE

        # Now appending orientation bias information to the MAF
        python /usr/local/orientationBiasFilter/AppendOrientationBiasFields2MAF.py \
        -i ${pairName} -m "${pairName}.${stub}_Q.maf.annotated" -b ${tumorBam} \
        -f $OUTPUT_INFO_FILE

        REF_ALLELE_COMP=$( cat "${pairName}.ref_allele_compliment.${stub}.txt" )
        ARTIFACT_ALLELE_COMP=$( cat "${pairName}.artifact_allele_compliment.${stub}.txt")

        bash -c "source /matlab_source_file_2012a.sh && /usr/local/orientationBiasFilter/orientationBiasFilter ${pairName}.OrientationBiasInfo.maf \
        ${pairName}.OrientationBiasFilter.maf . '0' '1' '0.96' '0.01' '-1' '30' '1.5' \
        $REF_ALLELE_COMP $ARTIFACT_ALLELE_COMP i_${stub}" ;

        #########################################

        #merge back indels into OBF output
        python /usr/local/bin/tsvConcatFiles.py $INDEL_MAF "${pairName}.OrientationBiasFilter.maf" \
        --outputFilename="${pairName}.OrientationBiasFilter.${stub}.indel_snp_merged.filtered.maf"

        python /usr/local/bin/tsvConcatFiles.py $INDEL_MAF "${pairName}.OrientationBiasFilter.unfiltered.maf" \
        --outputFilename="${pairName}.OrientationBiasFilter.${stub}.indel_snp_merged.unfiltered.maf"

        python /usr/local/bin/add_judgement_column.py \
        --input "${pairName}.OrientationBiasFilter.${stub}.indel_snp_merged.unfiltered.maf" \
        --output "${pairName}.OrientationBiasFilter.${stub}.indel_snp_merged.unfiltered.with_annotations.maf" \
        --column "i_${stub}_cut" \
        --pass_flag "0"
        
        zip -r ${pairName}.${stub}_OBF_figures.zip ./figures

    >>>

    runtime {
        docker: "gcr.io/broad-getzlab-workflows/cga_production_pipeline:v0.2.ccle"
        bootDiskSizeGb: if diskGB_boot != "" then diskGB_boot else default_diskGB_boot
        preemptible: if preemptible != "" then preemptible else default_preemptible
        cpu: if cpu != "" then cpu else default_cpu
        disks: "local-disk ${diskGB} HDD"
        memory: machine_memoryGB + "GB"
    }

    output {      
        Float q_val=read_float("${pairName}.orientation_BiasQ.txt")
        File OBF_figures="${pairName}.${stub}_OBF_figures.zip"
        Int num_passed_mutations=read_int("${pairName}.OrientationBiasFilter.maf.pass_count.txt")
        Int num_rejected_mutations=read_int("${pairName}.OrientationBiasFilter.maf.reject_count.txt")
        File WXS_Mutation_OBF_filtered_maf="${pairName}.OrientationBiasFilter.${stub}.indel_snp_merged.filtered.maf"
        File WXS_Mutation_OBF_unfiltered_maf="${pairName}.OrientationBiasFilter.${stub}.indel_snp_merged.unfiltered.maf"
        File WXS_Mutation_OBF_unfiltered_maf_with_annotations="${pairName}.OrientationBiasFilter.${stub}.indel_snp_merged.unfiltered.with_annotations.maf"
    }
}


task MAFPonFilter {

    # TASK INPUT PARAMS
    File MAFFile
    String pairName
    File PONFile
    File cytoBandFile
    String stub

    String TOTNStr
    String NMIN
    String THRESH
    String CODING_ONLY
    String MIN_ALT_COUNT
    String public_release

    # RUNTIME INPUT PARAMS
    String preemptible
    String diskGB_boot
    String diskGB_buffer
    String memoryGB
    String cpu

    # DEFAULT VALUES
    String default_cpu = "1"
    String default_memoryGB = "10"
    String default_preemptible = "1"
    String default_diskGB_boot = "15"
    String default_diskGB_buffer = "20"
    String default_NMIN = "1"
    String default_THRESH = "-2.5"
    String default_CODING_ONLY = "0"
    String default_MIN_ALT_COUNT = "3"
    String default_TOTNStr = "TN"
    String default_public_release = "false"

    # COMPUTE MEMORY SIZE
    Int machine_memoryGB = if memoryGB != "" then memoryGB else default_memoryGB

    # COMPUTE DISK SIZE
    Int machine_diskGB_buffer = if diskGB_buffer != "" then diskGB_buffer else default_diskGB_buffer
    Int diskGB = ceil(size(MAFFile, "G") + size(PONFile, "G") + size(cytoBandFile, "G") + machine_diskGB_buffer)

    String input_TOTNStr = if TOTNStr != "" then TOTNStr else default_TOTNStr
    String input_NMIN = if NMIN != "" then NMIN else default_NMIN
    String input_THRESH = if THRESH != "" then THRESH else default_THRESH
    String input_CODING_ONLY = if CODING_ONLY != "" then CODING_ONLY else default_CODING_ONLY
    String input_MIN_ALT_COUNT  = if MIN_ALT_COUNT != "" then MIN_ALT_COUNT else default_MIN_ALT_COUNT
    String input_public_release = if public_release != "" then public_release else default_public_release

    parameter_meta {
        MAFFile : "input MAF for PON filter analysis"
        PONFile : "formatted panel-of-normals file"
        cytoBandFile : "TSV file of chromosomal annotation: chr, start, end, band, stain"
        pairName : "used to name the output files and title other output"
        TOTNStr : "indicating pair status : can be 'TO' or 'TN'"
        NMIN : "minimum count of samples in a token PoN bin that are used to estimate the PoN log-likelihood (to reduce the effect of tumor contamination in the PoN)"
        thresh : "threshold for pass/fail with log likelihood"
        WCUT : "threshold for pass/fail with pon-computed weight"
        CODING_ONLY : "analyze coding regions only?"
        MIN_ALT_COUNT : "the minimum t_alt_count for mutations in the maf that pass the maf_pon_filter"
    }

    command <<<

        set -euxo pipefail
 
        #the cytoBand file should be in the expected location
        mkdir -pv /xchip/cga/reference/annotation/db/ucsc/hg19/
        cp -v  ${cytoBandFile} /xchip/cga/reference/annotation/db/ucsc/hg19/cytoBand.txt
        ls -alht /xchip/cga/reference/annotation/db/ucsc/hg19/cytoBand.txt

        #run the filter
        #Note "." is used for the parameter file whose usage is *not* exposed/functional in this WDL
        bash -c "source /matlab_source_file_2013a.sh && /usr/local/bin/maf_pon_filter \
        ${MAFFile} ${PONFile} ${pairName} ${input_TOTNStr} . ${input_NMIN} ${input_THRESH} \
        0.5 . ${input_CODING_ONLY} ${input_MIN_ALT_COUNT}" ;

        # Count number of passed and rejected mutations
        python /usr/local/bin/count_variants.py "${pairName}.pon_annotated.pass.maf" "${pairName}.count_passed_mutations.txt"
        python /usr/local/bin/count_variants.py "${pairName}.pon.blacklist.txt" "${pairName}.count_rejected_mutations.txt"

        if [ "${input_public_release}" = true ] ; 
        then
            python /usr/local/bin/remove_columns.py \
            --IN_MAF "${pairName}.pon_annotated.pass.maf" \
            --OUT_MAF "${pairName}.pon_annotated.pass.public.maf"
            python /usr/local/bin/remove_columns.py \
            --IN_MAF "${pairName}.pon_annotated.maf" \
            --OUT_MAF "${pairName}.pon_annotated.public.maf"
            mv "${pairName}.pon_annotated.pass.public.maf" "${pairName}.pon_annotated.pass.maf"
            mv "${pairName}.pon_annotated.public.maf" "${pairName}.pon_annotated.maf"
        fi

python <<CODE
# Get the header from the passed maf (header is the same in the passed maf and the full maf)
header = open("${pairName}.pon_annotated.pass.maf","r").readline()
header = header.split("\t")
new_header = []
for column in header:
    if column.startswith("pon"):
        new_header.append("${stub}_" + column)
    else:
        new_header.append(column)
header = "\t".join(new_header)
with open("${pairName}.pon_annotated.renamed.pass.maf", "w") as f:
    f.write(header)

with open("${pairName}.pon_annotated.renamed.maf", "w") as f:
    f.write(header)
CODE

        tail -n +2 "${pairName}.pon_annotated.pass.maf" | cat >> "${pairName}.pon_annotated.renamed.pass.maf"
        mv "${pairName}.pon_annotated.renamed.pass.maf" "${pairName}.pon_annotated.pass.maf"
        tail -n +2 "${pairName}.pon_annotated.maf" | cat >> "${pairName}.pon_annotated.renamed.maf"
        mv "${pairName}.pon_annotated.renamed.maf" "${pairName}.pon_annotated.maf"

        mv "${pairName}.pon_annotated.pass.maf" "${pairName}.${stub}.pon_annotated.pass.maf"
        mv "${pairName}.pon_annotated.maf" "${pairName}.${stub}.pon_annotated.maf"

        python /usr/local/bin/add_judgement_column.py \
        --input "${pairName}.${stub}.pon_annotated.maf" \
        --output "${pairName}.${stub}.pon_annotated.with_annotations.maf" \
        --column "${stub}_pon_pass_loglike" \
        --pass_flag "1"
    >>>

    runtime {
        docker: "gcr.io/broad-getzlab-workflows/cga_production_pipeline:v0.2.ccle"
        bootDiskSizeGb: if diskGB_boot != "" then diskGB_boot else default_diskGB_boot
        preemptible: if preemptible != "" then preemptible else default_preemptible
        cpu: if cpu != "" then cpu else default_cpu
        disks: "local-disk ${diskGB} HDD"
        memory: machine_memoryGB + "GB"
    }

    output {
        Int num_rejected_mutations=read_int("${pairName}.count_rejected_mutations.txt")
        Int num_passed_mutations=read_int("${pairName}.count_passed_mutations.txt")
        File allMaf="${pairName}.${stub}.pon_annotated.maf"
        File passMaf="${pairName}.${stub}.pon_annotated.pass.maf"
        File allMaf_with_annotations="${pairName}.${stub}.pon_annotated.with_annotations.maf"
    }
}


task blat {

    # TASK INPUT PARAMS
    File tumorBam
    File tumorBamIdx
    File MAF
    File hg19_bit
    String pairName

    # FILE SIZE
    Int tumorBam_size

    # RUNTIME INPUT PARAMS
    String preemptible
    String diskGB_boot
    String diskGB_buffer
    String memoryGB
    String cpu

    # DEFAULT VALUES
    String default_cpu = "1"
    String default_memoryGB = "10"
    String default_preemptible = "1"
    String default_diskGB_boot = "15"
    String default_diskGB_buffer = "20"

    # COMPUTE MEMORY SIZE
    Int machine_memoryGB = if memoryGB != "" then memoryGB else default_memoryGB

    # COMPUTE DISK SIZE
    Int machine_diskGB_buffer = if diskGB_buffer != "" then diskGB_buffer else default_diskGB_buffer
    Int diskGB = ceil(tumorBam_size + size(MAF, "G") + size(hg19_bit, "G") + machine_diskGB_buffer)

    parameter_meta {
        tumorBam : "sample tumor BAM file"
        tumorBamIdx : "sample tumor BAI file (indexed BAM)"
        MAF : "filename pointing to a mutation annotation format (MAF) file (data for somatic point mutations)"
        pairName : "tumor sample name, string prefix of the output"
    }

    command <<<

        set -euxo pipefail        

        cp -v ${hg19_bit} /opt/hg19.2bit
        python /opt/realign.py ${tumorBam} ${MAF} ${pairName}

        # Count number of passed and rejected mutations
        python /usr/local/bin/count_variants.py "${pairName}.blat.maf" "${pairName}.count_passed_mutations.txt"
        python /usr/local/bin/count_variants.py "${pairName}.blat.rejected.maf" "${pairName}.count_rejected_mutations.txt"

        python /usr/local/bin/add_judgement_column.py \
        --input "${pairName}.blat.all.maf" \
        --output "${pairName}.blat.all.with_annotations.maf" \
        --column "realign_judgment" \
        --pass_flag "KEEP"

    >>>

    runtime {
        docker: "gcr.io/broad-getzlab-workflows/cga_production_pipeline:v0.2.ccle"
        bootDiskSizeGb: if diskGB_boot != "" then diskGB_boot else default_diskGB_boot
        preemptible: if preemptible != "" then preemptible else default_preemptible
        cpu: if cpu != "" then cpu else default_cpu
        disks: "local-disk ${diskGB} HDD"
        memory: machine_memoryGB + "GB"
    }

    output {
        Int num_rejected_mutations=read_int("${pairName}.count_rejected_mutations.txt")
        Int num_passed_mutations=read_int("${pairName}.count_passed_mutations.txt")
        File blat_results="${pairName}.blat.maf"
        File debug_results="${pairName}.blat.rejected.maf"
        File blat_all_maf="${pairName}.blat.all.maf"
        File blat_all_maf_with_annotations="${pairName}.blat.all.with_annotations.maf"
    }
}

task merge_mafs_task {

    # TASK INPUT PARAMS
    File oxoGOBF_maf
    File ffpeOBF_maf
    Array[File] pon_filtered_mafs   
    File blat_maf
    String pairName

    # RUNTIME INPUT PARAMS
    String preemptible
    String diskGB_boot
    String diskGB_buffer
    String memoryGB
    String cpu

    # DEFAULT VALUES
    String default_cpu = "1"
    String default_memoryGB = "10"
    String default_preemptible = "1"
    String default_diskGB_boot = "15"
    String default_diskGB_buffer = "20"

    # COMPUTE MEMORY SIZE
    Int machine_memoryGB = if memoryGB != "" then memoryGB else default_memoryGB

    # COMPUTE DISK SIZE
    Int machine_diskGB_buffer = if diskGB_buffer != "" then diskGB_buffer else default_diskGB_buffer
    Int diskGB = ceil(size(oxoGOBF_maf, "G") + size(ffpeOBF_maf, "G") + size(blat_maf, "G") + machine_diskGB_buffer)

    parameter_meta {
        oxoGOBF_maf       : "MAF with annotations of passing or failing OxoG Orientation Bias filter"
        ffpeOBF_maf       : "MAF with annotations of passing or failing FFPE Orientation Bias filter"        
        blat_maf          : "MAF with annotations of passing or failing BLAT Re-aligner filter"
        pon_filtered_mafs : "MAF(s) with annotations of passing or failing MAFPONFilter"
    }

    command <<<

        set -euxo pipefail

python <<CODE

pon_filtered_mafs_array = '${sep="," pon_filtered_mafs}'.split(",")

command_string_file = open('command_string.txt', 'w')
command_string = ""
for filename in pon_filtered_mafs_array: 
    pon_name = filename.split(".")[1]
    command_string = command_string + "--{0}={1} ".format(pon_name, filename)

command_string_file.write(command_string)
command_string_file.close()
CODE

        COMMAND_STRING=`cat command_string.txt`

        python /usr/local/bin/merge_mafs_after_filters.py \
        --pair_name=${pairName} \
        --oxogOBF=${oxoGOBF_maf} \
        --ffpeOBF=${ffpeOBF_maf} \
        --blat=${blat_maf} \
        $COMMAND_STRING

        python /usr/local/bin/check_for_snvs.py \
        --IN_MAF="${pairName}.passed_all_filters.maf" \
        --OUT_FILE="${pairName}.found_snv"

    >>>

    runtime {
        docker         : "gcr.io/broad-getzlab-workflows/cga_production_pipeline:v0.2.ccle"
        bootDiskSizeGb : if diskGB_boot != "" then diskGB_boot else default_diskGB_boot
        preemptible    : if preemptible != "" then preemptible else default_preemptible
        cpu            : if cpu != "" then cpu else default_cpu
        disks          : "local-disk ${diskGB} HDD"
        memory         : machine_memoryGB + "GB"
    }

    output {
        File merged_intersection_maf="${pairName}.passed_all_filters.maf"
        File merged_union_maf="${pairName}.merged_union.maf"
        Boolean found_snv=read_boolean("${pairName}.found_snv")
    }
}


task mutation_validator {

    # TASK INPUT PARAMS
    String pairName
    File MAF
    File tumorBam
    File tumorBamIdx
    File normalBam
    File normalBamIdx
    File? tumorRNABam
    File? tumorRNABamIdx

    String maf_type

    # FILE SIZE
    Int tumorBam_size
    Int normalBam_size

    # RUNTIME INPUT PARAMS
    String preemptible
    String diskGB_boot
    String diskGB_buffer
    String memoryGB
    String cpu

    # DEFAULT VALUES
    String default_cpu = "1"
    String default_memoryGB = "10"
    String default_preemptible = "1"
    String default_diskGB_boot = "15"
    String default_diskGB_buffer = "20"
    String default_maf_type = "wex"

    # COMPUTE MEMORY SIZE
    Int machine_memoryGB = if memoryGB != "" then memoryGB else default_memoryGB

    # COMPUTE DISK SIZE
    Int machine_diskGB_buffer = if diskGB_buffer != "" then diskGB_buffer else default_diskGB_buffer
    Int diskGB = ceil(tumorBam_size + normalBam_size + size(MAF, "G") 
                        + size(tumorRNABam, "G") + size(tumorRNABamIdx, "G")
                        + machine_diskGB_buffer)

    String input_maf_type = if maf_type != "" then maf_type else default_maf_type

    parameter_meta {
        tumorBam : "sample tumor BAM file"
        tumorBamIdx : "sample tumor BAI file (indexed BAM file)"
        normalBam : "sample normal BAM file"
        normalBamIdx : "sample normal BAI file (indexed BAM file)"
        pairName : "a string for the name of the pair under analysis used for naming output files"

    }

    command <<<

        set -euxo pipefail

        SNV_MAF="${pairName}.snv.maf"
        INDEL_MAF="${pairName}.indel.maf"
        PREPROCESSED_FILE="${pairName}.pileup_preprocessing.txt"
        VALIDATED_SNV_MAF="${pairName}.snp.validated.maf"
        VALIDATED_INDEL_MAF="${pairName}.indel.validated.maf"
        VALIDATED_MAF="${pairName}.validated.maf"

        RNATYPE="hg19-chr"

        # split MAF into INDELs and SNPs
        python /usr/local/bin/split_maf_indel_snp.py -i ${MAF} -o $SNV_MAF -f Variant_Type -v "SNP|DNP|TNP|MNP"
        python /usr/local/bin/split_maf_indel_snp.py -i ${MAF} -o $INDEL_MAF -f Variant_Type -v "INS|DEL"

        python /opt/src/fh_MutationValidatorPreprocess/validation_wrapper_firehose_library_hack.py \
        --mafsnp $SNV_MAF \
        --mafindel $INDEL_MAF \
        --wextumor ${tumorBam} \
        --wexnormal ${normalBam} \
        --rnatumor ${default="None" tumorRNABam} \
        --rnatype $RNATYPE \
        --out ${pairName}

        bash -c "source /matlab_source_file_2012a.sh && /opt/src/fh_MutationValidator/mutation_validator_wrapper $PREPROCESSED_FILE $SNV_MAF ${input_maf_type} ${pairName}.snp"
        bash -c "source /matlab_source_file_2012a.sh && /opt/src/fh_MutationValidator/mutation_validator_wrapper $PREPROCESSED_FILE $INDEL_MAF ${input_maf_type} ${pairName}.indel"
        
        python /usr/local/bin/tsvConcatFiles.py $VALIDATED_SNV_MAF $VALIDATED_INDEL_MAF \
        --outputFilename=$VALIDATED_MAF

    >>>

    runtime {
        docker: "gcr.io/broad-getzlab-workflows/cga_production_pipeline:v0.2.ccle"
        bootDiskSizeGb: if diskGB_boot != "" then diskGB_boot else default_diskGB_boot
        preemptible: if preemptible != "" then preemptible else default_preemptible
        cpu: if cpu != "" then cpu else default_cpu
        disks: "local-disk ${diskGB} HDD"
        memory: machine_memoryGB + "GB"
    }

    output {
        #File pileup_preprocessing_txt="${pairName}.pileup_preprocessing.txt"
        File validated_maf="${pairName}.validated.maf"
    }
}


task gatk_acnv_only{

    # TASK INPUT PARAMS
    File tumorBam
    File tumorBamIdx
    File normalBam
    File normalBamIdx
    String pairName
    File refFasta
    File refFastaDict
    File refFastaIdx
    File common_snp_list
    File PoN
    File gatk_protected    
    #File het_pull_down_from_call_stats_file_tumor_only_diff
    File one_thousand_genomes_common_variants_minor_allele_freq_five
    File mutect1_call_stats

    # FILE SIZE
    Int tumorBam_size
    Int normalBam_size
    Int refFasta_size

    # RUNTIME INPUT PARAMS
    String preemptible
    String diskGB_boot
    String diskGB_buffer
    String memoryGB
    String cpu

    # OPTIONAL PARAMS DEFAULT VALUES
    String default_cpu = "1"
    String default_memoryGB = "10"
    String default_preemptible = "1"
    String default_diskGB_boot = "15"
    String default_diskGB_buffer = "20"

    # COMPUTE MEMORY SIZE
    Int machine_memoryGB = if memoryGB != "" then memoryGB else default_memoryGB
    Int command_memoryGB = machine_memoryGB - 1

    # COMPUTE DISK SIZE
    Int machine_diskGB_buffer = if diskGB_buffer != "" then diskGB_buffer else default_diskGB_buffer
    Int diskGB = ceil(tumorBam_size + normalBam_size + refFasta_size + size(gatk_protected, "G") 
                    + size(PoN, "G") + size(common_snp_list, "G") + size(mutect1_call_stats, "G")
                    + machine_diskGB_buffer)

    parameter_meta {
        tumorBam: "sample tumor BAM file "
        tumorBamIdx : "sample tumor BAI file (indexed BAM)"
        normalBam : "sample normal BAM file"
        normalBamIdx: "sample normal BAI file (indexed BAM)"
        pairName: "a string for the name of the pair under analysis used for naming output files"
        refFasta : "FASTA file for reference genome"
        refFastaIdx : "FASTA file index for the reference genome"
        refFastaDict : "FASTA file dictionary for the reference genome"
        recapseg_bed : "a BED of genomic targets indicating loci for analysis"
        common_snp_list: "a picard interval-list file indicating common heterozygous sites"
        PoN: "a panel-of-normals file; generated for example perhaps by this workflow http://gatkforums.broadinstitute.org/gatk/discussion/comment/31332/"
        padTargets: "a boolean flag indicating whether PadTargets is invoked on the targets before analysis"
    }

    command <<<

        set -euxo pipefail
         
        EFFECTIVE_TARGETS="padded_bed.bed"
        python /usr/local/bin/get_targets.py ${PoN} $EFFECTIVE_TARGETS

        #Coverage Collection
        /usr/local/jre1.8.0_73/bin/java "-Xmx${command_memoryGB}g" -Djava.library.path=/usr/lib/jni/ \
        -jar ${gatk_protected} CalculateTargetCoverage \
        --output ${pairName}.coverage.tsv \
        --groupBy SAMPLE \
        --transform PCOV \
        --targets $EFFECTIVE_TARGETS \
        --targetInformationColumns FULL \
        --input ${tumorBam} \
        --reference ${refFasta} \
        --interval_set_rule UNION \
        --interval_padding 0 \
        --secondsBetweenProgressUpdates 10.0 \
        --disableSequenceDictionaryValidation false \
        --createOutputBamIndex true \
        --help false \
        --version false \
        --verbosity INFO \
        --QUIET false

        #Normalization
        /usr/local/jre1.8.0_73/bin/java "-Xmx${command_memoryGB}g" -Djava.library.path=/usr/lib/jni/ \
        -jar ${gatk_protected} NormalizeSomaticReadCounts \
        --input ${pairName}.coverage.tsv \
        --targets $EFFECTIVE_TARGETS \
        --panelOfNormals ${PoN} \
        --factorNormalizedOutput ${pairName}.fnt.tsv \
        --tangentNormalized ${pairName}.tn.tsv \
        --betaHatsOutput  ${pairName}.betaHats.tsv \
        --preTangentNormalized ${pairName}.preTN.tsv \
        --help false \
        --version false \
        --verbosity INFO \
        --QUIET false
      
        # Segmentation by Tangent-normalized Coverage
        /usr/local/jre1.8.0_73/bin/java "-Xmx${command_memoryGB}g" -Djava.library.path=/usr/lib/jni/ \
        -jar ${gatk_protected} PerformSegmentation \
        --tangentNormalized ${pairName}.tn.tsv \
        --output ${pairName}.seg \
        --log2Input true \
        --alpha 0.01 \
        --nperm 10000 \
        --pmethod HYBRID \
        --minWidth 2 \
        --kmax 25 \
        --nmin 200 \
        --eta 0.05 \
        --trim 0.025 \
        --undoSplits NONE \
        --undoPrune 0.05 \
        --undoSD 3 \
        --help false \
        --version false \
        --verbosity INFO \
        --QUIET false

        /usr/local/jre1.8.0_73/bin/java "-Xmx${command_memoryGB}g" -Djava.library.path=/usr/lib/jni/ \
        -jar ${gatk_protected} GetHetCoverage \
        --reference ${refFasta} \
        --normal ${normalBam} \
        --tumor ${tumorBam} \
        --snpIntervals ${common_snp_list} \
        --normalHets ${pairName}.normal.hets.tsv \
        --tumorHets ${pairName}.tumor.hets.tsv \
        --pvalueThreshold 0.05 \
        --help false \
        --version false \
        --verbosity INFO \
        --QUIET false \
        --VALIDATION_STRINGENCY LENIENT

        mkdir -v plotting
        echo "plotting segmented copy ratio"

        /usr/local/jre1.8.0_73/bin/java "-Xmx${command_memoryGB}g" -Djava.library.path=/usr/lib/jni/ \
        -jar ${gatk_protected} PlotSegmentedCopyRatio \
        --tangentNormalized ${pairName}.tn.tsv \
        --preTangentNormalized ${pairName}.preTN.tsv \
        --sequenceDictionaryFile ${refFastaDict} \
        --outputPrefix ${pairName} \
        --segments  ${pairName}.seg \
        --output plotting/ \
        --log2Input true \
        --help false \
        --version false \
        --verbosity INFO \
        --QUIET false

        echo "allelic cnv"

# Remove header from 
python <<CODE
with open('${pairName}.tn.no_header.tsv', 'w') as writer:
    with open('${pairName}.tn.tsv', 'rb') as reader:
        for line in reader:
            if not line.startswith('#'):
                line = line.strip('\n').split('\t')
                writer.write('\t'.join([line[3]]+line[:3]+[line[-1]])+'\n')
CODE

# Insert Reference Allele and Tumor_Seq_Allele1
python <<CODE
header = ['Chromosome', 'Start_position', 'Reference_Allele', 'Tumor_Seq_Allele1', 'i_t_ref_count', 'i_t_alt_count']
with open("${pairName}.tumor.hets.reformatted.tsv", 'w') as writer:
    with open("${pairName}.tumor.hets.tsv", 'rb') as reader:
        writer.write('\t'.join(header) + '\n')
        for line in reader:
            if not line.startswith('CONTIG'):
                line = line.strip('\n').split('\t')
                writer.write('\t'.join(line[:2]+['A']+['T']+line[2:])+'\n')
CODE

# Re-scale segment mean in seg file
python <<CODE
import math
header = ['Sample', 'Chromosome', 'Start', 'End', 'Num_Probes', 'Segment_Mean']
with open("${pairName}.seg.re_scaled.tsv", 'w') as writer:
    with open("${pairName}.seg", 'rb') as reader:
        writer.write('\t'.join(header) + '\n')
        for line in reader:
            if not line.startswith('Sample'):
                line = line.strip('\n').split('\t')
                segment_mean = math.log(float(line[-1]))/math.log(2)
                writer.write('\t'.join(line[:-1]+[str(segment_mean)])+'\n')
CODE

        python /usr/local/bin/het_pull_down_from_call_stats_file_tumor_only_diff.py \
        --input-call-stats-filename ${mutect1_call_stats} \
        --db-variant-key-filename ${one_thousand_genomes_common_variants_minor_allele_freq_five} \
        --no-normal true \
        --output-filename "${pairName}.tumor.hets.reformatted.tsv"

        Rscript /opt/Allelic_CapSeg/AllelicCapseg_cli.R --SID=${pairName} \
        --capseg.probe.fn="${pairName}.tn.no_header.tsv" \
        --capseg.seg.fn="${pairName}.seg.re_scaled.tsv" \
        --germline.het.fn="${pairName}.tumor.hets.reformatted.tsv" \
        --drop.x=FALSE \
        --drop.y=TRUE \
        --seg.merge.thresh=0.5 \
        --min.seg.size=3 \
        --verbose=TRUE \
        --base.output.dir=. \
        --initial.merge=TRUE \
        --split.merge=TRUE \
        --outlier.thresh=0.005 \
        --working.dir=/opt/Allelic_CapSeg/

    >>>

    runtime {
        docker: "gcr.io/broad-getzlab-workflows/cga_production_pipeline:v0.2.ccle"
        bootDiskSizeGb: if diskGB_boot != "" then diskGB_boot else default_diskGB_boot
        preemptible: if preemptible != "" then preemptible else default_preemptible
        cpu: if cpu != "" then cpu else default_cpu
        disks: "local-disk ${diskGB} HDD"
        memory: machine_memoryGB + "GB"
    }

    output {
        File gatk_cnv_coverage_file="${pairName}.coverage.tsv"
        File gatk_cnv_seg_file="${pairName}.seg"
        File gatk_cnv_tn_coverage="${pairName}.tn.tsv"
        File gatk_cnv_pre_tn_coverage="${pairName}.preTN.tsv"
        File gatk_het_ad_normal="${pairName}.normal.hets.tsv"
        File gatk_het_ad_tumor="${pairName}.tumor.hets.tsv"
        Array[File] gatk_cnv_all_plots=glob("plotting/*.png")
        File alleliccapseg_plot="plots/${pairName}/${pairName}_FullGenome.png"
        File alleliccapseg_tsv="results/${pairName}.tsv"
        String alleliccapseg_skew=read_string("results/${pairName}.skew")
        File hets_low_1000g_vaf_added_in="${pairName}.tumor.hets.reformatted.tsv"
    }
}


task absolute {
    
    # TASK INPUT PARAMS
    File seg_file
    File maf
    String skew
    String pairName

    # RUNTIME INPUT PARAMS
    String preemptible
    String diskGB_boot
    String diskGB_buffer
    String memoryGB
    String cpu

    # OPTIONAL PARAMS DEFAULT VALUES
    String default_cpu = "1"
    String default_memoryGB = "7"
    String default_preemptible = "1"
    String default_diskGB_boot = "15"
    String default_diskGB_buffer = "20"

    # COMPUTE MEMORY SIZE
    Int machine_memoryGB = if memoryGB != "" then memoryGB else default_memoryGB

    # COMPUTE DISK SIZE
    Int machine_diskGB_buffer = if diskGB_buffer != "" then diskGB_buffer else default_diskGB_buffer
    Int diskGB = ceil(size(seg_file, "G") + size(maf, "G") + machine_diskGB_buffer)

    parameter_meta {
        seg_file : "filename pointing to the input - either a HAPSEG file or a segmentation file"
        maf : "filename pointing to a mutation annotation format (MAF) file. This specifies the data for somatic point mutations to be used by ABSOLUTE."
        skew : ""
        pairName : "a string for the name of the pair under analysis used for naming output files"
    }

    command <<<

        set -euxo pipefail

        SNV_MAF="${pairName}.snv.maf"
        INDEL_MAF="${pairName}.indel.maf"
        python /usr/local/bin/split_maf_indel_snp.py -i ${maf} -o $SNV_MAF -f Variant_Type -v "SNP|DNP|TNP|MNP"
        python /usr/local/bin/split_maf_indel_snp.py -i ${maf} -o $INDEL_MAF -f Variant_Type -v "INS|DEL"

        # TODO: replace with python script
        grep -v "NA" ${seg_file} > no_nan_segs.tsv

        Rscript /xchip/tcga/Tools/absolute/releases/v1.5/run/ABSOLUTE_cli_start.R \
        --seg_dat_fn no_nan_segs.tsv \
        --maf_fn $SNV_MAF \
        --indelmaf_fn $INDEL_MAF \
        --sample_name ${pairName} \
        --results_dir . \
        --ssnv_skew ${skew} \
        --abs_lib_dir /xchip/tcga/Tools/absolute/releases/v1.5/

    >>>

    runtime {
        docker         : "gcr.io/broad-getzlab-workflows/cga_production_pipeline:v0.2.ccle"
        bootDiskSizeGb : if diskGB_boot != "" then diskGB_boot else default_diskGB_boot
        preemptible    : if preemptible != "" then preemptible else default_preemptible
        cpu            : if cpu != "" then cpu else default_cpu
        disks          : "local-disk ${diskGB} HDD"
        memory         : machine_memoryGB + "GB"
    }

    output {
        # Plot showing the Purity/Ploidy values and the solutions
        File absolute_highres_plot="${pairName}.ABSOLUTE_plot.pdf"
        # An R file containing an object seg.dat which provides all of the information used to generate the plot.
        File absolute_rdata="${pairName}.PP-modes.data.RData"
    }

}


task lego_plotter_task {

    # TASK INPUT PARAMS
    File maf
    File mut_categs
    String pairName
    String covString

    # RUNTIME INPUT PARAMS
    String preemptible
    String diskGB_boot
    String diskGB_buffer
    String memoryGB
    String cpu

    # OPTIONAL PARAMS DEFAULT VALUES
    String default_cpu = "1"
    String default_memoryGB = "10"
    String default_preemptible = "1"
    String default_diskGB_boot = "15"
    String default_diskGB_buffer = "20"
    String default_covString = "exome"

    # COMPUTE MEMORY SIZE
    Int machine_memoryGB = if memoryGB != "" then memoryGB else default_memoryGB
    Int command_memoryGB = machine_memoryGB - 1

    # COMPUTE DISK SIZE
    Int machine_diskGB_buffer = if diskGB_buffer != "" then diskGB_buffer else default_diskGB_buffer
    Int diskGB = ceil(size(maf, "G") + size(mut_categs, "G") + machine_diskGB_buffer)

    String coverage = if covString != "" then covString else default_covString

    parameter_meta {
        maf : "maf file for plotting"
        mut_categs : "reference data"
        pairName : "name of the pair to use"
        covString : "'exome', 'genome', or 'unit' to pick typical exome coverage, genome coverage, or 'unit' as every base counts equally in the normalization."
    }

    command <<<

        #increase verbosity to inform of commands run
        set -x 

        #create directory where category data is expected
        mkdir -pv /xchip/cga/reference/hg19/annotation/db/tracks/hg19/c65e/

        #copy the category data
        cp -vf ${mut_categs} /xchip/cga/reference/hg19/annotation/db/tracks/hg19/c65e/

        #link for the name
        ln -vs ${maf} ${pairName}.maf

        #run the plotter and acquire the exit code (plot rates and zscale)
        /usr/local/bin/run_call_lego_plotter.sh ${pairName}.maf  . ${coverage} yes yes
        PLOTTER_EXIT_CODE=$? ;

        #encode the PNGS
        for PNG in `find *.png`; do 
            uuencode  -m $PNG  /dev/stdout | grep -Pv '^begin'|grep -Pv '^====$'|tr -d "\n" > $PNG.enc ; 
            UUE=$? ; 
            if [ "$UUE" -ne "0" ] ; then exit 2 ; fi ;
            done 

        #update the HTML in memory/place
python_cmd="
import glob
import re
for H in glob.glob('*_legos.html'):
    reader=open(H,'r')
    #store the HTML into a single string
    all_lines=''
    for line in reader:
        all_lines=all_lines+line
    reader.close()
    #iterate over the encs
    for ENC in glob.glob('*.enc'):
        print 'got enc = ',ENC
        PNG=re.sub('\.enc$','',ENC)
        enc_reader=open(ENC,'r')
        enc_data=''
        for enc_line in enc_reader:
            enc_data=enc_data+enc_line.strip()
        enc_reader.close()
        #embed the img, close the src, and add a title
        enc_replacement='data:image/png;base64,'+str(enc_data)+'\" title=\"'+str(PNG)+'\" '
        all_lines=all_lines.replace(PNG,enc_replacement)
    writer=open(H,'w')
    writer.write(all_lines)
    writer.close()
"
        python -c "$python_cmd" 
        PY_EXIT=$? ;
        if [ "$PY_EXIT" -ne "0" ] ; then exit 3 ; fi ;

        #propagate exit code
        bash -c "exit $PLOTTER_EXIT_CODE"

    >>>

    runtime {
        docker         : "gcr.io/broad-getzlab-workflows/cga_production_pipeline:v0.2.ccle"
        bootDiskSizeGb : if diskGB_boot != "" then diskGB_boot else default_diskGB_boot
        preemptible    : if preemptible != "" then preemptible else default_preemptible
        cpu            : if cpu != "" then cpu else default_cpu
        disks          : "local-disk ${diskGB} HDD"
        memory         : machine_memoryGB + "GB"
    }

    output {
        Array[File] ais=glob("*.ai")
        Array[File] pngs=glob("*.png")
        Array[File] figs=glob("*.fig")
        Array[File] pss=glob("*.ps")
        File mut_legos_html="${pairName}.maf.mutation_legos.html"
    }
}
