workflow CommonVariantFilter {
    String sampleId
    File maf

    Int min_exac_ac = 10
    Int min_filter_depth = 0
    Boolean? disable_whitelist = false
    Boolean? filter_noncoding = false

    Int RAM = 8
    Int SSD = 60
    Int preemptible = 3

    String docker_tag = "1.0.3"

    meta {
        author: "Brendan Reardon"
        email: "breardon@broadinstitute.org"
        laboratory: "Van Allen Lab"
        institution: "Dana-Farber Cancer Institute, Broad Institute of MIT & Harvard"
        github: "https://github.com/vanallenlab/common_variant_filter"
        license: "MIT License"
    }

    call commonfilterTask {
        input: sampleId=sampleId,
            maf=maf,
            min_exac_ac=min_exac_ac,
            min_filter_depth=min_filter_depth,
            disable_whitelist=disable_whitelist,
            filter_noncoding=filter_noncoding,
            RAM=RAM,
            SSD=SSD,
            preemptible=preemptible,
            docker_tag=docker_tag
    }

    output  {
        File commonfilter_annotated_maf = commonfilterTask.annotatedMAF
        File commonfilter_passed_maf = commonfilterTask.passedMAF
        File commonfilter_rejected_maf = commonfilterTask.rejectedMAF
        Int commonfilter_considered_count = commonfilterTask.consideredCount
        Int commonfilter_pass_count = commonfilterTask.passCount
        Int commonfilter_reject_count = commonfilterTask.rejectCount
    }
}

task commonfilterTask {
    String sampleId
    File maf

    Int min_exac_ac
    Int min_filter_depth
    Boolean? disable_whitelist
    Boolean? filter_noncoding

    Int? RAM
    Int? SSD
    Int? preemptible

    String? docker_tag

    command {
        args="--min_exac_ac "${min_exac_ac}" "
        args+="--min_filter_depth "${min_filter_depth}" "
        args+=${true="--filter_noncoding" false="" filter_noncoding}" "
        args+=${true="--disable_wl" false="" disable_whitelist}" "

        sed '/^#/ d' < ${maf} > uncommented_maf.maf
        python /common_variant_filter.py --id ${sampleId} --maf uncommented_maf.maf $args
    }

    output  {
        File annotatedMAF="${sampleId}.common_variant_filter.annotated.maf"
        File passedMAF="${sampleId}.common_variant_filter.pass.maf"
        File rejectedMAF="${sampleId}.common_variant_filter.reject.maf"
        String consideredCount=read_string("considered.txt")
		String passCount=read_string("passed.txt")
		String rejectCount=read_string("rejected.txt")
    }

    runtime {
        disks: "local-disk " + SSD + " SSD"
        docker: "vanallenlab/common_variant_filter:" + docker_tag
        memory: RAM + " GB"
        preemptible: preemptible
    }
}
