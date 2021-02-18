workflow filterMAF_for_CGA_workflow {
    call filterMaf {}
}

task filterMaf {
  File inMAFfn
  Float minAF
  Float maxExAC_AF
  Float maxExAC_AF_COSMIC
  Float maxPon_loglike
  Boolean onlyCoding
  Float minCoverage
  File filtermaf_script
  Float minAltReads
  File TCGAhotspotMutFN
  Float TCGAhotspotMinCnt
  Float COSMIChotspotMinCnt
  File blackListFN
  File? intersectWithMafFN
  String sample_id

  File addExAC_column

  command <<<
    set -euo pipefail
    touch ${sample_id}_outMAFfn.txt
    touch ${sample_id}_outMAFannotatedFN.txt

    Rscript ${addExAC_column} ${inMAFfn}

    Rscript ${filtermaf_script} CLEANED.maf ${sample_id}"_outMAFfn.txt" \
    ${sample_id}"_outMAFannotatedFN.txt" ${minAF} ${maxExAC_AF} \
    ${maxExAC_AF_COSMIC} ${maxPon_loglike} ${onlyCoding} \
    ${minCoverage} ${minAltReads} ${TCGAhotspotMutFN} ${TCGAhotspotMinCnt} \
    ${COSMIChotspotMinCnt} ${blackListFN} ${default="NULL" intersectWithMafFN}
  >>>

  output {
    File outMAFfn = "${sample_id}_outMAFfn.txt"
    File outMAFannotatedFN = "${sample_id}_outMAFannotatedFN.txt"
  }

  runtime {
    docker: "mghandi/ccle_docker_r"
    memory: "10GB"
    disks: "local-disk 50 HDD"
    preemptible: "5"
  }
}
