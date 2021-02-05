task aggregateMAFs_selectFields {


    Array[File] inMafFNs
    String outputFN_prfx
    String SelectFields
    File aggregate_selectedFields_MAF_Script
    
    Int memory
    Int disk_space
    Int num_preempt

    
    command {
        Rscript ${aggregate_selectedFields_MAF_Script} "${outputFN_prfx}.mergedMAF.txt" ${write_lines(inMafFNs)} ${SelectFields}
    }

    output {
        File mergedMAFfn="${outputFN_prfx}.mergedMAF.txt"
    }

    runtime {
        docker: "gcr.io/ccle-docker-images/qwbroad/ccle-firecloud:1"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        preemptible: "${num_preempt}"
    }

    meta {
        author: "Mahmoud Ghandi"
    }
}


workflow aggregateMAFs_selectFields_workflow {
    call aggregateMAFs_selectFields
}
