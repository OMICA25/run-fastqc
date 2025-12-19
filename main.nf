#!/usr/bin/env nextflow

raw_reads = params.rawReads
out_dir = file(params.outDir)

out_dir.mkdir()

read_pair = Channel.fromFilePairs("${raw_reads}/*R[1,2].fastq", type: 'file')

process runFastQC{
    tag { "${params.projectName}.rFQC.${sample}" }
    cpus { 2 }
    publishDir "${out_dir}/qc/raw/${sample}", mode: 'copy', overwrite: false

    input:
        set sample, file(in_fastq) from read_pair

    output:
        file("${sample}_fastqc/*.zip") into fastqc_files

    """
    mkdir ${sample}_fastqc
    fastqc --outdir ${sample}_fastqc \
    -t ${task.cpus} \
    ${in_fastq.get(0)} \
    ${in_fastq.get(1)}
    """
}

process runMultiQC{
    tag { "${params.projectName}.rMQC" }
    publishDir "${out_dir}/qc/raw", mode: 'copy', overwrite: false

    input:
        file('*') from fastqc_files.collect()

    output:
        file('multiqc_report.html')

    """
    multiqc .
    """
}

workflow.onComplete {

    println ( workflow.success ? """
        Pipeline execution summary
        ---------------------------
        Completed at: ${workflow.complete}
        Duration    : ${workflow.duration}
        Success     : ${workflow.success}
        workDir     : ${workflow.workDir}
        exit status : ${workflow.exitStatus}
        """ : """
        Failed: ${workflow.errorReport}
        exit status : ${workflow.exitStatus}
        """
    )
}
