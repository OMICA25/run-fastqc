#!/usr/bin/env nextflow

raw_reads = params.rawReads
out_dir  = file(params.outDir)

out_dir.mkdir()

/*
  Single-end mode:
  - read all fastq / fastq.gz files from the input directory
  - create a channel of tuples: (sampleName, file)
  - sampleName uses the file base name (without extension)
*/
reads = Channel.fromPath("${raw_reads}/*.{fastq,fastq.gz}").map { f -> tuple(f.baseName, f) }

process runFastQC{
    tag { "${params.projectName}.rFQC.${sample}" }
    cpus { 2 }
    publishDir "${out_dir}/qc/raw/${sample}", mode: 'copy', overwrite: false

    input:
        tuple val(sample), file(in_fastq) from reads

    output:
        file("${sample}_fastqc/*.zip") into fastqc_files

    """
    mkdir -p ${sample}_fastqc
    fastqc --outdir ${sample}_fastqc \
    -t ${task.cpus} \
    ${in_fastq}
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
