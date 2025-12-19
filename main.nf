#!/usr/bin/env nextflow

raw_reads = params.rawReads
out_dir = file(params.outDir)

out_dir.mkdir()

// Single-end reads: accept .fastq and .fastq.gz files.
// Map each input file to a tuple: [sampleName, file]
read_files = Channel.fromPath("${raw_reads}/*fastq*", checkIfExists: true)

reads_channel = read_files.map { f ->
    // derive a sample name by removing common suffixes and extensions
    def name = f.getName()
        .replaceAll(/(\.fastq|\.fq)(\.gz)?$/, '')   // remove .fastq, .fq and optional .gz
        .replaceAll(/_R?1(_\d+)?$/, '')             // remove _R1, R1_001, _1, etc.
    return [name, f]
}

process runFastQC{
    tag { "${params.projectName}.rFQC.${sample}" }
    cpus { 2 }
    publishDir "${out_dir}/qc/raw/${sample}", mode: 'copy', overwrite: false

    input:
        set sample, file(in_fastq) from reads_channel

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
    publishDir "${out_dir}/qc

