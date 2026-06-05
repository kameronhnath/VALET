#!/usr/bin/env nextflow

/*
 * Processes for the mate pair process
 */
process RUN_BWA {

    container 'quay.io/biocontainers/bwa:0.7.19--h577a1d6_1'

    input:
    tuple val(bin_id), path(bin_path)

    output:
    tuple val(bin_id), path("align.sam")

    script:
    """
    bwa index ${bin_path}/contigs.fasta
    bwa mem ${bin_path}/contigs.fasta ${bin_path}/lib_1.fq ${bin_path}/lib_2.fq > align.sam
    """
}

process SAM_TO_BAM {

    container 'biocontainers/samtools:v1.9-4-deb_cv1'

    input:
    tuple val(bin_id), path(samfile)

    output:
    tuple val(bin_id), path("align.bam")

    script:
    """
    samtools sort -o align.bam ${samfile}
    """
}

process MATE_PAIR_CHECKER {

    // https://hub.docker.com/r/biocontainers/python3-pysam/tags

    input:
    path abc

    output:
    path abc

    script:
    """
    """

}