#!/usr/bin/env nextflow

/*
 * Align the binned reads (by coverage)
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

/*
 * Convert the alignment file to bam
 */
process SAM_TO_BAM {

    container 'kameronhn/samtools:latest'

    input:
    tuple val(bin_id), path(samfile)

    output:
    tuple val(bin_id), path("align.bam")

    script:
    """
    samtools sort -o align.bam ${samfile}
    """
}

/*
 * Run the mate pair checker
 */
process MATE_PAIR_CHECKER {

    container "docker.io/kameronhn/valet-python:latest"

    input:
    tuple val(bin_id), path(bin_path), path(bamfile)
    path script_file

    output:
    path "${bin_id}.mate_error.bed"

    script:
    def adjusted_support = Math.ceil(bin_id.toInteger() / 1.5) as int
    """
    python3 ${script_file} \
        -1 ${bin_path}/lib_1.fq \
        -2 ${bin_path}/lib_2.fq \
        -bam ${bamfile} \
        -o ${bin_id}.mate_error.bed \
        -s ${adjusted_support}
    """

}

/*
 * Merge the resulting bed files
 */
process MERGE_BEDS_MATE_ERROR {

    container 'quay.io/biocontainers/bedtools:2.31.1--hf5e1c6e_0'

    input:
    path beds

    output:
    path "mate.errors.merged.bed"

    script:
    """
    cat ${beds.join(' ')} \
        | bedtools sort -i - \
        > mate.errors.merged.bed
    """

}