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

    container 'quay.io/biocontainers/mulled-v2-480c331443a1d7f4cb82aa41315ac8ea4c9c0b45:3e0fc1ebdf2007459f18c33c65d38d2b031b0052-0'

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