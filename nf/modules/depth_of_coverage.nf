#!/usr/bin/env nextflow

/*
 * Determine the depth of coverage
 */
process DEPTH_OF_COVERAGE {

    container "docker.io/kameronhn/valet-python:latest"

    input:
    path chunk
    val window_size
    path scriptfile

    output:
    path "${chunk}_errors_coverage.bed"   , emit:coverage_errrors

    script:
    """
    python ${scriptfile} \
        --mpileup-file ${chunk} \
        -o ${chunk}_errors_coverage.bed \
        -g -e \
        -w ${window_size}
    """
}

/*
 * Merge the bed outputs if split into chunks
 */
process MERGE_BEDS_COVERAGE {

    container 'quay.io/biocontainers/bedtools:2.31.1--hf5e1c6e_0'

    input:
    path beds

    output:
    path "coverage.errors.merged.bed"   , emit:merged_bed

    script:
    """
    cat ${beds.join(' ')} \
        | bedtools sort -i - \
        > coverage.errors.merged.bed
    """
}