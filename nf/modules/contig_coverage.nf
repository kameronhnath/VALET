#!/usr/bin/env nextflow

/*
 * Calculate contig coverage.
 */
process CONTIG_COVERAGE {

    container "docker.io/kameronhn/valet-python:latest"

    input:
    path assembly
    path pileup
    path scriptfile

    output:
    path "coverage"   , emit:coverage

    script:
    """
    python ${scriptfile} \
        --pileup_file ${pileup} \
        --out_filename coverage
    """
}