#!/usr/bin/env nextflow

/*
 * Calculate contig coverage.
 */
process CONTIG_COVERAGE {

    input:
    path assembly
    path pileup

    output:
    path "coverage"   , emit:coverage

    script:
    """
    python ${projectDir}/../src/py/calculate_contig_coverage.py \
        --pileup_file ${pileup} \
        --out_filename coverage
    """
}