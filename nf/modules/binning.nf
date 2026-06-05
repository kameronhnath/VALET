#!/usr/bin/env nextflow

/*
 * Bin reads and contigs by coverage
 */
process BIN_READS_AND_CONTIGS {

    input:
    path samfile
    path coverage_file
    path assembly
    val threads
    val min_coverage

    output:
    path "out/bins",    emit:bins
    path "out/bins/*",  emit:bin_dirs

    script:
    """
    python ${projectDir}/../src/py/nf/bin_reads.py \
        --samfile ${samfile} \
        --coverage_file ${coverage_file} \
        --threads ${threads} \
        --min_coverage ${min_coverage}

    python ${projectDir}/../src/py/nf/bin_contigs.py \
        --assembly ${assembly} \
        --coverage_file ${coverage_file} \
        --threads ${threads} \
        --min_coverage ${min_coverage}
    """
}