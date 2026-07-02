#!/usr/bin/env nextflow

/*
 * Bin reads and contigs by coverage
 */
process BIN_READS_AND_CONTIGS {

    container "docker.io/kameronhn/valet-python:latest"

    input:
    path samfile
    path coverage_file
    path assembly
    val threads
    val min_coverage
    path reads_scriptfile
    path contigs_scriptfile

    output:
    path "out/bins",    emit:bins
    path "out/bins/*",  emit:bin_dirs

    script:
    """
    python ${reads_scriptfile} \
        --samfile ${samfile} \
        --coverage_file ${coverage_file} \
        --threads ${threads} \
        --min_coverage ${min_coverage}

    python ${contigs_scriptfile} \
        --assembly ${assembly} \
        --coverage_file ${coverage_file} \
        --threads ${threads} \
        --min_coverage ${min_coverage}
    """
}