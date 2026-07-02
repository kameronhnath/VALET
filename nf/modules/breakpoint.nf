#!/usr/bin/env nextflow

/*
 * Artificially split up the unaligned reads.
 */
process BREAKPOINT_SPLITTER {

    container "docker.io/kameronhn/valet-python:latest"

    input:
    path unaligned
    path script_file

    output:
    path "breakpoint_reads",    emit:breakpoint_reads

    script:
    """
    mkdir breakpoint_reads
    python ${script_file} \
        --unaligned ${unaligned} \
        --output breakpoint_reads/
    """
}

process BREAKPOINT_FINDER {

    container 'docker.io/kameronhn/breakpoint:latest'

    input:
    path script_file
    path assembly
    path reads_dir
    val breakpoint_bins
    val threads
    path coverage_file

    output:
    path "breakpoint/",                         emit:breakpoint_dir
    path "breakpoint/interesting_bins.bed",     emit:breakpoint_bins


    script:
    """
    python ${script_file} \
        --assembly-file ${assembly} \
        --reads-dir ${reads_dir} \
        --bin-size ${breakpoint_bins} \
        --coverage ${coverage_file} \
        --threads ${threads} \
        --output breakpoint/
    """
}

process BREAKPOINT_BED_SORT {

    container 'quay.io/biocontainers/bedtools:2.31.1--hf5e1c6e_0'

    input:
    path bedfile

    output:
    path "breakpoints.bed"

    script:
    """
    bedtools sort -i ${bedfile} > breakpoints.bed
    """
}