#!/usr/bin/env nextflow

/*
 * Artificially split up the unaligned reads.
 */
process BREAKPOINT_SPLITTER {

    input:
    path unaligned

    output:
    path "breakpoint_reads",    emit:breakpoint_reads

    script:
    """
    mkdir breakpoint_reads
    python ${projectDir}/../src/py/nf/breakpoint_splitter.py \
        --unaligned ${unaligned} \
        --output breakpoint_reads/
    """
}

process BREAKPOINT_FINDER {

    container 'quay.io/biocontainers/mulled-v2-da271c8774be2b8dbd760259a085347c47897e8b:3eba2b66b7b62c358a394d58655730a65ba3b4c8-0'

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