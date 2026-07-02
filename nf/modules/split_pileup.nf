#!/usr/bin/env nextflow

/*
 * Split the pileup file based on the number of threads (input).
 */
process SPLIT_PILEUP {

    container "docker.io/kameronhn/valet-python:latest"

    input:
    path pileup
    val threads
    path scriptfile

    output:
    path "${pileup}.*"   , emit:files

    script:
    """
    python ${scriptfile} \
        --pileup_file ${pileup} \
        --chunks ${threads}
    """
}