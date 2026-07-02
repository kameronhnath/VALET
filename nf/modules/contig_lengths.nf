#!/usr/bin/env nextflow

/*
 * Dict of contig lengths.
 */
process CONTIG_LENGTHS {

    container "docker.io/kameronhn/valet-python:latest"

    input:
    path samfile
    path scriptfile

    output:
    path "contig.lengths"   , emit:contig_lengths

    script:
    """
    python ${scriptfile} \
        --sam_filename ${samfile} \
        --out_filename contig.lengths
    """
}