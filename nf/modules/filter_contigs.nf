#!/usr/bin/env nextflow

/*
 * Filter out short contigs
 */
process FILTER_CONTIGS {

    container "docker.io/kameronhn/valet-python:latest"

    input:
    path assembly
    path scriptfile

    output:
    path "filtered.fasta"       , emit:filtered_fasta
    path "contig_lengths.tsv"   , emit:contig_lengths

    script:
    """
    python ${scriptfile} \
        --fasta_filename $assembly \
        --filtered_fasta_filename filtered.fasta \
        --length_filename contig_lengths.tsv
    """
}