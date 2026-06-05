#!/usr/bin/env nextflow

/*
 * Filter out short contigs
 */
process FILTER_CONTIGS {

    input:
    path assembly

    output:
    path "filtered.fasta"       , emit:filtered_fasta
    path "contig_lengths.tsv"   , emit:contig_lengths

    script:
    """
    python ${projectDir}/../src/py/filter_short_contigs.py \
        --fasta_filename $assembly \
        --filtered_fasta_filename filtered.fasta \
        --length_filename contig_lengths.tsv
    """
}