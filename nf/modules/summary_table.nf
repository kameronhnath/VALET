#!/usr/bin/env nextflow

/*
 * Generate a summary file and summary table
 */
process GENERATE_SUMMARY {

    input:
    path coverage_bed
    path breakpoint_bed
    path matepair_bed
    path contig_lengths_file
    path contig_coverage
    path filtered_contigs

    output:
    path "summary.bed"      , emit:summary_bed
    path "summary.tsv"      , emit:summary_table

    script:
    """
    python ${projectDir}/../src/py/nf/summary_table.py \
        --coverage_bed ${coverage_bed} \
        --breakpoint_bed ${breakpoint_bed} \
        --matepair_bed ${matepair_bed} \
        --contig_lengths_file ${contig_lengths_file} \
        --summary_file summary.bed \
        --abundance_filename ${contig_coverage} \
        --filtered_contigs ${filtered_contigs} \
        --table_filename summary.tsv
    """
}