#!/usr/bin/env nextflow

/*
 * Dict of contig lengths.
 */
process CONTIG_LENGTHS {

    input:
    path samfile

    output:
    path "contig.lengths"   , emit:contig_lengths

    script:
    """
    python ${projectDir}/../src/py/get_contig_lengths.py \
        --sam_filename ${samfile} \
        --out_filename contig.lengths
    """
}