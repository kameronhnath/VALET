#!/usr/bin/env nextflow

/*
 * Split the pileup file based on the number of threads (input).
 */
process SPLIT_PILEUP {

    input:
    path pileup
    val threads

    output:
    path "${pileup}.*"   , emit:files

    script:
    """
    python ${projectDir}/../src/py/split_pileup.py \
        --pileup_file ${pileup} \
        --chunks ${threads}
    """
}