#!/usr/bin/env nextflow

process FIND_SUSPICIOUS_REGIONS {

    container 'quay.io/biocontainers/bedtools:2.31.1--hf5e1c6e_0'

    input:
    path bedfile
    val suspicious_flank_size
    val min_suspicious_regions

    output:
    path "suspicious.bed"

    script:
    """
    sort -k1,1 -k2,2n ${bedfile} \
        | cut -f1-4 \
        | bedtools merge \
            -i - \
            -d ${suspicious_flank_size} \
            -c 4 \
            -o distinct \
        | awk -F'\\t' '
            {
                n=split(\$4,a,",")
                if (n >= ${min_suspicious_regions})
                    print
            }
        ' \
        > suspicious.bed
    """
}