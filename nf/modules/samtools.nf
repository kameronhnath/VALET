#!/usr/bin/env nextflow

/*
 * Runs the samtools commands as part of the pipeline
 */
process RUN_SAMTOOLS {

    container 'kameronhn/samtools:latest'

    input:
    path assembly
    path samfile

    output:
    path "library.bam"              , emit: bam
    path "sorted_library.bam"       , emit: sorted_bam
    path "sorted_library.bam.bai"   , emit: bam_index
    path "mpileup_output.out"       , emit: pileup

    script:
    """
    # Convert SAM -> BAM
    samtools view -F 0x100 -bS $samfile > library.bam

    # Sort BAM
    samtools sort -@ 16 -o sorted_library.bam library.bam

    # Index BAM
    samtools index sorted_library.bam

    # Generate pileup
    samtools mpileup -C50 -A -f $assembly sorted_library.bam > mpileup_output.out
    """
}