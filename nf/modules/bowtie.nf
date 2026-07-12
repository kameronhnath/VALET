#!/usr/bin/env nextflow

/*
 * Build index and map reads
 */
process RUN_BOWTIE {

    container 'quay.io/biocontainers/bowtie2:2.5.4--he20e202_2'

    input:
    path assembly
    path filtered_fasta
    path read1
    path read2

    output:
    path "${assembly.baseName}.index*.bt2"       , emit:index
    path "${assembly.baseName}.sam"              , emit:sam
    path "${assembly.baseName}.paired.sam"       , emit:paired_sam
    path "unaligned.fastq"                       , emit:unaligned_reads

    script:
    def name = assembly.baseName
    """
    bowtie2-build ${filtered_fasta} ${name}.index
    bowtie2 -x ${name}.index -U ${read1} ${read2} --reorder -p 4 --un unaligned.fastq -S ${name}.sam
    bowtie2 -x ${name}.index -1 ${read1} -2 ${read2} --reorder -p 4 -S ${name}.paired.sam
    """
}