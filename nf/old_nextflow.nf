#!/usr/bin/env nextflow

// Module INCLUDE statements
include { FILTER_CONTIGS } from './modules/filter_contigs.nf'
include { RUN_BOWTIE } from './modules/bowtie.nf'
include { CONTIG_LENGTHS } from './modules/contig_lengths.nf'
include { RUN_SAMTOOLS } from './modules/samtools.nf'
include { CONTIG_COVERAGE } from './modules/contig_coverage.nf'

/*
 * Pipeline parameters
 */

// Primary input
params {
    assemblies: String
    reads1: Path
    reads2: Path
}

workflow {

    

    main:

    assemblies_ch = channel.from(
        params.assemblies
            .split(',')
            .collect { f -> file(f.trim()) }
    ).map { fasta ->
        tuple(fasta.baseName,fasta)
    }

    // Call processes

    FILTER_CONTIGS(assemblies_ch)

    RUN_BOWTIE(assemblies_ch, FILTER_CONTIGS.out.filtered_fasta, params.reads1, params.reads2)

    CONTIG_LENGTHS(assemblies_ch, RUN_BOWTIE.out.sam)

    RUN_SAMTOOLS(assemblies_ch, RUN_BOWTIE.out.sam)

    CONTIG_COVERAGE(assemblies_ch, RUN_SAMTOOLS.out.pileup)

    publish:
    filtered_fasta = FILTER_CONTIGS.out.filtered_fasta
    index = RUN_BOWTIE.out.index
    sam = RUN_BOWTIE.out.sam
    bam = RUN_SAMTOOLS.out.bam
}

output {
    filtered_fasta {
        path 'filtered_fasta'
    }
    index {
        path 'index'
    }
    sam {
        path 'sam'
    }
    bam {
        path 'bam'
    }
}
