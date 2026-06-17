#!/usr/bin/env nextflow

// Module INCLUDE statements
include { FILTER_CONTIGS } from './modules/filter_contigs.nf'
include { RUN_BOWTIE } from './modules/bowtie.nf'
include { CONTIG_LENGTHS } from './modules/contig_lengths.nf'
include { RUN_SAMTOOLS } from './modules/samtools.nf'
include { CONTIG_COVERAGE } from './modules/contig_coverage.nf'
include { SPLIT_PILEUP } from './modules/split_pileup.nf'
include { DEPTH_OF_COVERAGE ; MERGE_BEDS_COVERAGE } from './modules/depth_of_coverage.nf'
include { BREAKPOINT_SPLITTER ; BREAKPOINT_FINDER ; BREAKPOINT_BED_SORT } from './modules/breakpoint.nf'
include { BIN_READS_AND_CONTIGS } from './modules/binning.nf'
include { RUN_BWA ; SAM_TO_BAM ; MATE_PAIR_CHECKER ; MERGE_BEDS_MATE_ERROR } from './modules/mate_pair.nf'
include { GENERATE_SUMMARY } from './modules/summary_table.nf'
include { FIND_SUSPICIOUS_REGIONS } from './modules/suspicious_regions.nf'

/*
 * Pipeline parameters
 */

// Primary input
params {
    assembly: Path
    reads1: Path
    reads2: Path
    threads: Integer
    window_size: Integer
    breakpoint_bins: Integer
    min_suspicious_regions: Integer
    suspicious_flank_size: Integer
}

workflow {

    main:

    // Call processes

    assembly = params.assembly

    FILTER_CONTIGS(assembly)

    RUN_BOWTIE(assembly, FILTER_CONTIGS.out.filtered_fasta, params.reads1, params.reads2)

    CONTIG_LENGTHS(RUN_BOWTIE.out.sam)

    RUN_SAMTOOLS(assembly, RUN_BOWTIE.out.sam)

    CONTIG_COVERAGE(assembly, RUN_SAMTOOLS.out.pileup)

    SPLIT_PILEUP(RUN_SAMTOOLS.out.pileup, params.threads)

    chunks_ch = SPLIT_PILEUP.out.files.flatten()

    DEPTH_OF_COVERAGE(chunks_ch, params.window_size)

    merged_coverage_error_bed = DEPTH_OF_COVERAGE.out.coverage_errrors.collect()

    MERGE_BEDS_COVERAGE(merged_coverage_error_bed)

    BREAKPOINT_SPLITTER(RUN_BOWTIE.out.unaligned_reads)

    BREAKPOINT_FINDER(
        "${projectDir}/../src/py/nf/breakpoint_finder.py",
        FILTER_CONTIGS.out.filtered_fasta,
        BREAKPOINT_SPLITTER.out.breakpoint_reads,
        params.breakpoint_bins,
        params.threads,
        CONTIG_COVERAGE.out.coverage
    )

    BREAKPOINT_BED_SORT(BREAKPOINT_FINDER.out.breakpoint_bins)

    BIN_READS_AND_CONTIGS(
        RUN_BOWTIE.out.sam,
        CONTIG_COVERAGE.out.coverage,
        FILTER_CONTIGS.out.filtered_fasta,
        params.threads,
        10
    )

    bin_channel = BIN_READS_AND_CONTIGS.out.bin_dirs.flatten().filter{ path -> path.isDirectory() }.map { dir -> tuple(dir.name, dir) }

    RUN_BWA(bin_channel)

    SAM_TO_BAM(RUN_BWA.out)
    
    bin_channel_merged = bin_channel.join(SAM_TO_BAM.out)

    MATE_PAIR_CHECKER(bin_channel_merged, "${projectDir}/../src/py/mate_pairs.py")

    MERGE_BEDS_MATE_ERROR(MATE_PAIR_CHECKER.out.collect())

    GENERATE_SUMMARY(
        MERGE_BEDS_COVERAGE.out.merged_bed,
        BREAKPOINT_BED_SORT.out,
        MERGE_BEDS_MATE_ERROR.out,
        CONTIG_LENGTHS.out.contig_lengths,
        CONTIG_COVERAGE.out.coverage,
        FILTER_CONTIGS.out.contig_lengths
    )

    FIND_SUSPICIOUS_REGIONS(GENERATE_SUMMARY.out.summary_bed, params.suspicious_flank_size, params.min_suspicious_regions)

    

    publish:
    index = RUN_BOWTIE.out.index
    splits = BREAKPOINT_SPLITTER.out.breakpoint_reads
    summary = GENERATE_SUMMARY.out.summary_bed

}

output {

    index {
        path 'index'
    }
    splits {
        path 'splits'
    }
    summary {
        path 'summary'
    }
}
