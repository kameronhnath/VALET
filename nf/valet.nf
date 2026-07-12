#!/usr/bin/env nextflow

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

params.assembly = null
params.reads1 = null
params.reads2 = null

params.threads = 4
params.window_size = 501
params.breakpoint_bins = 50
params.min_suspicious_regions = 2
params.suspicious_flank_size = 0 //1000
params.min_coverage = 10 // 10
params.ignore_end_distances = 0 //150

workflow {

    main:

    // Print params
    params.each { key, value ->
        println "${key}: ${value}"
    }

    // Grab assembly file
    assembly = file(params.assembly)

    FILTER_CONTIGS(assembly, "${projectDir}/../src/py/filter_short_contigs.py")

    RUN_BOWTIE(assembly, FILTER_CONTIGS.out.filtered_fasta, file(params.reads1), file(params.reads2))

    CONTIG_LENGTHS(RUN_BOWTIE.out.sam, "${projectDir}/../src/py/get_contig_lengths.py")

    RUN_SAMTOOLS(assembly, RUN_BOWTIE.out.sam)

    CONTIG_COVERAGE(assembly, RUN_SAMTOOLS.out.pileup, "${projectDir}/../src/py/calculate_contig_coverage.py")

    SPLIT_PILEUP(RUN_SAMTOOLS.out.pileup, params.threads, "${projectDir}/../src/py/split_pileup.py")

    chunks_ch = SPLIT_PILEUP.out.files.flatten()

    DEPTH_OF_COVERAGE(chunks_ch, params.window_size, "${projectDir}/../src/py/nf/depth_of_coverage.py")

    merged_coverage_error_bed = DEPTH_OF_COVERAGE.out.coverage_errrors.collect()

    MERGE_BEDS_COVERAGE(merged_coverage_error_bed)

    BREAKPOINT_SPLITTER(RUN_BOWTIE.out.unaligned_reads, "${projectDir}/../src/py/nf/breakpoint_splitter.py")

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
        RUN_BOWTIE.out.paired_sam,
        CONTIG_COVERAGE.out.coverage,
        FILTER_CONTIGS.out.filtered_fasta,
        params.threads,
        params.min_coverage,
        "${projectDir}/../src/py/nf/bin_reads.py",
        "${projectDir}/../src/py/nf/bin_contigs.py"
    )

    // Convert bins into a channel
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
        FILTER_CONTIGS.out.contig_lengths,
        params.ignore_end_distances,
        "${projectDir}/../src/py/nf/summary_table.py"
    )

    FIND_SUSPICIOUS_REGIONS(GENERATE_SUMMARY.out.summary_bed, params.suspicious_flank_size, params.min_suspicious_regions)

    
    publish:
    suspicious = FIND_SUSPICIOUS_REGIONS.out
    summary = GENERATE_SUMMARY.out.summary_bed
    coverage = MERGE_BEDS_COVERAGE.out.merged_bed
    breakpoint = BREAKPOINT_BED_SORT.out
    mate_error = MERGE_BEDS_MATE_ERROR.out
    summary_table = GENERATE_SUMMARY.out.summary_table
    alignment = RUN_SAMTOOLS.out.sorted_bam

}

output {

    suspicious {
        path 'suspicious'
    }
    summary {
        path 'summary'
    }
    coverage {
        path 'results'
    }
    breakpoint {
        path 'results'
    }
    mate_error {
        path 'results'
    }
    summary_table {
        path 'summary_table'
    }
    alignment {
        path 'results'
    }
    
}
