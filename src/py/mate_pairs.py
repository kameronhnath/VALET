import pysam
import argparse
import numpy as np
from collections import defaultdict

#python src/py/mate_pairs.py -1 test/lib1.1.fastq -2 test/lib1.2.fastq -bam sequences/aln.bam

def main():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-1",
        dest = "first_mates",
        help = "First mates"
    )

    parser.add_argument(
        "-2",
        dest = "second_mates",
        help = "Second mates"
    )

    parser.add_argument(
        "-bam", "--bam",
        dest = "bam_file",
        help = "BAM file"
    )

    parser.add_argument(
        "-o", "--output",
        dest = "out_file",
        help = "Output file",
        default = "matepair_errors.bed"
    )

    parser.add_argument(
        "--binsize",
        dest = "bin_size",
        help = "How many base pairs the bins that contain the errors spans",
        default = 300
    )

    parser.add_argument(
        "-s","support",
        dest = "support",
        help = "Number of errors in a bin to constitute a suspicious region",
        default = 10
    )

    args = parser.parse_args()

    pairs = {}
    unruly_reads = []

    bam = pysam.AlignmentFile(args.bam_file, "rb")

    # bam.lengths calls a function each time - pulling it out of any loop to avoid runtime issues
    lengths = bam.lengths
    references = bam.references

    # For each read
    for read in bam.fetch(until_eof=True):

        if read.is_paired and not read.is_unmapped:

            # Flag any mate pairs on different contigs, or mate pairs with the same orientation
            if read.is_reverse == read.mate_is_reverse or read.reference_name != read.next_reference_name:
                unruly_reads.append([read.reference_name, read.reference_start])

            # Get read name
            name = read.query_name

            # Maintain a pair record of the following data
            # 0 - Contig id
            # 1 - Start position of forward pair
            # 2 - End position of forward pair
            # 3 - Start position of reverse pair
            # 4 - End position of reverse pair
            pair_record = [-1, -1, -1, -1, -1]

            if name in pairs:
                pair_record = pairs[name]

            # Fill in the pair record
            if read.is_reverse:
                pair_record[4] = read.reference_end
                pair_record[3] = read.reference_start
            else:
                pair_record[1] = read.reference_start
                pair_record[2] = read.reference_end
                pair_record[0] = read.reference_id

            # Save the pair record
            pairs[name] = pair_record

    # QC - Remove any pairs that dont have both mates
    pairs_to_delete = []
    for name, pair in pairs.items():
        if not (pair[1] != -1 and pair[3] != -1):
            pairs_to_delete.append(name)  
    for name in pairs_to_delete:
            del pairs[name]

    # Remove reads that are too close to the ends of the contigs - recalculate statistics and repeat until stable
    removed_reads = True
    while removed_reads:

        removed_reads = False

        distances = []
        for pair in pairs.values():
            distances.append(pair[4] - pair[1])

        # Get statistics
        dist_data = np.array(distances)
        mean = np.mean(dist_data)
        sd = np.std(dist_data)
        lower_cutoff = mean

        pairs_to_delete = []
        for name, pair in pairs.items():
            contig_length = lengths[pair[0]]
            upper_cutoff = contig_length - mean
            if not (pair[4] > lower_cutoff and pair[1] < upper_cutoff):
                pairs_to_delete.append(name)
                removed_reads = True
        for name in pairs_to_delete:
            del pairs[name]
    
    # Flag any outliers of library size
    upper_bound = mean + sd
    lower_bound = mean + sd
    median = np.median(dist_data)
    q1 = np.percentile(dist_data, 25)
    q3 = np.percentile(dist_data, 75)
    for name, pair in pairs.items():
        distance = pair[4] - pair[1]
        if distance < lower_bound or distance > upper_bound:
            unruly_reads.append([references[pair[0]], pair[1]])
            unruly_reads.append([references[pair[0]], pair[3]])
    

    # Bin outliers
    bins = defaultdict(int)
    for read in unruly_reads:
        contig = read[0]
        bin_start = (read[1] // args.bin_size) * args.bin_size
        bins[(contig,bin_start)] += 1

    # Write results to bed file
    with open(args.outfile, "w") as out_file:
        for (contig, start), support in bins.items():
            if support > args.support:
                out_file.write(
                    f"{contig}\t{start}\t{start + args.bin_size}\tMatePairError\n"
                )

if __name__=="__main__":
    main()