import argparse

""" Generate the summary files from the individual error files -> Then generate a summary table

Args:
    coverage_bed: Output from the depth of coverage checker step.
    breakpoint_bed: Output from the breakpoint checker.
    matepair_bed: Output from the mate pair error checker.
    contig_lengths_file: Tab separated file of contig lengths.
    summary_file: Summary file.
    abundance_filename: Contig abundances file.
    filtered_contigs: Tab separated file of filtered contig lengths.
    table_filename: Table file output.

Returns:
    Final mis-assemblies in tuple format.
"""

# Helper function to load in the contig lengths file
def load_contig_lengths(filename):
    contig_lengths = {}

    with open(filename) as infile:
        for line in infile:
            contig, length = line.rstrip().split('\t')
            contig_lengths[contig] = int(length)

    return contig_lengths

# Helper function for contig abundances
def get_contig_abundances(abundance_filename):
    """
    Return a dictionary of contig names => contig abundances from the '/coverage/temp.cvg'.

    Returns:
        A dictionary mapping from contig names to abundances.
    """

    abundance_file = open(abundance_filename, 'r')

    # Build a dictionary of contig abundances.
    contig_abundances = {}

    for line in abundance_file:
        contig_and_abundance = line.strip().split()
        contig_abundances[contig_and_abundance[0]] = int(round(float(contig_and_abundance[1])))

    return contig_abundances

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-c", "--coverage_bed",
        dest = "coverage_bed"
    )

    parser.add_argument(
        "-b", "--breakpoint_bed",
        dest = "breakpoint_bed"
    )

    parser.add_argument(
        "-m", "--matepair_bed",
        dest = "matepair_bed"
    )

    parser.add_argument(
        "-l", "--contig_lengths_file",
        dest = "contig_lengths_file"
    )

    parser.add_argument(
        "-s", "--summary_file",
        dest = "summary_file"
    )

    parser.add_argument(
        "-i", "--ignore_end_distances",
        dest = "ignore_end_distances",
        default=100, type=int
    )

    parser.add_argument(
        "-a", "--abundance_filename",
        dest = "abundance_filename"
    )

    parser.add_argument(
        "-f", "--filtered_contigs",
        dest = "filtered_contigs"
    )

    parser.add_argument(
        "-t", "--table_filename",
        dest = "table_filename"
    )

    # Parse args and set up related variables using helper functions
    args = parser.parse_args()
    contig_lengths = load_contig_lengths(args.contig_lengths_file)
    all_contig_lengths = contig_lengths
    results_filenames = [args.coverage_bed, args.breakpoint_bed, args.matepair_bed]
    contig_abundances = get_contig_abundances(args.abundance_filename)
    filtered_contig_lengths = load_contig_lengths(args.filtered_contigs)

    table_file = open(args.table_filename, 'w')

    summary_file = open(args.summary_file, 'w')

    misassemblies = []
    for results_file in results_filenames:
        if results_file:
            for line in open(results_file, 'r'):
                misassemblies.append(line.strip().split('\t'))

    # Sort misassemblies by start site.
    misassemblies.sort(key = lambda misassembly: (misassembly[0], int(misassembly[1]), int(misassembly[2])))
    final_misassemblies = []
    for misassembly in misassemblies:

        # Truncate starting/ending region if it is near the end of the contigs.
        if int(misassembly[1]) <= args.ignore_end_distances and \
            int(misassembly[2]) > args.ignore_end_distances:
          misassembly[1] = str(args.ignore_end_distances + 1)

        if int(misassembly[2]) >= (contig_lengths[misassembly[0]] - args.ignore_end_distances) and \
            int(misassembly[1]) < (contig_lengths[misassembly[0]] - args.ignore_end_distances):
          misassembly[2] = str(contig_lengths[misassembly[0]] - args.ignore_end_distances - 1)

        # Don't print a flagged region if it occurs near the ends of the contig.
        if int(misassembly[1]) > args.ignore_end_distances and \
                int(misassembly[2]) < (contig_lengths[misassembly[0]] - args.ignore_end_distances):
            summary_file.write('\t'.join(misassembly) + '\n')

            final_misassemblies.append(misassembly)

    summary_file.close()

    table_file.write("contig_name\tcontig_length\tabundance\tlow_cov\tlow_cov_bps\thigh_cov\thigh_cov_bps\tmate_error\tmate_error_bps\tbreakpoints\tbreakpoints_bps\n")

    prev_contig = None
    curr_contig = None

    # Misassembly signatures
    low_coverage = 0
    low_coverage_bps = 0
    high_coverage = 0
    high_coverage_bps = 0
    mate_error = 0
    mate_error_bps = 0
    breakpoints = 0
    breakpoints_bps = 0

    processed_contigs = set()

    for misassembly in final_misassemblies:

        curr_contig = misassembly[0]

        if prev_contig is None:
            prev_contig = curr_contig

        if curr_contig != prev_contig:
            # Output previous contig stats.
            table_file.write(prev_contig + '\t' + str(filtered_contig_lengths[prev_contig]) + '\t' + str(contig_abundances[prev_contig]) + '\t' + \
                str(low_coverage) + '\t' + str(low_coverage_bps) + '\t' + str(high_coverage) + '\t' + \
                str(high_coverage_bps) + '\t' + str(mate_error) + '\t' + str(mate_error_bps) + '\t' + str(breakpoints) + '\t' + \
                str(breakpoints_bps) + '\n')

            processed_contigs.add(prev_contig)

            # Reset misassembly signature counts.
            low_coverage = 0
            low_coverage_bps = 0
            high_coverage = 0
            high_coverage_bps = 0
            mate_error = 0
            mate_error_bps = 0
            breakpoints = 0
            breakpoints_bps = 0

            prev_contig = curr_contig

        # Process the current contig misassembly. Mate_error_bps is largely meaningless since mate pair errors are established with a specific bin width.
        if misassembly[3] == 'MatePairError':
            mate_error += 1
            mate_error_bps += (int(misassembly[2]) - int(misassembly[1]) + 1)

        elif misassembly[3] == 'Low_coverage':
            low_coverage += 1
            low_coverage_bps += (int(misassembly[2]) - int(misassembly[1]) + 1)

        elif misassembly[3] == 'High_coverage':
            high_coverage += 1
            high_coverage_bps += (int(misassembly[2]) - int(misassembly[1]) + 1)

        elif misassembly[3] == 'Breakpoint_finder':
            breakpoints += 1
            breakpoints_bps += (int(misassembly[2]) - int(misassembly[1]) + 1)

        else:
            print("Unhandled error: " + misassembly[3])

    if prev_contig and prev_contig in contig_abundances.keys():
        # Output previous contig stats.
        table_file.write(prev_contig + '\t' + str(filtered_contig_lengths[prev_contig]) + '\t' + str(contig_abundances[prev_contig]) + '\t' + \
            str(low_coverage) + '\t' + str(low_coverage_bps) + '\t' + str(high_coverage) + '\t' + \
            str(high_coverage_bps) + '\t' + str(mate_error) + '\t' + str(mate_error_bps) + '\t' + str(breakpoints) + '\t' + \
            str(breakpoints_bps) + '\n')

        processed_contigs.add(prev_contig)

    # We need to add the remaining, error-free contigs.
    for contig in filtered_contig_lengths:
        if contig not in processed_contigs and contig not in filtered_contig_lengths and contig in contig_abundances.keys():
            table_file.write(contig + '\t' + str(filtered_contig_lengths[contig]) + '\t' + str(contig_abundances[contig]) + '\t' + \
                '0\t0\t0\t0\t0\t0\t0\t0\n')
            processed_contigs.add(contig)


    # Finally, add the contigs that were filtered out prior to evaluation.
    for contig in all_contig_lengths:
        if contig not in processed_contigs:
            table_file.write(contig + '\t' + str(all_contig_lengths[contig]) + '\t' + 'NA\t' + \
                'NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n')
            processed_contigs.add(contig)

    table_file.close()

if __name__=="__main__":
    main()