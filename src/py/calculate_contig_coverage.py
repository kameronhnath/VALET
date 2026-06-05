import argparse

"""
Calculate contig coverage.  The coverage of a contig is the mean per-bp coverage.

Args:
    pileup_file: filename of the samtools formatted pileup file.
Returns:
    Filename of the coverage file.
"""

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-p", "--pileup_file",
        dest = "pileup_file"
    )

    parser.add_argument(
        "-o", "--out_filename",
        dest = "out_filename"
    )

    args = parser.parse_args()

    coverage_file = open(args.out_filename, 'w')

    prev_contig = None
    curr_contig = None

    length = 0
    curr_coverage = 0

    for record in open(args.pileup_file, 'r'):
        fields = record.strip().split()

        if prev_contig != fields[0]:
            if prev_contig:
                coverage_file.write(prev_contig + '\t' + str(float(curr_coverage) / length) + '\n')

            prev_contig = fields[0]
            length = 0
            curr_coverage = 0

        curr_coverage += int(fields[3])
        length += 1
    if prev_contig:
        coverage_file.write(prev_contig + '\t' + str(float(curr_coverage) / length) + '\n')
    coverage_file.close()

if __name__=="__main__":
    main()