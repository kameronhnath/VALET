import argparse

"""
Return a dictionary of contig names => contig lengths from a SAM file.

Args:
    sam_filename: input SAM filename
    out_filename: output file

Returns:
    Dictionary of contig names => contig lengths.
"""

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-s", "--sam_filename",
        dest = "sam_filename"
    )

    parser.add_argument(
        "-o", "--out_filename",
        dest = "out_filename"
    )

    args = parser.parse_args()

    outfile = open(args.out_filename, 'w')\
    
    sam_file = open(args.sam_filename, 'r')

    # Build dictionary of contig lengths.
    contig_lengths = {}
    line = sam_file.readline()
    while line.startswith("@"):

        if line.startswith("@SQ"):

            row = line.split()
            contig_name = row[1].split(':')[1]
            length = int(row[2].split(':')[1])
            contig_lengths[contig_name] = length
            outfile.write(str(contig_lengths[contig_name]) + '\t' + str(length) + '\n')

        line = sam_file.readline()

if __name__=="__main__":
    main()