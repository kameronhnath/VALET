import argparse

def main():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-f", "--fasta_filename",
        dest = "fasta_filename"
    )

    parser.add_argument(
        "-ff", "--filtered_fasta_filename",
        dest = "filtered_fasta_filename"
    )

    parser.add_argument(
        "-lf", "--length_filename",
        dest = "length_filename"
    )

    parser.add_argument(
        "-m", "--min_contig_length",
        dest = "min_contig_length",
        type=int,
        default = 1000
    )

    args = parser.parse_args()

    filtered_assembly_file = open(args.filtered_fasta_filename, 'w')
    length_file = open(args.length_filename, 'w')
    curr_length = 0
    with open(args.fasta_filename,'r') as assembly:
        for contig in contig_reader(assembly):
            curr_length = len(''.join(contig['sequence']))

            if curr_length >= args.min_contig_length:
                filtered_assembly_file.write(contig['name'])
                filtered_assembly_file.writelines(contig['sequence'])
                filtered_assembly_file.write('\n')

            length_file.write(contig['name'].strip()[1:] + "\t" + str(curr_length) + "\n")

    filtered_assembly_file.close()
    length_file.close()

def contig_reader(fasta_file):
    """ Iterator that takes a fasta file and returns entries.
    IMPORTANT NOTE: names are stripped at first space.
    """
    save_line = ""
    contig = {}
    in_contig = False
    for line in fasta_file:
        if line[0] == '>' and in_contig:
            save_line = line
            ret_contig = contig
            contig = {}
            contig['sequence'] = []
            contig['name'] = line.split()[0].strip() + "\n"
            yield ret_contig
        elif line[0] == '>':
            contig['name'] = line.split()[0].strip() + "\n"
            contig['sequence'] = []
            in_contig = True
        else:
            contig['sequence'].append(line.strip())
    yield contig

if __name__=="__main__":
    main()