import os
import argparse
import math
import subprocess

""" Bin assembly by their coverages.

Args:
    options: Command line options.
    assembly_filename: Assembly FASTA filename.
    contig_abundances: Dictionary containing contig_name => abundance.
    output_dir: Current assembly output directory.
"""

# Helper functions

def ensure_dir(f):
    d = os.path.dirname(f)
    if not os.path.exists(d):
        os.makedirs(d)
    assert os.path.exists(d)

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

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-a", "--assembly",
        dest = "assembly"
    )

    parser.add_argument(
        "-c", "--coverage_file",
        dest = "coverage_file"
    )

    parser.add_argument(
        "-t", "--threads",
        dest = "threads",
        type=int
    )

    parser.add_argument(
        "-m", "--min_coverage",
        dest = "min_coverage",
        type=int
    )

    args = parser.parse_args()

    output_dir = "out"

    contig_abundances = get_contig_abundances(args.coverage_file)

    # First create a file in the format (abundance, header, seq)
    abundance_contig_filename = output_dir + '/bins/abun_contig'
    ensure_dir(abundance_contig_filename)
    abundance_contig_file = open(abundance_contig_filename, 'w')


    # Contig abundances log file
    abundance_log_filename = output_dir + '/bins/abun_contig.log'
    ensure_dir(abundance_log_filename)
    abundance_log_file = open(abundance_log_filename, 'w')

    with open(args.assembly, 'r') as assembly:
        for contig in contig_reader(assembly):
            '''
            COMMENTED NATE's ISSUE
            '''
            # if contig['name'] not in contig_abundances.keys():
            #     abundance_log_file.write(contig['name'] + "not in contig abundance file, excluded read pair analysis")
            #     continue
            abundance_contig_file.write(str(int(math.ceil(contig_abundances[contig['name'][1:].strip()]))) + '\t' +\
                    contig['name'][1:].strip() + '\t' + ''.join(contig['sequence']).strip() + '\n')

    abundance_contig_file.close()
    abundance_log_file.close()

    # Sort the contigs by abundance.
    try:
        call_arr = ['sort', '-nk1,1', '-k2,2', '-T', './', '--parallel=' + str(args.threads), abundance_contig_filename, '-o', abundance_contig_filename + '.sorted']
        subprocess.check_call(call_arr)
    except:
        call_arr = ['sort', '-nk1,1', '-k2,2', '-T', './', abundance_contig_filename, '-o', abundance_contig_filename + '.sorted']
        subprocess.call(call_arr)


    prev_abun = None
    curr_abun = None

    abundance_contig_file = open(abundance_contig_filename + '.sorted', 'r')
    for line in abundance_contig_file:
        tuple = line.split('\t')
        curr_abun = tuple[0]

        if int(curr_abun) < args.min_coverage:
            continue

        if prev_abun is None or curr_abun != prev_abun:
            # Setup the writer.
            contig_writer = open(output_dir + '/bins/' + curr_abun + '/contigs.fasta', 'w')

        contig_writer.write('>' + tuple[1] + '\n' + tuple[2])

        prev_abun = curr_abun

    abundance_contig_file.close()
    abundance_log_file.close()

if __name__=="__main__":
    main()