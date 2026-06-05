import os
import argparse
import math
import subprocess

"""
Bins reads by coverage.

Args:
    samfile: input SAM filename
    coverage_file
    threads
    min_coverage

Returns:
    Directory of bins.
"""

# http://stackoverflow.com/questions/19570800/reverse-complement-dna
revcompl = lambda x: ''.join([{'A':'T','C':'G','G':'C','T':'A','N':'N','R':'N','M':'N','Y':'N','S':'N','W':'N','K':'N'}[B] for B in x][::-1])

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

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-s", "--samfile",
        dest = "samfile"
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

    # First create a file in the format (abundance, header, seq, quality)
    abundance_read_filename = output_dir + '/bins/abun_read'
    ensure_dir(abundance_read_filename)
    abundance_read_file = open(abundance_read_filename, 'w')
    
    # log file for read coverage binning
    abundance_log_filename = output_dir + '/bins/abun_read.log'
    ensure_dir(abundance_log_filename)
    abundance_log_file = open(abundance_log_filename, 'w')    

    # Skip the header sequence.
    sam_file = open(args.samfile, 'r')
    line = sam_file.readline()
    while line.startswith("@"):
        line = sam_file.readline()

    seq = None
    quals = None
    while line:
        tuple = line.split('\t')

        ## check for missing contigs not present on abundance file
        ## This is potentially due to only a few reads mapping to the contig        
    ##if tuple[2] not in contig_abundances.keys():
      #  abundance_log_file.write(tuple[2] + "not in contig abundance file read " +  tuple[0] +  " excluded from coverage and read pair analysis")
            #line = sam_file.readline()
        #continue

        if tuple[2] != '*':
            if tuple[2] not in contig_abundances.keys():
                abundance_log_file.write(tuple[2] + "not in contig abundance file read " +  tuple[0] +  " excluded from coverage and read pair analysis\n")
                line = sam_file.readline()
                continue
            elif int(tuple[1]) & 0x10 == 0:
                seq = tuple[9]
                quals = tuple[10]
            else:
                seq = revcompl(tuple[9])
                quals = tuple[10][::-1]

            abundance_read_file.write(str(int(math.ceil(contig_abundances[tuple[2]]))) + '\t' + tuple[0] + '\t' + seq + '\t' + quals + '\n')

        line = sam_file.readline()
    abundance_read_file.close()

    # Sort the abundance file by (1) abundance, (2) first mate, (3) second mate.
    try:
        call_arr = ['sort', '-nk1,1', '-k2,2', '-T', './', '--parallel=' + str(args.threads), abundance_read_filename, '-o', abundance_read_filename + '.sorted']
        subprocess.check_call(call_arr)
    except:
        call_arr = ['sort', '-nk1,1', '-k2,2', '-T', './', abundance_read_filename, '-o', abundance_read_filename + '.sorted']
        subprocess.call(call_arr)

    # Write out each read to there correct bin folder.
    path_to_bins = []

    abundance_read_file = open(abundance_read_filename + '.sorted', 'r')
    """
    5       HWUSI-EAS626_102891784:2:90:12172:7648/1        GGTCGTGTTTCATTGGTTAAATCACCAAATCCTTCCATATCCACGATACCACGCATATCTTTTTTCTTTAAAAATTCACGCAACCCATCTTTAACTTCAG    ??90AA:6@?4A??<@???@4.1+64=<63@:??@?@??5EE>EEEA@FFCEE=AE=B@CDDDEBBFFBFFFFFEDFEFFGGFGGGGGGGGGGGGGGGGG
    5       HWUSI-EAS626_102891784:2:90:12172:7648/2        CTGTAAAATCATCATCTACATTAAAGTAAACAGGATTACCGTCGGCCTTAATACCATATATACCAATCTGATTACCAACACTTGAAATATCGTCTCACTC    GGFGGGGGGGGGGGGGGGGGGGGGEGCGGGGGGGFEFFFD@?@@@EEEAEEEEEE?DDDDEEEEECCB@@@CAACEEEEECA;A9AAA??%%%%%%%%%%
    5       HWUSI-EAS626_102891784:2:94:12848:7954/1        TTTAAAAATTCACGCAACCCATCTTTAACTTCAGGTTTTACCTCATCATATGTAAAAATCTCATTATTTTCAACAGTAATATTTGAATAATCTTTAGGAC    A?AACE@EEBDAEDDCADAECDECDBDCF@DEEDEEBEFCEGEGGFGGGGFFGGGFFEFFAC@3@GFGGGGGGGGGGGGGGGGFGGGGGGFGEGGGGGGG
    5       HWUSI-EAS626_102891784:2:94:12848:7954/2        AATGTCTTTTTTAATCGGTTGGCCATTTTCATCTGTTTCAGAGCCATAAACTGGTCGTGTTTCATTGGTTAAATCACCAAATCCTTCCATATCCACGATA    FGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGFGGCGGEGGFGDEFADGDDGFDEEEEFEEDGEEEEEEEBDDGDEEEEEDCEBD@?@?@
    7       HWUSI-EAS626_102891784:1:100:14551:5313/1       TTCTTGCCTGTCTTTAATCCGAANAGCGCCCCTGACCAATCAGAAAAAGCACGCCTGTTGCTTGCCGCAAACTGACTATATAAGTCTTGTG     GGGGGFGGGGGGGGGGGGGGCCC!CCDDDDGGG:GGGDGDGFFFFGGGEGGGEGEFACDFDDEEFEEGEEEECEEEGAEEEEBEBDACBB#
    7       HWUSI-EAS626_102891784:1:102:10765:10655/1      CAGCAATCCAGTCTTTAACTTCTGGGTGCCATGCAGGATGCGGTATATAAACCTGTCCAGCTTCCCACATTGGAGACACTGACGCCGCACG     GGGGGGGGGGGGGGGGGGGGGGGGAECADEFFFFEGGGGGFGGEEAFFFDFGGGFEGEGEGBEGEGEEEDEE?EBEDEEEBCDBCECB=A#
    7       HWUSI-EAS626_102891784:1:102:10765:10655/2      GCTTTGCAGTAGCGTCAGGATACATGCGGGACATGGCTCTAATAGCGTCTAGCGTCTCAGTAAAGCTTAAACGCTTGTGGCACCAGTTAGGGCGCAGGTA    >>>?=,?:D?D?CCCBCDEF5BEEEFGDGEFCD?EGGEBGFEGGGGFFGGEGGGFGGGDGFGFGGEEEE?AFFBGFFFFEECAFFGGG=FFFEGGDFGGG
    """

    
    pairs = {}
    bin_files = {}

    # Sometimes the sequence files are missing mates, we need to remove any unpaired sequence.
    for line in abundance_read_file:

        # Get all elements of the line
        line_parts = line.strip().split("\t")
        abundance = line_parts[0]

        # Remove any spaces or characters that come before the read number
        readname = line_parts[1]
        readname = readname.split()[0]
        readname = readname.split('/')[0]
        
        sequence = line_parts[2]
        quality_score = line_parts[3]

        # Skip the read if it's below the min coverage
        if int(abundance) < args.min_coverage: 
            continue

        # Initialize the bin files and directory if needed
        if abundance not in bin_files:
            bin_directory = os.path.join(output_dir, 'bins', abundance)
            os.makedirs(bin_directory, exist_ok=True)
            path_to_bins.append(bin_directory + '/')

            f1 = open(os.path.join(bin_directory, 'lib_1.fq'), 'w')
            f2 = open(os.path.join(bin_directory, 'lib_2.fq'), 'w')
            bin_files[abundance] = (f1, f2)
        
        # Add the read to the pairs dict if it is the first one found, otherwise add the pair to the bin files
        if readname not in pairs:
            pairs[readname] = (abundance, sequence, quality_score)
        else:
            prev_abundance, prev_sequence, prev_quality = pairs.pop(readname)

            f1, f2 = bin_files[prev_abundance]

            f1.write(f"@{readname}/1\n{prev_sequence}\n+\n{prev_quality}\n")
            f2.write(f"@{readname}/2\n{sequence}\n+\n{quality_score}\n")

    # Close file writers
    for f1, f2 in bin_files.values():
        f1.close()
        f2.close()

if __name__=="__main__":
    main()