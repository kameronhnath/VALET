import random
from pathlib import Path
import math
import argparse

# Python script to generate test data for valet. Generates a genome and draws reads from it given some parameters. Then adds errors to the genome and tracks those errors.
##  ARGS:
##      -o (--outfile)          -> output directory
##      -s (--seed)             -> seed for RNG to ensure consistent outputs
##      -e (--volume_errors)    -> number of errors added per contig
##      -c (--num_contigs)      -> number of contigs
##      -g (--size)             -> length of each contig

# Generate a random string of nucleotides
def random_genome(length):
    return ''.join(random.choices('ACGT', k=length))

def reverse_complement(seq):
    comp = str.maketrans("ACGT", "TGCA")
    return seq.translate(comp)[::-1]

# Create non-overlapping regions within a genome
def create_regions(length, num, min_size, max_size):

    # Impossible when regions are too large
    if num*min_size >= length:
        return []
    
    regions = []
    while len(regions) < num:
        size = random.randint(min_size,max_size)
        gen_start = random.randint(0, length - size)
        gen_end = gen_start + size
        overlap = False

        # Ensure regions do not overlap
        for (start,end) in regions:
            if gen_start < end and gen_end > start:
                overlap = True
                break

        if not overlap:
            regions.append((gen_start,gen_end))

    return sorted(regions, reverse=True)

# Given a 'genome' randomly generate paired reads. 
##  Reads are drawn from random intervals within the genome.
##  Coverage deterimines how many reads are generated.
##  Mean and sd insert are used to determine the random size of the reads.
##  Low-coverage regions can be added that have a lower chance of reads being drawn.
def reads_from_genome(
        genome,
        r1_file,
        r2_file,
        contig_name,
        coverage = 10,
        mean_insert = 300,
        sd_insert = 50,
        low_coverage_regions = []):
    
    # Calculate the number of reads needed to approximately reach the coverage threshold
    genome_length = len(genome)
    num_fragments = int((coverage * genome_length) / mean_insert)

    for i in range(num_fragments):
        # Generate a random insert size using a gaussian distribution about the mean
        insert_size = int(random.gauss(mean_insert, sd_insert))
        read_length = math.floor(insert_size / 2 - 20)

        # Determine the start and end of the read
        start = random.randint(0, genome_length - insert_size)
        end = start + insert_size

        fragment = genome[start:end]

        forward_read = fragment[:read_length]
        reverse_read = reverse_complement(fragment[-read_length:])

        # Dummy perfect quality scores
        quality_scores = "I" * read_length

        read_id = f"@SIM:{contig_name}:{i}:start={start}:insert={insert_size}"

        # If the random read is in low coverage region there is a 90% chance to not include it
        include = True
        for (reg_start, reg_end) in low_coverage_regions:
            if (start > reg_start and start < reg_end) or (end > reg_start and end < reg_end):
                if random.randint(0,100) < 90:
                    include = False
                break

        if include:
            # Write to FASTQ files
            r1_file.write(f"{read_id}/1\n{forward_read}\n+\n{quality_scores}\n")
            r2_file.write(f"{read_id}/2\n{reverse_read}\n+\n{quality_scores}\n")
    
    return low_coverage_regions
    
# Write genome to file with fasta header
def write_genome_to_file(genome,filepath):
    outfile = open(filepath + "contigs.fasta", "w")
    for name,sequence in genome:
        outfile.write(">" + name + "\n")
        outfile.write(sequence + "\n")

# Write error info to a bed file
def errors_to_bed(errors, outfile):
    with open(outfile, "w") as f:
        for err in errors:
            f.write(
                f"{err['contig']}\t"
                f"{err['start']}\t"
                f"{err['end']}\t"
                f"{err['type']}\n"
            )

def main():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-o", "--outfile",
        dest="file_path",
        default="./"
    )

    parser.add_argument(
        "-s", "--seed",
        dest="seed",
        default=0
    )

    parser.add_argument(
        "-e", "--volume_errors",
        dest="volume_errors",
        default=10
    )

    parser.add_argument(
        "-c", "--num_contigs",
        dest="num_contigs",
        default=1, type=int
    )

    parser.add_argument(
        "-g", "--size",
        dest="size",
        default=1000000, type=int
    )

    args = parser.parse_args()
    e = int(args.volume_errors)

    # Open read files
    r1_file = open(args.file_path + "R1.fastq", "w")
    r2_file = open(args.file_path + "R2.fastq", "w")
    
    random.seed(args.seed)

    contigs = []
    errors = []
    for i in range(1,args.num_contigs+1):

        # Generate genome (for the contig)
        genome = random_genome(args.size)

        # Determine a random coverage value
        coverage = random.randint(10,90)

        # Generate error regions
        error_regions = create_regions(len(genome),e,1000,1500)

        # Draw low-coverage regions from the error regions
        num_low_coverage_regions = math.floor(e / 5)
        low_coverage_regions = error_regions[:num_low_coverage_regions]

        # Generate the reads using the random coverage value, label the contig with the current index
        low_coverage_regions = reads_from_genome(
            genome=genome,
            r1_file=r1_file,
            r2_file=r2_file,
            contig_name=i,
            coverage=coverage,
            mean_insert=300,
            sd_insert=50,
            low_coverage_regions=low_coverage_regions)
        
        # Add errors records for the low coverage regions
        contig_errors = []
        for start, end in low_coverage_regions:
            contig_errors.append({
                "type":"low_coverage",
                "contig":"contig" + str(i),
                "start":start,
                "end":end
            })

        # Add the remaining erros
        for start, end in error_regions:

            # Pick an error type
            error_type = random.choice([
                "deletion",
                "inversion",
                "duplication"
            ])

            match error_type:
                case "deletion":

                    genome = genome[:start] + genome[end:]

                    # Adjust position of existing errors
                    for err in contig_errors:
                        if err['start'] > start:
                            err['start'] = err['start'] - (end - start)
                            err['end'] = err['end'] - (end - start)

                    # Add the new error record
                    contig_errors.append({
                        "type":"deletion",
                        "contig":"contig" + str(i),
                        "start":start,
                        "end":end
                    })
                
                case "inversion":

                    genome = genome[:start] + reverse_complement(genome[start:end]) + genome[end:]

                    contig_errors.append({
                        "type":"inversion",
                        "contig":"contig" + str(i),
                        "start":start,
                        "end":end
                    })

                case "duplication":

                    genome = genome[:end] + genome[start:end] + genome[end:]

                    err_length = end - start

                    for err in contig_errors:
                        if err['start'] > start:
                            err['start'] = err['start'] + err_length
                            err['end'] = err['end'] + err_length

                    contig_errors.append({
                        "type":"duplication",
                        "contig":"contig" + str(i),
                        "start":end,
                        "end":end + err_length
                    })

                case "nearby_duplication":

                    location = random.randint(start,min(len(genome),end + 500))
                    genome = genome[:location] + genome[start:end] + genome[location:]

                    err_length = end - start

                    contig_errors.append({
                        "type":"nearby_duplication",
                        "contig":"contig" + str(i),
                        "start":location,
                        "end":location + err_length
                    })

                case "relocation":

                    distance = random.randint(5000,10000)

                    err_length = end - start

                    direction = random.randint(0,1)
                    if direction == 1:
                        relocation = min(len(genome), end + distance)
                        genome = genome[:start] + genome[end:relocation] + genome[start:end] + genome[relocation:]
                    else:
                        relocation = max(0, start - distance)
                        genome = genome[:relocation] + genome[start:end] + genome[relocation:start] + genome[end:]

                    contig_errors.append({
                        "type":"relocation_removed",
                        "contig":"contig" + str(i),
                        "start":start,
                        "end":start
                    })

                    contig_errors.append({
                        "type":"relocation_placed",
                        "contig":"contig" + str(i),
                        "start":relocation,
                        "end":relocation + err_length
                    })
            
        contigs.append(["contig" + str(i), genome])
        errors = errors + contig_errors
    
    write_genome_to_file(contigs,args.file_path)

    errors_to_bed(errors, args.file_path + "errors.bed")
    

if __name__=="__main__":
    main()