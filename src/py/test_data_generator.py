import random
from pathlib import Path
import math
import argparse

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
            if (gen_start > start and gen_start < end) or (gen_end > start and gen_end < end):
                overlap = True
                break

        if not overlap:
            regions.append((gen_start,gen_end))

    return regions

# Given a 'genome' randomly generate paired reads. 
##  Reads are drawn from random intervals within the genome.
##  Coverage deterimines how many reads are generated.
##  The seed parameter is used to seed the rng.
##  Mean and sd insert are used to determine the random size of the reads.
##  Reads are written to fasta files in the outfile path.
##  Low-coverage regions can be added that have a lower chance of reads being drawn.
def reads_from_genome(
        genome,
        coverage = 10,
        mean_insert = 300,
        sd_insert = 50,
        seed = 42,
        outfile = "./",
        num_low_coverage_regions = 0):
    
    random.seed(seed)

    # Calculate the number of reads needed to approximately reach the coverage threshold
    genome_length = len(genome)
    num_fragments = int((coverage * genome_length) / mean_insert)

    # Paired read files
    r1_file = open(outfile + "R1.fastq", "w")
    r2_file = open(outfile + "R2.fastq", "w")

    # Create low coverage regions
    low_coverage_regions = create_regions(genome_length, num_low_coverage_regions, 3000,3000)

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

        read_id = f"@SIM:{i}:start={start}:insert={insert_size}"

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


### FUNCTIONS TO ADD ERRORS TO THE INITIAL GENOME ###

def add_deletions(genome,num,min_size,max_size):
    # Get random non-overlapping regions
    regions = create_regions(len(genome),num,min_size,max_size)

    # Remove them in reverse sorted order
    for (start,end) in sorted(regions, reverse=True):
        genome = genome[:start] + genome[end:]

    # Return the subsetted genome
    return genome

def add_tandem_insertions(genome,num,min_size,max_size):
    # Get random non-overlapping regions
    regions = create_regions(len(genome),num,min_size,max_size)

    # Duplicate them in reverse sorted order
    for (start,end) in sorted(regions, reverse=True):
        genome = genome[:end] + genome[start:end] + genome[end:]

    # Return the subsetted genome
    return genome

def add_nearby_insertions(genome,num,min_size,max_size,range):
    regions = create_regions(len(genome),num,min_size,max_size)

    for (start,end) in sorted(regions, reverse=True):
        location = random.randint(start,min(len(genome),end + range + 200))
        genome = genome[:location] + genome[start:end] + genome[location:]

    return genome

def add_inversions(genome,num,min_size,max_size):
    regions = create_regions(len(genome),num,min_size,max_size)

    for (start,end) in sorted(regions):
        genome = genome[:start] + reverse_complement(genome[start:end]) + genome[end:]
    
    return genome

def add_duplications(genome,num,min_size,max_size):
    regions = create_regions(len(genome),num,min_size,max_size)

    for (start,end) in sorted(regions):
        genome = genome[:end] + genome[start:end] + genome[end:]
    
    return genome

def add_relocation(genome,size,distance):
    start = random.randint(distance+1, len(genome) - size - distance)
    end = start + size
    direction = random.randint(0,1)
    if direction == 1:
        relocation = end + distance
        genome = genome[:start] + genome[end:relocation] + genome[start:end] + genome[relocation:]
    else:
        relocation = start - distance
        genome = genome[:relocation] + genome[start:end] + genome[relocation:start] + genome[end:]
    return genome

def add_relocations(genome,num):
    for i in range(num):
        size = random.randint(700,1500)
        distance = random.randint(5000,10000)
        genome = add_relocation(genome,size,distance)
    return genome
    
        
# Write genome to file with fasta header
def write_genome_to_file(genome,filepath):
    outfile = open(filepath + "contigs.fasta", "w")
    for i,contig in enumerate(genome):
        outfile.write(">contig" + str(i+1) + "\n")
        outfile.write(contig + "\n")

# Split up the genome into a set number of contigs
def split_to_contigs(genome, num_contigs):
    if num_contigs <= 1:
        return [genome]
    len_per_contig = len(genome) // num_contigs
    contigs = []
    for i in range(0,num_contigs - 1):
        contigs.append(genome[i*len_per_contig:(i+1)*len_per_contig])
    contigs.append(genome[(num_contigs-1)*len_per_contig:])
    return contigs

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

    # The volume errors argument lets the user tweak the amount of structural errors added to the genome after reads are drawn
    parser.add_argument(
        "-e", "--volume_errors",
        dest="volume_errors",
        default=1
    )

    parser.add_argument(
        "-c", "--num_contigs",
        dest="num_contigs",
        default=0, type=int
    )

    parser.add_argument(
        "-g", "--genome_size",
        dest="genome_size",
        default=5000000, type=int
    )

    args = parser.parse_args()
    e = int(args.volume_errors)

    # Create a random genome
    genome = random_genome(args.genome_size)

    # Draw reads from it
    reads_from_genome(
        genome=genome,
        coverage=20,
        mean_insert=300,
        sd_insert=50,
        seed=args.seed,
        outfile=args.file_path,
        num_low_coverage_regions=e*3)
    
    # Add errors to the genome so that when the reads are aligned, structural errors are detected
    genome = add_deletions(genome,e*10,500,1000)
    genome = add_inversions(genome,e*10,500,1000)
    genome = add_nearby_insertions(genome,e*8,500,1000,500)
    genome = add_tandem_insertions(genome,e*8,500,1000)
    genome = add_duplications(genome,e*10,500,1000)
    genome = add_relocations(genome,e*3)

    genome = split_to_contigs(genome, args.num_contigs)
    
    write_genome_to_file(genome,args.file_path)
    

if __name__=="__main__":
    main()