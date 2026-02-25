import random
import random
from pathlib import Path


def random_genome(length):
    return ''.join(random.choices('ACGT', k=length))

def reverse_complement(seq):
    comp = str.maketrans("ACGT", "TGCA")
    return seq.translate(comp)[::-1]


def reads_from_genome(genome,coverage,mean_insert,sd_insert,seed,outfile,
                      mate_error_regions,mate_error_region_size,mate_error_support):
    
    random.seed(seed)

    genome_length = len(genome)
    num_fragments = int((coverage * genome_length) / mean_insert)

    r1_file = open(outfile + "R1.fastq", "w")
    r2_file = open(outfile + "R2.fastq", "w")
    truth_file = open(outfile + "_truth.tsv", "w")

    for i in range(num_fragments):
        # Generate a random insert size using a gaussian distribution about the mean
        insert_size = int(random.gauss(mean_insert, sd_insert))
        read_length = insert_size / 2 - 20

        # Determine the start and end of the read
        start = random.randint(0, genome_length - insert_size)
        end = start + insert_size

        fragment = genome[start:end]

        forward_read = fragment[:read_length]
        reverse_read = reverse_complement(fragment[-read_length:])

        # Quality scores
        quality_scores = "I" * read_length

        read_id = f"@SIM:{i}:start={start}:insert={insert_size}"

        # Write to FASTQ files
        r1_file.write(f"{read_id}/1\n{forward_read}\n+\n{quality_scores}\n")
        r2_file.write(f"{read_id}/2\n{reverse_read}\n+\n{quality_scores}\n")

        # Write ground truth
        truth_file.write(f"{i}\t{start}\t{end}\t{insert_size}\n")

    # Get mate error regions
    regions = []
    for i in range(mate_error_regions):
        start = random.randint(0, genome_length - mate_error_region_size)
        end = start + mate_error_region_size
        regions.append((start,end))

    # Populate mate error regions
    for (start, end) in regions:
        for i in range(mate_error_support):
            gene_start = random.randint(start, end - 50)
            gene = genome[gene_start:gene_start + 50]
            quality_scores = "I" * 50
            read_id = f"@Error:{i}:start={gene_start}"
            r1_file.write(f"{read_id}/1\n{forward_read}\n+\n{quality_scores}\n")

    r1_file.close()
    r2_file.close()
    truth_file.close()


genome = random_genome(5000000)
reads_from_genome(genome)