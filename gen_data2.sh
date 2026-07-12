#!/bin/bash

#!/bin/bash

### Script to generate data that compares results between valet and frcBAM -> generates sample data with many small contigs

### Example usage: chmod +x gen_data.sh & ./gen_data 1,2,3,4,5

if [ -z "$1" ]; then
    echo "Usage: $0 1,2,3,4"
    exit 1
fi

NUMBERS=$(echo "$1" | tr ',' ' ')

OUTDIR="$(realpath ./test_data/contigs)"
RESULTS="$(realpath ./test_data/contigs/results.tsv)"
RESULTS_FRC="$(realpath ./test_data/contigs/frc_results.tsv)"

for n in $NUMBERS; do

    TEST_DIR="${OUTDIR}/${n}"

    echo "=============================="
    echo "Running test: ${n}"
    echo "=============================="

    mkdir -p "${TEST_DIR}"

    python src/py/test_data_generator_v2.py \
        -s "${n}" \
        -e 15 \
        -c 18 \
        -o "${TEST_DIR}/" \
        -g 300000
    
    bwa index "${TEST_DIR}/contigs.fasta"

    bwa mem "${TEST_DIR}/contigs.fasta" "${TEST_DIR}/R1.fastq" "${TEST_DIR}/R2.fastq" > "${TEST_DIR}/alignment.sam"

    samtools view -bS "${TEST_DIR}/alignment.sam" > "${TEST_DIR}/alignment.bam"

    samtools sort "${TEST_DIR}/alignment.bam" -o "${TEST_DIR}/alignment.sorted.bam"

    samtools index "${TEST_DIR}/alignment.sorted.bam"

    mkdir -p "${TEST_DIR}/frc"

    ~/FRC_align/bin/FRC --pe-sam "${TEST_DIR}/alignment.sorted.bam" --output "${TEST_DIR}/frc/t"

    nextflow nf/valet.nf \
        --assembly "${TEST_DIR}/contigs.fasta" \
        --reads1 "${TEST_DIR}/R1.fastq" \
        --reads2 "${TEST_DIR}/R2.fastq" \
        -output-dir "${TEST_DIR}"

    LABEL="r${n}"

    echo "VALET RESULTS:"

    python test_detection.py \
        -g "${TEST_DIR}/errors.bed" \
        -v "${TEST_DIR}/summary/summary.bed" \
        -l "${LABEL}" \
        -f "${RESULTS}"

    echo "FRC RESULTS:"

    python frc_detection.py \
        -g "${TEST_DIR}/errors.bed" \
        -v "${TEST_DIR}/frc/t_Features.txt" \
        -l "${LABEL}" \
        -f "${RESULTS_FRC}"

done