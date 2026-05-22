#!/bin/bash
#SBATCH --job-name=rDNA_filter_subsample
#SBATCH --time=04:00:00                 
#SBATCH --ntasks-per-node=1                       
#SBATCH --cpus-per-task=48               
#SBATCH --mem=168G                      
#SBATCH --output=rDNA_filter_%j.out     
#SBATCH --error=rDNA_filter_%j.err      

set -euo pipefail

# ==========================================
# ENVIRONMENT & VARIABLES
# ==========================================
module purge
module load GCC/12.3.0 Winnowmap/2.03 SAMtools/1.18 seqtk/1.4

REF_FASTA="A3rDNAarray.fasta"
READS_IN="../A3HiFireads.fastq.gz"
OUTPUT_PREFIX="clean_rDNA"

THREADS=$SLURM_CPUS_PER_TASK

# Filtering parameters
MIN_LENGTH=10000

# Subsampling parameters
COVERAGES=(5 15 30)      # Array of target coverages
TARGET_REGION_SIZE=1196040     # Size of the locus in base pairs (~1.2 Mbp). Change to match individual array. 

echo "=========================================="
echo " Starting Read Filtering & Subsampling Pipeline"
echo " Target Coverages: ${COVERAGES[*]}X across ${TARGET_REGION_SIZE} bp"
echo "=========================================="

# 1. Map reads to the main array
echo " [1/5] Mapping reads to reference units using winnowmap2..."
winnowmap -ax map-pb -t "$THREADS" "$REF_FASTA" "$READS_IN" > "${OUTPUT_PREFIX}_mapped.sam"

# 2. Filter strictly for Primary Hits
echo " [2/5] Filtering for primary alignments only..."
samtools view -bS -F 2308 -@ "$THREADS" "${OUTPUT_PREFIX}_mapped.sam" > "${OUTPUT_PREFIX}_primary.bam"

# 3. Convert back to FASTQ
echo " [3/5] Extracting mapped reads to FASTQ..."
samtools fastq -@ "$THREADS" "${OUTPUT_PREFIX}_primary.bam" > "${OUTPUT_PREFIX}_primary.fastq"

# 4. Filter by minimum read length
echo " [4/5] Filtering out reads shorter than ${MIN_LENGTH} bp..."
seqtk seq -L "$MIN_LENGTH" "${OUTPUT_PREFIX}_primary.fastq" > "${OUTPUT_PREFIX}_long_reads.fastq"

# 5. Dynamically subsample to Target Coverages
echo " [5/5] Counting total available bases in filtered reads..."
TOTAL_BASES=$(awk 'NR%4==2{sum+=length($0)}END{print sum}' "${OUTPUT_PREFIX}_long_reads.fastq")
echo "       Available bases in long reads: $TOTAL_BASES"

echo "=========================================="
echo " Subsampling Loop Started"
echo "=========================================="

for TARGET_COVERAGE in "${COVERAGES[@]}"; do
    echo " -> Processing ${TARGET_COVERAGE}X coverage..."
    
    # Calculate how many total bases we need for this specific coverage
    TARGET_BASES=$(echo "$TARGET_COVERAGE * $TARGET_REGION_SIZE" | bc)
    echo "    Target bases needed: $TARGET_BASES"
    
    # Calculate the fraction of reads to keep
    FRACTION=$(echo "scale=4; $TARGET_BASES / $TOTAL_BASES" | bc)
    FINAL_OUT="${OUTPUT_PREFIX}_${TARGET_COVERAGE}x_subsampled.fastq"
    
    # Check if we already have less data than the target
    if (( $(echo "$FRACTION >= 1.0" | bc -l) )); then
        echo "    Available coverage is less than or equal to target. Copying all long reads to this dataset."
        cp "${OUTPUT_PREFIX}_long_reads.fastq" "$FINAL_OUT"
    else
        echo "    Subsampling to ${FRACTION} fraction of reads (Seed: 42)..."
        seqtk sample -s 42 "${OUTPUT_PREFIX}_long_reads.fastq" "$FRACTION" > "$FINAL_OUT"
    fi
    echo "    Saved to: $FINAL_OUT"
    echo "------------------------------------------"
done

# Cleanup intermediate mapping files to save space
rm "${OUTPUT_PREFIX}_mapped.sam"
rm "${OUTPUT_PREFIX}_primary.bam"
rm "${OUTPUT_PREFIX}_primary.fastq"

echo "=========================================="
echo " Pipeline Complete!"
echo " Your coverage-specific FASTQ files have been generated."
echo "=========================================="
