#!/bin/bash
#SBATCH --job-name=mask_primary
#SBATCH --output=mask_primary_%j.log
#SBATCH --error=mask_primary_%j.err
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=03:00:00

# --- Stop script on any error ---
set -euo pipefail

# --- Variables ---
INPUT_FASTA="DelgadoYA3.fasta"
INPUT_BAM="mappedreads.bam"         # The original raw BAM
FILTERED_BAM="mappedreads_primary.bam" # The new clean BAM

# Output Files
OUTPUT_MASKED_FASTA="DelgadoY_validated_masked.fasta"
COVERAGE_BED="primary_coverage_windows.bed"
THREADS=$SLURM_CPUS_PER_TASK

# --- Load Modules ---
module purge
module load GCC/13.3.0 BEDTools/2.31.1 SAMtools/1.21

echo "----------------------------------------------------"
echo "Step 0: Filtering BAM for Primary Alignments Only"
echo "----------------------------------------------------"

# -F 2308 excludes: Unmapped (4), Secondary (256), Supplementary (2048)
echo "Filtering secondary alignments from $INPUT_BAM..."
samtools view -@ "$THREADS" -b -F 2308 "$INPUT_BAM" > "$FILTERED_BAM"

echo "Indexing filtered BAM..."
samtools index "$FILTERED_BAM"

echo "Filtering complete. Using $FILTERED_BAM for all subsequent steps."

echo "----------------------------------------------------"
echo "Step 1: Calculating Zero Coverage & Masking Assembly"
echo "----------------------------------------------------"

# 1. Identify 0 coverage regions using the FILTERED BAM
# -bga gives genome regions. awk filters for 0 count.
bedtools genomecov -ibam "$FILTERED_BAM" -bga | awk 'BEGIN{OFS="\t"} $4 == 0 {print $1, $2, $3}' > zero_coverage.bed

# 2. Calculate stats for the log
MASK_COUNT=$(awk '{sum += $3-$2} END {print sum}' zero_coverage.bed)
echo "Masking $MASK_COUNT bases due to lack of primary read support."

# 3. Mask the FASTA
bedtools maskfasta -fi "$INPUT_FASTA" -bed zero_coverage.bed -fo "$OUTPUT_MASKED_FASTA"
echo "Masked assembly saved to: $OUTPUT_MASKED_FASTA"
