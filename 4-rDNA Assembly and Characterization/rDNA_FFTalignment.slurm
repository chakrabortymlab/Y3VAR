#!/bin/bash
#SBATCH --job-name=rdna_fftns2_backup
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=360G
#SBATCH --time=96:00:00
#SBATCH --output=rdna_fftns2_%j.out
#SBATCH --error=rdna_fftns2_%j.err

set -euo pipefail

# ==============================================================================
# CONFIGURATION
# ==============================================================================
INPUT_FILES=(
    "ISO1_X:./ISO1Xmorphs.fasta"
    "ISO1_Y:./ISO1Ymorphs.fasta"
    "A3_X:./A3Xmorphs.fasta"
    "A3_Y:./A3Ymorphs.fasta"
    "A4_X:./A4Xmorphs.fasta"
    "A4_Y:./A4Ymorphs.fasta"
)

OUT_DIR="Homogenization_FFTNS2_GAPSFIXED_MAXITERATE"
mkdir -p "$OUT_DIR"
THREADS="${SLURM_CPUS_PER_TASK:-8}"

COMBINED_FA="${OUT_DIR}/All_Strains_Combined_Units.fasta"
ALIGNED_FA="${OUT_DIR}/All_Strains_Aligned_FFTNS2.fasta"

module purge
module load GCC/13.2.0 MAFFT/7.526-with-extensions


echo ">>> Step 1: Merging units and standardizing headers via AWK..."
for item in "${INPUT_FILES[@]}"; do
    TAG="${item%%:*}"
    FILEPATH="${item##*:}"
    
    if [[ ! -f "$FILEPATH" ]]; then
        echo "    - Warning: $FILEPATH not found. Skipping."
        continue
    fi
    
    awk -v tag="$TAG" '
    BEGIN { count = 1 }
    /^>/ {
        raw_header = substr($0, 2)
        split(raw_header, parts, "::")
        num_splits = split(parts[1], subparts, "_Unit_")
        unit_type = (num_splits > 1) ? subparts[2] : "UnknownType"
        printf ">%s_%s_%d\n", tag, unit_type, count
        count++
        next
    }
    { print $0 }
    ' "$FILEPATH" >> "$COMBINED_FA"
done

echo ">>> Step 2: Running FFT-NS-2 Alignment..."
# --retree 2: builds the tree twice for better accuracy
# --maxiterate 1000: adds a layer of iterative refinement without the massive L-INS-i overhead
mafft --retree 2 --maxiterate 1000 --ep 0.0 --thread "$THREADS" "$COMBINED_FA" > "$ALIGNED_FA"

echo "Done. Alignment saved to $ALIGNED_FA"
