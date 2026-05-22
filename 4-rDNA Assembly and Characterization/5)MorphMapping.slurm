#!/bin/bash
#SBATCH --job-name=rdna_morph_verification
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=168G
#SBATCH --time=02:00:00
#SBATCH --output=rdna_verify_%j.out
#SBATCH --error=rdna_verify_%j.err

set -euo pipefail

# ==============================================================================
# CONFIGURATION
# ==============================================================================
# The reference should be the whole Y chromosome assembly
REF_FA="DelgadoYA3.fasta"

# Array of your successful ribotin-ref output directories for different coverages
COV_DIRS=("5xRef_output" "15xRef_output" "30xRef_output")

OUT_DIR="Morph_mapping_Results"
mkdir -p "$OUT_DIR"
TARGET_SCAFFOLD="Y_scaffold2"   #Main rDNA-containing scaffold
THREADS="${SLURM_CPUS_PER_TASK:-32}"

SUMMARY="${OUT_DIR}/Verification_Summary.txt"
echo -e "Dataset\tGold_Mapped(>98%)\tRelaxed_Mapped(>80%)\tTotal_Mapped" > "$SUMMARY"

module purge
module load GCC/13.2.0  minimap2/2.29 Biopython/1.84

echo "=============================================================================="
echo " Starting Dual-Pass (asm5 + asm20) Morph Verification"
echo " Target: $TARGET_SCAFFOLD"
echo "=============================================================================="

for DIR in "${COV_DIRS[@]}"; do
    DATASET=$(basename "$DIR")
    echo ">>> Processing Dataset: $DATASET"
    
    MORPH_FA="${DIR}/morphs.fa"
    
    if [[ ! -f "$MORPH_FA" ]]; then 
        echo ">>> WARNING: $MORPH_FA not found. Skipping dataset $DATASET."
        continue
    fi

    PAF_ASM5="${OUT_DIR}/${DATASET}_asm5.paf"
    PAF_ASM20="${OUT_DIR}/${DATASET}_asm20.paf"
    BED_FINAL="${OUT_DIR}/${DATASET}_verified_dual_pass.bed"

    # 1. Mapping - Pass 1 (Strict)
    minimap2 -c -t "$THREADS" -x asm5 -N 100000 "$REF_FA" "$MORPH_FA" > "$PAF_ASM5"

    # 2. Mapping - Pass 2 (Relaxed)
    minimap2 -c -t "$THREADS" -x asm20 -N 100000 "$REF_FA" "$MORPH_FA" > "$PAF_ASM20"

    # 3. Two-Pass Resolution via Python
    python3 -c "
import sys

def parse_paf(filename, target_filter, min_cov):
    alns = []
    with open(filename) as f:
        for line in f:
            p = line.strip().split('\t')
            if len(p) < 12 or p[5] != target_filter: continue
            qlen, qstart, qend = int(p[1]), int(p[2]), int(p[3])
            if (qend - qstart) / qlen < min_cov: continue
            # Handle missing AS tag safely
            score = 0
            for tag in p[12:]:
                if tag.startswith('AS:i:'):
                    score = int(tag.split(':')[2])
                    break
            alns.append({'score': score, 'ts': int(p[7]), 'te': int(p[8]), 'q': p[0], 's': p[4]})
    return alns

# Load both: Pass 1 is strict, Pass 2 is slightly more relaxed on coverage
strict_alns = parse_paf(sys.argv[1], sys.argv[4], 0.98)
relaxed_alns = parse_paf(sys.argv[2], sys.argv[4], 0.80)

strict_alns.sort(key=lambda x: x['score'], reverse=True)
relaxed_alns.sort(key=lambda x: x['score'], reverse=True)

accepted_coords = []
final_output = []

# Step A: Lock in Gold hits
for aln in strict_alns:
    if not any(max(aln['ts'], a[0]) < min(aln['te'], a[1]) for a in accepted_coords):
        accepted_coords.append((aln['ts'], aln['te']))
        final_output.append(f\"{sys.argv[4]}\t{aln['ts']}\t{aln['te']}\t{aln['q']}_Gold\t{aln['score']}\t{aln['s']}\")

# Step B: Fill gaps with Relaxed hits
for aln in relaxed_alns:
    if not any(max(aln['ts'], a[0]) < min(aln['te'], a[1]) for a in accepted_coords):
        accepted_coords.append((aln['ts'], aln['te']))
        final_output.append(f\"{sys.argv[4]}\t{aln['ts']}\t{aln['te']}\t{aln['q']}_Relaxed\t{aln['score']}\t{aln['s']}\")

with open(sys.argv[3], 'w') as out:
    for line in sorted(final_output, key=lambda x: int(x.split('\t')[1])):
        out.write(line + '\n')
" "$PAF_ASM5" "$PAF_ASM20" "$BED_FINAL" "$TARGET_SCAFFOLD"

    # 4. Append Results to Summary
    GOLD_COUNT=$(grep "_Gold" "$BED_FINAL" | wc -l || echo 0)
    RELAX_COUNT=$(grep "_Relaxed" "$BED_FINAL" | wc -l || echo 0)
    TOTAL_COUNT=$((GOLD_COUNT + RELAX_COUNT))
    
    echo -e "${DATASET}\t${GOLD_COUNT}\t${RELAX_COUNT}\t${TOTAL_COUNT}" >> "$SUMMARY"
    echo "  -> Finished $DATASET: $GOLD_COUNT Gold, $RELAX_COUNT Relaxed. Total: $TOTAL_COUNT"
done

echo "=============================================================================="
echo " Verification Complete. Summary file: $SUMMARY"
echo "=============================================================================="
