#!/bin/bash
#SBATCH --job-name=morph_abundance_estimation
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=morph_abundance_%j.out

set -euo pipefail

# ==============================================================================
# CONFIGURATION
# ==============================================================================

# Your single, unified long-read dataset, not subsampled to a specific coverage
READS_FQ="clean_rDNA_long_reads.fastq"

# Define your datasets: Prefix : Coverage : Morph FASTA : RepeatMasker .out
DATASETS=(
    "5x:30:./5xRef_output/5xmorphs.fa:./5xRef_output/5xmorphs.fa.out"
    "15x:30:./15xRef_output/15xmorphs.fa:./15xRef_output/15xmorphs.fa.out"
    "30x:30:./30xRef_output/30xmorphs.fa:./30xRef_output/30xmorphs.fa.out"
)

# Output directory
OUT_DIR="Array_Estimates"
mkdir -p "$OUT_DIR"
THREADS="${SLURM_CPUS_PER_TASK:-32}"

module purge
module load GCC/13.2.0 minimap2/2.29 Biopython/1.84 SAMtools/1.21

echo "=============================================================================="
echo " Starting Unique-Sequence Anchored Morph Abundance Estimation"
echo " Single Readset: $READS_FQ"
echo "=============================================================================="

if [[ ! -f "$READS_FQ" ]]; then
    echo ">>> FATAL ERROR: $READS_FQ not found in the current directory."
    exit 1
fi

for DATASET in "${DATASETS[@]}"; do
    IFS=':' read -r PREFIX COVERAGE MORPH_FA RM_OUT <<< "$DATASET"
    
    if [[ ! -f "$MORPH_FA" ]] || [[ ! -f "$RM_OUT" ]]; then 
        echo ">>> WARNING: Missing FASTA or RM .out for $PREFIX. Skipping."
        continue
    fi

    echo ">>> Processing $PREFIX dataset (Coverage Baseline: ${COVERAGE}x)"
    
    PAF_OUT="${OUT_DIR}/${PREFIX}_reads_to_morphs.paf"
    SAM_RAW="${OUT_DIR}/${PREFIX}_raw.sam"
    SAM_CLEAN="${OUT_DIR}/${PREFIX}_winning_reads.sam"
    BAM_FINAL="${OUT_DIR}/${PREFIX}_winning_reads.bam"
    TSV_OUT="${OUT_DIR}/${PREFIX}_Array_Abundance.tsv"

    # 1. Map to PAF (Used for the math and logic)
    minimap2 -c -t "$THREADS" -x map-hifi "$MORPH_FA" "$READS_FQ" > "$PAF_OUT"
    
    # 2. Map to SAM (Used to generate the IGV visualizer)
    minimap2 -a -t "$THREADS" -x map-hifi "$MORPH_FA" "$READS_FQ" > "$SAM_RAW"

    # 3. Python Engine: Enforce Rules & Filter SAM
    python3 -c "
import sys
import csv
from collections import defaultdict

paf_file = sys.argv[1]
coverage = float(sys.argv[2])
tsv_out = sys.argv[3]
rm_file = sys.argv[4]
sam_in = sys.argv[5]
sam_out = sys.argv[6]

# Step A: Parse RepeatMasker to establish the 'Unique Window' (28S start -> ETS start)
morph_unique_windows = {}
try:
    with open(rm_file, 'r') as f:
        lines = f.readlines()[3:]
        coords = defaultdict(lambda: {'28S_starts': [], 'ETS_starts': []})
        
        for line in lines:
            parts = line.strip().split()
            if len(parts) < 10: continue
            
            m_name = parts[4]
            start_coord = int(parts[5])
            rep_name = parts[9]
            
            if '28S' in rep_name:
                coords[m_name]['28S_starts'].append(start_coord)
            elif 'ETS' in rep_name:
                coords[m_name]['ETS_starts'].append(start_coord)
                
        for m_name, data in coords.items():
            if not data['28S_starts'] or not data['ETS_starts']:
                morph_unique_windows[m_name] = (0, float('inf'))
                continue
                
            first_28s = min(data['28S_starts'])
            first_ets = min(data['ETS_starts'])
            
            var_start = min(first_28s, first_ets)
            var_end = max(first_28s, first_ets)
                
            morph_unique_windows[m_name] = (var_start, var_end)
            
except Exception as e:
    print(f'Error parsing RM file: {e}')
    sys.exit(1)

# =====================================================================
# Step A.5: Export the Unique Variable Regions to an IGV-ready BED file
# =====================================================================
variable_bed_path = tsv_out.replace('_Array_Abundance.tsv', '_Variable_Regions.bed')
with open(variable_bed_path, 'w') as bf:
    bf.write('track name=Hyper-Variable_Regions description=28S_to_ETS_Window color=128,0,128\n')
    for m_name, (v_start, v_end) in morph_unique_windows.items():
        if v_end != float('inf'):
            clean_m = m_name.split('_')[0]
            if 'morphconsensus' in m_name:
                clean_m = 'Morph_' + m_name.split('morphconsensus')[1].split('_')[0]
            
            bf.write(f'{m_name}\t{v_start}\t{v_end}\t{clean_m}_Variable_Region\t0\t+\n')

# Step B: Parse PAF & Enforce Unique Window + Identity Rules
read_best_match = {}
morph_lengths = {}

with open(paf_file, 'r') as f:
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) < 12: continue
        
        q_name = parts[0]           
        t_name = parts[5]           
        t_len = int(parts[6])       
        t_start = int(parts[7])
        t_end = int(parts[8])
        
        matching_bases = int(parts[9])
        align_block_len = int(parts[10])
        
        morph_lengths[t_name] = t_len
        
        # Rule Check 1: Calculate Sequence Identity
        seq_identity = matching_bases / align_block_len
        if seq_identity < 0.95:
            continue
        
        # Rule Check 2: Calculate Physical Variable Region Coverage
        var_start, var_end = morph_unique_windows.get(t_name, (0, t_len))
        if var_end == float('inf'): var_end = t_len
        
        var_len = var_end - var_start
        if var_len <= 0: var_len = t_len 
        
        overlap_start = max(t_start, var_start)
        overlap_end = min(t_end, var_end)
        overlap = max(0, overlap_end - overlap_start)
        
        unique_cov = overlap / var_len
        
        if unique_cov < 0.95:
            continue
            
        score = int(parts[9])
        for tag in parts[12:]:
            if tag.startswith('AS:i:'):
                score = int(tag.split(':')[2])
                break
        
        # Highlander Rule (Best Score Wins)
        if q_name not in read_best_match:
            read_best_match[q_name] = {'morph': t_name, 'score': score}
        else:
            if score > read_best_match[q_name]['score']:
                read_best_match[q_name] = {'morph': t_name, 'score': score}

# Step C: Create a valid (Read, Morph) Set and Filter the SAM File
valid_pairs = set()
for q_name, data in read_best_match.items():
    valid_pairs.add((q_name, data['morph']))

written_reads = set()
with open(sam_in, 'r') as s_in, open(sam_out, 'w') as s_out:
    for line in s_in:
        # Keep all headers
        if line.startswith('@'):
            s_out.write(line)
            continue
            
        parts = line.split('\t')
        if len(parts) < 3: continue
        
        q_name = parts[0]
        flag = int(parts[1])
        t_name = parts[2]
        
        # Drop secondary alignments (256) and supplementary alignments (2048)
        if (flag & 256) or (flag & 2048):
            continue
        
        # Write ONLY the single primary alignment for the verified winning morph
        if (q_name, t_name) in valid_pairs and q_name not in written_reads:
            s_out.write(line)
            written_reads.add(q_name)

# Step D: Tally and Math for TSV
morph_read_counts = defaultdict(int)
for q_name, data in read_best_match.items():
    morph_read_counts[data['morph']] += 1

output_data = []
total_array_bp = 0
total_copy_number = 0

for morph, length in morph_lengths.items():
    count = morph_read_counts.get(morph, 0)
    copy_number = count / coverage
    estimated_bp = copy_number * length
    
    total_copy_number += copy_number
    total_array_bp += estimated_bp
    
    clean_morph = morph.split('_')[0]
    if 'morphconsensus' in morph:
        clean_morph = 'Morph_' + morph.split('morphconsensus')[1].split('_')[0]
        
    output_data.append([clean_morph, length, count, round(copy_number, 2), round(estimated_bp, 2)])

output_data.sort(key=lambda x: x[3], reverse=True)

with open(tsv_out, 'w', newline='') as f:
    writer = csv.writer(f, delimiter='\t')
    writer.writerow(['Morph_ID', 'Length_bp', 'Reads_Mapped_Perfectly', 'Estimated_Copy_Number', 'Estimated_Total_bp'])
    for row in output_data:
        writer.writerow(row)
        
    writer.writerow([])
    writer.writerow(['TOTAL_ESTIMATED_ARRAY', '-', sum(row[2] for row in output_data), round(total_copy_number, 2), round(total_array_bp, 2)])

" "$PAF_OUT" "$COVERAGE" "$TSV_OUT" "$RM_OUT" "$SAM_RAW" "$SAM_CLEAN"

    # 4. Convert Clean SAM to BAM, Sort, and Index
    echo "    - Compiling IGV BAM file..."
    samtools view -bS "$SAM_CLEAN" | samtools sort -@ "$THREADS" -o "$BAM_FINAL"
    samtools index "$BAM_FINAL"
    
    # 5. Cleanup temporary mapping files to save storage space
    rm "$PAF_OUT" "$SAM_RAW" "$SAM_CLEAN"

    echo "  -> Completed $PREFIX. Results saved to $TSV_OUT & $BAM_FINAL"
done

echo "=============================================================================="
echo " Array Estimation Complete!"
echo "=============================================================================="
