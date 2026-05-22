#!/bin/bash
#SBATCH --job-name=winnowmap_XY_cov
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=164G
#SBATCH --time=48:00:00
#SBATCH --output=winnowmap_XY_cov_%j.out
#SBATCH --error=winnowmap_XY_cov_%j.err

set -euo pipefail

# ==============================================================================
# CONFIGURATION
# ==============================================================================

declare -A READS=(
    ["ISO1"]="/scratch/group/chakraborty_lab/adolfo/fertilitygenes/readconfirmation4codingsequences/ISO1HiFireads.fastq.gz"
    ["A3"]="/scratch/group/chakraborty_lab/adolfo/fertilitygenes/readconfirmation4codingsequences/A3HiFireads.fastq.gz"
    ["A4"]="/scratch/group/chakraborty_lab/adolfo/fertilitygenes/readconfirmation4codingsequences/A4HiFireads.fastq.gz"
)

declare -A Y_REFS=( ["ISO1"]="dmelISO1Y.fasta" ["A3"]="dmelA3Y.fasta" ["A4"]="dmelA4Y.fasta" )
declare -A Y_SCAFFOLDS=( ["ISO1"]="Scaffold_2" ["A3"]="scaffold_2" ["A4"]="scaffold_2" )
declare -A X_REFS=( ["ISO1"]="ISO1HarshX.fasta" ["A3"]="A3HarshX.fasta" ["A4"]="A4HarshX.fasta" )

THREADS=48
KMER_SIZE=21
HASH_SIZE="8G"
SUMMARY="XY_array_coverage_summary.tsv"

echo -e "Strain\tTarget_BAM\tContig\tArray_Size(bp)\tMean_Depth\tMedian_Depth\tMin_Depth\tMax_Depth\tBreadth(>0X)" > "$SUMMARY"

echo "=============================================================================="
echo " Starting Winnowmap XY & Coverage Pipeline on $(hostname) at $(date)"
echo "=============================================================================="

# ==============================================================================
# COVERAGE FUNCTION
# ==============================================================================
# We define your coverage logic as a function so we can easily call it 
# for the X, Y, and rDNA BAM files without repeating code.

calculate_coverage() {
    local STRAIN=$1
    local BAM_FILE=$2
    local TARGET_DESC=$3
    
    local DEPTH_OUT="${STRAIN}_${TARGET_DESC}_raw_depth.tsv"
    local TRACK_OUT="${STRAIN}_${TARGET_DESC}_10kb_coverage.bedGraph"

    if [ ! -f "$BAM_FILE" ]; then
        echo ">>> Warning: Cannot find $BAM_FILE. Skipping coverage for $TARGET_DESC..."
        return
    fi

    echo "  -> Generating raw depth file: $DEPTH_OUT..."
    samtools depth -a "$BAM_FILE" > "$DEPTH_OUT"

    echo "  -> Calculating summary statistics..."
    awk -v strain="$STRAIN" -v target="$TARGET_DESC" '
    BEGIN { OFS="\t"; }
    {
        contig = $1; d = $3;
        
        if (!(contig in count)) {
            count[contig] = 0; sum[contig] = 0;
            min[contig] = 9999999; max[contig] = 0;
            cov_gt_0[contig] = 0;
        }
        
        # Store depth for median calculation (grouped by contig)
        depths[contig, count[contig]] = d;
        count[contig]++;
        sum[contig] += d;
        if (d < min[contig]) min[contig] = d;
        if (d > max[contig]) max[contig] = d;
        if (d > 0) cov_gt_0[contig]++;
    }
    END {
        for (c in count) {
            n = count[c];
            if (n > 0) {
                mean = sum[c] / n;
                breadth = (cov_gt_0[c] / n) * 100;
                
                # Extract and sort depths for this specific contig
                delete temp_d;
                for (i = 0; i < n; i++) temp_d[i+1] = depths[c, i];
                asort(temp_d);
                
                if (n % 2 == 1) { median = temp_d[int(n/2) + 1]; } 
                else { median = (temp_d[n/2] + temp_d[n/2 + 1]) / 2; }
                
                printf "%s\t%s\t%s\t%d\t%.2f\t%.2f\t%d\t%d\t%.2f%%\n", strain, target, c, n, mean, median, min[c], max[c], breadth;
            }
        }
    }' "$DEPTH_OUT" >> "$SUMMARY"

    echo "  -> Generating 10kb rolling window BedGraph for IGV..."
    awk -v w=10000 '
    BEGIN { sum=0; count=0; start=1; prev_contig=""; }
    {
        curr_contig = $1;
        # Reset window if we switch to a new contig
        if (curr_contig != prev_contig && prev_contig != "") {
            if (count > 0) print prev_contig"\t"start-1"\t"prev_pos"\t"sum/count;
            sum=0; count=0; start=1;
        }
        
        sum += $3; count++; prev_pos = $2; prev_contig = curr_contig;
        
        if (count == w) {
            print curr_contig"\t"start-1"\t"$2"\t"sum/w;
            sum=0; count=0; start=$2+1;
        }
    }
    END {
        if (count > 0) print prev_contig"\t"start-1"\t"prev_pos"\t"sum/count;
    }' "$DEPTH_OUT" > "$TRACK_OUT"
}

# ==============================================================================
# MAIN PIPELINE LOOP
# ==============================================================================

for STRAIN in "ISO1" "A3" "A4"; do
    echo ""
    echo "#########################################################"
    echo " PROCESSING STRAIN: $STRAIN"
    echo "#########################################################"

    HIFI_READS="${READS[$STRAIN]}"
    
    # --- Y CHROMOSOME CONFIG ---
    Y_REF="${Y_REFS[$STRAIN]}"
    TARGET_SCAFFOLD="${Y_SCAFFOLDS[$STRAIN]}"
    Y_MERYL_DB="${STRAIN}_Y_merylDB"
    Y_REP_KMERS="${STRAIN}_Y_repetitive_k15.txt"
    Y_BAM_ALL="${STRAIN}_Y_WholeChr_Mapped.bam"
    Y_BAM_RDNA="${STRAIN}_Y_rDNA_Only.bam"
    Y_FQ_OUT="${STRAIN}_Y_rDNA_reads.fastq"
    Y_JF_OUT="${STRAIN}_Y_rDNA_k${KMER_SIZE}.jf"

    # --- X CHROMOSOME CONFIG ---
    X_REF="${X_REFS[$STRAIN]}"
    X_MERYL_DB="${STRAIN}_X_merylDB"
    X_REP_KMERS="${STRAIN}_X_repetitive_k15.txt"
    X_BAM_ALL="${STRAIN}_X_WholeChr_Mapped.bam"

    # --------------------------------------------------------------------------
    # PART 1: Y CHROMOSOME MAPPING & rDNA ISOLATION
    # --------------------------------------------------------------------------
    echo "--- PART 1: Y Chromosome ($STRAIN) ---"
    
    module purge
    module load GCC/13.2.0 Meryl/1.4.1
    echo "[Y-1] Building Meryl database..." 
    meryl count k=15 output "$Y_MERYL_DB" "$Y_REF"
    meryl print greater-than distinct=0.9998 "$Y_MERYL_DB" > "$Y_REP_KMERS"

    module purge
    module load GCC/12.3.0 Winnowmap/2.03 SAMtools/1.18
    echo "[Y-2] Mapping HiFi reads to $Y_REF..."
    winnowmap -W "$Y_REP_KMERS" -t "$THREADS" -ax map-pb --secondary=no "$Y_REF" "$HIFI_READS" | \
        samtools view -@ "$THREADS" -bS -F 2308 - | \
        samtools sort -@ "$THREADS" -o "$Y_BAM_ALL"
    samtools index "$Y_BAM_ALL"

    echo "[Y-3] Isolating reads mapping specifically to $TARGET_SCAFFOLD..."
    samtools view -@ "$THREADS" -b "$Y_BAM_ALL" "$TARGET_SCAFFOLD" > "$Y_BAM_RDNA"
    samtools index "$Y_BAM_RDNA"

    echo "[Y-4] Extracting mapped rDNA reads and running Jellyfish..."
    samtools fastq -@ "$THREADS" "$Y_BAM_RDNA" > "$Y_FQ_OUT"
    
    module purge
    module load GCC/13.2.0 Jellyfish/2.3.1
    jellyfish count -C -m "$KMER_SIZE" -s "$HASH_SIZE" -t "$THREADS" -o "$Y_JF_OUT" "$Y_FQ_OUT"
    jellyfish histo -t "$THREADS" "$Y_JF_OUT" > "${STRAIN}_Y_rDNA_k${KMER_SIZE}.histo"

    # --------------------------------------------------------------------------
    # PART 2: X CHROMOSOME MAPPING
    # --------------------------------------------------------------------------
    echo "--- PART 2: X Chromosome ($STRAIN) ---"
    
    module purge
    module load GCC/13.2.0 Meryl/1.4.1
    echo "[X-1] Building Meryl database..." 
    meryl count k=15 output "$X_MERYL_DB" "$X_REF"
    meryl print greater-than distinct=0.9998 "$X_MERYL_DB" > "$X_REP_KMERS"

    module purge
    module load GCC/12.3.0 Winnowmap/2.03 SAMtools/1.18
    echo "[X-2] Mapping HiFi reads to $X_REF..."
    winnowmap -W "$X_REP_KMERS" -t "$THREADS" -ax map-pb --secondary=no "$X_REF" "$HIFI_READS" | \
        samtools view -@ "$THREADS" -bS -F 2308 - | \
        samtools sort -@ "$THREADS" -o "$X_BAM_ALL"
    samtools index "$X_BAM_ALL"

    # --------------------------------------------------------------------------
    # PART 3: COVERAGE ANALYSIS
    # --------------------------------------------------------------------------
    echo "--- PART 3: Coverage Analysis ($STRAIN) ---"
    
    # Reload the specific SAMtools version you requested for coverage
    module purge
    module load GCC/13.2.0 SAMtools/1.21

    echo "[Cov-1] Analyzing Y Chromosome (Whole)..."
    calculate_coverage "$STRAIN" "$Y_BAM_ALL" "Y_WholeChr"

    echo "[Cov-2] Analyzing Y Chromosome (rDNA Scaffold Only)..."
    calculate_coverage "$STRAIN" "$Y_BAM_RDNA" "Y_rDNA_Target"

    echo "[Cov-3] Analyzing X Chromosome..."
    calculate_coverage "$STRAIN" "$X_BAM_ALL" "X_WholeChr"

done

echo "=============================================================================="
echo "Pipeline complete at $(date)!"
echo "Master coverage summary saved to: $SUMMARY"
echo "=============================================================================="
