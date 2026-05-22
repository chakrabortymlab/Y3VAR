#!/bin/bash
#SBATCH --job-name=CGcontig_map
#SBATCH --output=CGcontig_map_%j.out
#SBATCH --error=CGcontig_map_%j.err
#SBATCH --time=02:00:00           # 2 hours should be plenty for a 1kb contig
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16         # 8 threads for minimap2 and samtools
#SBATCH --mem=128G                 # 16 GB RAM

module purge
module load GCC/13.2.0 minimap2/2.29 SAMtools/1.21

# --- 2. Define Input Files ---
CONTIG="CGcontig.fasta"
MALE_READS="/scratch/user/aadelgado6/Delgado/A4/reads/SRR29479667.fastq"     # Update with your actual male reads file
FEMALE_READS="/scratch/group/chakraborty_lab/alex/dmel/single_fly/hifi_reads/A4B_hifi.fq" # Update with your actual female reads file

# --- 3. Define Output Files ---
MALE_BAM="male_CGcontig_strict.bam"
FEMALE_BAM="female_CGcontig_strict.bam"

MALE_DEPTH="male_CGcontig_depth.txt"
FEMALE_DEPTH="female_CGcontig_depth.txt"
SUMMARY="CGcontig_Strict_Coverage_Summary.txt"

echo "Starting STRICT HiFi Mapping to CGcontig..." > $SUMMARY
echo "Date: $(date)" >> $SUMMARY
echo "------------------------------------------------" >> $SUMMARY

# --- 4. STRICT Mapping for Male Reads ---
echo "Mapping Male HiFi reads (STRICT)..."
# minimap2 restricts to primary only. Samtools filters out unmapped (-F 4), 
# secondary/supplementary (-F 2304), and requires perfect mapping quality (-q 60).
minimap2 -t 8 -ax map-hifi --secondary=no $CONTIG $MALE_READS | \
samtools view -@ 8 -b -h -F 2308 -q 60 - | \
samtools sort -@ 8 -o $MALE_BAM

samtools index $MALE_BAM
samtools depth -a $MALE_BAM > $MALE_DEPTH

# --- 5. STRICT Mapping for Female Reads ---
echo "Mapping Female HiFi reads (STRICT)..."
minimap2 -t 8 -ax map-hifi --secondary=no $CONTIG $FEMALE_READS | \
samtools view -@ 8 -b -h -F 2308 -q 60 - | \
samtools sort -@ 8 -o $FEMALE_BAM

samtools index $FEMALE_BAM
samtools depth -a $FEMALE_BAM > $FEMALE_DEPTH

# --- 6. Quick Summary Stats ---
echo -e "\n--- STRICT Mapping Summary ---" >> $SUMMARY

MALE_AVG=$(awk '{sum+=$3} END { if (NR > 0) print sum/NR; else print 0 }' $MALE_DEPTH)
FEMALE_AVG=$(awk '{sum+=$3} END { if (NR > 0) print sum/NR; else print 0 }' $FEMALE_DEPTH)

echo "Average Male Coverage (Strict): ${MALE_AVG}X" >> $SUMMARY
echo "Average Female Coverage (Strict): ${FEMALE_AVG}X" >> $SUMMARY

echo "------------------------------------------------" >> $SUMMARY
echo "Strict Mapping Complete!" >> $SUMMARY

cat $SUMMARY
