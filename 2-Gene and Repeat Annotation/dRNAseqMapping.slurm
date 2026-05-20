#!/bin/bash
#SBATCH --job-name=A4_Seq
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH -o A4.%j.out
#SBATCH -e A4.%j.err

module purge
module load GCC/13.2.0 minimap2/2.29 SAMtools/1.21

# Index your assembly first (optional but faster)
minimap2 -d ISO1_assembly.mmi DelgadoYISO-1v11.fasta

# Run the alignment
# -ax splice: Spliced alignment (RNA to Genome)
minimap2 -ax splice -uf -k14 -t 32 ISO1_assembly.mmi A4_F_wholebodyRNA_ONT.fq.gz > A4F_vs_ISO1.sam

samtools view -bS A4F_vs_ISO1.sam | samtools sort -o A4F_vs_ISO1_sorted.bam
samtools index A4F_vs_ISO1_sorted.bam
