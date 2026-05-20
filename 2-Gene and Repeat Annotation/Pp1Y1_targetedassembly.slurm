#!/bin/bash
#SBATCH --job-name=Y1reconstruction
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=2:00:00
#SBATCH -o Y1fix.%j.out
#SBATCH -e Y1fix.%j.err

module purge
module load GCC/13.2.0 minimap2/2.29 SAMtools/1.21 seqtk/1.4

echo "Starting mapping step"

minimap2 -t 8 -ax map-hifi Pp1-Y1.fasta ISO1HiFireads.fastq.gz | \
samtools view -@ 8 -bS -F 2308 -q 60 | \
samtools sort -@ 8 -o targeted_Pp1Y1.bam

samtools index targeted_Pp1Y1.bam

echo "Mapping complete"

echo "Extracting mapped reads"

samtools view targeted_Pp1Y1.bam | cut -f1 | sort | uniq > mapped_reads.txt

seqtk subseq ISO1HiFireads.fastq.gz mapped_reads.txt > pp1y1_baited.fastq

echo "Read extraction complete"

echo "Starting Flye assembly"
module purge
module load GCC/13.3.0 Flye/2.9.5
flye \
  --pacbio-hifi pp1y1_baited.fastq \
  --out-dir pp1y1_assembly \
  --genome-size 5k \
  --threads 8

echo "Assembly complete"
