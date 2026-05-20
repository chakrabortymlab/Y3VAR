#!/bin/bash
#SBATCH --job-name=map_winnow
#SBATCH --output=map_winnow_%j.log
#SBATCH --error=map_winnow_%j.err
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=24:00:00

# --- variables ---
REF="DelgadoYA3.fasta"
READS="SRR29479666.fastq.gz"
THREADS=32

KMER=28
MERYL_DB="${READS%.*}.k${KMER}.meryl"
REP_KMERS="${READS%.*}.rep.k${KMER}.txt"
DISTINCT="0.9998"

# --- build meryl DB from reads and produce repetitive k-mer list ---
module purge
module load GCC/13.2.0 Meryl/1.4.1
meryl threads="$THREADS" count k="$KMER" "$REF" output "$MERYL_DB"
meryl print greater-than distinct="$DISTINCT" "$MERYL_DB" > "$REP_KMERS"

module purge

# --- winnowmap mapping to SAM ---
module load GCC/12.2.0 Winnowmap/2.03
winnowmap -k 28 -t "$THREADS" -W "$REP_KMERS" -ax map-pb "$REF" "$READS" > alignment.sam

module purge

# --- convert, sort, index with samtools ---
module load GCC/13.2.0 SAMtools/1.21
samtools view -@ "$THREADS" -Sb alignment.sam | samtools sort -@ "$THREADS" -o mappedreads.bam
samtools index mappedreads.bam

# --- cleanup ---
rm -f alignment.sam
