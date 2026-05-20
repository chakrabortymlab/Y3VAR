#!/bin/bash
#SBATCH --job-name=Y_depth
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH -o y_inclusive.%j.out
#SBATCH -e y_inclusive.%j.err

module purge
module load GCC/13.3.0 SAMtools/1.21

samtools depth -a mappedreads.bam > raw_coverage_depth.txt
