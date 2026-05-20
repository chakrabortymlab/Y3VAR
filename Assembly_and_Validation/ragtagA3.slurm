#!/bin/bash
#SBATCH --job-name=ragtag
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH -o ragtag.%j.out
#SBATCH -e ragtag.%j.err

module load GCC/13.2.0 OpenMPI/4.1.6 RagTag/2.1.0

ragtag.py scaffold DelgadoY.fasta YA3.fasta -o ragtag --mm2-params "-x asm20"
