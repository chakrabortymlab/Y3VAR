#!/bin/bash
#SBATCH --job-name=BLAT_Y_CDS
#SBATCH --output=BLAT_CDS_%j.log
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --ntasks-per-node=1

module purge
module load GCC/13.3.0 BLAT/3.7

# 1. DEFINE REFERENCE
# We use the new Carvalho et al. 2026 CDS reference
REF_CDS="CarvalhoY_CDS.fasta"

# 2. RUN BLAT FOR EACH ASSEMBLY
# -t=dna: Target (Assembly) is DNA
# -q=dna: Query (CDS) is DNA
# -maxIntron=2000000: Essential for bridging large Y-linked introns (kl-3/kl-5)
# -fine: (Optional but recommended) Increases sensitivity for small exons in CDS mapping
ASSEMBLIES=("DelgadoYA3.fasta" "DelgadoYA4.fasta" "DelgadoYISO-1v11.fasta")

for ASM in "${ASSEMBLIES[@]}"; do
    NAME=$(basename $ASM .fasta)
    echo "Processing $NAME with reference $REF_CDS..."
    
    # Run BLAT
    blat "$ASM" "$REF_CDS" \
         -t=dna -q=dna \
         -maxIntron=2000000 \
         -fine \
         "results_${NAME}_CDS.psl"
done
