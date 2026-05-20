from Bio import SeqIO
from Bio.Seq import Seq

def load_exon_sequence(filename):
    """
    Reads a FASTA file and returns ONLY the nucleotide sequence as a string.
    This ignores the header (>...) automatically.
    """
    try:
        record = SeqIO.read(filename, "fasta")
        return str(record.seq).upper()
    except ValueError:
        # If your exon files have multiple records, this captures them all as one string
        records = SeqIO.parse(filename, "fasta")
        return "".join(str(r.seq) for r in records).upper()

def patch_assembly():
    # 1. Load the assembly (v1.1 includes your Pp1Y1 contig)
    # Using a list of records is safer for Biopython 1.84+
    print("Loading assembly...")
    assembly_records = list(SeqIO.parse("DelgadoYISO-1v2.fasta", "fasta"))
    
    # Map IDs to their index in the list for easy access
    idx_map = {rec.id: i for i, rec in enumerate(assembly_records)}
    
    # 2. Load your Exon sequences (Update filenames as needed)
    print("Loading exon patches...")
    # NOTE: Replace these with your actual filenames
    kl3_ex_14_15 = load_exon_sequence("kl3_exons_14and15.fasta")
    kl3_ex_2_3   = load_exon_sequence("kl3_exons_2and3.fasta")
    kl5_ex_16_17 = load_exon_sequence("kl5_exons16and17.fasta")
    kl5_ex_7_15  = load_exon_sequence("kl5_exons7_through_15.fasta")
    
    gap = "N" * 100
    
    # --- SURGERY ON Y_ptg000011l (kl-3) ---
    if "Y_ptg000011l" in idx_map:
        print("Patching kl-3 on Y_ptg000011l...")
        # Get the sequence as a string
        target_idx = idx_map["Y_ptg000011l"]
        scaff_str = str(assembly_records[target_idx].seq)
        
        # Patch highest coordinate first: 6,251,425
        pos1 = 6251425
        scaff_str = scaff_str[:pos1] + gap + kl3_ex_14_15 + gap + scaff_str[pos1:]
        
        # Patch lower coordinate: 5,743,472
        pos2 = 5743472
        scaff_str = scaff_str[:pos2] + gap + kl3_ex_2_3 + gap + scaff_str[pos2:]
        
        # Update the record with a proper Seq object
        assembly_records[target_idx].seq = Seq(scaff_str)
    else:
        print("Error: Y_ptg000011l not found in assembly!")

    # --- SURGERY ON Y_scaffold6 (kl-5) ---
    if "Y_scaffold6" in idx_map:
        print("Patching kl-5 on Y_scaffold6...")
        target_idx = idx_map["Y_scaffold6"]
        scaff_str = str(assembly_records[target_idx].seq)
        
        # Patch highest coordinate first: 260,568
        pos3 = 260568
        scaff_str = scaff_str[:pos3] + gap + kl5_ex_16_17 + gap + scaff_str[pos3:]
        
        # Patch lower coordinate: 256,712
        pos4 = 256712
        scaff_str = scaff_str[:pos4] + gap + kl5_ex_7_15 + gap + scaff_str[pos4:]
        
        assembly_records[target_idx].seq = Seq(scaff_str)
    else:
        print("Error: Y_scaffold6 not found in assembly!")

    # 3. Write the finalized Curated Assembly
    print("Writing final assembly...")
    SeqIO.write(assembly_records, "Delgado_ISO1_Curated_FINAL.fasta", "fasta")
    print("Success! Final assembly created: Delgado_ISO1_Curated_FINAL.fasta")

if __name__ == "__main__":
    patch_assembly()
