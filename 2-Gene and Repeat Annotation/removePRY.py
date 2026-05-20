import sys

# ==========================================
# CONFIGURATION
# ==========================================
# Make sure these filenames match exactly what you have in your folder
FASTA_IN = "DelgadoYISO-1v10.fasta"       # Your Patched WDY Assembly
GFF_IN   = "DelgadoY_Corrected_Final.gff3" # Your validated annotation

FASTA_OUT = "DelgadoY_ISO1_v11_Final.fasta"
GFF_OUT   = "DelgadoY_ISO1_v11_Final.gff3"

TARGET_ID = "PRY_Copy2"  # The ID of the artifact to remove
SCAFFOLD  = "Y_ptg000011l" # The scaffold it lives on

def run_surgery():
    print(f"--- PREPARING SURGERY ON {TARGET_ID} ---")
    
    # 1. LOCATE THE ARTIFACT
    # We need the exact Start and End to know what to cut.
    cut_start = None
    cut_end = None
    
    with open(GFF_IN, 'r') as f:
        for line in f:
            if TARGET_ID in line and "gene" in line.split('\t')[2]:
                cols = line.strip().split('\t')
                cut_start = int(cols[3])
                cut_end = int(cols[4])
                break
    
    if not cut_start:
        print(f"Error: Could not find gene coordinates for {TARGET_ID} in {GFF_IN}")
        return

    deletion_size = cut_end - cut_start + 1
    print(f"Target Found: {TARGET_ID}")
    print(f"Location: {SCAFFOLD}:{cut_start}-{cut_end}")
    print(f"Deletion Size: {deletion_size} bp")

    # 2. MODIFY THE FASTA (PHYSICAL DELETION)
    print(f"\nReading {FASTA_IN}...")
    
    new_fasta_lines = []
    current_header = None
    seq_fragments = []
    
    # Helper to process a scaffold once fully read
    def process_sequence(header, sequence):
        scaff_name = header.split()[0]
        if scaff_name == SCAFFOLD:
            print(f"  Processing {scaff_name} (Original Len: {len(sequence)} bp)...")
            
            # Python is 0-based, GFF is 1-based.
            # We want to keep 0 to (start-1-1) AND (end) to end
            # indices: 
            #   Keep 0 to (cut_start - 2)  [Prefix]
            #   Skip (cut_start - 1) to (cut_end - 1)
            #   Keep (cut_end) to end      [Suffix]
            
            # Example: Delete bases 3-4 (Length 2). 1-based.
            # Seq: A B C D E
            # 0:A, 1:B, 2:C, 3:D, 4:E
            # Cut start=3, end=4.
            # Prefix: Seq[:2] -> A B (indices 0,1)
            # Suffix: Seq[4:] -> E   (index 4+)
            
            p_idx = cut_start - 1
            s_idx = cut_end
            
            prefix = sequence[:p_idx]
            suffix = sequence[s_idx:]
            
            new_seq = prefix + suffix
            print(f"  Cut complete. New Len: {len(new_seq)} bp.")
            return new_seq
        else:
            return sequence

    # Read FASTA
    with open(FASTA_IN, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_header:
                    full_seq = "".join(seq_fragments)
                    mod_seq = process_sequence(current_header, full_seq)
                    new_fasta_lines.append(f">{current_header}\n")
                    # Write in 60-char blocks
                    for i in range(0, len(mod_seq), 60):
                        new_fasta_lines.append(mod_seq[i:i+60] + "\n")
                
                current_header = line[1:]
                seq_fragments = []
            else:
                seq_fragments.append(line)
        # Process last entry
        if current_header:
            full_seq = "".join(seq_fragments)
            mod_seq = process_sequence(current_header, full_seq)
            new_fasta_lines.append(f">{current_header}\n")
            for i in range(0, len(mod_seq), 60):
                new_fasta_lines.append(mod_seq[i:i+60] + "\n")

    print(f"Writing {FASTA_OUT}...")
    with open(FASTA_OUT, 'w') as f:
        f.writelines(new_fasta_lines)

    # 3. MODIFY THE GFF (COORDINATE SHIFT)
    print(f"\nUpdating Annotation Coordinates -> {GFF_OUT}...")
    
    with open(GFF_IN, 'r') as fin, open(GFF_OUT, 'w') as fout:
        for line in fin:
            if line.startswith("#"):
                fout.write(line)
                continue
            
            cols = line.strip().split('\t')
            if len(cols) < 9: continue
            
            gff_scaff = cols[0]
            gff_start = int(cols[3])
            gff_end   = int(cols[4])
            attributes = cols[8]
            
            # A. SKIP the deleted artifact
            if TARGET_ID in attributes:
                continue # Delete this line
                
            # B. Check if we need to shift
            if gff_scaff == SCAFFOLD:
                # If feature is AFTER the cut, shift it left
                if gff_start > cut_end:
                    new_start = gff_start - deletion_size
                    new_end   = gff_end - deletion_size
                    
                    # Reconstruct line
                    cols[3] = str(new_start)
                    cols[4] = str(new_end)
                    fout.write("\t".join(cols) + "\n")
                    
                # If feature is BEFORE the cut, keep it as is
                elif gff_end < cut_start:
                    fout.write(line)
                    
                # If feature OVERLAPS the cut (shouldn't happen if we cut a gene), warn
                else:
                    print(f"Warning: Feature overlaps deletion zone! Dropping: {attributes}")
            else:
                # Other scaffolds untouched
                fout.write(line)

    print("\n--- SURGERY SUCCESSFUL ---")
    print(f"1. Removed {deletion_size} bp from {SCAFFOLD}.")
    print(f"2. Deleted {TARGET_ID} from GFF.")
    print(f"3. Shifted downstream genes to match new coordinates.")
    print(f"Final files: {FASTA_OUT}, {GFF_OUT}")

if __name__ == "__main__":
    run_surgery()
