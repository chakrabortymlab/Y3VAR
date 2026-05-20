import sys

def patch_genome(fasta_file, scaffold_name):
    print(f"Reading {fasta_file}...")
    
    # 1. Read the Genome into Memory
    # (Simple parser for single-line or wrapped FASTA)
    scaffolds = {}
    current_scaff = None
    seq_fragments = []
    
    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_scaff:
                    scaffolds[current_scaff] = "".join(seq_fragments)
                current_scaff = line[1:].split()[0] # Get ID
                seq_fragments = []
            else:
                seq_fragments.append(line)
        if current_scaff:
            scaffolds[current_scaff] = "".join(seq_fragments)

    if scaffold_name not in scaffolds:
        print(f"Error: Scaffold {scaffold_name} not found!")
        return

    seq = scaffolds[scaffold_name]
    print(f"Original Length: {len(seq)} bp")

    # ==========================================
    # 2. DEFINE COORDINATES (1-based input -> 0-based python)
    # ==========================================
    # Orphan Exon (Exon 1): ~2,480,175 to 2,480,551
    # We add some padding (50bp) to be safe
    orphan_start = 2480175 - 1
    orphan_end   = 2480551
    
    # Target Location: Just upstream of Body (Body starts ~279,109 on Minus strand)
    # We insert at 280,000 to leave a small ~1kb intron
    insert_point = 280000 

    # ==========================================
    # 3. PERFORM SURGERY
    # ==========================================
    
    # A. Extract the Orphan Sequence
    orphan_seq = seq[orphan_start:orphan_end]
    print(f"Extracted Orphan Sequence ({len(orphan_seq)} bp)")
    
    # B. Mask the old spot with 'N's (preserve length of that region)
    # We do this BEFORE insertion to keep indices simple for now, 
    # but we are creating a NEW sequence, so exact indices matter.
    
    prefix = seq[:insert_point]
    suffix = seq[insert_point:]
    
    # C. Reconstruct Sequence: 
    # Prefix + [ORPHAN] + [Spacer] + Suffix (with Orphan masked out)
    
    # Wait, if we insert, the suffix indices shift. 
    # We need to handle the masking in the suffix carefully.
    
    # Let's mask the original sequence first
    seq_list = list(seq)
    for i in range(orphan_start, orphan_end):
        seq_list[i] = 'N'
    masked_seq = "".join(seq_list)
    
    # Now cut and split for insertion
    new_prefix = masked_seq[:insert_point]
    new_suffix = masked_seq[insert_point:]
    
    # Insert the Orphan (Alive) at the new spot
    # We add 100 'N's as a fake intron between Orphan and Body
    fake_intron = "N" * 100 
    
    final_seq = new_prefix + orphan_seq + fake_intron + new_suffix
    
    scaffolds[scaffold_name] = final_seq
    
    print(f"Patched Length: {len(final_seq)} bp")
    print(f"  (Added {len(orphan_seq) + 100} bp to the scaffold)")

    # ==========================================
    # 4. SAVE OUTPUT
    # ==========================================
    out_name = "DelgadoY_Patched.fasta"
    with open(out_name, 'w') as f:
        for name, s in scaffolds.items():
            f.write(f">{name}\n")
            # Write in 60-char lines
            for i in range(0, len(s), 60):
                f.write(s[i:i+60] + "\n")
                
    print(f"Done! Saved to {out_name}")

# RUN IT
# Replace with your actual FASTA filename
patch_genome("DelgadoYISO-1v9.fasta", "Y_ptg000011l")
