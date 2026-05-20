import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def run_insertion(target_fasta, source_fasta, scaffold_name, insert_pos, output_file):
    print(f"--- Manual Insertion Tool ---")
    print(f"Target Genome: {target_fasta}")
    print(f"Source Seq:    {source_fasta}")
    print(f"Insertion loc: {scaffold_name} at base {insert_pos}")

    # 1. Load the Source Sequence (Chang PRY/FDY)
    # We take the first record found in the file, ignoring the header.
    source_record = next(SeqIO.parse(source_fasta, "fasta"))
    insert_seq_str = str(source_record.seq)
    print(f"Loaded source sequence length: {len(insert_seq_str)} bp")

    # 2. Load the Genome Assembly
    # We use to_dict to make it mutable
    genome_dict = SeqIO.to_dict(SeqIO.parse(target_fasta, "fasta"))

    if scaffold_name not in genome_dict:
        print(f"ERROR: Scaffold '{scaffold_name}' not found in genome!")
        return

    target_record = genome_dict[scaffold_name]
    original_seq = str(target_record.seq)
    original_len = len(original_seq)

    # 3. Perform Insertion
    # Convert 1-based coordinate to 0-based index
    # Example: Insert at 5 means mapped between base 5 and 6. 
    # Python index 5 starts at the 6th base. 
    split_index = insert_pos 
    
    upstream = original_seq[:split_index]
    downstream = original_seq[split_index:]
    
    # Construct new sequence
    new_seq_str = upstream + "NNNNNNNNNN" + insert_seq_str + "NNNNNNNNNN" + downstream
    
    # Update the record
    genome_dict[scaffold_name].seq = Seq(new_seq_str)
    
    print(f"Insertion successful.")
    print(f"Old Scaffold Length: {original_len}")
    print(f"New Scaffold Length: {len(genome_dict[scaffold_name].seq)}")
    print(f"Delta: +{len(genome_dict[scaffold_name].seq) - original_len} bp")

    # 4. Write Output
    print(f"Writing to {output_file}...")
    with open(output_file, "w") as out_f:
        SeqIO.write(genome_dict.values(), out_f, "fasta")
    print("Done.")

if __name__ == "__main__":
    # You can run this directly or via command line
    # Default values based on your prompt:
    
    TARGET = "DelgadoYISO-1-v4.fasta"
    SOURCE = "Reference_Chang_PRY.fasta" # Containing the PRY sequence
    SCAFFOLD = "Y_ptg000011l"
    POS = 6466047 # 1-based coordinate
    OUTPUT = "DelgadoYISO-1v5.fasta"

    run_insertion(TARGET, SOURCE, SCAFFOLD, POS, OUTPUT)
