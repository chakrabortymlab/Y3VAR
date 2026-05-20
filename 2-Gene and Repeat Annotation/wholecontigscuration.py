import sys
import os
from Bio import SeqIO
from Bio.Seq import Seq

def run_curation():
    # --- Configuration ---
    input_fasta = "DelgadoYISO-1v5.fasta"  
    chang_kl5_fasta = "Reference_Chang_kl-5.fasta"
    # chang_ppr_fasta removed - no longer needed
    ppr_exon3_fasta = "Ppr-yexon3.fasta"
    output_fasta = "DelgadoYISO-1v8.fasta"
    
    print(f"--- Master Curation Pipeline (v9: Exon3 Replaces Chang Frag) ---")
    
    # 1. Load Assemblies
    print(f"Loading {input_fasta}...")
    genome_dict = SeqIO.to_dict(SeqIO.parse(input_fasta, "fasta"))
    
    print(f"Loading {chang_kl5_fasta}...")
    chang_kl5_dict = SeqIO.to_dict(SeqIO.parse(chang_kl5_fasta, "fasta"))
    
    print(f"Loading {ppr_exon3_fasta}...")
    ppr_exon3_record = SeqIO.read(ppr_exon3_fasta, "fasta")

    # 2. Extract All Fragments
    sequences = {}

    # -- A. New Exon 3 --
    sequences["Ppr_Exon3"] = ppr_exon3_record.seq
    print(f"  Captured Ppr_Exon3 ({len(sequences['Ppr_Exon3'])} bp)")

    # -- B. Extract Y_unlocalizedcontigs --
    target_unlocalized = [
        "Y_unlocalizedcontig1", 
        "Y_unlocalizedcontig3", 
        "Y_unlocalizedcontig4", 
        "Y_unlocalizedcontig6"
    ]
    
    for ctg in target_unlocalized:
        if ctg in genome_dict:
            sequences[ctg] = genome_dict[ctg].seq
            print(f"  Captured {ctg} ({len(sequences[ctg])} bp)")
        else:
            print(f"  ERROR: {ctg} missing from input genome!")
            return

    # -- C. Extract Chang Fragments (kl-5 only) --
    def get_slice(record, start, end):
        return record.seq[start-1:end]

    kl5_contig = "Y_scaffold6:31929-426170"
    if kl5_contig in chang_kl5_dict:
        sequences["Chang_kl5_Frag1"] = get_slice(chang_kl5_dict[kl5_contig], 224669, 228645)
        sequences["Chang_kl5_Frag2"] = get_slice(chang_kl5_dict[kl5_contig], 362856, 368028)
        sequences["Chang_kl5_Frag3"] = get_slice(chang_kl5_dict[kl5_contig], 373760, 374242)
        print(f"  Captured Chang kl-5 fragments from {kl5_contig}")
    else:
        print(f"  ERROR: Contig {kl5_contig} not found in {chang_kl5_fasta}")
        return

    # 3. Define Operations
    # Format: { "scaf": name, "pos": 1-based-int, "id": key, "prio": int }
    # Logic: Higher Priority = Processed First = Pushed Downstream by Lower Priority inserts at same site.
    ops = []

    # --- Y_scaffold5 ---
    # 1. Y_unlocalizedcontig6 at 672,257
    ops.append({"scaf": "Y_scaffold5", "pos": 672257, "id": "Y_unlocalizedcontig6", "prio": 1})

    # --- Y_ptg000011l ---
    # 1. Y_unlocalizedcontig4 at 2,371,175
    # 2. Ppr_Exon3 RIGHT AFTER Y_unlocalizedcontig4.
    # Desired Physical Order: [Contig4] -> [Exon3]
    # Execution Order: Exon3 (Prio 2), then Contig4 (Prio 1)
    ops.append({"scaf": "Y_ptg000011l", "pos": 2371175, "id": "Y_unlocalizedcontig4", "prio": 1})         
    ops.append({"scaf": "Y_ptg000011l", "pos": 2371175, "id": "Ppr_Exon3", "prio": 2}) 

    # --- Y_scaffold6 ---
    # 1. Y_unlocalizedcontig3 at 266,192
    # 2. Chang_kl5_Frag1 RIGHT BEFORE Y_unlocalizedcontig3.
    # Desired Physical Order: [Chang] -> [Contig3]
    # Execution Order: Contig3 (Prio 2), then Chang (Prio 1)
    ops.append({"scaf": "Y_scaffold6", "pos": 266192, "id": "Y_unlocalizedcontig3", "prio": 2})          
    ops.append({"scaf": "Y_scaffold6", "pos": 266192, "id": "Chang_kl5_Frag1", "prio": 1}) 

    # 3. Y_unlocalizedcontig1 at 392,716
    # 4. Chang_kl5_Frag2 RIGHT BEFORE Y_unlocalizedcontig1.
    # Desired Physical Order: [Chang] -> [Contig1]
    ops.append({"scaf": "Y_scaffold6", "pos": 392716, "id": "Y_unlocalizedcontig1", "prio": 2})
    ops.append({"scaf": "Y_scaffold6", "pos": 392716, "id": "Chang_kl5_Frag2", "prio": 1})

    # 5. Chang_kl5_Frag3 at 401,613
    ops.append({"scaf": "Y_scaffold6", "pos": 401613, "id": "Chang_kl5_Frag3", "prio": 1})


    # 4. Execute Operations
    print("\n--- Executing Insertions ---")
    
    scaf_ops = {}
    for op in ops:
        if op["scaf"] not in scaf_ops: scaf_ops[op["scaf"]] = []
        scaf_ops[op["scaf"]].append(op)

    for scaf, tasks in scaf_ops.items():
        if scaf not in genome_dict:
            print(f"Skipping {scaf} (not found)...")
            continue
            
        print(f"Processing {scaf}...")
        
        # Sort: Pos Descending, then Prio Descending
        # Ex: Prio 2 (Exon3) -> Prio 1 (Contig4)
        tasks.sort(key=lambda x: (x["pos"], x["prio"]), reverse=True)
        
        mutable_seq = genome_dict[scaf].seq
        
        for t in tasks:
            pos = t["pos"]
            seq_id = t["id"]
            insert_seq = sequences[seq_id]
            
            upstream = mutable_seq[:pos]
            downstream = mutable_seq[pos:]
            mutable_seq = upstream + insert_seq + downstream
            
            print(f"  -> Inserted {seq_id} at orig pos {pos} (Prio {t['prio']})")
            
        genome_dict[scaf].seq = mutable_seq
        print(f"  -> New length: {len(mutable_seq)}")

    # 5. Remove Redundant Contigs
    print("\n--- Cleaning Up ---")
    to_delete = [
        "Y_unlocalizedcontig1", 
        "Y_unlocalizedcontig2", 
        "Y_unlocalizedcontig3", 
        "Y_unlocalizedcontig4", 
        "Y_unlocalizedcontig5", 
        "Y_unlocalizedcontig6"
    ]
    
    final_records = []
    for rec_id, record in genome_dict.items():
        if rec_id in to_delete:
            print(f"  Removing {rec_id}")
            continue
        final_records.append(record)

    # 6. Write Output
    print(f"\nWriting {output_fasta}...")
    with open(output_fasta, "w") as out_f:
        SeqIO.write(final_records, out_f, "fasta")
    print("Done.")

if __name__ == "__main__":
    run_curation()
