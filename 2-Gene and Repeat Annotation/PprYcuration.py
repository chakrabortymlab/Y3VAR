import sys
import os
import re
from Bio import SeqIO
from Bio.Seq import Seq

# ==========================================
# CONFIGURATION
# ==========================================
ASSEMBLY_FILE = "DelgadoYISO-1v8.fasta"
EXON_FILE = "PprYexons.fasta"
SAM_FILE = "ppry_normalization.sam"
OUTPUT_FILE = "DelgadoYISO-v10_Normalized.fasta"

# Alignment Thresholds
MIN_IDENTITY = 0.80      # Detect anything resembling an exon
MAX_GENE_SPAN = 2000000  # 2 Mb max distance for "Same Gene" logic

def parse_sam(sam_file):
    """Parses SAM to find all potential hits for each exon."""
    hits = []
    print(f"Reading {sam_file}...")
    
    with open(sam_file, 'r') as f:
        for line in f:
            if line.startswith("@"): continue
            cols = line.split("\t")
            if len(cols) < 11: continue
            
            flag = int(cols[1])
            if flag & 4: continue # Unmapped
            
            qname = cols[0]
            chrom = cols[2]
            start = int(cols[3])
            cigar = cols[5]
            seq_len = len(cols[9])
            
            # Calculate Identity using NM tag (Edit Distance)
            nm_tag = next((x for x in cols[11:] if x.startswith("NM:i:")), None)
            nm = int(nm_tag.split(":")[2]) if nm_tag else 0
            matches = seq_len - nm
            identity = matches / seq_len
            
            # Calculate Reference End Pos (approximate based on CIGAR)
            ref_len = 0
            for l_str, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar):
                if op in ['M', '=', 'X', 'D', 'N']:
                    ref_len += int(l_str)
            end = start + ref_len - 1
            
            # Parse Exon Number from header (e.g., ...:5)
            ex_match = re.search(r":(\d+)", qname)
            exon_num = int(ex_match.group(1)) if ex_match else 0
            
            if identity >= MIN_IDENTITY and exon_num > 0:
                hits.append({
                    "chrom": chrom,
                    "start": start,
                    "end": end,
                    "exon": exon_num,
                    "id": identity,
                    "qname": qname
                })
    return hits

def determine_main_locus(hits):
    """
    Finds the cluster containing the most unique exons.
    Returns: (scaffold_name, start_range, end_range)
    """
    if not hits: return None, 0, 0
    
    # Sort by pos
    hits.sort(key=lambda x: (x["chrom"], x["start"]))
    
    clusters = []
    current_cluster = [hits[0]]
    
    for i in range(1, len(hits)):
        prev = current_cluster[-1]
        curr = hits[i]
        
        # Check distance
        dist = curr["start"] - prev["end"]
        if prev["chrom"] == curr["chrom"] and dist < MAX_GENE_SPAN:
            current_cluster.append(curr)
        else:
            clusters.append(current_cluster)
            current_cluster = [curr]
    clusters.append(current_cluster)
    
    # Select best cluster (most unique exons found)
    best_cluster = []
    best_score = 0
    
    for cl in clusters:
        unique_exons = set(h["exon"] for h in cl)
        score = len(unique_exons)
        if score > best_score:
            best_score = score
            best_cluster = cl
            
    if not best_cluster: return None, 0, 0
    
    chrom = best_cluster[0]["chrom"]
    l_start = min(h["start"] for h in best_cluster)
    l_end = max(h["end"] for h in best_cluster)
    
    return chrom, l_start, l_end

def run_normalization():
    print("--- Starting Assembly Normalization (v10) ---")
    
    # 0. Load Exon Sequences
    exon_seqs = {}
    for rec in SeqIO.parse(EXON_FILE, "fasta"):
        match = re.search(r":(\d+)", rec.id)
        if match:
            exon_seqs[int(match.group(1))] = str(rec.seq)
            
    # 1. Parse SAM
    if not os.path.exists(SAM_FILE):
        print(f"Error: {SAM_FILE} not found. Please run minimap2 first.")
        return

    hits = parse_sam(SAM_FILE)
    if not hits:
        print("No valid exon matches found.")
        return
    
    # 2. Identify Main Locus
    locus_chrom, locus_start, locus_end = determine_main_locus(hits)
    print(f"Identified Main Gene Locus on {locus_chrom} ({locus_start}-{locus_end})")
    
    # 3. Plan Edits
    # Format: {"chrom": str, "start": int, "end": int, "type": "DEL" or "REPLACE" or "INS", "seq": str, "prio": int}
    edits = []
    fixed_exons = set()
    
    # Sort hits for processing
    hits.sort(key=lambda x: (x["chrom"], x["start"]))
    
    # Group by Exon Number
    hits_by_exon = {i: [] for i in range(1, 7)}
    for h in hits:
        hits_by_exon[h["exon"]].append(h)
        
    for exon_num in range(1, 7):
        exon_hits = hits_by_exon[exon_num]
        
        # Define Inside vs Outside Locus
        inside = [h for h in exon_hits if h["chrom"] == locus_chrom and h["start"] >= locus_start and h["start"] <= locus_end + 500000]
        outside = [h for h in exon_hits if h not in inside]
        
        # A. DELETE OUTSIDERS
        for h in outside:
            print(f"  [Plan] Remove Duplicate Exon {exon_num} on {h['chrom']}:{h['start']}")
            edits.append({
                "chrom": h["chrom"], "start": h["start"], "end": h["end"], 
                "type": "DEL", "prio": 1
            })
            
        # B. HANDLE INSIDERS
        if inside:
            # Pick best match to serve as the anchor for replacement
            best_hit = max(inside, key=lambda x: x["id"])
            
            # Force perfect sequence (REPLACE)
            print(f"  [Plan] Enforce Perfect Exon {exon_num} at {locus_chrom}:{best_hit['start']}")
            edits.append({
                "chrom": locus_chrom, "start": best_hit["start"], "end": best_hit["end"],
                "type": "REPLACE", "seq": exon_seqs[exon_num], "prio": 2
            })
            fixed_exons.add(exon_num)
            
            # Delete local duplicates
            for h in inside:
                if h != best_hit:
                    print(f"  [Plan] Remove Local Duplicate Exon {exon_num} at {h['start']}")
                    edits.append({
                        "chrom": locus_chrom, "start": h["start"], "end": h["end"], 
                        "type": "DEL", "prio": 1
                    })
        else:
            print(f"  [WARN] Exon {exon_num} MISSING from main locus.")

    # C. Handle Missing Exons (Insertions)
    missing = [e for e in range(1, 7) if e not in fixed_exons]
    
    # Map fixed exons to their locations to use as anchors
    anchor_map = {} 
    for op in edits:
        if op["type"] == "REPLACE":
            for e_num, seq in exon_seqs.items():
                if seq == op["seq"]:
                    anchor_map[e_num] = op["start"]
                    break
    
    for m in missing:
        insert_pos = -1
        
        # Look upstream
        prev_ex = m - 1
        while prev_ex > 0:
            if prev_ex in anchor_map:
                insert_pos = anchor_map[prev_ex] + 2000 # 2kb downstream
                break
            prev_ex -= 1
            
        if insert_pos == -1:
            # Look downstream
            next_ex = m + 1
            while next_ex <= 6:
                if next_ex in anchor_map:
                    insert_pos = max(1, anchor_map[next_ex] - 2000) # 2kb upstream
                    break
                next_ex += 1
                
        if insert_pos != -1:
            print(f"  [Plan] Insert Missing Exon {m} at {locus_chrom}:{insert_pos}")
            edits.append({
                "chrom": locus_chrom, "start": insert_pos, "end": insert_pos, 
                "type": "INS", "seq": exon_seqs[m], "prio": 3
            })
        else:
            print(f"  [CRITICAL] Could not place Exon {m} (No neighbors found).")

    # 4. Execute Edits
    # Sort: Chromosome -> Pos Descending -> Priority
    # Processing backwards (Descending Pos) ensures coordinates don't shift for remaining edits.
    edits.sort(key=lambda x: (x["chrom"], x["start"], x["prio"]), reverse=True)
    
    print("\n--- Applying Edits ---")
    print(f"Loading {ASSEMBLY_FILE}...")
    genome = SeqIO.to_dict(SeqIO.parse(ASSEMBLY_FILE, "fasta"))
    
    for op in edits:
        scaf = op["chrom"]
        if scaf not in genome: continue
        
        mutable_seq = genome[scaf].seq
        start = op["start"] - 1 # 0-based
        end = op["end"]         # 0-based exclusive
        
        if op["type"] == "DEL":
            genome[scaf].seq = mutable_seq[:start] + mutable_seq[end:]
            
        elif op["type"] == "REPLACE":
            new_seq = Seq(op["seq"])
            genome[scaf].seq = mutable_seq[:start] + new_seq + mutable_seq[end:]
            
        elif op["type"] == "INS":
            new_seq = Seq(op["seq"])
            genome[scaf].seq = mutable_seq[:start] + new_seq + mutable_seq[start:]
            
    print(f"Writing {OUTPUT_FILE}...")
    with open(OUTPUT_FILE, "w") as out_f:
        SeqIO.write(genome.values(), out_f, "fasta")
    print("Done.")

if __name__ == "__main__":
    if not os.path.exists(ASSEMBLY_FILE) or not os.path.exists(EXON_FILE):
        print("Missing FASTA files.")
    else:
        run_normalization()
