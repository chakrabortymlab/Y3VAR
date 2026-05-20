import argparse
import os
import sys
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import BiopythonDeprecationWarning
import warnings

#python CDSreconstruction.py --genome DelgadoYA3.fasta --psl results_A3_CDS.psl --refs CarvalhoY_CDS.fasta --strain A3 --outdir A3
# =============================================================================
# Y-CHROMOSOME PIPELINE v11 (DNA CDS MODE)
# Reconstructs genes using a DNA CDS Reference query instead of Protein.
# =============================================================================

# Silence Biopython warnings
warnings.simplefilter('ignore', BiopythonDeprecationWarning)

# GLOBAL LOGS
AUDIT_LOG = []
SUMMARY_LOG = []

def log_decision(strain, gene, copy, step, action, reason, details, 
                 scaffold="N/A", dna_coords="N/A", strand="N/A"):
    """Helper to append structured data to the audit log."""
    AUDIT_LOG.append({
        "Strain": strain,
        "Gene": gene,
        "Copy_ID": copy,
        "Step": step,
        "Action": action,
        "Logic/Math": reason,
        "Details": details,
        "Scaffold": scaffold,
        "DNA_Coordinates": dna_coords,
        "Strand": strand
    })

def parse_psl_and_extract(psl_path, genome_dict, gene_name_filter=None):
    """
    Step 1: Parse PSL (DNA Query -> DNA Target).
    Extracts the matching target sequence for every block.
    """
    gene_hits = {}
    print(f"  - Parsing PSL file: {psl_path}...")
    
    with open(psl_path, 'r') as f:
        for line in f:
            cols = line.strip().split('\t')
            # PSL format usually has 21 columns. Valid lines start with numbers.
            if len(cols) < 21 or not cols[0].isdigit(): continue
            
            # Col 9 = qName (Query Gene Name, e.g., "Pp1Y1")
            gene_name = cols[9]
            if gene_name_filter and gene_name != gene_name_filter: continue

            matches = int(cols[0]) # Nucleotide matches
            t_name = cols[13]      # Target Scaffold Name
            
            if t_name not in genome_dict: continue

            # DNA Extraction
            scaffold_seq = genome_dict[t_name].seq
            strand = cols[8] # e.g., "++" or "+-"
            
            # qStart/qEnd are now Nucleotide coordinates (0-based)
            q_start = int(cols[11])
            q_end = int(cols[12])
            
            block_count = int(cols[17])
            # Remove trailing commas and convert to ints
            block_sizes = [int(x) for x in cols[18].strip(',').split(',')]
            t_starts = [int(x) for x in cols[20].strip(',').split(',')]
            t_size = int(cols[14])
            
            full_hit_seq = ""
            dna_span_starts = []
            dna_span_ends = []

            for i in range(block_count):
                size_nt = block_sizes[i] 
                # Note: No "* 3" multiplication needed here anymore!
                
                start_on_scaffold = t_starts[i]
                end_on_scaffold = start_on_scaffold + size_nt
                
                dna_span_starts.append(start_on_scaffold)
                dna_span_ends.append(end_on_scaffold)
                
                # Strand handling:
                # If target strand is '-' (e.g. "+-" or "--"), BLAT coordinates 
                # are usually relative to the forward strand, but the biological 
                # sequence is reverse complement.
                if strand.endswith('-'):
                    # Biopython and BLAT coordinate math for reverse strand:
                    # BLAT tStart is usually from the start of the REVERSE strand if the target is reversed?
                    # Actually, standard PSL is always + strand coordinates for Target.
                    # We just extract + strand and then RevComp it.
                    
                    chunk = str(scaffold_seq[start_on_scaffold:end_on_scaffold])
                    full_hit_seq += str(Seq(chunk).reverse_complement())
                else:
                    full_hit_seq += str(scaffold_seq[start_on_scaffold:end_on_scaffold])

            # Store the hit
            hit = {
                'gene': gene_name,
                'qStart': q_start,   # Nucleotide Index in Reference CDS
                'qEnd': q_end,       # Nucleotide Index in Reference CDS
                'dna_seq': full_hit_seq,
                'scaffold': t_name,
                't_start_min': min(dna_span_starts),
                't_end_max': max(dna_span_ends),
                'dna_range_str': f"{min(dna_span_starts)}-{max(dna_span_ends)}",
                'cds_range_str': f"nt {q_start}-{q_end}",
                'strand': strand,
                'matches': matches,
                'status': 'active' 
            }
            gene_hits.setdefault(gene_name, []).append(hit)
            
    return gene_hits

def resolve_conflicts(all_hits, strain_id):
    """
    Step 2: Resolve Genomic Overlaps.
    If two hits claim the exact same genomic locus, keep the one with better score.
    """
    print("\n  - Resolving Genomic Conflicts...")
    SCORE_MARGIN = 1.1 # 10% margin
    
    flat_hits = []
    for gene, hits in all_hits.items():
        flat_hits.extend(hits)
    
    # Sort by genomic position
    flat_hits.sort(key=lambda x: (x['scaffold'], x['t_start_min']))
    
    scaffold_groups = {}
    for hit in flat_hits:
        scaffold_groups.setdefault(hit['scaffold'], []).append(hit)
        
    for scaf, hits in scaffold_groups.items():
        if not hits: continue
        
        for i in range(len(hits)):
            if hits[i]['status'] == 'rejected': continue
            
            for j in range(i + 1, len(hits)):
                if hits[j]['status'] == 'rejected': continue
                
                # Check for physical overlap on the scaffold
                if hits[j]['t_start_min'] < hits[i]['t_end_max']:
                    # Overlap detected
                    h1 = hits[i]
                    h2 = hits[j]
                    
                    if h1['gene'] == h2['gene']: continue # Self-overlaps handled in chaining
                        
                    # Score Comparison
                    if h1['matches'] >= h2['matches']:
                        winner, loser = h1, h2
                    else:
                        winner, loser = h2, h1
                    
                    if winner['matches'] > (loser['matches'] * SCORE_MARGIN):
                        loser['status'] = 'rejected'
                        log_decision(strain_id, loser['gene'], "Rejected", "Conflict Resolution", 
                                     "Genomic Collision", "Lower Match Score", 
                                     f"Overlapped by {winner['gene']}", 
                                     scaffold=loser['scaffold'], dna_coords=loser['dna_range_str'])
                else:
                    break
                    
    cleaned_hits = {}
    for hit in flat_hits:
        if hit['status'] == 'active':
            cleaned_hits.setdefault(hit['gene'], []).append(hit)
            
    return cleaned_hits

def find_single_best_copy(hits, ref_cds_len, strain_id, gene_name):
    """
    Step 3: Greedy Linear Stitching (DNA Mode)
    """
    
    # Sort by Query (CDS) Start Position
    hits.sort(key=lambda x: (x['qStart'], -x['matches']))
    
    chain = []
    
    for hit in hits:
        if not chain:
            chain.append(hit)
            log_decision(strain_id, gene_name, "Best_Copy", "Chaining", "Anchor", 
                         "First Fragment", "Started Chain", 
                         scaffold=hit['scaffold'], dna_coords=hit['cds_range_str'])
            continue
            
        last_hit = chain[-1]
        
        # LOGIC:
        # Check if this hit continues the CDS sequence.
        # Allow small overlap (e.g. 30bp) to account for alignment edge fuzziness
        
        overlap_tolerance = 30 # bases
        
        if hit['qStart'] >= (last_hit['qEnd'] - overlap_tolerance): 
            chain.append(hit)
            
            # Check for Scaffold Jump
            if hit['scaffold'] != last_hit['scaffold']:
                log_decision(strain_id, gene_name, "Best_Copy", "Chaining", "Scaffold Jump", 
                             "Fragmentation Detected", f"Jumped from {last_hit['scaffold']} to {hit['scaffold']}", 
                             scaffold=hit['scaffold'], dna_coords=hit['cds_range_str'])
            else:
                log_decision(strain_id, gene_name, "Best_Copy", "Chaining", "Extension", 
                             "Linear Extension", "Added to chain", 
                             scaffold=hit['scaffold'], dna_coords=hit['cds_range_str'])
        else:
            # Overlap rejection
            log_decision(strain_id, gene_name, "Best_Copy", "Chaining", "Skipped Overlap", 
                         "Redundant / Duplicate", f"Hit {hit['cds_range_str']} overlaps existing coverage", 
                         scaffold=hit['scaffold'])

    # --- STITCHING ---
    copy_id = f"{strain_id}_{gene_name}_CDS"
    
    stitched_seq = Seq("")
    current_cds_end = chain[0]['qStart'] 
    gap_n_count = 0
    
    for hit in chain:
        raw_seq = Seq(hit['dna_seq'])
        gap_nt = hit['qStart'] - current_cds_end
        
        if len(stitched_seq) == 0:
            stitched_seq += raw_seq
            current_cds_end = hit['qEnd']
        else:
            if gap_nt > 0:
                # Fill Gap with Ns
                gap_n_count += gap_nt
                stitched_seq += Seq("N" * gap_nt)
                stitched_seq += raw_seq
                current_cds_end = hit['qEnd']
                log_decision(strain_id, gene_name, "Best_Copy", "Stitching", "Gap Filled", 
                             f"{gap_nt} bp Gap", "Inserted Ns")
                
            elif gap_nt < 0:
                # Trim overlap
                trim_bp = abs(gap_nt)
                if trim_bp >= len(raw_seq): 
                    continue
                trimmed_seq = raw_seq[trim_bp:]
                stitched_seq += trimmed_seq
                current_cds_end = hit['qEnd']
            else:
                # Perfect abutment
                stitched_seq += raw_seq
                current_cds_end = hit['qEnd']

    final_dna = str(stitched_seq).upper()
    
    # QC: Translate to check for internal stops (Quality Check Only)
    aa_seq = Seq(final_dna).translate()
    internal_stops = 0
    if '*' in aa_seq[:-1]: # Ignore stop at very end
        internal_stops = aa_seq.count('*')
        log_decision(strain_id, gene_name, "Best_Copy", "QC", "Internal Stop Detected", 
                     f"{internal_stops} Stops found", "Potential Frameshift or Pseudogene")

    return [{
        "header": copy_id,
        "sequence": final_dna,
        "stats": {"gap_N_bp": gap_n_count, "internal_stops": internal_stops}
    }]

def main():
    parser = argparse.ArgumentParser(description="Y-Chromosome Pipeline v11 (DNA CDS Mode)")
    parser.add_argument("--genome", required=True, help="Genome FASTA (Assembly)")
    parser.add_argument("--psl", required=True, help="BLAT PSL (DNA-DNA)")
    parser.add_argument("--refs", required=True, help="Reference CDS FASTA (CarvalhoY_CDS.fasta)")
    parser.add_argument("--strain", required=True, help="Strain Name (e.g., A3, ISO1)")
    parser.add_argument("--outdir", default="Reconstructed_CDS_v11", help="Output Directory")
    
    args = parser.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    
    print(f"Loading genome: {args.genome}...")
    genome = SeqIO.to_dict(SeqIO.parse(args.genome, "fasta"))
    
    print(f"Loading CDS references: {args.refs}...")
    # Load CDS lengths for completeness calc
    refs = {r.id: len(r.seq) for r in SeqIO.parse(args.refs, "fasta")}
    
    # 1. Parse Everything
    all_hits = parse_psl_and_extract(args.psl, genome)
    
    # 2. Resolve Conflicts
    active_hits = resolve_conflicts(all_hits, args.strain)
    
    # 3. Process Survivors
    for gene, hit_list in active_hits.items():
        # Match gene name to reference (fuzzy match if needed, or exact)
        # We assume BLAT qName matches Reference FASTA header ID exactly
        if gene not in refs:
            # Try fuzzy matching if exact match fails
            candidates = [k for k in refs.keys() if gene in k]
            if not candidates:
                print(f"Warning: Gene {gene} found in PSL but not in Ref FASTA keys. Skipping.")
                continue
            ref_len = refs[candidates[0]]
        else:
            ref_len = refs[gene]
            
        print(f"  - Reconstructing CDS for {gene}...")
        results = find_single_best_copy(hit_list, ref_len, args.strain, gene)
        
        for res in results:
            out_file = os.path.join(args.outdir, f"{res['header']}.fasta")
            with open(out_file, "w") as f:
                f.write(f">{res['header']}\n{res['sequence']}\n")
            
            total_len = len(res['sequence'])
            completeness = 0
            if ref_len > 0:
                completeness = ((total_len - res['stats']['gap_N_bp']) / ref_len) * 100
            
            SUMMARY_LOG.append({
                "Filename": f"{res['header']}.fasta",
                "Strain": args.strain,
                "Gene": gene,
                "Ref_CDS_Length": ref_len,
                "Reconstructed_Length": total_len,
                "Gap_N_bp": res['stats']['gap_N_bp'],
                "Internal_Stops": res['stats']['internal_stops'],
                "Completeness_Pct": round(completeness, 2)
            })

    excel_path = os.path.join(args.outdir, f"{args.strain}_Pipeline_v11_Audit.xlsx")
    with pd.ExcelWriter(excel_path) as writer:
        pd.DataFrame(SUMMARY_LOG).to_excel(writer, sheet_name='Gene_Stats', index=False)
        pd.DataFrame(AUDIT_LOG).to_excel(writer, sheet_name='Detailed_Log', index=False)
        
    print(f"\nDone. Validation Report saved to {excel_path}")

if __name__ == "__main__":
    main()
