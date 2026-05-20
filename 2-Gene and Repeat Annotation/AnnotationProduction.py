import sys
import os
import csv
from collections import defaultdict

# ==========================================
# CONFIGURATION
# ==========================================
PSL_FILE = "results_DelgadoYISO-1v9_CDS.psl"
GFF_FILE = "DelgadoY_BLATPSL_Stitch.gff3"

# ONLY stitch Ppr-Y. All others are treated as distinct copies.
STITCH_TARGETS = ["PprY", "Ppr-Y"]

# Colors
GENE_COLORS = {
    "kl5": "#FF0000", "kl-5": "#FF0000",
    "kl3": "#00FF00", "kl-3": "#00FF00",
    "kl2": "#0000FF", "kl-2": "#0000FF",
    "PprY": "#FFA500", "Ppr-Y": "#FFA500", # Orange for Ppr-Y
    "WDY": "#800080",
    "ORY": "#00FFFF",
    "PRY": "#FF00FF",
    "CCY": "#A52A2A",
    "Pp1Y1": "#808080", "Pp1-Y1": "#808080",
    "Pp1Y2": "#C0C0C0", "Pp1-Y2": "#C0C0C0",
    "FDY": "#FFD700",
    "ARY": "#008080",
    "CG41561_PB": "#4B0082"
}

def get_color(name):
    for key in GENE_COLORS:
        if key in name:
            return GENE_COLORS[key]
    return "#000000"

def parse_psl_line(cols):
    """Parses a PSL line into a structured dictionary."""
    return {
        "gene": cols[9],
        "chrom": cols[13],
        "strand": cols[8],
        "g_start": int(cols[15]) + 1,
        "g_end": int(cols[16]),
        "matches": int(cols[0]),
        "block_count": int(cols[17]),
        "block_sizes": [int(x) for x in cols[18].split(',') if x],
        "t_starts": [int(x) for x in cols[20].split(',') if x]
    }

def run_conversion():
    if not os.path.exists(PSL_FILE):
        print(f"Error: {PSL_FILE} not found.")
        return

    print(f"Reading {PSL_FILE}...")
    
    # Store hits
    grouped_hits = defaultdict(list)
    
    with open(PSL_FILE, 'r') as f:
        for line in f:
            if not line[0].isdigit(): continue 
            cols = line.strip().split('\t')
            if len(cols) < 21: continue
            
            hit = parse_psl_line(cols)
            
            # Filter Noise
            if hit['matches'] < 50: continue

            # Group by Gene, Chrom, Strand
            key = (hit['gene'], hit['chrom'], hit['strand'])
            grouped_hits[key].append(hit)

    print(f"Processing {len(grouped_hits)} gene loci...")
    
    final_entries = []

    for (gene, chrom, strand), hits in grouped_hits.items():
        
        # CHECK: Is this Ppr-Y?
        is_ppry = any(target in gene for target in STITCH_TARGETS)
        
        if is_ppry:
            # --- STITCHING LOGIC (Ppr-Y ONLY) ---
            # Merge all fragments into ONE gene model
            print(f"  > Stitching {len(hits)} fragments for {gene} on {chrom}...")
            
            all_exons = []
            total_score = 0
            
            for h in hits:
                total_score += h['matches']
                for i in range(h['block_count']):
                    b_start = h['t_starts'][i] + 1
                    b_len = h['block_sizes'][i]
                    b_end = b_start + b_len - 1
                    all_exons.append((b_start, b_end))
            
            # Sort exons genomically
            all_exons.sort(key=lambda x: x[0])
            
            if not all_exons: continue
            
            # Create ONE entry
            final_entries.append({
                "type": "stitched",
                "gene": gene,
                "chrom": chrom,
                "strand": strand,
                "start": all_exons[0][0],
                "end": all_exons[-1][1],
                "exons": all_exons,
                "score": total_score
            })
            
        else:
            # --- STANDARD LOGIC (Everything Else) ---
            # Treat every hit as a distinct entity/copy
            
            # Sort by score (best first)
            hits.sort(key=lambda x: x['matches'], reverse=True)
            
            for h in hits:
                exons = []
                for j in range(h['block_count']):
                    b_start = h['t_starts'][j] + 1
                    b_len = h['block_sizes'][j]
                    b_end = b_start + b_len - 1
                    exons.append((b_start, b_end))
                
                final_entries.append({
                    "type": "single",
                    "gene": gene,
                    "chrom": chrom,
                    "strand": strand,
                    "start": h['g_start'],
                    "end": h['g_end'],
                    "exons": exons,
                    "score": h['matches']
                })

    # Sort final list genomically
    final_entries.sort(key=lambda x: (x['chrom'], x['start']))

    print(f"Writing GFF3 to {GFF_FILE}...")
    
    # Counter for ID generation
    gene_counts = defaultdict(int)

    with open(GFF_FILE, 'w') as out:
        out.write("##gff-version 3\n")
        
        for e in final_entries:
            gene_base = e['gene']
            gene_counts[gene_base] += 1
            copy_num = gene_counts[gene_base]
            
            # Naming
            if e['type'] == "stitched":
                # Ppr-Y gets the "Complete" label
                display_name = f"{gene_base} (Complete)"
                note = "Note=Stitched_PprY_Model;"
            else:
                # Others get Copy numbers
                if copy_num == 1:
                    display_name = f"{gene_base}"
                    suffix = ""
                else:
                    display_name = f"{gene_base}_Copy{copy_num}"
                    suffix = f"_copy{copy_num}"
                note = ""

            gene_id = f"{gene_base}_id{copy_num}"
            mrna_id = f"{gene_id}.t1"
            color = get_color(gene_base)
            
            # Write Gene & mRNA
            out.write(f"{e['chrom']}\tBLAT_Specific\tgene\t{e['start']}\t{e['end']}\t.\t{e['strand']}\t.\tID={gene_id};Name={display_name};color={color};{note}\n")
            out.write(f"{e['chrom']}\tBLAT_Specific\tmRNA\t{e['start']}\t{e['end']}\t.\t{e['strand']}\t.\tID={mrna_id};Parent={gene_id};Name={display_name};color={color}\n")
            
            # Write Exons
            # Sort exons by genomic position
            sorted_exons = sorted(e['exons'], key=lambda x: x[0])
            num_exons = len(sorted_exons)
            
            for i, (estart, eend) in enumerate(sorted_exons):
                # Numbering based on Strand
                if e['strand'] == '-':
                    ex_num = num_exons - i
                else:
                    ex_num = i + 1
                    
                ex_id = f"{mrna_id}.exon{ex_num}"
                attrs = f"ID={ex_id};Parent={mrna_id};Name=Exon {ex_num};color={color}"
                
                out.write(f"{e['chrom']}\tBLAT_Specific\texon\t{estart}\t{eend}\t.\t{e['strand']}\t.\t{attrs}\n")
                out.write(f"{e['chrom']}\tBLAT_Specific\tCDS\t{estart}\t{eend}\t.\t{e['strand']}\t0\t{attrs}\n")

    print(f"Done. {GFF_FILE} created.")
    print("Only Ppr-Y was stitched. All other genes are preserved as distinct BLAT hits.")

if __name__ == "__main__":
    run_conversion()
