import sys
import os
from collections import defaultdict

# ==========================================
# CONFIGURATION
# ==========================================
PSL_FILE = "results_DelgadoYA4_CDS.psl"
GFF_FILE = "DelgadoYA4_Corrected_Final.gff3"

# Colors for IGV
GENE_COLORS = {
    "kl5": "#FF0000", "kl-5": "#FF0000",
    "kl3": "#00FF00", "kl-3": "#00FF00",
    "kl2": "#0000FF", "kl-2": "#0000FF",
    "PprY": "#FFA500", "Ppr-Y": "#FFA500",
    "WDY": "#800080",
    "ORY": "#00FFFF",
    "PRY": "#FF00FF",
    "CCY": "#A52A2A",
    "Pp1Y1": "#808080",
    "Pp1Y2": "#C0C0C0",
    "FDY": "#FFD700",
    "ARY": "#008080",
    "CG41561_PB": "#4B0082"
}

def get_color(name):
    for key in GENE_COLORS:
        if key in name:
            return GENE_COLORS[key]
    return "#000000"

def parse_psl(filename):
    """Parses PSL file into a list of hit dictionaries."""
    hits = []
    with open(filename, 'r') as f:
        for line in f:
            if not line[0].isdigit(): continue
            cols = line.strip().split('\t')
            if len(cols) < 21: continue
            
            # PSL Columns (0-based)
            # 0:matches, 8:strand, 9:qName, 13:tName, 15:tStart, 16:tEnd
            # 17:blockCount, 18:blockSizes, 20:tStarts
            
            hits.append({
                "gene": cols[9],
                "chrom": cols[13],
                "strand": cols[8],
                "g_start": int(cols[15]), # 0-based
                "g_end": int(cols[16]),
                "matches": int(cols[0]),
                "block_count": int(cols[17]),
                "block_sizes": [int(x) for x in cols[18].split(',') if x],
                "t_starts": [int(x) for x in cols[20].split(',') if x]
            })
    return hits

def overlaps(hit1, hit2):
    """Checks if two hits overlap significantly."""
    if hit1['chrom'] != hit2['chrom']: return False
    
    start = max(hit1['g_start'], hit2['g_start'])
    end = min(hit1['g_end'], hit2['g_end'])
    overlap = end - start
    
    if overlap <= 0: return False
    
    # Check if overlap is substantial (e.g., > 20% of the smaller hit)
    min_len = min(hit1['g_end']-hit1['g_start'], hit2['g_end']-hit2['g_start'])
    return overlap > (0.2 * min_len)

def process_genes(all_hits):
    """Applies gene-specific logic to filter and stitch hits."""
    
    # Group by Gene Name (e.g., all 'kl-5' hits together)
    # We normalize names: 'kl5' and 'kl-5' go to 'kl-5' group
    groups = defaultdict(list)
    for h in all_hits:
        base_name = h['gene']
        # Normalize common names
        if "kl5" in base_name or "kl-5" in base_name: key = "kl-5"
        elif "Ppr" in base_name: key = "Ppr-Y"
        elif "WDY" in base_name: key = "WDY"
        else: key = base_name
        groups[key].append(h)
        
    final_models = []
    
    for gene_key, hits in groups.items():
        # --- LOGIC PER GENE ---
        
        if gene_key == "Ppr-Y":
            # STRATEGY: Stitch, but enforced STRAND CONSISTENCY.
            # 1. Count score per strand
            plus_score = sum(h['matches'] for h in hits if h['strand'] == '+')
            minus_score = sum(h['matches'] for h in hits if h['strand'] == '-')
            
            # 2. Pick winner strand
            winner_strand = '+' if plus_score > minus_score else '-'
            
            # 3. Filter hits for winner strand ONLY
            valid_hits = [h for h in hits if h['strand'] == winner_strand]
            
            # 4. Stitch them all into ONE model
            # (Ppr-Y is naturally fragmented in BLAT, so we merge everything)
            merged_model = stitch_hits(valid_hits, gene_key, winner_strand)
            final_models.append(merged_model)
            
        elif gene_key == "kl-5":
            # STRATEGY: "The Highlander" (There can be only one... per locus)
            # 1. Sort by Score
            hits.sort(key=lambda x: x['matches'], reverse=True)
            
            # 2. Pick the Champion (Best Hit)
            champion = hits[0]
            
            # 3. Keep Champion. Check others.
            # If another hit overlaps the Champion, KILL IT.
            # If it doesn't overlap (e.g., a real copy 2Mb away), keep it.
            kept_hits = [champion]
            
            for candidate in hits[1:]:
                if not overlaps(champion, candidate):
                    kept_hits.append(candidate)
            
            # Convert kept hits to models
            for i, h in enumerate(kept_hits):
                model = hit_to_model(h, gene_key, i+1)
                final_models.append(model)
                
        elif gene_key == "WDY":
            # STRATEGY: Prioritize the 6-exon "Super Hit". Kill overlapping noise.
            hits.sort(key=lambda x: x['matches'], reverse=True)
            
            kept_hits = []
            for h in hits:
                # If hit is tiny (<500 score) AND overlaps a bigger kept hit, skip it
                is_noise = False
                if h['matches'] < 500:
                    for k in kept_hits:
                        if overlaps(h, k):
                            is_noise = True
                            break
                if not is_noise:
                    kept_hits.append(h)
            
            # WDY doesn't need stitching anymore (it's one hit in v10)
            for i, h in enumerate(kept_hits):
                model = hit_to_model(h, gene_key, i+1)
                final_models.append(model)
        
        else:
            # DEFAULT STRATEGY (Other genes)
            # Sort by score, keep all unique/non-overlapping hits
            hits.sort(key=lambda x: x['matches'], reverse=True)
            for i, h in enumerate(hits):
                model = hit_to_model(h, gene_key, i+1)
                final_models.append(model)

    return final_models

def stitch_hits(hits, gene_name, strand):
    """Merges multiple PSL hits into a single Gene Model."""
    if not hits: return None
    
    # Collect all exons from all hits
    all_exons = []
    total_score = 0
    chrom = hits[0]['chrom']
    
    for h in hits:
        total_score += h['matches']
        for i in range(h['block_count']):
            # PSL tStarts are 0-based. GFF needs 1-based.
            start = h['t_starts'][i] + 1
            size = h['block_sizes'][i]
            end = start + size - 1
            all_exons.append((start, end))
            
    # Sort exons
    all_exons.sort(key=lambda x: x[0])
    
    # Global coordinates
    g_start = all_exons[0][0]
    g_end = all_exons[-1][1]
    
    return {
        "gene": gene_name,
        "chrom": chrom,
        "strand": strand,
        "start": g_start,
        "end": g_end,
        "exons": all_exons,
        "score": total_score,
        "type": "stitched"
    }

def hit_to_model(h, gene_name, copy_num):
    """Converts a single PSL hit to a Model object."""
    exons = []
    for i in range(h['block_count']):
        start = h['t_starts'][i] + 1
        size = h['block_sizes'][i]
        end = start + size - 1
        exons.append((start, end))
        
    return {
        "gene": gene_name,
        "chrom": h['chrom'],
        "strand": h['strand'],
        "start": h['g_start'] + 1,
        "end": h['g_end'], # PSL end is exclusive? No, PSL tEnd is exclusive. Check parsing.
                           # Actually standard PSL `tEnd` is exclusive.
        "exons": exons,
        "score": h['matches'],
        "type": "single",
        "copy": copy_num
    }

def write_gff(models, filename):
    """Writes the models to GFF3."""
    models.sort(key=lambda x: (x['chrom'], x['start']))
    
    with open(filename, 'w') as out:
        out.write("##gff-version 3\n")
        
        # Track counts for naming
        counts = defaultdict(int)
        
        for m in models:
            if not m: continue
            
            gene = m['gene']
            counts[gene] += 1
            count = counts[gene]
            
            # Name
            if m['type'] == 'stitched':
                name = gene
                suffix = ""
            else:
                # If it's the main large hit, call it "Gene", else "Gene_CopyX"
                # Simple heuristic: If score > 2000, it's the main gene
                if m['score'] > 2000 and count == 1:
                    name = gene
                    suffix = ""
                else:
                    name = f"{gene}_Copy{count}"
                    suffix = f"_copy{count}"
            
            gene_id = f"{gene}{suffix}"
            mrna_id = f"{gene_id}.t1"
            color = get_color(gene)
            
            # Write Gene Line
            out.write(f"{m['chrom']}\tClean_GFF\tgene\t{m['start']}\t{m['end']}\t.\t{m['strand']}\t.\tID={gene_id};Name={name};color={color}\n")
            out.write(f"{m['chrom']}\tClean_GFF\tmRNA\t{m['start']}\t{m['end']}\t.\t{m['strand']}\t.\tID={mrna_id};Parent={gene_id};Name={name};color={color}\n")
            
            # Write Exons
            sorted_exons = sorted(m['exons'], key=lambda x: x[0])
            total_exons = len(sorted_exons)
            
            for i, (estart, eend) in enumerate(sorted_exons):
                # Strand numbering
                if m['strand'] == '-':
                    ex_num = total_exons - i
                else:
                    ex_num = i + 1
                    
                ex_id = f"{mrna_id}.exon{ex_num}"
                attrs = f"ID={ex_id};Parent={mrna_id};Name=Exon {ex_num};color={color}"
                
                out.write(f"{m['chrom']}\tClean_GFF\texon\t{estart}\t{eend}\t.\t{m['strand']}\t.\t{attrs}\n")
                out.write(f"{m['chrom']}\tClean_GFF\tCDS\t{estart}\t{eend}\t.\t{m['strand']}\t0\t{attrs}\n")

if __name__ == "__main__":
    hits = parse_psl(PSL_FILE)
    cleaned_models = process_genes(hits)
    write_gff(cleaned_models, GFF_FILE)
    print(f"Success! Corrected annotation saved to {GFF_FILE}")
