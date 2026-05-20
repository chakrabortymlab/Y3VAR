import pandas as pd
import pysam
import re
import os

# ==========================================
# CONFIGURATION - CHANGE THESE TO YOUR FILES
# ==========================================
GFF_FILE = "DelgadoYISO1_Corrected_Final.gff3" # Your complete ISO1 reference

# Dictionary of your BAM files mapping to the Strain name
# Add your exact file paths here!
BAM_FILES = {
    "A4_Male": "A4M_vs_ISO1_sorted.bam",
    "A4_Female": "A4F_vs_ISO1_sorted.bam",
    "BL156_Male": "BL156_vs_ISO1_sorted.bam"
}

OUTPUT_TSV = "Supplementary_Table_Exon_Support_with_Females.tsv"

SUPPORT_THRESHOLD = 1.0  # Average depth > 1x

def get_primary_transcripts(gff_path):
    """Pass 1: Identifies the transcript with the most complete exon set for each gene."""
    gene_transcripts = {}
    
    with open(gff_path, 'r') as f:
        for line in f:
            if line.startswith("#"): continue
            cols = line.strip().split('\t')
            if len(cols) < 9 or cols[2] != 'exon': continue
            
            attr = cols[8]
            name_match = re.search(r'Name=Exon\s*(\d+)', attr, re.IGNORECASE)
            parent_match = re.search(r'Parent=([^;]+)', attr)
            
            if name_match and parent_match:
                t_id = parent_match.group(1)
                gene_name = t_id.split('.t')[0]
                
                # Standardize gene families to group fragments together
                if 'kl-5' in gene_name or 'kl5' in gene_name: base_gene = 'kl-5'
                elif 'kl-3' in gene_name or 'kl3' in gene_name: base_gene = 'kl-3'
                elif 'kl-2' in gene_name or 'kl2' in gene_name: base_gene = 'kl-2'
                elif 'Ppr' in gene_name: base_gene = 'Ppr-Y'
                elif 'ORY' in gene_name: base_gene = 'ORY'
                elif 'PRY' in gene_name: base_gene = 'PRY'
                elif 'WDY' in gene_name: base_gene = 'WDY'
                elif 'Pp1Y1' in gene_name: base_gene = 'Pp1-Y1'
                elif 'Pp1Y2' in gene_name: base_gene = 'Pp1-Y2'
                else:
                    # Strip out _copyX for others (like CCY) to find the main one
                    base_gene = re.sub(r'_copy\d+', '', gene_name)
                
                exon_num = int(name_match.group(1))
                
                if base_gene not in gene_transcripts:
                    gene_transcripts[base_gene] = {}
                if t_id not in gene_transcripts[base_gene]:
                    gene_transcripts[base_gene][t_id] = set()
                    
                gene_transcripts[base_gene][t_id].add(exon_num)
                
    primary_transcripts = {}
    for base_gene, transcripts in gene_transcripts.items():
        # Select the transcript ID that has the highest number of unique exons
        best_t_id = max(transcripts.keys(), key=lambda k: len(transcripts[k]))
        primary_transcripts[best_t_id] = base_gene
        
    return primary_transcripts

def parse_clean_exons(gff_path, primary_transcripts):
    """Pass 2: Extracts only the exons belonging to the primary transcripts."""
    exons = []
    
    with open(gff_path, 'r') as f:
        for line in f:
            if line.startswith("#"): continue
            cols = line.strip().split('\t')
            if len(cols) < 9 or cols[2] != 'exon': continue
            
            scaffold = cols[0]
            start = int(cols[3])
            end = int(cols[4])
            attr = cols[8]
            
            name_match = re.search(r'Name=Exon\s*(\d+)', attr, re.IGNORECASE)
            parent_match = re.search(r'Parent=([^;]+)', attr)
            
            if name_match and parent_match:
                t_id = parent_match.group(1)
                
                # ONLY KEEP IT if it belongs to one of our chosen primary models
                if t_id in primary_transcripts:
                    base_gene = primary_transcripts[t_id]
                    exon_num = int(name_match.group(1))
                    
                    exons.append({
                        'Gene': base_gene,
                        'Exon': exon_num,
                        'Scaffold': scaffold,
                        'Start': start,
                        'End': end,
                        'Length': end - start + 1
                    })
                    
    # Drop any weird duplicates just in case the same primary transcript listed an exon twice
    df = pd.DataFrame(exons).drop_duplicates(subset=['Gene', 'Exon'])
    return df

def add_bam_support(exons_df, strain_name, bam_path):
    print(f"Processing {strain_name} BAM file: {bam_path}...")
    
    if not os.path.exists(bam_path):
        print(f"  -> ERROR: {bam_path} not found! Skipping {strain_name}.")
        exons_df[f'{strain_name}_Supported'] = "ERROR"
        return exons_df

    if not os.path.exists(bam_path + ".bai"):
        print("  -> BAM index not found. Creating index...")
        pysam.index(bam_path)
        
    bam = pysam.AlignmentFile(bam_path, "rb")
    avg_cov_list, support_list = [], []
    
    for _, row in exons_df.iterrows():
        scaffold = row['Scaffold']
        start = row['Start'] - 1 
        end = row['End']
        
        try:
            cov_stats = bam.count_coverage(scaffold, start, end)
            total_bases = sum([sum(base_cov) for base_cov in cov_stats])
            avg_cov = total_bases / row['Length']
        except ValueError:
            avg_cov = 0.0
            
        avg_cov_list.append(round(avg_cov, 2))
        support_list.append("Yes" if avg_cov >= SUPPORT_THRESHOLD else "No")
            
    bam.close()
    
    exons_df[f'{strain_name}_Avg_Depth'] = avg_cov_list
    exons_df[f'{strain_name}_Supported'] = support_list
    return exons_df

if __name__ == "__main__":
    print("Pass 1: Identifying primary gene models...")
    primary_t_ids = get_primary_transcripts(GFF_FILE)
    print(f"  -> Found complete models for {len(primary_t_ids)} genes.")
    
    print("Pass 2: Extracting clean exons...")
    results_df = parse_clean_exons(GFF_FILE, primary_t_ids)
    
    # Sort logically
    results_df = results_df.sort_values(by=['Gene', 'Exon'])
    print(f"  -> Extracted {len(results_df)} unique exons.")
    
    for strain, bam_file in BAM_FILES.items():
        results_df = add_bam_support(results_df, strain, bam_file)
    
    cols = ['Gene', 'Exon', 'Scaffold', 'Start', 'End', 'Length']
    for strain in BAM_FILES.keys():
        cols.extend([f'{strain}_Avg_Depth', f'{strain}_Supported'])
    
    results_df = results_df[cols]
    
    print(f"\nSaving clean supplementary table to {OUTPUT_TSV}...")
    results_df.to_csv(OUTPUT_TSV, sep='\t', index=False)
    print("Done!")
