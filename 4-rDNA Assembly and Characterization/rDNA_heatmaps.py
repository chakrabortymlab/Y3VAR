import numpy as np
import pandas as pd
import re
import os
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
import matplotlib.patches as patches
from matplotlib.colors import LinearSegmentedColormap

# ==========================================
# 0. CONFIGURATION & CUSTOM GRADIENT
# ==========================================
min_identity_scale = 0.50  # Matches your 50% request

# Custom Okabe-Ito Gradient: Blue -> Sky -> Green -> Yellow -> Orange -> Deep Rust
okabe_colors = ["#0072B2", "#56B4E9", "#009E73", "#F0E442", "#E69F00", "#8C3E00"]
okabe_cmap = LinearSegmentedColormap.from_list("OkabeItoContinuous", okabe_colors, N=256)

# Expanded to distinguish X vs Y linked insertions
UNIT_TYPE_COLORS = {
    'Canonical': '#009E73',           # Bluish Green (Shared)
    
    # X-Linked (Cool Colors)
    'X_R1_inserted': '#56B4E9',       # Sky Blue
    'X_R2_inserted': '#0072B2',       # Blue
    'X_R1_R2_inserted': '#CC79A7',     # Reddish Purple
    
    # Y-Linked (Warm Colors)
    'Y_R1_inserted': '#E69F00',       # Orange
    'Y_R2_inserted': '#D55E00',       # Vermillion
    'Y_R1_R2_inserted': '#F0E442',     # Yellow
}

# ==========================================
# 1. PARSING FUNCTIONS
# ==========================================
def parse_fasta(file_path):
    print(f"Parsing FASTA: {file_path}")
    sequences = {}
    curr_id, curr_seq = None, []
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line: continue
            if line.startswith(">"):
                if curr_id: sequences[curr_id] = "".join(curr_seq).upper()
                curr_id, curr_seq = line[1:], []
            else:
                curr_seq.append(line)
    if curr_id: sequences[curr_id] = "".join(curr_seq).upper()
    return sequences

def parse_structural_bed(bed_path, target_contig, chrom_label):
    bed_entries = []
    try:
        with open(bed_path, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < 4 or "track" in parts[0]: continue
                chrom, start, end, raw_name = parts[0], int(parts[1]), int(parts[2]), parts[3]
                if target_contig.lower() not in chrom.lower(): continue
                
                clean_type = re.sub(r'(?i)^rDNA_Unit_', '', raw_name).capitalize()
                
                # Append Chromosome to inserted types to distinguish X vs Y
                if 'inserted' in clean_type.lower():
                    clean_type = f"{chrom_label}_{clean_type}"
                    
                bed_entries.append({'start': start, 'end': end, 'unit_type': clean_type})
    except FileNotFoundError:
        return []
    return sorted(bed_entries, key=lambda x: x['start'])

def sort_key_type(m):
    return (0, m) if 'Canonical' in m else (1, m)

# ==========================================
# 2. MAIN EXECUTION
# ==========================================
if __name__ == "__main__":
    fasta_file = "All_Strains_Aligned_FFTNS2_GAPs_1000xiterate.fasta"
    datasets = {
        'ISO1_X': ('ISO1XPhaseMatched_CLEAN.bed', 'CM090845.1'),
        'ISO1_Y': ('ISO1YPhaseMatched_CLEAN.bed', 'Y_Scaffold2'),
        'A3_X':   ('A3XPhaseMatched_CLEAN.bed', 'CM090833.1'),
        'A3_Y':   ('A3YPhaseMatched_CLEAN.bed', 'Y_scaffold2'),
        'A4_X':   ('A4XPhaseMatched_CLEAN.bed', 'CM090839.1'),
        'A4_Y':   ('A4YPhaseMatched_CLEAN.bed', 'Y_scaffold2')
    }

    out_dir = "Structural_Homogenization_Results"
    os.makedirs(out_dir, exist_ok=True)

    sequences = parse_fasta(fasta_file)
    seq_arrays = {sid: np.frombuffer(seq.encode('ascii'), dtype='S1') for sid, seq in sequences.items()}

    fasta_pool = {ds: {} for ds in datasets}
    for sid in sequences:
        match = re.search(r'^(ISO1|A3|A4)_([XY])_(.+)_(\d+)$', sid, re.IGNORECASE)
        if match:
            strain, chrom, raw_type = match.group(1).upper(), match.group(2).upper(), match.group(3).capitalize()
            
            # Append Chromosome to inserted types to distinguish X vs Y
            if 'inserted' in raw_type.lower():
                unit_type = f"{chrom}_{raw_type}"
            else:
                unit_type = raw_type
                
            dataset_key = f"{strain}_{chrom}"
            if dataset_key in fasta_pool:
                if unit_type not in fasta_pool[dataset_key]: fasta_pool[dataset_key][unit_type] = []
                fasta_pool[dataset_key][unit_type].append(sid)

    dataset_data = {}
    for dataset_key, (bed_file, contig) in datasets.items():
        # Pass the chromosome label ('X' or 'Y') to the BED parser
        chrom_label = dataset_key.split('_')[1]
        ordered_entries = parse_structural_bed(bed_file, target_contig=contig, chrom_label=chrom_label)
        if not ordered_entries: continue
        
        local_pool = {k: list(v) for k, v in fasta_pool[dataset_key].items()}
        labels, seqs, types = [], [], []
        for i, entry in enumerate(ordered_entries):
            u_type = entry['unit_type']
            if u_type in local_pool and len(local_pool[u_type]) > 0:
                matched_sid = local_pool[u_type].pop(0)
                labels.append(f"U{i+1}_{u_type}")
                seqs.append(seq_arrays[matched_sid])
                types.append(u_type)
        dataset_data[dataset_key] = {'labels': labels, 'seqs': seqs, 'types': types}

    keys = list(dataset_data.keys())
    for idx1 in range(len(keys)):
        for idx2 in range(idx1, len(keys)):
            key1, key2 = keys[idx1], keys[idx2]
            l1, s1, t1 = dataset_data[key1]['labels'], dataset_data[key1]['seqs'], dataset_data[key1]['types']
            l2, s2, t2 = dataset_data[key2]['labels'], dataset_data[key2]['seqs'], dataset_data[key2]['types']
            n1, n2 = len(l1), len(l2)
            if n1 < 2 or n2 < 2: continue

            print(f"Comparing {key1} vs {key2}...")
            sim_matrix = np.zeros((n1, n2))
            for i in range(n1):
                start_j = i if key1 == key2 else 0
                for j in range(start_j, n2):
                    v1, v2 = s1[i], s2[j]
                    valid = (v1 != b'-') & (v2 != b'-')
                    overlap = np.sum(valid)
                    ident = np.sum(v1[valid] == v2[valid]) / overlap if overlap > 100 else 0.0
                    sim_matrix[i, j] = ident
                    if key1 == key2: sim_matrix[j, i] = ident

             # --- PLOTTING ---
            if key1 == key2:
                fig, ax = plt.subplots(figsize=(18, 6))
                verts, facecolors = [], []
                for i in range(n1):
                    for j in range(i, n1):
                        x_c, y_c = (i+j)/2.0, (j-i)/2.0
                        verts.append([(x_c, y_c-0.5), (x_c-0.5, y_c), (x_c, y_c+0.5), (x_c+0.5, y_c)])
                        facecolors.append(sim_matrix[i, j])
                
                # Removed rasterized=True so SVGs are true vector format
                collection = PolyCollection(verts, cmap=okabe_cmap, edgecolors='white', linewidths=0.2, rasterized=False)
                collection.set_array(np.array(facecolors))
                collection.set_clim(vmin=min_identity_scale, vmax=1.0)
                ax.add_collection(collection)

                # Annotation Track
                track_colors = [UNIT_TYPE_COLORS.get(t, '#999999') for t in t1]
                ax.add_collection(PolyCollection([[(i-0.5,-1.2),(i-0.5,-0.5),(i+0.5,-0.5),(i+0.5,-1.2)] for i in range(n1)], 
                                               facecolors=track_colors, edgecolors='white', linewidths=0.5, rasterized=False))
                
                ax.autoscale_view(); ax.set_aspect('equal')
                ax.set_title(f"{key1} Intra-Array (Min Scale: {int(min_identity_scale*100)}%)", fontweight='bold')
                fig.colorbar(collection, ax=ax, shrink=0.6).set_label('Identity')
                
                # SAVE BOTH FORMATS
                plt.savefig(os.path.join(out_dir, f"Heatmap_INTRA_{key1}.png"), dpi=600, bbox_inches='tight')
                plt.savefig(os.path.join(out_dir, f"Heatmap_INTRA_{key1}.svg"), format='svg', bbox_inches='tight')
                plt.close()

            else:
                fig, ax = plt.subplots(figsize=(max(8, n1*0.2), max(8, n2*0.2)))
                verts = [[(i-0.5, j-0.5), (i-0.5, j+0.5), (i+0.5, j+0.5), (i+0.5, j-0.5)] for i in range(n1) for j in range(n2)]
                
                # Removed rasterized=True so SVGs are true vector format
                collection = PolyCollection(verts, cmap=okabe_cmap, edgecolors='white', linewidths=0.1, rasterized=False)
                collection.set_array(sim_matrix.flatten())
                collection.set_clim(vmin=min_identity_scale, vmax=1.0)
                ax.add_collection(collection)
                
                ax.autoscale_view(); ax.set_aspect('equal')
                ax.set_title(f"{key1} vs {key2}", fontweight='bold')
                
                # SAVE BOTH FORMATS
                plt.savefig(os.path.join(out_dir, f"Heatmap_INTER_{key1}_vs_{key2}.png"), dpi=600, bbox_inches='tight')
                plt.savefig(os.path.join(out_dir, f"Heatmap_INTER_{key1}_vs_{key2}.svg"), format='svg', bbox_inches='tight')
                plt.close()

    # ==========================================
    # 3. GENERATE UNIVERSAL UNIT TYPE LEGEND
    # ==========================================
    if UNIT_TYPE_COLORS:
        print("\nGenerating Universal Unit Type Legend...")
        legend_fig, legend_ax = plt.subplots(figsize=(8, len(UNIT_TYPE_COLORS) * 0.4))
        legend_ax.axis('off')
        
        # Only add to legend if it's not 'Unknown' or alternatives just to keep it clean
        legend_items = {k: v for k, v in UNIT_TYPE_COLORS.items() if k not in ['X_R1_r2_inserted', 'Y_R1_r2_inserted']}
        legend_patches = [patches.Patch(color=color, label=label) for label, color in legend_items.items()]
        
        ncols = 1 if len(legend_items) <= 6 else 2
        legend_ax.legend(handles=legend_patches, loc='center', title="Chromosome-Specific Unit Types", 
                         fontsize=12, title_fontsize=14, ncol=8, frameon=True)
        
        legend_out_png = os.path.join(out_dir, "Universal_Unit_Type_Legend.png")
        legend_out_svg = os.path.join(out_dir, "Universal_Unit_Type_Legend.svg")
        
        legend_fig.savefig(legend_out_png, dpi=600, bbox_inches='tight')
        legend_fig.savefig(legend_out_svg, format='svg', bbox_inches='tight')
        plt.close(legend_fig)
        print(f"Saved Universal Unit Type Legend to {legend_out_svg}")

    print("\nAnalysis complete with custom Okabe-Ito gradient and pure vector SVGs.")
