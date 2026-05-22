import numpy as np
import pandas as pd
import re
import os
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
import matplotlib.patches as patches

# ==========================================
# 0. CONFIGURATION & STRUCTURAL COLORS
# ==========================================
UNIT_TYPE_COLORS = {
    'Canonical': '#2ca02c',        # Green
    'R1_inserted': '#1f77b4',      # Blue
    'R2_inserted': '#d62728',      # Red
    'R1r2_inserted': '#9467bd',    # Purple
    'R1_r2_inserted': '#9467bd',   # Purple alternative
    'Unknown': '#7f7f7f'           # Grey for anything else
}

# ==========================================
# 1. PARSING FUNCTIONS
# ==========================================
def parse_fasta(file_path):
    print(f"Parsing FASTA: {file_path}")
    sequences = {}
    curr_id = None
    curr_seq = []
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line: continue
            if line.startswith(">"):
                if curr_id is not None:
                    sequences[curr_id] = "".join(curr_seq).upper()
                curr_id = line[1:]
                curr_seq = []
            else:
                curr_seq.append(line)
    if curr_id is not None:
        sequences[curr_id] = "".join(curr_seq).upper()
    return sequences

def parse_structural_bed(bed_path, target_contig):
    bed_entries = []
    try:
        with open(bed_path, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < 4 or "track" in parts[0]: continue
                
                chrom, start, end, raw_name = parts[0], int(parts[1]), int(parts[2]), parts[3]
                if target_contig.lower() not in chrom.lower(): continue
                
                clean_type = re.sub(r'(?i)^rDNA_Unit_', '', raw_name).capitalize()
                bed_entries.append({'start': start, 'end': end, 'unit_type': clean_type})
    except FileNotFoundError:
        print(f"[!] BED file not found: {bed_path}")
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

    print("Loading Sequence Data...")
    sequences = parse_fasta(fasta_file)
    seq_arrays = {sid: np.frombuffer(seq.encode('ascii'), dtype='S1') for sid, seq in sequences.items()}

    fasta_pool = { ds: {} for ds in datasets.keys() }
    for sid in sequences.keys():
        match = re.search(r'^(ISO1|A3|A4)_([XY])_(.+)_(\d+)$', sid, re.IGNORECASE)
        if match:
            strain, chrom, unit_type = match.group(1).upper(), match.group(2).upper(), match.group(3).capitalize()
            dataset_key = f"{strain}_{chrom}"
            if dataset_key in fasta_pool:
                if unit_type not in fasta_pool[dataset_key]: fasta_pool[dataset_key][unit_type] = []
                fasta_pool[dataset_key][unit_type].append(sid)

    # --- 2A. PRE-PROCESS ALL DATASETS ---
    print("\nPre-processing BED coordinates for all datasets...")
    dataset_data = {}
    
    for dataset_key, (bed_file, contig) in datasets.items():
        ordered_entries = parse_structural_bed(bed_file, target_contig=contig)
        if not ordered_entries: continue
            
        local_pool = {k: list(v) for k, v in fasta_pool[dataset_key].items()}
        labels, seqs, types = [], [], []
        
        for i, entry in enumerate(ordered_entries):
            u_type, start_coord = entry['unit_type'], entry['start']
            if u_type in local_pool and len(local_pool[u_type]) > 0:
                matched_sid = local_pool[u_type].pop(0)
                labels.append(f"Unit{i+1}_{u_type}_Pos{start_coord}")
                seqs.append(seq_arrays[matched_sid])
                types.append(u_type)
                
        dataset_data[dataset_key] = {'labels': labels, 'seqs': seqs, 'types': types}

    # --- 2B. CALCULATE & PLOT ALL COMBINATIONS ---
    keys = list(dataset_data.keys())
    
    for idx1 in range(len(keys)):
        for idx2 in range(idx1, len(keys)):
            key1, key2 = keys[idx1], keys[idx2]
            
            labels1, seqs1, types1 = dataset_data[key1]['labels'], dataset_data[key1]['seqs'], dataset_data[key1]['types']
            labels2, seqs2, types2 = dataset_data[key2]['labels'], dataset_data[key2]['seqs'], dataset_data[key2]['types']
            n1, n2 = len(labels1), len(labels2)
            
            if n1 < 2 or n2 < 2: continue
            
            print(f"\n==========================================")
            print(f"Comparing {key1} vs {key2}...")
            
            is_intra = (key1 == key2)
            sim_matrix = np.zeros((n1, n2))
            
            for i in range(n1):
                # For INTRA, we only calculate upper triangle to save time. For INTER, we calculate full matrix.
                start_j = i if is_intra else 0
                for j in range(start_j, n2):
                    s1, s2 = seqs1[i], seqs2[j]
                    valid = (s1 != b'-') | (s2 != b'-')
                    overlap_length = np.sum(valid)
                    ident = np.sum(s1[valid] == s2[valid]) / overlap_length if overlap_length > 100 else 0.0
                    sim_matrix[i, j] = ident
                    if is_intra: sim_matrix[j, i] = ident
            
            # Save TSV
            df_matrix = pd.DataFrame(sim_matrix, index=labels1, columns=labels2)
            tsv_path = os.path.join(out_dir, f"Matrix_{key1}_vs_{key2}.tsv")
            df_matrix.to_csv(tsv_path, sep='\t')
            
            # --- PLOTTING LOGIC ---
            unique_types = sorted(list(set(types1 + types2)), key=sort_key_type)
            
            if is_intra:
                # ========================================
                # INTRA-ARRAY PLOTS (Rotated Diamond + Line)
                # ========================================
                fig_heat, ax_heat = plt.subplots(figsize=(18, 6))
                verts, facecolors = [], []
                
                for i in range(n1):
                    for j in range(i, n1):
                        ident = sim_matrix[i, j]
                        x_center, y_center = (i + j) / 2.0, (j - i) / 2.0
                        verts.append([(x_center, y_center - 0.5), (x_center - 0.5, y_center), 
                                      (x_center, y_center + 0.5), (x_center + 0.5, y_center)])
                        facecolors.append(ident)

                facecolors = np.array(facecolors)
                min_ident = facecolors[facecolors > 0].min() if np.any(facecolors > 0) else 0.40
                collection = PolyCollection(verts, cmap='YlOrRd', edgecolors='white', linewidths=0.3)
                collection.set_array(facecolors)
                collection.set_clim(vmin=max(0.40, min_ident), vmax=1.0)
                ax_heat.add_collection(collection)
                
                # 1D Track
                track_verts, track_colors = [], []
                for i in range(n1):
                    c = UNIT_TYPE_COLORS.get(types1[i], UNIT_TYPE_COLORS['Unknown'])
                    track_verts.append([(i - 0.5, -1.2), (i - 0.5, -0.5), (i + 0.5, -0.5), (i + 0.5, -1.2)])
                    track_colors.append(c)
                ax_heat.add_collection(PolyCollection(track_verts, facecolors=track_colors, edgecolors='white', linewidths=0.5))
                
                ax_heat.autoscale_view()
                ax_heat.set_aspect('equal')
                for spine in ax_heat.spines.values(): spine.set_visible(False)
                ax_heat.spines['bottom'].set_visible(True)
                ax_heat.spines['bottom'].set_position(('data', -1.2))
                
                ticks = np.arange(0, n1, 10 if n1 >= 20 else 5)
                ax_heat.set_xticks(ticks)
                ax_heat.set_xticklabels(ticks + 1, fontsize=10)
                ax_heat.set_yticks([]) 
                ax_heat.set_xlabel(f"{key1} Array Position", fontsize=12, labelpad=10)
                ax_heat.set_title(f"{key1} Intra-Array Homogenization", fontsize=18, pad=20, fontweight='bold')
                fig_heat.colorbar(collection, ax=ax_heat, shrink=0.6, pad=0.02, aspect=15).set_label('Sequence Identity', rotation=270, labelpad=20)

                plt.savefig(os.path.join(out_dir, f"Heatmap_INTRA_{key1}.png"), dpi=600, bbox_inches='tight')
                plt.close(fig_heat)

                # LINE PLOT
                fig_line, ax_line = plt.subplots(figsize=(16, 4))
                adj_ident = np.diag(sim_matrix, k=1) * 100 
                x_line = np.arange(n1 - 1) + 0.5
                
                ax_line.plot(x_line, adj_ident, color='black', linewidth=1.5, zorder=3)
                ax_line.scatter(x_line, adj_ident, color='black', s=15, zorder=4)
                ax_line.axhline(99, color='red', linestyle='--', linewidth=1, zorder=2, alpha=0.7)
                
                for i in range(n1):
                    c = UNIT_TYPE_COLORS.get(types1[i], UNIT_TYPE_COLORS['Unknown'])
                    ax_line.add_patch(patches.Rectangle((i - 0.5, 0), 1, 10, linewidth=0, facecolor=c, zorder=1))
                    ax_line.add_patch(patches.Rectangle((i - 0.5, 10), 1, 95, linewidth=0, facecolor=c, alpha=0.15, zorder=0))
                    
                ax_line.set_ylim(0, 105)
                ax_line.set_xlim(-0.5, n1 - 0.5) 
                ax_line.set_ylabel("Identity (%)", fontsize=12, fontweight='bold')
                ax_line.set_xticks(ticks)
                ax_line.set_xticklabels([f"U{t+1}" for t in ticks], fontsize=10)
                ax_line.set_title(f"{key1} Adjacent Unit Identity", fontsize=14, fontweight='bold')

                plt.savefig(os.path.join(out_dir, f"LinePlot_{key1}.png"), dpi=600, bbox_inches='tight')
                plt.close(fig_line)

            else:
                # ========================================
                # INTER-ARRAY PLOTS (Rectangular Grid)
                # ========================================
                # Scale figure size relative to array lengths to keep squares strictly square
                fig_inter, ax_inter = plt.subplots(figsize=(max(8, n1*0.2), max(8, n2*0.2)))
                verts, facecolors = [], []
                
                for i in range(n1):
                    for j in range(n2):
                        verts.append([(i-0.5, j-0.5), (i-0.5, j+0.5), (i+0.5, j+0.5), (i+0.5, j-0.5)])
                        facecolors.append(sim_matrix[i, j])

                facecolors = np.array(facecolors)
                min_ident = facecolors[facecolors > 0].min() if np.any(facecolors > 0) else 0.40
                collection = PolyCollection(verts, cmap='YlOrRd', edgecolors='white', linewidths=0.1)
                collection.set_array(facecolors)
                collection.set_clim(vmin=max(0.40, min_ident), vmax=1.0)
                ax_inter.add_collection(collection)
                
                # Bottom Track (Key1 Types)
                track_verts_bot, track_colors_bot = [], []
                for i in range(n1):
                    c = UNIT_TYPE_COLORS.get(types1[i], UNIT_TYPE_COLORS['Unknown'])
                    track_verts_bot.append([(i - 0.5, -1.2), (i - 0.5, -0.6), (i + 0.5, -0.6), (i + 0.5, -1.2)])
                    track_colors_bot.append(c)
                ax_inter.add_collection(PolyCollection(track_verts_bot, facecolors=track_colors_bot, edgecolors='white', linewidths=0.5))

                # Left Track (Key2 Types)
                track_verts_left, track_colors_left = [], []
                for j in range(n2):
                    c = UNIT_TYPE_COLORS.get(types2[j], UNIT_TYPE_COLORS['Unknown'])
                    track_verts_left.append([(-1.2, j - 0.5), (-1.2, j + 0.5), (-0.6, j + 0.5), (-0.6, j - 0.5)])
                    track_colors_left.append(c)
                ax_inter.add_collection(PolyCollection(track_verts_left, facecolors=track_colors_left, edgecolors='white', linewidths=0.5))

                ax_inter.autoscale_view()
                ax_inter.set_aspect('equal')
                for spine in ax_inter.spines.values(): spine.set_visible(False)
                
                ticks_x = np.arange(0, n1, 10 if n1 >= 20 else 5)
                ax_inter.set_xticks(ticks_x)
                ax_inter.set_xticklabels(ticks_x + 1, fontsize=10)
                ax_inter.set_xlabel(f"{key1} Array Position", fontsize=14, labelpad=15)
                
                ticks_y = np.arange(0, n2, 10 if n2 >= 20 else 5)
                ax_inter.set_yticks(ticks_y)
                ax_inter.set_yticklabels(ticks_y + 1, fontsize=10)
                ax_inter.set_ylabel(f"{key2} Array Position", fontsize=14, labelpad=15)
                
                ax_inter.set_title(f"Inter-Array Homogenization: {key1} vs {key2}", fontsize=18, pad=20, fontweight='bold')
                fig_inter.colorbar(collection, ax=ax_inter, shrink=0.5, pad=0.03, aspect=20).set_label('Sequence Identity', rotation=270, labelpad=20)

                plt.savefig(os.path.join(out_dir, f"Heatmap_INTER_{key1}_vs_{key2}.png"), dpi=600, bbox_inches='tight')
                plt.close(fig_inter)

    # PANEL C: LEGEND (Drawn once at the end)
    all_types_in_run = []
    for d in dataset_data.values(): all_types_in_run.extend(d['types'])
    unique_types = sorted(list(set(all_types_in_run)), key=sort_key_type)
    
    if unique_types:
        fig_leg = plt.figure(figsize=(10, 1.5))
        ax_leg = fig_leg.add_subplot(111)
        ax_leg.axis('off') 
        ax_leg.legend(handles=[patches.Patch(facecolor=UNIT_TYPE_COLORS.get(m, UNIT_TYPE_COLORS['Unknown']), label=m) for m in unique_types], 
                      loc='center', title=f"Structural Variants", title_fontsize=14, fontsize=12, ncol=min(8, len(unique_types)), frameon=False)
        plt.savefig(os.path.join(out_dir, f"Structural_Legend_Universal.png"), dpi=600, bbox_inches='tight')
        plt.close(fig_leg)

print("\nAll structural analysis complete.")
