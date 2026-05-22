import os
import re
import csv
import glob
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from collections import defaultdict

# ==============================================================================
# CONFIGURATION
# ==============================================================================

# Target prefixes for the files in your current directory
TARGET_PREFIXES = ["5x", "15x", "30x"]

# Output Master Files
MASTER_CSV_OUTPUT = "All_Morphs_Structural_Summary.csv"
ABUNDANCE_TSV_OUTPUT = "Morph_Abundance_Summary.tsv"

# Sleek, Modern Color Palette for Subunits
COLOR_MAP = {
    '18S': '#2b83ba',       'ITS': '#abdda4',       '5.8S': '#fdae61',      
    '2S': '#ffd92f',        '28S': '#1a9641',       'ETS': '#a6d96a',       
    'NTS_240': '#d7191c',   'NTS_330': '#ea8285',   'NTS_95': '#6a3d9a',    
    'IGS': '#cab2d6',       'R1': '#8c510a',        'R2': '#d8b365',        
    'Jockey': '#f1b6da',    'Other': '#d3d3d3'      
}

# Pre-defined high-contrast RGB colors for IGV Tracks (Allows up to 20 distinct morphs)
# Morph 0 will be red, Morph 1 will be green, etc.
IGV_MORPH_COLORS = [
    "230,25,75",   "60,180,75",   "255,225,25",  "0,130,200", 
    "245,130,48",  "145,30,180",  "70,240,240",  "240,50,230", 
    "210,245,60",  "250,190,212", "0,128,128",   "220,190,255", 
    "154,99,36",   "255,250,200", "128,0,0",     "170,255,195", 
    "128,128,0",   "255,216,177", "0,0,117",     "128,128,128"
]

# ==============================================================================
# DATA PARSING & EXTRACTION FUNCTIONS
# ==============================================================================

def clean_morph_name(raw_name):
    """Extracts the Morph number to perfectly match the internal annotations."""
    match = re.search(r'morphconsensus(\d+)', raw_name)
    if match: return f"Morph {match.group(1)}"
    return raw_name.split('_')[0]

def parse_repeatmasker_morphs(filepath):
    morphs_data = defaultdict(lambda: {'length': 0, 'annotations': []})
    if not os.path.exists(filepath): return {}
    try:
        with open(filepath, 'r') as f:
            for _ in range(3): next(f) # Skip header
            for line in f:
                parts = line.strip().split()
                if not parts: continue
                try:
                    q_name, q_start, q_end, strand, rep_name = parts[4], int(parts[5])-1, int(parts[6]), parts[8], parts[9]
                except (IndexError, ValueError): continue
                
                if q_end > morphs_data[q_name]['length']: morphs_data[q_name]['length'] = q_end
                
                label = "Other"
                if "18S" in rep_name: label = "18S"
                elif "28S" in rep_name: label = "28S"
                elif "5_8S" in rep_name: label = "5.8S"
                elif "2s" in rep_name.lower(): label = "2S"
                elif "ITS" in rep_name: label = "ITS"
                elif "ETS" in rep_name: label = "ETS"
                elif "NTS_240" in rep_name: label = "NTS_240"
                elif "NTS_330" in rep_name: label = "NTS_330"
                elif "NTS_95" in rep_name: label = "NTS_95"
                elif "NTS" in rep_name or "IGS" in rep_name: label = "IGS"
                elif "R1" in rep_name: label = "R1"
                elif "R2" in rep_name: label = "R2"
                elif "Jockey" in rep_name: label = "Jockey"
                
                morphs_data[q_name]['annotations'].append({
                    'name': label, 'start': q_start, 'end': q_end,
                    'strand': strand, 'color': COLOR_MAP.get(label, COLOR_MAP['Other'])
                })
    except Exception as e:
        print(f"Error parsing {filepath}: {e}")
        return {}
    return morphs_data

def get_non_overlapping_length(annotations, target_name):
    feats = [a for a in annotations if a['name'] == target_name]
    if not feats: return 0
    feats.sort(key=lambda x: x['start'])
    merged = []
    for current in feats:
        if not merged: merged.append([current['start'], current['end']])
        else:
            prev = merged[-1]
            if current['start'] <= prev[1]: prev[1] = max(prev[1], current['end'])
            else: merged.append([current['start'], current['end']])
    return sum(end - start for start, end in merged)

def extract_te_insertions(annotations):
    tes = [ann for ann in annotations if ann['name'] in ['R1', 'R2', 'Jockey']]
    if not tes: return "None"
    tes_sorted = sorted(tes, key=lambda x: x['start'])
    return " + ".join([f"{t['name']}({t['end'] - t['start']}bp)" for t in tes_sorted])

def orient_and_phase_annotations(annotations, length):
    if not annotations: return annotations
    def get_true_start(feat_name):
        feats = sorted([a for a in annotations if a['name'] == feat_name], key=lambda x: x['start'])
        if not feats: return None
        if len(feats) == 1: return feats[0]['start']
        max_gap = -1
        true_start = feats[0]['start']
        for i in range(len(feats)):
            prev, curr = feats[i-1], feats[i]
            gap = (curr['start'] - prev['end']) % length
            if gap > max_gap: max_gap, true_start = gap, curr['start']
        return true_start

    pos_18 = get_true_start('18S')
    anchor_pos = get_true_start('5.8S') or get_true_start('2S') or get_true_start('28S')

    is_forward = True
    if pos_18 is not None and anchor_pos is not None:
        if (pos_18 - anchor_pos) % length < (anchor_pos - pos_18) % length:
            is_forward = False

    if not is_forward:
        annotations = [{'name': a['name'], 'start': length - a['end'], 'end': length - a['start'], 'strand': a['strand'], 'color': a['color']} for a in annotations]

    pos_18 = get_true_start('18S')
    shift_amount = pos_18 if pos_18 is not None else 0

    aligned_anns = []
    for ann in annotations:
        new_start = (ann['start'] - shift_amount) % length
        size = ann['end'] - ann['start']
        new_end = new_start + size
        if new_end > length:
            ann1, ann2 = ann.copy(), ann.copy()
            ann1['start'], ann1['end'] = new_start, length
            ann2['start'], ann2['end'] = 0, new_end % length
            aligned_anns.extend([ann1, ann2])
        else:
            ann_shifted = ann.copy()
            ann_shifted['start'], ann_shifted['end'] = new_start, new_end
            aligned_anns.append(ann_shifted)
    return aligned_anns

# ==============================================================================
# PLOTTING FUNCTION
# ==============================================================================

def plot_dense_morphs(morphs_data, condition_name, output_dir):
    if not morphs_data: return
    def get_sort_key(item):
        match = re.search(r'\d+', clean_morph_name(item[0]))
        return int(match.group()) if match else float('inf')

    sorted_morphs = sorted(morphs_data.items(), key=get_sort_key)
    num_morphs = len(sorted_morphs)
    
    fig, ax = plt.subplots(figsize=(14, max(4, num_morphs * 0.25)))
    legend_elements, max_x, y_labels, y_ticks = {}, 0, [], []

    for idx, (morph_id, data) in enumerate(sorted_morphs):
        length, anns = data['length'], data['annotations']
        y_pos = num_morphs - 1 - idx
        max_x = max(max_x, length)
            
        ax.plot([0, length], [y_pos, y_pos], color='#e0e0e0', linewidth=1.0, zorder=1)
        
        for ann in orient_and_phase_annotations(anns, length):
            ax.add_patch(patches.Rectangle((ann['start'], y_pos - 0.25), ann['end'] - ann['start'], 0.5, facecolor=ann['color'], edgecolor='none', zorder=2, alpha=0.9))
            if ann['name'] not in legend_elements:
                legend_elements[ann['name']] = patches.Patch(facecolor=ann['color'], label=ann['name'])

        y_ticks.append(y_pos)
        y_labels.append(clean_morph_name(morph_id))

    ax.set_ylim(-1, num_morphs)
    ax.set_xlim(0, max_x + 500)
    ax.set_yticks(y_ticks)
    ax.set_yticklabels(y_labels, fontsize=max(6, 12 - (num_morphs / 10)))
    ax.set_xlabel('Position (bp)', fontsize=12, fontweight='bold', color='#333333')
    for spine in ['top', 'right', 'left']: ax.spines[spine].set_visible(False)
    ax.spines['bottom'].set_color('#cccccc')
    ax.tick_params(axis='y', length=0)
    ax.set_title(f"rDNA Morphs: {condition_name} ({num_morphs} variants)", fontsize=14, fontweight='bold', pad=20, color='#333333')

    order = ['18S', 'ITS', '5.8S', '2S', '28S', 'ETS', 'NTS_240', 'NTS_330', 'NTS_95', 'IGS', 'R1', 'R2', 'Jockey', 'Other']
    ax.legend(handles=[legend_elements[k] for k in order if k in legend_elements], loc='center left', bbox_to_anchor=(1.02, 0.5), title="Subunits", frameon=False)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"{condition_name}_Morph_Structures.png"), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(output_dir, f"{condition_name}_Morph_Structures.svg"), format='svg', bbox_inches='tight')
    plt.close(fig)

# ==============================================================================
# MAIN EXECUTION
# ==============================================================================

if __name__ == "__main__":
    print("=======================================================================")
    print(" Starting Comprehensive Morph Analysis Pipeline")
    print("=======================================================================\n")

    # Create a single directory to hold all output txts, plots, and IGV tracks
    master_out_dir = "Annotations_Output"
    os.makedirs(master_out_dir, exist_ok=True)

    with open(MASTER_CSV_OUTPUT, 'w', newline='') as csvfile, \
         open(ABUNDANCE_TSV_OUTPUT, 'w') as tsvfile:
        
        # 1. Initialize CSV structural writer
        csv_writer = csv.writer(csvfile)
        csv_writer.writerow([
            "Condition", "Morph_ID", "Total_Length_bp", "TE_Insertions", 
            "18S_bp", "28S_bp", "NTS_240_bp", "NTS_330_bp", "NTS_95_bp", "IGS_bp"
        ])
        
        # 2. Initialize TSV abundance writer
        tsvfile.write("Condition\tMorph_ID\tFound_in_Assembly\tCopies_Mapped\n")
        
        for prefix in TARGET_PREFIXES:
            print(f">>> Processing {prefix} dataset...")
            
            # Predict RM out filename
            rm_file = f"{prefix}morphs.fa.out"
            if not os.path.exists(rm_file):
                rm_file = f"{prefix}morphs.out"
                if not os.path.exists(rm_file):
                    print(f"  -> WARNING: Could not find {prefix}morphs.fa.out. Skipping.")
                    continue
            
            # Predict Mapped BED filename
            mapped_bed_files = glob.glob(f"*{prefix}*verified_dual_pass.bed")
            
            morphs_data = parse_repeatmasker_morphs(rm_file)
            if not morphs_data:
                print(f"  -> No valid annotations found in {rm_file}.")
                continue
                
            txt_output = os.path.join(master_out_dir, f"{prefix}_Feature_Alignments.txt")
            
            # Track all unique morph IDs assembled by Ribotin
            assembled_morph_ids = []
            
            with open(txt_output, 'w') as txtfile:
                for morph_id, data in morphs_data.items():
                    length, raw_anns = data['length'], data['annotations']
                    aligned_anns = orient_and_phase_annotations(raw_anns, length)
                    clean_id = clean_morph_name(morph_id)
                    assembled_morph_ids.append(clean_id)
                    
                    txtfile.write(f"=== {prefix} | {clean_id} ({length} bp) ===\n")
                    for ann in sorted(aligned_anns, key=lambda x: x['start']):
                        txtfile.write(f"[{ann['start']}..{ann['end']}] {ann['name']} ({ann['end'] - ann['start']} bp)\n")
                    txtfile.write("\n")
                    
                    csv_writer.writerow([
                        prefix, clean_id, length, extract_te_insertions(aligned_anns), 
                        get_non_overlapping_length(aligned_anns, '18S'), get_non_overlapping_length(aligned_anns, '28S'),
                        get_non_overlapping_length(aligned_anns, 'NTS_240'), get_non_overlapping_length(aligned_anns, 'NTS_330'),
                        get_non_overlapping_length(aligned_anns, 'NTS_95'), get_non_overlapping_length(aligned_anns, 'IGS')
                    ])
            
            # Process Mapped BED file for Abundance & IGV Colorizing
            mapped_counts = defaultdict(int)
            
            if mapped_bed_files:
                original_bed = mapped_bed_files[0]
                igv_out_bed = os.path.join(master_out_dir, f"{prefix}_IGV_Colored_Track.bed")
                
                with open(original_bed, 'r') as in_bed, open(igv_out_bed, 'w') as out_bed:
                    # Write IGV Header for rendering 9-column itemRgb colors
                    out_bed.write(f'track name="{prefix} Morphs" description="Ribotin Morph Mappings" itemRgb="On" visibility="pack"\n')
                    
                    for line in in_bed:
                        parts = line.strip().split()
                        if len(parts) < 6: continue
                        chrom, start, end, raw_name, score, strand = parts[:6]
                        
                        clean_id = clean_morph_name(raw_name)
                        mapped_counts[clean_id] += 1
                        
                        # Extract Mapping Quality (Gold vs Relaxed)
                        map_status = "Gold" if "Gold" in raw_name else "Relaxed" if "Relaxed" in raw_name else ""
                        display_name = f"{clean_id.replace(' ', '_')}"
                        if map_status: display_name += f"_({map_status})"
                        
                        # Assign fixed RGB color based on Morph ID
                        match = re.search(r'\d+', clean_id)
                        m_num = int(match.group()) if match else 0
                        rgb_color = IGV_MORPH_COLORS[m_num % len(IGV_MORPH_COLORS)]
                        
                        # Write standard 9-column BED for IGV
                        # Columns: chrom, start, end, name, score, strand, thickStart, thickEnd, itemRgb
                        out_bed.write(f"{chrom}\t{start}\t{end}\t{display_name}\t{score}\t{strand}\t{start}\t{end}\t{rgb_color}\n")
                
                print(f"  -> Generated IGV-ready Colored BED file: {igv_out_bed}")
            else:
                print(f"  -> WARNING: Mapped BED file not found for {prefix}. Abundance counts will be 0.")
            
            # Write Abundance TSV
            for clean_id in assembled_morph_ids:
                copies = mapped_counts.get(clean_id, 0)
                found = "Yes" if copies > 0 else "No"
                tsvfile.write(f"{prefix}\t{clean_id}\t{found}\t{copies}\n")
            
            plot_dense_morphs(morphs_data, prefix, master_out_dir)
            print(f"  -> Successfully generated plots and summaries for {prefix}")

    print("\n=======================================================================")
    print(f" Pipeline Complete!")
    print(f"  - Structural Summary CSV:  {MASTER_CSV_OUTPUT}")
    print(f"  - Morph Abundance TSV:     {ABUNDANCE_TSV_OUTPUT}")
    print(f"  - IGV Colored BED Tracks:  Found in Annotations_Output/")
    print("=======================================================================")
