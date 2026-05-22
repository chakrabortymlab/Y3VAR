import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import seaborn as sns
import os

# Set global font sizes for axis ticks
plt.rcParams['xtick.labelsize'] = 21
plt.rcParams['ytick.labelsize'] = 21

# Extended Okabe-Ito Colorblind Palette for our baseline categories
c_sky = '#56B4E9'
c_blue = '#0072B2'
c_verm = '#D55E00'
c_green = '#009E73'
c_orange = '#E69F00'
c_purple = '#CC79A7'
c_yellow = '#F0E442'
c_black = '#000000'

def get_unit_type(name):
    # Hierarchy to grab specific structural variants
    if 'R1_r2_inserted' in name: return 'R1_R2'
    if 'R1_inserted' in name: return 'R1'
    if 'R2_inserted' in name: return 'R2'
    if 'I_dm_inserted' in name: return 'I_dm'
    if 'Canonical' in name: return 'Canonical'
    
    # If it is missing those keywords, return the actual string identifier so nothing is vaguely labeled
    parts = name.split('_')
    return parts[0] if len(parts) > 0 else name

def plot_strain_histograms(strain_name):
    matrix_file = f'Matrix_{strain_name}_Y_vs_{strain_name}_Y.tsv'
    if not os.path.exists(matrix_file):
        print(f"Warning: {matrix_file} not found. Skipping {strain_name}.")
        return

    df = pd.read_csv(matrix_file, sep='\t', index_col=0)
    rows = df.index.tolist()
    cols = df.columns.tolist()
    
    out_dir = 'Histograms'
    os.makedirs(out_dir, exist_ok=True)
    
    # 1. COMPILE ALL DATA AND DYNAMICALLY LABEL PAIRS
    all_data = []
    for i in range(len(rows)):
        for j in range(i + 1, len(cols)):
            val = df.iloc[i, j] * 100
            t1, t2 = get_unit_type(rows[i]), get_unit_type(cols[j])
            
            # Alphabetize the pair so "R1 vs Canonical" and "Canonical vs R1" are grouped exactly the same
            pair = sorted([t1, t2])
            combo_name = f"{pair[0]} vs {pair[1]}"
                
            all_data.append({
                'val': val,
                'Comparison': combo_name
            })
            
    df_all = pd.DataFrame(all_data)
    
    if df_all.empty:
        return

    # Extract every single unique combination found in the data (No hardcoded omissions)
    unique_combos = df_all['Comparison'].unique()
    
    # Map the primary combinations to our preferred palette
    base_palette = {
        'Canonical vs Canonical': c_sky,
        'R1 vs R1': c_green,
        'R2 vs R2': c_purple,
        'Canonical vs R1': c_orange,
        'Canonical vs R2': c_yellow,
        'R1 vs R2': c_verm
    }
    
    # Dynamically assign colors so any new combination (like I_dm vs R1) gets its own unique label and color
    palette_dict = {}
    extra_colors = sns.color_palette("tab20", len(unique_combos))
    color_idx = 0
    
    for combo in unique_combos:
        if combo in base_palette:
            palette_dict[combo] = base_palette[combo]
        else:
            palette_dict[combo] = extra_colors[color_idx]
            color_idx += 1

    # OVERALL COMBINED HISTOGRAM WITH EXTERNAL LEGEND
    plt.figure(figsize=(11, 6))
    
    sns.histplot(
        data=df_all, 
        x='val', 
        hue='Comparison', 
        bins=40, 
        kde=False, 
        palette=palette_dict, 
        multiple="stack", 
        edgecolor=c_black, 
        alpha=0.8,
        legend=False
    )
    plt.yscale('log')
    
    mean_val = np.mean(df_all['val'])
    median_val = np.median(df_all['val'])
    min_val, max_val = df_all['val'].min(), df_all['val'].max()
    
    plt.axvline(mean_val, color=c_blue, linestyle='solid', linewidth=2, label=f'Mean: {mean_val:.2f}%')
    plt.axvline(median_val, color=c_verm, linestyle='dashed', linewidth=2, label=f'Median: {median_val:.2f}%')
    
    padding = (max_val - min_val) * 0.05 if max_val != min_val else 1.0
    plt.xlim(min_val - padding, max_val + padding)
    
    plt.title(f'Overall Intra-Y Sequence Identity ({strain_name})\nAll Units Combined (n={len(df_all)})', fontsize=22, fontweight='bold')
    plt.xlabel('Sequence Identity (%)', fontsize=24)
    plt.ylabel('Count (Log)', fontsize=24)
    
    handles, labels = plt.gca().get_legend_handles_labels()
    # Add every dynamic combination that exists in the dataset to the legend
    for comp_name in sorted(unique_combos):
        handles.append(Patch(facecolor=palette_dict[comp_name], edgecolor=c_black, alpha=0.8, label=comp_name))
        labels.append(comp_name)
    
    # Push legend outside the plot box to the right
    plt.legend(handles=handles, labels=labels, fontsize=18, bbox_to_anchor=(1.05, 1), loc='upper left')
    
    overall_base = os.path.join(out_dir, f"{strain_name}_Intra_Y_Overall_Histogram")
    plt.savefig(f"{overall_base}.png", dpi=300, bbox_inches='tight')
    plt.savefig(f"{overall_base}.svg", format='svg', bbox_inches='tight')
    plt.close()

    # 2. SEPARATE INDIVIDUAL COMBO HISTOGRAMS (Generated for EVERY unique pairing)
    for combo_name in unique_combos:
        # Filter dataframe for just this specific combination
        combo_df = df_all[df_all['Comparison'] == combo_name]
        vals = combo_df['val'].tolist()
        
        if vals:
            plt.figure(figsize=(11, 6))
            sns.histplot(vals, bins=30, kde=False, color=palette_dict[combo_name], edgecolor=c_black, alpha=0.8)
            plt.yscale('log')
            
            mean_val = np.mean(vals)
            median_val = np.median(vals)
            min_val, max_val = min(vals), max(vals)
            
            plt.axvline(mean_val, color=c_blue, linestyle='solid', linewidth=2, label=f'Mean: {mean_val:.2f}%')
            plt.axvline(median_val, color=c_verm, linestyle='dashed', linewidth=2, label=f'Median: {median_val:.2f}%')
            
            padding = (max_val - min_val) * 0.05 if max_val != min_val else 1.0
            plt.xlim(min_val - padding, max_val + padding)
            
            plt.title(f"Intra-Y {combo_name} ({strain_name})\n(n={len(vals)})", fontsize=22, fontweight='bold')
            plt.xlabel('Sequence Identity (%)', fontsize=24)
            plt.ylabel('Count (Log)', fontsize=24)
            
            handles, labels = plt.gca().get_legend_handles_labels()
            handles.append(Patch(facecolor=palette_dict[combo_name], edgecolor=c_black, alpha=0.8, label=combo_name))
            labels.append(combo_name)
            
            plt.legend(handles=handles, labels=labels, fontsize=18, bbox_to_anchor=(1.05, 1), loc='upper left')
            
            # Replace spaces for safer file naming
            safe_combo_name = combo_name.replace(" ", "_")
            combo_base = os.path.join(out_dir, f"{strain_name}_Intra_Y_{safe_combo_name}_Histogram")
            plt.savefig(f"{combo_base}.png", dpi=300, bbox_inches='tight')
            plt.savefig(f"{combo_base}.svg", format='svg', bbox_inches='tight')
            plt.close()

# Run the plot generation for all three strains
for strain in ['ISO1', 'A3', 'A4']:
    plot_strain_histograms(strain)

print("All plots generated successfully with explicit dynamic labeling.")
