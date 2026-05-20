import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import sys
import os

# --- CONFIGURATION ---
INPUT_DEPTH_FILE = "raw_coverage_depth.txt"
OUTPUT_FILENAME = "HiFi_Coverage_Atlas_Fixed_YA3_Axis.pdf"

# Visualization Settings
WINDOW_SIZE = 2000       # 2kb windows (smooths the noise)
SCAFFOLDS_PER_PAGE = 4
PAGE_WIDTH, PAGE_HEIGHT = 20, 12

# Y-AXIS FIX FOR REPETITIVE GENOMES
# Sets the y-axis limit to 2.5x the mean coverage. 
# This focuses the visualization on the uniform regions while clipping extreme outliers.
Y_AXIS_MULTIPLIER = 2.5

# Filtering
MIN_SCAFFOLD_LEN = 20000 # Ignore contigs/scaffolds < 20kb

def load_and_bin_data(filepath, window_size):
    """
    Reads the massive samtools depth file, splits it by scaffold, 
    and averages the depth into windows (bins).
    """
    print(f"Reading {filepath} (this may take a moment)...")
    
    # Read columns: Chromosome, Position, Depth
    try:
        # We read the file in one go.
        df = pd.read_csv(filepath, sep='\t', names=['Chrom', 'Pos', 'Depth'])
    except Exception as e:
        print(f"Error reading file: {e}")
        sys.exit(1)

    print("Binning data (calculating sliding window averages)...")
    
    binned_data = []
    
    # Process each scaffold
    for chrom, group in df.groupby('Chrom'):
        max_pos = group['Pos'].max()
        if max_pos < MIN_SCAFFOLD_LEN:
            continue
            
        # Create bins
        group['Bin_ID'] = (group['Pos'] // window_size) * window_size
        
        # Calculate mean depth per bin
        binned = group.groupby('Bin_ID')['Depth'].mean().reset_index()
        binned['Chrom'] = chrom
        
        binned_data.append(binned)
        
    if not binned_data:
        print("No scaffolds passed the length filter.")
        sys.exit(1)

    return pd.concat(binned_data, ignore_index=True)

def plot_coverage(pdf, df):
    scaff_sizes = df.groupby('Chrom')['Bin_ID'].max()
    sorted_scaffolds = scaff_sizes.sort_values(ascending=False).index.tolist()
    
    # Calculate Global Average Coverage
    global_mean_cov = df['Depth'].mean()
    print(f"Global Mean Coverage: {global_mean_cov:.2f}x")

    for i in range(0, len(sorted_scaffolds), SCAFFOLDS_PER_PAGE):
        chunk = sorted_scaffolds[i : i + SCAFFOLDS_PER_PAGE]
        
        fig, axes = plt.subplots(nrows=len(chunk), ncols=1, figsize=(PAGE_WIDTH, PAGE_HEIGHT), constrained_layout=True)
        if len(chunk) == 1: axes = [axes]
        
        fig.suptitle(f"HiFi Read Coverage Atlas - Page {(i//SCAFFOLDS_PER_PAGE)+1}", fontsize=18, fontweight='bold')

        for ax, scaff_id in zip(axes, chunk):
            subset = df[df['Chrom'] == scaff_id]
            
            # Stats for this specific scaffold
            local_mean = subset['Depth'].mean()
            local_median = subset['Depth'].median()
            max_pos = subset['Bin_ID'].max()
            
            # --- PLOTTING ---
            ax.set_title(f"{scaff_id} (Length: ~{max_pos:,} bp)", fontsize=14, fontweight='bold', loc='left')
            
            # Plot the filled area
            ax.fill_between(subset['Bin_ID'], subset['Depth'], color='#1f77b4', alpha=0.6, label='Read Depth')
            ax.plot(subset['Bin_ID'], subset['Depth'], color='#1f77b4', linewidth=0.5)
            
            # Add Baseline (Global Mean)
            ax.axhline(global_mean_cov, color='red', linestyle='--', linewidth=1, label=f'Global Mean ({global_mean_cov:.0f}x)')
            
            # Formatting
            ax.set_ylabel("Depth (x)", fontsize=10)
            ax.set_xlabel("Position (bp)", fontsize=10)
            
            # DYNAMIC Y-AXIS LIMIT FIX: Cap the Y-axis based on the Y_AXIS_MULTIPLIER
            # Ensures normal regions are clearly visible by clipping extreme outliers
            y_limit = max(global_mean_cov * Y_AXIS_MULTIPLIER, local_mean * Y_AXIS_MULTIPLIER)
            ax.set_ylim(0, y_limit + 5) # Add 5x headroom for clarity
            ax.set_xlim(0, max_pos)
            
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.grid(axis='y', linestyle=':', alpha=0.5)
            
            # Add Stats Box
            stats_text = (
                f"Mean Depth: {local_mean:.1f}x\n"
                f"Median Depth: {local_median:.1f}x\n"
                f"Zero Cov: {(subset['Depth'] == 0).sum()} bins"
            )
            ax.text(0.98, 0.85, stats_text, transform=ax.transAxes, 
                    fontsize=10, verticalalignment='top', horizontalalignment='right',
                    bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="gray", alpha=0.9))

            ax.legend(loc='upper left', fontsize=8)

        pdf.savefig(fig)
        plt.close(fig)
        print(f"Saved Page {(i//SCAFFOLDS_PER_PAGE)+1}")

def main():
    if not os.path.exists(INPUT_DEPTH_FILE):
        print(f"Error: Could not find {INPUT_DEPTH_FILE}")
        print("ACTION REQUIRED: Run this command first:")
        print(f"samtools depth -a mapped.bam > {INPUT_DEPTH_FILE}")
        return

    # 1. Load and Bin
    binned_df = load_and_bin_data(INPUT_DEPTH_FILE, WINDOW_SIZE)
    
    # 2. Plot
    with PdfPages(OUTPUT_FILENAME) as pdf:
        plot_coverage(pdf, binned_df)
        
    print(f"\n[SUCCESS] Coverage Atlas saved to: {OUTPUT_FILENAME}")

if __name__ == "__main__":
    main()
