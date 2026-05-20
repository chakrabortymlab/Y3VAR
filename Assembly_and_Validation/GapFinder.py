import sys
import os
import matplotlib.pyplot as plt

# --- CONFIGURATION ---
INPUT_FILE = "raw_coverage_depth.txt"
OUTPUT_CSV = "A3Supplementary_Table_ZeroCoverageGaps.csv"
OUTPUT_PDF = "A3Supplementary_Figure_GapDistribution.pdf"

def process_gaps(filepath):
    print(f"Scanning {filepath} for zero-coverage regions...")
    
    gaps = []
    in_gap = False
    chrom = ""
    start = 0
    end = 0

    # 1. Read line-by-line (Memory safe for massive files)
    with open(filepath, 'r') as fin:
        for line in fin:
            parts = line.strip().split('\t')
            if len(parts) != 3: continue
                
            c, p, d = parts[0], int(parts[1]), int(parts[2])
            
            if d == 0:
                if not in_gap:
                    chrom, start, in_gap = c, p, True
                elif c != chrom:
                    gaps.append((chrom, start, end, end - start + 1))
                    chrom, start = c, p
                end = p
            else:
                if in_gap:
                    gaps.append((chrom, start, end, end - start + 1))
                    in_gap = False

    if in_gap:
        gaps.append((chrom, start, end, end - start + 1))

    return gaps

def write_reviewer_table(gaps, out_csv):
    """Writes a CSV with a metadata header."""
    total_gaps = len(gaps)
    if total_gaps == 0:
        print("No gaps found. Assembly has 100% coverage support!")
        sys.exit(0)

    sizes = [g[3] for g in gaps]
    total_gap_bases = sum(sizes)
    max_gap = max(sizes)
    avg_gap = total_gap_bases / total_gaps

    print(f"Writing documented data to {out_csv}...")
    
    with open(out_csv, 'w') as fout:
        # --- REVIEWER METADATA HEADER ---
        fout.write("# SUPPLEMENTARY DATA: Zero-Coverage Regions\n")
        fout.write("# Description: Genomic coordinates lacking primary read support (0x coverage).\n")
        fout.write(f"# Total Gaps Identified: {total_gaps:,}\n")
        fout.write(f"# Total Bases in Gaps: {total_gap_bases:,} bp\n")
        fout.write(f"# Average Gap Size: {avg_gap:.1f} bp\n")
        fout.write(f"# Maximum Gap Size: {max_gap:,} bp\n")
        fout.write("# \n")
        fout.write("# Column Definitions:\n")
        fout.write("# Scaffold_ID: Name of the contig/scaffold\n")
        fout.write("# Start_bp: 1-based start coordinate of the gap\n")
        fout.write("# End_bp: 1-based end coordinate of the gap\n")
        fout.write("# Gap_Length_bp: Total length of the zero-coverage region\n")
        fout.write("# --------------------------------------------------\n")
        
        # --- TABULAR DATA ---
        fout.write("Scaffold_ID,Start_bp,End_bp,Gap_Length_bp\n")
        for g in gaps:
            fout.write(f"{g[0]},{g[1]},{g[2]},{g[3]}\n")

def plot_gap_distribution(gaps, out_pdf):
    """Generates a histogram of gap sizes."""
    sizes = [g[3] for g in gaps]
    
    print(f"Generating supplementary figure: {out_pdf}...")
    
    plt.figure(figsize=(8, 6))
    
    # We use a log-scale for bins because gap sizes often vary from 1bp to 10,000+bp
    # This prevents massive gaps from squishing the entire plot
    bins = plt.hist(sizes, bins=50, color='darkred', edgecolor='black', alpha=0.7, log=True)
    
    plt.title("Distribution of Zero-Coverage Gap Sizes", fontsize=14, fontweight='bold')
    plt.xlabel("Gap Size (bp)", fontsize=12)
    plt.ylabel("Frequency (Log Scale)", fontsize=12)
    
    # Add summary text box to the plot
    stats_text = (
        f"Total Gaps: {len(sizes):,}\n"
        f"Total Missing bp: {sum(sizes):,}\n"
        f"Mean Size: {sum(sizes)/len(sizes):.1f} bp"
    )
    plt.text(0.95, 0.95, stats_text, transform=plt.gca().transAxes, 
             fontsize=10, verticalalignment='top', horizontalalignment='right',
             bbox=dict(boxstyle="round,pad=0.4", fc="white", ec="gray", alpha=0.9))
    
    plt.grid(axis='y', linestyle='--', alpha=0.5)
    plt.tight_layout()
    plt.savefig(out_pdf, dpi=300)
    plt.close()

if __name__ == "__main__":
    if not os.path.exists(INPUT_FILE):
        print(f"Error: Could not find {INPUT_FILE}")
        sys.exit(1)
        
    found_gaps = process_gaps(INPUT_FILE)
    write_reviewer_table(found_gaps, OUTPUT_CSV)
    plot_gap_distribution(found_gaps, OUTPUT_PDF)
    print("\n✅ Success! Both supplementary files are ready.")
