import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import os
import re  


OKABE_ITO = [
    "#009E73",  
    "#56B4E9",  
    "#E69F00",  
    "#F0E442", 
    "#0072B2",  
    "#D55E00", 
    "#CC79A7",  
    "#999999"   
]

def create_depth_map(depth_df, assembled_array_size):
    """Creates a fast NumPy array for base-pair level math."""
    depth_map = np.zeros(assembled_array_size + 1)
    valid_pos = depth_df['Pos'][depth_df['Pos'] <= assembled_array_size]
    valid_depth = depth_df['Depth'][depth_df['Pos'] <= assembled_array_size]
    depth_map[valid_pos.values] = valid_depth.values
    return depth_map

def run_robust_rdna_analysis(depth_map, bed_df, assembled_array_size):
    """Calculates coverage stats strictly across mapped rDNA variants (Unit level & Aggregated)."""
    all_variant_depths = []
    unit_stats = []
    
    for idx, row in bed_df.iterrows():
        s_idx = max(0, int(row['Array_Start']) - 1)
        e_idx = min(assembled_array_size, int(row['Array_End']))
        
        unit_depths = depth_map[s_idx:e_idx]
        
        # Calculate coverage using ONLY sequence > 0X to avoid skewing from assembly gaps
        valid_depths = unit_depths[unit_depths > 0]
        
        if len(valid_depths) > 0:
            all_variant_depths.append(valid_depths)
            unit_stats.append({
                'Variant_ID': row['Name'],
                'Start': int(row['Array_Start']),
                'End': int(row['Array_End']),
                'Mean_Depth': np.mean(valid_depths),
                'Median_Depth': np.median(valid_depths)
            })

    if not all_variant_depths:
        return None, pd.DataFrame()

    # Pool all variant base-pairs to get the true global sequence coverage stats
    concatenated_depths = np.concatenate(all_variant_depths)
    unit_df = pd.DataFrame(unit_stats)
    
    stats = {
        'assembled_copies': len(bed_df),
        'assembled_bp': sum(len(d) for d in all_variant_depths),
        'mean_depth': np.mean(concatenated_depths),
        'median_depth': np.median(concatenated_depths)
    }

    return stats, unit_df

def main():
    # ==============================================================================
    # CONFIGURATION & HARDCODED DATA
    # ==============================================================================
    targets_config = [
        {
            "target_id": "A3_X_CM090833.1", "strain": "A3", "contig": "CM090833.1",
            "hifi_depth": "./HiFi/A3_X_CM090833.1_raw_depth.tsv",  
            "bed_file": "../morphs/mappedmorphs/A3X.bed", "offset": 0 
        },
        {
            "target_id": "A4_X_CM090839.1", "strain": "A4", "contig": "CM090839.1",
            "hifi_depth": "./HiFi/A4_X_CM090839.1_raw_depth.tsv",
            "bed_file": "../morphs/mappedmorphs/A4X.bed", "offset": 0
        },
        {
            "target_id": "ISO1_X_CM090845.1", "strain": "ISO1", "contig": "CM090845.1",
            "hifi_depth": "./HiFi/ISO1_X_CM090845.1_raw_depth.tsv",
            "bed_file": "../morphs/mappedmorphs/ISO1X.bed", "offset": 0
        },
        {
            "target_id": "ISO1_Y_Y_Scaffold2", "strain": "ISO1", "contig": "Y_Scaffold2",
            "hifi_depth": "./HiFi/ISO1_Y_Y_Scaffold2_raw_depth.tsv",
            "bed_file": "../morphs/mappedmorphs/ISO1Y.bed", "offset": 0
        },
        {
            "target_id": "A3_Y_Y_scaffold2", "strain": "A3", "contig": "Y_scaffold2",
            "hifi_depth": "./HiFi/A3_Y_Y_scaffold2_raw_depth.tsv",
            "bed_file": "../morphs/mappedmorphs/A3Y.bed", "offset": 0
        },
        {
            "target_id": "A4_Y_Y_scaffold2", "strain": "A4", "contig": "Y_scaffold2",
            "hifi_depth": "./HiFi/A4_Y_Y_scaffold2_raw_depth.tsv",
            "bed_file": "../morphs/mappedmorphs/A4Y.bed", "offset": 0
        },
    ]

    # Hardcoded Assembly Gaps
    hardcoded_gaps = pd.DataFrame([
        {"Strain": "A3", "Contig": "Y_scaffold2", "Gap_Start": 923608, "Gap_End": 923707},
        {"Strain": "A3", "Contig": "Y_scaffold2", "Gap_Start": 1345899, "Gap_End": 1345998},
        {"Strain": "ISO1", "Contig": "Y_Contig33", "Gap_Start": 244658, "Gap_End": 246911},
        {"Strain": "ISO1", "Contig": "Y_Contig33", "Gap_Start": 247471, "Gap_End": 249802},
        {"Strain": "ISO1", "Contig": "Y_Contig33", "Gap_Start": 258301, "Gap_End": 260816},
        {"Strain": "ISO1", "Contig": "CM090845.1", "Gap_Start": 208215, "Gap_End": 209214},
        {"Strain": "A4", "Contig": "Y_scaffold2", "Gap_Start": 1093550, "Gap_End": 1093649}
    ])

    # Hardcoded Pp1Y2 Gene Coordinates
    hardcoded_genes = pd.DataFrame([
        {"Target": "ISO1_Y_Y_Scaffold2", "Start": 1239696, "End": 1240625, "Name": "Pp1Y2"},
        {"Target": "A3_Y_Y_scaffold2", "Start": 1337939, "End": 1338868, "Name": "Pp1Y2"},
        {"Target": "A4_Y_Y_scaffold2", "Start": 394106, "End": 395035, "Name": "Pp1Y2"}
    ])

    out_dir = "Stacked_rDNA_Figures"
    os.makedirs(out_dir, exist_ok=True)
    print(f"Outputs will be saved to: ./{out_dir}/")

    # --- OKABE-ITO SYSTEM COLORS ---
    COLOR_UNIQUE = '#009E73'      # Bluish Green for active unit sequence
    COLOR_INTERGENIC = '#CC79A7'  # Reddish Purple for intergenic sequence
    COLOR_GAP = '#999999'         # Grey for assembly gaps

    for config in targets_config:
        tid = config["target_id"]
        strain = config["strain"]
        contig = config["contig"]
        
        print(f"\n=========================================================")
        print(f"[{tid}] Starting rDNA Coverage Mapping...")

        if not os.path.exists(config["hifi_depth"]) or not os.path.exists(config["bed_file"]):
            print(f"[{tid}] Warning: Missing required files. Skipping...")
            continue

        hifi_df = pd.read_csv(config["hifi_depth"], sep='\t', header=None, names=['Contig', 'Pos', 'Depth'])
        
        assembled_array_size = hifi_df['Pos'].max()
        hifi_map = create_depth_map(hifi_df, assembled_array_size)
        
        bed_cols = ['Chr', 'Start', 'End', 'Name', 'Score', 'Strand']
        # Skip the 'track' header line and strictly grab only the first 6 columns
        bed_df = pd.read_csv(config["bed_file"], sep='\t', header=None, names=bed_cols, usecols=[0, 1, 2, 3, 4, 5], skiprows=1)
        
        # Filter for the correct contig to ensure architecture plot aligns with target
        bed_df = bed_df[bed_df['Chr'] == contig].copy()
        
        # Extract just the identifier (e.g., "Morph_8" from "Morph_8_(Gold)")
        bed_df['Name'] = bed_df['Name'].str.extract(r'(Morph_\d+)', expand=False).fillna(bed_df['Name'])
        # Rename "Morph" to "Variant" so it propagates to all figures and reports
        bed_df['Name'] = bed_df['Name'].str.replace('Morph', 'Variant')
        
        gap_df = hardcoded_gaps[(hardcoded_gaps['Strain'] == strain) & (hardcoded_gaps['Contig'] == contig)]
        gene_df = hardcoded_genes[hardcoded_genes['Target'] == tid]

        offset = config["offset"]
        array_end_coord = offset + assembled_array_size
        bed_df = bed_df[(bed_df['Start'] >= offset) & (bed_df['End'] <= array_end_coord)].copy()
        bed_df['Array_Start'] = bed_df['Start'] - offset + 1
        bed_df['Array_End'] = bed_df['End'] - offset + 1

        unit_mask = np.zeros(assembled_array_size + 1, dtype=bool)
        for _, row in bed_df.iterrows():
            s_idx = max(0, int(row['Array_Start']) - 1)
            e_idx = min(assembled_array_size, int(row['Array_End']))
            unit_mask[s_idx:e_idx] = True

        intergenic_mask = ~unit_mask
        ig_diff = np.diff(intergenic_mask.astype(int))
        ig_starts = np.where(ig_diff == 1)[0] + 1
        if intergenic_mask[0]: ig_starts = np.insert(ig_starts, 0, 0)
        ig_ends = np.where(ig_diff == -1)[0] + 1
        if intergenic_mask[-1]: ig_ends = np.append(ig_ends, len(intergenic_mask))
        intergenic_blocks = list(zip(ig_starts, ig_ends))

        print(f"[{tid}] Running array collapse statistics...")
        hifi_stats, hifi_unit_df = run_robust_rdna_analysis(hifi_map, bed_df, assembled_array_size)

        # ==============================================================================
        # WRITE EXPANDED SUMMARY REPORT
        # ==============================================================================
        report_file = os.path.join(out_dir, f"{tid}_collapse_report.txt")
        with open(report_file, 'w') as f:
            f.write(f"=== {tid} rDNA Array Collapse Report ===\n\n")
            
            if hifi_stats is not None and not hifi_unit_df.empty:
                f.write(f"--- HiFi Long-Read Platform ---\n")
                f.write(f"  Mapped Variant Copies: {hifi_stats['assembled_copies']}\n")
                f.write(f"  Valid Mapped Sequence: {hifi_stats['assembled_bp']:,.0f} bp\n")
                f.write(f"  Overall Mean rDNA Depth: {hifi_stats['mean_depth']:,.1f}X\n")
                f.write(f"  Overall Median rDNA Depth: {hifi_stats['median_depth']:,.1f}X\n\n")
                
                f.write(f"  [ Aggregated Variant Type Coverage ]\n")
                variant_summary = hifi_unit_df.groupby('Variant_ID').agg(
                    Count=('Variant_ID', 'count'),
                    Mean_of_Means=('Mean_Depth', 'mean'),
                    Mean_of_Medians=('Median_Depth', 'mean')
                ).reset_index().sort_values(by='Variant_ID')
                
                for _, m_row in variant_summary.iterrows():
                    f.write(f"    {m_row['Variant_ID']} (n={m_row['Count']} units): Mean = {m_row['Mean_of_Means']:,.1f}X | Median = {m_row['Mean_of_Medians']:,.1f}X\n")
                
                f.write(f"\n  [ Individual Variant Unit Coverage ]\n")
                for _, u_row in hifi_unit_df.iterrows():
                    f.write(f"    {u_row['Variant_ID']} (Coords: {u_row['Start']}-{u_row['End']}): Mean = {u_row['Mean_Depth']:,.1f}X | Median = {u_row['Median_Depth']:,.1f}X\n")
                f.write("\n")
                
        # High Res Binning
        bin_size_hifi = 500       
        
        def bin_data(df, unit_mask_array, current_bin_size):
            if df.empty: return pd.DataFrame()
            df_plot = df[df['Depth'] > 0].copy()
            df_plot['Window'] = (df_plot['Pos'] // current_bin_size) * current_bin_size
            binned = df_plot.groupby('Window')['Depth'].median().reset_index()
            binned['is_unit'] = [np.mean(unit_mask_array[w:w+current_bin_size]) > 0.50 for w in binned['Window']]
            return binned

        binned_hifi = bin_data(hifi_df, unit_mask, bin_size_hifi)

        # ==============================================================================
        # GENERATE SHARED VARIANT COLORMAP (OKABE-ITO)
        # ==============================================================================
        unique_variants = bed_df['Name'].unique() if not bed_df.empty else []
        
        # Sort variants naturally (e.g., Variant_1, Variant_2, Variant_10)
        unique_variants = sorted(unique_variants, key=lambda x: [int(c) if c.isdigit() else c.lower() for c in re.split(r'(\d+)', str(x))])
        
        # Cycle through Okabe-Ito colors for the variants
        variant_colors = {var: OKABE_ITO[i % len(OKABE_ITO)] for i, var in enumerate(unique_variants)}
        
        print(f"[{tid}] Variant color codes: " + ", ".join(f"{var}={variant_colors[var]}" for var in unique_variants))
        
        # ==============================================================================
        # 1. GENERATE SEPARATE VARIANT LEGEND
        # ==============================================================================
        if len(unique_variants) > 0:
            legend_fig, legend_ax = plt.subplots(figsize=(8, max(2, len(unique_variants) * 0.3)))
            legend_ax.axis('off')
            legend_patches = [patches.Patch(color=variant_colors[var], label=var) for var in unique_variants]
            
            ncols = 1 if len(unique_variants) <= 10 else 2
            
            legend_ax.legend(handles=legend_patches, loc='center', title=f"{tid} Variant Colors", 
                             fontsize=12, title_fontsize=14, ncol=ncols, frameon=True)
            
            legend_out_png = os.path.join(out_dir, f"{tid}_variant_legend.png")
            legend_fig.savefig(legend_out_png, dpi=300, bbox_inches='tight')
            plt.close(legend_fig)
            print(f"[{tid}] Saved separate variant legend to {legend_out_png}")

        # ==============================================================================
        # 2. GENERATE COVERAGE LANDSCAPE PLOT
        # ==============================================================================
        fig, ax = plt.subplots(figsize=(24, 6))
        fig.suptitle(f"rDNA Scaffold Coverage Landscape", fontsize=18, fontweight='bold', y=0.95)

        x_hf, y_hf = binned_hifi['Window'] / 1e6, binned_hifi['Depth']
        hifi_colors = np.where(binned_hifi['is_unit'], COLOR_UNIQUE, COLOR_INTERGENIC)

        ax.plot(x_hf, y_hf, color='dimgrey', linewidth=0.5, zorder=1, alpha=0.6)
        ax.scatter(x_hf, y_hf, c=hifi_colors, s=8, zorder=2, alpha=0.9, edgecolors='none')
        
        if hifi_stats:
            # Okabe-Ito Vermillion and Blue
            ax.axhline(hifi_stats['mean_depth'], color='#D55E00', linestyle='--', linewidth=2.5, alpha=0.9, label=f"Mean Depth ({hifi_stats['mean_depth']:,.1f}X)", zorder=5)
            ax.axhline(hifi_stats['median_depth'], color='#0072B2', linestyle='-.', linewidth=2.5, alpha=0.9, label=f"Median Depth ({hifi_stats['median_depth']:,.1f}X)", zorder=5)
        
        ax.set_ylabel("HiFi Depth (X)", fontsize=14, fontweight='bold')
        ax.set_yscale('symlog', linthresh=50)
        max_hifi = y_hf.max()
        ax.set_ylim(0, max_hifi * 1.5 if (pd.notna(max_hifi) and max_hifi > 50) else 150)
        ax.set_xlabel("Scaffold Position (Mbp)", fontsize=14, fontweight='bold')

        # SHADE ASSEMBLY GAPS
        if not gap_df.empty:
            for _, row in gap_df.iterrows():
                start_mb, end_mb = row['Gap_Start'] / 1e6, row['Gap_End'] / 1e6
                min_vis_width = 0.01 
                if (end_mb - start_mb) < min_vis_width:
                    midpoint = (start_mb + end_mb) / 2.0
                    start_mb = midpoint - (min_vis_width / 2.0)
                    end_mb = midpoint + (min_vis_width / 2.0)
                # Okabe-Ito Grey for gaps
                ax.axvspan(start_mb, end_mb, color=COLOR_GAP, alpha=0.4, lw=0, zorder=0)

        # X-AXIS ANNOTATION: INTERGENIC SEQUENCE
        for start, end in intergenic_blocks:
            ax.add_patch(patches.Rectangle((start / 1e6, 0), (end - start) / 1e6, 0.04, 
                                           color=COLOR_INTERGENIC, transform=ax.get_xaxis_transform(), 
                                           zorder=5))

        # X-AXIS ANNOTATION: RDNA VARIANTS
        for _, row in bed_df.iterrows():
            variant_name = row['Name']
            color = variant_colors[variant_name]
            unit_width = (row['Array_End'] - row['Array_Start']) / 1e6
            rect = patches.Rectangle((row['Array_Start'] / 1e6, 0), unit_width, 0.04, 
                                     color=color, transform=ax.get_xaxis_transform(), 
                                     zorder=10)
            ax.add_patch(rect)
            
        # GENE ANNOTATION (Pp1Y2)
        if not gene_df.empty:
            for _, row in gene_df.iterrows():
                g_mid = ((row['Start'] + row['End']) / 2.0) / 1e6
                y_max = ax.get_ylim()[1]
                ax.axvline(g_mid, color='black', linestyle='-', linewidth=2, zorder=12)
                ax.plot(g_mid, y_max * 0.88, marker='v', color='black', markersize=10, zorder=15)
                bbox_props = dict(boxstyle="round,pad=0.3", fc="white", ec="black", lw=1.5, alpha=0.9)
                ax.text(g_mid, y_max * 0.92, row['Name'], ha='center', va='bottom', 
                        fontsize=12, fontweight='bold', color='black', bbox=bbox_props, zorder=20)

        # Clean up the legend
        ax.plot([], [], marker='o', color=COLOR_UNIQUE, linestyle='None', markersize=6, label='Mapped Variant')
        ax.plot([], [], marker='o', color=COLOR_INTERGENIC, linestyle='None', markersize=6, label='Intergenic Sequence')
        if not gap_df.empty:
            ax.fill_between([], [], color=COLOR_GAP, alpha=0.4, label='Assembly Gap')
        if not gene_df.empty:
            ax.plot([], [], color='black', linestyle='-', linewidth=2, marker='v', markersize=8, label='Pp1Y2 Gene')
        
        handles, labels = ax.get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        ax.legend(by_label.values(), by_label.keys(), loc='upper left', bbox_to_anchor=(1.02, 1), framealpha=0.9, fontsize=11)

        # X-AXIS LIMITS & INVERSION
        left_limit = 0
        right_limit = assembled_array_size / 1e6

        if "_X_" in tid and not bed_df.empty:
            array_start_mb = bed_df['Array_Start'].min() / 1e6
            array_end_mb = bed_df['Array_End'].max() / 1e6
            buffer_mb = max((array_end_mb - array_start_mb) * 0.05, 0.05)
            left_limit = max(0, array_start_mb - buffer_mb)
            right_limit = array_end_mb + buffer_mb

        ax.set_xlim(right_limit, left_limit)

        plot_out_png = os.path.join(out_dir, f"{tid}_stacked_coverage.png")
        plot_out_svg = os.path.join(out_dir, f"{tid}_stacked_coverage.svg")
        plt.savefig(plot_out_png, dpi=300, bbox_inches='tight')
        plt.savefig(plot_out_svg, format='svg', bbox_inches='tight')
        print(f"[{tid}] Saved coverage plots to {plot_out_png} and .svg")
        plt.close(fig)

        # ==============================================================================
        # 3. GENERATE VARIANT ARCHITECTURE PLOT
        # ==============================================================================
        if not bed_df.empty:
            arch_fig, arch_ax = plt.subplots(figsize=(20, 3))
            arch_fig.suptitle(f"{tid} rDNA Variant Architecture", fontsize=18, fontweight='bold', y=1.05)

            orig_start_min = bed_df['Start'].min()
            orig_end_max = bed_df['End'].max()
            local_relative_end = (orig_end_max + 5000) - orig_start_min

            # Draw Backbone Line
            arch_ax.plot([-5000, local_relative_end], [0.5, 0.5], color='black', lw=1, zorder=1)
            
            # Plot Variants perfectly matched to Okabe-Ito colors
            for _, row in bed_df.iterrows():
                color = variant_colors[row['Name']]
                rel_start = row['Start'] - orig_start_min
                width = row['End'] - row['Start']
                
                arch_ax.add_patch(patches.Rectangle((rel_start, 0.5 - 0.15), width, 0.3, 
                                                     color=color, lw=0.5, ec='black', zorder=2))

            arch_ax.set_ylim(0.1, 0.9)
            arch_ax.set_yticks([])
            arch_ax.set_ylabel(strain, fontsize=14, fontweight='bold', rotation=0, labelpad=30, va='center')
            arch_ax.spines['top'].set_visible(False)
            arch_ax.spines['right'].set_visible(False)
            arch_ax.spines['left'].set_visible(False)
            
            # Invert X-axis
            arch_ax.set_xlim(local_relative_end, -5000)
            
            def kb_formatter(x, pos):
                if local_relative_end > 1_000_000:
                    return f"{x / 1e6:.2f} Mb"
                return f"{x / 1000:.0f} kb"
                
            arch_ax.xaxis.set_major_formatter(plt.FuncFormatter(kb_formatter))
            arch_ax.tick_params(axis='x', labelsize=11, length=5)
            arch_ax.set_xlabel("Relative Array Coordinates", fontsize=13, fontweight='bold', labelpad=10)

            # Strain-Specific Individual Legend using the synchronized colors
            legend_elements = [patches.Rectangle((0,0),1,1, facecolor=variant_colors[var], edgecolor='black', label=var) for var in unique_variants]
            arch_ax.legend(handles=legend_elements, loc='center left', bbox_to_anchor=(1.01, 0.5), 
                           title=f"{strain} Variants", title_fontproperties={'weight':'bold'}, 
                           fontsize=10, frameon=False)

            arch_out_png = os.path.join(out_dir, f"{tid}_architecture.png")
            arch_out_svg = os.path.join(out_dir, f"{tid}_architecture.svg")
            plt.savefig(arch_out_png, dpi=300, bbox_inches='tight')
            plt.savefig(arch_out_svg, format='svg', bbox_inches='tight')
            plt.close(arch_fig)
            print(f"[{tid}] Saved architecture plots to {arch_out_png} and .svg")
        else:
            print(f"[{tid}] No bed entries found. Skipping architecture plot.")

if __name__ == "__main__":
    main()
