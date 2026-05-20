import matplotlib.pyplot as plt
import matplotlib.collections as mc
import numpy as np
import os
from statistics import median

# ==========================================
# MAXIMUM CONTRAST STRAND COLORS
# ==========================================
FWD_COLOR = "#004488"  # Dark Cobalt Blue
REV_COLOR = "#EE6677"  # Bright Rose-Red

def parse_paf_data(paf_file):
    """Parses PAF and records reference scaffold order."""
    alignments = []
    q_lens, t_lens = {}, {}
    t_order_in_file = []
    with open(paf_file, 'r') as f:
        for line in f:
            cols = line.strip().split('\t')
            if len(cols) < 12: continue
            q_name, q_len, t_name, t_len = cols[0], int(cols[1]), cols[5], int(cols[6])
            q_lens[q_name], t_lens[t_name] = q_len, t_len
            if t_name not in t_order_in_file: 
                t_order_in_file.append(t_name)
            alignments.append({
                'q': q_name, 'qs': int(cols[2]), 'qe': int(cols[3]), 'strand': cols[4],
                't': t_name, 'ts': int(cols[7]), 'te': int(cols[8]), 
                'alen': int(cols[10])
            })
    return alignments, q_lens, t_lens, t_order_in_file

def get_ref_ordered_query(alns, t_order, q_lens):
    """Orders Query scaffolds by their median position on the Reference."""
    q_to_t_stats = {}
    for a in alns:
        q, t, l = a['q'], a['t'], a['alen']
        q_to_t_stats.setdefault(q, {}).setdefault(t, 0)
        q_to_t_stats[q][t] += l
    
    q_best_t = {q: max(targets, key=targets.get) for q, targets in q_to_t_stats.items()}
    t_name_to_idx = {name: i for i, name in enumerate(t_order)}
    
    q_sort_vals = {}
    for q in q_lens:
        if q not in q_best_t:
            q_sort_vals[q] = (len(t_order), 0)
            continue
        best_t = q_best_t[q]
        starts = [a['ts'] for a in alns if a['q'] == q and a['t'] == best_t]
        q_sort_vals[q] = (t_name_to_idx[best_t], median(starts))
    return sorted(q_lens.keys(), key=lambda x: q_sort_vals[x])

def create_publication_dotplot(paf_file, y_label, output_base):
    print(f"Generating Publication Plot: {output_base}")
    alns, q_lens, t_lens, t_order = parse_paf_data(paf_file)
    q_order = get_ref_ordered_query(alns, t_order, q_lens)
    
    t_offsets, curr_t = {}, 0
    for name in t_order: t_offsets[name] = curr_t; curr_t += t_lens[name]
    q_offsets, curr_q = {}, 0
    for name in q_order: q_offsets[name] = curr_q; curr_q += q_lens[name]
    
    fig, ax = plt.subplots(figsize=(24, 24))
    
    # 1. GRID (Fine and light in the background)
    for name in t_order: ax.axvline(t_offsets[name], color='#333333', linewidth=0.4, alpha=0.4)
    for name in q_order: ax.axhline(q_offsets[name], color='#333333', linewidth=0.4, alpha=0.4)
        
    # 2. ALIGNMENTS (Binary high-contrast coloring)
    lines, colors = [], []
    for a in alns:
        if a['alen'] < 1: continue 
        tx1, tx2 = t_offsets[a['t']] + a['ts'], t_offsets[a['t']] + a['te']
        qy1, qy2 = q_offsets[a['q']] + a['qs'], q_offsets[a['q']] + a['qe']
        
        lines.append([(tx1, qy1), (tx2, qy2)] if a['strand'] == '+' else [(tx1, qy2), (tx2, qy1)])
        colors.append(FWD_COLOR if a['strand'] == '+' else REV_COLOR)

    lc = mc.LineCollection(lines, colors=colors, linewidths=9, alpha=1.0)
    ax.add_collection(lc)
    
    # 3. AXIS FORMATTING (Hide tick labels)
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.set_xticks([])
    ax.set_yticks([])
    
    ax.set_xlabel("ISO-1 (This study)", fontweight='bold', fontsize=60, labelpad=25)
    ax.set_ylabel(y_label, fontweight='bold', fontsize=60, labelpad=25)
    
    # 4. LEGEND
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], color=FWD_COLOR, lw=8, label='Forward Orientation (+)'),
        Line2D([0], [0], color=REV_COLOR, lw=8, label='Reverse Orientation (-)')
    ]
    ax.legend(handles=legend_elements, loc='upper left', fontsize=44, 
              frameon=True, framealpha=1.0, edgecolor='black', borderpad=1)

    ax.set_xlim(0, curr_t); ax.set_ylim(0, curr_q)
    plt.tight_layout()
    plt.savefig(f"{output_base}.png", dpi=600)
    plt.savefig(f"{output_base}.svg")
    plt.close()

if __name__ == "__main__":
    tasks = [
        ("map_dmelA3YDelgado_to_dmelISO1YDelgado.paf", "A3", "A3_publication"),
        ("map_dmelA4YDelgado_to_dmelISO1YDelgado.paf", "A4", "A4_publication"),
        ("map_YChang_to_dmelISO1YDelgado.paf", "Chang et al (ISO-1)", "YChang_publication")
    ]
    for paf, y_label, out in tasks:
        if os.path.exists(paf): 
            create_publication_dotplot(paf, y_label, out)
        else:
            print(f"File not found: {paf}")
