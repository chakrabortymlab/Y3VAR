import os
from collections import defaultdict 

# ==============================================================================
# RIBOTIN-AWARE ARRAY SEGMENTER (ALL-IN-ONE PIPELINE)
# ==============================================================================

# Ensure these match the exact contig names in your RepeatMasker file!
ALLOWED_CONTIGS = set(["Y_scaffold2"])

def parse_and_generate_raw_units(rm_bed):
    """
    Scans the RepeatMasker baseline and extracts phase-matched rDNA units.
    Returns a list of raw unit dictionaries.
    """
    features = defaultdict(list)
    raw_units = []

    if not os.path.exists(rm_bed):
        print(f">>> ERROR: {rm_bed} not found!")
        return raw_units

    # 1. Parse RepeatMasker BED
    with open(rm_bed) as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 4: continue
            chrom = parts[0]
            if chrom not in ALLOWED_CONTIGS: continue
            s, e, name = int(parts[1]), int(parts[2]), parts[3]
            features[chrom].append((s, e, name))
            
    for chrom, feats in features.items():
        feats.sort(key=lambda x: x[0])
        
        # Helper: Merge fragmented RM annotations of the same gene (within 3kb)
        def merge_feature(target_string, max_gap=3000):
            merged = []
            for s, e, name in feats:
                if target_string in name:
                    if not merged:
                        merged.append([s, e])
                    else:
                        if s - merged[-1][1] <= max_gap:
                            merged[-1][1] = max(merged[-1][1], e)
                        else:
                            merged.append([s, e])
            return merged

        m_18s = merge_feature('18S')
        m_28s = merge_feature('28S')
        m_ets = merge_feature('ETS')

        # Helper: Find the closest feature in a specific direction
        def get_closest(target_pos, candidates, direction='downstream'):
            valid = []
            for cs, ce in candidates:
                if direction == 'downstream' and cs > target_pos:
                    valid.append((cs, ce))
                elif direction == 'upstream' and ce < target_pos:
                    valid.append((cs, ce))
            if not valid: return None
            if direction == 'downstream': return min(valid, key=lambda x: x[0])
            if direction == 'upstream': return max(valid, key=lambda x: x[1])

        # Build Phase-Matched Units anchored on 18S
        for s_18, e_18 in m_18s:
            closest_28_down = get_closest(e_18, m_28s, 'downstream')
            closest_28_up = get_closest(s_18, m_28s, 'upstream')

            dist_down = (closest_28_down[0] - e_18) if closest_28_down else float('inf')
            dist_up = (s_18 - closest_28_up[1]) if closest_28_up else float('inf')

            if dist_down == float('inf') and dist_up == float('inf'): continue 

            strand = '+' if dist_down < dist_up else '-'

            if strand == '+':
                term_ets = get_closest(closest_28_down[1], m_ets, 'downstream')
                if not term_ets: continue
                unit_start, unit_end = s_18, term_ets[1]
            else:
                term_ets = get_closest(closest_28_up[0], m_ets, 'upstream')
                if not term_ets: continue
                unit_start, unit_end = term_ets[0], e_18

            # Scan for TEs perfectly inside this specific unit's boundaries
            inserted_tes = set()
            for fs, fe, fname in feats:
                if fs >= unit_start and fe <= unit_end:
                    base_name = fname.split('#')[0]
                    if base_name in ['Doc', 'I_DM']: inserted_tes.add(base_name)
                    elif 'R1' in base_name: inserted_tes.add('R1')
                    elif 'R2' in base_name: inserted_tes.add('R2')

            if not inserted_tes:
                label = "rDNA_Unit_Canonical"
            else:
                label = "rDNA_Unit_" + "_".join(sorted(list(inserted_tes))) + "_Inserted"

            raw_units.append({
                'chrom': chrom, 
                'start': unit_start, 
                'end': unit_end, 
                'name': label, 
                'score': '0', 
                'strand': strand, 
                'len': unit_end - unit_start
            })

    return raw_units

def clean_and_write_bed(raw_units, out_bed):
    """
    Filters out stunted units and overlapping paradoxes, then writes the final BED6.
    """
    print(f"Total raw units generated: {len(raw_units)}")
    
    # 1. Strict Length Filter (9kb to 30kb)
    length_filtered = [u for u in raw_units if 9000 <= u['len'] <= 30000]
    print(f"Units remaining after length filter: {len(length_filtered)}")

    # Sort strictly by chromosome, then start coordinate
    length_filtered.sort(key=lambda x: (x['chrom'], x['start']))

    # 2. Strict Overlap Filter (Fixes branching and Möbius strip collisions)
    clean_units = []
    for u in length_filtered:
        if not clean_units:
            clean_units.append(u)
        else:
            prev = clean_units[-1]
            
            # Only check overlap if they are on the same chromosome
            if u['chrom'] == prev['chrom']:
                overlap = max(0, min(u['end'], prev['end']) - max(u['start'], prev['start']))
                
                # If they overlap by more than 1000bp, keep the more complete (longer) one
                if overlap > 1000:
                    if u['len'] > prev['len']:
                        clean_units[-1] = u
                else:
                    clean_units.append(u)
            else:
                clean_units.append(u)

    print(f"Final pristine units saved to BED: {len(clean_units)}\n")

    # 3. Write standard BED6
    with open(out_bed, 'w') as f:
        for u in clean_units:
            f.write(f"{u['chrom']}\t{u['start']}\t{u['end']}\t{u['name']}\t{u['score']}\t{u['strand']}\n")

# ==============================================================================
# EXECUTION
# ==============================================================================
if __name__ == "__main__":
    input_rm_file = "A3RM.bed"                   # must be BED file
    output_clean_bed = "PhaseMatched_CLEAN.bed"    # The final Ribotin-ready BED

    raw = parse_and_generate_raw_units(input_rm_file)
    if raw:
        clean_and_write_bed(raw, output_clean_bed)
