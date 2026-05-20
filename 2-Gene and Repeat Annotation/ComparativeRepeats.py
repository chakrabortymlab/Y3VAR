import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import argparse
import sys
import numpy as np

#Usage: python ComparativeRepeats.py --iso1_te ISO-1tes.fasta.out --iso1_sat ISO-1sats.fasta.out --a3_te A3tes.fasta.out --a3_sat A3sats.fasta.out --a4_te A4tes.fasta.out --a4_sat A4sats.fasta.out --output Comparative_Y_Extensive_Data.xlsx

# =============================================================================
# COLOR PALETTES (Okabe-Ito)
# =============================================================================
# Standard categorical Okabe-Ito palette for discrete data (Figure 2)
OKABE_ITO = ["#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000"]

# =============================================================================
# USER DEFINED LISTS
# =============================================================================

# 1. Simple Repeats for Cytological/Array Analysis
SIMPLE_REPEAT_TARGETS = [
    "AAGAG", "AAGAGAG", "AAGAC", "AATAT", "AATAG", 
    "AATAC", "AATAAAC", "AATAGAC"
]

# 2. Specific Satellites for Figure 2 (Includes Simple + Complex)
TARGET_SATELLITES = [
    "HETRP_DM", "15mer_SAT", "353bp_SAT", "356bp_SAT", "260bp_SAT", 
    "359bp_SAT", "360bp_SAT", "90U_SAT", "AAGAG", "AAGAGAG", 
    "AAGAC", "AATAT", "AATAG", "AATAC", "AATAAAC", "AATAGAC"
]

# 3. Full Library List for Detailed Excel Reporting
LIBRARY_ELEMENTS = {
    # --- SATELLITES ---
    "HETRP_DM": "Satellite/Satellite", "15mer_SAT": "Satellite/Satellite",
    "353bp_SAT": "Satellite/1pt688", "356bp_SAT": "Satellite/1pt688",
    "260bp_SAT": "Satellite/1pt688", "359bp_SAT": "Satellite/1pt688",
    "360bp_SAT": "Satellite/Satellite", "90U_SAT": "Satellite/Satellite",
    "AAGAG": "Satellite", "AAGAGAG": "Satellite", "AAGAC": "Satellite",
    "AATAT": "Satellite", "AATAG": "Satellite", "AATAC": "Satellite",
    "AATAAAC": "Satellite", "AATAGAC": "Satellite",
    # --- TEs ---
    "MARINA": "DNA/TcMar-Mariner", "NOF_FB": "DNA/MULE-NOF", "PLACW_DM": "DNA/P",
    "PROTOP": "DNA/P", "PROTOP_A": "DNA/P", "PROTOP_B": "DNA/P",
    "LOOPER1_DM": "DNA/PiggyBac", "POGO": "DNA/TcMar-Pogo", "POGON1": "DNA/TcMar-Pogo",
    "BARI1": "DNA/TcMar-Tc1", "BARI_DM": "DNA/TcMar-Tc1", "FB4_DM": "DNA/TcMar-Tc1",
    "MINOS": "DNA/TcMar-Tc1", "Mariner2_DM": "DNA/TcMar-Tc1", "PARIS": "DNA/TcMar-Tc1",
    "S2_DM": "DNA/TcMar-Tc1", "S_DM": "DNA/TcMar-Tc1", "TC1-2_DM": "DNA/TcMar-Tc1",
    "TC1_DM": "DNA/TcMar-Tc1", "UHU": "DNA/TcMar-Tc1", "M4DM": "DNA/CMC-Transib",
    "TRANSIB1": "DNA/CMC-Transib", "TRANSIB2": "DNA/CMC-Transib", "TRANSIB3": "DNA/CMC-Transib",
    "TRANSIB4": "DNA/CMC-Transib", "Transib-N1_DM": "DNA/CMC-Transib", "HOBO": "DNA/hAT-hobo",
    "DMCR1A": "Non-LTR_retrotransposon/CR1", "IVK_DM": "Non-LTR_retrotransposon/I",
    "I_DM": "Non-LTR_retrotransposon/I", "BS": "Non-LTR_retrotransposon/Jockey",
    "BS2": "Non-LTR_retrotransposon/Jockey", "BS3_DM": "Non-LTR_retrotransposon/Jockey",
    "BS4_DM": "Non-LTR_retrotransposon/Jockey", "DOC": "Non-LTR_retrotransposon/Jockey",
    "DOC2_DM": "Non-LTR_retrotransposon/Jockey", "DOC3_DM": "Non-LTR_retrotransposon/Jockey",
    "DOC4_DM": "Non-LTR_retrotransposon/Jockey", "DOC5_DM": "Non-LTR_retrotransposon/Jockey",
    "DOC6_DM": "Non-LTR_retrotransposon/Jockey", "FW2_DM": "Non-LTR_retrotransposon/Jockey",
    "FW3_DM": "Non-LTR_retrotransposon/Jockey", "FW_DM": "Non-LTR_retrotransposon/Jockey",
    "G3_DM": "Non-LTR_retrotransposon/Jockey", "G4_DM": "Non-LTR_retrotransposon/Jockey",
    "G5A_DM": "Non-LTR_retrotransposon/Jockey", "G5_DM": "Non-LTR_retrotransposon/Jockey",
    "G6_DM": "Non-LTR_retrotransposon/Jockey", "G7_DM": "Non-LTR_retrotransposon/Jockey",
    "G_DM": "Non-LTR_retrotransposon/Jockey", "HELENA": "Non-LTR_retrotransposon/Jockey",
    "HELENA_RT": "Non-LTR_retrotransposon/Jockey", "Jockey2": "Non-LTR_retrotransposon/Jockey",
    "LINEJ1_DM": "Non-LTR_retrotransposon/Jockey", "Baggins1": "Non-LTR_retrotransposon/LOA",
    "Bilbo": "Non-LTR_retrotransposon/LOA", "TRIM": "Non-LTR_retrotransposon/LOA",
    "DNAREP1_DM": "RC/Helitron", "ISFUN1": "RC/Helitron", "PENELOPE": "Non-LTR_retrotransposon/Penelope",
    "DMRT1A": "Non-LTR_retrotransposon/R1", "DMRT1B": "Non-LTR_retrotransposon/R1",
    "DMRT1C": "Non-LTR_retrotransposon/R1", "R1-2_DM": "Non-LTR_retrotransposon/R1",
    "R1_DM": "Non-LTR_retrotransposon/R1", "R2B_DM": "Non-LTR_retrotransposon/R2",
    "R2_DM": "Non-LTR_retrotransposon/R2", "AT_rich": "Low_complexity/Low_complexity",
    "GC_rich": "Low_complexity/Low_complexity", "Helitron1_DM": "RC/Helitron",
    "DMRP1": "Other/Other", "DMRPR": "Other/Other", "Transib5": "DNA/CMC-Transib",
    "MARWOLEN1": "DNA/TcMar-Tc1", "Galileo_DB": "DNA/P", "Galileo_DM": "DNA/P",
    "LmeSINE1c": "SINE/tRNA-Deu-L2", "Hoyak3": "DNA/hAT-Pegasus", "Hoana8": "DNA/hAT-hobo",
    "Hoana3": "DNA/hAT-hobo", "Hoana5": "DNA/hAT-Ac", "Hoana7": "DNA/hAT-Ac",
    "Hoyak1": "DNA/hAT-hobo", "Hoana2": "DNA/hAT-Ac", "Hoyak2": "DNA/hAT-Ac",
    "Hoana4": "DNA/hAT-hobo", "Hoana1": "DNA/hAT-hobo", "Hoana6": "DNA/hAT-Ac",
    "Hopseu2": "DNA/hAT-Pegasus", "Homo1": "DNA/hAT-hobo", "Howilli1": "DNA/hAT-Ac",
    "Homo11": "DNA/hAT-hobo", "Hovi1": "DNA/hAT-hobo", "Homo4": "DNA/hAT-hobo",
    "Howilli2": "DNA/hAT-hobo", "Homo3": "DNA/hAT-Ac", "Howilli4": "DNA/hAT-hobo",
    "Homo6": "DNA/hAT-Pegasus", "Hovi2": "DNA/hAT-hobo", "Homo8": "DNA/hAT-Ac",
    "Hogri1": "DNA/hAT-hobo", "Howilli3": "DNA/hAT-Ac", "Hopers3": "DNA/hAT-Pegasus",
    "Homo5": "DNA/hAT-Ac", "Homo9": "DNA/hAT-Ac", "Howilli5": "DNA/hAT-hobo",
    "Homo10": "DNA/hAT-Ac", "Homo7": "DNA/hAT-Ac", "Hopers2": "DNA/hAT-hobo",
    "Homo2": "DNA/hAT-Ac", "LOA": "Non-LTR_retrotransposon/LOA",
    "dmel_rDNA_18S_rDNA": "rDNA/rDNA", "dmel_rDNA_28S_rDNA": "rDNA/rDNA",
    "dmel_rDNA_5_8S_rDNA": "rDNA/rDNA", "dmel_rDNA_ETS": "rDNA/rDNA",
    "dmel_rDNA_ITS1": "rDNA/rDNA", "dmel_rDNA_ITS2": "rDNA/rDNA",
    "Jockey-3_Dmel_08212020": "Non-LTR_retrotransposon/Jockey", "NTS_Uniq_Dmel": "Other/IGS",
    "NTS_95_Dmel": "Other/IGS", "NTS_330_Dmel": "Other/IGS", "NTS_240_Dmel": "Other/IGS",
    "HETA": "Non-LTR_retrotransposon/Jockey", "Heta-1_D": "Non-LTR_retrotransposon/Jockey",
    "Heta-2": "Non-LTR_retrotransposon/Jockey", "Heta-3": "Non-LTR_retrotransposon/Jockey",
    "Heta-5": "Non-LTR_retrotransposon/Jockey", "TART-A": "Non-LTR_retrotransposon/Jockey",
    "TAHRE": "Non-LTR_retrotransposon/Jockey", "TART_B1": "Non-LTR_retrotransposon/Jockey",
    "TART-C_NTPR": "Non-LTR_retrotransposon/Jockey", "Jockey-1_DSi": "Non-LTR_retrotransposon/Jockey",
    "Jockey-4_DSim": "Non-LTR_retrotransposon/Jockey", "Jockey-5_DYa": "Non-LTR_retrotransposon/Jockey",
    "Helena_Ds": "Non-LTR_retrotransposon/Jockey", "DMR_DV": "Other/Other",
    "Polinton-1_DY": "DNA/Maverick", "hAT-1_DSi": "DNA/hAT-Ac", "P-1_DY": "DNA/P",
    "DNA4-1_DK": "RC/Helitron", "Jockey-1_DEr": "Non-LTR_retrotransposon/Jockey",
    "Transib-1_DRh": "DNA/CMC-Transib", "Jockey-3_DGri": "Non-LTR_retrotransposon/Jockey",
    "Jockey-6_DK": "Non-LTR_retrotransposon/Jockey", "Helitron-1_DEu": "RC/Helitron",
    "Helitron-N1_DBi": "RC/Helitron", "TART_DV": "Non-LTR_retrotransposon/Jockey",
    "hAT-2_DY": "DNA/hAT-hobo", "dmel_rDNA_2s_rDNA": "Other/rDNA", "CIRCE": "LTR/Gypsy",
    "Copia2_I": "LTR/Copia", "Copia2_LTR_DM": "LTR/Copia", "Copia_I": "LTR/Copia",
    "Copia_LTR": "LTR/Copia", "DM1731_I": "LTR/Copia", "DM1731_LTR": "LTR/Copia",
    "ACCORD2_I": "LTR/Gypsy", "ACCORD2_LTR": "LTR/Gypsy", "ACCORD_I": "LTR/Gypsy",
    "ACCORD_LTR": "LTR/Gypsy", "BLASTOPIA_I": "LTR/Gypsy", "BLASTOPIA_LTR": "LTR/Gypsy",
    "BLOOD_I": "LTR/Gypsy", "BLOOD_LTR": "LTR/Gypsy", "BURDOCK_I": "LTR/Gypsy",
    "BURDOCK_LTR": "LTR/Gypsy", "DM176_I": "LTR/Gypsy", "DM176_LTR": "LTR/Gypsy",
    "DM297_I": "LTR/Gypsy", "DM297_LTR": "LTR/Gypsy", "DM412": "LTR/Gypsy",
    "DM412B_LTR": "LTR/Gypsy", "DMLTR5": "LTR/Gypsy", "DMTOM1_LTR": "LTR/Gypsy",
    "GTWIN_I": "LTR/Gypsy", "GTWIN_LTR": "LTR/Gypsy", "Gypsy2_I": "LTR/Gypsy",
    "Gypsy2_LTR": "LTR/Gypsy", "Gypsy3_I": "LTR/Gypsy", "Gypsy3_LTR": "LTR/Gypsy",
    "Gypsy4_I": "LTR/Gypsy", "Gypsy4_LTR": "LTR/Gypsy", "Gypsy5_I": "LTR/Gypsy",
    "Gypsy5_LTR": "LTR/Gypsy", "Gypsy6A_LTR": "LTR/Gypsy", "Gypsy6_I": "LTR/Gypsy",
    "Gypsy6_LTR": "LTR/Gypsy", "Gypsy7_I": "LTR/Gypsy", "Gypsy7_LTR": "LTR/Gypsy",
    "Gypsy8_I": "LTR/Gypsy", "Gypsy8_LTR": "LTR/Gypsy", "Gypsy9_I": "LTR/Gypsy",
    "Gypsy9_LTR": "LTR/Gypsy", "Gypsy10_I": "LTR/Gypsy", "Gypsy10_LTR": "LTR/Gypsy",
    "Gypsy11_I": "LTR/Gypsy", "Gypsy11_LTR": "LTR/Gypsy", "Gypsy12A_LTR": "LTR/Gypsy",
    "Gypsy12_I": "LTR/Gypsy", "Gypsy12_LTR": "LTR/Gypsy", "Gypsy_I": "LTR/Gypsy",
    "Gypsy_LTR": "LTR/Gypsy", "HMSBEAGLE_I": "LTR/Gypsy", "IDEFIX_I": "LTR/Gypsy",
    "IDEFIX_LTR": "LTR/Gypsy", "Invader1_I": "LTR/Gypsy", "Invader1_LTR": "LTR/Gypsy",
    "Invader2_I": "LTR/Gypsy", "Invader2_LTR": "LTR/Gypsy", "Invader3_I": "LTR/Gypsy",
    "Invader3_LTR": "LTR/Gypsy", "Invader4_I": "LTR/Gypsy", "Invader4_LTR": "LTR/Gypsy",
    "Invader5_I": "LTR/Gypsy", "Invader5_LTR": "LTR/Gypsy", "Invader6_I": "LTR/Gypsy",
    "Invader6_LTR": "LTR/Gypsy", "MDG1_I": "LTR/Gypsy", "MDG1_LTR": "LTR/Gypsy",
    "MDG3_I": "LTR/Gypsy", "MDG3_LTR": "LTR/Gypsy", "MICROPIA_I": "LTR/Gypsy",
    "MICROPIA_LTR": "LTR/Gypsy", "NOMAD_I": "LTR/Gypsy", "NOMAD_LTR": "LTR/Gypsy",
    "OSVALDO_I": "LTR/Gypsy", "OSVALDO_LTR": "LTR/Gypsy", "QUASIMODO_I": "LTR/Gypsy",
    "QUASIMODO_LTR": "LTR/Gypsy", "ROVER-I_DM": "LTR/Gypsy", "ROVER-LTR_DM": "LTR/Gypsy",
    "STALKER4_I": "LTR/Gypsy", "STALKER4_LTR": "LTR/Gypsy", "Stalker2_I": "LTR/Gypsy",
    "Stalker2_LTR": "LTR/Gypsy", "Stalker3_LTR": "LTR/Gypsy", "TABOR_I": "LTR/Gypsy",
    "TABOR_LTR": "LTR/Gypsy", "TIRANT_I": "LTR/Gypsy", "TIRANT_LTR": "LTR/Gypsy",
    "TOM_I": "LTR/Gypsy", "TOM_LTR": "LTR/Gypsy", "TRANSPAC_I": "LTR/Gypsy",
    "TRANSPAC_LTR": "LTR/Gypsy", "TV1I": "LTR/Gypsy", "TV1LTR": "LTR/Gypsy",
    "Ulysses_I": "LTR/Gypsy", "Ulysses_LTR": "LTR/Gypsy", "ZAM_I": "LTR/Gypsy",
    "ZAM_LTR": "LTR/Gypsy", "BATUMI_I": "LTR/Pao", "BATUMI_LTR": "LTR/Pao",
    "BEL_I": "LTR/Pao", "BEL_LTR": "LTR/Pao", "DIVER2_I": "LTR/Pao",
    "DIVER2_LTR": "LTR/Pao", "DIVER_I": "LTR/Pao", "DIVER_LTR": "LTR/Pao",
    "MAX_I": "LTR/Pao", "MAX_LTR": "LTR/Pao", "NINJA_I": "LTR/Pao",
    "NINJA_LTR": "LTR/Pao", "ROOA_I": "LTR/Pao", "ROOA_LTR": "LTR/Pao",
    "ROO_I": "LTR/Pao", "ROO_LTR": "LTR/Pao", "TRAM_I": "LTR/Pao",
    "TRAM_LTR": "LTR/Pao", "TLD1": "LTR/Gypsy", "TLD2": "LTR/Gypsy",
    "FROGGER_I": "LTR/Copia", "FROGGER_LTR": "LTR/Copia", "QUASIMODO2-I_DM": "LTR/Gypsy",
    "QUASIMODO2-LTR_DM": "LTR/Gypsy", "Chimpo_LTR": "LTR/Gypsy", "Chimpo_I": "LTR/Gypsy",
    "Copia1-I_DM": "LTR/Copia", "Copia1-LTR_DM": "LTR/Copia", "Gypsy2-I_DM": "LTR/Gypsy",
    "Gypsy_6B": "LTR/Gypsy", "Pifo_LTR": "LTR/Gypsy", "Pifo_I": "LTR/Gypsy",
    "Chouto_LTR": "LTR/Gypsy", "Chouto_I": "LTR/Gypsy", "Bica_LTR": "LTR/Gypsy",
    "Bica_I": "LTR/Gypsy", "Nobel_I": "LTR/Pao", "Nobel_LTR": "LTR/Pao",
    "BEL3_DM-LTR": "LTR/Pao", "BEL3_DM-I": "LTR/Pao", "Gypsy1-LTR_DM": "LTR/Gypsy",
    "Gypsy1-I_DM": "LTR/Gypsy", "Gypsy20-I_Dya": "LTR/Gypsy", "Gypsy20-LTR_Dya": "LTR/Gypsy",
    "Gypsy-24_DYa-LTR": "LTR/Gypsy", "Gypsy-24_DYa-I": "LTR/Gypsy",
    "Gypsy-8_DEl-I": "LTR/Gypsy", "Gypsy-8_DEl-LTR": "LTR/Gypsy",
    "Gypsy-26_DYa-LTR": "LTR/Gypsy", "Gypsy-26_DYa-I": "LTR/Gypsy",
    "Gypsy-1_DSim-LTR": "LTR/Gypsy", "Gypsy-1_DSim-I": "LTR/Gypsy",
    "Gypsy-6_DSe-LTR": "LTR/Gypsy", "Gypsy-6_DSe-I": "LTR/Gypsy",
    "Gypsy-4_DSim-I": "LTR/Gypsy", "Gypsy-4_DSim-LTR": "LTR/Gypsy",
}

# =============================================================================
# HELPER LOGIC
# =============================================================================

def clean_family_name(class_str):
    if pd.isna(class_str): return "Unknown"
    class_str = str(class_str)
    if '/' in class_str:
        return class_str.split('/')[-1]
    return class_str

def get_rev_comp(seq):
    complement = str.maketrans('ACGTacgt', 'TGCAtgca')
    return seq.translate(complement)[::-1]

def get_rotations(seq):
    rotations = set()
    n = len(seq)
    for i in range(n):
        rotations.add(seq[i:] + seq[:i])
    return rotations

def build_variant_map(targets):
    mapping = {}
    for t in targets:
        if all(c in 'ACGTacgt' for c in t):
            for rot in get_rotations(t):
                mapping[rot] = t
                mapping[rot.upper()] = t
                mapping[rot.lower()] = t
            rc = get_rev_comp(t)
            for rot in get_rotations(rc):
                mapping[rot] = t
                mapping[rot.upper()] = t
                mapping[rot.lower()] = t
        else:
            mapping[t] = t
    return mapping

def clean_repeat_name(name):
    if pd.isna(name): return name
    name = str(name)
    if name.startswith('(') and ')n' in name:
        return name.replace('(', '').replace(')n', '')
    return name

SAT_MAP = build_variant_map(TARGET_SATELLITES)

def merge_intervals(df, gap_tol=500):
    """
    Merges adjacent repeats of the same name to calculate true array length.
    df must be sorted by Contig and Start.
    Returns a dataframe with ['Contig', 'Repeat_Name', 'Array_Length'].
    """
    merged_rows = []
    
    if df.empty: return pd.DataFrame(columns=['Contig', 'Repeat_Name', 'Array_Length'])
    
    df = df.sort_values(['Contig', 'Start'])
    
    # Iteration
    for (contig, rep), group in df.groupby(['Contig', 'Repeat_Name']):
        group = group.sort_values('Start')
        
        current_start = None
        current_end = None
        
        for _, row in group.iterrows():
            start = row['Start']
            end = row['End']
            
            if current_start is None:
                current_start = start
                current_end = end
                continue
                
            # Check overlap or proximity
            if start <= current_end + gap_tol:
                # Merge
                current_end = max(current_end, end)
            else:
                # Close current array, start new
                merged_rows.append({
                    'Contig': contig,
                    'Repeat_Name': rep,
                    'Array_Length': current_end - current_start + 1
                })
                current_start = start
                current_end = end
                
        # Append last one
        if current_start is not None:
            merged_rows.append({
                'Contig': contig,
                'Repeat_Name': rep,
                'Array_Length': current_end - current_start + 1
            })
            
    return pd.DataFrame(merged_rows)

# =============================================================================
# PARSING
# =============================================================================

def parse_rm_out(filepath, sample_name, repeat_type):
    if not os.path.exists(filepath):
        print(f"Warning: {filepath} not found.")
        return pd.DataFrame()

    try:
        column_names = [
            'score', 'div', 'del', 'ins', 'query', 'begin', 'end', 'left', 
            'strand', 'repeat', 'class', 'r_begin', 'r_end', 'r_left', 'ID', 'star'
        ]
        
        df = pd.read_csv(
            filepath, 
            sep=r'\s+',  # Fixed to regex raw string to prevent Pandas FutureWarnings
            skiprows=3, 
            header=None, 
            names=column_names, 
            engine='python'
        )
        
        if df.empty: return pd.DataFrame()

        df = df[['query', 'begin', 'end', 'repeat', 'class']].copy()
        df.columns = ['Contig', 'Start', 'End', 'Repeat_Name', 'Repeat_Class']
        
        df['Start'] = pd.to_numeric(df['Start'], errors='coerce')
        df['End'] = pd.to_numeric(df['End'], errors='coerce')
        df['Length_bp'] = df['End'] - df['Start'] + 1
        df['Sample'] = sample_name
        df['Type_Category'] = repeat_type 
        
        df['Repeat_Name'] = df['Repeat_Name'].apply(clean_repeat_name)
        if repeat_type == 'Satellite':
            df['Repeat_Name'] = df['Repeat_Name'].apply(lambda x: SAT_MAP.get(x, x))
        
        df['Clean_Family'] = df['Repeat_Class'].apply(clean_family_name)
        
        return df

    except Exception as e:
        print(f"Error parsing {filepath}: {e}")
        return pd.DataFrame()

def summarize_abundance(full_df, group_col):
    summary = full_df.groupby(['Sample', group_col])['Length_bp'].sum().reset_index()
    pivot = summary.pivot(index=group_col, columns='Sample', values='Length_bp').fillna(0)
    return pivot

# =============================================================================
# MAIN PIPELINE
# =============================================================================

def run_pipeline(args):
    print("--- Starting Comparative Pipeline ---")
    
    samples = {
        'ISO-1': {'TE': args.iso1_te, 'Sat': args.iso1_sat},
        'A3':    {'TE': args.a3_te,   'Sat': args.a3_sat},
        'A4':    {'TE': args.a4_te,   'Sat': args.a4_sat}
    }
    
    all_data = []
    
    for strain, files in samples.items():
        print(f"Processing {strain}...")
        te_df = parse_rm_out(files['TE'], strain, 'TE')
        sat_df = parse_rm_out(files['Sat'], strain, 'Satellite')
        combined = pd.concat([te_df, sat_df], ignore_index=True)
        if not combined.empty:
            all_data.append(combined)

    if not all_data:
        print("Error: No data parsed from any files.")
        sys.exit(1)

    master_df = pd.concat(all_data, ignore_index=True)
    
    if master_df.empty:
        print("Error: The final merged dataframe is empty.")
        sys.exit(1)

    print(f"Total elements parsed: {len(master_df)}")

    # 3. Generate Statistics & Figures
    # AESTHETIC UPGRADE: Set publication-quality context
    sns.set_theme(style="ticks", font_scale=1.3, rc={"axes.grid": True, "grid.alpha": 0.3})
    plt.rcParams['font.family'] = 'sans-serif'
    # Fallback fonts if Arial is not installed on the server
    plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans', 'Liberation Sans'] 

    te_only = master_df[master_df['Type_Category'] == 'TE']
    sat_only = master_df[master_df['Type_Category'] == 'Satellite']

    # --- 1. FIGURES ---
    
    # Figure 1: TE Heatmap (Abundance Density)
    print("Generating Fig 1: TE Heatmap...")
    
    # Filter out Satellites and Simple Repeats from TE data to ensure purity
    ignored_families = ['Satellite', 'Simple_repeat', 'Low_complexity']
    te_clean = te_only[~te_only['Clean_Family'].isin(ignored_families)]
    
    te_family_stats = te_clean.groupby(['Clean_Family', 'Sample'])['Length_bp'].sum().reset_index()
    te_pivot = te_family_stats.pivot(index='Clean_Family', columns='Sample', values='Length_bp').fillna(0)
    top_families = te_pivot.sum(axis=1).nlargest(25).index
    
    plt.figure(figsize=(10, 12))
    # Using Okabe-Ito Blue (#0072B2) to generate a sequential gradient for the Heatmap 
    cmap_okabe = sns.light_palette("#0072B2", as_cmap=True)
    sns.heatmap(te_pivot.loc[top_families], cmap=cmap_okabe, annot=True, fmt=".0f", 
                linewidths=0.5, linecolor='lightgray', cbar_kws={'label': 'Total Base Pairs'})
    plt.title('Genomic Abundance of Top 25 TE Families', fontsize=16, weight='bold')
    plt.ylabel('TE Family', fontsize=14, weight='bold')
    plt.xlabel('Strain', fontsize=14, weight='bold')
    plt.tight_layout()
    plt.savefig('Figure_1_TE_Heatmap.svg', format='svg', bbox_inches='tight') # bbox tight ensures labels are safe
    plt.close()
    
    # Figure 2: Satellite Bar (Targeted)
    print("Generating Fig 2: Satellite Bar...")
    filtered_sat = sat_only[sat_only['Repeat_Name'].isin(TARGET_SATELLITES)]
    
    if not filtered_sat.empty:
        sat_agg = filtered_sat.groupby(['Repeat_Name', 'Sample'])['Length_bp'].sum().reset_index()
        plt.figure(figsize=(14, 8))
        # Explicitly passing the categorical Okabe-Ito color palette
        sns.barplot(data=sat_agg, x='Repeat_Name', y='Length_bp', hue='Sample', palette=OKABE_ITO, edgecolor='black', linewidth=0.8)
        plt.title('Abundance of Target Satellite Repeats', fontsize=16, weight='bold')
        plt.ylabel('Total Base Pairs', fontsize=14, weight='bold')
        plt.xlabel('Satellite Name', fontsize=14, weight='bold')
        plt.xticks(rotation=45, ha='right')
        plt.grid(axis='y', linestyle='--', alpha=0.5)
        plt.legend(title='Strain', fontsize=12)
        plt.tight_layout()
        plt.savefig('Figure_2_Satellite_Targeted.svg', format='svg', bbox_inches='tight')
        plt.close()
    else:
        print("Warning: No matching targeted satellites found for Figure 2.")
    
    # --- 2. REPORTS ---
    print(f"Writing reports to {args.output}...")
    # Explicitly calling 'openpyxl' to prevent missing engine crashes
    with pd.ExcelWriter(args.output, engine='openpyxl') as writer:
        
        # A. Array Statistics for SIMPLE REPEATS (Hotspots)
        simple_sats = sat_only[sat_only['Repeat_Name'].isin(SIMPLE_REPEAT_TARGETS)].copy()
        
        if not simple_sats.empty:
            # 1. Total Abundance per Contig
            simple_contig_stats = simple_sats.groupby(['Sample', 'Contig', 'Repeat_Name'])['Length_bp'].sum().reset_index()
            simple_contig_stats.rename(columns={'Length_bp': 'Total_Abundance_bp'}, inplace=True)
            
            # 2. Longest Array Calculation (Merging)
            array_stats_list = []
            for sample_name in simple_sats['Sample'].unique():
                sample_subset = simple_sats[simple_sats['Sample'] == sample_name]
                merged = merge_intervals(sample_subset)
                merged['Sample'] = sample_name
                if not merged.empty:
                    max_arrays = merged.groupby(['Contig', 'Repeat_Name'])['Array_Length'].max().reset_index()
                    max_arrays['Sample'] = sample_name
                    array_stats_list.append(max_arrays)
                
            if array_stats_list:
                all_max_arrays = pd.concat(array_stats_list)
                final_simple_stats = pd.merge(simple_contig_stats, all_max_arrays, on=['Sample', 'Contig', 'Repeat_Name'], how='outer')
                final_simple_stats.sort_values(['Sample', 'Repeat_Name', 'Total_Abundance_bp'], ascending=[True, True, False], inplace=True)
                final_simple_stats.to_excel(writer, sheet_name='Simple_Satellite_Hotspots', index=False)
        
        # B. Library Specific Abundance
        lib_df = pd.DataFrame(list(LIBRARY_ELEMENTS.keys()), columns=['Repeat_Name'])
        lib_df['Class'] = list(LIBRARY_ELEMENTS.values())
        for strain in ['ISO-1', 'A3', 'A4']:
            strain_data = master_df[master_df['Sample'] == strain]
            strain_counts = strain_data.groupby('Repeat_Name')['Length_bp'].sum()
            lib_df[f'{strain}_bp'] = lib_df['Repeat_Name'].map(strain_counts).fillna(0)
        lib_df.to_excel(writer, sheet_name='Library_Specific_Abundance', index=False)

        # C. General Summaries
        summarize_abundance(te_only, 'Repeat_Class').to_excel(writer, sheet_name='Summary_TE_Family_bp')
        summarize_abundance(sat_only, 'Repeat_Name').to_excel(writer, sheet_name='Summary_All_Sat_bp')
        
        # D. Detailed Locations
        for strain in ['ISO-1', 'A3', 'A4']:
            det = master_df[master_df['Sample'] == strain].sort_values(['Contig', 'Start'])
            det.to_excel(writer, sheet_name=f'{strain}_Detailed_Locations', index=False)
            
        # E. Top Repetitive Contigs
        contig_stats = master_df.groupby(['Sample', 'Contig', 'Type_Category'])['Length_bp'].sum().unstack(fill_value=0)
        contig_stats['Total_Repeat_bp'] = contig_stats.sum(axis=1)
        contig_stats.sort_values('Total_Repeat_bp', ascending=False).head(100).to_excel(writer, sheet_name='Top_Repetitive_Contigs')

    print("Pipeline Complete.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--iso1_te", required=True)
    parser.add_argument("--iso1_sat", required=True)
    parser.add_argument("--a3_te", required=True)
    parser.add_argument("--a3_sat", required=True)
    parser.add_argument("--a4_te", required=True)
    parser.add_argument("--a4_sat", required=True)
    parser.add_argument("--output", default="Comparative_Y_Extensive_Data.xlsx")
    args = parser.parse_args()
    run_pipeline(args)
