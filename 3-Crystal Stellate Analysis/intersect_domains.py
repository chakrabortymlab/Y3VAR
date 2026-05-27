import cooler
import pandas as pd
import numpy as np
import itertools

# ==========================================
# 1. CONFIGURATION PARAMETERS
# ==========================================
COOL_FILE = 'CantonS_inter_25kb_transformed.cool'  
ISLAND_FILE = 'Island_Summary_ISO-1.tsv'
CHROM = 'Y_ptg000011l'
OUTPUT_FILE = 'Domain_Contact_Frequencies_CantonS.tsv'

print(f"Loading annotation file: {ISLAND_FILE}")
# ==========================================
# 2. LOAD ISLAND ANNOTATIONS
# ==========================================
islands_df = pd.read_csv(ISLAND_FILE, sep='\t')

# Isolate just the sequential SuSte islands for neighbor and network analysis
suste_islands = [f"SuSte_{i}" for i in range(1, 8)]

domains = []
for idx, row in islands_df.iterrows():
    domains.append({
        'name': row['Structural_Domain'],
        'start': int(row['Start_Coordinate']),
        'end': int(row['End_Coordinate'])
    })

print(f"Loading Hi-C matrix: {COOL_FILE}")
# ==========================================
# 3. LOAD HI-C DATA
# ==========================================
c = cooler.Cooler(COOL_FILE)

# ==========================================
# 4. COMPUTE PAIRWISE CONTACTS
# ==========================================
print("Calculating contacts between all domain pairs...")
results = []

for dom_a, dom_b in itertools.combinations_with_replacement(domains, 2):
    region1 = (CHROM, dom_a['start'], dom_a['end'])
    region2 = (CHROM, dom_b['start'], dom_b['end'])
    
    try:
        sub_matrix = c.matrix(balance=False).fetch(region1, region2)
        if sub_matrix.size > 0:
            mean_contact = np.nanmean(sub_matrix)
            sum_contact = np.nansum(sub_matrix)
            valid_bins = np.count_nonzero(~np.isnan(sub_matrix))
        else:
            mean_contact, sum_contact, valid_bins = 0, 0, 0
    except Exception as e:
        mean_contact, sum_contact, valid_bins = 0, 0, 0

    interaction_type = "Intra-Domain" if dom_a['name'] == dom_b['name'] else "Inter-Domain"

    results.append({
        'Domain_A': dom_a['name'],
        'Domain_B': dom_b['name'],
        'Interaction_Type': interaction_type,
        'Mean_Obs_Exp_Ratio': round(mean_contact, 4),
        'Total_Contact_Sum': round(sum_contact, 4),
        'Number_of_Interacting_Bins': valid_bins
    })

results_df = pd.DataFrame(results)
# Sort baseline pairwise data by interaction strength
results_df = results_df.sort_values(by='Mean_Obs_Exp_Ratio', ascending=False)

# ==========================================
# 5. GENERATE NETWORK AND NEIGHBOR METRICS
# ==========================================
print("Generating pattern network and linear neighbor rows (O/E and Absolute Sums)...")

# Dictionaries to track both metrics [Mean_Obs_Exp, Total_Sum]
intra_oe, intra_sum = [], []
neighbor_oe, neighbor_sum = [], []
odd_oe, odd_sum = [], []
even_oe, even_sum = [], []

odds = ['SuSte_1', 'SuSte_3', 'SuSte_5', 'SuSte_7']
evens = ['SuSte_2', 'SuSte_4', 'SuSte_6']

# Helper maps to find values easily
oe_map = {}
sum_map = {}
for idx, row in results_df.iterrows():
    pair = tuple(sorted([row['Domain_A'], row['Domain_B']]))
    oe_map[pair] = row['Mean_Obs_Exp_Ratio']
    sum_map[pair] = row['Total_Contact_Sum']

# Parse through our target SuSte list to categorize interactions
for i in range(len(suste_islands)):
    d1 = suste_islands[i]
    
    # 1. Collect Intra-Island entries
    intra_oe.append(oe_map.get((d1, d1), 0))
    intra_sum.append(sum_map.get((d1, d1), 0))
    
    for j in range(i + 1, len(suste_islands)):
        d2 = suste_islands[j]
        val_oe = oe_map.get((d1, d2), 0)
        val_sum = sum_map.get((d1, d2), 0)
        
        # 2. Collect Linear Neighbors
        if j == i + 1:
            neighbor_oe.append(val_oe)
            neighbor_sum.append(val_sum)
        
        # 3. Collect Odd vs Even Networks
        if d1 in odds and d2 in odds:
            odd_oe.append(val_oe)
            odd_sum.append(val_sum)
        elif d1 in evens and d2 in evens:
            even_oe.append(val_oe)
            even_sum.append(val_sum)

# Compute the global averages/sums
summary_rows = [
    {
        'Domain_A': 'SUMMARY_METRIC',
        'Domain_B': 'OVERALL_INTRA_ISLAND_METRIC',
        'Interaction_Type': 'Summary',
        'Mean_Obs_Exp_Ratio': round(np.mean(intra_oe), 4) if intra_oe else 0,
        'Total_Contact_Sum': round(np.sum(intra_sum), 4) if intra_sum else 0,
        'Number_of_Interacting_Bins': 0
    },
    {
        'Domain_A': 'SUMMARY_METRIC',
        'Domain_B': 'OVERALL_LINEAR_NEIGHBOR_METRIC',
        'Interaction_Type': 'Summary',
        'Mean_Obs_Exp_Ratio': round(np.mean(neighbor_oe), 4) if neighbor_oe else 0,
        'Total_Contact_Sum': round(np.sum(neighbor_sum), 4) if neighbor_sum else 0,
        'Number_of_Interacting_Bins': 0
    },
    {
        'Domain_A': 'SUMMARY_METRIC',
        'Domain_B': 'ODD_NETWORK_METRIC_(1,3,5,7)',
        'Interaction_Type': 'Summary',
        'Mean_Obs_Exp_Ratio': round(np.mean(odd_oe), 4) if odd_oe else 0,
        'Total_Contact_Sum': round(np.sum(odd_sum), 4) if odd_sum else 0,
        'Number_of_Interacting_Bins': 0
    },
    {
        'Domain_A': 'SUMMARY_METRIC',
        'Domain_B': 'EVEN_NETWORK_METRIC_(2,4,6)',
        'Interaction_Type': 'Summary',
        'Mean_Obs_Exp_Ratio': round(np.mean(even_oe), 4) if even_oe else 0,
        'Total_Contact_Sum': round(np.sum(even_sum), 4) if even_sum else 0,
        'Number_of_Interacting_Bins': 0
    }
]

# Convert summaries to DataFrame and append to bottom of our true grid
summary_df = pd.DataFrame(summary_rows)
final_df = pd.concat([results_df, summary_df], ignore_index=True)

# ==========================================
# 6. EXPORT TO TRUE TSV SPREADSHEET
# ==========================================
final_df.to_csv(OUTPUT_FILE, sep='\t', index=False)
print(f"\nSuccess! Complete TSV spreadsheet with statistical summaries saved to: {OUTPUT_FILE}")
