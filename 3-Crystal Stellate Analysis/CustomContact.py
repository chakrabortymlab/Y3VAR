import cooler
import pandas as pd
import numpy as np
import plotly.graph_objects as go

# ==========================================
# 1. CONFIGURATION PARAMETERS
# ==========================================
COOL_FILE = 'CantonS_inter_25kb_transformed.cool'
CHROM = 'Y_ptg000011l'
ROI_START = 2852210
ROI_END = 5180375

print("Loading Hi-C data...")
# ==========================================
# 2. EXTRACT DIRECTLY FROM .COOL 
# ==========================================
c = cooler.Cooler(COOL_FILE)

# Fetch the pixels and bin coordinates, then join them together
pixels_df = c.pixels().fetch(CHROM)
bins_df = c.bins().fetch(CHROM)
df = cooler.annotate(pixels_df, bins_df)

print("Filtering to your specific region...")
# ==========================================
# 3. FILTER TO REGION OF INTEREST
# ==========================================
df = df[(df['start1'] >= ROI_START) & (df['start1'] <= ROI_END)]
df = df[(df['start2'] >= ROI_START) & (df['start2'] <= ROI_END)]

# Keep only the upper triangle
df = df[df['start2'] >= df['start1']]

# Convert coordinates to Megabases (Mb)
df['bin1_mb'] = df['start1'] / 1_000_000
df['bin2_mb'] = df['start2'] / 1_000_000

print("Transforming into a triangular peak...")
# ==========================================
# 4. TRIANGULAR TRANSFORMATION
# ==========================================
df['midpoint_mb'] = (df['bin1_mb'] + df['bin2_mb']) / 2
df['distance_mb'] = df['bin2_mb'] - df['bin1_mb']

# Reshape into the peak grid
pivot_df = df.pivot(index='distance_mb', columns='midpoint_mb', values='count')

# Sort axes logically
pivot_df = pivot_df.sort_index(ascending=True)
pivot_df = pivot_df[sorted(pivot_df.columns)]

# Extract matrix array and fill empty background spaces with 0
dense_matrix = pivot_df.values
dense_matrix = np.nan_to_num(dense_matrix) 

print("Calculating optimal color contrast...")
# ==========================================
# 5. DYNAMIC SCALING & PLOTTING
# ==========================================
# Extract true contacts to calculate the 95th percentile color saturation limit
true_contacts = dense_matrix[dense_matrix > 0]
vmax_threshold = np.percentile(true_contacts, 95) if len(true_contacts) > 0 else 1.5

# --- BRIGHTENED OKABE-ITO DIVERGENT COLORSCALE ---
okabe_ito_bright = [
    [0.0, '#0072B2'],  # Okabe-Ito Blue (Dark canvas background)
    [0.5, '#78C9FF'],  # Brightened Sky Blue (Makes baseline contacts pop)
    [1.0, '#FF6600']   # Brightened Vermilion (Blazing peak enrichment)
]

fig = go.Figure(data=go.Heatmap(
    z=dense_matrix,
    x=pivot_df.columns,
    y=pivot_df.index,
    colorscale=okabe_ito_bright,  
    zmin=0,                 # 0 maps to Dark Blue
    zmid=1.0,               # Anchors 1.0 (Expected baseline) to Bright Sky Blue
    zmax=vmax_threshold,    # Enriched peaks map to Bright Vermilion
    colorbar=dict(
        title=dict(text="Obs/Exp Ratio", side="right"),
        thickness=20,
        len=0.7
    ),
    hovertemplate=(
        "Midpoint: %{x:.3f} Mb<br>" +
        "Interaction Distance: %{y:.3f} Mb<br>" +
        "Transformation Value: %{z:.4f}<extra></extra>"
    )
))

# Apply layout optimizations
fig.update_layout(
    title={
        'text': f"<b>Su(Ste) / PCKR Region Contact Map (1 kb Res)</b><br><sup>Coordinates: {ROI_START/1_000_000:.3f} Mb - {ROI_END/1_000_000:.3f} Mb</sup>",
        'y': 0.95,
        'x': 0.5,
        'xanchor': 'center',
        'yanchor': 'top'
    },
    xaxis_title="Genomic Coordinate (Mb)",
    yaxis_title="Interaction Distance (Mb)",
    width=1000,
    height=600,
    xaxis=dict(showgrid=False, zeroline=False),
    yaxis=dict(showgrid=False, zeroline=False, autorange=True),
    # Match the background to the Okabe-Ito Blue for seamless edges
    plot_bgcolor='#0072B2',   
    paper_bgcolor='white'
)

print("Saving plots to file...")
# ==========================================
# 6. EXPORT AND SHOW
# ==========================================
# Save as a high-resolution PNG (scale=2 doubles the pixel density for publications)
fig.write_image("SuSte_PCKR_ContactMap_OkabeIto_Bright.png", scale=2)

# Save as a vector SVG
fig.write_image("SuSte_PCKR_ContactMap_OkabeIto_Bright.svg")

# Launch the interactive browser version
fig.show()

print("Success! Interactive plot opened and PNG/SVG files saved in your directory.")
