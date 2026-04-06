import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind
from pathlib import Path
import numpy as np

# --- CONFIGURATION ---
EXCEL_PATH = list(Path('../results/modelV2').glob('Pax7_Genotype_Analysis_*.xlsx'))[0]
OUTPUT_DIR = Path('../results/modelV2/plots_stats')
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

AGES_ORDER = ['p4', 'p7', 'p10', 'p15', 'p20', 'p35']
GENOTYPE_ORDER = ['WT', 'KIKO']  # Forces WT left, KIKO right
PALETTE = {'WT': '#0072B2', 'KIKO': '#D55E00'} # Blue and Orange

# --- LOAD DATA ---
df = pd.read_excel(EXCEL_PATH)
df['age'] = pd.Categorical(df['age'], categories=AGES_ORDER, ordered=True)
df['genotype'] = pd.Categorical(df['genotype'], categories=GENOTYPE_ORDER, ordered=True)

# Automatically identify all metrics (Density, Entropy, and Normalized Entropy)
density_metrics = ['whole muscle nuclei/um2'] + [c for c in df.columns if 'nuclei/um2 g' in c]
entropy_metrics = [c for c in df.columns if 'normalized_entropy' in c]
all_metrics = density_metrics + entropy_metrics

def add_stat_annotation(ax, x1, x2, y, p_val):
    """Adds brackets and significance stars to the plot."""
    if p_val < 0.001: stars = "***"
    elif p_val < 0.01: stars = "**"
    elif p_val < 0.05: stars = "*"
    else: return 

    h = y * 0.02
    ax.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.2, c='k')
    ax.text((x1+x2)*.5, y+h, stars, ha='center', va='bottom', color='k', fontsize=11)

# --- PLOTTING LOOP ---
for metric in all_metrics:
    if metric not in df.columns: continue
    
    # Check if the metric has any data
    if df[metric].dropna().empty: continue

    plt.figure(figsize=(10, 6), dpi=300)
    
    # 1. Plot raw points with jitter
    # hue_order ensures WT is left, KIKO is right
    sns.stripplot(data=df, x='age', y=metric, hue='genotype', 
                  hue_order=GENOTYPE_ORDER, dodge=True, 
                  alpha=0.3, palette=PALETTE, jitter=0.2, size=4)
    
    # 2. Plot Mean and SD
    sns.pointplot(data=df, x='age', y=metric, hue='genotype', 
                  hue_order=GENOTYPE_ORDER, dodge=0.4, join=False, 
                  palette=PALETTE, errorbar='sd', capsize=0.1, 
                  markers="_", scale=1.2)

    ax = plt.gca()
    
    # 3. Statistical Comparison per Age Group
    for i, age in enumerate(AGES_ORDER):
        sub = df[df['age'] == age]
        wt_vals = sub[sub['genotype'] == 'WT'][metric].dropna()
        ki_vals = sub[sub['genotype'] == 'KIKO'][metric].dropna()
        
        if len(wt_vals) >= 3 and len(ki_vals) >= 3:
            _, p_val = ttest_ind(wt_vals, ki_vals)
            
            # Dynamically position the star above the local max
            local_max = sub[metric].max()
            if not np.isnan(local_max):
                y_pos = local_max * 1.05
                # Adjusting x-offsets to match the dodged points
                add_stat_annotation(ax, i - 0.2, i + 0.2, y_pos, p_val)

    # Styling based on metric type
    clean_title = metric.replace("_", " ").title()
    plt.title(f'Genotype Comparison: {clean_title}', fontsize=14, pad=20)
    
    if 'entropy' in metric.lower():
        plt.ylabel('Spatial Entropy (Disorder Index)', fontsize=12)
    else:
        plt.ylabel(r'Density ($\mu m^{-2}$)', fontsize=12)

    plt.xlabel('Postnatal Age', fontsize=12)
    plt.legend(title='Genotype', frameon=False, loc='upper right')
    sns.despine()
    
    # Save with safe filename
    safe_name = metric.replace('/', '_').replace(' ', '_')
    plt.savefig(OUTPUT_DIR / f"{safe_name}.png", bbox_inches='tight')
    plt.close()

print(f"Successfully generated {len(all_metrics)} plots in {OUTPUT_DIR}")