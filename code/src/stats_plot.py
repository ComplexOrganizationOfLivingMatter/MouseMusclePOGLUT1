import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind
from pathlib import Path
import numpy as np

# --- CONFIGURATION ---
MODEL_FOLDER = 'modelV2'
excel_files = list(Path(f"../../results/{MODEL_FOLDER}").glob('Pax7_Genotype_Analysis_*.xlsx'))
EXCEL_PATH = excel_files[0]
OUTPUT_DIR = Path(f"../../results/{MODEL_FOLDER}/plots_stats")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

AGES_ORDER = ['p0','p4', 'p7', 'p10', 'p15', 'p20', 'p35']
GENOTYPE_ORDER = ['WT', 'KIKO']
PALETTE = {'WT': '#0072B2', 'KIKO': '#D55E00'}

# --- LOAD DATA ---
df = pd.read_excel(EXCEL_PATH)
df['age'] = pd.Categorical(df['age'], categories=AGES_ORDER, ordered=True)
df['genotype'] = pd.Categorical(df['genotype'], categories=GENOTYPE_ORDER, ordered=True)

# Sort to align Dataframe rows with Plot points
df = df.sort_values(['age', 'genotype']).reset_index(drop=True)

density_metrics = ['whole muscle nuclei/um2'] + [c for c in df.columns if 'nuclei/um2 g' in c]
all_metrics = [m for m in density_metrics if 'entropy' not in m.lower()]

def add_stat_annotation(ax, x1, x2, y, p_val):
    if p_val < 0.001: stars = "***"
    elif p_val < 0.01: stars = "**"
    elif p_val < 0.05: stars = "*"
    else: return 
    h = y * 0.02
    ax.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.2, c='k')
    ax.text((x1+x2)*.5, y+h, stars, ha='center', va='bottom', color='k', fontsize=11)

# --- PLOTTING LOOP ---
for metric in all_metrics:
    if metric not in df.columns or df[metric].dropna().empty: continue

    plt.figure(figsize=(14, 8), dpi=300)
    
    # 1. Plot raw points (Jitter + Dodge)
    strip = sns.stripplot(data=df, x='age', y=metric, hue='genotype', 
                         hue_order=GENOTYPE_ORDER, dodge=True, 
                         alpha=0.4, palette=PALETTE, jitter=0.2, size=5)
    
    # 2. Plot Mean and SD
    sns.pointplot(data=df, x='age', y=metric, hue='genotype', 
                  hue_order=GENOTYPE_ORDER, dodge=0.4, linestyle='none', 
                  palette=PALETTE, errorbar='sd', capsize=0.1, 
                  markers="_", markersize=10)

    ax = plt.gca()

    # --- UPDATED LABELING: FIRST 2 WORDS ---
    for i, genotype in enumerate(GENOTYPE_ORDER):
        offsets = strip.collections[i].get_offsets()
        mask = (df['genotype'] == genotype) & (df[metric].notna())
        subset = df[mask].copy()
        
        if len(offsets) == len(subset):
            for j in range(len(offsets)):
                full_name = Path(str(subset.iloc[j]['image'])).stem
                # Split by space and take first 2 words
                words = full_name.split(' ')
                label = " ".join(words[:2])
                
                ax.text(offsets[j, 0] + 0.04, offsets[j, 1], label, 
                        fontsize=5.5, alpha=0.8, rotation=35, 
                        ha='left', va='bottom')
        else:
            # Manual fallback if point counts diverge
            for age_idx, age in enumerate(AGES_ORDER):
                age_subset = subset[subset['age'] == age]
                x_base = age_idx - 0.2 if genotype == 'WT' else age_idx + 0.2
                
                for k in range(len(age_subset)):
                    val = age_subset.iloc[k][metric]
                    full_name = Path(str(age_subset.iloc[k]['image'])).stem
                    label = " ".join(full_name.split(' ')[:2])
                    
                    ax.text(x_base + 0.05, val, label, 
                            fontsize=5, alpha=0.7, rotation=35)

    # 3. Statistical Comparison
    for i, age in enumerate(AGES_ORDER):
        sub = df[df['age'] == age]
        wt_vals = sub[sub['genotype'] == 'WT'][metric].dropna()
        ki_vals = sub[sub['genotype'] == 'KIKO'][metric].dropna()
        
        if len(wt_vals) >= 2 and len(ki_vals) >= 2:
            _, p_val = ttest_ind(wt_vals, ki_vals)
            local_max = sub[metric].max()
            if not np.isnan(local_max):
                y_pos = local_max * 1.25 # Increased to clear labels
                add_stat_annotation(ax, i - 0.2, i + 0.2, y_pos, p_val)

    # Styling
    plt.title(f'Comparison: {metric.replace("_", " ").title()}', fontsize=14, pad=35)
    plt.ylabel(r'Density ($\mu m^{-2}$)', fontsize=12)
    plt.xlabel('Postnatal Age', fontsize=12)
    
    handles, labels = ax.get_legend_handles_labels()
    plt.legend(handles[0:2], labels[0:2], title='Genotype', frameon=False, loc='upper right')
    sns.despine()
    
    safe_name = metric.replace('/', '_').replace(' ', '_')
    plt.savefig(OUTPUT_DIR / f"{safe_name}.png", bbox_inches='tight')
    plt.close()

print(f"Success! Density plots with 2-word labels generated in {OUTPUT_DIR}")