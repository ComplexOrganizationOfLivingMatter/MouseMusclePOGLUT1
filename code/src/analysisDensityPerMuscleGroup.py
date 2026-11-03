"""
Comprehensive Plotting Pipeline for Multi-Feature Muscle Analysis
Reads pre-calculated features from the multi-sheet Excel report.
No calculations are performed from scratch.
"""

import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# ============================================================
# CONFIGURATION & CONSTANTS
# ============================================================
BASE_DIR = Path("../../")
RESULT_DIR = BASE_DIR / "results" / "largeModel_dotDetection"

# Set plotting style & consistent colors/order
sns.set_theme(style="whitegrid")
plt.rcParams.update({'font.sans-serif': 'Arial', 'font.size': 11})

GENOTYPE_ORDER = ["WT", "KIKO"]
GENOTYPE_BOX_PALETTE = {"WT": "#a6c8e0", "KIKO": "#fba375"}
GENOTYPE_DOT_PALETTE = {"WT": "#1a426e", "KIKO": "#b33d00"}

AGE_SAMPLE_PERCENTAGE_THRESHOLD = 0.20


def load_latest_excel_report(result_dir: Path) -> dict[str, pd.DataFrame]:
    """Finds and loads all sheets from the most recent multi-feature Excel report."""
    excel_files = sorted(result_dir.glob("MultiFeature_Analysis_*.xlsx"))
    if not excel_files:
        excel_files = sorted(result_dir.glob("*.xlsx"))
        excel_files = [f for f in excel_files if not f.name.startswith("~") and "outlier" not in f.name.lower()]
        
    if not excel_files:
        raise FileNotFoundError(f"No analysis Excel report found in {result_dir}")
        
    latest_file = excel_files[-1]
    print(f"[INFO] Loading data from: {latest_file.name}")
    
    xls = pd.ExcelFile(latest_file)
    return {sheet_name: pd.read_excel(xls, sheet_name=sheet_name) for sheet_name in xls.sheet_names}


# ============================================================
# OUTLIER DETECTION & FILTERING HELPER (Mean ± 2*SD)
# ============================================================
def filter_and_log_outliers(df: pd.DataFrame, value_col: str, metric_name: str) -> tuple[pd.DataFrame, list]:
    """
    Identifies and logs data points further than 2 standard deviations from the mean 
    for each combination of age, genotype, and muscle_group.
    Returns the cleaned dataframe and a list of logged outlier records.
    """
    filtered_groups = []
    outlier_records = []
    group_cols = ["age", "genotype", "muscle_group"]
    
    for _, group in df.groupby(group_cols):
        mean_val = group[value_col].mean()
        std_val = group[value_col].std()
        
        if pd.isna(std_val) or std_val == 0:
            filtered_groups.append(group)
        else:
            lower_bound = mean_val - 2 * std_val
            upper_bound = mean_val + 2 * std_val
            
            clean_group = group[(group[value_col] >= lower_bound) & (group[value_col] <= upper_bound)]
            outliers = group[(group[value_col] < lower_bound) | (group[value_col] > upper_bound)]
            
            filtered_groups.append(clean_group)
            
            for _, row in outliers.iterrows():
                outlier_records.append({
                    "metric": metric_name,
                    "image": row.get("image", "unknown"),
                    "age": row["age"],
                    "genotype": row["genotype"],
                    "muscle_group": row["muscle_group"],
                    "value": row[value_col],
                    "group_mean": mean_val,
                    "group_std": std_val,
                    "lower_bound": lower_bound,
                    "upper_bound": upper_bound
                })
            
    if filtered_groups:
        return pd.concat(filtered_groups, ignore_index=True), outlier_records
    return df, outlier_records


# ============================================================
# PLOTTING HELPER WITH ORGANIZED DIRECTORIES
# ============================================================
def generate_feature_plots(
    df_raw: pd.DataFrame, 
    feature_prefix: str, 
    metric_label: str, 
    subfolder_name: str
) -> list:
    """
    Melts wide per-group columns for a specific feature, cleans up muscle group labels,
    filters outliers (Mean ± 2*SD), and generates shared and independent Y-axis plots.
    """
    feature_cols = [c for c in df_raw.columns if c.startswith(feature_prefix)]
    if not feature_cols:
        print(f"[WARN] No columns found matching prefix: '{feature_prefix}'")
        return []

    print(f"\n[STAGE] Generating plots for: '{metric_label}' ({subfolder_name})")

    melted_df = pd.melt(
        df_raw, 
        id_vars=["image", "genotype", "age"], 
        value_vars=feature_cols, 
        var_name="muscle_group_raw", 
        value_name="value"
    )

    melted_df["muscle_group"] = melted_df["muscle_group_raw"].apply(
        lambda x: f"G{re.search(r'g(\d+)', str(x)).group(1)}" if re.search(r'g(\d+)', str(x)) else str(x)
    )
    
    melted_df = melted_df[~melted_df["muscle_group"].isin(["G0", "G1", "G6", "G13"])].copy()
    
    melted_df["group_num"] = melted_df["muscle_group"].str.extract(r'(\d+)').astype(int)
    melted_df = melted_df.sort_values("group_num")
    melted_df["age"] = melted_df["age"].astype(str).str.strip()

    filtered_df, outlier_records = filter_and_log_outliers(melted_df, "value", subfolder_name)

    global_max = filtered_df["value"].max()
    if np.isnan(global_max) or global_max == 0:
        global_max = 1.0
    shared_y_limit = global_max * 1.15

    base_plot_dir = RESULT_DIR / "plots_density_by_age" / subfolder_name
    shared_dir = base_plot_dir / "shared_yaxis"
    independent_dir = base_plot_dir / "independent_yaxis"
    
    shared_dir.mkdir(parents=True, exist_ok=True)
    independent_dir.mkdir(parents=True, exist_ok=True)

    ages = filtered_df["age"].unique()

    for age in ages:
        df_age = filtered_df[filtered_df["age"] == age].copy()
        if df_age.empty:
            continue

        total_samples = df_age["image"].nunique()
        group_counts = df_age.dropna(subset=["value"]).groupby("muscle_group")["image"].nunique()
        
        valid_groups = [
            grp for grp, count in group_counts.items()
            if (count / total_samples) >= AGE_SAMPLE_PERCENTAGE_THRESHOLD
        ]

        df_age_filtered = df_age[df_age["muscle_group"].isin(valid_groups)]
        if df_age_filtered.empty:
            continue

        present_genotypes = [g for g in GENOTYPE_ORDER if g in df_age_filtered["genotype"].unique()]
        safe_age = str(age).replace(" ", "_").replace("/", "_")

        # --- VARIANT A: SHARED Y-AXIS ---
        try:
            plt.figure(figsize=(14, 6))
            sns.boxplot(
                data=df_age_filtered, x="muscle_group", y="value", hue="genotype",
                hue_order=GENOTYPE_ORDER, palette=GENOTYPE_BOX_PALETTE, boxprops=dict(alpha=0.6), 
                showcaps=True, showfliers=False
            )
            sns.stripplot(
                data=df_age_filtered, x="muscle_group", y="value", hue="genotype",
                hue_order=GENOTYPE_ORDER, palette=GENOTYPE_DOT_PALETTE, dodge=True, jitter=0.2, size=6, linewidth=0.8, edgecolor="black"
            )
            plt.ylim(0, shared_y_limit)
            plt.title(f"{metric_label} per Muscle Group (Shared Scale, Mean±2SD Filtered) — Age: {age}", fontsize=13, fontweight="bold", pad=15)
            plt.xlabel("Muscle Group", fontweight="bold")
            plt.ylabel(metric_label, fontweight="bold")
            
            handles, labels = plt.gca().get_legend_handles_labels()
            if handles:
                by_label = dict(zip(labels[:len(present_genotypes)], handles[:len(present_genotypes)]))
                plt.legend(by_label.values(), by_label.keys(), title="Genotype", frameon=True)

            plt.tight_layout()
            plt.savefig(shared_dir / f"{subfolder_name}_{safe_age}_shared.png", dpi=300)
            plt.close()
        except Exception:
            plt.close()

        # --- VARIANT B: INDEPENDENT Y-AXIS ---
        try:
            plt.figure(figsize=(14, 6))
            sns.boxplot(
                data=df_age_filtered, x="muscle_group", y="value", hue="genotype",
                hue_order=GENOTYPE_ORDER, palette=GENOTYPE_BOX_PALETTE, boxprops=dict(alpha=0.6), 
                showcaps=True, showfliers=False
            )
            sns.stripplot(
                data=df_age_filtered, x="muscle_group", y="value", hue="genotype",
                hue_order=GENOTYPE_ORDER, palette=GENOTYPE_DOT_PALETTE, dodge=True, jitter=0.2, size=6, linewidth=0.8, edgecolor="black"
            )
            
            age_max = df_age_filtered["value"].max()
            if not np.isnan(age_max) and age_max > 0:
                plt.ylim(0, age_max * 1.15)

            plt.title(f"{metric_label} per Muscle Group (Independent Scale, Mean±2SD Filtered) — Age: {age}", fontsize=13, fontweight="bold", pad=15)
            plt.xlabel("Muscle Group", fontweight="bold")
            plt.ylabel(metric_label, fontweight="bold")
            
            handles, labels = plt.gca().get_legend_handles_labels()
            if handles:
                by_label = dict(zip(labels[:len(present_genotypes)], handles[:len(present_genotypes)]))
                plt.legend(by_label.values(), by_label.keys(), title="Genotype", frameon=True)

            plt.tight_layout()
            plt.savefig(independent_dir / f"{subfolder_name}_{safe_age}_independent.png", dpi=300)
            plt.close()
        except Exception:
            plt.close()

    return outlier_records


# ============================================================
# MAIN EXECUTION
# ============================================================
if __name__ == "__main__":
    print("=" * 60)
    print("MULTI-FEATURE PLOTTING PIPELINE")
    print("=" * 60)

    sheets_dict = load_latest_excel_report(RESULT_DIR)
    all_outliers = []

    # Map features to generate plots for across sheets: (sheet_name, feature_column_prefix, plot_title_label, subfolder_name)
    feature_map = [
        # --- Pax7 Features ---
        ("Pax7_Features", "pax7_cell_density g", "Pax7 Density (cells / um2)", "pax7_density"),
        ("Pax7_Features", "pax7_per_muscle_fibre g", "Pax7 per Muscle Fibre", "pax7_per_muscle_fibre"),
        ("Pax7_Features", "pax7_per_total_nuclei g", "Pax7 per Total Nuclei", "pax7_per_total_nuclei"),
        ("Pax7_Features", "voronoi_cell_area_mean g", "Mean Voronoi Cell Area (um2)", "voronoi_cell_area_mean"),
        ("Pax7_Features", "voronoi_cell_area_std g", "SD of Voronoi Cell Area (um2)", "voronoi_cell_area_std"),
        
        # --- Muscle Fibres ---
        ("Muscle_Fibres", "muscle_fibre_density g", "Muscle Fibre Density (fibres / um2)", "muscle_fibre_density"),
        ("Muscle_Fibres", "muscle_fibre_area_mean g", "Mean Muscle Fibre Area (um2)", "muscle_fibre_area_mean"),
        ("Muscle_Fibres", "muscle_fibre_area_std g", "SD of Muscle Fibre Area (um2)", "muscle_fibre_area_std"),

        # --- Nuclei ---
        ("Nuclei", "nuclei_density g", "Nuclei Density (nuclei / um2)", "nuclei_density"),
        ("Nuclei", "nuclei_per_muscle_fibre g", "Nuclei per Muscle Fibre", "nuclei_per_muscle_fibre"),
    ]

    for sheet_name, prefix, label, subfolder in feature_map:
        if sheet_name in sheets_dict:
            df_sheet = sheets_dict[sheet_name]
            outliers = generate_feature_plots(df_sheet, prefix, label, subfolder)
            all_outliers.extend(outliers)

    # ========================================================
    # EXPORT OUTLIER REPORT TO EXCEL (.XLSX)
    # ========================================================
    outlier_excel_path = RESULT_DIR / "detected_outliers_report.xlsx"
    if all_outliers:
        df_outliers = pd.DataFrame(all_outliers)
        with pd.ExcelWriter(outlier_excel_path, engine="openpyxl") as writer:
            df_outliers.to_excel(writer, sheet_name="Outliers", index=False)
        print(f"\n[REPORT] Detected {len(df_outliers)} outlier data points.")
        print(f"-> Outlier tracking report saved to Excel:\n   {outlier_excel_path}")
    else:
        with pd.ExcelWriter(outlier_excel_path, engine="openpyxl") as writer:
            pd.DataFrame(columns=["metric", "image", "age", "genotype", "muscle_group", "value"]).to_excel(writer, sheet_name="Outliers", index=False)
        print(f"\n[REPORT] No outliers exceeding 2*SD were detected. Blank report saved to:\n   {outlier_excel_path}")

    print("\n" + "=" * 60)
    print("[SUCCESS] ALL FEATURE PLOTS AND EXCEL OUTLIER REPORTS GENERATED")
    print(f"Results located in: {RESULT_DIR / 'plots_density_by_age'}")
    print("=" * 60)