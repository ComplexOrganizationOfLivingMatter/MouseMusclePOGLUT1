import numpy as np
import pandas as pd
import cv2
import matplotlib.pyplot as plt
from pathlib import Path
from skimage import io, measure

# --- CONFIGURATION ---
AGES = ['p4', 'p7', 'p10', 'p15', 'p20', 'p35']

UM_PER_PIXEL = 0.5055443
PIXEL_AREA_UM2 = UM_PER_PIXEL ** 2

PATH_VORONOI = Path('../results/modelV2/voronoiImages')
PATH_MASKS = Path('../data/annotatedMuscleMasks')
RESULT_DIR = Path('../results/modelV2/heatmaps_comparison')
RESULT_DIR.mkdir(parents=True, exist_ok=True)

# Separate heatmaps folder
SEPARATE_HEATMAPS_DIR = RESULT_DIR / "separate_heatmaps"
SEPARATE_HEATMAPS_DIR.mkdir(parents=True, exist_ok=True)

CMAP = 'RdBu_r'
V_LIMIT = 0.6

def get_wt_averages_um2():
    avg_areas_um2 = {}
    print("--- Calculating WT Baselines (um2) ---")
    for age in AGES:
        age_folder = PATH_VORONOI / age
        if not age_folder.exists(): continue
        wt_files = [f for f in age_folder.glob('*.tif') if 'KI' not in f.name.upper() and 'POGLU' not in f.name.upper()]
        all_areas_px = []
        for f in wt_files:
            img = io.imread(str(f))
            props = measure.regionprops(img)
            all_areas_px.extend([p.area for p in props if p.label > 0])
        avg_areas_um2[age] = (np.mean(all_areas_px) * PIXEL_AREA_UM2) if all_areas_px else 1.0
        print(f"  {age} WT Mean: {avg_areas_um2[age]:.2f} um2")
    return avg_areas_um2

def extract_detailed_stats(wt_stats):
    """Extracts Detailed Stats per Muscle Region/Image.
       Calculates Mean and Std Area per muscle group per image.
       The baseline (avg_wt) is in UM2 from wt_stats.
    """
    detailed_data = []
    print("\n--- Extracting Detailed Voronoi Stats per Zone ---")

    for age in AGES:
        age_folder = PATH_VORONOI / age
        if not age_folder.exists(): continue
        avg_wt_um2 = wt_stats.get(age, 1.0)
        # Convert baseline from UM2 back to PIXELS to calculate log-ratios during extraction
        avg_wt_px = avg_wt_um2 / PIXEL_AREA_UM2

        # Create age folders within the separate heatmaps directory
        age_heatmap_dir = SEPARATE_HEATMAPS_DIR / age
        age_heatmap_dir.mkdir(parents=True, exist_ok=True)

        for vor_path in age_folder.glob('*.tif'):
            img_vor = io.imread(str(vor_path))
            mask_path = PATH_MASKS / age / vor_path.name
            mask_zones = io.imread(str(mask_path)) if mask_path.exists() else None
            if mask_zones is not None and mask_zones.ndim == 3: mask_zones = mask_zones[:,:,0]

            if mask_zones is not None and mask_zones.shape != img_vor.shape:
                mask_zones = cv2.resize(mask_zones, (img_vor.shape[1], img_vor.shape[0]), interpolation=cv2.INTER_NEAREST)

            row = {
                'image': vor_path.name,
                'age': age,
                'genotype': 'KIKO' if 'KI' in vor_path.name.upper() else 'WT'
            }

            props = measure.regionprops(img_vor)
            label_to_area = {p.label: p.area * PIXEL_AREA_UM2 for p in props}

            # Lookup array for heatmap log-ratio: log10(Area/avg_wt_px)
            max_label = img_vor.max()
            lookup = np.zeros(max_label + 1)
            for p in props:
                val = np.log10(p.area / avg_wt_px)
                lookup[p.label] = np.clip(val, -V_LIMIT, V_LIMIT)
            heatmap = lookup[img_vor]

            # Detailed stats calculation per zone
            for z in range(1, 16):
                if z == 6: continue
                if mask_zones is not None:
                    zone_ids = np.unique(img_vor[mask_zones == z])
                    zone_ids = zone_ids[zone_ids > 0]
                    areas = [label_to_area[lid] for lid in zone_ids if lid in label_to_area]
                    if areas:
                        row[f'mean_area_g{z}_um2'] = np.mean(areas)
                        row[f'std_area_g{z}_um2'] = np.std(areas)
                    else:
                        row[f'mean_area_g{z}_um2'] = np.nan
                        row[f'std_area_g{z}_um2'] = np.nan

            detailed_data.append(row)

            # --- SEPARATE HEATMAP SAVING ---
            heatmap_display = np.copy(heatmap)
            if mask_zones is not None:
                heatmap_display[(mask_zones == 0) | (mask_zones == 6)] = np.nan

            fig_sep, ax_sep = plt.subplots(figsize=(10, 10), facecolor='black')
            im_sep = ax_sep.imshow(heatmap_display, cmap=CMAP, vmin=-V_LIMIT, vmax=V_LIMIT)

            if mask_zones is not None:
                mask_bin = (mask_zones > 0).astype(np.uint8)
                contours, _ = cv2.findContours(mask_bin, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
                for cnt in contours:
                    poly = cnt.reshape(-1, 2)
                    ax_sep.plot(poly[:,0], poly[:,1], color='white', lw=0.5, alpha=0.5)

            ax_sep.set_title(f"Organization Map: {vor_path.stem}\n(Age: {age.upper()})", color='white', fontsize=12)
            ax_sep.axis('off')

            cbar_ax_sep = fig_sep.add_axes([0.3, 0.08, 0.4, 0.02])
            cb_sep = fig_sep.colorbar(im_sep, cax=cbar_ax_sep, orientation='horizontal')
            cb_sep.set_label(r'$\log_{10}(Area / Mean_{WT})$', color='white')
            cb_sep.ax.xaxis.set_tick_params(color='white', labelcolor='white')

            sep_heatmap_path = age_heatmap_dir / f"{vor_path.stem}_heatmap.png"
            plt.savefig(sep_heatmap_path, facecolor='black', bbox_inches='tight', dpi=150)
            plt.close(fig_sep)

    df_details = pd.DataFrame(detailed_data)
    csv_path = RESULT_DIR / "Detailed_Voronoi_Stats_um2.csv"
    df_details.to_csv(csv_path, index=False)
    print(f"  [SAVED] Detailed stats to: {csv_path.name}")
    print(f"  [SAVED] Separate heatmaps to age folders in {SEPARATE_HEATMAPS_DIR.absolute()}")
    return df_details

def create_comparison_montage(age, avg_wt_um2):
    age_folder = PATH_VORONOI / age
    wt_files = sorted([f for f in age_folder.glob('*.tif') if 'KI' not in f.name.upper()])
    ki_files = sorted([f for f in age_folder.glob('*.tif') if 'KI' in f.name.upper()])
    
    n_cols = max(len(wt_files), len(ki_files))
    if n_cols == 0: return

    # Baseline conversion
    avg_wt_px = avg_wt_um2 / PIXEL_AREA_UM2

    # CRITICAL: Reduce figure size and DPI to save memory
    fig, axes = plt.subplots(2, n_cols, figsize=(n_cols * 3, 6), facecolor='black', dpi=100)
    if n_cols == 1: axes = axes.reshape(2, 1)

    last_im = None
    for row_idx, (files, label) in enumerate(zip([wt_files, ki_files], ["WT", "KIKO"])):
        for col_idx in range(n_cols):
            ax = axes[row_idx, col_idx]
            if col_idx < len(files):
                f_path = files[col_idx]
                m_path = PATH_MASKS / age / f_path.name
                
                h_map, m_zones = process_heatmap_data(f_path, m_path, avg_wt_px)
                
                # --- MEMORY PROTECTION: DOWNSAMPLE FOR DISPLAY ---
                # If image is huge, resize it to max 2000px width for the montage
                max_dim = 2000
                h, w = h_map.shape
                if w > max_dim or h > max_dim:
                    scale = max_dim / max(h, w)
                    new_size = (int(w * scale), int(h * scale))
                    # Use NEAREST to preserve the Voronoi logic
                    h_map = cv2.resize(h_map, new_size, interpolation=cv2.INTER_NEAREST)
                    if m_zones is not None:
                        m_zones = cv2.resize(m_zones, new_size, interpolation=cv2.INTER_NEAREST)

                last_im = ax.imshow(h_map, cmap=CMAP, vmin=-V_LIMIT, vmax=V_LIMIT)
                
                # Plot Boundaries (thinner for downsampled view)
                if m_zones is not None:
                    mask_bin = (m_zones > 0).astype(np.uint8)
                    contours, _ = cv2.findContours(mask_bin, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
                    for cnt in contours:
                        poly = cnt.reshape(-1, 2)
                        ax.plot(poly[:,0], poly[:,1], color='white', lw=0.3, alpha=0.5)
                
                ax.set_title(f"{label}", color='white', fontsize=8)
            ax.axis('off')

    plt.suptitle(f"Comparison: {age.upper()}", color='white', fontsize=14, y=0.98)
    
    if last_im:
        cbar_ax = fig.add_axes([0.3, 0.08, 0.4, 0.02])
        cb = fig.colorbar(last_im, cax=cbar_ax, orientation='horizontal')
        cb.ax.xaxis.set_tick_params(color='white', labelcolor='white', labelsize=8)

    save_name = RESULT_DIR / f"Montage_{age}_Comparison.png"
    # Use explicit clear to free memory immediately after saving
    plt.savefig(save_name, facecolor='black', bbox_inches='tight', dpi=150)
    plt.close(fig)
    plt.clf() # Clear current figure
    del h_map # Force delete large array from memory
    print(f"  [SAVED] Montage for {age}")
    
def process_heatmap_data(vor_path, mask_path, avg_wt_px):
    """
    Loads Voronoi and Mask images, calculates log-normalized areas,
    and returns the heatmap array and the mask zones.
    """
    # 1. Load images
    img_vor = io.imread(str(vor_path))
    mask_zones = io.imread(str(mask_path)) if mask_path.exists() else None
    
    # Handle multi-channel masks if necessary
    if mask_zones is not None and mask_zones.ndim == 3:
        mask_zones = mask_zones[:, :, 0]

    # 2. Calculate properties and map to log-space
    props = measure.regionprops(img_vor)
    max_label = img_vor.max()
    lookup = np.zeros(max_label + 1)
    
    for p in props:
        # log10(Area / WT_Mean)
        val = np.log10(p.area / avg_wt_px)
        lookup[p.label] = np.clip(val, -V_LIMIT, V_LIMIT)
    
    heatmap = lookup[img_vor]
    
    # 3. Apply Anatomical Masking (Background=0, Bone=6)
    if mask_zones is not None:
        if mask_zones.shape != heatmap.shape:
            mask_zones = cv2.resize(mask_zones, (heatmap.shape[1], heatmap.shape[0]), 
                                    interpolation=cv2.INTER_NEAREST)
        
        # Set non-muscle areas to NaN so they appear black/transparent in the plot
        heatmap[(mask_zones == 0) | (mask_zones == 6)] = np.nan
        
    return heatmap, mask_zones
    

if __name__ == "__main__":
    wt_stats = get_wt_averages_um2()
    
    # Detailed stats will also create separate heatmaps in age folders
    detailed_vor_df = extract_detailed_stats(wt_stats)
    
    # Still create the main comparison montages
    for age in AGES:
        if age in wt_stats:
            create_comparison_montage(age, wt_stats[age])