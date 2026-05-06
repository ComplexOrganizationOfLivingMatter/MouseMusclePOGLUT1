import numpy as np
import pandas as pd
import cv2
import matplotlib.pyplot as plt
from pathlib import Path
from skimage import io, measure
import matplotlib.image as mpimg
import gc

# --- CONFIGURATION ---
# Ages to process based on folder names
AGES = ['p0', 'p4', 'p7', 'p10', 'p15', 'p20', 'p35']
UM_PER_PIXEL = 0.5055443
PIXEL_AREA_UM2 = UM_PER_PIXEL ** 2
MODEL_FOLDER = 'modelV2'

# Path Setup
PATH_VORONOI = Path(f"../../results/{MODEL_FOLDER}/voronoiImages")
PATH_MASKS = Path('../../data/annotatedMuscleMasks')
RESULT_DIR = Path(f"../../results/{MODEL_FOLDER}/heatmaps_comparison")
RESULT_DIR.mkdir(parents=True, exist_ok=True)
SEPARATE_HEATMAPS_DIR = RESULT_DIR / "separate_heatmaps"
SEPARATE_HEATMAPS_DIR.mkdir(parents=True, exist_ok=True)

# Visualization Settings
CMAP = 'RdBu_r'  # Blue: Small/Dense (High log ratio) | Red: Big/Sparse (Low log ratio)
V_LIMIT = 0.6    # Limits for the log10 color scale
RENDER_MAX_DIM = 2000  # Resolution cap for PNGs to prevent Out Of Memory (OOM) crashes

ONLY_COMPOSITES = True  # Set to True to skip TIF processing and use existing PNGs

# =========================================================
# UTILITIES
# =========================================================

def get_downsampled_label_image(label_img, max_dim=RENDER_MAX_DIM):
    """
    Resizes uint32 label images safely. 
    OpenCV resize doesn't support uint32, so we cast to float32 and use NEAREST interpolation
    to preserve integer label IDs.
    """
    h, w = label_img.shape
    if max(h, w) <= max_dim:
        return label_img
    scale = max_dim / max(h, w)
    new_size = (int(w * scale), int(h * scale))
    
    img_float = label_img.astype(np.float32)
    resized_float = cv2.resize(img_float, new_size, interpolation=cv2.INTER_NEAREST)
    return resized_float.astype(np.uint32)

def create_thick_outlines(label_img):
    """
    Detects boundaries between Voronoi cells. 
    Dilates the edges so they are visible even after downsampling.
    Returns a masked array where only the edges are non-transparent.
    """
    edges = np.zeros(label_img.shape, dtype=bool)
    # Compare each pixel with neighbor below and to the right
    edges[:-1, :] |= (label_img[:-1, :] != label_img[1:, :])
    edges[:, :-1] |= (label_img[:, :-1] != label_img[:, 1:])
    
    kernel = np.ones((2, 2), np.uint8)
    thick_edges = cv2.dilate(edges.astype(np.uint8), kernel)
    return np.ma.masked_where(thick_edges == 0, thick_edges)

def clean_muscle_name(filename, is_ki):
    """
    Formats titles for the plots. 
    Example: 'Muscle_WT_561' -> 'WT 561'
    Example: 'Muscle_KI_119' -> 'KIKO 119'
    """
    parts = filename.replace('_', ' ').split(' ')
    val_id = parts[-1] # Take the last part (usually the ID)
    prefix = "KIKO" if is_ki else "WT"
    return f"{prefix} {val_id}"

# =========================================================
# CORE PROCESSING
# =========================================================

def get_wt_averages_um2():
    """Calculates age-specific mean areas for WT baseline to normalize the heatmaps."""
    avg_areas_um2 = {}
    print("Step 1: Calculating WT Baselines...")
    for age in AGES:
        age_folder = PATH_VORONOI / age
        if not age_folder.exists(): continue
        
        # Filter files that are Wild Type (not KI or POGLUT)
        wt_files = [f for f in age_folder.glob('*.tif') if 'KI' not in f.name.upper() and 'POGLU' not in f.name.upper()]
        all_areas_px = []
        
        for f in wt_files:
            img = io.imread(str(f))
            props = measure.regionprops(img)
            all_areas_px.extend([p.area for p in props if p.label > 0])
            del img, props
            gc.collect() # Force memory cleanup
            
        avg_areas_um2[age] = (np.mean(all_areas_px) * PIXEL_AREA_UM2 if all_areas_px else 1.0)
    return avg_areas_um2

def extract_detailed_stats(wt_stats):
    """Generates individual heatmaps with black backgrounds and cell outlines."""
    detailed_data = []
    print("Step 2: Processing Individual Heatmaps...")
    for age in AGES:
        age_folder = PATH_VORONOI / age
        if not age_folder.exists(): continue
        avg_wt_px = wt_stats.get(age, 1.0) / PIXEL_AREA_UM2
        age_heatmap_dir = SEPARATE_HEATMAPS_DIR / age
        age_heatmap_dir.mkdir(parents=True, exist_ok=True)

        for vor_path in age_folder.glob('*.tif'):
            img_vor = io.imread(str(vor_path))
            
            # --- 1. FULL RESOLUTION CALCULATION (For Stats) ---
            props = measure.regionprops(img_vor)
            lookup = np.zeros(img_vor.max() + 1, dtype=np.float32)
            label_to_area = {}
            for p in props:
                area_um2 = p.area * PIXEL_AREA_UM2
                label_to_area[p.label] = area_um2
                # Calculate normalized log10 ratio
                lookup[p.label] = np.clip(np.log10(p.area / avg_wt_px), -V_LIMIT, V_LIMIT)

            # Muscle Mask logic
            mask_path = PATH_MASKS / age / vor_path.name
            mask_zones = None
            if mask_path.exists():
                mask_zones = io.imread(str(mask_path))
                if mask_zones.ndim == 3: mask_zones = mask_zones[:, :, 0]
                if mask_zones.shape != img_vor.shape:
                    mask_zones = cv2.resize(mask_zones, (img_vor.shape[1], img_vor.shape[0]), interpolation=cv2.INTER_NEAREST)
            
            # Save data to CSV row
            row = {'image': vor_path.name, 'age': age, 'genotype': 'KIKO' if 'KI' in vor_path.name.upper() else 'WT'}
            if mask_zones is not None:
                for z in range(1, 16):
                    if z == 6: continue
                    zone_ids = np.unique(img_vor[mask_zones == z])
                    areas = [label_to_area[lid] for lid in zone_ids if lid in label_to_area]
                    row[f'mean_area_g{z}_um2'] = np.mean(areas) if areas else np.nan
            detailed_data.append(row)

            # --- 2. DOWNSAMPLED PLOTTING (For Visuals) ---
            img_vor_small = get_downsampled_label_image(img_vor)
            heatmap_small = lookup[img_vor_small]
            edge_overlay = create_thick_outlines(img_vor_small)

            # Handle background: set values outside mask to NaN (will show as black)
            if mask_zones is not None:
                mask_small = cv2.resize(mask_zones, (img_vor_small.shape[1], img_vor_small.shape[0]), interpolation=cv2.INTER_NEAREST)
                heatmap_small[(mask_small == 0) | (mask_small == 6)] = np.nan

            # Setup Figure with strictly black background
            fig, ax = plt.subplots(figsize=(10, 10), facecolor='black')
            ax.set_facecolor('black') 
            
            clean_name = clean_muscle_name(vor_path.stem, 'KI' in vor_path.name.upper())
            ax.set_title(clean_name, color='white', fontsize=14)
            
            # Plot Heatmap and Outlines
            ax.imshow(heatmap_small, cmap=CMAP, vmin=-V_LIMIT, vmax=V_LIMIT)
            ax.imshow(edge_overlay, cmap='binary_r', interpolation='nearest', alpha=0.9)

            # Draw white muscle contours
            if mask_zones is not None:
                mask_bin = (mask_small > 0).astype(np.uint8)
                contours, _ = cv2.findContours(mask_bin, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
                for cnt in contours:
                    poly = cnt.reshape(-1, 2)
                    ax.plot(poly[:, 0], poly[:, 1], color='white', lw=1.0, alpha=0.5)

            ax.axis('off')
            plt.savefig(age_heatmap_dir / f"{vor_path.stem}_heatmap.png", facecolor='black', edgecolor='none', bbox_inches='tight', dpi=120)
            plt.close(fig)
            
            del img_vor, img_vor_small, heatmap_small, props, lookup
            gc.collect()

    pd.DataFrame(detailed_data).to_csv(RESULT_DIR / "Detailed_Voronoi_Stats_um2.csv", index=False)


def create_comparison_montage(age, avg_wt_um2):
    """Creates a side-by-side comparison with cleaned names and black background."""
    print(f"Step 3: Creating Montage for {age}...")
    age_folder = PATH_VORONOI / age
    wt_files = sorted([f for f in age_folder.glob('*.tif') if 'KI' not in f.name.upper()])
    ki_files = sorted([f for f in age_folder.glob('*.tif') if 'KI' in f.name.upper()])
    n_cols = max(len(wt_files), len(ki_files))
    if n_cols == 0: return

    avg_wt_px = avg_wt_um2 / PIXEL_AREA_UM2
    
    # 1. FAIL-SAFE: Initialize figure with black facecolor
    fig, axes = plt.subplots(2, n_cols, figsize=(n_cols * 5, 10), facecolor='black')
    
    # 2. FAIL-SAFE: Manually set the figure background color to black
    fig.patch.set_facecolor('black')
    
    if n_cols == 1: axes = axes.reshape(2, 1)

    last_im = None
    for row_idx, files in enumerate([wt_files, ki_files]):
        is_ki = (row_idx == 1)
        for col_idx in range(n_cols):
            ax = axes[row_idx, col_idx]
            
            # 3. FAIL-SAFE: Force every subplot background to black
            ax.set_facecolor('black')
            
            if col_idx < len(files):
                f_path = files[col_idx]
                img_vor = io.imread(str(f_path))
                props = measure.regionprops(img_vor)
                lookup = np.zeros(img_vor.max() + 1, dtype=np.float32)
                for p in props:
                    lookup[p.label] = np.clip(np.log10(p.area / avg_wt_px), -V_LIMIT, V_LIMIT)
                
                img_small = get_downsampled_label_image(img_vor, max_dim=1200)
                h_map_small = lookup[img_small]
                edges = create_thick_outlines(img_small)
                
                # Plotting
                last_im = ax.imshow(h_map_small, cmap=CMAP, vmin=-V_LIMIT, vmax=V_LIMIT)
                ax.imshow(edges, cmap='binary_r', interpolation='nearest', alpha=0.8)
                
                # Title clean up
                ax.set_title(clean_muscle_name(f_path.stem, is_ki), color='white', fontsize=12)
                
                del img_vor, img_small, h_map_small, props, lookup
            else:
                # Ensure even empty slots are black
                ax.set_visible(False) 
            
            ax.axis('off')

    # Add normalized scale bar (Colorbar)
    if last_im:
        # Create a new axis for the colorbar to prevent layout issues
        cax = fig.add_axes([0.3, 0.05, 0.4, 0.02]) 
        cbar = fig.colorbar(last_im, cax=cax, orientation='horizontal')
        cbar.set_label(r'$\log_{10}(\frac{Area}{Mean\_WT})$', color='white', size=14)
        cbar.ax.xaxis.set_tick_params(color='white', labelcolor='white')
        cbar.outline.set_edgecolor('white')

    # 4. FAIL-SAFE: Force save parameters
    plt.savefig(
        RESULT_DIR / f"Montage_{age}_Comparison.png", 
        facecolor=fig.get_facecolor(), # Use the figure's black color
        edgecolor='none', 
        bbox_inches='tight', 
        transparent=False, # Ensure it's not saved as transparent/white
        dpi=150
    )
    plt.close(fig)
    gc.collect()
    


def create_montage_from_existing_pngs(age):
    """
    Loads pre-existing PNG heatmaps and arranges them into a black montage
    with a global title and normalized scale bar.
    """
    age_folder = SEPARATE_HEATMAPS_DIR / age
    if not age_folder.exists():
        print(f"!!! Error: Folder not found: {age_folder}")
        return

    # Find WT and KI files
    all_pngs = list(age_folder.glob('*.png'))
    wt_files = sorted([f for f in all_pngs if 'KI' not in f.name.upper()])
    ki_files = sorted([f for f in all_pngs if 'KI' in f.name.upper()])
    
    n_cols = max(len(wt_files), len(ki_files))
    if n_cols == 0:
        print(f"--- No PNGs found in {age_folder}")
        return

    print(f"Generating Composite for {age} ({n_cols} columns)...")

    # Initialize Figure
    fig, axes = plt.subplots(2, n_cols, figsize=(n_cols * 5, 11), facecolor='black')
    fig.patch.set_facecolor('black')
    
    # Global Title: "p20 Age Comparison"
    fig.suptitle(f"{age.upper()} Age Comparison", color='white', fontsize=22, fontweight='bold', y=0.95)
    
    if n_cols == 1:
        axes = axes.reshape(2, 1)

    for row_idx, files in enumerate([wt_files, ki_files]):
        for col_idx in range(n_cols):
            ax = axes[row_idx, col_idx]
            ax.set_facecolor('black')
            
            if col_idx < len(files):
                file_to_load = files[col_idx]
                try:
                    img = mpimg.imread(str(file_to_load))
                    ax.imshow(img)
                except Exception as e:
                    print(f"Error loading {file_to_load.name}: {e}")
            
            ax.axis('off')

    # --- ADD NORMALIZED SCALE BAR ---
    # We create a "dummy" image mapping to generate the colorbar logic
    # so it matches the original log10(Area/WT) scale
    sm = plt.cm.ScalarMappable(cmap=CMAP, norm=plt.Normalize(vmin=-V_LIMIT, vmax=V_LIMIT))
    sm._A = [] # Required for matplotlib internal logic
    
    # Create colorbar at the bottom
    cax = fig.add_axes([0.3, 0.07, 0.4, 0.02]) # [left, bottom, width, height]
    cbar = fig.colorbar(sm, cax=cax, orientation='horizontal')
    cbar.set_label(r'$\log_{10}(\frac{Area}{Mean\_WT})$', color='white', size=16, labelpad=10)
    cbar.ax.xaxis.set_tick_params(color='white', labelcolor='white')
    cbar.outline.set_edgecolor('white')

    # Save with forced black background
    output_path = RESULT_DIR / f"Montage_{age}_Comparison.png"
    plt.savefig(
        output_path, 
        facecolor='black', 
        edgecolor='none', 
        bbox_inches='tight', 
        transparent=False, 
        dpi=150
    )
    plt.close(fig)
    gc.collect()
    print(f"Saved: {output_path.name}")

# =========================================================
# EXECUTION
# =========================================================

if __name__ == "__main__":
    if ONLY_COMPOSITES:
        print("Mode: Montage generation only (skipping TIF processing).")
        for age in AGES:
            create_montage_from_existing_pngs(age)
    else:
        wt_stats = get_wt_averages_um2()
        extract_detailed_stats(wt_stats)
        for age in AGES:
            if age in wt_stats:
                create_comparison_montage(age, wt_stats[age])