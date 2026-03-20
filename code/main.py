import numpy as np
import pandas as pd
import cv2
import tifffile
from scipy.stats import entropy
from skimage import io, measure, filters, morphology
from scipy.ndimage import distance_transform_edt
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from datetime import datetime

# --- BIOLOGICAL CONSTANTS ---
BONE_ZONE = 6            
MAX_MUSCLE_ZONE = 15     
MIN_AREA_UM2 = 10        
MAX_NUCLEI_AREA_UM2 = 70 
UM_PER_PIXEL = 0.5055443 

def process_single_image(mask_path, path_preds, output_mask_root, output_vor_root):
    img_name = mask_path.name
    age_folder = mask_path.parent.name 
    genotype = "KIKO" if "KI" in img_name.upper() else "WT"
    
    mask_file_path = output_mask_root / age_folder / img_name
    vor_file_path = output_vor_root / age_folder / img_name
    mask_file_path.parent.mkdir(parents=True, exist_ok=True)
    vor_file_path.parent.mkdir(parents=True, exist_ok=True)
    
    pixel_area_um2 = UM_PER_PIXEL ** 2

    try:
        
        # CHECK IF PREVIOUSLY PROCESSED
        if vor_file_path.exists() and mask_file_path.exists():
                try:
                    # 1. Load the 32-bit Labeled Voronoi Map
                    global_voronoi = io.imread(str(vor_file_path))
                    
                    # 2. Load the Binary Nuclei Mask (The one that looks white)
                    binary_nuclei = io.imread(str(mask_file_path)) > 0
                    
                    # 3. THE MAGIC STEP: 
                    # Multiply the labeled territories by the binary nuclei positions.
                    # This 'crops' the Voronoi labels back down to just the nuclei shapes.
                    final_nuclei_mask = (global_voronoi * binary_nuclei).astype(np.uint32)
                    
                    # 4. Verify the count matches the unique IDs in the Voronoi
                    # (Subtract 1 to ignore the background/0)
                    actual_count = len(np.unique(global_voronoi)) - 1
                    print(f"  [CACHE] Recovered {actual_count} nuclei for: {img_name}")
                    
                    # Load mask_zones for the rest of the extraction
                    mask_zones_raw = io.imread(str(mask_path))
                    mask_zones = mask_zones_raw[:,:,0] if mask_zones_raw.ndim == 3 else mask_zones_raw
                    
                    goto_extraction = True 
                except Exception as e:
                    print(f"  [WARN] Reload failed, reprocessing... {e}")
                    goto_extraction = False
            
        if not goto_extraction:
            mask_zones_raw = io.imread(str(mask_path))
            mask_zones = mask_zones_raw[:,:,0] if mask_zones_raw.ndim == 3 else mask_zones_raw
            h, w = mask_zones.shape
            
            # --- STAGE 1: LOAD & NORMALIZE PREDICTION ---
            pred_path = next((f for f in path_preds if f.name == img_name), None)
            if not pred_path: return None
            
            img_pred_raw = io.imread(str(pred_path))
            if img_pred_raw.ndim == 3: img_pred_raw = img_pred_raw[:,:,0]
            
            # Force to 8-bit scale regardless of input depth
            if img_pred_raw.dtype == np.float32 or img_pred_raw.dtype == np.float64:
                if img_pred_raw.max() <= 1.01:
                    img_pred = (img_pred_raw * 255).astype(np.uint8)
                else:
                    img_pred = img_pred_raw.astype(np.uint8)
            else:
                img_pred = img_pred_raw.astype(np.uint8)
    
            if img_pred.shape != (h, w):
                img_pred = cv2.resize(img_pred, (w, h), interpolation=cv2.INTER_NEAREST)
    
            # --- STAGE 2: BINARIZATION & FILTERING ---
            if "modelv2" in str(pred_path).lower():
                # Use a high-confidence threshold for Model V2
                bw = (img_pred >= 250).astype(np.uint8)
            else:
                # Otsu for Model V1
                try:
                    thresh = filters.threshold_otsu(img_pred)
                    bw = (img_pred > thresh).astype(np.uint8)
                except ValueError: # If image is uniform
                    bw = np.zeros_like(img_pred, dtype=np.uint8)
    
            # Morphological clean and size filter
            labels = measure.label(morphology.closing(bw, morphology.disk(1)))
            final_nuclei_mask = morphology.remove_small_objects(
                labels, min_size=MIN_AREA_UM2/pixel_area_um2
            ).astype(np.uint32)
            
            # Exclude Bone and Outside regions
            final_nuclei_mask[(mask_zones == 0) | (mask_zones == BONE_ZONE)] = 0
            
            # Save segmented nuclei
            tifffile.imwrite(str(mask_file_path), (final_nuclei_mask > 0).astype(np.uint8)*255)
    
            # --- STAGE 2.5: BOUNDING-BOX LOCAL VORONOI ---
            global_voronoi = np.zeros((h, w), dtype=np.uint32)
            
            # Pre-group nuclei by zone to avoid searching the full mask 15 times
            all_props = measure.regionprops(final_nuclei_mask)
            zone_to_nuclei = {z: [] for z in range(1, MAX_MUSCLE_ZONE + 1)}
            for p in all_props:
                cy, cx = map(int, p.centroid)
                z_id = mask_zones[cy, cx]
                if 0 < z_id <= MAX_MUSCLE_ZONE and z_id != BONE_ZONE:
                    zone_to_nuclei[z_id].append(p)
    
            for z, nuclei_list in zone_to_nuclei.items():
                if not nuclei_list: continue
                
                # 1. Get the mask for this specific muscle
                z_mask = (mask_zones == z)
                
                # 2. Find the BOUNDING BOX (The smallest rectangle containing the muscle)
                coords = np.argwhere(z_mask)
                y_min, x_min = coords.min(axis=0)
                y_max, x_max = coords.max(axis=0)
                
                # 3. CROP the workspace to this box
                # Adding +1 to max because slicing is exclusive at the end
                z_crop = z_mask[y_min:y_max+1, x_min:x_max+1]
                
                # 4. Create SEEDS only for nuclei belonging to this zone inside the crop
                seed_crop = np.zeros(z_crop.shape, dtype=np.uint32)
                for p in nuclei_list:
                    # Get nucleus pixels relative to the bounding box origin
                    p_mask = (final_nuclei_mask == p.label)
                    # Slice the global nucleus mask to fit our local crop
                    p_crop = p_mask[y_min:y_max+1, x_min:x_max+1]
                    seed_crop[p_crop] = p.label
    
                # 5. LOCAL VORONOI: Grow the seeds within the crop
                # 'indices' finds the local coordinates of the nearest nucleus
                _, indices = distance_transform_edt(seed_crop == 0, return_indices=True)
                local_vor = seed_crop[indices[0], indices[1]]
                
                # 6. CLIP & PASTE: Only keep pixels inside the actual muscle shape
                # Paste the result back into the global 32-bit image
                global_voronoi[y_min:y_max+1, x_min:x_max+1][z_crop] = local_vor[z_crop]
    
            tifffile.imwrite(str(vor_file_path), global_voronoi, compression='zlib')

        # --- STAGE 4: EXTRACTION ---
        props = measure.regionprops(final_nuclei_mask)
        if not props: return {"image": img_name, "age": age_folder, "genotype": genotype, "status": "no_nuclei"}

        centroids = np.array([p.centroid for p in props]).astype(int)
        areas_um2 = np.array([p.area * pixel_area_um2 for p in props])
        nuclei_counts = np.ceil(areas_um2 / MAX_NUCLEI_AREA_UM2)
        zone_ids = mask_zones[centroids[:, 0], centroids[:, 1]]
        
        res = {"image": img_name, "age": age_folder, "genotype": genotype}
        total_a, total_n = 0, 0

        for z in range(1, MAX_MUSCLE_ZONE + 1):
            if z == BONE_ZONE: continue
            z_m = (mask_zones == z)
            if np.any(z_m):
                z_area = np.sum(z_m) * pixel_area_um2
                z_n_count = np.sum(nuclei_counts[zone_ids == z])
                total_a += z_area
                total_n += z_n_count
                
                res[f'area g{z} (um2)'] = z_area
                res[f'nuclei/um2 g{z}'] = z_n_count / z_area if z_area > 0 else 0
                
                # We need at least 2 cells to calculate a meaningful distribution
                v_cells_z = global_voronoi[z_m]
                u_ids, counts = np.unique(v_cells_z[v_cells_z > 0], return_counts=True)
                
                # S = Richness (Number of unique Voronoi territories)
                S = len(u_ids) 
                
                if S >= 3:
                    # 1. Calculate Probabilities
                    probs = counts / np.sum(counts)
                    
                    # 2. Raw Shannon Entropy (H) - using base e
                    raw_ent = entropy(probs)
                    res[f'entropy g{z}'] = raw_ent
                    
                    # 3. Normalized Entropy (Pielou's Evenness J')
                    # J' = H / ln(S)
                    # This scales the result between 0 (perfect order) and 1 (max disorder)
                    norm_ent = raw_ent / np.log(S)
                    res[f'normalized_entropy g{z}'] = np.clip(norm_ent, 0, 1)
                else:
                    # Not enough cells to calculate a meaningful distribution
                    res[f'entropy g{z}'] = np.nan
                    res[f'normalized_entropy g{z}'] = np.nan

        res['total nuclei detected'] = total_n
        res['whole muscle area (um2)'] = total_a
        res['whole muscle nuclei/um2'] = total_n / total_a if total_a > 0 else 0
        return res

    except Exception as e:
        print(f"  [ERR] {img_name}: {e}")
        return None

# (Keep the __main__ block as is from the previous version)
# --- EXECUTION ---
if __name__ == "__main__":
    BASE_DIR = Path('..')
    PATH_MASKS = list((BASE_DIR / 'data/annotatedMuscleMasks').rglob('*.tif'))
    PATH_PREDS = list((BASE_DIR / 'data/predictions/all_predictions_18032026_modelv2/input20260318_154449/results/input20260318_154449_1/per_image').rglob('*.tif'))
    
    RESULT_DIR = BASE_DIR / 'results' / 'modelV2'
    QC_ROOT = RESULT_DIR / 'binaryNucleiStemCells_filtered'
    VOR_ROOT = RESULT_DIR / 'voronoiImages'
    for d in [QC_ROOT, VOR_ROOT]: d.mkdir(parents=True, exist_ok=True)

    print(f"\n{'='*50}\nSTARTING PIPELINE: {len(PATH_MASKS)} Files Found\n{'='*50}")

    results_list = []
    with ProcessPoolExecutor() as executor:
        futures = {executor.submit(process_single_image, mp, PATH_PREDS, QC_ROOT, VOR_ROOT): mp for mp in PATH_MASKS}
        for i, future in enumerate(futures):
            res = future.result()
            if res: results_list.append(res)
            if (i + 1) % 5 == 0 or (i+1) == len(PATH_MASKS):
                print(f"--- Global Progress: {i+1}/{len(PATH_MASKS)} complete ---")

    print("\n[STAGE] Reordering Columns and Deduplicating...")
    df = pd.DataFrame(results_list)
    df = df.loc[:, ~df.columns.duplicated()]

    base_cols = ['image', 'age', 'genotype', 'total nuclei detected', 'whole muscle area (um2)', 'whole muscle nuclei/um2']
    metrics = ['area', 'nuclei/um2', 'normalized_entropy']
    group_cols = []
    for z in range(1, 16):
        if z == BONE_ZONE: continue
        for m in metrics:
            matches = [c for c in df.columns if (c.startswith(m) and (c.endswith(f' g{z}') or c.endswith(f' g{z} (um2)')))]
            group_cols.extend(matches)
    
    df = df.reindex(columns=base_cols + [c for c in group_cols if c in df.columns])

    print("[STAGE] Generating Summary Tables...")
    grouped = df.groupby(['genotype', 'age'], sort=False)
    summary_mean = grouped.mean(numeric_only=True)
    summary_std = grouped.std(numeric_only=True)

    summary_final = pd.DataFrame(index=summary_mean.index)
    for col in df.columns:
        if col in summary_mean.columns:
            summary_final[f"{col}_mean"] = summary_mean[col]
            summary_final[f"{col}_std"] = summary_std[col]

    timestamp = datetime.now().strftime("%Y%m%d_%H%M")
    output_path = RESULT_DIR / f"Pax7_Genotype_Analysis_{timestamp}.xlsx"
    with pd.ExcelWriter(output_path, engine='openpyxl') as writer:
        df.to_excel(writer, sheet_name='Raw_Data', index=False)
        summary_final.to_excel(writer, sheet_name='Summary_Stats')

    print(f"\n{'='*50}\n[SUCCESS] Final Report Saved: {output_path.name}\n{'='*50}")