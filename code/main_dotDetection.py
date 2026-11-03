"""
Created on Mon Sep 7 15:46:58 2026
@author: pedro
"""

import cv2
import tifffile
import numpy as np
import pandas as pd
from pathlib import Path
from datetime import datetime
from scipy.ndimage import distance_transform_edt
from skimage import io, measure
from concurrent.futures import ProcessPoolExecutor, as_completed

# ============================================================
# CONSTANTS & CONFIGURATION
# ============================================================
BONE_ZONE = 6
MAX_MUSCLE_ZONE = 15
UM_PER_PIXEL = 0.5055443
PIXEL_AREA_UM2 = UM_PER_PIXEL**2
EXCLUDED_ZONES = {0, BONE_ZONE, 13}
VALID_ZONES = [z for z in range(1, MAX_MUSCLE_ZONE + 1) if z not in EXCLUDED_ZONES]


def load_single_channel(path: Path) -> np.ndarray:
    """Reads an image and extracts the first channel if multi-channel, preserving high bit depths (uint16/uint32)."""
    img = io.imread(str(path))
    if img.ndim == 3:
        img = img[..., 0]
    return img


def process_single_image(
    mask_path: Path, 
    pred_dict: dict, 
    fibres_dict: dict, 
    nuclei_dict: dict, 
    output_mask_root: Path, 
    output_vor_root: Path
) -> dict | None:
    img_name = mask_path.name
    age_folder = mask_path.parent.name
    genotype = "KIKO" if "KI" in img_name.upper() else "WT"

    mask_file_path = output_mask_root / age_folder / img_name
    vor_file_path = output_vor_root / age_folder / img_name

    mask_file_path.parent.mkdir(parents=True, exist_ok=True)
    vor_file_path.parent.mkdir(parents=True, exist_ok=True)

    try:
        # 1. Load Pax7 prediction
        pred_path = pred_dict.get(img_name)
        if pred_path is None:
            print(f"[WARN] Prediction not found: {img_name}")
            return None

        img_pred = load_single_channel(pred_path)
        pred_h, pred_w = img_pred.shape

        mask_zones = load_single_channel(mask_path)

        # Dimension Check & Nearest Neighbor Transformation for Zone Mask
        if mask_zones.shape != (pred_h, pred_w):
            mask_zones = cv2.resize(mask_zones, (pred_w, pred_h), interpolation=cv2.INTER_NEAREST)

        h, w = mask_zones.shape
        
        # Create strict valid mask based on annotated zones (excluding background, bone, and excluded zones)
        valid_mask = np.isin(mask_zones, list(VALID_ZONES))

        use_cache = False

        # --- Cache Check for Pax7 / Voronoi ---
        if mask_file_path.exists() and vor_file_path.exists():
            try:
                binary_nuclei = load_single_channel(mask_file_path) > 0
                global_voronoi = io.imread(str(vor_file_path)).astype(np.uint32)
                
                final_nuclei_mask = measure.label(
                    (global_voronoi * binary_nuclei) > 0, connectivity=1
                ).astype(np.uint32)
                use_cache = True
            except Exception as e:
                print(f"[WARN] Cache failed for {img_name}. Reprocessing. Reason: {e}")
                
        # --- Pax7 Prediction & Voronoi Generation ---
        if not use_cache:
            raw_mask = measure.label(img_pred > 0, connectivity=1)
            raw_mask[~valid_mask] = 0
            final_nuclei_mask = measure.label(raw_mask > 0, connectivity=1).astype(np.uint32)

            tifffile.imwrite(
                str(mask_file_path), (final_nuclei_mask > 0).astype(np.uint8) * 255
            )

            global_voronoi = np.zeros((h, w), dtype=np.uint32)
            all_props = measure.regionprops(final_nuclei_mask)
            zone_to_nuclei = {z: [] for z in VALID_ZONES}

            for p in all_props:
                cy, cx = map(int, p.centroid)
                if 0 <= cy < h and 0 <= cx < w:
                    z_id = mask_zones[cy, cx]
                    if z_id in zone_to_nuclei:
                        zone_to_nuclei[z_id].append(p)

            for z, nuclei_list in zone_to_nuclei.items():
                if not nuclei_list:
                    continue
                z_mask = mask_zones == z
                if not np.any(z_mask):
                    continue

                coords = np.argwhere(z_mask)
                y_min, x_min = coords.min(axis=0)
                y_max, x_max = coords.max(axis=0)

                z_crop = z_mask[y_min : y_max + 1, x_min : x_max + 1]
                seed_crop = np.zeros(z_crop.shape, dtype=np.uint32)

                for p in nuclei_list:
                    cy, cx = map(int, p.centroid)
                    ly, lx = cy - y_min, cx - x_min
                    if 0 <= ly < seed_crop.shape[0] and 0 <= lx < seed_crop.shape[1]:
                        seed_crop[ly, lx] = p.label

                _, indices = distance_transform_edt(seed_crop == 0, return_indices=True)
                local_voronoi = seed_crop[indices[0], indices[1]]

                global_crop = global_voronoi[y_min : y_max + 1, x_min : x_max + 1]
                global_crop[z_crop] = local_voronoi[z_crop]
                global_voronoi[y_min : y_max + 1, x_min : x_max + 1] = global_crop

            tifffile.imwrite(str(vor_file_path), global_voronoi, compression="zlib")

        # ========================================================
        # LOAD & MASK CELLPOSE SEGMENTATIONS (FIBRES & NUCLEI)
        # ========================================================
        fibre_path = fibres_dict.get(img_name)
        nuclei_path = nuclei_dict.get(img_name)

        fibre_props = []
        if fibre_path and Path(fibre_path).exists():
            fibre_img = load_single_channel(Path(fibre_path))
            if fibre_img.shape != (h, w):
                fibre_img = cv2.resize(fibre_img, (w, h), interpolation=cv2.INTER_NEAREST)
            
            fibre_masked = fibre_img.copy()
            fibre_masked[~valid_mask] = 0
            fibre_labeled = measure.label(fibre_masked, connectivity=1)
            fibre_props = measure.regionprops(fibre_labeled)

        total_nuclei_props = []
        if nuclei_path and Path(nuclei_path).exists():
            nuclei_img = load_single_channel(Path(nuclei_path))
            if nuclei_img.shape != (h, w):
                nuclei_img = cv2.resize(nuclei_img, (w, h), interpolation=cv2.INTER_NEAREST)
            
            nuclei_masked = nuclei_img.copy()
            nuclei_masked[~valid_mask] = 0
            # Integer label pass handles recycled distant IDs while preserving touching boundaries
            nuclei_labeled = measure.label(nuclei_masked, connectivity=1)
            total_nuclei_props = measure.regionprops(nuclei_labeled)

        pax7_props = measure.regionprops(final_nuclei_mask)
        
        res_pax7 = {"image": img_name, "age": age_folder, "genotype": genotype}
        res_fibres = {"image": img_name, "age": age_folder, "genotype": genotype}
        res_nuclei = {"image": img_name, "age": age_folder, "genotype": genotype}

        tot_area = 0.0
        tot_pax7 = 0
        tot_fibres = len(fibre_props)
        tot_total_nuclei = len(total_nuclei_props)

        zone_pax7 = {z: 0 for z in VALID_ZONES}
        for p in pax7_props:
            cy, cx = map(int, p.centroid)
            if 0 <= cy < h and 0 <= cx < w:
                z_id = mask_zones[cy, cx]
                if z_id in zone_pax7:
                    zone_pax7[z_id] += 1

        zone_fibres = {z: [] for z in VALID_ZONES}
        for p in fibre_props:
            cy, cx = map(int, p.centroid)
            if 0 <= cy < h and 0 <= cx < w:
                z_id = mask_zones[cy, cx]
                if z_id in zone_fibres:
                    zone_fibres[z_id].append(p.area * PIXEL_AREA_UM2)

        zone_nuclei_count = {z: 0 for z in VALID_ZONES}
        for p in total_nuclei_props:
            cy, cx = map(int, p.centroid)
            if 0 <= cy < h and 0 <= cx < w:
                z_id = mask_zones[cy, cx]
                if z_id in zone_nuclei_count:
                    zone_nuclei_count[z_id] += 1

        # ========================================================
        # PER-ZONE METRIC CALCULATIONS
        # ========================================================
        for z in VALID_ZONES:
            z_mask = mask_zones == z
            if not np.any(z_mask):
                continue

            z_area = np.sum(z_mask) * PIXEL_AREA_UM2
            tot_area += z_area

            # 1. Muscle Fibres Metrics
            f_areas = zone_fibres[z]
            f_count = len(f_areas)
            res_fibres[f"muscle_zone_area g{z} (um2)"] = z_area
            res_fibres[f"muscle_fibre_count g{z}"] = f_count
            res_fibres[f"muscle_fibre_density g{z} (fibres/um2)"] = f_count / z_area if z_area > 0 else 0
            res_fibres[f"muscle_fibre_area_mean g{z} (um2)"] = float(np.mean(f_areas)) if f_count > 0 else 0.0
            res_fibres[f"muscle_fibre_area_std g{z} (um2)"] = float(np.std(f_areas)) if f_count > 0 else 0.0

            # 2. Nuclei Metrics
            n_count = zone_nuclei_count[z]
            res_nuclei[f"muscle_zone_area g{z} (um2)"] = z_area
            res_nuclei[f"nuclei_count g{z}"] = n_count
            res_nuclei[f"nuclei_density g{z} (nuclei/um2)"] = n_count / z_area if z_area > 0 else 0
            res_nuclei[f"_nuclei_per_muscle_fibre g{z}"] = n_count / f_count if f_count > 0 else 0.0

            # 3. Pax7 Features Metrics
            pax7_count = zone_pax7[z]
            tot_pax7 += pax7_count
            res_pax7[f"muscle_zone_area g{z} (um2)"] = z_area
            res_pax7[f"pax7_cell_count g{z}"] = pax7_count
            res_pax7[f"pax7_cell_density g{z} (pax7_cells/um2)"] = pax7_count / z_area if z_area > 0 else 0
            res_pax7[f"pax7_per_muscle_fibre g{z}"] = pax7_count / f_count if f_count > 0 else 0.0
            res_pax7[f"pax7_per_nuclei g{z}"] = pax7_count / n_count if n_count > 0 else 0.0

            # Voronoi cell size statistics for Pax7
            voronoi_cells = global_voronoi[z_mask]
            voronoi_cells = voronoi_cells[voronoi_cells > 0]
            if voronoi_cells.size > 0:
                _, territory_sizes = np.unique(voronoi_cells, return_counts=True)
                cell_areas = territory_sizes * PIXEL_AREA_UM2
                res_pax7[f"voronoi_cell_area_mean g{z} (um2)"] = float(np.mean(cell_areas))
                res_pax7[f"voronoi_cell_area_std g{z} (um2)"] = float(np.std(cell_areas))
            else:
                res_pax7[f"voronoi_cell_area_mean g{z} (um2)"] = np.nan
                res_pax7[f"voronoi_cell_area_std g{z} (um2)"] = np.nan

        # Whole muscle summary values
        all_fibre_areas = [area for areas in zone_fibres.values() for area in areas]
        
        res_fibres["total_muscle_fibres"] = tot_fibres
        res_fibres["whole_muscle_area (um2)"] = tot_area
        res_fibres["whole_muscle_fibre_density (fibres/um2)"] = tot_fibres / tot_area if tot_area > 0 else 0
        res_fibres["whole_muscle_fibre_area_mean (um2)"] = float(np.mean(all_fibre_areas)) if all_fibre_areas else 0.0
        res_fibres["whole_muscle_fibre_area_std (um2)"] = float(np.std(all_fibre_areas)) if all_fibre_areas else 0.0

        res_nuclei["total_nuclei_count"] = tot_total_nuclei
        res_nuclei["whole_muscle_area (um2)"] = tot_area
        res_nuclei["whole_muscle_nuclei_density (nuclei/um2)"] = tot_total_nuclei / tot_area if tot_area > 0 else 0
        res_nuclei["whole_muscle_nuclei_per_muscle_fibre"] = tot_total_nuclei / tot_fibres if tot_fibres > 0 else 0.0

        res_pax7["total_pax7_cells_detected"] = tot_pax7
        res_pax7["whole_muscle_area (um2)"] = tot_area
        res_pax7["whole_muscle_pax7_density (pax7_cells/um2)"] = tot_pax7 / tot_area if tot_area > 0 else 0
        res_pax7["whole_muscle_pax7_per_muscle_fibre"] = tot_pax7 / tot_fibres if tot_fibres > 0 else 0.0
        res_pax7["whole_muscle_pax7_per_total_nuclei"] = tot_pax7 / tot_total_nuclei if tot_total_nuclei > 0 else 0.0

        return {"pax7": res_pax7, "fibres": res_fibres, "nuclei": res_nuclei}

    except Exception as e:
        print(f"[ERROR] {img_name}: {e}")
        return None


# ============================================================
# MAIN EXECUTION
# ============================================================
if __name__ == "__main__":
    BASE_DIR = Path("..")
    PATH_MASKS = list((BASE_DIR / "data" / "annotatedMuscleMasks").rglob("*.tif"))
    PATH_PREDS = list((BASE_DIR / "data" / "predictions" / "all_largeModel_detection" / "per_image_local_max_check_post_processing").rglob("*.tif"))
    
    # Cellpose Segmentation Paths (Relative and robust)
    SEG_BASE = BASE_DIR / "results" / "cellposeSegmentation"
    PATH_FIBRES = list((SEG_BASE / "muscleFibres").rglob("*.tif"))
    PATH_NUCLEI = list((SEG_BASE / "nuclei").rglob("*.tif"))

    PRED_DICT = {p.name: p for p in PATH_PREDS}
    FIBRES_DICT = {p.name: p for p in PATH_FIBRES}
    NUCLEI_DICT = {p.name: p for p in PATH_NUCLEI}

    RESULT_DIR = BASE_DIR / "results" / "largeModel_dotDetection"
    QC_ROOT = RESULT_DIR / "binaryNucleiDots"
    VOR_ROOT = RESULT_DIR / "voronoiImages"

    QC_ROOT.mkdir(parents=True, exist_ok=True)
    VOR_ROOT.mkdir(parents=True, exist_ok=True)

    print("\n" + "=" * 60)
    print("STARTING MULTI-FEATURE PIPELINE (INTUITIVE FEATURE NAMES)")
    print(f"Muscle masks found: {len(PATH_MASKS)}")
    print(f"BiaPy predictions found: {len(PATH_PREDS)}")
    print(f"Muscle fibre segmentations found: {len(PATH_FIBRES)}")
    print(f"Nuclei segmentations found: {len(PATH_NUCLEI)}")
    print("=" * 60 + "\n")

    pax7_list, fibres_list, nuclei_list = [], [], []

    with ProcessPoolExecutor() as executor:
        futures = {
            executor.submit(
                process_single_image, 
                mask_path, 
                PRED_DICT, 
                FIBRES_DICT, 
                NUCLEI_DICT, 
                QC_ROOT, 
                VOR_ROOT
            ): mask_path
            for mask_path in PATH_MASKS
        }

        for i, future in enumerate(as_completed(futures), start=1):
            res = future.result()
            if res is not None:
                pax7_list.append(res["pax7"])
                fibres_list.append(res["fibres"])
                nuclei_list.append(res["nuclei"])

            if i % 5 == 0 or i == len(PATH_MASKS):
                print(f"[PROGRESS] {i}/{len(PATH_MASKS)} images complete")

    # --- Build DataFrames ---
    print("\n[STAGE] Creating dataframes...")
    df_pax7 = pd.DataFrame(pax7_list)
    df_fibres = pd.DataFrame(fibres_list)
    df_nuclei = pd.DataFrame(nuclei_list)

    if df_pax7.empty:
        raise RuntimeError("No images were successfully processed.")

    def reorder_df(df, base_cols, metric_templates):
        df = df.loc[:, ~df.columns.duplicated()]
        group_cols = [
            col
            for z in VALID_ZONES
            for tmpl in metric_templates
            for col in [tmpl.format(z=z)]
            if col in df.columns
        ]
        return df.reindex(columns=base_cols + group_cols)

    df_pax7 = reorder_df(
        df_pax7,
        [
            "image", "age", "genotype", 
            "total_pax7_cells_detected", "whole_muscle_area (um2)", 
            "whole_muscle_pax7_density (pax7_cells/um2)", 
            "whole_muscle_pax7_per_muscle_fibre", 
            "whole_muscle_pax7_per_nuclei"
        ],
        [
            "muscle_zone_area g{z} (um2)", 
            "pax7_cell_count g{z}", 
            "pax7_cell_density g{z} (pax7_cells/um2)", 
            "pax7_per_muscle_fibre g{z}", 
            "pax7_per_nuclei g{z}", 
            "voronoi_cell_area_mean g{z} (um2)", 
            "voronoi_cell_area_std g{z} (um2)"
        ]
    )

    df_fibres = reorder_df(
        df_fibres,
        [
            "image", "age", "genotype", 
            "total_muscle_fibres", "whole_muscle_area (um2)", 
            "whole_muscle_fibre_density (fibres/um2)", 
            "whole_muscle_fibre_area_mean (um2)", 
            "whole_muscle_fibre_area_std (um2)"
        ],
        [
            "muscle_zone_area g{z} (um2)", 
            "muscle_fibre_count g{z}", 
            "muscle_fibre_density g{z} (fibres/um2)", 
            "muscle_fibre_area_mean g{z} (um2)", 
            "muscle_fibre_area_std g{z} (um2)"
        ]
    )

    df_nuclei = reorder_df(
        df_nuclei,
        [
            "image", "age", "genotype", 
            "total_nuclei_count", "whole_muscle_area (um2)", 
            "whole_muscle_nuclei_density (nuclei/um2)", 
            "whole_muscle_nuclei_per_muscle_fibre"
        ],
        [
            "muscle_zone_area g{z} (um2)", 
            "nuclei_count g{z}", 
            "nuclei_density g{z} (nuclei/um2)", 
            "nuclei_per_muscle_fibre g{z}"
        ]
    )

    # ========================================================
    # SAVE MULTI-SHEET EXCEL REPORT
    # ========================================================
    timestamp = datetime.now().strftime("%Y%m%d_%H%M")
    output_path = RESULT_DIR / f"MultiFeature_Analysis_{timestamp}.xlsx"

    print(f"[STAGE] Saving report to Excel: {output_path}")
    with pd.ExcelWriter(output_path, engine="openpyxl") as writer:
        df_pax7.to_excel(writer, sheet_name="Pax7_Features", index=False)
        df_fibres.to_excel(writer, sheet_name="Muscle_Fibres", index=False)
        df_nuclei.to_excel(writer, sheet_name="Nuclei", index=False)

    print("\n" + "=" * 60)
    print("[SUCCESS] PIPELINE FINISHED")
    print(f"Report saved as:\n{output_path}")
    print("=" * 60)