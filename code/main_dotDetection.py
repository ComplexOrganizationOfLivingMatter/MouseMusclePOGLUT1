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
    """Reads an image and extracts the first channel if multi-channel."""
    img = io.imread(str(path))
    return img[..., 0] if img.ndim == 3 else img


def process_single_image(mask_path: Path, pred_dict: dict, output_mask_root: Path, output_vor_root: Path) -> dict | None:
    img_name = mask_path.name
    age_folder = mask_path.parent.name
    genotype = "KIKO" if "KI" in img_name.upper() else "WT"

    mask_file_path = output_mask_root / age_folder / img_name
    vor_file_path = output_vor_root / age_folder / img_name

    mask_file_path.parent.mkdir(parents=True, exist_ok=True)
    vor_file_path.parent.mkdir(parents=True, exist_ok=True)

    try:
        # Load prediction first to know its dimensions
        pred_path = pred_dict.get(img_name)
        if pred_path is None:
            print(f"[WARN] Prediction not found: {img_name}")
            return None

        img_pred = load_single_channel(pred_path)
        pred_h, pred_w = img_pred.shape

        mask_zones = load_single_channel(mask_path)

        # --- Dimension Check & Nearest Neighbor Transformation for Mask ---
        if mask_zones.shape != (pred_h, pred_w):
            print(f"[INFO] Dimension mismatch for {img_name}: Mask {mask_zones.shape} vs Prediction {(pred_h, pred_w)}. Resizing mask using nearest neighbor.")
            mask_zones = cv2.resize(mask_zones, (pred_w, pred_h), interpolation=cv2.INTER_NEAREST)

        h, w = mask_zones.shape
        use_cache = False

        # --- Cache Check ---
        if mask_file_path.exists() and vor_file_path.exists():
            try:
                binary_nuclei = load_single_channel(mask_file_path) > 0
                global_voronoi = io.imread(str(vor_file_path)).astype(np.uint32)
                
                final_nuclei_mask = measure.label(
                    (global_voronoi * binary_nuclei) > 0, connectivity=1
                ).astype(np.uint32)
                
                cached_count = len(measure.regionprops(final_nuclei_mask))
                print(f"[CACHE] {img_name}: {cached_count} nuclei recovered")
                use_cache = True
            except Exception as e:
                print(f"[WARN] Cache failed for {img_name}. Reprocessing. Reason: {e}")
                
        # --- Prediction & Voronoi Generation ---
        if not use_cache:
            # Label components & remove excluded zones
            raw_mask = measure.label(img_pred > 0, connectivity=1)
            raw_mask[np.isin(mask_zones, list(EXCLUDED_ZONES))] = 0
            final_nuclei_mask = measure.label(raw_mask > 0, connectivity=1).astype(np.uint32)

            # Save binary nuclei mask
            tifffile.imwrite(
                str(mask_file_path), (final_nuclei_mask > 0).astype(np.uint8) * 255
            )

            # Group nuclei centroids per zone
            global_voronoi = np.zeros((h, w), dtype=np.uint32)
            all_props = measure.regionprops(final_nuclei_mask)
            zone_to_nuclei = {z: [] for z in VALID_ZONES}

            for p in all_props:
                cy, cx = map(int, p.centroid)
                z_id = mask_zones[cy, cx]
                if z_id in zone_to_nuclei:
                    zone_to_nuclei[z_id].append(p)

            # Local Voronoi per zone
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
            print(f"[DONE] {img_name}: {len(all_props)} nuclei detected")

        # --- Data Extraction ---
        props = measure.regionprops(final_nuclei_mask)
        if not props:
            return {"image": img_name, "age": age_folder, "genotype": genotype, "status": "no_nuclei"}

        centroids = np.array([p.centroid for p in props]).astype(int)
        zone_ids = mask_zones[centroids[:, 0], centroids[:, 1]]

        res = {"image": img_name, "age": age_folder, "genotype": genotype}
        total_area, total_nuclei = 0.0, 0

        for z in VALID_ZONES:
            z_mask = mask_zones == z
            if not np.any(z_mask):
                continue

            z_area = np.sum(z_mask) * PIXEL_AREA_UM2
            z_nuclei_count = int(np.sum(zone_ids == z))

            total_area += z_area
            total_nuclei += z_nuclei_count

            res[f"area g{z} (um2)"] = z_area
            res[f"nuclei/um2 g{z}"] = z_nuclei_count / z_area if z_area > 0 else 0

            # Voronoi cell size statistics (mean & std in um2)
            voronoi_cells = global_voronoi[z_mask]
            voronoi_cells = voronoi_cells[voronoi_cells > 0]
            if voronoi_cells.size > 0:
                _, territory_sizes = np.unique(voronoi_cells, return_counts=True)
                cell_areas = territory_sizes * PIXEL_AREA_UM2
                res[f"voronoi_mean g{z} (um2)"] = float(np.mean(cell_areas))
                res[f"voronoi_std g{z} (um2)"] = float(np.std(cell_areas))
            else:
                res[f"voronoi_mean g{z} (um2)"] = np.nan
                res[f"voronoi_std g{z} (um2)"] = np.nan

        res["total nuclei detected"] = total_nuclei
        res["whole muscle area (um2)"] = total_area
        res["whole muscle nuclei/um2"] = total_nuclei / total_area if total_area > 0 else 0

        return res

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
    PRED_DICT = {p.name: p for p in PATH_PREDS}

    RESULT_DIR = BASE_DIR / "results" / "largeModel_dotDetection"
    QC_ROOT = RESULT_DIR / "binaryNucleiDots"
    VOR_ROOT = RESULT_DIR / "voronoiImages"

    QC_ROOT.mkdir(parents=True, exist_ok=True)
    VOR_ROOT.mkdir(parents=True, exist_ok=True)

    print("\n" + "=" * 60)
    print("STARTING DOT-BASED NUCLEI PIPELINE")
    print(f"Muscle masks found: {len(PATH_MASKS)}")
    print(f"BiaPy predictions found: {len(PATH_PREDS)}")
    print("=" * 60 + "\n")

    results_list = []
    with ProcessPoolExecutor() as executor:
        futures = {
            executor.submit(process_single_image, mask_path, PRED_DICT, QC_ROOT, VOR_ROOT): mask_path
            for mask_path in PATH_MASKS
        }

        for i, future in enumerate(as_completed(futures), start=1):
            res = future.result()
            if res is not None:
                results_list.append(res)

            if i % 5 == 0 or i == len(PATH_MASKS):
                print(f"[PROGRESS] {i}/{len(PATH_MASKS)} images complete")

    # --- DataFrame Formatting & Ordering ---
    print("\n[STAGE] Creating dataframe...")
    df = pd.DataFrame(results_list)
    if df.empty:
        raise RuntimeError("No images were successfully processed.")

    df = df.loc[:, ~df.columns.duplicated()]

    base_cols = [
        "image", "age", "genotype",
        "total nuclei detected", "whole muscle area (um2)", "whole muscle nuclei/um2"
    ]
    group_cols = [
        col
        for z in VALID_ZONES
        for col in [
            f"area g{z} (um2)", 
            f"pax7/um2 g{z}", 
            f"voronoi_mean g{z} (um2)", 
            f"voronoi_std g{z} (um2)"
        ]
        if col in df.columns
    ]

    df = df.reindex(columns=base_cols + group_cols)

    # ========================================================
    # GENERATE SUMMARY TABLES
    # ========================================================
    print("[STAGE] Generating summary statistics...")

    numeric_cols = df.select_dtypes(include="number").columns

    summary_final = df.groupby(["genotype", "age"], sort=False)[numeric_cols].agg(
        ["mean", "std"]
    )

    summary_final.columns = [f"{col}_{stat}" for col, stat in summary_final.columns]

    # ========================================================
    # SAVE EXCEL REPORT
    # ========================================================
    timestamp = datetime.now().strftime("%Y%m%d_%H%M")
    output_path = RESULT_DIR / f"Pax7_DotDetection_Analysis_{timestamp}.xlsx"

    with pd.ExcelWriter(output_path, engine="openpyxl") as writer:
        df.to_excel(writer, sheet_name="Raw_Data", index=False)
        summary_final.to_excel(writer, sheet_name="Summary_Stats")

    # ========================================================
    # FINISHED
    # ========================================================
    print("\n" + "=" * 60)
    print("[SUCCESS] PIPELINE FINISHED")
    print(f"Report saved as:\n{output_path}")
    print("=" * 60)