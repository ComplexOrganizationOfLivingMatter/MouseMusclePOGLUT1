from datetime import datetime
import gc
from pathlib import Path
import cv2
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from skimage import io, measure

# =========================================================
# CONFIGURATION
# =========================================================
AGES = ["p0", "p4", "p7", "p10", "p15", "p20", "p35"]
UM_PER_PIXEL = 0.5055443
PIXEL_AREA_UM2 = UM_PER_PIXEL**2
MODEL_FOLDER = "largeModel_dotDetection"

BONE_ZONE = 6
EXCLUDED_ZONES = {0, BONE_ZONE, 13}
VALID_ZONES = [z for z in range(1, 16) if z not in EXCLUDED_ZONES]

# Paths Setup
BASE_DIR = Path("../../")
PATH_MASKS = BASE_DIR / "data" / "annotatedMuscleMasks"
PATH_VORONOI = BASE_DIR / "results" / MODEL_FOLDER / "voronoiImages"
PATH_CELLPOSE_FIBRES = (
    BASE_DIR / "results" / "cellposeSegmentation" / "muscleFibres"
)

# Output Setup
RESULT_DIR = BASE_DIR / "results" / MODEL_FOLDER / "heatmaps_combined_comparison"
RESULT_DIR.mkdir(parents=True, exist_ok=True)

# Separate Output Directories for Cellpose modes & Voronoi
HEATMAPS_CELLPOSE_REGION_DIR = RESULT_DIR / "heatmaps_cellpose_region_by_region"
HEATMAPS_CELLPOSE_WHOLELEG_DIR = RESULT_DIR / "heatmaps_cellpose_whole_leg"
HEATMAPS_VORONOI_DIR = RESULT_DIR / "heatmaps_voronoi_region_by_region"

for d in [HEATMAPS_CELLPOSE_REGION_DIR, HEATMAPS_CELLPOSE_WHOLELEG_DIR, HEATMAPS_VORONOI_DIR]:
  d.mkdir(parents=True, exist_ok=True)

# Visualization Settings
CMAP = "RdBu_r"  # Blue: Small/Dense | Red: Big/Sparse
V_LIMIT = 0.6  # Color scale limit for log10

# Toggle resolution for separate individual heatmaps:
SAVE_FULL_RESOLUTION_HEATMAPS = True
RENDER_MAX_DIM = 2000


# =========================================================
# UTILITIES & DEBUG LOGGING
# =========================================================
def log_step(message, level="INFO"):
  time_str = datetime.now().strftime("%H:%M:%S")
  prefix = {
      "INFO": "[INFO]",
      "DEBUG": "[DEBUG]",
      "WARN": "[WARNING]",
      "ERROR": "[ERROR]",
  }.get(level, "[INFO]")
  print(f"{time_str} {prefix} {message}")


def clean_muscle_name(filename, is_ki):
  stem = Path(filename).stem
  clean_str = stem.replace("Muscle_", "").replace("_", " ")
  parts = clean_str.split()
  filtered_parts = [p for p in parts if p.upper() not in ["WT", "KI", "KIKO"]]
  prefix = "KIKO" if is_ki else "WT"
  return f"{prefix} {' '.join(filtered_parts)}" if filtered_parts else f"{prefix} {stem}"


def resize_image_safely(img, max_dim=RENDER_MAX_DIM):
  h, w = img.shape[:2]
  if max(h, w) <= max_dim:
    return img
  scale = max_dim / max(h, w)
  new_size = (int(w * scale), int(h * scale))
  img_float = img.astype(np.float32)
  resized = cv2.resize(img_float, new_size, interpolation=cv2.INTER_NEAREST)
  return resized.astype(img.dtype)


def create_thick_outlines(label_img):
  edges = np.zeros(label_img.shape, dtype=bool)
  edges[:-1, :] |= label_img[:-1, :] != label_img[1:, :]
  edges[:, :-1] |= label_img[:, :-1] != label_img[:, 1:]
  kernel = np.ones((2, 2), np.uint8)
  thick_edges = cv2.dilate(edges.astype(np.uint8), kernel)
  return np.ma.masked_where(thick_edges == 0, thick_edges)


# =========================================================
# DYNAMIC WT BASELINE CALCULATIONS
# =========================================================
def compute_wt_baselines_per_zone(source_dir: Path, label_type: str):
  """Computes per-region WT dynamic baselines (zone by zone) for each age."""
  log_step(f"Computing per-region WT dynamic baselines for {label_type}...")
  zone_areas_by_age_zone = {age: {z: [] for z in VALID_ZONES} for age in AGES}

  for age in AGES:
    age_src_folder = source_dir / age
    if not age_src_folder.exists():
      continue

    wt_files = [
        f
        for f in age_src_folder.glob("*.tif")
        if "KI" not in f.name.upper() and "POGLU" not in f.name.upper()
    ]
    log_step(f"[{label_type}] Age {age}: Found {len(wt_files)} WT baseline files.")

    for tif_path in wt_files:
      img_name = tif_path.name
      img_label = io.imread(str(tif_path))
      if img_label.ndim == 3:
        img_label = img_label[..., 0]

      h, w = img_label.shape
      mask_path = PATH_MASKS / age / img_name
      if not mask_path.exists():
        continue

      mask_zones = io.imread(str(mask_path))
      if mask_zones.ndim == 3:
        mask_zones = mask_zones[..., 0]
      if mask_zones.shape != (h, w):
        mask_zones = cv2.resize(mask_zones, (w, h), interpolation=cv2.INTER_NEAREST)

      props = measure.regionprops(img_label.astype(np.int32))
      for p in props:
        if p.label == 0:
          continue
        cy, cx = map(int, p.centroid)
        if 0 <= cy < h and 0 <= cx < w:
          z_id = mask_zones[cy, cx]
          if z_id in VALID_ZONES:
            cell_area_um2 = p.area * PIXEL_AREA_UM2
            zone_areas_by_age_zone[age][z_id].append(cell_area_um2)

  wt_baselines = {}
  for age in AGES:
    wt_baselines[age] = {}
    for z in VALID_ZONES:
      areas = zone_areas_by_age_zone[age][z]
      wt_baselines[age][z] = np.mean(areas) if len(areas) > 0 else np.nan

  log_step(f"Per-region WT baselines computed successfully for {label_type}.")
  return wt_baselines


def compute_wt_baselines_whole_leg(source_dir: Path, label_type: str):
  """Computes whole-leg WT dynamic baselines (global average across all valid muscle regions) for each age."""
  log_step(f"Computing whole-leg WT dynamic baselines for {label_type}...")
  whole_leg_areas_by_age = {age: [] for age in AGES}

  for age in AGES:
    age_src_folder = source_dir / age
    if not age_src_folder.exists():
      continue

    wt_files = [
        f
        for f in age_src_folder.glob("*.tif")
        if "KI" not in f.name.upper() and "POGLU" not in f.name.upper()
    ]

    for tif_path in wt_files:
      img_name = tif_path.name
      img_label = io.imread(str(tif_path))
      if img_label.ndim == 3:
        img_label = img_label[..., 0]

      h, w = img_label.shape
      mask_path = PATH_MASKS / age / img_name
      if not mask_path.exists():
        continue

      mask_zones = io.imread(str(mask_path))
      if mask_zones.ndim == 3:
        mask_zones = mask_zones[..., 0]
      if mask_zones.shape != (h, w):
        mask_zones = cv2.resize(mask_zones, (w, h), interpolation=cv2.INTER_NEAREST)

      props = measure.regionprops(img_label.astype(np.int32))
      for p in props:
        if p.label == 0:
          continue
        cy, cx = map(int, p.centroid)
        if 0 <= cy < h and 0 <= cx < w:
          z_id = mask_zones[cy, cx]
          if z_id in VALID_ZONES:
            cell_area_um2 = p.area * PIXEL_AREA_UM2
            whole_leg_areas_by_age[age].append(cell_area_um2)

  wt_whole_leg_baselines = {}
  for age in AGES:
    areas = whole_leg_areas_by_age[age]
    wt_whole_leg_baselines[age] = np.mean(areas) if len(areas) > 0 else np.nan

  log_step(f"Whole-leg WT baselines computed successfully for {label_type}.")
  return wt_whole_leg_baselines


# =========================================================
# SEQUENTIAL PROCESSING WITH BOUNDARY CONTOURS & DIRECTORIES
# =========================================================
def process_heatmap_set_sequential(
    source_dir: Path,
    output_base_dir: Path,
    wt_baselines_dict: dict,
    label_type: str,
    mode: str,
):
  log_step(f"Starting sequential processing for: {label_type} (Mode: {mode}) from {source_dir}")

  if not source_dir.exists():
    log_step(f"Source directory does not exist: {source_dir}", level="ERROR")
    return

  total_processed = 0

  for age in AGES:
    age_src_folder = source_dir / age
    if not age_src_folder.exists():
      continue

    age_baselines = wt_baselines_dict.get(age, {})
    age_out_dir = output_base_dir / age
    age_out_dir.mkdir(parents=True, exist_ok=True)

    tif_files = list(age_src_folder.glob("*.tif"))
    log_step(f"[{label_type} - {mode}] Processing age {age}: {len(tif_files)} files found.")

    for idx, tif_path in enumerate(tif_files, start=1):
      img_name = tif_path.name
      is_ki = "KI" in img_name.upper()

      try:
        img_label = io.imread(str(tif_path))
        if img_label.ndim == 3:
          img_label = img_label[..., 0]
        img_label = img_label.astype(np.int32)

        h, w = img_label.shape
        mask_path = PATH_MASKS / age / img_name
        if not mask_path.exists():
          continue

        mask_zones = io.imread(str(mask_path))
        if mask_zones.ndim == 3:
          mask_zones = mask_zones[..., 0]
        if mask_zones.shape != (h, w):
          mask_zones = cv2.resize(mask_zones, (w, h), interpolation=cv2.INTER_NEAREST)

        valid_mask = np.isin(mask_zones, VALID_ZONES)
        img_label_masked = img_label.copy()
        img_label_masked[~valid_mask] = 0

        props = measure.regionprops(img_label_masked)
        if len(props) == 0:
          continue

        lookup = np.zeros(img_label.max() + 1, dtype=np.float32)

        for p in props:
          if p.label == 0:
            continue
          cy, cx = map(int, p.centroid)
          if 0 <= cy < h and 0 <= cx < w:
            z_id = mask_zones[cy, cx]
            if z_id in VALID_ZONES:
              cell_area_um2 = p.area * PIXEL_AREA_UM2

              if mode == "region_by_region":
                wt_mean_um2 = age_baselines.get(z_id, np.nan)
              else:  # whole_leg
                wt_mean_um2 = age_baselines

              if not np.isnan(wt_mean_um2) and wt_mean_um2 > 0:
                ratio = np.log10(cell_area_um2 / wt_mean_um2)
                lookup[p.label] = np.clip(ratio, -V_LIMIT, V_LIMIT)

        if SAVE_FULL_RESOLUTION_HEATMAPS:
          heatmap = lookup[img_label_masked]
          edge_overlay = create_thick_outlines(img_label_masked)
          mask_render = mask_zones
          dpi = 300
          fig_w, fig_h = w / dpi, h / dpi
        else:
          img_render = resize_image_safely(img_label_masked)
          heatmap = lookup[img_render]
          edge_overlay = create_thick_outlines(img_render)
          mask_render = resize_image_safely(mask_zones)
          dpi = 120
          fig_w, fig_h = 10, 10

        heatmap[~np.isin(mask_render, VALID_ZONES)] = np.nan

        fig, ax = plt.subplots(figsize=(fig_w, fig_h), facecolor="black")
        ax.set_facecolor("black")

        clean_name = clean_muscle_name(tif_path.stem, is_ki)
        ax.set_title(f"{clean_name} ({label_type} - {mode})", color="white", fontsize=12 if SAVE_FULL_RESOLUTION_HEATMAPS else 14)

        ax.imshow(heatmap, cmap=CMAP, vmin=-V_LIMIT, vmax=V_LIMIT)
        ax.imshow(edge_overlay, cmap="binary_r", interpolation="nearest", alpha=0.9)

        # 1. Outer boundaries (White)
        mask_bin = (np.isin(mask_render, VALID_ZONES)).astype(np.uint8)
        contours, _ = cv2.findContours(mask_bin, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
        for cnt in contours:
          poly = cnt.reshape(-1, 2)
          ax.plot(poly[:, 0], poly[:, 1], color="white", lw=0.8 if SAVE_FULL_RESOLUTION_HEATMAPS else 1.0, alpha=0.6)

        # 2. Inner muscle group / zone boundaries (Black)
        for z in VALID_ZONES:
          zone_bin = (mask_render == z).astype(np.uint8)
          if np.any(zone_bin):
            z_contours, _ = cv2.findContours(zone_bin, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
            for z_cnt in z_contours:
              z_poly = z_cnt.reshape(-1, 2)
              ax.plot(z_poly[:, 0], z_poly[:, 1], color="black", lw=0.8 if SAVE_FULL_RESOLUTION_HEATMAPS else 1.0, alpha=0.8)

        ax.axis("off")
        output_heatmap_path = age_out_dir / f"{tif_path.stem}_heatmap.png"
        plt.savefig(
            output_heatmap_path,
            facecolor="black",
            edgecolor="none",
            bbox_inches="tight",
            pad_inches=0 if SAVE_FULL_RESOLUTION_HEATMAPS else 0.1,
            dpi=dpi,
        )
        plt.close(fig)
        plt.close('all')
        total_processed += 1
        gc.collect()

      except Exception as e:
        log_step(f"Error processing {img_name}: {str(e)}", level="ERROR")

  log_step(f"Completed {label_type} ({mode}) processing. Total heatmaps generated: {total_processed}")


# =========================================================
# MONTAGE GENERATION
# =========================================================
def create_comparison_montage(age, heatmap_subfolder: Path, title_suffix: str, mode: str):
  age_folder = heatmap_subfolder / age
  if not age_folder.exists():
    return

  all_pngs = list(age_folder.glob("*.png"))
  wt_files = sorted([f for f in all_pngs if "KI" not in f.name.upper()])
  ki_files = sorted([f for f in all_pngs if "KI" in f.name.upper()])

  n_cols = max(len(wt_files), len(ki_files))
  if n_cols == 0:
    return

  fig, axes = plt.subplots(2, n_cols, figsize=(n_cols * 5, 11), facecolor="black")
  fig.patch.set_facecolor("black")
  fig.suptitle(
      f"{age.upper()} {title_suffix} Comparison ({mode})",
      color="white",
      fontsize=20,
      fontweight="bold",
      y=0.95,
  )

  if n_cols == 1:
    axes = axes.reshape(2, 1)

  for row_idx, files in enumerate([wt_files, ki_files]):
    for col_idx in range(n_cols):
      ax = axes[row_idx, col_idx]
      ax.set_facecolor("black")

      if col_idx < len(files):
        file_to_load = files[col_idx]
        try:
          img = mpimg.imread(str(file_to_load))
          ax.imshow(img)
        except Exception as e:
          log_step(f"Failed to load PNG for montage {file_to_load.name}: {e}", level="ERROR")

      ax.axis("off")

  sm = plt.cm.ScalarMappable(cmap=CMAP, norm=plt.Normalize(vmin=-V_LIMIT, vmax=V_LIMIT))
  sm._A = []

  cax = fig.add_axes([0.3, 0.07, 0.4, 0.02])
  cbar = fig.colorbar(sm, cax=cax, orientation="horizontal")
  cbar.set_label(r"$\log_{10}(\frac{Cell\_Area}{WT\_Mean\_Area})$", color="white", size=16, labelpad=10)
  cbar.ax.xaxis.set_tick_params(color="white", labelcolor="white")
  cbar.outline.set_edgecolor("white")

  output_filename = f"Montage_{age}_{title_suffix}_{mode}_Comparison.png"
  output_path = RESULT_DIR / output_filename
  plt.savefig(
      output_path,
      facecolor="black",
      edgecolor="none",
      bbox_inches="tight",
      transparent=False,
      dpi=150,
  )
  plt.close(fig)
  plt.close('all')
  gc.collect()
  log_step(f"Successfully saved and cleared montage: {output_filename}")


# =========================================================
# EXECUTION (CELLPOSE FIRST, THEN VORONOI)
# =========================================================
if __name__ == "__main__":
  log_step("=== WORKFLOW INITIALIZED (CELLPOSE DUAL-MODE & VORONOI REGION-WISE) ===")

  # ---------------------------------------------------------
  # 1. PROCESS CELLPOSE (Region-by-Region)
  # ---------------------------------------------------------
  log_step("--- PROCESSING CELLPOSE (Region-by-Region) ---")
  wt_cp_region = compute_wt_baselines_per_zone(PATH_CELLPOSE_FIBRES, "Cellpose")
  process_heatmap_set_sequential(
      source_dir=PATH_CELLPOSE_FIBRES,
      output_base_dir=HEATMAPS_CELLPOSE_REGION_DIR,
      wt_baselines_dict=wt_cp_region,
      label_type="Cellpose",
      mode="region_by_region",
  )

  # ---------------------------------------------------------
  # 2. PROCESS CELLPOSE (Whole-Leg)
  # ---------------------------------------------------------
  log_step("--- PROCESSING CELLPOSE (Whole-Leg) ---")
  wt_cp_whole = compute_wt_baselines_whole_leg(PATH_CELLPOSE_FIBRES, "Cellpose")
  process_heatmap_set_sequential(
      source_dir=PATH_CELLPOSE_FIBRES,
      output_base_dir=HEATMAPS_CELLPOSE_WHOLELEG_DIR,
      wt_baselines_dict=wt_cp_whole,
      label_type="Cellpose",
      mode="whole_leg",
  )

  # ---------------------------------------------------------
  # 3. PROCESS VORONOI (Region-by-Region ONLY)
  # ---------------------------------------------------------
  log_step("--- PROCESSING VORONOI (Region-by-Region) ---")
  wt_vor_region = compute_wt_baselines_per_zone(PATH_VORONOI, "Voronoi")
  process_heatmap_set_sequential(
      source_dir=PATH_VORONOI,
      output_base_dir=HEATMAPS_VORONOI_DIR,
      wt_baselines_dict=wt_vor_region,
      label_type="Voronoi",
      mode="region_by_region",
  )

  # ---------------------------------------------------------
  # 4. MONTAGE COMPILATION PHASE
  # ---------------------------------------------------------
  log_step("Starting Montage Compilation Phase...")
  for age in AGES:
    log_step(f"Compiling Cellpose region-wise montages for age {age}...")
    create_comparison_montage(age, HEATMAPS_CELLPOSE_REGION_DIR, "Cellpose", mode="region_by_region")

    log_step(f"Compiling Cellpose whole-leg montages for age {age}...")
    create_comparison_montage(age, HEATMAPS_CELLPOSE_WHOLELEG_DIR, "Cellpose", mode="whole_leg")

    log_step(f"Compiling Voronoi region-wise montages for age {age}...")
    create_comparison_montage(age, HEATMAPS_VORONOI_DIR, "Voronoi", mode="region_by_region")

  log_step("=== WORKFLOW FINISHED SUCCESSFULLY ===")