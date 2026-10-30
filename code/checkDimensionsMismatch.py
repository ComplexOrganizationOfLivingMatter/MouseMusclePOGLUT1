# -*- coding: utf-8 -*-
"""
Created on Thu Oct 29 14:40:46 2026

@author: pedro
"""

"""
Standalone script to check and find dimension mismatches between 
annotated muscle masks and prediction images.
"""

from pathlib import Path
from skimage import io
import pandas as pd


def load_single_channel_shape(path: Path) -> tuple:
    """Reads an image shape and extracts spatial dimensions if multi-channel."""
    img = io.imread(str(path))
    # Return (height, width) regardless of whether it's 2D or has color/channel dimensions
    return img.shape[:2] if img.ndim >= 3 else img.shape


def main():
    BASE_DIR = Path("..")
    
    # Define paths
    masks_dir = BASE_DIR / "data" / "annotatedMuscleMasks"
    preds_dir = BASE_DIR / "data" / "predictions" / "all_largeModel_detection" / "per_image_local_max_check_post_processing"

    PATH_MASKS = list(masks_dir.rglob("*.tif"))
    PATH_PREDS = list(preds_dir.rglob("*.tif"))

    # Map prediction paths to filenames for O(1) lookups
    PRED_DICT = {p.name: p for p in PATH_PREDS}

    print("=" * 60)
    print("CHECKING IMAGE PAIR DIMENSIONS")
    print(f"Total masks found: {len(PATH_MASKS)}")
    print(f"Total predictions found: {len(PATH_PREDS)}")
    print("=" * 60)

    mismatches = []
    missing_preds = []
    checked_count = 0

    for mask_path in PATH_MASKS:
        img_name = mask_path.name
        
        if img_name not in PRED_DICT:
            missing_preds.append(img_name)
            continue

        pred_path = PRED_DICT[img_name]
        checked_count += 1

        # Get spatial dimensions (height, width)
        mask_shape = load_single_channel_shape(mask_path)
        pred_shape = load_single_channel_shape(pred_path)

        if mask_shape != pred_shape:
            mismatches.append({
                "image": img_name,
                "mask_shape": mask_shape,
                "pred_shape": pred_shape,
                "mask_path": str(mask_path),
                "pred_path": str(pred_path)
            })

    # ========================================================
    # REPORT FINDINGS
    # ========================================================
    print(f"\n[SUMMARY]")
    print(f"  - Successfully paired and checked: {checked_count}")
    print(f"  - Missing predictions: {len(missing_preds)}")
    print(f"  - Dimension mismatches found: {len(mismatches)}")

    if missing_preds:
        print(f"\n[WARNING] The following {len(missing_preds)} masks had no matching prediction file:")
        for name in missing_preds[:10]:  # Show first 10
            print(f"    * {name}")
        if len(missing_preds) > 10:
            print(f"    ... and {len(missing_preds) - 10} more.")

    if mismatches:
        print(f"\n[MISMATCHES] Found {len(mismatches)} pairs with dimension differences:")
        for m in mismatches:
            print(f"  - {m['image']}")
            print(f"      Mask Shape:       {m['mask_shape']}")
            print(f"      Prediction Shape: {m['pred_shape']}")

        # Save mismatch details to a CSV report
        output_csv = BASE_DIR / "results" / "dimension_mismatches_report.csv"
        output_csv.parent.mkdir(parents=True, exist_ok=True)
        
        df_mismatches = pd.DataFrame(mismatches)
        df_mismatches.to_csv(output_csv, index=False)
        print(f"\n[SAVED] Detailed mismatch report saved to:\n  {output_csv}")
    else:
        print("\n[SUCCESS] All paired images have identical dimensions. No mismatches found!")

    print("=" * 60)


if __name__ == "__main__":
    main()