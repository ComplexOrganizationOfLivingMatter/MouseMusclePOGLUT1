# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 13:05:55 2026

@author: pedro
"""

import os
import numpy as np
import tifffile as tiff
import json

# ============================================================
# TRAINING DATA PATH
# ============================================================

train_folder = r"F:\Lab\MouseMusclePOGLUT1\data\trainingDataset\individual label\Pax7_v2_Mar\after_curation_cleanArtifacts\train\raw"


# ============================================================
# CLEAN SIGNAL (5–99 percentile)
# ============================================================

def clean_signal(x):
    low = np.percentile(x, 5)
    high = np.percentile(x, 99)
    return x[(x > low) & (x < high)]


# ============================================================
# LOAD IMAGE
# ============================================================

def load_img(path):
    return tiff.imread(path)


# ============================================================
# COMPUTE STATS
# ============================================================

def compute_stats(folder):

    pax7_vals = []
    col_vals  = []
    nuc_vals  = []

    files = [f for f in os.listdir(folder) if f.endswith(".tif")]

    print(f"Found {len(files)} training images")

    for f in files:

        img = load_img(os.path.join(folder, f))

        if img.ndim != 3 or img.shape[-1] != 3:
            continue

        pax7_vals.append(clean_signal(img[..., 0]))
        col_vals.append(clean_signal(img[..., 1]))
        nuc_vals.append(clean_signal(img[..., 2]))

    pax7 = np.concatenate(pax7_vals)
    col  = np.concatenate(col_vals)
    nuc  = np.concatenate(nuc_vals)

    stats = {
        "pax7_mean": float(np.mean(pax7)),
        "pax7_std": float(np.std(pax7)),

        "collagen_mean": float(np.mean(col)),
        "collagen_std": float(np.std(col)),

        "nuclei_mean": float(np.mean(nuc)),
        "nuclei_std": float(np.std(nuc)),
    }

    return stats


# ============================================================
# SAVE JSON
# ============================================================

stats = compute_stats(train_folder)

out_file = os.path.join(train_folder, "reference_stats.json")

with open(out_file, "w") as f:
    json.dump(stats, f, indent=4)

print("Saved:", out_file)
print(stats)