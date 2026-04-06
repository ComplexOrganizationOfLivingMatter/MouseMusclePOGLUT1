import os
import glob
import math
import numpy as np
import pandas as pd
from datetime import datetime
from skimage import io, measure, morphology, filters
from scipy import ndimage
from scipy.spatial.distance import cdist

# %% ================= USER INPUT =================
# Ground truth folder containing labeled masks
root_path = r'F:\Lab\MouseMusclePOGLUT1\results\detectionMetrics'

gt_folder = os.path.join(root_path, 'imagesTest', 'GT emi_vicky combined')

# Folder containing predicted masks
name_pred_images = 'modelV1'
pred_folder = os.path.join(root_path, 'imagesTest', name_pred_images)

dist_thresh = 10       # Distance threshold (pixels) to consider a centroid matched
out_excel = os.path.join(root_path, f'detection_evaluation_{name_pred_images}_autobinary_160326.xlsx')
img_ext = '*.tif'      # Image file extension
binary_thr = 250       # Threshold to binarize predicted masks (unused in active logic, matching MATLAB)

# Min nuclei size
microns_per_pixel = 0.5055443
min_area_um2 = 10 # minimum nuclei size in um^2
min_area_pixels = math.ceil(min_area_um2 / (microns_per_pixel**2))

# %% ==============================================

# List all ground truth files
search_pattern = os.path.join(gt_folder, img_ext)
gt_files = glob.glob(search_pattern)

n_images = len(gt_files) # Total number of images
eps = np.finfo(float).eps # Machine epsilon to prevent division by zero

# Initialize lists to store per-image metrics (lists are more Pythonic for building DataFrames than pre-allocating zeros)
image_names = []
TPs, FPs, FNs, GTs, PRs, Precisions, Recalls, F1s = [], [], [], [], [], [], [], []

# %% ================= MAIN LOOP =================
for img_path in gt_files:
    file_name = os.path.basename(img_path)
    image_names.append(file_name)
    
    # --- Load images ---
    gt = io.imread(img_path)
    pr_path = os.path.join(pred_folder, file_name)
    pr = io.imread(pr_path)
    
    if pr.dtype != np.uint8:
        BW = (pr * 255).astype(np.uint8)
    else:
        BW = pr
        
    # --- Binarize images ---
    gt_bin = gt > 0
    
    # Predicted mask: threshold (Otsu) + fill holes + morphological close
    # Handle the case where the image might be completely black/constant (Otsu would fail)
    if BW.min() == BW.max():
        pr_bin = BW > 0
    else:
        thresh = filters.threshold_otsu(BW)
        pr_bin = BW > thresh
        
    # Fill holes
    pr_bin = ndimage.binary_fill_holes(pr_bin)
    
    # Morphological close with disk of radius 1
    disk_1 = morphology.disk(1)
    pr_bin = morphology.binary_closing(pr_bin, disk_1)
    
    # Remove small nuclei
    pr_bin = morphology.remove_small_objects(pr_bin, min_size=min_area_pixels)
    
    # --- Label connected components ---
    # connectivity=2 in 2D is equivalent to 8-connectivity in MATLAB
    gt_L = measure.label(gt_bin, connectivity=2)
    pr_L = measure.label(pr_bin, connectivity=2)
    
    # --- Extract centroids ---
    gt_props = measure.regionprops(gt_L)
    pr_props = measure.regionprops(pr_L)
    
    n_GT = len(gt_props)
    n_PR = len(pr_props)
    
    GTs.append(n_GT)
    PRs.append(n_PR)
    
    # --- Edge cases: zero objects ---
    if n_GT == 0 and n_PR == 0:
        TPs.append(0); FPs.append(0); FNs.append(0)
        Precisions.append(1.0); Recalls.append(1.0); F1s.append(1.0)
        continue
        
    if n_GT == 0:
        TPs.append(0); FPs.append(n_PR); FNs.append(0)
        Precisions.append(0.0); Recalls.append(1.0); F1s.append(0.0)
        continue
        
    if n_PR == 0:
        TPs.append(0); FPs.append(0); FNs.append(n_GT)
        Precisions.append(1.0); Recalls.append(0.0); F1s.append(0.0)
        continue
        
    # --- Build centroid arrays ---
    # skimage centroids are (row, col) representing (y, x). Order doesn't matter for distance.
    gt_C = np.array([prop.centroid for prop in gt_props]) # n_GT x 2
    pr_C = np.array([prop.centroid for prop in pr_props]) # n_PR x 2
    
    # --- Compute distance matrix between GT and predicted centroids ---
    D = cdist(gt_C, pr_C) # n_GT x n_PR matrix
    
    matched_GT = np.zeros(n_GT, dtype=bool)
    matched_PR = np.zeros(n_PR, dtype=bool)
    
    # --- Greedy one-to-one matching ---
    TP_count = 0
    while True:
        min_D = np.min(D)
        if min_D > dist_thresh:
            break
            
        # Get 2D indices of the minimum value
        g, p = np.unravel_index(np.argmin(D), D.shape)
        
        # Mark as matched
        TP_count += 1
        matched_GT[g] = True
        matched_PR[p] = True
        
        # Remove matched row and column from distance matrix by setting to infinity
        D[g, :] = np.inf
        D[:, p] = np.inf
        
    # --- Store per-image results ---
    TPs.append(TP_count)
    FNs.append(np.sum(~matched_GT))
    FPs.append(np.sum(~matched_PR))
    
    # Compute Precision / Recall / F1 per image
    p_val = TP_count / (TP_count + FPs[-1] + eps)
    r_val = TP_count / (TP_count + FNs[-1] + eps)
    f1_val = 2 * p_val * r_val / (p_val + r_val + eps)
    
    Precisions.append(p_val)
    Recalls.append(r_val)
    F1s.append(f1_val)

# %% ================= CREATE TABLE =================

# Per-image table
df = pd.DataFrame({
    'ImageName': image_names,
    'TP': TPs,
    'FP': FPs,
    'FN': FNs,
    'Total_GT': GTs,
    'Total_Pred': PRs,
    'Precision': Precisions,
    'Recall': Recalls,
    'F1': F1s,
    'DistanceThreshold_px': [dist_thresh] * n_images
})

# --- Add a final row with averages across all images ---
avg_row = pd.DataFrame({
    'ImageName': ["Average"],
    'TP': [np.mean(TPs)],
    'FP': [np.mean(FPs)],
    'FN': [np.mean(FNs)],
    'Total_GT': [np.mean(GTs)],
    'Total_Pred': [np.mean(PRs)],
    'Precision': [np.mean(Precisions)],
    'Recall': [np.mean(Recalls)],
    'F1': [np.mean(F1s)],
    'DistanceThreshold_px': [dist_thresh]
})

df = pd.concat([df, avg_row], ignore_index=True)

# %% ================= SAVE TO EXCEL =================
if os.path.isfile(out_excel):
    # Read existing table
    df_old = pd.read_excel(out_excel)
    
    # Concatenate old and new tables
    df = pd.concat([df_old, df], ignore_index=True)

# Write table to Excel (index=False prevents pandas from writing row numbers)
df.to_excel(out_excel, index=False)

# %% ================= REPORT =================
print('\n=== Detection results ===')
print(f'Number of images: {n_images}')
print(f'Saved per-image metrics and average to: {out_excel}')