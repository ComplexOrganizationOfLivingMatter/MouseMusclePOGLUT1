import os
import glob
import random
import numpy as np
import warnings
from skimage import io

# %% Parameters
crop_size = 1020
n_pairs_to_crop = 300
data_path = os.path.join('..', '..', 'data')

# %% Paths
# Recursive globbing (**) in Python requires the recursive=True flag
annotated_masks_pattern = os.path.join(data_path, 'annotatedMuscleMasks', '**', '*.tif')
path_tifs = glob.glob(annotated_masks_pattern, recursive=True)

path_filtered_mask = os.path.join('..', '..', 'results', 'binaryNucleiStemCells_filtered')
path_raw = os.path.join(data_path, 'rawImages', 'tifsHistMatched')

# Load files
mask_files = glob.glob(os.path.join(path_filtered_mask, '*.tif'))
raw_files = glob.glob(os.path.join(path_raw, '**', '*.tif'), recursive=True)

path2save_dataset = os.path.join(data_path, 'trainingDataset', 'individual label', 'Pax7_v2', 'backup')
output_train_raw = os.path.join(path2save_dataset, 'train', 'raw')
output_train_mask = os.path.join(path2save_dataset, 'train', 'mask')
output_test_raw = os.path.join(path2save_dataset, 'test', 'raw')
output_test_mask = os.path.join(path2save_dataset, 'test', 'mask')

# Create output folders if they don't exist
folders = [output_train_raw, output_train_mask, output_test_raw, output_test_mask]
for f in folders:
    os.makedirs(f, exist_ok=True)

# %% Prepare list of matching mask–raw pairs
# We map filenames to their full paths using dictionaries. This avoids the O(N) 
# linear search MATLAB does with `strcmp`, making lookups instantaneous.
raw_dict = {os.path.basename(f): f for f in raw_files}
mask_dict = {os.path.basename(f): f for f in mask_files}
selected_names = [os.path.basename(f) for f in path_tifs]

valid_pairs = []
for name in selected_names:
    if name in mask_dict and name in raw_dict:
        valid_pairs.append(name)

print(f"Found {len(valid_pairs)} valid image pairs.")
if not valid_pairs:
    raise ValueError("No matching mask–raw image pairs found.")

# Checks if two rectangles [x, y, w, h] overlap
def rects_overlap(r1, r2):
    return not (r1[0] >= r2[0] + r2[2] or  # r1 is to the right of r2
                r1[0] + r1[2] <= r2[0] or  # r1 is to the left of r2
                r1[1] >= r2[1] + r2[3] or  # r1 is below r2
                r1[1] + r1[3] <= r2[1])    # r1 is above r2

# Random cropping and saving
n_crops_saved = 0
max_attempts = n_pairs_to_crop * 10
attempt = 0

# crop_history maps filenames to a list of bounding boxes: [x, y, width, height]
crop_history = {} 

while n_crops_saved < n_pairs_to_crop and attempt < max_attempts:
    attempt += 1
    
    # Randomly pick a valid pair
    file_name = random.choice(valid_pairs)
    
    # Load mask and raw
    mask_path = mask_dict[file_name]
    raw_path = raw_dict[file_name]
    
    mask = io.imread(mask_path)
    raw = io.imread(raw_path)
    
    # Check image size
    h, w = mask.shape[:2]
    if h < crop_size or w < crop_size:
        continue
        
    # Random crop coordinates (Python is 0-indexed)
    # random.randint(a, b) includes both a and b.
    x = random.randint(0, w - crop_size)
    y = random.randint(0, h - crop_size)
    new_rect = [x, y, crop_size, crop_size]
    
    # Check overlap with previous crops of the same image
    if file_name in crop_history:
        prev_rects = crop_history[file_name]
        # If any previous rectangle overlaps with the new one, skip
        if any(rects_overlap(pr, new_rect) for pr in prev_rects):
            continue
    else:
        crop_history[file_name] = []
        
    # Crop the images
    mask_crop = mask[y:y+crop_size, x:x+crop_size]
    
    # Handle both grayscale (2D) and multi-channel (3D) raw arrays properly
    if raw.ndim > 2:
        raw_crop = raw[y:y+crop_size, x:x+crop_size, :]
    else:
        raw_crop = raw[y:y+crop_size, x:x+crop_size]
        
    # Skip empty mask crops (checks if all values are 0)
    if not np.any(mask_crop):
        continue
        
    # Save crop
    n_crops_saved += 1
    
    if n_crops_saved <= 250:
        out_mask = os.path.join(output_train_mask, f'mask_{n_crops_saved:03d}.tif')
        out_raw = os.path.join(output_train_raw, f'raw_{n_crops_saved:03d}.tif')
    else:
        out_mask = os.path.join(output_test_mask, f'mask_{n_crops_saved-250:03d}.tif')
        out_raw = os.path.join(output_test_raw, f'raw_{n_crops_saved-250:03d}.tif')
    
    # Save images (suppressing low contrast warnings for empty/sparse masks)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        io.imsave(out_mask, mask_crop)
        io.imsave(out_raw, raw_crop)
        
    # Update crop history
    crop_history[file_name].append(new_rect)
    print(f"Saved crop {n_crops_saved}/{n_pairs_to_crop} from {file_name}")

if n_crops_saved < n_pairs_to_crop:
    print(f"⚠️ Warning: Only {n_crops_saved} crops were generated (insufficient non-empty regions).")
else:
    print("✅ Dataset generation completed successfully.")