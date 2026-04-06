import os
import glob
import numpy as np
import matplotlib.pyplot as plt
from skimage import io
from skimage.exposure import match_histograms

# Histogram Matching function
# as reference image we will consider KIKO_1_7um

dir_tifs = os.path.join('..', 'data', 'rawImages', 'tifs')
days_folder = 'p4'

# Read the reference image
ref_path = os.path.join(dir_tifs, 'Initial dataset', 'KIKO 1 p20.tif')
img_reference = io.imread(ref_path)

# Find all images matching the '*.tif' pattern
search_pattern = os.path.join(dir_tifs, days_folder, '*.tif')
list_images = glob.glob(search_pattern)

# Define and create the target directory if it doesn't exist
target_dir = os.path.join(dir_tifs, '..', 'tifsHistMatched', days_folder)
os.makedirs(target_dir, exist_ok=True) 

# Loop through the matching images
for img_path in list_images:
    print(img_path)
    
    # Extract just the file name (equivalent to listImages(i).name)
    img_name = os.path.basename(img_path)
    
    # Read the target image
    img2HM = io.imread(img_path)
    
    # Perform histogram matching
    # channel_axis is set conditionally to handle both RGB/multichannel and grayscale images
    is_multichannel = img2HM.ndim > 2
    matched_image = match_histograms(img2HM, img_reference, 
                                     channel_axis=-1 if is_multichannel else None)
    
    # Display the original and matched image side-by-side 
    combined_image = np.hstack((img2HM, matched_image))
    
    plt.imshow(combined_image, cmap='gray' if not is_multichannel else None)
    plt.axis('off')
    plt.show(block=False)
    plt.pause(0.5) # Pauses for half a second to display before closing
    plt.close()
    
    # Save the matched image
    save_path = os.path.join(target_dir, img_name)
    io.imsave(save_path, matched_image)