################################################################################

# Translocation Reporters analysis script


################################################################################
# import libraries

import numpy as np
import matplotlib.pyplot as plt
from PIL import Image
import tifffile as tiff
from skimage.filters import threshold_otsu
from skimage.morphology import remove_small_objects, remove_small_holes, binary_opening, disk, binary_dilation
import os
from glob import glob
import csv
import sys

################################################################################
# import custom libs

# sys.path.append("/Volumes/sils-mc/13776452/Python_scripts")
sys.path.append('/Users/m.wehrens/Documents/git_repos/_UVA/2024_small-projects/2025_analyze_translocation_reporters_Julian/')

# Segment
from Functions.Segmentation import (segment_nucleus, create_cytoplasm_roi)

# Tracking
from Functions.Discarded.Cell_tracker import (extract_centroids, segment_and_extract_centroids, visualize_tracked_centroids)

# Measure
from Functions.Intensity_measurements import (measure_intensities_for_all_timepoints, save_intensities_to_csv)

# More measure ?! (including multiple frames?)
from Functions.Individual_measurements import (segment_and_extract_centroids, measure_cell_intensities, save_individual_intensities_to_csv)

################################################################################


# Enable interactive mode
plt.ion()
# Exit interactive mode
plt.ioff()

# Input and output folders
input_folder = "/Users/m.wehrens/Data_UVA/2024_10_Sebastian-KTR/202503_DATA_julian/Forskolin/"
output_folder = "/Users/m.wehrens/Data_UVA/2024_10_Sebastian-KTR/202503_OUTPUT-testmw/"
os.makedirs(output_folder, exist_ok=True)

# loop over tif files in input directory
# for each file, separately analyze and create a csv output file
for file_path in glob(os.path.join(input_folder, "*.tif")):
    # file_path = glob(os.path.join(input_folder, "*.tif"))[0]
    
    # read file
    image_stack = tiff.imread(file_path)
    print(f"Processing file: {file_path}")

    num_timepoints, num_channels, height, width = image_stack.shape
    time_index = 0
    original_image = image_stack[:, 0][time_index]
    segmented_mask = segment_nucleus(original_image)
        # plt.imshow(original_image); plt.show(); 

    nucleus_channel = image_stack[:, 1]
    segmented_masks = [segment_nucleus(nucleus_channel[time_index]) for time_index in range(nucleus_channel.shape[0])]
    cytoplasm_rois = [create_cytoplasm_roi(mask, dilation_radius=5) for mask in segmented_masks]

    cyan_channel = image_stack[:, 0]
    nucleus_intensities_cyan, cytoplasm_intensities_cyan = measure_intensities_for_all_timepoints(cyan_channel, segmented_masks, cytoplasm_rois)

    green_channel = image_stack[:, 2]
    nucleus_intensities_green, cytoplasm_intensities_green = measure_intensities_for_all_timepoints(green_channel, segmented_masks, cytoplasm_rois)

    timepoints = range(len(nucleus_intensities_cyan))
    output_csv = os.path.join(output_folder, os.path.splitext(os.path.basename(file_path))[0] + "_results.csv")

    save_intensities_to_csv(nucleus_intensities_green, cytoplasm_intensities_green, timepoints, output_csv)

    image_stack = tiff.imread(file_path)
    print(f"Processing file: {file_path}")

    tracking_data = segment_and_extract_centroids(image_stack)
    # print(f"Tracking data for {file_path}:", tracking_data)

    #visualize_tracked_centroids(image_stack, tracking_data)
    cell_intensity_data = measure_cell_intensities(image_stack, tracking_data)
    # print(f"Intensity data for {file_path}:", cell_intensity_data)

    for label, intensities in cell_intensity_data.items():
        for frame, (nucleus_intensity, cytoplasm_intensity) in enumerate(zip(intensities["nucleus"], intensities["cytoplasm"])):
            # print(f"Label {label}, Frame {frame}: Nucleus Intensity = {nucleus_intensity}, Cytoplasm Intensity = {cytoplasm_intensity}")  # Commented out
            pass

    output_csv = os.path.join(output_folder, os.path.splitext(os.path.basename(file_path))[0] + "_intensities.csv")
    save_individual_intensities_to_csv(cell_intensity_data, output_csv)
    print(f"Intensity data saved to {output_csv}")
    print(f"Results saved to {output_csv}")
