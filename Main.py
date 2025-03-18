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
# from Functions.Segmentation import (segment_nucleus, create_cytoplasm_roi)
import Functions.Segmentation as TRseg
    # import importlib; importlib.reload(TRseg)

# Tracking
# from Functions.Discarded.Cell_tracker import (extract_centroids, segment_and_extract_centroids, visualize_tracked_centroids)

# Measure
# from Functions.Intensity_measurements import (measure_intensities_for_all_timepoints, save_intensities_to_csv)
import Functions.Intensity_measurements as TRmeas
    # import importlib; importlib.reload(TRmeas)

import Functions.Plotting as TRplt
    # import importlib; importlib.reload(TRplt)

# More measure ?! (including multiple frames?)
# from Functions.Individual_measurements import (segment_and_extract_centroids, measure_cell_intensities, save_individual_intensities_to_csv)

################################################################################


# Enable interactive mode
plt.ion()
# Exit interactive mode
plt.ioff()

# Input and output folders
input_folder = "/Users/m.wehrens/Data_UVA/2024_10_Sebastian-KTR/202503_DATA_julian/Forskolin/"
output_folder = "/Users/m.wehrens/Data_UVA/2024_10_Sebastian-KTR/202503_OUTPUT-testmw/"
os.makedirs(output_folder, exist_ok=True)

# Channel information
NUCLEAR_CHANNEL = 1

# loop over tif files in input directory
# for each file, separately analyze and create a csv output file
for file_path in glob(os.path.join(input_folder, "*.tif")):
    # file_path = glob(os.path.join(input_folder, "*.tif"))[0]
    
    # get the filename, without the extension
    file_name = os.path.splitext(os.path.basename(file_path))[0]
    
    # read file
    image_stack = tiff.imread(file_path)
    print(f"Processing file: {file_path}")

    # num_timepoints, num_channels, height, width = image_stack.shape

    # segment the nuclei
    imgstack_nucleus = image_stack[:, NUCLEAR_CHANNEL]
    segmented_masks_preliminary = [TRseg.segment_nucleus(imgstack_nucleus[time_index]) for time_index in range(imgstack_nucleus.shape[0])]

    # For frames t>0, make the labeling consistent with frame t=0
    # The updated labeling is stored in segmented_masks_tracked.
    # To create new labels for a frame, updated information is necessary,
    # therefor, each frm+1 frame is constructed based segmented_masks_tracked[frm]
    # and segmented_masks_preliminary[frm+1].
    segmented_masks_tracked = [segmented_masks_preliminary[0]]
    for frm in range(len(segmented_masks_preliminary)-1):
        label_maskplus1 = TRseg.track_nuclei(segmented_masks_tracked[frm], segmented_masks_preliminary[frm+1])
        segmented_masks_tracked.extend([label_maskplus1])
    
    # Plot the tracking of the first N frames
    TRplt.plot_labels_framesX(segmented_masks_tracked, range_start=0, range_end=12, text_xoffset=50, output_folder=output_folder, file_name=file_name)
    
    # create the cytoplasmic regions (regions of interest, ROI)
    cytoplasm_rois = [TRseg.create_cytoplasm_roi(mask, dilation_radius=5) for mask in segmented_masks_tracked]

    # analyze the data
    # TO DO: MAKE THIS A LOOP OVER CUSTOM COLORS
    # TO DO: MAKE IT OUTPUT BOTH THE MEAN SIGNAL AND SINGLE CELLS
    cyan_channel = image_stack[:, 0]
    nucleus_intensities_cyan, cytoplasm_intensities_cyan = measure_intensities_for_all_timepoints(cyan_channel, segmented_masks_tracked, cytoplasm_rois)

    green_channel = image_stack[:, 2]
    nucleus_intensities_green, cytoplasm_intensities_green = measure_intensities_for_all_timepoints(green_channel, segmented_masks_tracked, cytoplasm_rois)

    timepoints = range(len(nucleus_intensities_cyan))
    output_csv = os.path.join(output_folder, os.path.splitext(os.path.basename(file_path))[0] + "_results.csv")

    # TO DO: MAKE THIS SAVE BOTH COLORS
    save_intensities_to_csv(nucleus_intensities_green, cytoplasm_intensities_green, timepoints, output_csv)


    
