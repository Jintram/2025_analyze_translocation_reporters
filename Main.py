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

import pandas as pd

import seaborn as sns

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
MAPPING_CHANNELS = {'nucleus':0, 'ERK':1, 'PKA':2}
nuclear_channel = MAPPING_CHANNELS['nucleus']

# TO DO: THE LOOP BELOW IS TOO LONG, PERHAPS WE SHOULD SPLIT THE CODE INTO FUNCTIONS

# loop over tif files in input directory
# for each file, separately analyze and create a csv output file
df_list=[]
for file_path in glob(os.path.join(input_folder, "*.tif")):
    # file_path = glob(os.path.join(input_folder, "*.tif"))[0]
    
    # get the filename, without the extension
    file_name = os.path.splitext(os.path.basename(file_path))[0]
    
    # read file
    image_stack = tiff.imread(file_path)
    print(f"Processing file: {file_path}")

    # num_timepoints, num_channels, height, width = image_stack.shape

    # segment the nuclei
    imgstack_nucleus = image_stack[:, nuclear_channel]
    nucleus_masks_preliminary = [TRseg.segment_nucleus(imgstack_nucleus[time_index]) for time_index in range(imgstack_nucleus.shape[0])]

    # For frames t>0, make the labeling consistent with frame t=0
    # The updated labeling is stored in nucleus_masks_tracked.
    # To create new labels for a frame, updated information is necessary,
    # therefor, each frm+1 frame is constructed based nucleus_masks_tracked[frm]
    # and nucleus_masks_preliminary[frm+1].
    nucleus_masks_tracked = [nucleus_masks_preliminary[0]]
    for frm in range(len(nucleus_masks_preliminary)-1):
        label_maskplus1 = TRseg.track_nuclei(nucleus_masks_tracked[frm], nucleus_masks_preliminary[frm+1])
        nucleus_masks_tracked.extend([label_maskplus1])
    # Plot the nuclear segmentation and tracking of the first N frames
    TRplt.plot_labels_framesX(nucleus_masks_tracked, range_start=0, range_end=12, text_xoffset=50, output_folder=output_folder, file_name=file_name, suffix='_nuclei')
    
    # Create the cytoplasmic regions (regions of interest, ROI)
    cytoplasm_masks_tracked = [TRseg.create_cytoplasm_roi(mask, dilation_radius=5) for mask in nucleus_masks_tracked] # tracked in is tracked out :)
    # Plot the rings of the first N frames
    TRplt.plot_labels_framesX(cytoplasm_masks_tracked, range_start=0, range_end=12, text_xoffset=50, output_folder=output_folder, file_name=file_name, suffix='_cytorings')

    # TO DO: DO WE WANT SOME KIND OF BACKGROUND CORRECTION?

    # Now combine the segmentation with the intensity signals to calculate nuclear and cytoplasmic signals per cell
    # nucleus_masks_tracked, cytoplasm_masks_tracked, image_stack_current    
    for thekey in list(MAPPING_CHANNELS.keys()):        
        # thekey=list(MAPPING_CHANNELS.keys())[1]
        
        if thekey == 'nucleus': continue
        
        image_stack_intensity = image_stack[:, MAPPING_CHANNELS[thekey]]
        
        df_current = \
            TRmeas.measure_intensities_for_all_timepoints(image_stack_intensity, nucleus_masks_tracked, cytoplasm_masks_tracked)
        
        # Calculate nucleus/cyto ratio
        df_current['Ratio_nucleus_div_cytoplasm'] = df_current['Intensity_nucleus']/df_current['Intensity_cytoplasm']
        # Add key to the df
        df_current['Key'] = thekey
        # Add sample filename to df
        df_current['Sample'] = file_name
        
        # export current df
        df_current.to_csv(os.path.join(output_folder, f"{file_name}_{thekey}_results.csv"), index=False)
        df_current.to_excel(os.path.join(output_folder, f"{file_name}_{thekey}_results.xlsx"), index=False)
        
        # save current df to list
        df_list.append(df_current)
        
# concatenate all dfs
df_data_all = pd.concat(df_list, ignore_index=True)
# Save those too
df_data_all.to_csv(os.path.join(output_folder, f"{file_name}_ALL_results.csv"), index=False)
df_data_all.to_excel(os.path.join(output_folder, f"{file_name}_ALL_results.xlsx"), index=False)

################################################################################
# Now create a plot of the signals

# Optionally, simply load data
# (This might be usefull if we split these two 
# df_data_all = pd.read_csv(os.path.join(output_folder, f"{file_name}_ALL_results.csv"), index=False)

# TO DO: CONVERT THE BELOW CODE TO FUNCTIONS

# plot the data
# df_data_all.columns
sns.set_theme(style="whitegrid")
g = sns.FacetGrid(df_data_all, col="Key", col_wrap=2, height=4, sharey=False)
# Plot individual cells with hue='Cell'
_=g.map(sns.lineplot, "Frame", "Intensity_nucleus", data=df_data_all.loc[df_data_all['Cell'] != 'all'], hue='Cell', legend=False)
_=g.map(sns.lineplot, "Frame", "Intensity_cytoplasm", linestyle='--', data=df_data_all.loc[df_data_all['Cell'] != 'all'], hue='Cell', legend=False)
# Plot the average (black lines)
_=g.map(sns.lineplot, "Frame", "Intensity_nucleus", data=df_data_all.loc[df_data_all['Cell'] == 'all'], color='black', units='Cell', estimator=None, linewidth=2, label='Average Nucleus')
_=g.map(sns.lineplot, "Frame", "Intensity_cytoplasm", data=df_data_all.loc[df_data_all['Cell'] == 'all'], color='black', units='Cell', estimator=None, linewidth=2, linestyle='--', label='Average Cytoplasm')
# Add legend 
g.add_legend()
# Save
# plt.tight_layout()
g.figure.savefig(os.path.join(output_folder, f"PLOT_{file_name}_Intensity_plot_nuc-cyto-separate.pdf"), dpi=300, bbox_inches='tight')
plt.close(g.figure)


# plot the data
# df_data_all.columns
sns.set_theme(style="whitegrid")
g = sns.FacetGrid(df_data_all, col="Key", col_wrap=2, height=4, sharey=False)
g.map(sns.lineplot, "Frame", "Ratio_nucleus_div_cytoplasm", data=df_data_all.loc[df_data_all['Cell']=='all'], color='black', units='Cell', estimator=None, linewidth=2)
g.map(sns.lineplot, "Frame", "Ratio_nucleus_div_cytoplasm", data=df_data_all.loc[df_data_all['Cell']!='all'], hue='Cell')
g.set_axis_labels("Time", "Intensity")
g.add_legend()
# Save
# plt.tight_layout()
g.figure.savefig(os.path.join(output_folder, f"PLOT_{file_name}_Intensity_plot_nuc-cyto-ratio.pdf"), dpi=300, bbox_inches='tight')
plt.close(g.figure)









    
