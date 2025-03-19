################################################################################

# Translocation Reporters analysis script

# Expected tiff file input such that:
# image_stack.shape = num_timepoints, num_channels, height, width 

################################################################################
# import libraries
################################################################################

import os
import sys
from glob import glob

import numpy as np
import tifffile as tiff
import pandas as pd

################################################################################
# Import functions from other files in this repo
################################################################################

# sys.path.append("/Volumes/sils-mc/13776452/Python_scripts")
sys.path.append('/Users/m.wehrens/Documents/git_repos/_UVA/2024_small-projects/2025_analyze_translocation_reporters_Julian/')

import Functions.Segmentation as TRseg # Segmentation & tracking functions
    # import importlib; importlib.reload(TRseg)
import Functions.Intensity_measurements as TRmeas
    # import importlib; importlib.reload(TRmeas)
import Functions.Image_corrections as TRcorrect
    # import importlib; importlib.reload(TRcorrect)
import Functions.Plotting as TRplt
    # import importlib; importlib.reload(TRplt)

################################################################################
# Settings
################################################################################

# Enter/exit interactive mode
# plt.ion(); plt.ioff()
# plt.style.use("default")

# Read in settings from command
if (len(sys.argv) > 1):
    
    input_folder  = sys.argv[1]
    output_folder = sys.argv[2]
    AUTO_BACKGROUND_CORRECTION = bool(int(sys.argv[3]))
    
    # loop over remaining arguments
    MAPPING_CHANNELS = {}
    for idx in range(4, len(sys.argv)-1, 2):
        # print(sys.argv[idx], sys.argv[idx+1])
        MAPPING_CHANNELS[sys.argv[idx]] = int(sys.argv[idx+1])
    
    print(MAPPING_CHANNELS)

else:

    print('='*80)
    print('Please call this script as follows: \n')
    print('python analyze_transl_rep.py /input/folder/path/ /output/folder/path/ 0|1 nucleus 0 name1 1 name2 2\n')
    print('Where respectively folders can be customized, 0 or 1 is chosen to indicate auto background correction, ')
    print("and 'nucleus 0 ..' indicates in which channel nucleus and custom named channels to analyze can be found.\n")
    print('Exiting')
    print('='*80)
    sys.exit()

os.makedirs(output_folder, exist_ok=True)
    
if False:
    
    # Input and output folders
    input_folder = "/Users/m.wehrens/Data_UVA/2024_10_Sebastian-KTR/202503_DATA_julian/Forskolin/"
    # input_folder = "/Users/m.wehrens/Data_UVA/2024_10_Sebastian-KTR/202503_DATA_julian/testdata/"
    output_folder = "/Users/m.wehrens/Data_UVA/2024_10_Sebastian-KTR/202503_OUTPUT-testmw/"
    
    # SETTINGS
    MAPPING_CHANNELS = {'nucleus':0, 'ERK':1, 'PKA':2}
    nuclear_channel = MAPPING_CHANNELS['nucleus']
    AUTO_BACKGROUND_CORRECTION = False # only use this if there are areas in the picture with no signal


######################################################################
# Functions that constitute the loop below
# (More sophisticated functions are in other files)
######################################################################


def segment_and_track_nuclei(imgstack_nucleus, output_folder, file_name):
    # imgstack_nucleus = image_stack[:, nuclear_channel]
    # Note that some parameters below are defined implicitly by global values

    # segment the nuclei
    nucleus_masks_preliminary = [TRseg.segment_nucleus(imgstack_nucleus[time_index]) for time_index in range(imgstack_nucleus.shape[0])]

    # For frames t>0, make the labeling consistent with frame t=0
    # The updated labeling is stored in nucleus_masks_tracked.
    # To create new labels for a frame, updated information is necessary,
    # therefor, each frm+1 frame is constructed based nucleus_masks_tracked[frm]
    # and nucleus_masks_preliminary[frm+1].
    nucleus_masks_tracked = [nucleus_masks_preliminary[0]]
    for frm in range(len(nucleus_masks_preliminary)-1):
        label_maskplus1, _ = TRseg.track_nuclei(nucleus_masks_tracked[frm], nucleus_masks_preliminary[frm+1])
        nucleus_masks_tracked.extend([label_maskplus1])
    # Plot the nuclear segmentation and tracking of the first N frames
    TRplt.plot_labels_framesX(nucleus_masks_tracked, range_start=0, range_end=12, text_xoffset=50, output_folder=output_folder, file_name=file_name, suffix='_nuclei')
    
    return np.array(nucleus_masks_tracked), np.array(nucleus_masks_preliminary)


def create_cytoplasm_masks(nucleus_masks_tracked, output_folder, file_name, dilation_radius=5, margin_radius=0):
    
    cytoplasm_masks_tracked = [TRseg.create_cytoplasm_roi(mask, dilation_radius=dilation_radius, margin_radius=margin_radius) for mask in nucleus_masks_tracked] # tracked in is tracked out :)
    TRplt.plot_labels_framesX(cytoplasm_masks_tracked, range_start=0, range_end=12, text_xoffset=50, output_folder=output_folder, file_name=file_name, suffix='_cytorings')

    return np.array(cytoplasm_masks_tracked)


def calculate_intensity_values_to_df(MAPPING_CHANNELS, thekey, image_stack, nucleus_masks_tracked, cytoplasm_masks_tracked):
        
    image_stack_intensity = image_stack[:, MAPPING_CHANNELS[thekey]]
    
    # Optional, background correction
    if AUTO_BACKGROUND_CORRECTION:
        image_stack_intensity_corrected = np.array([TRcorrect.correct_background(img) for img in image_stack_intensity])
            # plt.imshow(image_stack_intensity[0]); plt.title('Not corrected'); plt.show(); plt.close()
            # plt.imshow(image_stack_intensity_corrected[0]); plt.title('Corrected'); plt.show(); plt.close()
    else:
        image_stack_intensity_corrected = image_stack_intensity
    
    # Determine the intensity alues for all timepoints
    df_current = \
        TRmeas.measure_intensities_for_all_timepoints(image_stack_intensity_corrected, nucleus_masks_tracked, cytoplasm_masks_tracked)
    
    # Calculate nucleus/cyto ratio
    df_current['Ratio_nucleus_div_cytoplasm'] = df_current['Intensity_nucleus']/df_current['Intensity_cytoplasm']
    # Add key to the df
    df_current['Key'] = thekey
    # Add sample filename to df
    df_current['Sample'] = file_name
        
    return df_current
        
        
######################################################################
# Main loop
################################################################################


# loop over tif files in input directory
# for each file, separately analyze and create a csv output file
df_list=[]
for file_path in glob(os.path.join(input_folder, "*.tif")):
    # file_path = glob(os.path.join(input_folder, "*.tif"))[0]
    
    # Read current file        
    print(f"Processing file: {file_path}")
    image_stack = tiff.imread(file_path)    
    file_name = os.path.splitext(os.path.basename(file_path))[0] # used further down

    # Segment the nuclei and track them such that labels are consistent throughout segmentation
    nucleus_masks_tracked, nucleus_masks_preliminary = segment_and_track_nuclei(image_stack[:, nuclear_channel], output_folder, file_name)

    # Create the cytoplasmic regions (regions of interest, ROI), and plot the rings of the first N frames
    cytoplasm_masks_tracked = create_cytoplasm_masks(nucleus_masks_tracked, output_folder, file_name, dilation_radius=5, margin_radius=0)

    # Now go over the channels and calculate the intensities (and ratios) for both
    keys_to_plot = [key for key in list(MAPPING_CHANNELS.keys()) if not (key=='nucleus')]
    for thekey in keys_to_plot: # thekey = keys_to_plot[1]
    
        df_current = calculate_intensity_values_to_df(MAPPING_CHANNELS, thekey, image_stack, nucleus_masks_tracked, cytoplasm_masks_tracked)
                
        # export current df
        df_current.to_csv(os.path.join(output_folder, f"{file_name}_{thekey}_results.csv"), index=False)
        df_current.to_excel(os.path.join(output_folder, f"{file_name}_{thekey}_results.xlsx"), index=False)
        # save current df to list
        df_list.append(df_current)


# concatenate all dfs
df_data_all = pd.concat(df_list, ignore_index=True)
# df_data_all = pd.concat(df_list[:2], ignore_index=True)
# Save those too
df_data_all.to_csv(os.path.join(output_folder, f"ALL_results.csv"), index=False)
df_data_all.to_excel(os.path.join(output_folder, f"ALL_results.xlsx"), index=False)

################################################################################
# Now create a plot of the signals
################################################################################

# Optionally, load data
# df_data_all = pd.read_csv(os.path.join(output_folder, f"ALL_results.csv"), index=False)

for CURRENT_SAMPLE in np.unique(df_data_all['Sample']):
    
    TRplt.plot_intensity_nuc_cyto(df_data_all.loc[df_data_all['Sample']==CURRENT_SAMPLE, ], output_folder, CURRENT_SAMPLE)
    TRplt.plot_intensity_ratio(df_data_all.loc[df_data_all['Sample']==CURRENT_SAMPLE, ], output_folder, CURRENT_SAMPLE)


################################################################################



