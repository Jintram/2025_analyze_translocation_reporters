

import numpy as np
import matplotlib.pyplot as plt
from PIL import Image
import tifffile as tiff
from skimage.filters import threshold_otsu
from skimage.morphology import remove_small_objects, remove_small_holes, binary_opening, disk, binary_dilation
import os
from glob import glob
import csv
from skimage.measure import label

import pandas as pd


# Function to visualize a specific time point and channel
def visualize_timepoint(stack, time_index, channel_index):
    plt.figure(figsize=(6, 6))
    plt.title(f"Time Point: {time_index}, Channel: {channel_index}")
    plt.imshow(stack[time_index, channel_index], cmap="gray")
    plt.axis("off")
    plt.show(block=False)
    plt.waitforbuttonpress()  # Wait for a button press to continue


# Function to measure mean intensity in a given mask (nucleus or cytoplasm)
def measure_intensity(image, mask):
    return np.mean(image[mask])

def thedata_rowasdict(time_index, cell_lbl, nucleus_intensity, cytoplasm_intensity):
    
    return {
        "Frame": time_index,
        "Cell": str(cell_lbl),
        "Intensity_nucleus": nucleus_intensity,
        "Intensity_cytoplasm": cytoplasm_intensity,
        }
    

# Function to measure intensities for both nucleus and cytoplasm for all time points
def measure_intensities_for_all_timepoints(image_stack_intensity, nucleus_masks_tracked, cytoplasm_masks_tracked):
    
    # A list to store the data
    # Each entry will be a dictionary (see function thedata_rowasdict above)
    # {'Frame': .., 'Cell': .., 'Intensity_nucleus': .., 'Intensity_cytoplasm': ..}
    thedata = []
    
    # Go over frames
    for time_index in range(image_stack_intensity.shape[0]):
        
        current_image = image_stack_intensity[time_index]
        
        # Calculate the mean intensities per cell                
        for cell_lbl in range(1, np.max(nucleus_masks_tracked)+1):
            
            # Calculate mean intensity for this cell
            nucleus_intensity = measure_intensity(current_image, nucleus_masks_tracked[time_index]==cell_lbl)
            cytoplasm_intensity = measure_intensity(current_image, cytoplasm_masks_tracked[time_index]==cell_lbl)        
            
            # Add information to list
            thedata.append(thedata_rowasdict(time_index, cell_lbl, nucleus_intensity, cytoplasm_intensity))
            
        # Also calculate the overall mean intensities, add to list
        thedata.append(thedata_rowasdict(time_index, 'all', 
                                         measure_intensity(current_image, nucleus_masks_tracked[time_index]>0), 
                                         measure_intensity(current_image, cytoplasm_masks_tracked[time_index]>0)))

    # Convert to dataframe    
    df_intensities = pd.DataFrame(thedata)
    
    return df_intensities



# Function to save intensities to CSV
def save_intensities_to_csv(nucleus_intensities, cytoplasm_intensities, timepoints, output_csv):
    with open(output_csv, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(["Timepoint", "Nucleus", "Cytoplasm"])
        for timepoint, nucleus_intensity, cytoplasm_intensity in zip(timepoints, nucleus_intensities, cytoplasm_intensities):
            writer.writerow([timepoint, nucleus_intensity, cytoplasm_intensity])
