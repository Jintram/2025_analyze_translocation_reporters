

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


def segment_nucleus(image):
    thresh = threshold_otsu(image)
    binary_mask = image > thresh
    clean_mask = remove_small_objects(binary_mask, min_size=30)
    clean_mask = remove_small_holes(clean_mask, area_threshold=50)
    opened_mask = binary_opening(clean_mask, footprint=disk(2))
    
    # Label the mask to assign unique labels to each nucleus
    labeled_mask = label(opened_mask)
    
    return labeled_mask

def create_cytoplasm_roi(nucleus_mask, dilation_radius=5):
    # Create a binary mask for the nuclei
    binary_nucleus_mask = nucleus_mask > 0

    # Dilate the binary nucleus mask
    dilated_mask = binary_dilation(binary_nucleus_mask, footprint=disk(dilation_radius))

    # Create the cytoplasm ring by subtracting the nucleus mask from the dilated mask
    cytoplasm_ring = dilated_mask ^ binary_nucleus_mask

    # Assign the same labels as the nucleus mask to the cytoplasm ring
    cytoplasm_mask = np.where(cytoplasm_ring, nucleus_mask, 0)

    return cytoplasm_mask
