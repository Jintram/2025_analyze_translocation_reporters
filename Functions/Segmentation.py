

import numpy as np
import matplotlib.pyplot as plt
from PIL import Image
import tifffile as tiff
from skimage.filters import threshold_otsu
from skimage.morphology import remove_small_objects, remove_small_holes, binary_opening, disk, binary_dilation, dilation
import os
from glob import glob
import csv
from skimage.measure import label


def segment_nucleus(image, min_size_objects=30,  area_threshold_holes=50, footprint_opening = 2):
    '''
    Segment an image with nuclei based on simple Otsu thresholding.
    Make small improvements to the image by removing small objects and holes, and opening.
    Returns labeled mask.
    '''
    
    # Perform thresholding
    thresh = threshold_otsu(image)
    binary_mask = image > thresh
    
    # Remove small objects and holes
    clean_mask = remove_small_objects(binary_mask, min_size=min_size_objects)
    clean_mask = remove_small_holes(clean_mask, area_threshold=area_threshold_holes)
    
    # Slight opening to ±smooth image
    opened_mask = binary_opening(clean_mask, footprint=disk(footprint_opening))
    
    # Label the mask to assign unique labels to each nucleus
    labeled_mask = label(opened_mask)
    
    return labeled_mask

def create_cytoplasm_roi(nucleus_mask, dilation_radius=5, margin_radius = 0):
    '''
    Create cytoplasm ring by dilating nuclear mask
    and removing the original nuclei again.
    
    If desired, the part that's removed can also be dilated, such that 
    a margin is introduced between the original nuclei and the cytoplasm ring.
    
    Note that if two nuclei are very close, the dilated mask pixels close to 
    both nuclei will have the value of the highest label. This is arbitrary
    and might introduce small artifacts. This is not expected to have
    big effects on measurements.
    '''
    # nucleus_mask = segmented_masks[0]; dilation_radius=5; margin_radius = 0
    
    # plt.imshow(nucleus_mask); plt.show(); plt.close()

    # Dilate the binary nucleus mask
    dilated_mask = dilation(nucleus_mask, footprint=disk(dilation_radius+margin_radius))
        # plt.imshow(dilated_mask); plt.show(); plt.close()

    # now remove the areas of the original nuclei, if desired including a margin
    dilated_nucleus_mask = dilation(nucleus_mask, footprint=disk(margin_radius))
    dilated_mask[dilated_nucleus_mask > 0] = 0
        # plt.imshow(dilated_mask); plt.show(); plt.close()

    return dilated_mask
