

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


def track_nuclei(mask_t, mask_tplus1):
    # mask_t = nucleus_masks_preliminary[0]; mask_tplus1 = nucleus_masks_preliminary[1]
    # problem case:
    # mask_t = nucleus_masks_preliminary[2]; mask_tplus1 = nucleus_masks_preliminary[3]
    
    ''''
    label propagation
    
    Given two labeled masks, which contain similar segmentations, but one at timepoints t, 
    and one at timepoint t+1, the goal of this function is to make the labeling consistent 
    between the masks of the two timepoints.

    The function will take as input two labeled masks, mask_t and mask_tplus1, go over 
    each of the labels within mask_t, select the corresponding pixels in mask_tplus1, 
    and determine which label is most often found in those pixels of mask_tplus1.
    
    Note 1: this strategy might lead to lineages ending in 0, leaving that object untracked
    in all subsequent frames.
    
    Note 2: a more sophisticated approach would be to calculate an overlap matrix ij between
    labels i from mask_t and labels j from mask_tplus1, but this would be more computationally
    expensive, and the result would be the same. (The advantage of that would be to more 
    easily identify problem cases with multiple overlaps.)
    '''
    
    # obtain non-zero labels from mask_t
    mask_t_labels = np.unique(mask_t)
    mask_t_labels = mask_t_labels[mask_t_labels>0]
    
    # 
    mask_tplus1_corrected = np.zeros_like(mask_tplus1)
    
    # now for each label check with which label in t+1 they overlap most
    the_mapping = {}
    for lbl in mask_t_labels:
        # lbl=2
        
        overlapping_labels = mask_tplus1[mask_t == lbl]
            # test case
            # overlapping_labels = np.array([0, 10, 10, 10, 2, 10, 10, 10, 3, 3, 0 , 0, 0, 0, 0, 0, 0])        
            
        # get the label that's most frequent non-zero value (mode)            
        # first, determine the frequency of the labels in the corresponding area in the previous frame
        label_frequency = np.bincount(overlapping_labels)
        # if non-zero values found, take the most frequent value (mode)
        if len(label_frequency)>1:
            mode_label = label_frequency[1:].argmax()+1 # remove 0; hence also +1
        # if no non-zero labels are found, make it zero (will be ignored later)
        else:
            mode_label = 0
                
        # now assing that at the correct position in the updated mask
        # we have determined mode_label to be the matching label in the t+1 frame
        # this matches lbl in the frame t. 
        # so we need to set all positions of the mode_label to the label
        # consistent with the frame t, ie lbl.
        if (not mode_label==0):
            mask_tplus1_corrected[mask_tplus1==mode_label] = lbl
            
            # save how labels are updated
            the_mapping[mode_label] = lbl # currently for debugging, could be used also
    
    return mask_tplus1_corrected, the_mapping