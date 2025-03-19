

# import ndimage
from scipy import ndimage
import numpy as np

def get_background_mask(img_int, ESTIMATED_OBJECT_RADIUS=30):
    # img_int= image_stack_intensity[0];ESTIMATED_OBJECT_RADIUS=30
    '''
    Automatically determine the location of a background based on 
    a dilation and search for minimum.
    Both the size of the background box and the dilatino are based
    on the estimated object radius.
    This should be the objects that appear in the data.
    '''
    
    background_mask_size_odd = ESTIMATED_OBJECT_RADIUS + 1 if ESTIMATED_OBJECT_RADIUS % 2 == 0 else ESTIMATED_OBJECT_RADIUS
    background_mask_halfsize = background_mask_size_odd//2
    
    # apply a large min filter, size of nucleus
    img_int_max = ndimage.maximum_filter(img_int, size=background_mask_size_odd)
        # plt.imshow(img_int_max); plt.show(); plt.close()
    
    # mask with lowest values
    low_mask = img_int_max==np.min(img_int_max)
        
    # distance transform the mask
    mask_distance = ndimage.distance_transform_edt(low_mask)
    
    # get the location of the mask_distance's maximum
    # cut some margins using ESTIMATED_OBJECT_RADIUS, such that the final box doesn't go outside the image
    mask_distance_margincut = mask_distance[background_mask_halfsize:-background_mask_halfsize, background_mask_halfsize:-background_mask_halfsize]
    max_loc = np.array(np.unravel_index(np.argmax(mask_distance_margincut), mask_distance_margincut.shape))+background_mask_halfsize
        
    # produce the boundaries of a box around the maximum (±ESTIMATED_OBJECT_RADIUS), similarly formatted as the bbox in regionprops
    background_bbox_coords = (max_loc[0]-background_mask_halfsize, max_loc[1]-background_mask_halfsize,
                                max_loc[0]+background_mask_halfsize, max_loc[1]+background_mask_halfsize)

    # produce a mask according to the bbox
    background_mask = np.zeros(img_int.shape, dtype=bool)
    background_mask[background_bbox_coords[0]:background_bbox_coords[2], background_bbox_coords[1]:background_bbox_coords[3]] = True
      
    if False:
        # Show 
        plt.imshow(img_int_max)
        plt.imshow(low_mask, alpha=0.5, cmap='gray')
        plt.imshow(mask_distance, alpha=0.5, cmap='gray')    
        plt.scatter(max_loc[1], max_loc[0], color='r', s=100)
        # draw the bbox
        plt.plot([background_bbox_coords[1], background_bbox_coords[3], background_bbox_coords[3], background_bbox_coords[1], background_bbox_coords[1]],\
                    [background_bbox_coords[0], background_bbox_coords[0], background_bbox_coords[2], background_bbox_coords[2], background_bbox_coords[0]], color='r')             
        # draw the mask
        # plt.contour(background_mask, levels=[0.5], colors='lightblue')
        plt.grid(False)
        plt.show(); plt.close()
    
    return background_mask

def correct_background(img_int, ESTIMATED_OBJECT_RADIUS=30):
    
    # Determine which part of the image can be considered background
    background_mask = get_background_mask(img_int, ESTIMATED_OBJECT_RADIUS=30)
    
    # Determine a background value
    median_background = np.median(img_int[background_mask])
    
    # Return an image with the background subtracted
    img_int_corrected = img_int.copy()
    img_int_corrected[img_int_corrected<median_background] = median_background # make sure no negative values result next line
    img_int_corrected = img_int_corrected - median_background
        # img_int_corrected[background_mask]
    
    return img_int_corrected
    
    
    