

import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

from skimage.measure import regionprops

import numpy as np

# generate a jet colormap, but make the first color white
my_jet_colors = np.concatenate([[[1, 1, 1, 1]], plt.cm.jet(np.linspace(0, 1, 256))])
jet_custom = ListedColormap(my_jet_colors)

def plot_nuclear_seg(segmented_masks, imgstack_nucleus):
    '''
    Plot the segmentation of the first 4 frames in a 
    2x2 panel plot.
    '''
    
    # plot the first three segmentation masks
    fig, ax = plt.subplots(2, 2, figsize=(10/2.54, 10/2.54))
    axf = ax.flatten()
    
    # now plot first three
    for idx, mask in enumerate(segmented_masks[:4]):
        #row, col = divmod(idx, 2)
        axf[idx].imshow(imgstack_nucleus[idx, ], cmap='viridis')
        axf[idx].contour(mask, bins=2, colors='white')
        axf[idx].set_title(f"Mask {idx+1}")
        axf[idx].axis('off')  
        
    plt.show()
    
def plot_nuclear_segmove(segmented_masks, imgstack_nucleus):
    '''
    Plots the segmentation mask outlines on top of each other
    for the starting time frames. Color coded for time using
    viridis.
    '''
    
    # plot the first three segmentation masks
    fig, ax = plt.subplots(1, 1, figsize=(10/2.54, 10/2.54))
    
    # create a series of hex colors based on viridis of length 10
    color_series = plt.cm.viridis(np.linspace(0, 1, 10))
    # now plot first three
    for idx, mask in enumerate(segmented_masks[:10]):
        ax.contour(mask, bins=2, colors=[color_series[idx]])
    
    # add a color bar
    sm = plt.cm.ScalarMappable(cmap='viridis', norm=plt.Normalize(vmin=0, vmax=9))
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, orientation='vertical', fraction=0.046, pad=0.04)
    cbar.set_label('Time Frame')  
    
    plt.show()
      
def plot_illustration_tracking(frame_t, frame_tplus1):
    # frame_t = segmented_masks_preliminary[2]; frame_tplus1 = segmented_masks_preliminary[3]
    
    fig, ax = plt.subplots(1, 2, figsize=(10/2.54, 10/2.54))
    plt.rcParams.update({'font.size': 12})
    
    _=ax[0].set_title('Frame t')
    _=ax[0].imshow(frame_t, cmap=jet_custom)
    _=ax[0].contour(frame_tplus1>0, bins=2, colors='black')
    
    for region in regionprops(frame_t):
        y0, x0 = region.centroid
        # center aligned
        _=ax[0].text(x0, y0, region.label, color='black',                      
                     bbox=dict(facecolor='white', alpha=0.3, edgecolor='none'), ha='center', va='center')  
    _=ax[1].set_title('Frame t+1')    
    _=ax[1].imshow(frame_tplus1, cmap=jet_custom)
    
    for region in regionprops(frame_tplus1):
        y0, x0 = region.centroid
        _=ax[1].text(x0, y0, region.label, color='black', 
                        bbox=dict(facecolor='white', alpha=0.3, edgecolor='none'), ha='center', va='center')
        
    # remove ticks and tick labels
    for a in ax:
        a.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
    
    plt.show()
    plt.close('all')
    

def plot_labels_framesX(labeled_masks, range_start=0, range_end=10, text_xoffset=20, output_folder=None, file_name=None):
    # labeled_masks=segmented_masks_tracked; range_start=1; range_end=3; text_xoffset=50
    '''
    Use jet-color coded imshow and text labels to 
    create a multipanel plot to show the labels
    of first 10 frames frames (0..9)
    '''
    
    # determine number of panels based on given range
    num_panels = range_end - range_start
    panel_rows = int(np.ceil(np.sqrt(num_panels)))
    panel_cols = int(np.ceil(num_panels / panel_rows))
    
    fig, ax = plt.subplots(panel_cols, panel_rows, figsize=(5*panel_rows/2.54, 5*panel_cols/2.54))
    plt.rcParams.update({'font.size': 6})
    axf = ax.flatten()
    
    for idx, frm in enumerate(range(range_start, range_end)):
        
        _=axf[idx].imshow(labeled_masks[frm], cmap=jet_custom)
        _=axf[idx].set_title(f"Frame {frm}")
        _=axf[idx].axis('off')
        
        for region in regionprops(labeled_masks[frm]):
            y0, x0 = region.centroid
            _=axf[idx].text(x0+text_xoffset, y0, region.label, color='black', 
                            bbox=dict(facecolor='white', alpha=0.3, edgecolor='none'))  
    # remove left-over panels
    for panel_idx in range(range_end, panel_rows*panel_cols):
        axf[panel_idx].axis('off')
    
    plt.tight_layout()
    
    if (output_folder is None) or (file_name is None):
        plt.show()
        plt.close(fig)
    else:
        plt.savefig(output_folder+file_name+'.pdf', dpi=300, bbox_inches='tight')
        plt.close(fig)
    
