

import matplotlib.pyplot as plt



def plot_nuclear_seg(segmented_masks, imgstack_nucleus):
    
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