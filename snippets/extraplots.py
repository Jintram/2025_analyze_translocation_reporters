''' Additional plots to check a few things. '''




import matplotlib.pyplot as plt




def plotaframe(image_stack,nucleus_masks_tracked,t,ch):
    # t = 0; ch = 2
    
    plt.imshow(image_stack[t,ch,:,:])
    plt.contour(nucleus_masks_tracked[t, :,:], colors='r', linewidths=0.5)
    plt.show()
    plt.close()