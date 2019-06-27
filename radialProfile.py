import numpy as np

def azimuthalAverage(img,center=None):
    """
    Calculate the azimuthally averaged radial profile.

    img - The 2D image
    center - The [x,y] pixel coordinates used as the center. The default is 
             None, which then uses the center of the image (including 
             fracitonal pixels).
    
    """


    # Calculate the indices from the image
    y, x = np.indices(img.shape)

    if not center:
        center = np.array([(x.max()-x.min())/2.0, (y.max()-y.min())/2.0])

    r = np.hypot(x-center[0],y-center[1])

    # flat arrays
    img_flat = img.flatten()
    r_flat = r.flatten() 

    # drop nan's
    ind_valid = ~np.isnan(img_flat)
    img_flat = img_flat[ind_valid]
    r_flat = r_flat[ind_valid]

    # Get sorted radii
    ind = np.argsort(r_flat)
    r_sorted = r_flat[ind]
    i_sorted = img_flat[ind]

    # Get the integer part of the radii (bin size = 1)
    r_int = r_sorted.astype(int)

    # Find all pixels that fall within each radial bin.
    deltar = r_int[1:] - r_int[:-1]  # Assumes all radii represented
    rind = np.where(deltar)[0]       # location of changed radius
    nr = rind[1:] - rind[:-1]        # number of radius bin
    
    # Cumulative sum to figure out sums for each radius bin
    csim = np.cumsum(i_sorted,dtype=float)
    tbin = csim[rind[1:]] - csim[rind[:-1]]

    radial_prof = tbin / nr

    return radial_prof























