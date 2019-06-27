''' About the module:
Purpose: To subtract baselines in spectral lines.

Functions:
        baseline
'''


import numpy as np
import matplotlib.pyplot as plt


def baseline(v,img,vrange,order,model=False):
    '''
    Purpose: To subtract baselines in spectral lines.
    
    Inputs:
        v: 1D array, freq. axis array of img, can be in freq., chan., or velo.
        img: nD array, the spectral cube, must have shape (v,dec,ra)
        vrange: in 'v' unit, line-free ranges, [[v1,v2],[v3,v4],...]
        order: 0 or 1, constant or linear baseline
        model: bool, whether to output the baseline model

    Outputs:
        img_l: img with baseline subtracted
        img_b: baseline model data
    '''

    v = np.array(v)
    img = np.array(img)
    vrange = np.array(vrange)

    # line-free range
    s = img.shape
    x = np.arange(s[0])
    xrange = (vrange-vrange[0])/(vrange[-1]-vrange[0])*(s[0]-1)
    xrange = xrange.round().astype(np.int)

    # line-free data
    x_free = []
    img_free = []
    for xmin, xmax in xrange:
        x_free.append(np.arange(xmin,xmax+1))
        img_free.append(img[xmin:xmax+1])
    x_free = np.concatenate(tuple(x_free))    
    img_free = np.concatenate(tuple(img_free))

    # baselining
    img_b = {}
    if order==0:
        img_b = np.nanmean(img_free,axis=(0))
        img_b = img_b[np.newaxis,...]
    elif order==1:
        x_free_av = np.mean(x_free)
        img_free_av = np.nanmean(img_free,axis=0)
        b = np.nanmean(x_free[...,np.newaxis,np.newaxis]*img_free,axis=0)
        b -= x_free_av*img_free_av
        b /= (np.mean(x_free**2)-x_free_av**2)
        a = img_free_av-b*x_free_av
        img_b= b[np.newaxis,...]*x[...,np.newaxis,np.newaxis]+a[np.newaxis,...]
    else:
        print('Wrong order.')
    img_l = img - img_b

    # output
    if model:
        return img_l,img_b
    else:
        return img_l











