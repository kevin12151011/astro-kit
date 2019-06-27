'''
mask module:
    To create mask arrays.

Features:
    1. x,y are 2D array, in case of non-Euclidean mapping.
'''

import numpy as np
from matplotlib.path import Path



def converter(a,yes,no):
    '''
    To modify the values in array a.

    Inputs
        a: bool array
        yes: set a[a] = yes
        no: set a[~a] = no

    Outputs
        a1: the modified array
    '''

    a1 = {}
    if type(yes)==bool and type(no)==bool:
        a1 = a.copy()
    else:
        a1 = np.zeros(a.shape,dtype=float)
        a1[a] = yes
        a1[~a] = no 

    return a1


def circle(X,Y,x0,y0,r,yes=True,no=False):
    '''
    To create round mask array.

    Inputs:
        X, Y: 2D array, coordinates of the points, shape=(ny,nx).
        x0, y0: center of the circle
        r: radius of the circle
        yes: bool or float, value to fill pixels which meet the condition
        no: bool or float, value to fill pixels which don't meet the condition
    
    Outputs:
        mask
    '''

    R = ((X-x0)**2+(Y-y0)**2)**.5
    mask = np.full(R.shape,False)
    mask[R<=r] = True

    return converter(mask,yes,no)


def ellipse(X,Y,x0,y0,a,b,theta,yes=True,no=False):
    '''
    To create mask array.

    Inputs:
        X, Y: 2D array, coordinates of the points, shape=(ny,nx).
        x0, y0: center of the ellipse
        a, b: major and minor axis of the ellipse
        theta: tilt angle, between a- & x-axis, in deg
        yes: bool or float, value to fill pixels which meet the condition
        no: bool or float, value to fill pixels which don't meet the condition

    Outputs:
        mask
    '''
    
    theta *= np.pi/180
    X1 = (X-x0)*np.cos(theta) + (Y-y0)*np.sin(theta)
    Y1 = -(X-x0)*np.sin(theta)+ (Y-y0)*np.cos(theta)
    R = (X1/a)**2 + (Y1/b)**2
    mask = np.full(R.shape,False)
    mask[R<=1] = True
    
    return converter(mask,yes,no)


def polygon(X,Y,coords,yes=True,no=False):
    '''
    To create mask array.

    Inputs:
        X, Y: 2D array, coordinates of the points, shape=(ny,nx).
        coords: N*2 array, coords of the vertices of the polygon [x,y]
        yes: bool or float, value to fill pixels which meet the condition
        no: bool or float, value to fill pixels which don't meet the condition

    Outputs:
        mask
    '''

    coords = np.asarray(coords)
    s = X.shape
    posi = np.array([X,Y])  # matrix of positions [ra,dec] of the pixels
    posi = posi.transpose((1,2,0))
    posi = posi.reshape(-1,posi.shape[2])
    mask = Path(coords).contains_points(posi).reshape(s[0],s[1])

    return converter(mask,yes,no)


def sector_ring(X,Y,x0,y0,r1,r2,th1,th2,yes=True,no=False):
    '''
    To create mask array.

    Inputs:
        X, Y: 2D array, coordinates of the points, shape=(ny,nx).
        x0, y0: center of the sector ring 
        r1, r2: minor & major radius of the ring
        th1, th2: position angles of the sector, [0,360], in deg
        yes: bool or float, value to fill pixels which meet the condition
        no: bool or float, value to fill pixels which don't meet the condition
    
    Outputs:
        mask
    '''
    
    th1 *= np.pi/180
    th2 *= np.pi/180
    R = np.hypot(X-x0,Y-y0)
    Th = np.arctan2(Y-y0,X-x0)
    Th[Th<0] = Th[Th<0]+2*np.pi
    m1 = R>=r1
    m2 = R<=r2
    m3 = Th>=th1
    m4 = Th<=th2
    if th1 <= th2:
        mask = np.logical_and.reduce((m1,m2,m3,m4))
    else:
        mask = np.logical_and.reduce((m1,m2,np.logical_or(m3,m4)))

    return converter(mask,yes,no)


def elliptic_sector_ring(X,Y,x0,y0,a1,a2,b1,b2,theta,th1,th2,yes=True,
    no=False):
    '''
    To create mask array.

    Inputs:
        X, Y: 2D array, coordinates of the points, shape=(ny,nx).
        x0, y0: center of the ellipses
        a1, a2, b1, b2: major and minor axes of the ellipse
        theta: tilt angle, between a- & x-axis, in deg
        th1, th2: position angles of the sector, [0,360], in deg
        yes: bool or float, value to fill pixels which meet the condition
        no: bool or float, value to fill pixels which don't meet the condition

    Outputs:
        mask
    '''
    
    theta *= np.pi/180
    th1 *= np.pi/180
    th2 *= np.pi/180

    X1 = (X-x0)*np.cos(theta) + (Y-y0)*np.sin(theta)
    Y1 = -(X-x0)*np.sin(theta)+ (Y-y0)*np.cos(theta)
    R1 = (X1/a1)**2 + (Y1/b1)**2
    R2 = (X1/a2)**2 + (Y1/b2)**2
    Th = np.arctan2(Y-y0,X-x0)
    Th[Th<0] = Th[Th<0]+2*np.pi

    m1 = R1>=1
    m2 = R2<=1
    m3 = Th>=th1
    m4 = Th<=th2
    if th1 <= th2:
        mask = np.logical_and.reduce((m1,m2,m3,m4))
    else:
        mask = np.logical_and.reduce((m1,m2,np.logical_or(m3,m4)))

    return converter(mask,yes,no)


def test():

    import matplotlib.pyplot as plt

    x = np.linspace(0,100,100)**1
    y = np.linspace(0,100,100)**1
    X, Y = np.meshgrid(x,y)

    x0 = 50
    y0 = 40
    r = 50
    a = 60
    b = 10    
    theta = 30
    coords = [[10,10],[20,5],[130,50],[40,70],[70,20],[20,30],[10,60],[10,10]]

    r1 = 5
    r2 = 15
    th1 = 120
    th2 = 240
    a1, a2 = 25,30
    b1, b2 = 13, 27

    '''
    plt.figure()
    plt.imshow(circle(X,Y,x0,y0,r),origin='lower')
    plt.show()
    #'''

    '''
    plt.figure()
    plt.imshow(ellipse(X,Y,x0,y0,a,b,theta),origin='lower')
    plt.show()    
    #'''

    '''
    plt.figure()
    plt.imshow(polygon(X,Y,coords),origin='lower')
    plt.show()
    #'''

    '''
    plt.figure()
    plt.imshow(sector_ring(X,Y,x0,y0,r1,r2,th1,th2),origin='lower')
    plt.show()
    #'''

    '''
    plt.figure()
    plt.imshow(elliptic_sector_ring(X,Y,x0,y0,a1,a2,b1,b2,theta,th1,th2),
        origin='lower')
    plt.show()
    #'''






















