import numpy as np
from matplotlib.path import Path


def circle(x,y,x0,y0,r,boolean='union'):
    '''
    To check whether point (x,y) is in circles.

    Inputs:
        x,y         -- scalars or ndarrays, coords of points
        x0,y0       -- scalars or 1darrays, coords of centers of circles
        r           -- scalar or 1darray, radius (radii)
        *boolean    -- optional, 'intersection' or 'union', boolean operation on circles
        
    Outputs:
        ans         -- bools of the same shape as x, True for in the region
    '''

    # data type processing
    flag_scalar = False
    if np.isscalar(x):
        flag_scalar = True
        x, y = np.array([x]), np.array([y])
    else:
        x, y = np.array(x), np.array(y)
        
    if np.isscalar(x0):
        x0, y0, r = np.array([x0]), np.array([y0]), np.array([r])
    else:
        x0, y0, r = np.array(x0), np.array(y0), np.array(r)
        
    # main procedure
    if boolean=='union':
        ans = np.zeros_like(x,dtype=np.float)
        for i in range(len(x0)):
            ans += np.sign(r[i]-((x-x0[i])**2+(y-y0[i])**2)**.5) + 1
    elif boolean=='intersection':
        ans = np.ones_like(x,dtype=np.float)
        for i in range(len(x0)):
            ans *= np.sign(r[i]-((x-x0[i])**2+(y-y0[i])**2)**.5) + 1
    
    if flag_scalar:
        ans = ans[0]
    ans = ans.astype(bool)
    
    return ans


def ellipse(x,y,x0,y0,a,b,alpha):
    '''
    To check if the point (x,y) is in the eclipse.

    Inputs:
        x,y     -- scalar, coord of the point
        x0,y0   -- scalar, coord of the center
        a, b    -- the major and minor axis of the eclipse
        alpha   -- the angle between the major axis and x-axis
        
    Outputs:
        ans -- bool, True for in the eclipse
    '''

    x1 = (x-x0)*np.cos(alpha) + (y-y0)*np.sin(alpha)
    y1 = -(x-x0)*np.sin(alpha)+ (y-y0)*np.cos(alpha)
    if (x1/a)**2+(y1/b)**2 <=1:
        ans = True
    else:
        ans = False
    return ans
    
    
def polygon(x,y,coords):
    '''
    To check if the point (x,y) is in the polygon.
    
    Inputs:
        x, y    -- scalar, coord of the point
        coords  -- N*2 array, coords of the vertices of the polygon [x,y]

    Outputs:
        ans -- bool, True for in the polygon
    '''
        
    return Path(np.array(coords)).contains_point([x,y])
    

def test_circle():
    from matplotlib.patches import Circle
    import matplotlib.pyplot as plt
    
    #x, y = [1,2,3,4,5,4,3,2], [4,7,6,3,8,4,5,1]
    #x0, y0, r = [5,3], [5,4] ,[3,2]
    x, y = 4,7
    x0, y0, r = 5, 5, 3
    ans = circle(x,y,x0,y0,r,boolean='intersection')

    # display
    plt.figure()
    ax = plt.subplot()
    plt.axis('equal')
    '''
    for i in range(len(x)):
        if ans[i]:
            c = 'r'
        else:
            c = 'k'
        plt.scatter(x[i],y[i],c=c,edgecolor='none',s=100)
    '''
    if ans:
        c='r'
    else:
        c='k'
    plt.scatter(x,y,c=c,edgecolor='none',s=100)
    #for i in range(len(x0)):
    #   ax.add_patch(Circle((x0[i],y0[i]),r[i],facecolor='none',edgecolor='k'))
    ax.add_patch(Circle((x0,y0),r,facecolor='none',edgecolor='k'))
    plt.show()




























