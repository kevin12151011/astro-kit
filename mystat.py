import numpy as np 



def linear_regression(x,y):
    '''
    1D Linear regression analysis.

    Inputs:
        x, y: 1darray

    Outputs:
        (k,b): slop & intercept
    '''

    x = np.array(x)
    y = np.array(y)

    # drop nans
    ind = np.logical_and(~np.isnan(x),~np.isnan(y))
    x = x[ind]
    y = y[ind]

    # calculation
    x_a = np.mean(x)
    y_a = np.mean(y)
    k = np.sum((y-y_a)*(x-x_a))/np.sum((x-x_a)**2)
    b = y_a - k*x_a

    return k, b


