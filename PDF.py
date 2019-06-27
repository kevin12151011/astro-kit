import numpy as np 
import matplotlib.pyplot as plt 


def pdf2points(x,pdf,N):
    '''
    To generate random points given pdf.

    Inputs:
        x: 1darray
        pdf: 1darray, must > 0 and without nans.
        N: scalar, number of points.

    Outputs:
        x_sample: scalar or 1darray.
    '''

    # create CDF
    cdf = np.cumsum(pdf*np.gradient(x)) 
    cdf = (cdf-cdf[0])/(cdf[-1]-cdf[0])

    # generate random points
    x_sample = np.interp(np.random.rand(N),cdf,x)
    if N==1:
        x_sample = x_sample[0]

    return x_sample



def test():
    x = np.linspace(-2,2,1000)
    pdf = np.exp(-x**2/(2*.7**2))+5*np.exp(-(x-1)**2/(2*.1**2))
    N = 5000
    x_sample = pdf2points(x,pdf,N)

    plt.figure()
    plt.plot(x,5*pdf)
    plt.hist(x_sample,bins=int(N/5))
    plt.show()



# The function below has bugs, revise it before using!

# def pdf2points_2d(pdf,N):
#     '''
#     To generate random points given arbitrary pdf (ndarray).

#     Inputs:
#         pdf: 2darray, must > 0 and without nans.
#         N: scalar, number of points.

#     Outputs:
#         Coord: ndarray with shape (N,2).
#     '''

#     pdf = np.array(pdf)
#     Coord = np.full((N,2),np.nan)
#     s = pdf.shape

#     # create CDF
#     cdf_y = np.cumsum(np.sum(pdf,axis=1)) # marginal CDF of y
#     cdf_y = (cdf_y-cdf_y[0])/(cdf_y[-1]-cdf_y[0])
#     cdf_y_x = np.cumsum(pdf,axis=1) # conditional CDF of x given y
#     for i in range(s[0]):
#         cdf_y_x[i] = (cdf_y_x[i]-cdf_y_x[i,0])/(cdf_y_x[i,-1]-cdf_y_x[i,0])

#     # generate random points
#     Y = np.interp(np.random.rand(N),cdf_y,np.arange(s[0]))
#     X = np.full_like(Y,np.nan)
#     for i in range(N):
#         y = int(round(Y[i]))
#         X[i] = np.interp(np.random.rand(),cdf_y_x[y],np.arange(s[1]))
#     Coord = np.transpose(np.array([X,Y]))

#     return Coord


