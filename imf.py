import numpy as np
import matplotlib.pyplot as plt 
from numpy import sign, exp, log, log10, pi 
from scipy.special import erf


def imf(m,name,dlgM=False):
    '''
    To calculate IMF.

    Inputs:
        m: [Msun], scalar or 1darray
        name    -- type of IMF, can be 
                'Salpeter55'
                'Miller-Scalo'
                'Chabrier-i': for individual stars
                'Chabrier-s': for stellar systems
                'Kroupa'

        dlgM: bool, 'dM' or 'dlgM'

    Outputs:
        y: scalar or 1darray, the IMF

    Notes:
        1. The IMFs are not normed except the Kroupa IMF!
        2. If p(x) is normalized as int p(x)dx = 1, then int q(lgx)dlgx = 1.
    '''

    y = {}
    if name=='Salpeter55':
        y = m**-2.35

    elif name=='Miller-Scalo':
        y = (sign(1-m)+1)/2/m + (sign(m-1)+1)/2*m**-2.35

    elif name=='Chabrier-i':
        y = (1/.254*((sign(1-m)+1)/2* .158/log(10)/m*exp(-(log10(m)-
            (log10(.08))**2)/(2*.69**2))) + (sign(m-1)+1)/2* m**-2.35)

    elif name=='Chabrier-s':
        y = (1/.08*((sign(1-m)+1)/2* .086/log(10)/m*exp(-(log10(m)-
            (log10(.22))**2)/(2*.57**2))) + (sign(m-1)+1)/2* m**-2.35)

    elif name=='Kroupa':
        p1 = .3
        p2 = 1.3
        p3 = 2.3
        m1 = .08
        m2 = .5

        A = ((m1**(1-p1)/(1-p1) + (m1**(p2-p1)*m2**(1-p2)-m1**(1-p1))/(1-p2) -
             m1**(p2-p1)*m2**(1-p2)/(1-p3))**-1)
        B = A*m1**(p2-p1)
        C = B*m2**(p3-p2)
        y = ((sign(.08-m)+1)/2*A*m**-p1 + 
            .25*(sign(m-.08)+1)*(sign(.5-m)+1)* B*m**-p2 + 
            (sign(m-.5)+1)/2* C*m**-p3)


    elif name=='Chabrier05':
        # parameters
        mc = 1 # critical mass value
        alpha = 2.35 # slope
        m0 = .2 # peak value of lognormal
        sigma = (log10(m0)/log(10)/(1-alpha))**.5 # dispersion of lognormal

        A = ((pi/2)**.5*log(10)*sigma*(1-erf(log10(m0)/2**.5/sigma)) + 
            1/(alpha-1)*exp(-log10(m0)**2/2/sigma**2))**-1
        B = A*exp(-log10(m0)**2/2/sigma**2)

        y = ((sign(1-m)+1)/2* A/m*exp(-log10(m/m0)**2/2/sigma**2) + 
            (sign(m-1)+1)/2* B*m**-alpha)


    else:
        print('No IMF named as "%s".'%name)

    if dlgM:
        y *= (m*np.log(10)) 

    return y



def test():
    m = np.logspace(-3,3,1000)
    Name = ['Salpeter55','Miller-Scalo','Chabrier-i','Chabrier-s','Kroupa',
        'Chabrier05']

    plt.figure()
    plt.xscale('log')
    plt.yscale('log')
    for name in Name:
        plt.plot(m,imf(m,name))
    plt.legend(Name)
    plt.xlabel('m (Msun)')
    plt.ylabel('dp/dm')
    plt.show()



















