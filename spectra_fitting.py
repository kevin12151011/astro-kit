''' About the module:
Purpose: To find peaks & fit with Gaussians in the spectrum.

Functions:
        line_fit
        find_peaks
        fp_NH3
        smooth
        gaussian_fitter
        ft_NH3
        func_NH3
        gaussian
'''

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit



def line_fit(x,y,**lim): # ----------------------------------------------------
    '''
    Purpose: find and fit the spectral lines.
    
    Inputs:
        x   -- the independent variable
        y   -- the spectrum to be fitted
        **lim  -- upper and lower limits of a,b,c

    Outputs:
        flag: True for successful fitting
        para: [a,b,c,a,b,c,...], [] for unsuccessful fitting
    '''

    # find peaks & fitting
    para0 = find_peaks(x,y,sign=1,**lim)
    flag, para = gaussian_fitter(x,y,para0,**lim)
    
    ''' demo
    print('Results from find_peaks:')
    print('id    amp    pos    wid')
    for i in range(int(len(para0)/3)):
        print('%d\t%.2f\t%.2f\t%.2f\t%.2f'%(i,para0[3*i],para0[3*i+1],para0[3*i+2]))

    print('Results from line_fit:')
    print 'id    amp    pos    wid'
    for i in range(int(len(para)/3)):
        print('%d\t%.2f\t%.2f\t%.2f\t%.2f'%(i,para[3*i],para[3*i+1],para[3*i+2]))
    #'''

    return flag, para


def find_peaks(x,y,sign=1,**lim): # --------------------------------------------------
    '''
    Purpose: Find peaks in dataset (x,y).

    Inputs:
        x
        y       -- x & y are 1D array
        sign    -- sgin of peaks to be found, 1 for positive, -1 for negative, 0 for both
        **lim   -- 'a_min', 'a_max', 'b_min', 'b_max', 'c_min', 'c_max', 

    Outputs: 
        para0   -- (a,b,c,a,b,c,...)

    Algorithm:
        1. Calculate dy/dx and smooth it according to th_wid. According to the test, 
           peaks with widths < length/5 will be smooth out efficiently.
        2. Find peaks by identify zero-across in dy/dx.
        3. Check the amplitudes of the peaks.
    '''

    para0 = []

    # constants
    th_snr_dy = 1
    # th_snr_dy = 0.043
    smth_ratio = 1.4
        # time of smooth, in theory, the equivalent smooth length should be 
        # 2*sigma, so smth_ratio * N_sm**.5 = 2.
    N_sm = 2 

    for i in range(N_sm):
        if 'c_min' in lim:
            y = smooth(x,y,smth_ratio*lim['c_min'])

    dy = np.append(np.diff(y)/np.diff(x),0)
    if 'c_min' in lim:
        dy = smooth(x,dy,smth_ratio*lim['c_min'])
    std_dy = np.std(dy)

    mono = np.zeros_like(dy,dtype=np.int8)
    mono[dy > th_snr_dy*std_dy] = 1
    mono[dy < -th_snr_dy*std_dy] = -1
    
    # find peaks & check
    if sign==1 or sign==0:
        # find positive peaks
        i = 0
        while i < len(x)-1:
            find = False
            if mono[i]==1:
                j = i 
                while j < len(x)-1:
                    if mono[j+1] > mono[j]:
                        break
                    j += 1
                if mono[j]==-1:
                    amp = y[int((i+j)/2)]
                    pos = x[int((i+j)/2)]
                    wid = (x[j]-x[i])/2.*.5
                    # print(amp,pos,wid,x[i],x[j])
                    if lim['a_min'] < amp < lim['a_max'] and lim['b_min'] < pos < lim['b_max'] \
                        and lim['c_min'] < wid < lim['c_max']:
                        para0.extend([amp,pos,wid])
                        find = True
                i = j       
            if not find:
                i += 1

    if sign==-1 or sign==0:
        # find negative peaks 
        i = 0
        while i < len(x)-1:
            find = False
            if mono[i]==-1:
                j = i 
                while j < len(x)-1:
                    if mono[j] < mono[j-1]:
                        break
                    j += 1
                if mono[j-1]==1:
                    amp = y[int((i+j)/2)]
                    pos = x[int((i+j)/2)]
                    wid = (x[j]-x[i])/2.*.5
                    if -lim['a_max'] < amp < -lim['a_min'] and lim['b_min'] < pos < lim['b_max'] \
                        and lim['c_min'] < wid < lim['c_max']:
                        para0.extend([amp,pos,wid])
                        find = True
                i = j
            if not find:
                i += 1

    ''' display
    plt.figure()
    plt.plot(x,y)
    plt.figure()
    plt.plot(x,dy,x,dy.mean()+mono*(dy.max()-dy.mean()))
    plt.show()
    print(para0)
    #'''

    return para0


def fp_NH3(x,y,satellite,df,Spec,debug=False,**lim): # ------------------------------------------------------------
    '''
    Purpose: To find peaks in ammonia (n,n) lines.

    Inputs:
        x           -- ndarray, the independent variable
        y           -- ndarray, the spectrum to be fitted
        satellite   -- bool, True for harboring satellite lines
        Spec        -- dict for (1,1), containing the theoretical spectral info of NH3
                        sat : {'I', 'df'}, sat = 1,2,3,4,5
        debug       -- bool, whether to display the results
        **lim       -- 'a_min', 'a_max', 'b_min', 'b_max', 'c_min', 'c_max', c_min is mandatory

    Outputs:
        para0       -- (a,b,c) for one Gaussian or (Am,As1,As2,x0,sigma_x) for lines with satellites

    Notes:
        1. Only for finding one main peak in the spectrum (with satellite lines)!
        2. Algorithm:
            a. Smooth y to ys, using lim['c_min'];
            b. Find the Top5, regard them as the ammonia (1,1) lines;
            c. 3 kind of checking: 
                i. amplitude checking: for Top5
                ii. position checking: for the main line
                iii. relative position checking: for the satellites
        3. Three possible results: main, main+satellite1, main+satellite1+satellite2
    '''

    # some parameters
    smth_ratio = 8 # smooth length/lim['c_min']
    snr_m, snr_s = 1.5, 1.
    tol = 15 # [x unit], position tolerance for finding satellites

    dx = x[1]-x[0]
    dist_sate = [(Spec[3]['df'][0]-Spec[4]['df'][0])/df, (Spec[3]['df'][0]-Spec[5]['df'][0])/df]
    signal_wid = 1.1*dist_sate[1]  # for regions to calculate rms
    peak_wid = lim['c_max']
    I_weight = [np.sum(Spec[1]['I']),np.sum(Spec[2]['I']),np.sum(Spec[3]['I'])]

    rms = np.nan
    para0 = []

    # smooth
    x, y = np.array(x), np.array(y)
    ys = y.copy()
    if 'c_min' in lim:
        ys = smooth(x,ys,smth_ratio*lim['c_min'])
    ys1 = ys.copy() # only for finding maximum 

    # finding
    if satellite:

        # find the top 5
        I_agm = []
        for i in range(5):
            i_tmp = np.argmax(ys1)
            I_agm.append(i_tmp)
            ys1[int(i_tmp-peak_wid/dx):int(i_tmp+peak_wid/dx)] = -np.inf
        I_agm = np.sort(I_agm)
        Y_agm = x[I_agm]
        Y_max = ys[I_agm]

        # position checking for the main component
        if Y_agm[2]-signal_wid > 0 and Y_agm[2]+signal_wid < len(x):

            # calculate rms
            rms = np.std(np.append(y[x>Y_agm[2]+signal_wid],y[x<Y_agm[2]-signal_wid]))          

            # amplitude checking for the main component
            if Y_max[2] > snr_m*rms:
                para0.extend([Y_max[2]/I_weight[2],Y_agm[2],2*lim['c_min']])

                # relative position checking for satellite 1
                if Y_agm[3]+tol-Y_agm[2] >= dist_sate[0] >= Y_agm[3]-tol-Y_agm[2] \
                        and Y_agm[2]+tol-Y_agm[1] >= dist_sate[0] >= Y_agm[2]-tol-Y_agm[1]:

                    # amplitude checking for satellite 1    
                    if Y_max[1] > snr_s*rms and Y_max[3] > snr_s*rms:
                    
                        para0.insert(1,(Y_max[1]+Y_max[3])/2/I_weight[1])

                        # relative position checking for satellite 2
                        if Y_agm[4]+tol-Y_agm[2] >= dist_sate[1] >= Y_agm[4]-tol-Y_agm[2] \
                                and Y_agm[2]+tol-Y_agm[0] >= dist_sate[1] >= Y_agm[2]-tol-Y_agm[0]:

                            # amplitude checking for satellite 2
                            if Y_max[0] > snr_s*rms and Y_max[4] > snr_s*rms:
                                para0.insert(2,(Y_max[0]+Y_max[4])/2/I_weight[0])
    else:
        # find max
        i_agm = np.argmax(ys)
        y_agm = x[i_agm]
        y_max = ys[i_agm]

        # position checking
        if y_agm-signal_wid > 0 and y_agm+signal_wid < len(x):

            # calculate rms
            rms = np.std(np.append(y[x>y_agm+signal_wid],y[x<y_agm-signal_wid]))            

            # amplitude checking for the main component
            if y_max > snr_m*rms:
                para0.extend([y_max,y_agm,2*lim['c_min']])

    # debugging
    if debug:
        plt.figure()
        plt.plot(x,y,'k',alpha=.5,linewidth=1)
        plt.plot(x,ys,'r',linewidth=1)  
        plt.plot(x,ys1,'g',linewidth=4,alpha=.4)
        plt.plot(x,np.zeros_like(x)-snr_m*rms,'k--')
        plt.plot(x,np.zeros_like(x)+snr_m*rms,'k--')
        plt.plot(x,np.zeros_like(x)-snr_s*rms,'k:')
        plt.plot(x,np.zeros_like(x)+snr_s*rms,'k:')
        if para0:
            yf = func_NH3(x,df,Spec,*para0)
            plt.plot(x,yf,'b',linewidth=1)
        plt.show()
    
    return para0


def smooth(x,y,l): # ------------------------------------------------------
    '''
    Purpose: Smooth y in the dataset (x,y).

    Inputs:
        x -- must be equally spaced
        y -- x & y must be of the same shape
        l -- length of interval for smooth

    Output:
        ys
    '''

    r_pix = int(round(l/(x[1]-x[0])/2))
    l_pix = 2*r_pix
    ys = {}
    if l_pix==0:
        ys = y.copy()
    else:
        cdf = y.copy()
        cdf = np.append(np.zeros(r_pix),cdf)
        cdf = np.cumsum(np.append(cdf,np.zeros(r_pix)))
        ys = (cdf[l_pix:] - cdf[:-l_pix])/l_pix

    return ys


def gaussian_fitter(x,y,para0,**lim): # ------------------------------------------
    '''
    Purpose: multi-Gaussian fitting

    Inputs:
        x       -- independent variable
        y       -- the data
        para0   -- initial parameters, order:[a1,b1,c1,a2,b2,c2,...]
        **lim   

    Outputs:
        para: [a,b,c,a,b,c,...], [] for unsuccessful fitting
        flag: True for successful fitting
    ''' 

    n = int(len(para0)/3)

    # limits
    a_min = lim['a_min'] if 'a_min' in lim else -np.inf
    b_min = lim['b_min'] if 'b_min' in lim else -np.inf
    c_min = lim['c_min'] if 'c_min' in lim else -np.inf
    a_max = lim['a_max'] if 'a_max' in lim else np.inf
    b_max = lim['b_max'] if 'b_max' in lim else np.inf
    c_max = lim['c_max'] if 'c_max' in lim else np.inf
    bound=([a_min,b_min,c_min]*n, [a_max,b_max,c_max]*n)

    flag = True
    para = []
    para_orig = []
    try:
        para_orig, pcov = curve_fit(gaussian,x,y,para0,bounds=bound,
                            method='trf') # method = 'lm', 'trf', 'dogbox'
    except:
        flag = False

    else:
        for a,b,c in para_orig.reshape((-1,3)):
            if a_min<a<a_max and b_min<b<b_max and c_min<c<c_max:
                para.extend([a,b,c])
        para = np.array(para)
        
    finally:
        return flag, para


def ft_NH3(x,y,df,Spec,*para0,**lim): # ---------------------------------------------------------------
    '''
    Purpose: fitting with NH3 (1,1) model.

    Inputs:
        x       -- independent variable
        y       -- the data
        df      -- df of one pix
        Spec    -- dict, containing the theoretical spectral info of NH3
                        sat : {'I', 'df'}, sat = 1,2,3,4,5
        *para0  -- initial parameters, order:[Am,[As1,[As2]],x0,sigma_x]
        **lim

    Outputs:
        para     -- parameters after fitting
        flag_bad
    ''' 

    a_min = lim['a_min'] if 'a_min' in lim else -np.inf
    a_max = lim['a_max'] if 'a_max' in lim else np.inf
    b_min = lim['b_min'] if 'b_min' in lim else -np.inf
    b_max = lim['b_max'] if 'b_max' in lim else np.inf
    c_min = lim['c_min'] if 'c_min' in lim else -np.inf
    c_max = lim['c_max'] if 'c_max' in lim else np.inf
    if len(para0)==3:
        bounds = ([a_min,b_min,c_min], [a_max,b_max,c_max])
    elif len(para0)==4:
        bounds = ([a_min,a_min,b_min,c_min], [a_max,a_max,b_max,c_max])
    elif len(para0)==5:
        bounds = ([a_min,a_min,a_min,b_min,c_min], [a_max,a_max,a_max,b_max,c_max])

    try:
        para, pcov = curve_fit(lambda x,*para0: func_NH3(x,df,Spec,*para0),x,y,para0,bounds=bounds,method='trf') 
        # method = 'lm', 'trf', 'dogbox'
        flag_bad = 0
        return flag_bad, para

    except:
        para = np.zeros_like(para0)
        flag_bad = 1
        return flag_bad, para   
    

def func_NH3(x,df,Spec,*para): # -------------------------------------------
    '''
    To calculate the NH3 (1,1) line function.

    Inputs:
        x           -- the independent variable
        df          -- df of one pix
        Spec        -- dict, containing the theoretical spectral info of NH3
                        sat : {'I', 'df'}, sat = 1,2,3,4,5
        *para       -- Am,[As1,[As2]],x0,sig_x, depending on how many components are detected
                        A's are relative intensities for

    Outputs:
        y           -- the spectral function

    Notes:
        1. Only one velocity component is considered here.
    '''

    sat_ind = { 1:  2,\
                2:  1,\
                3:  0,\
                4:  1,\
                5:  2 }
    Sat = []
    if len(para)==3:
        Sat = [3]
    elif len(para)==4:
        Sat = [2,3,4]
    elif len(para)==5:
        Sat = [1,2,3,4,5]
    else:
        print('func_NH3: input error')

    y = np.zeros_like(x)
    para_sp = []
    for sat in Sat:
        ind = sat_ind[sat]
        for i in range(len(Spec[sat]['I'])):
            para_sp.extend([Spec[sat]['I'][i]*para[ind], para[-2]-Spec[sat]['df'][i]/df, para[-1]])
    y = gaussian(x,*para_sp)

    return y
    
    
def gaussian(x,*p): # ------------------------------------------------------------
    # p: [a1,b1,c1,a2,...]
    y = np.zeros_like(x)*1.
    for a,b,c in np.reshape(p,(-1,3)):
        y += a*np.exp(-.5*((x-b)/c)**2)

    return y


def test(): # ------------------------------------------------------------
    n = 1000
    x = np.linspace(-1,1,n)
    A = 2
    An = .6
    sigma = .05
    y0 = A*np.exp(-.5*(x/sigma)**2) + .7*A*np.exp(-.5*((x-4*sigma)/sigma)**2) 
    y = y0 + An*np.random.normal(size=len(y0))
    lim = { 'a_max':    2,\
            'a_min':    .1,\
            'b_max':    .5,\
            'b_min':    -.5,\
            'c_max':    .2,\
            'c_min':    .05 }

    para0 = find_peaks(x,y,sign=1,**lim)
    print(para0)

    fy = gaussian(x,*para0)


    plt.figure()
    plt.plot(x,y,x,fy)
    plt.show()












