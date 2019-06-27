from __future__ import division
import numpy as np
from scipy.optimize import fsolve

''' About the module:
Purpose: To derive physical properties using ammonia emission lines.
'''


def Temperature(I22,*I11):
	'''
	To calculate the kinematic temperature of gas using NH3 (1,1) & (2,2) lines.

	Inputs:
		I22		-- peak intensity of the NH3 (2,2) line
		I11		-- [in I22's unit], [[[Ip1], Ip2], Ip3, [Ip4, [Ip5]]]

	Outputs:
		tau	-- optical depth, = tau if I11 has no satellite data, ( =np.nan if not well-defined, the same below)
		T_r		-- rotational T
		T_k		-- kinematic T

	Notes:
		1. The same beam-filling factor is assumed.
		2. The results are robust when Tk = 1-40K. 
		3. For more info, see the doc in astro/formula/ammonia2T.pdf.
	'''

	# constants
	a = [.28, .22] # theoretical peak intensity ratio, 11m/11s

	I11_m = I11[int(len(I11)/2)]
	# tau
	tau = np.nan
	if len(I11)==3:
		ratio = I11[1]/((I11[0]+I11[2])/2)
		if 0 < ratio < 1/a[0]:
			try:
				tau = fsolve(lambda x : func_tau(x,ratio,a[0]),1e-9)
			except:
				pass

	elif len(I11)==5:
		Ratio = [I11[2]/((I11[1]+I11[3])/2), I11[2]/((I11[0]+I11[4])/2)]
		Tau = [np.nan]*2
		for i in range(len(Ratio)):
			if 0 < Ratio[i] < 1/a[i]:
				try:
					Tau[i] = fsolve(lambda x : func_tau(x,Ratio[i],a[i]),1e-9)
				except:
					pass

		if not Tau==[np.nan,np.nan]:
			tau = np.nanmean(Tau)
	else:
		print('Func Temperature: Input Error.')

	# T_r & T_k
	T_r = -41.5/(np.log(-.282/tau*np.log(1-I22/I11_m*(1-np.exp(-tau)))))
	T_k = T_r/(1-T_r/41.5*np.log(1+1.608*np.exp(-25.25/T_r)))

	return tau, T_r, T_k


def func_tau(x,ratio,a):
	return (1-x)/(1-x**a)-ratio





































