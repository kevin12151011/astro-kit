from __future__ import division

def dust(img,T_ex,lamda,kappa_niu,dist,th_p,th_maj,th_min):
	'''
	To calculate the mass using the dust emission data img.

	Inputs:
		img		-- [Jy/beam], 1D or 2D data of the dust emission lines
		T_ex		-- [K], excitation temperature
		lamda		-- [mm], wavelength of the emission
		kappa_niu	-- [cm^2/g], opacity coefficient
		dist		-- [pc], distance of the source
		th_p		-- [deg], the opening angle of one pixel
		th_maj		-- [deg], size of the major axis of the beam
		th_min		-- [deg], size of the minor axis of the beam

	Outputs:
		M		-- [M_sun], the mass

	Note:
		1. All assuming that the gas is optically thin
	'''
	
	import numpy as np

	A = (np.exp(1.44/(T_ex/10)/lamda)-1)*lamda**3
	M = .118*A/(kappa_niu/.01)*img*(dist/100)**2*th_p**2/th_maj/th_min

	return M



def CO(img,J,T_ex,dist,th_p,th_maj,th_min,tau,dv,H2_CO_ratio=1e4):
	'''	
	To calculate the mass using the CO rotational emission data img.

	Inputs:
		img		-- [Jy/beam], data of the CO rotational emission lines
		J		-- '1-0','2-1',...
		T_ex		-- [K], excitation temperature
		dist		-- [pc], distance of the source	
		th_p		-- [deg], the opening angle of one pixel
		th_maj		-- [deg], size of the major axis of the beam
		th_min		-- [deg], size of the minor axis of the beam		
		tau		-- optical depth	
		dv 		-- [km/s], delta v in one pixel along the v-axis

	Outputs:
		M		-- [M_sun], the mass

	Notes:
		1. The mass is proportional to the ratio n(H2)/n(CO), which is 1e4 by default.
		2. Only one of the emission lines is used in this function.
	'''

	import numpy as np

	_K2 = {	'1-0':	5.53,\
			'2-1':	16.6	}
	_K4 = {	'1-0':	2.15e-5,\
			'2-1':	1.34e-6}
	K2 = _K2[J]
	K4 = _K4[J]

	A = 4*np.log(2)/np.pi
	M = A*K4*np.exp(K2/T_ex)*(T_ex+.92)*(dist/1000)**2*th_p**2/th_maj/th_min*tau/(1-np.exp(-tau))*img*dv
	M *= H2_CO_ratio/1e4

	return M

	



















