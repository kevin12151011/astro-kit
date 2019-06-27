import numpy as np

def dust(img,T_ex,lamda,kappa_niu,dist,th_p,out,optical='thick'):
	'''
	To calculate the mass using the dust emission data img.

	Inputs:
		img			-- [Jy/sr], 1D or 2D data of the dust emission lines
		T_ex		-- [K], excitation temperature
		lamda		-- [mm], wavelength of the emission
		kappa_niu	-- [cm^2/g], opacity coefficient
		dist		-- [pc], distance of the source
		th_p		-- [deg], the opening angle of one pixel
		out			-- 'M' for mass, 'N' for column density output
		optical		-- thin or thick, thich by default

	Outputs:
		N_H2	-- [cm^-2], column density of H2
		M			-- [M_sun], the mass

	Note:
		CAUTION: If optically thick, you can't sum the img then calculate its N_H2, 'cause N_H2 is
				 no longer proportional to I_niu!
	'''
	
	
	if optical=='thin':
		N_H2 = 5.38e14*(np.exp(1.44/(T_ex/10)/lamda)-1)/(kappa_niu/.01)*img*lamda**3
	elif optical=='thick':
		N_H2 = -2.14e25/(kappa_niu/.01)* np.log(1-2.51e-11*img*lamda**3*(np.exp(1.44/(T_ex/10)/lamda)-1))
		
	M = 7.52e-24* dist**2*N_H2*th_p**2

	if out=='M':
		return M
	elif out=='N':
		return N_H2


def CO(trans,img,T_ex,tau,dv,dist,th_p,out):
	'''	
	To calculate the mass using the CO rotational emission data img.

	Inputs:
		trans	-- transition, '1-0', '2-1', etc.
		img		-- [Jy/sr], data of the CO rotational emission lines
		T_ex	-- [K], excitation temperature
		tau		-- optical depth
		dv 		-- [km/s], delta v in one pixel along the v-axis
		dist	-- [pc], distance of the source	
		th_p	-- [deg], the opening angle of one pixel
		out		-- 'M' for mass, 'N' for column density output

	Outputs:
		N_H2	-- [cm^-2], column density of H2, with the same shape as img
		M		-- [M_sun], the mass with the same shape as img

	Notes:
		1. The mass is proportional to the ratio n(H2)/n(CO), which is 1e4 by default.
		2. Only one of the emission lines is used in this function.
	'''
	
	
	coef = { '7-6': {'K1p': 44.55447028232269, 'K2': 154.88392948020572},
			 '5-4': {'K1p': 170.25107331205717, 'K2': 82.97560383930616},
			 '1-0': {'K1p': 105842.79235979453, 'K2': 5.53196565742015},
			 '6-5': {'K1p': 82.30745843685533, 'K2': 116.1643969144637},
			 '4-3': {'K1p': 414.8012287666839, 'K2': 55.31775085658724},
			 '2-1': {'K1p': 6619.307235136216, 'K2': 16.595720365843416},
			 '3-2': {'K1p': 1308.532469812296, 'K2': 33.191052485514604}}
	K1p, K2 = coef[trans]['K1p'], coef[trans]['K2']
	
	H2_CO = 1e4
	N_H2 = H2_CO* K1p*np.exp(K2/T_ex)*(T_ex+.92)* tau/(1-np.exp(-tau))* img*dv
	M = 7.52e-24* dist**2*N_H2*th_p**2
	
	if out=='M':
		return M
	elif out=='N':
		return N_H2



















