from __future__ import division

def box(x,lw,up):
	'''
	Box function, equals 1 in the interval (lw,up), .5 when x=lw or up, 0 for the others.

	Inputs:
		x	-- the independent variable
		lw,up 	-- 'none' for no boundaries
	
	Outputs:
		y
	'''

	import numpy as np
	from numpy import sign
	
	y = np.zeros_like(x)

	if not lw=='none':
		y += (sign(x-lw)+1)/2
	if not up=='none':
		y += (sign(up-x)+1)/2

	return y

