from __future__ import division

def coord_rot(x,y,z,axis,theta,rotate_coord):
	'''
	Purpose: rotation transformation of Cartesian coordinate system.

	Inputs:
		x,y,z -- scalar or ndarray, must be in the same size,original coordinates
		axis -- 'x','y' or 'z', the axis along which to rotate
		theta -- rotation angle, in rad
		rotate_coord -- 0 for rotating the points, 1 for rotating the coord system

	Outputs:
		x1,y1,z1 -- coordinates after rotation
	'''

	import numpy as np
	from numpy import sin,cos

	x = np.asarray(x)
	y = np.asarray(y)
	z = np.asarray(z)

	if rotate_coord==1:
		theta = -theta

	if axis=='x':
		ind = 0
	if axis=='y':
		ind = 1
	if axis=='z':
		ind = 2

	T = np.zeros([3,3])  # rotation matrix
	T[ind,ind] = 1.
	T[divmod(ind+1,3)[1],divmod(ind+1,3)[1]] = T[divmod(ind+2,3)[1],divmod(ind+2,3)[1]] = cos(theta)
	T[divmod(ind+1,3)[1],divmod(ind+2,3)[1]] = sin(theta)
	T[divmod(ind+2,3)[1],divmod(ind+1,3)[1]] = -sin(theta)

	x1 = np.zeros_like(x)
	y1 = np.zeros_like(x)
	z1 = np.zeros_like(x)

	s = x.shape
	x = x.flatten()
	y = y.flatten()
	z = z.flatten()

	coord1 = np.dot(T,[x,y,z])
	x1,y1,z1 = coord1[0],coord1[1],coord1[2]
	x1.reshape(s)
	y1.reshape(s)
	z1.reshape(s)

	return x1,y1,z1

