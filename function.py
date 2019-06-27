from __future__ import division
import numpy as np	
from numpy import sin,cos

def Gaussian(x,para):
	'''
	Purpose: calculate the Gaussian function.

	Inputs:	
		x -- 1D array
		para -- parameters, [h,x0,sigma,...]

	Outputs:
		z -- 	z = sum(  h*exp(-.5*(((x-x0)/sigma)^2 ))  )	
	'''
	
	x = np.asarray(x,dtype=np.float64)
	z = np.zeros_like(x)
	
	for i in np.arange(0,len(para),3,dtype=int):
		h = para[i]
		x0 = para[i+1]
		sigma = para[i+2]
		z += h*np.exp(-.5*((x-x0)/sigma)**2)
	
	return z


def Gaussian_2d(x,y,para):
	'''
	Purpose: calculate the 2D Gaussian function.

	Inputs:	
		x -- 1D array
		y -- 1D array
		para -- parameters, [h,x0,y0,theta,sigma_x,sigma_y,...]

	Outputs:
		z -- 	z = sum(  h*exp(-.5*((x1/sigma_x)^2 + (y1/sigma_y)^2))  )	
			where x1 = (x-x0)*cos(theta) + (y-y0)*sin(theta)
			      y1 = -(x-x0)*sin(theta)+ (y-y0)*cos(theta)
	'''

	x = np.asarray(x,dtype=np.float64)
	y = np.asarray(y,dtype=np.float64)

	X,Y = np.meshgrid(x,y)
	z = np.zeros_like(X)
	
	for i in np.arange(0,len(para),6,dtype=int):
		h = para[i]
		x0 = para[i+1]
		y0 = para[i+2]
		theta = para[i+3]
		sigma_x = para[i+4]
		sigma_y = para[i+5]
		
		x1 = (X-x0)*cos(theta) - (Y-y0)*sin(theta)
		y1 = (X-x0)*sin(theta) + (Y-y0)*cos(theta)
		z += h*np.exp(-.5*((x1/sigma_x)**2 + (y1/sigma_y)**2))

	return z


def Gaussian_2d_r(x,y,para):
	'''
	Purpose: calculate the round 2D Gaussian function.

	Inputs:	
		x -- 1D array
		y -- 1D array
		para -- parameters, [h,x0,y0,sigma,...]

	Outputs:
		z -- z = sum( h*exp(-.5*(r/sigma)^2 ), where r^2 = (x-x0)^2+(y-y0)^2	
	'''

	x = np.asarray(x,dtype=np.float64)
	y = np.asarray(y,dtype=np.float64)

	X,Y = np.meshgrid(x,y)
	z = np.zeros_like(X)

	for i in np.arange(0,len(para),4,dtype=int):
		h = para[i]
		x0 = para[i+1]
		y0 = para[i+2]
		sigma = para[i+3]
		
		r2 = (X-x0)**2 + (Y-y0)**2
		z += h*np.exp(-.5*r2/sigma**2)
	
	return z


def obj_Gaussian(para,x,y):
	'''
	Purpose: calculate the objective function of Gaussian fitting.

	Inputs:
		para
		x -- 1D array
		y -- the data to be fitted

	Outputs:
		fitness
	'''

	fitness = np.sum( (Gaussian(x,para)-y)**2 )
	return fitness


def obj_Gaussian_2d(para,x,y,z):
	'''
	Purpose: calculate the objective function of 2D Gaussian fitting.

	Inputs:
		para
		x -- 1D array
		y -- 1D array
		z -- the data to be fitted

	Outputs:
		fitness
	'''

	fitness = np.sum( (Gaussian_2d_r(x,y,para)-z)**2 )
	return fitness


def linear(x,k,b):
    return k * np.array(x) + b













