'''
To generate coefficients for mass calculation.

Output:
	coef -- K1p, K2, K3p
'''

import numpy as np

# constants in cgs
B = 57635.96e6
h = 6.6261e-27
k = 1.3807e-16
c = 2.9979e10

# load data
f = open('/home/kevin/astro/modules/mass/parameters.txt','r')
for i in range(3):
	f.readline()
data = f.readlines()
f.close()

# make para
para = {}
for line in data:
	line = line.split()
	para[line[0]] = {}
	para[line[0]]['niu'], para[line[0]]['A'], para[line[0]]['g_u'] = float(line[1])*1e9, 10**float(line[2]), float(line[3])
	
	
# make coef
coef = {}
for trans in para.keys():
	coef[trans] = {}
	J = float(trans.split('-')[1])
	niu = para[trans]['niu']
	g_u = para[trans]['g_u']
	A = para[trans]['A']
	
	coef[trans]['K1p'] = 1e-18* 4*np.pi*k/(B*h**2*c*g_u*A)
	coef[trans]['K2'] = h*(B*J*(J+1)+niu)/k
