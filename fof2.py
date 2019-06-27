'''
To find velocity-coherent structures.

Inputs:
	img

Outputs:
	grp -- list, groups of points
	iso -- list, isolated points
'''

import sys
sys.path.append('/home/hp/astro/fila/modules')
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pickle

import inner
from UnitConversion import *

# load data
img = pickle.load(open('variables/img.p','rb'))
s = img.shape

I_th = .3  # intensity threshold
r0 = .2 # in pc
dv0 = 2 # km/s/pc
r = l_ra2pix(r0)-l_ra2pix(0) # convert pc to pix
dv = (v2pix(dv0)-v2pix(0))/(l_ra2pix(1)-l_ra2pix(0))  # convert km/s/pc to pix/pix
region = 1

# change the type of img to list
lop = []  # list of points [ra,dec,v]
for i in range(s[2]):
	for j in range(s[1]):
		for k in range(s[0]):
			if img[k,j,i] > I_th:
				lop.append([i,j,k])

#''' create the range of indices i,j,k in the given region
ind = []
for i in range(-int(r),int(r)+1):
	for j in range(-int(r),int(r)+1):
		for k in range(-int(r*dv),int(r*dv)+1):
			if inner.inner(0,0,0,i,j,k,r,dv,1):
				ind.append([i,j,k])
ind = np.asarray(ind)
#'''


# ================================== fof method: find out the continuums ===================================

iso = []
grp = []
L = len(lop)
while not lop==[]:
	next = []
	next.append(lop[0])
	lop.pop(0)
	p = 0     
	while p < len(next):
		for i,j,k in [min(s[2]-1,next[p][0]+1),next[p][1],next[p][2]],[max(0,next[p][0]-1),next[p][1],next[p][2]],\
			     [next[p][0],min(s[1]-1,next[p][1]+1),next[p][2]],[next[p][0],max(0,next[p][1]-1),next[p][2]],\
			     [next[p][0],next[p][1],min(s[0]-1,next[p][2]+1)],[next[p][0],next[p][1],min(0,next[p][2]-1)]:
			if [i,j,k] in lop:
				next.append([i,j,k])
				lop.pop(lop.index([i,j,k]))
		p += 1

	if len(next) <= 10:
		iso.append(next)
	else:
		grp.append(next)
	print('%.3f%%'%(100*(1-len(lop)/L)))





# =========================================== analysis ===========================================

l_grp = len(grp)
n_grp = np.zeros(l_grp,dtype=np.int)  # number of points in each group
for i in range(l_grp):
	n_grp[i] = len(grp[i])
N_grp = sum(n_grp)  # number of points in 'points'
N_iso = 0
n_iso = np.zeros(len(iso),dtype=np.int)
for i in range(len(iso)):
	n_iso[i] = len(iso[i])
	N_iso += len(iso[i])



# ============================================ output ============================================

pickle.dump(grp,open('variables/grp.p','wb'))
pickle.dump(iso,open('variables/iso.p','wb'))


# =========================================== display ============================================

# display the initial parameters
print('region = %d'%region)
print('dv = %.2f km/s/pc'%dv0)
print('r = %.3f pc'%r0)

# components analysis
print('l_grp = %d'%l_grp)
print('N_grp = %d, ratio = %.1f%%'%(N_grp,N_grp/L*100))
print('N_iso = %d, ratio = %.1f%%'%(N_iso,N_iso/L*100))

''' histgram of n_grp
plt.figure()
plt.hist(n_grp,bins=500)
plt.title('histgram of n_grp')
plt.show()
#'''

''' histgram of n_iso
plt.figure()
plt.hist(n_iso)
plt.title('histgram of n_iso')
plt.show()
#'''

#------------------------------------- 3D graphs ---------------------------------

ra = []
dec = []
v = []
for i in range(l_grp):
	ra.append([])
	dec.append([])
	v.append([])
	for j in range(n_grp[i]):
		ra[i].append(grp[i][j][0])
		dec[i].append(grp[i][j][1])
		v[i].append(grp[i][j][2])

''' 3d groups
color_set = ['b','g','r','c','m','y']
fig = plt.figure()  
ax = fig.gca(projection='3d')   
plt.axis('equal')
ax.set_xlabel('l_ra[pc]')
ax.set_ylabel('l_dec[pc]')
ax.set_zlabel('v[km/s]')
for i in range(len(grp)):
	c = color_set[divmod(i,len(color_set))[1]]
	ax.scatter(ra[i],dec[i],v[i],s=2,edgecolors='none',c=c)
plt.title('3d groups, number of groups: %d'%l_grp)
plt.show()
#'''

''' 3d isolated points
fig = plt.figure()  
ax = fig.gca(projection='3d')   
plt.axis('equal')
ax.set_xlabel('l_ra[pc]')
ax.set_ylabel('l_dec[pc]')
ax.set_zlabel('v[km/s]')
for i in range(len(iso)):
	for j in range(len(iso[i])):
		ax.scatter(iso[i][j][2],iso[i][j][1],iso[i][j][0],s=6,edgecolors='none',c='k')
plt.title('3d isolated points')
plt.show()
#'''

#------------------------------------- 2D graphs ---------------------------------

''' 2d groups
color_set = ['b','g','r','c','m','y']
plt.figure()
plt.axis('equal')
plt.xlabel('l_ra[pc]')
plt.ylabel('l_dec[pc]')
for i in range(len(grp)):
	c = color_set[divmod(i,len(color_set))[1]]
	plt.scatter(ra[i],dec[i],s=2,edgecolors='none',c=c)
plt.title('2D groups, number of groups: %d, I_th = %.2f'%(l_grp,I_th))
plt.show()
#'''

''' 2d groups pv
color_set = ['b','g','r','c','m','y']
plt.figure()
for i in range(len(grp)):
	c = color_set[divmod(i,len(color_set))[1]]
	plt.scatter(dec[i],v[i],s=2,edgecolors='none',c=c)
plt.title('2D groups pv, number of groups: %d, I_th = %.2f'%(l_grp,I_th))
plt.show()
#'''

''' 2d isolated points
plt.figure()
plt.axis('equal')
plt.xlabel('l_ra[pc]')
plt.ylabel('l_dec[pc]')
for i in range(len(iso)):
	for j in range(len(iso[i])):
		plt.scatter(iso[i][j][0],iso[i][j][1],s=6,edgecolors='none',c='k')
plt.title('2d isolated points, number of iso: %d, I_th = %.2f'%(len(iso),I_th))
plt.show()
#'''


# ====================================== recycle bin =================================


'''

N = 0
while not lop==[]:
	next = []
	next.append(lop[0])
	lop.pop(0)
	p = 0     
	while p < len(next):
		for i,j,k in [[next[p][0]+ind[m,0],next[p][1]+ind[m,1],next[p][2]+ind[m,2]] for m in range(len(ind))]:
			if [i,j,k] in lop:
				next.append([i,j,k])
				lop.pop(lop.index([i,j,k]))
			N += 1
		p += 1
	if len(next) <= 1:
		iso.append(next)
	else:
		grp.append(next)
	print('%.3f%%,  N = %d'%(100*(1-len(lop)/L),N))
'''








