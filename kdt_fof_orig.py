'''
To find velocity-coherent structures with friends-of-friend method.

Outputs:
    grp -- list, groups of points
    iso -- list, isolated points
'''


from __future__ import division
import pickle
import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from mpl_toolkits.mplot3d import Axes3D
from scipy import spatial
import timeit



def inner(x0,y0,v0,x,y,v,r,dv):
    '''
    Purpose: Judge whether the given point is in the given velocity region
    
    Inputs:
        (x0,y0,v0) -- centeral point of the area
        (x,y,v)    -- the point to be judged
        r          -- radius of the sphere
        dv         -- gradient z

    Returns:
        flag -- True or false
    '''

    flag = False
    r1 = np.sqrt((x-x0)**2+(y-y0)**2)
    dv1 = np.abs(v-v0)/r1
    if dv1<=dv:
        flag = True

    return flag


'''

This is the best combination for identify the filamentary structure in our case: 
r0 = 0.045  &  II>0.04

'''


# control panel \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

dist = 1300 # [pc]
r0 = 0.045 # [pc], separation threshold, ~5.5 pixel
#dv0 = 3 # [km/s/pc], velocity gradient threshold
min_num_group = 20 # minimal number of points in a group
name_f = './fits/H13COp_vlsr_Int_multiv_rebin3_Nov03.fits' 

# /////////////////////////////////////////////////////////////////////////////

# load data
hdu = fits.open(name_f)
h = hdu[0].header
img = hdu[0].data
s = img.shape
th_p = h['CDELT2']



r = r0 /dist*180/np.pi/th_p  # [pix]
print('r = %.4f pc = %.4f pix'%(r0,r))
#dv = dv0 *dist/180*np.pi*th_p # [km/s/pix]
#print('dv = %.4f km/s/pc = %.4f km/s/pix'%(dv0,dv))




#====================================================================
#------------ First Step ----------------
# Find the seed points that have high intenisties
#-----------------------------------------

# change the type of img to list
vv = img[0:4,:,:]       # Vlsr
II = img[4:8,:,:]       # velocity integrated intensity
sigma = img[8:12,:,:]   # line width

nimg = np.where(II>0.04, vv, np.nan)
nI = np.where(II>0.04, II, np.nan)

# change the type of img to list
tmap = []
lop = []  # list of points [ra,dec,v]
rdlop = []
#Ilop = []
for i in range(s[1]):
    for j in range(s[2]):
        for v in nimg[:,i,j]:
            if not np.isnan(v):
                lop.append([i,j,v])
                rdlop.append([i,j])
                tmap.append([i,j,v])
                #Ilop.append([v])

Ilop = []
for i in range(s[1]):
    for j in range(s[2]):
        for I in nI[:,i,j]:
            if not np.isnan(I):
                Ilop.append([i,j,I])


sig = []  # list of points [ra,dec,sigma]
for i in range(s[1]):
    for j in range(s[2]):
        for si in sigma[:,i,j]:
            if not np.isnan(si):
                sig.append([si])


#nas = np.where(nI==np.nanmax(nI))
#startpoint = Ilop.index([nas[1],nas[2],np.nanmax(nI)])


starttime = timeit.default_timer()
iso = []
grp = []
L = len(lop)
start = 0
while lop:
    nxt = []
    
    #if start == 0:
     #   nxt.append(lop.pop(startpoint))
     #   rdlop.pop(startpoint)
     #   start = 1
        
    nxt.append(lop.pop(0))
    rdlop.pop(0)
    p = 0  
    #while p < len(nxt):
    while p < len(nxt):   
        tree = spatial.KDTree(rdlop)
        ind = tree.query_ball_point([nxt[p][0], nxt[p][1]], r)
        tlop = []
        for i in range(len(lop)):
            tlop.append([lop[i][0], lop[i][1], lop[i][2]])
        #print 'len(tlop)',len(tlop)
        if len(tlop) <= 1:
            break
        tind = tmap.index([ nxt[p][0], nxt[p][1], nxt[p][2] ])
        dv = sig[tind][0]
        #print tlop[n][0],tlop[n][1], tlop[n][2]
        #print 'tind', tind
        print 'dv', dv
        if len(ind)>0:
            for n in ind:
                if inner(nxt[p][0],nxt[p][1],nxt[p][2],tlop[n][0],tlop[n][1],tlop[n][2],r,dv):
                    tempind = lop.index([ tlop[n][0], tlop[n][1], tlop[n][2] ])
                    nxt.append(lop.pop(tempind))
                    rdlop.pop(tempind)
                    print 'len(lop)s1',len(lop)
                    #print 'len(tlop)1',len(tlop)
                    #print 'n', n
        p += 1
        #print 'p', p
        if len(lop) <= 1:
                break
    if len(nxt) < min_num_group:
        iso.append(nxt)
    else:
        grp.append(nxt)

    print('%.3f%%'%(100*(1-len(lop)/L)))

    if len(lop) <= 1:
        break

#stoptime = timeit.default_timer()
#print('Time: ', stoptime - starttime)

#print('Found %d groups.'%len(grp))
#print('Found %d isolated.'%len(iso))




#====================================================================
#------------ Second Step ----------------
# Find the seed points friends
#-----------------------------------------
#r0 = 0.12 # [pc], separation threshold, ~44 pixel
#dv0 = 6 # [km/s/pc], velocity gradient threshold
#min_num_group = 4 # minimal number of points in a group
#r = r0 /dist*180/np.pi/th_p  # [pix]
#dv = dv0 *dist/180*np.pi*th_p # [km/s/pix]
#print('dv = %.4f km/s/pc = %.4f km/s/pix'%(dv0,dv))
#print('r = %.4f pc = %.4f pix'%(r0,r))


# change the type of img to list
lop1 = []  # list of points [ra,dec,v]
rdlop1 = []
#Ilop1 = []
for i in range(s[1]):
    for j in range(s[2]):
        for v in img[0:4,i,j]:
            if not np.isnan(v):
                lop1.append([i,j,v])
                rdlop1.append([i,j])
                #Ilop1.append([v])

L = len(lop1)

for i in range(len(grp)):
    for j in range(len(grp[i])):
        tempind = lop1.index([grp[i][j][0],grp[i][j][1],grp[i][j][2]])
        lop1.pop(tempind)
        rdlop1.pop(tempind)



#start = timeit.default_timer()
grp1 = grp
for i in range(len(grp)):
    lop1.pop(0)
    rdlop1.pop(0)
    for j in range(len(grp[i])):
        tree = spatial.KDTree(rdlop1)
        ind = tree.query_ball_point([grp[i][j][0],grp[i][j][1]], r)
        tlop = []
        for k in range(len(lop1)):
            tlop.append([lop1[k][0], lop1[k][1], lop1[k][2]])
        #print 'len(tlop)',len(tlop)
        #if len(tlop) <= 1:
            #break
        tind = tmap.index([ grp[i][j][0], grp[i][j][1], grp[i][j][2] ])
        dv = sig[tind][0]
        if len(ind)>0:
            for n in ind:
                if inner(grp[i][j][0],grp[i][j][1],grp[i][j][2],tlop[n][0],tlop[n][1],tlop[n][2],r,dv):
                    tempind = lop1.index([tlop[n][0],tlop[n][1],tlop[n][2]])
                    grp1[i].append(lop1.pop(tempind))
                    rdlop1.pop(tempind)
                    print 'len(lop1)s2',len(lop1)
                    #print 'len(tlop)s2',len(tlop)
                    #print 'n', n
        

    print('%.3f%%'%(100*(1-len(lop1)/L)))

    #if len(lop1) <= 1:
    #    break

#stop = timeit.default_timer()
#print('Time: ', stop - start)




#====================================================================
#------------ Third Step ----------------
#-----------------------------------------
# Re-make friend for isolated points 
#===========================================


min_num_group = 30
iso2 = []
iso2 = lop1

rdlop2 = []
lop2 = []
for k in range(len(iso2)):
    lop2.append([iso2[k][0], iso2[k][1], iso2[k][2]])
    rdlop2.append([iso2[k][0], iso2[k][1]])


tmap = []
for i in range(s[1]):
    for j in range(s[2]):
        for v in img[0:4,i,j]:
            if not np.isnan(v):
                tmap.append([i,j,v])

sig = []  # list of points [ra,dec,sigma]
for i in range(s[1]):
    for j in range(s[2]):
        for si in img[8:12,i,j]:
            if not np.isnan(si):
                sig.append([si])

isoiso = []
isogrp = []
L = len(lop2)
while lop2:
    nxt = []
    nxt.append(lop2.pop(0))
    rdlop2.pop(0)
    p = 0  
    #while p < len(nxt):
    while p < len(nxt):   
        #tree = spatial.KDTree(rdlop)
        #print 'len(rdlop)',len(rdlop)
        #print 'len(lop)',len(lop)
        tree = spatial.KDTree(rdlop2)
        ind = tree.query_ball_point([nxt[p][0], nxt[p][1]], r)
        tlop = []
        for i in range(len(lop2)):
            tlop.append([lop2[i][0], lop2[i][1], lop2[i][2]])
        #print 'len(tlop)',len(tlop)
        if len(tlop) <= 1:
            break
        tind = tmap.index([ nxt[p][0], nxt[p][1], nxt[p][2] ])
        dv = sig[tind][0]
        if dv >=0.015:
            dv = 0.012
        if len(ind)>0:
            #print 'len(lop)0',len(lop)
            #print 'len(tlop)0',len(tlop)
            for n in ind:
                #if n > len(tlop):
                    #continue
                if inner(nxt[p][0],nxt[p][1],nxt[p][2],tlop[n][0],tlop[n][1],tlop[n][2],r,dv):
                    tempind = lop2.index([tlop[n][0],tlop[n][1],tlop[n][2]])
                    nxt.append(lop2.pop(tempind))
                    rdlop2.pop(tempind)
                    print 'len(lop2)s3',len(lop2)
                    #print 'len(tlop)1',len(tlop)
                    #print 'n', n
        p += 1
        #print 'p', p
        if len(lop2) <= 1:
                break
    if len(nxt) < min_num_group:
        isoiso.append(nxt)
    else:
        isogrp.append(nxt)

    print('%.3f%%'%(100*(1-len(lop2)/L)))

    if len(lop2) <= 1:
        break

stoptime = timeit.default_timer()
print('Time: ', stoptime - starttime)










# analysis ====================================================================

l_grp = len(grp1)
n_grp = np.zeros(l_grp,dtype=np.int)  # number of points in each group
for i in range(l_grp):
    n_grp[i] = len(grp1[i])
N_grp = sum(n_grp)  # number of points in 'points'
N_iso = 0
n_iso = np.zeros(len(iso),dtype=np.int)
for i in range(len(iso)):
    n_iso[i] = len(iso[i])
    N_iso += len(iso[i])



# output ======================================================================

# to screen -------------------------------------------------------------------

#print('dv = %.4f km/s/pc = %.4f km/s/pix'%(dv0,dv))
#print('r = %.4f pc = %.4f pix'%(r0,r))

loppp = []  # list of points [ra,dec,v]
for i in range(s[1]):
    for j in range(s[2]):
        for v in img[0:3,i,j]:
            if not np.isnan(v):
                loppp.append([i,j,v])
LL = len(loppp)
# components analysis
print('Found %d groups.'%l_grp)
print('%d points are in groups, ratio = %.1f%%'%(N_grp,N_grp/LL*100))
print('%d first step isolated points, ratio = %.1f%%'%(N_iso,N_iso/LL*100))
print('%d second step isolated points, ratio = %.1f%%'%(len(lop1),len(lop1)/LL*100))

print('%d second step isolated re-friend group points, ratio = %.1f%%'%(len(isoiso),len(isoiso)/LL*100))
print('%d third step isolated points, ratio = %.1f%%'%(len(lop2),len(lop2)/LL*100))

'''
# .p --------------------------------------------------------------------------
'''
output = '_kdt_adaptvel_rebin_nov03_0.045pcI0.04_final.p'
os.system('rm -rf ./variable/'+'*'+output)

pickle.dump(grp,open('./variable/grp'+output,'wb'))
pickle.dump(iso,open('./variable/iso'+output,'wb'))

pickle.dump(grp1,open('./variable/grp1'+output,'wb'))

pickle.dump(isogrp,open('./variable/isogrp'+output,'wb'))
pickle.dump(isoiso,open('./variable/isoiso'+output,'wb'))
#'''


'''
stop = timeit.default_timer()
print('Time: ', stop - start)








# .pdf ------------------------------------------------------------------------
'''
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
'''
ra = []
dec = []
v = []
for i in range(l_grp):
    ra.append([])
    dec.append([])
    v.append([])
    for j in range(n_grp[i]):
        dec[i].append(grp[i][j][0])
        ra[i].append(grp[i][j][1])
        v[i].append(grp[i][j][2])


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




''' #inspect specifi filament
ra = []
dec = []
v = []
i = 0
for j in range(len(grp[i])):
        ra.append(grp[i][j][0])
        dec.append(grp[i][j][1])
        v.append(grp[i][j][2])


color_set = ['b','g','r','c','m','y']
fig = plt.figure()  
ax = fig.gca(projection='3d')   
plt.axis('equal')
ax.set_xlabel('l_ra[pc]')
ax.set_ylabel('l_dec[pc]')
ax.set_zlabel('v[km/s]')
ax.scatter(ra,dec,v,s=2,edgecolors='none',c='blue')
#cset = ax.contourf(ra,dec,v, zdir='z', offset=-7.5, cmap='jet', vmin=-6, vmax=-0.7)
plt.title('3d groups, number of groups: %d'%l_grp)
plt.show()
'''



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
plt.title('2D groups, number of groups: %d'%l_grp)
plt.savefig('image/groups_pp.pdf')
plt.show()
#'''

''' 2d groups pv
color_set = ['b','g','r','c','m','y']
plt.figure()
for i in range(len(grp)):
    c = color_set[divmod(i,len(color_set))[1]]
    plt.scatter(dec[i],v[i],s=2,edgecolors='none',c=c)
plt.title('2D groups pv, number of groups: %d'%l_grp)
plt.savefig('image/groups_pv.pdf')
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
plt.title('2d isolated points, number of iso: %d'%len(iso))
plt.savefig('image/isolated.pdf')
plt.show()
#'''













