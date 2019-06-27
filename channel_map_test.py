import matplotlib.animation as animation
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy import wcs
import numpy as np
import aplpy

# channel map of CH3OH overlaid on continuum in contours.

plt.rcParams['axes.labelsize']=14
plt.rcParams['xtick.labelsize']=12
plt.rcParams['ytick.labelsize']=12
#plt.rcParams['font.family'] ='serif'

vmax = 8.e-3
string = 'CH3OH'

filename = "../data/G286_"+string+".test2.fits"
output = 'fch_map_'+string+'.pdf'

source=(159.63396,-58.31902)
size = 0.0012
rms_cont_high = 2.0e-5

data, hdr = fits.getdata('tmpcut.fits', 0, header = True) 
daplt1 = fits.PrimaryHDU(data=data,header=hdr)

def prep(name,vellim):

  data,hdr=fits.getdata(name,0,header=True)
  nxline=hdr['NAXIS1']
  nyline=hdr['NAXIS2']
  nvel=hdr['NAXIS3']
  delv=hdr['CDELT3']*1.e-3
  repro_data=np.ndarray(shape=(1,nyline,nxline))
  w=wcs.WCS(hdr)
  #w.wcs_world2pix(ra(in degree),dec(in degree),vel(in m/s depending the header),0) the ra, dec must be in the data range.
  temp = w.wcs_world2pix(source[0],source[1],vellim[0]*1.e3,0)[2]
  index0 = int(round(float(temp)))
  temp = w.wcs_world2pix(source[0],source[1],vellim[1]*1.e3,0)[2]
  index1 = int(round(float(temp)))
  temp=data[index0:index1+1,:,:]
  repro_data[0,:,:]=temp.sum(0)*abs(delv)
  hdr['NAXIS3']=1
  hdr['CDELT3']=1
  daplt=fits.PrimaryHDU(data=repro_data,header=hdr)
  return daplt


fig = plt.figure(figsize=(8,8)) 
vel_arr = np.arange(-28,-10,0.5)

def update(index):
  vel = vel_arr[index]
  daplt=prep(filename,(vel,vel+0.5))

  f = aplpy.FITSFigure(daplt,figure=fig,subplot=[0.18,0.13,.8,0.8])
  f.recenter(source[0],source[1],width=size,height=size)
  f.show_colorscale(aspect='auto',vmin=0,vmax=vmax,cmap='Blues')
  levs = rms_cont_high*np.array([4,6,10,15,20,30,50,100])
  f.show_contour(daplt1,levels=levs,colors='green',linewidths=0.6)

  f.ticks.set_xspacing(3/3600)
  f.add_scalebar(0.01*206265/2500/3600)
  f.scalebar.set_label('0.01pc')
  f.scalebar.set_color('k')
  f.set_tick_color('k') 

  channel_tag = float('%0.1f'%(vel+0.25))            
  f.add_label(0.2,0.9,str(channel_tag)+' km/s',color='k',relative=True,size='xx-large')

ani = animation.FuncAnimation(fig,update,interval=400,repeat_delay=1000,frames=34)
#ani.save('test.mp4')

# you may need to 'brew install imagemagick'...
ani.save('test.gif',writer='imagemagick', fps=2)
