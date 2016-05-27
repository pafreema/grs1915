from astropy.io import fits
import scipy.ndimage as nd
from spectral_cube import SpectralCube
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np

tmax=fits.getdata('H2CO_303_202_out.fits')
#get tmax map to find brightest pixel

mask=(tmax==0)
#mask where tmax is 0
#nd - multidimensional imaging. non-zero (true) elements of the mask are dilated, np.ones are used to dilate
data=~nd.morphology.binary_dilation(mask, np.ones((71, 71)))
# ~ reverses the true and false elements
#plt.imshow(data*tmax) shows tmax map with edges removed (got rid of noise at edges which would interfere with finding max of the signal. would have to adjust size of np.ones based off image)

np.max(data*tmax)
#finds max value
np.argmax(data*tmax)
#finds position of max - 1d
np.unravel_index(np.argmax(data*tmax), tmax.shape)
#gives back pixels in (x,y) of max, then assign ymax and xmax to those values
ymax, xmax=np.unravel_index(np.argmax(data*tmax), tmax.shape)

velocity=np.arange(35,95)
#make an array of velocity values as if you only have a single command in plt it starts x at 0

a=SpectralCube.read('H2CO_303_202.fits')
a[:, ymax, xmax]
#slice so you get entire spectrum at specific x and y, returns array of values
#plt.plot(a[:, ymax, xmax].value) gives you spectrum 
a1=a.to(u.K, equivalencies=a.beam.jtok_equiv(a.header['RESTFRQ']*u.Hz))
#gets a in K instead of Jy, need equivalencies as K~Jy/beam

plt.plot(velocity, a1[:, ymax, xmax].value)
plt.xlabel(r'$v$ (km s$^{-1}$)')
plt.ylabel(r'$T_b$ (K)')
plt.savefig('H2CO_303_202_brightspectrum.png')


#now find spectra for other lines of H2CO at same position
plt.clf()
b=SpectralCube.read('H2CO_322_221.fits')
b1=b.to(u.K, equivalencies=b.beam.jtok_equiv(b.header['RESTFRQ']*u.Hz))
plt.plot(velocity, b1[:, ymax, xmax].value)
plt.xlabel(r'$v$ (km s$^{-1}$)')
plt.ylabel(r'$T_b$ (K)')
plt.savefig('H2CO_322_221_brightspectrum.png')

plt.clf()
c=SpectralCube.read('H2CO_321_220.fits')
c1=c.to(u.K, equivalencies=c.beam.jtok_equiv(c.header['RESTFRQ']*u.Hz))
plt.plot(velocity, c1[:, ymax, xmax].value)
plt.xlabel(r'$v$ (km s$^{-1}$)')
plt.ylabel(r'$T_b$ (K)')
plt.savefig('H2CO_321_220_brightspectrum.png')


