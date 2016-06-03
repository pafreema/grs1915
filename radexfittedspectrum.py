import pyspeckit
import astropy.io.fits as pyfits
from pyspeckit.spectrum import models
from spectral_cube import SpectralCube
import astropy.units as u
import numpy as np
import scipy.ndimage as nd

#read in model grids to create H2CO radex fitter
#grids use the radex program to calculate brightness of H2CO under different conditions
texgrid1a=pyfits.getdata('/mnt/fastdata/pafreema/303-202_321-220_5kms_temperature_para_tex1.fits')
taugrid1a=pyfits.getdata('/mnt/fastdata/pafreema/303-202_321-220_5kms_temperature_para_tau1.fits')
texgrid2a=pyfits.getdata('/mnt/fastdata/pafreema/303-202_321-220_5kms_temperature_para_tex2.fits')
taugrid2a=pyfits.getdata('/mnt/fastdata/pafreema/303-202_321-220_5kms_temperature_para_tau2.fits')
hdra=pyfits.getheader('/mnt/fastdata/pafreema/303-202_321-220_5kms_temperature_para_tau2.fits')

texgrid1b=pyfits.getdata('/mnt/fastdata/pafreema/303-202_322-221_5kms_temperature_para_tex1.fits')
taugrid1b=pyfits.getdata('/mnt/fastdata/pafreema/303-202_322-221_5kms_temperature_para_tau1.fits')
texgrid2b=pyfits.getdata('/mnt/fastdata/pafreema/303-202_322-221_5kms_temperature_para_tex2.fits')
taugrid2b=pyfits.getdata('/mnt/fastdata/pafreema/303-202_322-221_5kms_temperature_para_tau2.fits')
hdrb=pyfits.getheader('/mnt/fastdata/pafreema/303-202_322-221_5kms_temperature_para_tau2.fits')

#models.formaldehyde.formaldehyde_radex is the model that we are going to fit 
#models.model.SpectralModel is a wrapper to deal with parinfo, multiple peaks, and annotations
#all parameters after the first are passed to the model function

#fit the 303-202 and 322-221 lines
#tex/tau grid sets frequency range (in GHz) over which frequency range is valid 
formaldehyde_radex_fitter_b=models.model.SpectralModel(models.formaldehyde_mm.formaldehyde_mm_radex, 5, parnames=['temperature', 'column', 'density', 'center', 'width'], parvalues=[50,12,4.5,0,1], parlimited=[(True, True), (True, True), (True, True), (False, False), (True, False)], parlimits=[(5,205), (10,17), (2,7), (0,0), (0,0)], parsteps=[0.01, 0.01, 0.1, 0, 0], fitunits='Hz', texgrid=((218.15, 218.25, texgrid1b), (218.4, 218.55, texgrid2b)), taugrid=((218.15, 218.25, taugrid1b), (218.4, 218.55, taugrid2b)), hdr=hdrb, shortvarnames=("T", "N", "n", "v", "\\sigma"), grid_vwidth=5.0)

#fit the 303-202 and 321-220 lines
formaldehyde_radex_fitter_a=models.model.SpectralModel(models.formaldehyde_mm.formaldehyde_mm_radex, 5, parnames=['temperature', 'column', 'density', 'center', 'width'], parvalues=[50,12,4.5,0,1], parlimited=[(True, True), (True, True), (True, True), (False, False), (True, False)], parlimits=[(5,205), (10,17), (2,7), (0,0), (0,0)], parsteps=[0.01, 0.01, 0.1, 0, 0], fitunits='Hz', texgrid=((218.15, 218.25, texgrid1a), (218.7, 218.8, texgrid2a)), taugrid=((218.15, 218.25, taugrid1a), (218.7, 218.8, taugrid2a)), hdr=hdra, shortvarnames=("T", "N", "n", "v", "\\sigma"), grid_vwidth=5.0)

#fit all three lines
formaldehyde_radex_fitter_both=models.model.SpectralModel(models.formaldehyde_mm.formaldehyde_mm_radex, 5, parnames=['temperature', 'column', 'density', 'center', 'width'], parvalues=[50,12,4.5,0,1], parlimited=[(True, True), (True, True), (True, True), (False, False), (True, False)], parlimits=[(5,205), (10,17), (2,7), (0,0), (0,0)], parsteps=[0.01, 0.01, 0.1, 0, 0], fitunits='Hz', texgrid=((218.15, 218.25, texgrid1b), (218.4, 218.55, texgrid2b), (218.7, 218.8, texgrid2a)), taugrid=((218.15, 218.25, taugrid1b), (218.4, 218.55, taugrid2b), (218.7, 218.8, taugrid2a)), hdr=hdrb, shortvarnames=("T", "N", "n", "v", "\\sigma"), grid_vwidth=5.0)

#?? if __name__="__main__":
datadir = '/mnt/fastdata/pafreema/'
#get spectra, into units of K and GHz
tmax=pyfits.getdata(datadir+'H2CO_303_202_out.fits')
mask=(tmax==0)
data=~nd.morphology.binary_dilation(mask, np.ones((71,71)))
ymax, xmax=np.unravel_index(np.argmax(data*tmax), tmax.shape)

a=SpectralCube.read(datadir+'H2CO_303_202.fits')
a1=a.to(u.K, equivalencies=a.beam.jtok_equiv(a.header['RESTFRQ']*u.Hz))
a2=a1.with_spectral_unit(u.Hz)

b=SpectralCube.read(datadir+'H2CO_322_221.fits')
b1=b.to(u.K, equivalencies=b.beam.jtok_equiv(b.header['RESTFRQ']*u.Hz))
b2=b1.with_spectral_unit(u.Hz)


c=SpectralCube.read(datadir+'H2CO_321_220.fits')
c1=c.to(u.K, equivalencies=c.beam.jtok_equiv(c.header['RESTFRQ']*u.Hz))
c2=c1.with_spectral_unit(u.Hz)

#combine three lines into one spectrum by concatenating numpy arrays, y-y axis,x-x axis 
y=np.concatenate((a2[:, ymax, xmax].value, b2[:, ymax, xmax].value, c2[:, ymax, xmax].value))
x=np.concatenate((a2.spectral_axis.value, b2.spectral_axis.value, c2.spectral_axis.value))
#need to take only values in concatenation and add units with xarrkwargs below
x /= 1e9
sp=pyspeckit.Spectrum(data=y, xarr=x, unit="$T_B$ (K)",xarrkwargs={'unit':'GHz'})

#add (register) fitters to the spectrum
sp.Registry.add_fitter('formaldehyde_mm_radex', formaldehyde_radex_fitter_both, 5)
sp.Registry.add_fitter('formaldehyde_mm_radex_a', formaldehyde_radex_fitter_a, 5)
sp.Registry.add_fitter('formaldehyde_mm_radex_b', formaldehyde_radex_fitter_b, 5)

#plot fit for all 3 ('both')
sp.plotter(figure=1)
sp.specfit(fittype='formaldehyde_mm_radex', guesses=[95, 13.2, 4.5, 67, 7.0], limits=[(20,200), (11,15), (3,5.5), (65,70), (2,15)], limited=[(True, True)]*5, fixed=[False, False, False, False, False])
#only change center parameter from example
sp.plotter.savefig('H2CO_all_radexfit.pdf')

#plot fit for 303-202 and 321-220 ('a')
sp.plotter(figure=2)
sp.specfit(fittype='formaldehyde_mm_radex_a', guesses=[95,13.2,4.5,67,7.0], limits=[(20,200), (11,15), (3,5.5), (65,70), (2,15)], limited=[(True, True)]*5, fixed=[False, False, False, False, False])
sp.plotter.savefig('H2CO_303_202_321_220_radexfit.pdf')

#plot fit for 303-202 and 322-221 ('b')
sp.plotter(figure=3)
sp.specfit(fittype='formaldehyde_mm_radex_b', guesses=[95,13.2,4.5,67,7.0], limits=[(20,200), (11,15), (3,5.5), (65,70), (2,15)], limited=[(True, True)]*5, fixed=[False, False, False, False, False])
sp.plotter.savefig('H2CO_303_202_322_221_radexfit.pdf')
