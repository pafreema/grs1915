from spectral_cube import SpectralCube
import astropy.io.fits as pyfits
import astropy.units as u

twelve=SpectralCube.read('12CO_combinewithtp.fits')
twelve_1=twelve.with_spectral_unit(u.km/u.s, velocity_convention='radio')
thirteen=SpectralCube.read('13CO_combinewithtp.fits')
thirteen_1=thirteen.with_spectral_unit(u.km/u.s, velocity_convention='radio')

#interpolate 12CO cube to match the 13CO spectral axis
twelve_interp=twelve_1.spectral_interpolate(spectral_grid=thirteen_1.spectral_axis)

#mask the 13CO initially using 5 * the rms signal
tmax=pyfits.open('13CO_combinewithtp.fits')
mask=tmax[0].data>0.05
badmask=tmax[0].data==0
keep=mask*(~badmask)

#mask the 13CO again, finding where the 12CO is 0.75 * the 12CO emission 
mask_thirteen=thirteen_1.with_mask(keep[0,:,:,:])
mask_thirteen_1=mask_thirteen.with_mask(thirteen_1>(twelve_interp*0.75))
mask_thirteen_1.write('13CO_combinewithtp_masked.fits')

#find the 0th moment of this masked cube
thirteen_0=mask_thirteen_1.moment(order=0)
thirteen_0.write('13CO_moment0_masked.fits')

#finding ratio of 13CO/12CO - shows nothing?
ratio=twelve_interp/thirteen_1
ratio.write('13CO_12CO_ratio.fits')

