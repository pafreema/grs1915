from spectral_cube import SpectralCube
import matplotlib.pyplot as plt
import astropy.units as u
import numpy as np

#importing fits files - a=12CO, b=13CO, c=C18O
#get a subcube corresponding to region of interest
#in a [800,800] 0.2 arcsec image (x)(y) - ridge (385,400)(555,575) - non thermal region (400,420)(570,595) - not-of-interest region (490,505)(490,505), high mass sfr (340, 360)(470, 485) 



#RIDGE
x0=385
x1=400
y0=560
y1=575

datadir='/mnt/fastdata/pafreema/'
a=SpectralCube.read(datadir+'12CO_combinewithtp.fits')
a1=a.to(u.K, equivalencies=a.beam.jtok_equiv(a.header['RESTFRQ']*u.Hz))
a2=a1.with_spectral_unit(u.km/u.s, velocity_convention='radio')
suba=a2[:, y0:y1, x0:x1]

b=SpectralCube.read(datadir+'13CO_combinewithtp.fits')
b1=b.to(u.K, equivalencies=b.beam.jtok_equiv(b.header['RESTFRQ']*u.Hz))
b2=b1.with_spectral_unit(u.km/u.s, velocity_convention='radio')
subb=b2[:, y0:y1, x0:x1]

c=SpectralCube.read(datadir+'C18O_combinewithtp.fits')
c1=c.to(u.K, equivalencies=c.beam.jtok_equiv(c.header['RESTFRQ']*u.Hz))
c2=c1.with_spectral_unit(u.km/u.s, velocity_convention='radio')
subc=c2[:, y0:y1, x0:x1]

vaxisa=suba.spectral_axis
vaxisb=subb.spectral_axis
vaxisc=subc.spectral_axis

meanspeca = np.mean(suba.filled_data[:],axis=(1,2))
meanspecb = np.mean(subb.filled_data[:],axis=(1,2))
meanspecc = np.mean(subc.filled_data[:],axis=(1,2))

fig=plt.figure(figsize=(6,6))
ax=fig.add_subplot(111)
ax.spines['top'].set_color('none')
ax.spines['bottom'].set_color('none')
ax.spines['left'].set_color('none')
ax.spines['right'].set_color('none')
ax.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')
ax.set_ylabel('$T_B$ (K)')
ax.set_xlabel('$v$ (km/s)')

ax1=fig.add_subplot(221)
plt.plot(vaxisa, meanspeca, label='$^{12}$CO', color='black', ls='-.')
plt.plot(vaxisb, meanspecb, label='$^{13}$CO', color='blue')
plt.plot(vaxisc, meanspecc, label='C$^{18}$O', color='red', ls='--')
plt.setp(ax1.get_xticklabels(), visible=False)
plt.axis([55,85, -5, 25])
plt.text(80, 22, '(a)')


#NON THERMAL
x0=400
x1=420
y0=570
y1=595

datadir='/mnt/fastdata/pafreema/'
a=SpectralCube.read(datadir+'12CO_combinewithtp.fits')
a1=a.to(u.K, equivalencies=a.beam.jtok_equiv(a.header['RESTFRQ']*u.Hz))
a2=a1.with_spectral_unit(u.km/u.s, velocity_convention='radio')
suba=a2[:, y0:y1, x0:x1]

b=SpectralCube.read(datadir+'13CO_combinewithtp.fits')
b1=b.to(u.K, equivalencies=b.beam.jtok_equiv(b.header['RESTFRQ']*u.Hz))
b2=b1.with_spectral_unit(u.km/u.s, velocity_convention='radio')
subb=b2[:, y0:y1, x0:x1]

c=SpectralCube.read(datadir+'C18O_combinewithtp.fits')
c1=c.to(u.K, equivalencies=c.beam.jtok_equiv(c.header['RESTFRQ']*u.Hz))
c2=c1.with_spectral_unit(u.km/u.s, velocity_convention='radio')
subc=c2[:, y0:y1, x0:x1]

vaxisa=suba.spectral_axis
vaxisb=subb.spectral_axis
vaxisc=subc.spectral_axis

meanspeca = np.mean(suba.filled_data[:],axis=(1,2))
meanspecb = np.mean(subb.filled_data[:],axis=(1,2))
meanspecc = np.mean(subc.filled_data[:],axis=(1,2))

ax2=fig.add_subplot(222, sharey=ax1)
plt.plot(vaxisa, meanspeca, label='$^{12}$CO', color='black', ls='-.')
plt.plot(vaxisb, meanspecb, label='$^{13}$CO', color='blue')
plt.plot(vaxisc, meanspecc, label='C$^{18}$O', color='red', ls='--')
plt.setp(ax2.get_xticklabels(), visible=False)
plt.setp(ax2.get_yticklabels(), visible=False)
plt.axis([55,85, -5, 25])
plt.legend(bbox_to_anchor=(-0.1, 1.0), loc=8, ncol=3, borderaxespad=0.)
plt.text(80, 22, '(b)')

#HMSFR
x0=340
x1=360
y0=470
y1=485

datadir='/mnt/fastdata/pafreema/'
a=SpectralCube.read(datadir+'12CO_combinewithtp.fits')
a1=a.to(u.K, equivalencies=a.beam.jtok_equiv(a.header['RESTFRQ']*u.Hz))
a2=a1.with_spectral_unit(u.km/u.s, velocity_convention='radio')
suba=a2[:, y0:y1, x0:x1]

b=SpectralCube.read(datadir+'13CO_combinewithtp.fits')
b1=b.to(u.K, equivalencies=b.beam.jtok_equiv(b.header['RESTFRQ']*u.Hz))
b2=b1.with_spectral_unit(u.km/u.s, velocity_convention='radio')
subb=b2[:, y0:y1, x0:x1]

c=SpectralCube.read(datadir+'C18O_combinewithtp.fits')
c1=c.to(u.K, equivalencies=c.beam.jtok_equiv(c.header['RESTFRQ']*u.Hz))
c2=c1.with_spectral_unit(u.km/u.s, velocity_convention='radio')
subc=c2[:, y0:y1, x0:x1]

vaxisa=suba.spectral_axis
vaxisb=subb.spectral_axis
vaxisc=subc.spectral_axis

meanspeca = np.mean(suba.filled_data[:],axis=(1,2))
meanspecb = np.mean(subb.filled_data[:],axis=(1,2))
meanspecc = np.mean(subc.filled_data[:],axis=(1,2))

ax3=fig.add_subplot(223)
plt.plot(vaxisa, meanspeca, label='$^{12}$CO', color='black', ls='-.')
plt.plot(vaxisb, meanspecb, label='$^{13}$CO', color='blue')
plt.plot(vaxisc, meanspecc, label='C$^{18}$O', color='red', ls='--')
plt.axis([55,85, -5, 25])
plt.text(80, 22, '(c)')


#NORMAL
x0=490
x1=505
y0=490
y1=505

datadir='/mnt/fastdata/pafreema/'
a=SpectralCube.read(datadir+'12CO_combinewithtp.fits')
a1=a.to(u.K, equivalencies=a.beam.jtok_equiv(a.header['RESTFRQ']*u.Hz))
a2=a1.with_spectral_unit(u.km/u.s, velocity_convention='radio')
suba=a2[:, y0:y1, x0:x1]

b=SpectralCube.read(datadir+'13CO_combinewithtp.fits')
b1=b.to(u.K, equivalencies=b.beam.jtok_equiv(b.header['RESTFRQ']*u.Hz))
b2=b1.with_spectral_unit(u.km/u.s, velocity_convention='radio')
subb=b2[:, y0:y1, x0:x1]

c=SpectralCube.read(datadir+'C18O_combinewithtp.fits')
c1=c.to(u.K, equivalencies=c.beam.jtok_equiv(c.header['RESTFRQ']*u.Hz))
c2=c1.with_spectral_unit(u.km/u.s, velocity_convention='radio')
subc=c2[:, y0:y1, x0:x1]

vaxisa=suba.spectral_axis
vaxisb=subb.spectral_axis
vaxisc=subc.spectral_axis

meanspeca = np.mean(suba.filled_data[:],axis=(1,2))
meanspecb = np.mean(subb.filled_data[:],axis=(1,2))
meanspecc = np.mean(subc.filled_data[:],axis=(1,2))

ax4=fig.add_subplot(224)
plt.plot(vaxisa, meanspeca, label='$^{12}$CO', color='black', ls='-.')
plt.plot(vaxisb, meanspecb, label='$^{13}$CO', color='blue')
plt.plot(vaxisc, meanspecc, label='C$^{18}$O', color='red', ls='--')
plt.setp(ax4.get_yticklabels(), visible=False)
plt.axis([55,85, -5, 25])
plt.text(80, 22, '(d)')
plt.tight_layout()
plt.savefig('COspectra.pdf')
plt.show()

