from astropy.io import fits
from spectral_cube import SpectralCube
from astropy.wcs.utils import add_stokes_axis_to_wcs
import radio_beam

#shortcuts to get the right directory
pa='/mnt/fastdata/pafreema/'
tpp='/mnt/fastdata/pafreema/totalpowerproduct/'
tppr='/mnt/fastdata/pafreema/totalpowerproduct/Radio_Peak_of_IRAS_19132+1035'

#remove all old files - beware of what line
os.system("rm -rf /mnt/fastdata/pafreema/13CO_fullchannel.*")
os.system("rm -rf /mnt/fastdata/pafreema/totalpowerproduct/Radio_Peak_of_IRAS_19132+1035_regrid9.spw33.I.sd.image")
os.system("rm -rf /mnt/fastdata/pafreema/13CO_fullchannelflux9.fits")
os.system("rm -rf /mnt/fastdata/pafreema/totalpowerproduct/Radio_Peak_of_IRAS_19132+1035_regrid9.spw33.I.sd.fits")
os.system("rm -rf /mnt/fastdata/pafreema/totalpowerproduct/fluxmultregrid9.spw33.I.sd.fits")
os.system("rm -rf /mnt/fastdata/pafreema/totalpowerproduct/fluxmultregrid9.spw33.I.sd.image")
os.system("rm -rf /mnt/fastdata/pafreema/totalpowerproduct/fluxmultregrid09.spw33.I.sd.image")
os.system("rm -rf /mnt/fastdata/pafreema/13CO_combinewithtp.*")

#shortcuts that will be specific for each spw/line
line='13CO'
spw='33' #single dish spw

myspw='8' #spw for array
myniter=1
mythreshold='8mJy'
mynchan=-1
mystart=''
mywidth=1
myimsize=800
mycell='0.2arcsec'
myrestfreq='220.39868GHz'

#create all channel template image

clean(vis=["/mnt/bigdata/pafreema/calibrated_final.ms", "/mnt/bigdata/pafreema/calib_MS_7m/calibrated_final.ms"],imagename=pa+line+'_fullchannel',outlierfile="",field=["4~36", "5~18"],spw=myspw, selectdata=True,timerange="",uvrange="",antenna="",scan="",observation="",intent="", mode="channel",resmooth=False,gridmode="",wprojplanes=-1,facets=1,cfcache="cfcache.dir", rotpainc=5.0,painc=360.0,aterm=True,psterm=False,mterm=True,wbawp=False,conjbeams=True, epjtable="",interpolation="linear", niter=myniter,gain=0.1,threshold=mythreshold,psfmode="clark",imagermode="csclean", ftmachine="mosaic",mosweight=False,scaletype="SAULT",multiscale=[], negcomponent=-1,smallscalebias=0.6, interactive=False,mask=[], nchan=mynchan,start=mystart,width=mywidth,outframe="lsrk",veltype="radio", imsize=[myimsize, myimsize],cell=mycell, phasecenter="J2000 19h15m38.305s 10d41m01.018s", restfreq=myrestfreq, stokes="I",weighting="natural",robust=0.0,uvtaper=False,outertaper=[''],innertaper=['1.0'],modelimage='', restoringbeam=[''], pbcor=True,minpb=0.2,usescratch=False,noise="1.0Jy", npixels=0,npercycle=100,cyclefactor=1.5,cyclespeedup=-1,nterms=1,reffreq="",chaniter=False, flatnoise=True,allowchunk=False)

os.system("rm -rf /mnt/fastdata/pafreema/13CO_fullchan9.residual /mnt/fastdata/pafreema/13CO_fullchan9.model /mnt/fastdata/pafreema/13CO_fullchan9.psf")

#using the already made [800,800] 0.2 arcsec _fullchannel.* templates
#regrid the sd image to the template image of the line
imregrid(imagename=tppr+'.spw'+spw+'.I.sd.fits', template=pa+line+'_fullchan9.image', output=tppr+'_regrid9.spw'+spw+'.I.sd.image', asvelocity=False, overwrite=True)


#export the template .flux file and the regridded sd file to fits

exportfits(imagename=pa+line+'_fullchan9.flux', fitsimage=pa+line+'_fullchannelflux9.fits', overwrite=True)
exportfits(imagename=tppr+'_regrid9.spw'+spw+'.I.sd.image', fitsimage=tppr+'_regrid9.spw'+spw+'.I.sd.fits', overwrite=True)

#with spectralcube multiply the flux and sd
flux=SpectralCube.read(pa+line+'_fullchannelflux9.fits')
flux.allow_huge_operations=True
sd=SpectralCube.read(tppr+'_regrid9.spw'+spw+'.I.sd.fits')
sd.allow_huge_operations=True
end=flux*sd

#need to add a 4th axis - Stokes and modify the headers to match the ones of the original sd file
hdu = fits.PrimaryHDU(end.filled_data[:].value.reshape((1,)+end.shape), header=add_stokes_axis_to_wcs(end.wcs, 3).to_header())

a=SpectralCube.read(tppr+'.spw'+spw+'.I.sd.fits')

hdu.header['BUNIT']='Jy/beam'
hdu.header['BMAJ']=a.header['BMAJ']
hdu.header['BMIN']=a.header['BMIN']
hdu.header['BPA']=a.header['BPA']
hdu.header['BZERO']=a.header['BZERO']
hdu.header['BSCALE']=a.header['BSCALE']
hdu.header['BTYPE']=a.header['BTYPE']
hdu.header['TELESCOP']=a.header['TELESCOP']
hdu.header['INSTRUME']=a.header['INSTRUME']
hdu.header['OBSERVER']=a.header['OBSERVER']
hdu.header['TIMESYS']=a.header['TIMESYS']
hdu.header['OBSDEC']=a.header['OBSDEC']
hdu.header['OBSRA']=a.header['OBSRA']
hdu.header['VELREF']=a.header['VELREF']

hdu.writeto(tpp+'fluxmultregrid9.spw'+spw+'.I.sd.fits', clobber=True)

#get back to a good casa image by imtrans twice
imtrans(imagename=tpp+'fluxmultregrid9.spw'+spw+'.I.sd.fits', outfile=tpp+'fluxmultregrid9.spw'+spw+'.I.sd.image', order='0132')

imtrans(imagename=tpp+'fluxmultregrid9.spw'+spw+'.I.sd.image', outfile=tpp+'fluxmultregrid09.spw'+spw+'.I.sd.image', order='0132')

#clean using this new file as a modelimage

myspw='8' #spw for array
myniter=20000
mythreshold='8mJy'
mynchan=60
mystart=40
mywidth=1
myimsize=800
mycell='0.2arcsec'
myrestfreq='220.36868GHz'

#12CO - spw 17 - myspw 0 - freq 230.538GHz - mynchan 80 mystart 250
#13CO - spw 33 - myspw 8 - freq 220.39868GHz - mynchan 60 mystart 50
#C18O - spw 31 - myspw 7 - freq 219.56035GHz - mynchan 60 mystart 160
#303-202 - spw 25 - myspw 4 - freq 218.22219GHz - mynchan 100 mystart 400 (full clean nchan 480 start 240)
#322-221 - spw 27 - myspw 5 - freq 218.47563GHz - mynchan 200 mystart 430 (full clean nchan 480 start 240)
#321-220 - spw 29 - myspw 6 - freq 218.76007GHz - mynchan 100 mystart 400 (full clean nchan 480 start 240)
#H30a - spw 21 - myspw 2 - freq 231.90093GHz - mynchan 80 mystart 80
#SiO - spw 23 - myspw 3 - freq 217.105 - (full clean nchan 480 start 240)


clean(vis=["/mnt/bigdata/pafreema/calibrated_final.ms", "/mnt/bigdata/pafreema/calib_MS_7m/calibrated_final.ms"],imagename=line+"_combinewithtp",outlierfile="",field=['4~36', '5~18'],spw=myspw, selectdata=True,timerange="",uvrange="",antenna="",scan="",observation="",intent="",  mode="channel",resmooth=False,gridmode="",wprojplanes=-1,facets=1,cfcache="cfcache.dir", rotpainc=5.0,painc=360.0,aterm=True,psterm=False,mterm=True,wbawp=False,conjbeams=True, epjtable="",interpolation="linear", niter=myniter,gain=0.1,threshold=mythreshold,psfmode="clark",imagermode="csclean", ftmachine="mosaic",mosweight=False,scaletype="SAULT",multiscale=[0, 1, 5, 10, 15], negcomponent=-1,smallscalebias=0.6, interactive=False,mask=[], nchan=mynchan,start=mystart,width=mywidth,outframe="lsrk",veltype="radio", imsize=[myimsize, myimsize],cell=mycell, phasecenter="J2000 19h15m38.305s 10d41m01.018s", restfreq=myrestfreq, stokes="I",weighting="natural",robust=0.0,uvtaper=False,outertaper=[''],innertaper=['1.0'],modelimage=tpp+'fluxmultregrid1.spw'+spw+'.I.sd.image', restoringbeam=[''], pbcor=True,minpb=0.2,usescratch=False,noise="1.0Jy", npixels=0,npercycle=100,cyclefactor=2.0,cyclespeedup=-1,nterms=1,reffreq="",chaniter=False, flatnoise=True,allowchunk=False)

