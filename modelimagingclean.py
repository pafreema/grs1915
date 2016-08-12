#get necessary packages
from astropy.io import fits
from spectral_cube import SpectralCube
from astropy.wcs.utils import add_stokes_axis_to_wcs
import radio_beam

#shortcuts to get the right directory
pa='/mnt/fastdata/pafreema/'
tpp='/mnt/fastdata/pafreema/totalpowerproduct/'
tppr='/mnt/fastdata/pafreema/totalpowerproduct/Radio_Peak_of_IRAS_19132+1035'

#shortcuts that will be specific for each spw/line
line='SiO'
spw='23' #single dish spw
myspw='3' #spw for array
myniter=1
mythreshold='8mJy'
mynchan=-1
mystart=''
mywidth=1
myimsize=800
mycell='0.2arcsec'
myrestfreq='217.105GHz'

#os.system("rm -rf 13CO_combinewithtp.*")

#beware of channels, with 2 vis may have changed
#12CO - spw 17 - myspw 0 - freq 230.538GHz - threshold 8mJy
#13CO - spw 33 - myspw 8 - freq 220.36868GHz - threshold 8mJy
#C18O - spw 31 - myspw 7 - freq 219.56035GHz - threshold 7mJy
#303-202 - spw 25 - myspw 4 - freq 218.22219GHz - threshold 7mJy
#322-221 - spw 27 - myspw 5 - freq 218.47563GHz - threshold 6mJy
#321-220 - spw 29 - myspw 6 - freq 218.76007GHz - threshold 6mJy
#H30a - spw 21 - myspw 2 - freq 231.90093GHz - threshold 4mJy


clean(vis=["/mnt/bigdata/pafreema/calibrated_final.ms","/mnt/bigdata/pafreema/calib_MS_7m/calibrated_final.ms"],imagename=line+"_try",outlierfile="",field=['4~36','5~18'],spw=myspw, selectdata=True,timerange="",uvrange="",antenna="",scan="",observation="",intent="",  mode="channel",resmooth=False,gridmode="",wprojplanes=-1,facets=1,cfcache="cfcache.dir", rotpainc=5.0,painc=360.0,aterm=True,psterm=False,mterm=True,wbawp=False,conjbeams=True, epjtable="",interpolation="linear", niter=myniter,gain=0.1,threshold=mythreshold,psfmode="clark",imagermode="csclean", ftmachine="mosaic",mosweight=False,scaletype="SAULT",multiscale=[], negcomponent=-1,smallscalebias=0.6, interactive=False,mask=[], nchan=mynchan,start=mystart,width=mywidth,outframe="lsrk",veltype="radio", imsize=[myimsize, myimsize],cell=mycell, phasecenter="J2000 19h15m38.305s 10d41m01.018s", restfreq=myrestfreq, stokes="I",weighting="natural",robust=0.0,uvtaper=False,outertaper=[''],innertaper=['1.0'],modelimage='', restoringbeam=[''], pbcor=True,minpb=0.2,usescratch=False,noise="1.0Jy", npixels=0,npercycle=100,cyclefactor=1.5,cyclespeedup=-1,nterms=1,reffreq="",chaniter=False, flatnoise=True,allowchunk=False)

