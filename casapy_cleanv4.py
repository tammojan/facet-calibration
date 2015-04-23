import os
import sys
import numpy

#casapy --nologger -c ~weeren/scripts/rx42_hba/casapy_cleanv2.py ms1 ms2 imout maskname 1mJy 1000 1 512 True

el=len(sys.argv)

print sys.argv[:]



msind = sys.argv.index('-c')
ms        = sys.argv[msind+2:el-7]
imageout  = sys.argv[el-7]
mask      = sys.argv[el-6]
threshold = sys.argv[el-5]
niter     = numpy.int(sys.argv[el-4])
nterms    = numpy.int(sys.argv[el-3])
imsizep   = numpy.int(sys.argv[el-2])
mscale    = sys.argv[el-1]
cycfactor = 2.5

if niter > 1500:
   cycfactor = 3.0


scales = []
if mscale == 'True':
   scales = [0,3,7,25,60,150]

#timer = '>15:32:00'
#timer = '<23:50:00' # to fix instabilities at the end of a 20 min block
timer = ""

nfacets = 1
wplanes = 1

if imsizep > 512:
   wplanes = 64

if imsizep > 799:
   wplanes = 96

if imsizep > 1023:
   wplanes = 128

if imsizep > 1599:
   wplanes = 256
   nfacets = 1

if imsizep > 2047:
  wplanes = 384
  nfacets = 1

if imsizep > 3000:
  wplanes = 448 
  nfacets = 1

if imsizep > 4095:
  wplanes = 512
  nfacets  = 1

print ms, imageout, mask, threshold, niter, nterms


imsize    = [imsizep, imsizep]
cell      = ['1.5arcsec', '1.5arcsec']
uvrange   =">80lambda"


masktmp = mask


if mask == "None":
   mask = ""


print masktmp.split(",")

if len(masktmp.split(",")) >= 2 :
  print 'Ok, we have more than one mask...'
  mask = []

  for idx in range(len(masktmp.split(","))):
    mask.append(masktmp.split(",")[idx])

print 'Using mask:', mask


clean(vis=ms,imagename=imageout,outlierfile="",field="",spw="",selectdata=True,timerange=timer,\
      uvrange=uvrange,antenna="",scan="",observation="",mode="mfs",gridmode="widefield",wprojplanes=wplanes,\
      facets=nfacets,cfcache="cfcache.dir",painc=360.0,epjtable="",interpolation="linear",        \
      niter=niter,gain=0.1,threshold=threshold,psfmode="clark",imagermode="csclean",        \
      ftmachine="mosaic",mosweight=False,scaletype="SAULT",multiscale=scales,negcomponent=-1,   \
      smallscalebias=0.6,interactive=False,mask=mask,nchan=-1,start=0,width=1,outframe="",  \
      veltype="radio",imsize=imsize,cell=cell,phasecenter="",restfreq="",stokes="I",        \
      weighting="briggs",robust=-0.25,uvtaper=False,outertaper=[''],innertaper=['1.0'],     \
      modelimage="",restoringbeam=[''],pbcor=False,minpb=0.2,usescratch=False,noise="1.0Jy",\
      npixels=0,npercycle=100,cyclefactor=cycfactor,cyclespeedup=-1,nterms=nterms,reffreq="",          \
      chaniter=False,flatnoise=True,allowchunk=False)
