import os
import sys
import numpy

#casapy --nologger -c ~weeren/scripts/rx42_hba/casapy_cleanv2.py ms1 ms2 imout maskname 1mJy 1000 1 512 True

el=len(sys.argv)

print sys.argv[:]


ms        = sys.argv[4:el-2]
imageout  = sys.argv[el-2]
imsizep   = numpy.int(sys.argv[el-1])


#timer = '>15:32:00'
timer = '<23:50:00' # to fix instabilities at the end of a 20 min block
timer = ''

nfacets = 1
wplanes = 1

if imsizep > 512:
   wplanes = 64

if imsizep > 799:
   wplanes = 72

if imsizep > 1023:
   wplanes = 96

if imsizep > 1599:
   wplanes = 128
   nfacets = 1

if imsizep > 2047:
  wplanes = 164
  nfacets = 1

if imsizep > 3000:
  wplanes = 256
  nfacets = 1

if imsizep > 4095:
  wplanes = 448
  nfacets  = 1

print ms, imageout


imsize    = [imsizep, imsizep]
cell      = ['30arcsec', '30arcsec']
uvrange   =">80lambda"









clean(vis=ms,imagename=imageout,outlierfile="",field="",spw="",selectdata=True,timerange=timer,\
      uvrange=uvrange,antenna="",scan="",observation="",mode="mfs",gridmode="widefield",wprojplanes=wplanes,\
      facets=1,cfcache="cfcache.dir",painc=360.0,epjtable="",interpolation="linear",        \
      niter=10,gain=0.1,threshold="0Jy",psfmode="clark",imagermode="csclean",        \
      ftmachine="mosaic",mosweight=False,scaletype="SAULT",multiscale=[],negcomponent=-1,   \
      smallscalebias=0.6,interactive=False,mask="",nchan=-1,start=0,width=1,outframe="",  \
      veltype="radio",imsize=imsize,cell=cell,phasecenter="",restfreq="",stokes="I",        \
      weighting="briggs",robust=0.0,uvtaper=False,outertaper=[''],innertaper=['1.0'],     \
      modelimage="",restoringbeam=[''],pbcor=False,minpb=0.2,usescratch=False,noise="1.0Jy",\
      npixels=0,npercycle=100,cyclefactor=1.5,cyclespeedup=-1,nterms=1,reffreq="",          \
      chaniter=False,flatnoise=True,allowchunk=False)
