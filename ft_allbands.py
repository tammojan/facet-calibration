import os
import sys
import numpy
import time

#casapy --nologger -c ft_allbands.py ms1 modimage nterms imsize

execfile('mytasks.py')

el=len(sys.argv)

print sys.argv[:]





mslist    = sys.argv[4:el-3]
modimage  = sys.argv[el-3]
ntermsi   = numpy.int(sys.argv[el-2])
imsizep   = numpy.int(sys.argv[el-1])
wplanes   = 1

if imsizep > 512:
   wplanes = 64

if imsizep > 799:
   wplanes = 96

if imsizep > 1023:
   wplanes = 128

if imsizep > 1599:
   wplanes = 256

if imsizep > 2047:
  wplanes = 384

if imsizep > 3000:
  wplanes = 448 

if imsizep > 4095:
  wplanes = 512



if ntermsi == 1:
  mod = [modimage]
if ntermsi == 2: 
  mod = [modimage+'.tt0',modimage+'.tt1'] 
if ntermsi == 3: 
  mod = [modimage+'.tt0',modimage+'.tt1',modimage+'.tt2'] 


for ms_id,ms in enumerate(mslist):
       
  if wplanes > 1:
     ftw(vis=ms,field="",spw="",model=mod,nterms=ntermsi,reffreq="",\
	 wprojplanes=wplanes,complist="",incremental=False,usescratch=True, async=False)
  else:
     ft(vis=ms,field="",spw="",model=mod,nterms=ntermsi,reffreq="",\
	  complist="",incremental=False,usescratch=True, async=False)

