import matplotlib
#matplotlib.use('GTK')
import numpy
import os
import sys
from scipy import interpolate
import time
from subprocess import Popen, PIPE
import pyrap.tables as pt
import pyrap.images
import lofar.parmdb
from coordinates_mode import *
import pwd

# location of this script - note this script should be in subdirectory 'use'
#SCRIPTPATH = os.path.dirname(sys.argv[0])
SCRIPTPATH = os.path.dirname(os.path.abspath(sys.argv[0]))


pi = numpy.pi

def create_phaseshift_parset_field(msin, msout):
  ndppp_parset = msin.split('.')[0] +'ndppp_avgphaseshift_check.parset' 
  os.system('rm -f ' + ndppp_parset)
  os.system('rm -rf ' + msout)
  f=open(ndppp_parset, 'w')
  f.write('msin ="%s"\n' % msin) 
  f.write('msin.datacolumn = SUBTRACTED_DATA_ALL\n')
  f.write('msin.autoweight = false\n')
  f.write('msout ="%s"\n' % msout)
  f.write('msout.writefullresflag=False\n')
  f.write('steps = [uv,avg1]\n')
  #f.write('shift.type        = phaseshift\n')
  #f.write('shift.phasecenter = [%s]\n' % direction)
  
  f.write('uv.type=uvwflagger\n')
  f.write('uv.uvmmax=2500.0\n')
  #f.write('uv.uvmmin=20.0\n')
  f.write('avg1.type = squash\n')
  f.write('avg1.freqstep = 20\n')
  f.write('avg1.timestep = 6\n')    
  f.close()
  return ndppp_parset
  
  
  
el=len(sys.argv)
print el

mslist    = sys.argv[1:el-2]
res_val   = numpy.float(str(sys.argv[el-2]))
source    = str(sys.argv[el-1])

msavglist = []
for ms_id, ms in enumerate(mslist):
  msavglist.append(ms.split('.')[0] + '.' + source + '.ms.avgcheck')
     

username = pwd.getpwuid(os.getuid())[0]





###########################################################################  
# NDPPP phase shift, less averaging (NEW: run 2 in parallel)
for ms_id, ms in enumerate(mslist):
     parset = create_phaseshift_parset_field(ms, msavglist[ms_id])

     cmd = "ps -u " + username + " | grep NDPPP | wc -l"
     output=numpy.int(Popen(cmd, shell=True, stdout=PIPE).communicate()[0])
     while output > 1 : # max 2 processes (max 2, instead of 3, because we also phaseshift) 
       time.sleep(10)
       output=numpy.int(Popen(cmd, shell=True, stdout=PIPE).communicate()[0])

     # START NDPPP BECAUSE LESS/EQ 2 PROCESSES ARE RUNNING	
     os.system('NDPPP ' + parset + '&')

# Check if all NDPPP processes are finished
output=numpy.int(Popen(cmd, shell=True, stdout=PIPE).communicate()[0]) 
while output > 0 :
      time.sleep(10)
      output=numpy.int(Popen(cmd, shell=True, stdout=PIPE).communicate()[0])
###########################################################################  


imsize = 2048


###########################################################################  
# IMAGE IN PARALLEL
for ms in msavglist:

     #cmd = "ps -fuww "+ username +" | grep 'casapy_cleanv4_checksubtract.py' | wc -l"
     cmd = "ps -fu "+ username +" | grep 'python -W ignore' | wc -l"
     output=numpy.int(Popen(cmd, shell=True, stdout=PIPE).communicate()[0])
     
     while output > ((8*2)+2) : # max 8 in parallel
       time.sleep(10)
       output=numpy.int(Popen(cmd, shell=True, stdout=PIPE).communicate()[0])
 
     imout = 'im'+ '_residual_' + source + '_' + ms.split('.')[0]
     os.system('casapy --nologger -c '+SCRIPTPATH+'/casapy_cleanv4_checksubtract.py ' +\
                ms + ' ' + imout + ' ' + str(imsize)+ '&')
     time.sleep(20)


# Check if all NDPPP processes are finished
time.sleep(20)
output=numpy.int(Popen(cmd, shell=True, stdout=PIPE).communicate()[0])
while output > 2 : # because this generates two procceses matching the "cmd", see line below
  #weeren    2855  1102  0 03:57 pts/1    00:00:00 /bin/sh -c ps -fu weeren | grep 'python -W ignore' 
  #weeren    2857  2855  0 03:57 pts/1    00:00:00 grep python -W ignore
  time.sleep(10)
  output=numpy.int(Popen(cmd, shell=True, stdout=PIPE).communicate()[0])
########################################################################### 


#conver the images to FITS format
for ms in msavglist:
  imout = 'im'+ '_residual_' + source + '_' + ms.split('.')[0]
  os.system('image2fits in=' + imout +'.image'     + ' ' + 'out='+ imout + '.fits')






stopcal = False
for ms in msavglist:
    
    
    # find the source which was done before the current one
    findim   = 'im'+ '_residual_' + '*' + '_' + ms.split('.')[0] + '.image'
    cmd      = 'ls -dt1 ' + findim
    output   = Popen(cmd, shell=True, stdout=PIPE).communicate()[0].split()
    pre_sourcename = 'empty'
    if len(output) > 1:
        pre_sourcename = output[1]
        #print 'Previous image was', pre_sourcename
    
    image = 'im'+ '_residual_' + source + '_' + ms.split('.')[0] + '.image'

    img    = pyrap.images.image(image)
    pixels = numpy.copy(img.getdata())
    maxval = numpy.copy(numpy.max(pixels))
    
    maxvalpre = 1e9
    if pre_sourcename != 'empty' :
      imgpre    = pyrap.images.image(pre_sourcename)
      pixelspre = numpy.copy(imgpre.getdata())
      maxvalpre = numpy.copy(numpy.max(pixelspre))
    
   
    print maxval, ' ' + image
    print maxvalpre, ' ' + pre_sourcename
    if  (maxval > res_val) or ((maxval*0.95) > maxvalpre) :
      stopcal = True
      print 'WARNING RESIDUAL TOO LARGE, STOPPING', maxval, res_val
      print 'WARNING RESIDUAL TOO LARGE, STOPPING, previous max in image', maxvalpre
      
while(stopcal):      
  time.sleep(100)
   
      
