import numpy
import os
import sys

#casapy --nologger -c make_empty_image.py in.ms out.image 1024 1.5arcsec

el=len(sys.argv)
print sys.argv[:]

ms        = sys.argv[el-4]
imageout  = sys.argv[el-3]
imsizep   = numpy.int(sys.argv[el-2]) 
cell      = str(sys.argv[el-1])


print ms, imageout, imsizep, cell

os.system('rm -rf ' + imageout)

im.open(ms)

imsize    = [imsizep, imsizep]
#cell      = '1.5arcsec'

im.defineimage(nx=imsize[0],ny=imsize[1],cellx=cell,celly=cell,stokes='I')


im.make(imageout)	

im.done()
im.close()
