import numpy
import os
import sys

#casapy --nologger -c make_empty_image.py in.ms out.image 1024 1.5arcsec

el=len(sys.argv)
print sys.argv[:]

mslistfile= sys.argv[el-4]
imageoutlistfile  = sys.argv[el-3]
imsizeplistfile   = sys.argv[el-2]
cell      = str(sys.argv[el-1])


print mslistfile, imageoutlistfile, imsizeplistfile, cell

mslist = numpy.load(mslistfile)
imageoutlist = numpy.load(imageoutlistfile)
imsizeplist = numpy.load(imsizeplistfile)


for i in range(len(mslist)):
    
    ms = mslist[i]
    imageout = imageoutlist[i]
    imsizep = numpy.int(numpy.float(imsizeplist[i]))
    
    os.system('rm -rf ' + imageout)
    
    print ms, imageout, imsizep

    im.open(ms)

    imsize    = [imsizep, imsizep]
    #cell      = '1.5arcsec'

    im.defineimage(nx=imsize[0],ny=imsize[1],cellx=cell,celly=cell,stokes='I')


    im.make(imageout)	

    im.done()
    im.close()
    
os.system('rm -rf '+mslistfile)
os.system('rm -rf '+imageoutlistfile)
os.system('rm -rf '+imsizeplistfile)

