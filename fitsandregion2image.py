import os
import sys
import numpy

el=len(sys.argv)
#print sys.argv[:]

inputim    = sys.argv[el-3]
imageout  = sys.argv[el-2]
region    = sys.argv[el-1]

print 'inputim', inputim
print 'imageout', imageout
print 'region', region

# don't importfits here to work round CASA bug
# conversion is done in calling code instead
#importfits(fitsimage=inputim,imagename=imageout,whichrep=0,\
#           whichhdu=-1,zeroblanks=True,overwrite=True,defaultaxes=False,\
#	   defaultaxesvalues=[],beam=[])


if region != 'None':
 print 'We have a region file', region
 ia.open(imageout)
 myRGN = rg.fromtextfile(filename=region, shape=ia.shape(), csys=ia.coordsys().torecord())
 im.regiontoimagemask(mask=imageout, region=myRGN)
 ia.close()

# avoid cleaning the edges of the image
edge = 25
tb.open(imageout, nomodify=False)
data = numpy.copy(tb.getcol('map'))

print numpy.shape(data)
sh = numpy.shape(data)
print sh


# mask edges
data[0:sh[0],0:edge,0,0,0] = 0
data[0:edge,0:sh[1],0,0,0] = 0
data[0:sh[0],sh[1]-edge:sh[1],0,0,0] = 0
data[sh[0]-edge:sh[0],0:sh[1],0,0,0] = 0



tb.putcol('map', data)
tb.flush()
tb.close()
