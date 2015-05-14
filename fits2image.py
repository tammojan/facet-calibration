import os
import sys
import numpy

el=len(sys.argv)
print sys.argv[:]

inputim     = sys.argv[el-2]
imageout  = sys.argv[el-1]


print 'inputim', inputim
print 'imageout', imageout

os.system('rm -rf ' + imageout)

importfits(fitsimage=inputim,imagename=imageout,whichrep=0,\
           whichhdu=-1,zeroblanks=True,overwrite=True,defaultaxes=False,\
           defaultaxesvalues=[],beam=[])
