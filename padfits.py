#!/usr/bin/python

# pad a FITS file

import pyfits
import sys
import numpy as np

def padfits(infile,outfile,scalefactor=1.2,verbose=False):

    hdu=pyfits.open(infile)
    imdata=hdu[0].data[0,0]

    (xsize,ysize)=imdata.shape
    assert(xsize==ysize)
    if verbose:
        print 'size is',xsize

    padsize=int(xsize*1.2)
    offset=(padsize-xsize)/2
    if verbose:
        print 'padding to',padsize
        print 'offset is',offset

    newdata=np.zeros((1,1,padsize,padsize))

    newdata[0,0,offset:offset+xsize,offset:offset+xsize]=imdata
    hdu[0].data=newdata
    hdu[0].header['CRPIX1']+=offset
    hdu[0].header['CRPIX2']+=offset
    hdu.writeto(outfile,clobber=True)
    return padsize

if __name__=='__main__':
    filename=sys.argv[1]
    outfile=sys.argv[2]
    padfits(filename,outfile,verbose=True)
