#!/usr/bin/python

# blank a facet FITS file using a mask image

import sys
import pyfits
import numpy as np
import pyrap.images

def blank_facet(imagename,maskname):

    imhdu=pyfits.open(imagename)
    maskim=pyrap.images.image(maskname)

    imdata=imhdu[0].data[0,0]
    maskdata=maskim.getdata()[0,0]

    assert(imdata.shape==maskdata.shape)

    nanmask=np.ones_like(imdata)*np.nan
    imdata=np.where(maskdata>0,imdata,nanmask)
    imhdu[0].data[0,0]=imdata
    outname=imagename.replace('.fits','.blanked.fits')
    imhdu.writeto(outname,clobber=True)
    return outname

if __name__=='__main__':
    imagename=sys.argv[1]
    maskname=sys.argv[2]
    blank_facet(imagename,maskname)
