#!/usr/bin/python

# check the rms of a file after removing the bright source mask. Can
# be used to see if a worse than expected solution has been obtained

import pyfits
import numpy as np
import sys
import pyrap.images as pi

def getrms(image,facetmask=None):
    """
    Return the rms of source-free regions in the image.
    image should be a FITS file of the form xxx-image.fits.
    xxxnm-fitsmask should exist.
    """
    mask=image.replace('-image.fits','nm.fitsmask')
    if facetmask is None:
        # guess facet mask name
        root=image[16:].replace('-image.fits','')
        facetmask='templatemask_'+root+'.masktmp'
   
    imhdu=pyfits.open(image)
    mkhdu=pyfits.open(mask)
    maskim=pi.image(facetmask)
    
    imdata=imhdu[0].data[0,0]
    maskdata=mkhdu[0].data[0,0]
    fmaskdata=maskim.getdata()[0,0]
    assert(imdata.shape==maskdata.shape)
    assert(imdata.shape==fmaskdata.shape)
    filter=(maskdata==0) & (fmaskdata==1)
    mdata=imdata[filter]
    return np.std(mdata)

# run main code from command line with list of imfield... FITS images.

if __name__=='__main__':
    d={}
    rms=[]
    for f in sys.argv[1:]:
        if not('nm' in f):
            r=getrms(f)
            rms.append(r)
            d[f]=r
            print f,r
    mean=np.mean(r)
    print 'mean rms is',r
    for k in d:
        if d[k]>2*mean:
            print 'image',k,'is bad'
