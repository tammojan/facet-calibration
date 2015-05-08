#! /usr/bin/python


#def convert_fits_to_image(fitsimage, force_stokes_I=True):
"""
Convert a fits image to a CASA image

For WSClean 1.7, use force_stokes_I = True to overide incorrect Stokes
keyword in model image headers. The restfreq is also set to avoid
problems with casapy2bbs.py
"""
import pyrap.images as pim
import pyrap.tables as pt
import numpy as np
import sys
 
fitsimage = sys.argv[1]
force_stokes_I = sys.argv[2]#True    

outfilename = fitsimage.split('.fits')[0] + '.image'
casaimage = pim.image(fitsimage)
casaimage.saveas(outfilename, overwrite=True)

if force_stokes_I:
    coords = casaimage.coordinates().dict()
    coords['stokes1']['stokes'] = ['I']
    freq = coords['spectral2']['wcs']['crval']
    coords['spectral2']['restfreqs'] = np.array([freq])
    outtable = pt.table(outfilename, readonly=False, ack=False)
    outtable.putkeywords({'coords': coords})
    outtable.done()
