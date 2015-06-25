#!/usr/bin/python

import os
import sys
import lofar.bdsm as bdsm
import logging

def do_makecleanmask_field_wsclean(image_name,threshpix,threshisl,atrousdo,ncores=8):
    mask_name  = image_name.split('-image')[0] + '.fitsmask'

    logging.info('makecleanmask_field_wsclean: Making mask: '+mask_name)

    os.system('rm -rf ' + mask_name)

    # DO THE SOURCE DETECTION
    img = bdsm.process_image(image_name, mean_map='zero',
                             rms_box=(300,60), thresh_pix=threshpix,
                             thresh_isl=threshisl,
                             atrous_do=atrousdo,
                             adaptive_rms_box=True,
                             adaptive_thresh=150,
                             rms_box_bright=(70,10),
                             atrous_jmax=3,ncores=ncores)

    #img.show_fit()

    # WRITE THE MASK FITS
    img.export_image(img_type='island_mask', img_format='fits', outfile=mask_name)


    # convert mask to casapy format

    #sys.exit()

    #img         = pyrap.images.image(mask_name)
    #pixels      = numpy.copy(img.getdata())
    #pixels_mask = 0.*numpy.copy(pixels)

    #hdulist   = pyfits.open(maskmodel + '.fits')
    #pixels_gs = hdulist[0].data

    #idx = numpy.where(pixels_gs > 0.0)
    #pixels_mask[idx] = 1.0

    #img.putdata(pixels_mask)


if __name__ == '__main__':
    #reproduce old command line behaviour
    
    argc=len(sys.argv)
    from optparse import OptionParser
    parser = OptionParser(usage='%prog [options] <Input image>')
    #parser.add_option('--outfile', dest='file_clusters', help='Name of the output skymodel file [default catalogue.clusters.skymodel]', \
    #                      metavar='VAL', default='catalog.clusters.skymodel')
    parser.add_option('--threshpix', dest='threshpix', help='Pixel threshold [default: 5]', metavar='VAL', default=5)
    parser.add_option('--threshisl', dest='threshisl', help='Island threshold [default: 3]', metavar='VAL', default=3)
    parser.add_option('--atrous_do', dest='atrous_do', help='Decompose in wavelet scales threshold [default: False]', metavar='VAL', default='False')
    parser.add_option('--casaregion',dest='casaregion',help='Casa region file [default: ""]', metavar='VAL', default='')
    (o, args) = parser.parse_args()


    # convert to boolean
    if o.atrous_do == 'False':
        o.atrous_do = False
    else:
        o.atrous_do = True
        o.threshisl = 3.0
        print 'Changing island threshold to 3 because atrous_do=True'


    image_name = args[0]

    do_makecleanmask_field_wsclean(image_name,float(o.threshpix),float(o.threshisl),o.atrous_do)
    
