#!/usr/bin/python

import os
import numpy
import pyrap.images
import pyfits
import sys
import lofar.bdsm as bdsm

def do_makecleanmask(image_name,threshpix,threshisl,atrousdo,ncores=8):
    mask_name  = image_name.split('.image')[0] + '.cleanmask'
    
    logging.info('makecleanmask: Making mask: '+mask_name)

    if atrousdo:
        threshisl = 4.0
        logging.info('Changing island threshold to %.1f because atrous_do=True' % threshisl)

    os.system('rm -rf ' + mask_name)

    # DO THE SOURCE DETECTION
    #img = bdsm.process_image(image_name, mean_map='zero', rms_box=(70,10), thresh_pix=numpy.float(o.threshpix), \
    #                         thresh_isl=numpy.float(o.threshisl), atrous_do=o.atrous_do,ini_method='curvature')
    
    img = bdsm.process_image(image_name, mean_map='zero',
                             rms_box=(80,20), thresh_pix=threshpix,
                             thresh_isl=threshisl,
                             atrous_do=atrousdo,ini_method='curvature',
                             adaptive_rms_box=True,
                             adaptive_thresh=150,
                             rms_box_bright=(35,7), rms_map=True, ncores=ncores)



    #img.show_fit()

    # WRITE THE ISLAND MASK
    img.export_image(img_type='island_mask',img_format='casa',outfile=mask_name, clobber=True)


if __name__ == '__main__':

    argc=len(sys.argv)
    from optparse import OptionParser
    parser = OptionParser(usage='%prog [options] <Input image>')
    #parser.add_option('--outfile', dest='file_clusters', help='Name of the output skymodel file [default catalogue.clusters.skymodel]', \
    #                      metavar='VAL', default='catalog.clusters.skymodel')
    parser.add_option('--threshpix', dest='threshpix', help='Pixel threshold [default: 5]', metavar='VAL', default=5)
    parser.add_option('--threshisl', dest='threshisl', help='Island threshold [default: 3]', metavar='VAL', default=3)
    parser.add_option('--atrous_do', dest='atrous_do', help='Decompose in wavelet scales threshold [default: False]', metavar='VAL', default='False')
    (o, args) = parser.parse_args()

        #if len(args) < 1: sys.exit("Missing BBS-format sky model.")


    # convert to boolean
    if o.atrous_do == 'False':
        o.atrous_do = False
    else:
        o.atrous_do = True

    image_name = args[0]

    do_makecleanmask(image_name,float(o.threshpix),float(o.threshisl),o.atrous_do)
    


