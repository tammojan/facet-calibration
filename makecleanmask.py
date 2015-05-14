#!/usr/bin/python

if __name__ == '__main__':

    import os
    import numpy
    import pyrap.images
    import pyfits
    import sys
    import lofar.bdsm as bdsm

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
        o.threshisl = 4.0
        print 'Changing island threshold to 4 because atrous_do=True'

    gs_cut = 1.25e-3
    image_name = args[0]

    mask_name  = image_name.split('.image')[0] + '.cleanmask'
    
    print '\n\n\n'
    print 'Making mask:', mask_name
    print '\n\n\n'

    os.system('rm -rf ' + mask_name)

    # DO THE SOURCE DETECTION
    #img = bdsm.process_image(image_name, mean_map='zero', rms_box=(70,10), thresh_pix=numpy.float(o.threshpix), \
    #                         thresh_isl=numpy.float(o.threshisl), atrous_do=o.atrous_do,ini_method='curvature')
    
    img = bdsm.process_image(image_name, mean_map='zero', rms_box=(80,20), thresh_pix=numpy.float(o.threshpix), \
                             thresh_isl=numpy.float(o.threshisl), atrous_do=o.atrous_do,ini_method='curvature', \
                             adaptive_rms_box=True, adaptive_thresh=150, rms_box_bright=(35,7), rms_map=True)



    #img.show_fit()

    # WRITE THE ISLAND MASK
    img.export_image(img_type='island_mask',img_format='casa',outfile=mask_name, clobber=True)


