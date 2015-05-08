#! /usr/bin/env python
import argparse
from argparse import RawTextHelpFormatter
from lofar import bdsm
import numpy
import sys


def main(image_name, mask_name, image_beam=False, atrous_do=False, threshisl=0.0, threshpix=0.0, rmsbox=None,
         iterate_threshold=False, adaptive_rmsbox=False, beam=None, img_format=None):
    convbox = rmsbox.rstrip(')').lstrip('(').split(',')
    rmsbox = (int(convbox[0]), int(convbox[1]))
    if image_beam:
        print 'do stuff for the beam'

    if atrous_do:
        threshisl = 4.0

    if iterate_threshold:
        # Start with high threshold and lower it until we get at least one island
        threshpix_orig = threshpix
        threshisl_orig = threshisl
        nisl = 0
        threshpix = 25
        threshisl = 15
        while nisl == 0:
            img = bdsm.process_image(image_name, mean_map='zero', rms_box=rmsbox,
                                     thresh_pix=numpy.float(threshpix), thresh_isl=numpy.float(threshisl),
                                     atrous_do=atrous_do, ini_method='curvature', beam=beam,
                                     adaptive_rms_box=adaptive_rmsbox, adaptive_thresh=20, quiet=True)
            nisl = img.nisl
            threshpix /= 1.2
            threshisl /= 1.2
        threshpix = threshpix_orig
        threshisl = threshisl_orig
    else:
        img = bdsm.process_image(image_name, mean_map='zero', rms_box=rmsbox,
                                 thresh_pix=numpy.float(threshpix), thresh_isl=numpy.float(threshisl),
                                 atrous_do=atrous_do, ini_method='curvature', beam=beam,
                                 adaptive_rms_box=adaptive_rmsbox, adaptive_thresh=20, quiet=True)

    img.export_image(img_type='island_mask', mask_dilation=0, outfile=mask_name,
                     img_format=img_format, clobber=True)

    log_file = mask_name + '.log'
    with open(log_file, 'wb') as f:
        f.write('# 5-sigma clipped rms (Jy/beam): {0}'.format(5.0 * img.clipped_rms))


if __name__ == '__main__':
    descriptiontext = "Make a clean mask.\n"

    parser = argparse.ArgumentParser(description=descriptiontext, formatter_class=RawTextHelpFormatter)
    parser.add_argument('image_name', help='Image name')
    parser.add_argument('mask_name', help='Mask name')
    parser.add_argument('-l', '--image_beam', help='Image beam', type=bool, default=False)
    parser.add_argument('-a', '--atrous_do', help='', type=bool, default=False)
    parser.add_argument('-i', '--threshisl', help='', type=float, default=0.0)
    parser.add_argument('-p', '--threshpix', help='', type=float, default=0.0)
    parser.add_argument('-x', '--rmsbox', help='', default=None)
    parser.add_argument('-t', '--iterate_threshold', help='', type=bool, default=False)
    parser.add_argument('-o', '--adaptive_rmsbox', help='', type=bool, default=False)
    parser.add_argument('-b', '--beam', help='', default=None)
    parser.add_argument('-f', '--img_format', help='', default=None)

    args = parser.parse_args()
    # convert rmsbox to tupel for bdsm
    #convbox = args.rmsbox.rstrip(')').lstrip('(').split(',')
    #tupbox = (int(convbox[0]), int(convbox[1]))
    main(args.image_name, args.mask_name, image_beam=args.image_beam, atrous_do=args.atrous_do,
         threshisl=args.threshisl, threshpix=args.threshpix, rmsbox=args.rmsbox, iterate_threshold=args.iterate_threshold,
         adaptive_rmsbox=args.adaptive_rmsbox, beam=args.beam, img_format=args.img_format)
