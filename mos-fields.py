#!/usr/bin/python

"""
A script to create a mosaic of LOFAR images
(The script was written to generate MSSS mosaics,
but it can be used in other cases as well)

This version hacked by Martin to use the masks and cluster images from
the wsclean facet code.

HISTORY
=======
v0.1	 G. Heald	       Created script
v0.2	 R. Breton	     Fixed avgpb behavior and weighting implementation
v0.3	 G. Heald	       Add weights, fix image naming, output sensitivity map
v0.4   S. van Velzen   Add beam & frequency info to header
v0.5   S. van Velzen   Use pyrap to save fits files
v0.6   J. Swinbank     Select a single Stokes parameter from input maps
v0.7	 G. Heald	       Fix RA behavior, add NCP option, pyfits tweak
v0.8	 G. Heald	       Fix behavior near RA=0
v0.9   G. Heald        Fix pyfits behavior (updating beam info)
v0.91  R. Breton       Can read the input image list as a file. Functional NCP flag (hardcoded to 20-degree angular size).
v0.92   G. Heald	Attempt to fix high dec overly-wide mosaics with a hardcoded adjustment
v0.93	G. Heald	Only load pylab if plotting is enabled
v0.94   G. Heald	Allow forcing maximum mosaic width (for high dec)

TO DO
=====

"""

version = '0.94 2014-11-17'

import pyrap.tables
import pyrap.images as pim
from pyrap import quanta
import numpy as np
import argparse
import pyfits
import os
import time
import glob
import re

def main(args):

    if args.plotimg:
        import pylab as plt
    # Generate lists of input images and check that they exist
    images=[]
    facets=[]
    psf_fwhm = [] # resolution
    frequency = [] # frequency of images (should be equal?)

    basestring=args.basestring
    imlist=glob.glob(basestring+'*.image')

    images=[i for i in imlist if not('nm' in i)]

    #construct image, facet number list
    images=[]
    fields=[]
    fnumbers=[]
    p=re.compile('imfield(\d)_cluster(.*)\.')
    for i in imlist:
        if 'nm' in i:
            continue
        m=p.match(i)
        if m is None:
            print 'failed to match',i
        assert(m is not None)
        images.append(i)
        fields.append(m.group(1))
        fnumbers.append(m.group(2))
    fnumberset=set(fnumbers)
    for f in fnumberset:
        fieldlist=[]
        for i,(field,facet) in enumerate(zip(fields,fnumbers)):
            if f==facet:
                fieldlist.append((field,i))
        while len(fieldlist)>1:
            # more than one field for the same facet...
            delfield,i=min(fieldlist)
            print 'de-duplicating',images[i]
            del(images[i])
            del(fields[i])
            del(fnumbers[i])
            del(fieldlist[fieldlist.index((delfield,i))])
    # now we have a non-redundant list
    for i in range(len(images)):
        print i,images[i],fields[i],fnumbers[i]

    # get the facet mask
    for fn in fnumbers:
        facets.append('templatemask_'+fn+'.masktmp')
        if not os.path.exists(facets[-1]):
            print "Error: facet image",facets[-1],"does not exist"
            return 1

    formstr = '{0:45s}  {1:45s} {2:s}  {3:s} {4:s} {5:s}'
    print formstr.format("-----","--------","------------","-------","-------","------")
    print formstr.format("Image", "FC image","Norm. weight", "Maj(ac)", "Min(ac)","PA(deg)")
    print formstr.format("-----","--------","------------","-------","-------","------")

    for i in range(len(images)):
        this_pim = pim.image(images[i])
        info_dict = this_pim.info()['imageinfo']['restoringbeam']
        # get beam info
        bpar_ma = quanta.quantity(info_dict['major']).get_value('deg')
        bpar_mi = quanta.quantity(info_dict['minor']).get_value('deg')
        bpar_pa = quanta.quantity(info_dict['positionangle']).get_value('deg')
        psf_fwhm.append([bpar_ma, bpar_mi, bpar_pa])
        frequency.append(this_pim.info()['coordinates']['spectral2']['restfreq'])
        print '{0:45.45s}  {1:45.45s} {2:0.2f}          {3:0.2f}    {4:0.2f}    {5:0.2f}'.format(images[i], facets[i], 0, bpar_ma*60, bpar_mi*60,bpar_pa)

    psf_fwhm = np.array(psf_fwhm)
    frequency = np.array(frequency)
    mean_psf_fwhm = np.mean(psf_fwhm, axis=0)
    mean_frequency = np.mean(frequency)
    print '\nmean Beam: {0:0.3f} maj (arcmin), {1:2.3f} min (arcmin), {2:0.2f} pa (deg)'.format(mean_psf_fwhm[0]*60, mean_psf_fwhm[1]*60, mean_psf_fwhm[2])
    print '(Frequency (MHz):', mean_frequency*1e-6

    if np.max(mean_frequency-frequency)/mean_frequency > 1e-6:
        print '\n\nWARNING.\nAre you using  images from different bands?'
        print 'Frequencies (Hz):', frequency
        time.sleep(2) # give user time to see this ...

    # Initialize some vectors
    declims = [] # store the limits of the declination axes
    #ralims = [] # store the limits of the r.a. axes
    raleft = []
    raright = []
    rainc = [] # store the r.a. increments in case they differ
    decinc = [] # store the dec increments in case they differ
    pims = [] # stores the pyrap images of the data
    pfcs = [] # stores the pyrap images of the facet images


# Get image frames for input images
    for im, fa in zip(images, facets):
        image = pim.image(im)
        sptcoords = image.coordinates().get_coordinate('spectral')
        nc = sptcoords.get_axis_size()
#        assert(sptcoords.get_image_axis() == 0)

        # Get Stokes axis. Ensure we are working with the Stokes parameter requested.
        stkcoords = image.coordinates().get_coordinate('stokes')
#        assert(stkcoords.get_image_axis() == 1)
        if stkcoords.get_axis_size() == 1:
            assert(stkcoords.get_stokes()[0] == args.stokes)
        else:
            stks = stkcoords.get_stokes().index(args.stokes)
            image = image.subimage(blc=(0, stks), trc=(nc-1, stks), dropdegenerate=False)
        ns = 1

        dircoords = image.coordinates().get_coordinate('direction')
        nx = dircoords.get_axis_size(axis=1)
        ny = dircoords.get_axis_size(axis=0)
        inc = dircoords.get_increment()
        ref = dircoords.get_referencepixel()
        val = dircoords.get_referencevalue()
        # wsclean image header is weird
        if val[1]<0:
            val[1]+=2*np.pi
        ra_axis = (range(nx)-ref[1])*inc[1]+val[1]
        dec_axis = (range(ny)-ref[0])*inc[0]+val[0]
        rainc.append(inc[1])
        decinc.append(inc[0])
        declims.append(min(dec_axis))
        declims.append(max(dec_axis))
        mean_ra = np.mean(ra_axis)
        print im,mean_ra
        #ralims.append((min(ra_axis)-mean_ra)*np.cos(val[0])+mean_ra)
        #ralims.append((max(ra_axis)-mean_ra)*np.cos(val[0])+mean_ra)
        raleft.append((ra_axis[0]-mean_ra)*np.cos(val[0])+mean_ra)
        raright.append((ra_axis[-1]-mean_ra)*np.cos(val[0])+mean_ra)
        pims.append(image)
        pfcs.append(pim.image(fa))


    # Generate the mosaic coordinate frame
    if not args.NCP:
        print('Using the regular mosaic mode.')
        master_dec = np.arange(min(declims),max(declims),min(decinc))
        if max(raleft)-min(raright) > 5.*np.pi/3.: # crossed RA=0
            print "Warning: I think the mosaic crosses RA=0, treating the coordinates as such."
            ##ralims[ralims>np.pi] -= 2.*np.pi
            #for i in range(len(ralims)):
            #    if ralims[i]>np.pi: ralims[i] = ralims[i]-2.*np.pi
            for i in range(len(raright)):
                raright[i] = raright[i]-2.*np.pi
        master_ra = np.arange(max(raleft),min(raright),max(rainc))
        lmra = len(master_ra)
        if args.maxwidth != 0:
            if lmra > args.maxwidth:
                xboundary = (lmra-args.maxwidth)/2
                master_ra = master_ra[xboundary:-xboundary]
        if args.verbose:
            print "Found ra,dec pixel increments (arcsec):"
            print np.array(rainc)*206265.,np.array(decinc)*206265.
        ma = pims[-1].coordinates()
        ma['direction'].set_referencepixel([len(master_dec)/2,len(master_ra)/2])
        ma['direction'].set_increment([decinc[np.argmin(np.abs(decinc))],rainc[np.argmin(np.abs(rainc))]])
        ma['direction'].set_referencevalue([master_dec[len(master_dec)/2],master_ra[len(master_ra)/2]])
    else:
        print('Using the special NCP mosaic mode.')
        ra_width = 20. / 180*np.pi
        dec_width = 20. / 180*np.pi
        rainc = rainc[np.argmin(np.abs(rainc))]
        decinc = decinc[np.argmin(np.abs(decinc))]
        ra_imsize = int(ra_width/np.abs(rainc))
        dec_imsize = int(dec_width/np.abs(decinc))
        master_ra = np.arange(ra_imsize, dtype=float)/ra_imsize*rainc-ra_width/2
        master_dec = np.arange(dec_imsize, dtype=float)/dec_imsize*decinc-dec_width/2
        ma = pims[-1].coordinates()
        ma['direction'].set_referencevalue([np.pi/2,0.])
        ma['direction'].set_increment([decinc,rainc])
        ma['direction'].set_referencepixel([dec_imsize/2.,ra_imsize/2.])

    # Initialize the arrays for the output image, sensitivity, and weights
    print 'making output image of size',len(master_dec),'x',len(master_ra)
    master_im = np.zeros((len(master_dec),len(master_ra)))
    master_mask = np.zeros((len(master_dec),len(master_ra)))

    # Reproject the images onto the master grid, weight and normalize
    for i in range(len(pims)):
        print 'doing image',i
        im = pims[i].regrid([2,3],ma,outshape=(nc,ns,len(master_dec),len(master_ra)))
        fa = pfcs[i].regrid([2,3],ma,outshape=(nc,ns,len(master_dec),len(master_ra)))
        imdata = np.squeeze(im.getdata())
        facmask = np.squeeze(fa.getdata())
        newim = imdata*facmask
#        newpb = pbdata
#        newwt = (weights[i]*newpb)**2
        master_im += newim
        master_mask += facmask
#        master_sens += newpb*newwt
#        master_weight += newwt

    print 'Blanking'
    blank=np.ones_like(im)*np.nan
    master_im=np.where(master_mask,master_im,blank)
    # Show image if requested
    if args.plotimg:
        plt.imshow(master_im,vmin=0.,vmax=0.5)
        plt.show()

    # Write fits files
    arrax = np.zeros( (1,1, len(master_im[:,0]), len(master_im[0,:])) )
    arrax[0,0,:,:] = master_im


    # Open new casa image for mosaic
    new_pim = pim.image('',shape=(1,1, len(master_dec),len(master_ra)), coordsys=ma)
    new_pim.putdata(arrax)
    # Write fits
    new_pim.tofits(args.outfits, overwrite=True)

    # need to add new beam info (not sure if this is possible with pyrap)
    hdu = pyfits.open(args.outfits,mode='update')
    header = hdu[0].header
    header.update('BMAJ',mean_psf_fwhm[0])
    header.update('BMIN',mean_psf_fwhm[1])
    header.update('BPA',mean_psf_fwhm[2])
    header.update('BUNIT',pims[-1].info()['unit'])
    header.update('RESTFRQ',mean_frequency)
    header.update('RESTFREQ',mean_frequency)
    newhdu = pyfits.PrimaryHDU(data=hdu[0].data, header=header)
    newhdu.writeto(args.outfits,clobber=True)

    return

print "LOFAR mosaic generator, v"+version+'\n'
parser = argparse.ArgumentParser(description="Mosaic MSSS images.")
parser.add_argument('-v','--verbose',help='Give some verbose output [default False]',action='store_true',default=False)
parser.add_argument('-o','--outfits',help='Output name of mosaic fits file [default mosaic.fits]',default='mosaic.fits')
parser.add_argument('-b','--basestring',help='Base string for image names, may include wild cards [default imfield]',default='imfield')
parser.add_argument('-p','--plotimg',help='Display image on screen? [default False]',action='store_true',default=False)
parser.add_argument('-S','--stokes',help='Stokes parameter to use?  [default I]',default='I')
parser.add_argument('-N','--NCP',help='Use NCP instead of SIN? This option does not work yet. [default False]',default=False,action='store_true')
parser.add_argument('-m','--maxwidth',help='Maximum number of pixels to consider for the width of the mosaic [default 0 = unlimited] This can be helpful at high declination.',default=0,type=int)
args = parser.parse_args()
main(args)
