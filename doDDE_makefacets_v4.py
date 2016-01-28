#!/usr/bin/env python

import matplotlib
import numpy
import os
import sys
from scipy import interpolate
import time
from subprocess import Popen, PIPE, STDOUT
import pyrap.tables as pt
import pwd
import pyrap.images
import matplotlib.pyplot as pl
from scipy.spatial import Voronoi
import matplotlib.path as mplPath
import matplotlib.transforms as mplTrans
import lofar.parmdb
from numpy import pi
from facet_utilities import run, bg, angsep, dec_to_str, ra_to_str

# Use scipy Voronoi tessellation to generate a set of masks
# original code from Martin Hardcastle
# Finite polygons code taken from https://gist.github.com/pv/8036995

# v3
# do the voronoi tesselation in image coordinates
# add crop to bounding box after voronoi tesselation
# general tidy up of code - set options at beginning
# note: possible matplotlib version issues with mpl.Path

# v2
# - added option flag useShifted to switch on shifting positions or not

# check high-DR
# check old selfcalv15

# check point source model


# TO DO
# - problem with e1 mask (sources not in clean mask)
# - problem with s3, diffuse source not in clean mask causes neg. bowl --> FIXED
# - check for mask gaps --> DONE
# - fit for clock in selfcal --> DONE
# - clean component overlap check -> avoid by giving up patches
# - reset phases to zero in selfcal part --> DONE
# - gain normalization --> DONE
# - HIGH-DYNAMIC RANGE (need to adjust merger/join parmdb, solution smoohting as well)
# - move m10 direction --> DONE
# - make normal template names files --> DONE
# - remove source list loading to prevent segfaults --> DONE

# - run the addback sources in seperate python script to avoid frequent Bus Errors DONE
# - run the addback sources FIELD in seperate python script to avoid frequent Bus Errors DONE
# - make sure parrallel FFT does not run out of memory (limit max. processes)/facets? DONE
# - automatic w-planes computation in imaging
# - how to deal with outlier sources low-freqs? # peel them first using low-freq bands only?
# - combine images into big image --> DONE
# - make fake avg pb images from mask to work with msss mosaic --> DONE
# - convolve images all to the same resolution (use casapy for that, easy to do)
# - What about second imaging and calibration cycle?
#    - 1. ADD back skymodel using "master solutions", then CORRECT using "master solutions"
#    - 2. image that again using the same setting
#    - 3. redo the subtract (will be slightly better....but solutions remained the same, just better noise) or just proceed to the next field?
# - how to deal with Toothbrush, how to subtract it from the data properly (make low-res images for that?)
# - memory control in full resolution subtract step --> DONE


def run_parallel(cmds, logs='none', Monitor=False, monitorname="run", nthreads=2):
    '''run a list of shell commands in run_parallel
optional - specify max number of processes (default 2)
           secify file to log output to for each command
    '''

    processes = set()
    max_processes = nthreads
    process_monitors = set()


    for icmd, cmd in enumerate(cmds):
        if logs not in ['auto', 'none']:
            cmdlog = logs[icmd]
        elif logs == 'auto':
            cmdlog = "log{icmd:d}".format(icmd=icmd)
        else:
            cmdlog = ''
        if logs != 'none':
            with open(cmdlog,'w') as f:
                print 'exec {icmd:d}: {cmd} > {log}'.format(cmd=cmd, icmd=icmd, log=cmdlog)
                p = Popen(cmd.split(), stdout=f, stderr=STDOUT)
                processes.add(p)
                pid=p.pid

                if Monitor:
                    monitor_cmd = 'python {scriptpath}/monitorjob.py {pid} {monitorname}'.format(pid=pid, monitorname=monitorname, scriptpath=scriptpath)
                    process_monitors.add(Popen(monitor_cmd.split()))
        else:
            print 'exec {icmd:d}: {cmd} > {log}'.format(cmd=cmd, icmd=icmd, log=cmdlog)
            p = Popen(cmd.split())
            processes.add(p)
            pid=p.pid

            if Monitor:
                monitor_cmd = 'python {scriptpath}/monitorjob.py {pid} {monitorname}'.format(pid=pid, monitorname=monitorname, scriptpath=scriptpath)
                process_monitors.add(Popen(monitor_cmd.split()))
        if len(processes) >= max_processes:
            os.wait()
            #if logs != 'none':
                #for ip, p in enumerate(processes):
                    #cmdlog = logs[ip]
                    #with open(cmdlog,'w') as f:
                        #ss = p.communicate()[0]
                        #try:
                            #f.write(ss)
                        #except:
                            #print ss
                            #print "log writing failed"
            processes.difference_update([p for p in processes if p.poll() is not None])
    for p in processes:
        if p.poll() is None:
            p.wait()
    # if we have submitted < max_processes, or on the last set #
    #os.wait()
    #if logs is not None:
        #for ip,p in enumerate(processes):
            #cmdlog = logs[ip]
            #with open(cmdlog,'w') as f:
                #f.write(p.communicate()[0])
    #processes.difference_update([p for p in processes if p.poll() is not None])

    return


def image_centre(image):
    '''
    return image centre
    '''
    im = pyrap.images.image(image)
    ic = im.coordinates()
    lon,lat = ic.get_referencevalue()[2]

    ra = lat*180./pi
    dec = lon*180./pi

    return ra,dec


def image_world_to_image(image,ra,dec):
    '''
    convert ra,dec arrays to x,y arrays for wcs in image
    '''
    im = pyrap.images.image(image)
    ic = im.coordinates()
    freq,stokes,pos = ic.get_referencevalue()


    x = []
    y = []
    for rai,deci in zip(ra,dec):
        f,c,yi,xi = im.topixel([freq,stokes[0],deci*pi/180,rai*pi/180.])
        x.append(xi)
        y.append(yi)
    x = numpy.array(x)
    y = numpy.array(y)

    return x,y


def image_image_to_world(image,x,y):
    '''
    convert x,y arrays to ra,dec arrays for wcs in image
    '''
    im = pyrap.images.image(image)
    ic = im.coordinates()
    freq,stokes,pos = ic.get_referencevalue()
    nf,nc,nx,ny = im.shape()

    ra = []
    dec = []
    for xi,yi in zip(x,y):
        f,c,deci,rai = im.toworld([0,0,yi,xi])
        ra.append(rai*180./pi)
        dec.append(deci*180./pi)
    ra = numpy.array(ra)
    dec = numpy.array(dec)

    return ra,dec


def ra_to_degrees(ra_str, delim=' '):
    '''
    converts array of strings or single string ra values to decimal degrees
    '''
    # string or array
    if isinstance(ra_str, str):
        if delim == 'h':
            ra_str = ra_str.replace('h',' ')
            ra_str = ra_str.replace('m',' ')
            t = ra_str.split()
        else:
            t = ra_str.split(delim)
        ra_deg = 15.*(float(t[0]) + float(t[1])/60. + float(t[2])/3600.)
        return ra_deg
    else:
        ra_deg = numpy.zeros(len(ra_str))
        for i,ra_s in enumerate(ra_str):
            if delim == 'h':
                ra_s = ra_s.replace('h',' ')
                ra_s = ra_s.replace('m',' ')
                t = ra_s.split()
            else:
                t = ra_s.split(delim)
            ra_deg[i] = 15.*(float(t[0]) + float(t[1])/60. + float(t[2])/3600.)
        return ra_deg


def dec_to_degrees(dec_str, delim=' '):
    '''
    converts array of strings or single string dec values to decimal degrees
    '''
    # string or array
    if isinstance(dec_str,str):
        if delim == 'd':
            dec_str = dec_str.replace('d',' ')
            dec_str = dec_str.replace('m',' ')
            t = dec_str.split()
        else:
            t = dec_str.split(delim)
        if '-' in dec_str:
        	dec_deg = (float(t[0]) - float(t[1])/60. - float(t[2])/3600.)
	else:
	        dec_deg = (float(t[0]) + float(t[1])/60. + float(t[2])/3600.)
        return dec_deg
    else:
        dec_deg = numpy.zeros(len(dec_str))
        for i,dec_s in enumerate(dec_str):
            if delim == 'd':
                dec_s = dec_s.replace('d',' ')
                dec_s = dec_s.replace('m',' ')
                t = dec_s.split()
            else:
                t = dec_s.split(delim)
        if '-' in dec_str:
        	dec_deg = (float(t[0]) - float(t[1])/60. - float(t[2])/3600.)
	else:
	        dec_deg = (float(t[0]) + float(t[1])/60. + float(t[2])/3600.)
        return dec_deg


def dir2pos(direction):

    def d2p(d):
        sra, sdec = d.split(',')
        rah = int(sra.split('h')[0])
        ram = int(sra.split('h')[1].split('m')[0])
        ras = float(sra.split('h')[1].split('m')[1])
        ra = 15.*(rah+ram/60.+ras/3600.)
        decd = int(sdec.split('d')[0])
        decm = int(sdec.split('d')[1].split('m')[0])
        decs = float(sdec.split('d')[1].split('m')[1])
 	if '-' in sdec:
        	dec = decd-decm/60.-decs/3600.
        else:
	        dec = decd+decm/60.+decs/3600.
	
        return ra,dec

    if isinstance(direction, numpy.ndarray) or  isinstance(direction, list):
        ra = numpy.zeros(len(direction))
        dec = numpy.zeros(len(direction))
        for i,d in enumerate(direction):
            ra[i], dec[i] = d2p(d)
    else:
        ra, dec = d2p(direction)
    return ra,dec


def create_dummyms_parset_formasks(msin, msout, numchanperms=20):
    ndppp_parset = (msin.split('.')[0]) +'_ndppp_dummy.parset'
    os.system('rm -f ' + ndppp_parset)

    f=open(ndppp_parset, 'w')
    f.write('msin ="%s"\n' % msin)
    f.write('msin.datacolumn = DATA\n')
    f.write('msin.autoweight = false\n')
    f.write('msout ="%s"\n' % msout)
    f.write('msout.writefullresflag=False\n')
    f.write('steps = [avg1]\n')
    f.write('avg1.type = squash\n')
    f.write('avg1.freqstep = %i\n' % numchanperms)
    f.write('avg1.timestep = 20\n')
    f.close()
    return ndppp_parset


def create_phaseshift_parset_formasks(msin, msout, source, direction):
    ndppp_parset = (msin.split('.')[0]) +'_ndppp_avgphaseshift.'+source+'.parset'
    os.system('rm -f ' + ndppp_parset)

    f=open(ndppp_parset, 'w')
    f.write('msin ="%s"\n' % msin)
    f.write('msin.datacolumn = DATA\n')
    f.write('msin.autoweight = false\n')
    f.write('msout ="%s"\n' % msout)
    f.write('msout.writefullresflag=False\n')
    f.write('steps = [shift]\n')
    f.write('shift.type        = phaseshift\n')
    f.write('shift.phasecenter = [%s]\n' % direction)
    f.close()
    return ndppp_parset


def voronoi_finite_polygons_2d(vor, radius=None):
    """
    Reconstruct infinite voronoi regions in a 2D diagram to finite
    regions.

    Parameters
    ----------
    vor : Voronoi
        Input diagram
    radius : float, optional
        Distance to 'points at infinity'.

    Returns
    -------
    regions : list of tuples
        Indices of vertices in each revised Voronoi regions.
    vertices : list of tuples
        Coordinates for revised Voronoi vertices. Same as coordinates
        of input vertices, with 'points at infinity' appended to the
        end.

    """


    # Use scipy Voronoi tessellation to generate a set of masks
    # code from Martin Hardcastle
    # Finite polygons code taken from https://gist.github.com/pv/8036995


    if vor.points.shape[1] != 2:
        raise ValueError("Requires 2D input")

    new_regions = []
    new_vertices = vor.vertices.tolist()

    center = vor.points.mean(axis=0)
    if radius is None:
        radius = vor.points.ptp().max()

    # Construct a map containing all ridges for a given point
    all_ridges = {}
    for (p1, p2), (v1, v2) in zip(vor.ridge_points, vor.ridge_vertices):
        all_ridges.setdefault(p1, []).append((p2, v1, v2))
        all_ridges.setdefault(p2, []).append((p1, v1, v2))

    # Reconstruct infinite regions
    for p1, region in enumerate(vor.point_region):
        vertices = vor.regions[region]

        if all(v >= 0 for v in vertices):
            # finite region
            new_regions.append(vertices)
            continue

        # reconstruct a non-finite region
        ridges = all_ridges[p1]
        new_region = [v for v in vertices if v >= 0]

        for p2, v1, v2 in ridges:
            if v2 < 0:
                v1, v2 = v2, v1
            if v1 >= 0:
                # finite ridge: already in the region
                continue

            # Compute the missing endpoint of an infinite ridge

            t = vor.points[p2] - vor.points[p1] # tangent
            t /= numpy.linalg.norm(t)
            n = numpy.array([-t[1], t[0]])  # normal

            midpoint = vor.points[[p1, p2]].mean(axis=0)
            direction = numpy.sign(numpy.dot(midpoint - center, n)) * n
            far_point = vor.vertices[v2] + direction * radius

            new_region.append(len(new_vertices))
            new_vertices.append(far_point.tolist())

        # sort region counterclockwise
        vs = numpy.asarray([new_vertices[v] for v in new_region])
        c = vs.mean(axis=0)
        angles = numpy.arctan2(vs[:,1] - c[1], vs[:,0] - c[0])
        new_region = numpy.array(new_region)[numpy.argsort(angles)]

        # finish
        new_regions.append(new_region.tolist())

    return new_regions, numpy.asarray(new_vertices)


def voronoi_finite_polygons_2d_box(vor, box):
    """
    Reconstruct infinite voronoi regions in a 2D diagram to finite
    box.

    Parameters
    ----------
    vor : Voronoi
        Input diagram
    box : (2,2) float array
        corners of bounding box
        numpy.array([[x1,y1],[x2,y2]])

    Returns
    -------
    poly : array of M (N,2) arrays
        polygon coordinates for M revised Voronoi regions.

    """

    if vor.points.shape[1] != 2:
        raise ValueError("Requires 2D input")
    if box.shape != (2,2):
        raise ValueError("Bounding box should be 2x2 array ((x1,y1),(x2,y2))")

    radius = numpy.max(box)

    # define the bounding box transform from the box extent - to be used to intersect with the regions
    bbox = mplTrans.Bbox(box)

    new_regions = []
    new_vertices = vor.vertices.tolist()

    center = vor.points.mean(axis=0)
    if radius is None:
        radius = vor.points.ptp().max()

    # Construct a map containing all ridges for a given point
    all_ridges = {}
    for (p1, p2), (v1, v2) in zip(vor.ridge_points, vor.ridge_vertices):
        all_ridges.setdefault(p1, []).append((p2, v1, v2))
        all_ridges.setdefault(p2, []).append((p1, v1, v2))

    # Reconstruct infinite regions
    for p1, region in enumerate(vor.point_region):
        vertices = vor.regions[region]

        if all(v >= 0 for v in vertices):
            # finite region
            new_regions.append(vertices)
            continue

        # reconstruct a non-finite region
        ridges = all_ridges[p1]
        new_region = [v for v in vertices if v >= 0]

        for p2, v1, v2 in ridges:
            if v2 < 0:
                v1, v2 = v2, v1
            if v1 >= 0:
                # finite ridge: already in the region
                continue

            # Compute the missing endpoint of an infinite ridge

            t = vor.points[p2] - vor.points[p1] # tangent
            t /= numpy.linalg.norm(t)
            n = numpy.array([-t[1], t[0]])  # normal

            midpoint = vor.points[[p1, p2]].mean(axis=0)
            direction = numpy.sign(numpy.dot(midpoint - center, n)) * n
            far_point = vor.vertices[v2] + direction * radius

            new_region.append(len(new_vertices))
            new_vertices.append(far_point.tolist())

        # sort region counterclockwise
        vs = numpy.asarray([new_vertices[v] for v in new_region])
        c = vs.mean(axis=0)
        angles = numpy.arctan2(vs[:,1] - c[1], vs[:,0] - c[0])
        new_region = numpy.array(new_region)[numpy.argsort(angles)]

        # finish
        new_regions.append(new_region.tolist())

    regions, imvertices = new_regions, numpy.asarray(new_vertices)
    #return new_regions, numpy.asarray(new_vertices)

    ## now force them to be in the bounding box
    poly = numpy.asarray([imvertices[v] for v in regions])

    newpoly = []

    for p in poly:
        polyPath = mplPath.Path(p)
        newpolyPath = polyPath.clip_to_bbox(bbox)
        pp = newpolyPath.vertices.transpose()
        newpoly.append(pp.transpose())

    return numpy.asarray(newpoly)


def make_facet_mask(imname, maskout, region, pad=True, edge=25):

    os.system('rm -rf ' + maskout)
    os.system('cp -r  ' + imname + ' ' + maskout)

    img    = pyrap.images.image(maskout)
    pixels = numpy.copy(img.getdata())
    sh     = numpy.shape(pixels)
    print "making facet mask " + maskout +" ("+str(sh[2])+")"

    facetmask    = (numpy.copy(pixels)*0.0) + 1.0

    region_pixel = []
    for p in region:
        p = p*pi/180.
        region_pixel.append(img.topixel([1,1,p[1], p[0]])[2:4])
    region_pixel = numpy.array(region_pixel)

    from PIL import Image, ImageDraw
    imgmask = Image.new('L', (sh[2], sh[3]), 0)
    lpoly=[tuple(v) for v in region_pixel]
    ImageDraw.Draw(imgmask).polygon(lpoly, outline=1, fill=1)
    mask = numpy.array(imgmask)
    mask = mask.transpose()  ## get the orientation right

    #pl.imshow(pixels[0,0,:,:])
    c = pl.imshow(mask)
    pl.colorbar()
    for p in region_pixel:
        pl.plot(p[0],p[1],'ko')

    if pad:
        # mask the edges of the image
        mask[0:edge,0:sh[3]] = 0.
        mask[0:sh[2],0:edge] = 0.

        mask[sh[2]-edge:sh[2],0:sh[3]] = 0.
        mask[0:sh[2],sh[3]-edge:sh[3]] = 0.
    facetmask[0,0,:,:] = mask

    img.putdata(facetmask)
    return


def add_facet_mask(maskout, region, value, direction, size,  lowres=15., actualres=1.5):
    print "adding facet to " + maskout +" s"+str(value)+ " ("+str(size)+")"

    #os.system('rm -rf ' + maskout)
    #os.system('cp -r  ' + imname + ' ' + maskout)

    img    = pyrap.images.image(maskout)
    pixels = numpy.copy(img.getdata())
    sh     = numpy.shape(pixels)

    facetmask    = numpy.copy(pixels)
    imagefacetmask = facetmask[0,0,:,:]


    ic = img.coordinates()
    freq,stokes,pos = ic.get_referencevalue()

    region_pixel = []
    for p in region:
        #print p
        p = p*pi/180.
        #print img.topixel([1,1,p[1], p[0]])[2:4]
        # handle negatives #
        #if p[0] >= 0 and p[1] >= 0:
        #print p[1], p[0]
        region_pixel.append(img.topixel([freq,stokes[0],p[1], p[0]])[2:4])
    region_pixel = numpy.array(region_pixel)

    from PIL import Image, ImageDraw
    imgmaskreg = Image.new('L', (sh[2], sh[3]), 0)
    lpoly=[tuple(v) for v in region_pixel]
    ImageDraw.Draw(imgmaskreg).polygon(lpoly, outline=1, fill=1)
    maskreg = numpy.array(imgmaskreg)
    maskreg = maskreg.transpose()  ## get the orientation right


    C = direction.split(',')
    cen_ra = ra_to_degrees(C[0],delim='h')
    cen_dec = dec_to_degrees(C[1],delim='d')
    cen_y, cen_x = img.topixel([1,1,cen_dec*pi/180., cen_ra*pi/180.])[2:4]

    D = int(size*actualres/lowres)/2
    #print D

    ## do the facet image mask
    imgmaskbox = Image.new('L', (sh[2], sh[3]), 0)
    facet_corners = numpy.array([[cen_y-D, cen_x-D],[cen_y-D, cen_x+D],[cen_y+D, cen_x+D], [cen_y+D, cen_x-D]])

    #print region_pixel
    #print facet_corners
    lpoly=[tuple(v) for v in facet_corners]
    ImageDraw.Draw(imgmaskbox).polygon(lpoly, outline=1, fill=1)
    maskbox = numpy.array(imgmaskbox)
    maskbox = maskbox.transpose()


    mask = value*maskreg*maskbox

    #import pylab as pl
    ##pl.imshow(pixels[0,0,:,:])
    #pl.imshow(mask)
    #for p in region_pixel:
        #pl.plot(p[0],p[1],'ko')

    #print mask
    #print numpy.sum(mask>0)
    #print imagefacetmask
    imagefacetmaskmask = (imagefacetmask == 0)
    thismask = (mask*imagefacetmaskmask)>0
    #thismask = numpy.where((mask*imagefacetmaskmask)>0)
    #pl.figure()
    #pl.imshow(mask)
    #pl.figure()
    #pl.imshow(imagefacetmaskmask)
    #pl.figure()
    #pl.imshow(mask*imagefacetmaskmask)
    #pl.show()
    #print thismask
    addmask = mask
    addmask[~thismask] *= 0

    #new_facetmask = facetmask == 0
    #addmask = numpy.where(facetmask[0,0,:,:]==0)
    #mask
    facetmask[0,0,:,:] += addmask
    #facetmask[0,0,:,:] += mask[thismask]

    img.putdata(facetmask)
    return


def show_facets(facetmap, directions, directions2=None, maxsize=6400, lowres=15., actualres=1.5, r=[1.,2.]):

    img    = pyrap.images.image(facetmap)
    pixels = numpy.copy(img.getdata())
    sh     = numpy.shape(pixels)
    NX = sh[2]
    NY = sh[3]

    pixels    = numpy.copy(pixels)
    mask = pixels[0,0,:,:]
    mask[mask==0] *= numpy.nan

    import matplotlib.pyplot as plt
    plt.figure()
    plt.title(facetmap)
    c = plt.imshow(mask, origin='lower', interpolation='none' )#, extent=[0,0,NX,NY])

    x0 = NX/2
    y0 = NY/2
    theta = numpy.arange(0,2*pi,pi/1000)
    for ri in r:
        r1 = ri*3600./lowres
        xi = x0 + r1*numpy.cos(theta)
        yi = y0 + r1*numpy.sin(theta)
        plt.plot(xi,yi,'k')
    plt.colorbar()
    plt.xlim(0,NX)
    plt.ylim(0,NY)


    for direction in directions:
        C = direction.split(',')
        cen_ra = ra_to_degrees(C[0],delim='h')
        cen_dec = dec_to_degrees(C[1],delim='d')
        #print 'original direction:', direction
        #print 'original direction world:', cen_ra, cen_dec
        cen_y, cen_x = img.topixel([1,1,cen_dec*pi/180., cen_ra*pi/180.])[2:4]
        #print 'original direction pixel:', cen_x, cen_y
        if directions2 is None:
            value = mask[cen_y, cen_x]
            try:
                svalue = str(int(value))
            except:
                svalue ='?'
            plt.text(cen_x, cen_y, svalue)

        plt.plot(cen_x, cen_y, 'ko')

    if directions2 is not None:
        for direction in directions2:
            C = direction.split(',')
            cen_ra = ra_to_degrees(C[0],delim='h')
            cen_dec = dec_to_degrees(C[1],delim='d')
            #print 'original direction:', direction
            #print 'original direction world:', cen_ra, cen_dec
            cen_y, cen_x = img.topixel([1,1,cen_dec*pi/180., cen_ra*pi/180.])[2:4]
            #print 'original direction pixel:', cen_x, cen_y
            value = mask[cen_y, cen_x]
            try:
                svalue = str(int(value))
            except:
                svalue ='?'

            plt.plot(cen_x, cen_y, 'k+')
            plt.text(cen_x, cen_y, svalue)


    plt.savefig(facetmap+'.png')

    return


def find_newsize_centre_lowresfacets(facetmap, region, direction, source, maxsize=6400, lowres=15., actualres=1.5, edge=25):

    value = int(source.replace('s',''))

    print "doing source " + source

    img    = pyrap.images.image(facetmap)
    pixels = numpy.copy(img.getdata())
    sh     = numpy.shape(pixels)
    NX = sh[2]
    NY = sh[3]

    pixels    = numpy.copy(pixels)
    mask = pixels[0,0,:,:]


    C = direction.split(',')
    cen_ra = ra_to_degrees(C[0],delim='h')
    cen_dec = dec_to_degrees(C[1],delim='d')
    #print 'original direction:', direction
    #print 'original direction world:', cen_ra, cen_dec
    cen_y, cen_x = img.topixel([1,1,cen_dec*pi/180., cen_ra*pi/180.])[2:4]
    #print 'original direction pixel:', cen_x, cen_y

    # get the facet pixels
    maskpix = numpy.where(mask == value)
    maskpix_y = maskpix[0]
    maskpix_x = maskpix[1]
    ## get the nonfacet pixels
    #nonmaskpix = numpy.where(mask != value)
    #nonmaskpix_x = nonmaskpix[0]
    #nonmaskpix_y = nonmaskpix[1]

    #import matplotlib.pyplot as plt
    #plt.imshow(mask == value)
    #plt.show()

    #sep = numpy.sqrt((nonmaskpix_x-cen_x)**2.+ (nonmaskpix_y-cen_y)**2.)
    #maxedgesep = numpy.max(sep)
    #minedgesep = numpy.min(sep)
    #print maxedgesep, minedgesep, maxedgesep/minedgesep

    X1, X2 = numpy.min(maskpix_x),numpy.max(maskpix_x)
    Y1, Y2 = numpy.min(maskpix_y),numpy.max(maskpix_y)
    #print 'X1,X2 = ', X1,X2
    #print 'Y1,Y2 = ', Y1,Y2
    #X01 = abs(numpy.max(maskpix_x)-cen_x)
    #X02 = abs(cen_x-numpy.min(maskpix_x))
    #Y01 = abs(numpy.max(maskpix_y)-cen_y)
    #Y02 = abs(cen_y-numpy.min(maskpix_y))
    #new_cen_x = (max(maskpix_x) + min(maskpix_x))/2.
    #new_cen_y = (max(maskpix_y) + min(maskpix_y))/2.
    #max_span_x = 2.*max(X01, X02)
    #max_span_y = 2.*max(Y01, Y02)
    #print 'max x dist to facet edge: ', max_span_x
    #print 'max y dist to facet edge: ', max_span_y

    maxrealsize = maxsize  #*lowres/actualres  # in actual pixels
    #print 'maxrealsize: ',maxrealsize

    new_cen_x = (X2+X1)/2.
    new_cen_y = (Y2+Y1)/2.
    span_x = X2-X1
    span_y = Y2-Y1
    span = max(span_x,span_y)
    real_span = span*lowres/actualres  # in actual pixels
    real_span += 2*edge  # add real edge pixels
    #print 'span (real)', span, real_span
    #trysizes = numpy.copy(sorted([6400,6144,5600,5400,5184,5000,4800,4608,4320,4096,3840,3600,3200,3072,2880,2560,2304,2048, 1600, 1536, 1200, 1024, 800, 512]))
    # for testing at lowres
    trysizes = numpy.copy(sorted([6400,6144,5600,5400,5184,5000,4800,4608,4320,4096,3840,3600,3200,3072,2880,2560,2304,2048, 1600, 1536, 1200, 1024, 800, 512, 256, 128, 64]))
    trysizes = trysizes*1.5/actualres  # set proper sizes when res is 1.5 arcsec
    trysizes = trysizes[trysizes<=maxrealsize]
    trysizes1 = numpy.copy(trysizes[::-1])
    #print trysizes1
    #print trysizes1>= real_span
    selind = numpy.sum(trysizes1>= real_span)-1
    if selind < 0:
        new_real_span = trysizes1[0]
    else:
        new_real_span = trysizes1[selind]
    new_span = new_real_span*actualres/lowres
    #print 'new_span (real)', new_span, new_real_span

    new_cen_dec, new_cen_ra = img.toworld([0,0,new_cen_x, new_cen_y])[2:4]
    new_cen_ra *= 180./pi
    new_cen_dec *= 180./pi
    new_position =  '{x:s},{y:s}'.format(x=ra_to_str(new_cen_ra,delim='h'),y=dec_to_str(new_cen_dec,delim='d'))

    #d  = new_span/2
    #this_mask = mask[max(0,new_cen_x-d):min(new_cen_x+d,NX), max(new_cen_y-d,0):min(new_cen_y+d,NY)]


    ##maskpix = numpy.where(mask == value)
    ##maskpix_x = maskpix[0]
    ##maskpix_y = maskpix[1]
    ##new_maxsize = max(numpy.sqrt((maskpix_x-new_cen_x)**2.+ (maskpix_y-new_cen_y)**2.))
    ##print maxsize, new_maxsize

    #outer = False
    #special = False
    #if (X1 == 0) or (Y1 == 0) or (X2 == NX) or (Y2 == NY):
        #print "***edge in mask***"
        #outer = True
        ##raw_input('ok')
    #if numpy.sum(this_mask == 0) > 0:
        #print "***boundary in mask***"
        #outer = True
        ##raw_input('ok')
    #if (maxedgesep/minedgesep) > 4.:
        #print "***elongated facet***"
        #special = True
    #axratio = 1.25
    #if (span_x/span_y > axratio) or (span_y/span_x > axratio):
        #print "***elongated facet***"
        #special = True

    #if outer or special:
        #print "using calibrator direction"
        ##maskpix = numpy.where(mask == value)
        ##maskpix_x = edgepix[0]
        ##maskpix_y = edgepix[1]
        ##maxsize = max(numpy.sqrt((maskpix_x-cen_x)**2.+ (maskpix_y-cen_y)**2.))
        ##new_maxsize = max(numpy.sqrt((maskpix_x-new_cen_x)**2.+ (maskpix_y-new_cen_y)**2.))
        ##print maxsize, new_maxsize
        #new_position =  direction

    #import matplotlib.pyplot as plt
    #plt.imshow(this_mask)
    #plt.show()

    ## update dist
    #nX01 = abs(numpy.max(maskpix_x)-new_cen_x)
    #nX02 = abs(new_cen_x-numpy.min(maskpix_x))
    #nY01 = abs(numpy.max(maskpix_y)-new_cen_y)
    #nY02 = abs(new_cen_y-numpy.min(maskpix_y))
    #new_min_boundary_dist = min(numpy.sqrt((edgepix_x-new_cen_x)**2.+ (edgepix_y-new_cen_y)**2.))

    ## if any x,y dist from centre is
    #if numpy.sum(numpy.array([X01, X02, Y01, Y02]) >= min_boundary_dist) > 0:
        #print "boundary hit on original coords"
        #boundary = True
    #if numpy.sum(numpy.array([nX01, nX02, nY01, nY02]) >= new_min_boundary_dist) > 0:
        #print "boundary hit on updated coords"
        #boundary = True

    #if not boundary:
        #new_cen_x = (max(maskpix_x) + min(maskpix_x))/2.
        #new_cen_y = (max(maskpix_y) + min(maskpix_y))/2.
        #new_span = max((max(maskpix_x) - min(maskpix_x)),(max(maskpix_y) - min(maskpix_y)))
    #else:
        #print "** boundary **"
        #new_cen_x = cen_x
        #new_cen_y = cen_y
        #new_span = min_boundary_dist

    #print new_cen_x, new_cen_y, new_span

    ##def find_imsize_from_region(region, resolution):
    #'''
    #get approx image size from region and image resolution - ignoring correct wcs
    #resolution - in arcsec
    #'''
    #edge = 25.

    #t = region.transpose()
    #ra = t[0]
    #dec = t[1]

    #span_x = (ra.max() - ra.min())*3600./resolution
    #span_y = (dec.max() - dec.min())*3600./resolution

    #max_span = max(span_x, span_y) + 2*edge

    #midra = (ra.max() + ra.min())/2.
    #middec = (dec.max() + dec.min())/2.
    #position =  '{x:s},{y:s}'.format(x=ra_to_str(midra,delim='h'),y=dec_to_str(middec,delim='d'))

    #trysizes = numpy.copy(sorted([6400,6144,5600,5400,5184,5000,4800,4608,4320,4096,3840,3600,3200,3072,2880,2560,2304,2048, 1600, 1536, 1200, 1024, 800, 512]))
    #trysizes1 = numpy.copy(trysizes[::-1])
    #this_size = trysizes1[numpy.sum(trysizes1>= max_span)]

    #return this_size, position


    #facetmask[0,0,:,:] += mask

    #img.putdata(facetmask)
    #print new_position, new_real_span

    return new_position, new_real_span


def find_newsize_lowresfacets(facetmap, region, direction, source, maxsize=6400, lowres=15., actualres=1.5, debug=True, findedge=None, edge=25):

#if 1:

    value = int(source.replace('s',''))

    if debug : print "doing source " + source

    #os.system('rm -rf ' + maskout)
    #os.system('cp -r  ' + imname + ' ' + maskout)

    img    = pyrap.images.image(facetmap)
    pixels = numpy.copy(img.getdata())
    sh     = numpy.shape(pixels)
    NX = sh[2]
    NY = sh[3]

    pixels    = numpy.copy(pixels)
    mask = pixels[0,0,:,:]


    C = direction.split(',')
    cen_ra = ra_to_degrees(C[0],delim='h')
    cen_dec = dec_to_degrees(C[1],delim='d')
    #print 'original direction:', direction
    #print 'original direction world:', cen_ra, cen_dec
    cen_y, cen_x = img.topixel([1,1,cen_dec*pi/180., cen_ra*pi/180.])[2:4]
    #print 'original direction pixel:', cen_x, cen_y

    # get the facet pixels
    maskpix = numpy.where(mask == value)
    maskpix_y = maskpix[0]
    maskpix_x = maskpix[1]
    ## get the nonfacet pixels
    #nonmaskpix = numpy.where(mask != value)
    #nonmaskpix_x = nonmaskpix[0]
    #nonmaskpix_y = nonmaskpix[1]

    #import matplotlib.pyplot as plt
    #plt.imshow(mask == value)
    #plt.show()

    #sep = numpy.sqrt((nonmaskpix_x-cen_x)**2.+ (nonmaskpix_y-cen_y)**2.)
    #maxedgesep = numpy.max(sep)
    #minedgesep = numpy.min(sep)
    #print maxedgesep, minedgesep, maxedgesep/minedgesep

    X1, X2 = numpy.min(maskpix_x),numpy.max(maskpix_x)
    Y1, Y2 = numpy.min(maskpix_y),numpy.max(maskpix_y)
    if debug : print 'X1,X2 = ', X1,X2
    if debug : print 'Y1,Y2 = ', Y1,Y2
    X01 = abs(X1-cen_x)
    X02 = abs(cen_x-X2)
    Y01 = abs(Y1-cen_y)
    Y02 = abs(cen_y-Y2)
    #new_cen_x = (max(maskpix_x) + min(maskpix_x))/2.
    #new_cen_y = (max(maskpix_y) + min(maskpix_y))/2.
    max_span_x = 2.*max(X01, X02)
    max_span_y = 2.*max(Y01, Y02)
    if debug : print 'max x dist to facet edge: ', max_span_x
    if debug : print 'max y dist to facet edge: ', max_span_y

    maxrealsize = maxsize  #*lowres/actualres  # in actual pixels
    #print 'maxrealsize: ',maxrealsize

    #new_cen_x = (X2+X1)/2.
    #new_cen_y = (Y2+Y1)/2.
    #span_x = X2-X1
    #span_y = Y2-Y1
    span = max(max_span_x,max_span_y)
    real_span = span*lowres/actualres  # in actual pixels
    if findedge is not None:
        edge=findedge(real_span)
        maxedge=findedge(maxsize)
        edge=int(min([maxedge,edge]))
        if debug: print 'Calculating edge to be:',edge
    real_span += 2*edge  # add real edge pixels
    if debug : print 'span (real)', span, real_span
    #trysizes = numpy.copy(sorted([6400,6144,5600,5400,5184,5000,4800,4608,4320,4096,3840,3600,3200,3072,2880,2560,2304,2048, 1600, 1536, 1200, 1024, 800, 512]))
    # for testing at lowres
    trysizes = numpy.copy(sorted([6400,6144,5600,5400,5184,5000,4800,4608,4320,4096,3840,3600,3200,3072,2880,2560,2304,2048, 1600, 1536, 1200, 1024, 800, 512, 256, 128, 64]))
    trysizes = trysizes*1.5/actualres  # set proper sizes when res is 1.5 arcsec
    trysizes = trysizes[trysizes<=maxrealsize]
    trysizes1 = numpy.copy(trysizes[::-1])
    #print trysizes1
    #print trysizes1>= real_span
    selind = numpy.sum(trysizes1>= real_span)-1
    if selind < 0:
        new_real_span = trysizes1[0]
    else:
        new_real_span = trysizes1[selind]
    new_span = new_real_span*actualres/lowres
    #print 'new_span (real)', new_span, new_real_span

    #new_cen_dec, new_cen_ra = img.toworld([0,0,new_cen_x, new_cen_y])[2:4]
    #new_cen_ra *= 180./pi
    #new_cen_dec *= 180./pi
    #new_position =  '{x:s},{y:s}'.format(x=ra_to_str(new_cen_ra,delim='h'),y=dec_to_str(new_cen_dec,delim='d'))

    if debug : print 'image size, edge:',new_real_span, edge

    return new_real_span, edge


def write_casapy_region(regfile, polygon, name):
    #poly[[x1, y1], [x2, y2], [x3, y3], ...]
    s= '''#CRTFv0 CASA Region Text Format version 0
'''
    s += 'poly['
    dirs = ['[{ra},{dec}]'.format(ra=ra_to_str(p[0]), dec=dec_to_str(p[1], delim='.')) for p in polygon]
    s += ', '.join(dirs)
    s += ']\n'
    #s+= '{name}'.format(name=name)
    #s+= r'} '+'\n'
    with open(regfile,'w') as f:
        f.write(s)
    return


def write_ds9_region(regfile, polygon, name):
    s= '''global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1
fk5
'''
    s += 'polygon('
    dirs = ['{ra},{dec}'.format(ra=ra_to_str(p[0]), dec=dec_to_str(p[1])) for p in polygon]
    s += ', '.join(dirs)
    s += r')# text={'
    s+= '{name}'.format(name=name)
    s+= r'} '+'\n'
    with open(regfile,'w') as f:
        f.write(s)
    return


def write_ds9_allregions(regfile, polygons, names):
    s= '''global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1
fk5
'''

    with open(regfile,'w') as f:
        f.write(s)
        for polygon, name in zip(polygons,names):
            s = 'polygon('
            dirs = ['{ra},{dec}'.format(ra=ra_to_str(p[0]), dec=dec_to_str(p[1])) for p in polygon]
            s += ', '.join(dirs)
            s += r')# text={'
            s+= '{name}'.format(name=name)
            s+= r'} '+'\n'
            f.write(s)
    return


def fwhm_freq(freq, alpha=1.02):
    c = 299.79  ## m MHz
    #D = 41.05  # for remote stations
    D = 30.75  # for core stations - DUAL_INNER
    lam = c/freq
    fwhm = 57.2957795*alpha*lam/D
    return fwhm


### main code starts here ###

if __name__=='__main__':

    username = pwd.getpwuid(os.getuid())[0]

    if len(sys.argv)<2:
        raise Exception('Give the path to the setup code')

    print 'Using',sys.argv[1],'as the setup code'
    execfile(sys.argv[1])
    print 'script path is',SCRIPTPATH

    sys.path.append(SCRIPTPATH)
    from coordinates_mode import *
    
    try:
        edge_scale
    except NameError:
        edge_scale=None

    source_info_rec = numpy.genfromtxt(peelsourceinfo, \
                                       dtype="S50,S25,S5,S5,i8,i8,i8,i8,S2,S255,S255,S255,S5", \
                                       names=["sourcelist","directions","atrous_do","mscale_field","imsizes",\
                                       "cellsizetime_p","cellsizetime_a","fieldsize","dynamicrange",\
                                       "regionselfc","regionfield","peelskymodel","outliersource"],comments='#')

    sourcelist = source_info_rec["sourcelist"]
    directions = source_info_rec["directions"]
    atrous_do = source_info_rec["atrous_do"]
    mscale_field = source_info_rec["mscale_field"]
    imsizes = source_info_rec["imsizes"]
    cellsizetime_p = source_info_rec["cellsizetime_p"]
    cellsizetime_a = source_info_rec["cellsizetime_a"]
    fieldsize = source_info_rec["fieldsize"]
    dynamicrange = source_info_rec["dynamicrange"]
    regionselfc = source_info_rec["regionselfc"]
    regionfield = source_info_rec["regionfield"]
    peelskymodel = source_info_rec["peelskymodel"]
    outliersource = source_info_rec["outliersource"]
    print sourcelist,directions


    sresolution = '{r:.1f}arcsec'.format(r=resolution)
    slowresolution = '{r:.1f}arcsec'.format(r=lowresolution)


    fwhm = fwhm_freq(FREQ)
    rad = [fwhm/2., numpy.sqrt(numpy.log(3.)/numpy.log(2.))*fwhm/2.]  # PB= 0.5, 0.3
    print 'FWHM: {ff:.2f} deg (draw radii at {ss} deg)'.format(ff=fwhm, ss=', '.join(['{rad:.2f}'.format(rad=radi) for radi in rad]))

    # Get the number of channels in the MS
    freq_tab     = pt.table(ms + '/SPECTRAL_WINDOW')
    numchanperms = freq_tab.getcol('NUM_CHAN')[0]
    logging.info('Number of channels per ms is {:d}'.format(numchanperms))
    freq_tab.close()

    #mslist = ['BOOTES24_SB190-199.2ch8s.ms']
    ### MAKE ALL THE MASKS AT ONCE ###



    # create low timeres, low freqres ms for faster shifting
    tmpms  = os.path.basename(ms).split('.')[0] + '.dummy.ms'
    parset = create_dummyms_parset_formasks(ms, tmpms, numchanperms=numchanperms)
    if not os.path.exists(tmpms):
        print "making dummy ms (freq/time averaged) set for image-making"
        os.system('NDPPP ' + parset)

    npixlow = int(maxsize_fieldimage*3600./lowresolution)
    npix = int(maxsize_fieldimage*3600./resolution)

    dummy_image = 'empty_wide.image'

    ## make low_res full image
    cmd = 'casapy --nogui --nologfile  -c '+SCRIPTPATH+'/make_empty_image.py '+ tmpms + ' ' + dummy_image + ' '  + str(npixlow) + ' ' +slowresolution
    if not os.path.exists(dummy_image):
        #os.system('rm -rf empty_wide.image')
        print "making image template: {n:d} pix at {s}/pix".format(n=npixlow, s=slowresolution)
        os.system(cmd)

    ra0,dec0 = image_centre(dummy_image)

    #tesselate#
    print "doing tesselation"
    ra,dec = dir2pos(directions)

    x,y = image_world_to_image(dummy_image,ra,dec)


    x1,y1 = 0.,0.
    x2,y2 = npixlow,npixlow
    box = numpy.array([[x1,y1],[x2,y2]])


    # use the image coordinates - as the non-linear sky is handled properly
    vor = Voronoi(numpy.array((x, y)).transpose())

    # convert the voronoi object into something useful (array of regions, truncated to imagesize)
    # impoly is an array of (n,2) arrays defining the polygon region
    impoly = voronoi_finite_polygons_2d_box(vor, box)

    # convert image coordinates of regions to
    worldpoly = []
    for p in impoly:
        pp = p.transpose()
        pra,pdec = image_image_to_world(dummy_image,pp[0],pp[1])
        worldpoly.append(numpy.array((pra,pdec)).transpose())
    worldpoly = numpy.asarray(worldpoly)


    plot_facets = True
    if plot_facets:
        pl.figure(figsize=(8,4))
        ax1 = pl.subplot(121)
        ax1.plot(x,y,'*')
        ax1.plot([x1,x1,x2,x2,x1],[y1,y2,y2,y1,y1])
        for p in impoly:
            pp = p.transpose()
            ax1.plot(pp[0],pp[1])


        ax1.set_xlabel('x')
        ax1.set_ylabel('y')


        ax2 = pl.subplot(122)
        ax2.plot(ra,dec,'*')
        for p in worldpoly:
            pp = p.transpose()
            ax2.plot(pp[0],pp[1])
        x1,x2 = ax2.get_xlim()
        ax2.set_xlabel('ra')
        ax2.set_ylabel('dec')
        pl.savefig('voronoi_facets.png')

    # use world coordinates now
    poly = worldpoly



    if 1:
        # initial tesselation with huge facets
        print "allow huge facets"
        facet_image_name_huge = 'facets_wide_huge.image'
        os.system('rm -rf {im2}'.format(im2=facet_image_name_huge))
        os.system('cp -r {im1} {im2}'.format(im1=dummy_image,im2=facet_image_name_huge))

        sizes = 2*npix*numpy.ones(len(sourcelist))

        for source_id,source in enumerate(sourcelist):
            value = int(source.replace('s',''))
            add_facet_mask(facet_image_name_huge, poly[source_id], value, directions[source_id], sizes[source_id], lowres=lowresolution, actualres=resolution)


        show_facets(facet_image_name_huge, directions, r=rad)



    if 1:
        # initial tesselation with huge facets
        print "limit facet sizes (inner and outer)"
        facet_image_name = 'facets_wide_maxsize.image'
        os.system('rm -rf {im2}'.format(im2=facet_image_name))
        os.system('cp -r {im1} {im2}'.format(im1=dummy_image,im2=facet_image_name))


        outlierdist = rad[0]  # make it at 0.5 power point
        # set max sizes
        sizes = []
        for source_id,source in enumerate(sourcelist):
            d = directions[source_id]

            C = d.split(',')
            ra = ra_to_degrees(C[0],delim='h')
            dec = dec_to_degrees(C[1],delim='d')
            dist_from_pointing = angsep(ra,dec,ra0,dec0)


            if dist_from_pointing > outlierdist:
                sizes.append(maxoutliersize)
            else:
                sizes.append(maxcentralsize)


        for source_id,source in enumerate(sourcelist):
            value = int(source.replace('s',''))
            add_facet_mask(facet_image_name, poly[source_id], value, directions[source_id], sizes[source_id], lowres=lowresolution, actualres=resolution)


        show_facets(facet_image_name, directions, r=rad)


    ## obsolete - never shift the facet centres
    ## get smallest size possible - with shifting the facet centre
    ## but use the max image size
    #newpositionsshift = []
    #newsizesshift = []
    #for source_id,source in enumerate(sourcelist):
        #new_pos, new_size = find_newsize_centre_lowresfacets(facet_image_name_huge, poly[source_id], directions[source_id], source, maxsize=max_fieldsize, lowres=lowresolution, actualres=resolution)
        #newpositionsshift.append(new_pos)
        #newsizesshift.append(new_size)


    #for source_id,source in enumerate(sourcelist):
        #print source, newsizesshift[source_id], newsizes[source_id], sizes[source_id]

    #facet_image_name = 'facets_wide_shifted.image'
    #os.system('rm -rf {im2}'.format(im2=facet_image_name))
    #os.system('cp -r {im1} {im2}'.format(im1=dummy_image,im2=facet_image_name))
    #for source_id,source in enumerate(sourcelist):
        #value = int(source.replace('s',''))
        #add_facet_mask(facet_image_name, poly[source_id], value, newpositionsshift[source_id], newsizesshift[source_id], lowres=lowresolution, actualres=resolution)
    #show_facets(facet_image_name, directions, directions2=newpositionsshift, r=rad)


    reduce_sizes = True
    if reduce_sizes:

        print "reducing image sizes where possible"

        # get smallest size possible - without shifting the facet centre
        newsizes = []
        edges = []
        for source_id,source in enumerate(sourcelist):
            new_size,new_edge = find_newsize_lowresfacets(facet_image_name_huge, poly[source_id], directions[source_id], source, maxsize=max_fieldsize, lowres=lowresolution, actualres=resolution, findedge=edge_scale)
            newsizes.append(new_size)
            edges.append(new_edge)

        facet_image_name = 'facets.image'
        os.system('rm -rf {im2}'.format(im2=facet_image_name))
        os.system('cp -r {im1} {im2}'.format(im1=dummy_image,im2=facet_image_name))
        for source_id,source in enumerate(sourcelist):
            value = int(source.replace('s',''))
	    print value, directions[source_id]
            add_facet_mask(facet_image_name, poly[source_id], value, directions[source_id], newsizes[source_id], lowres=lowresolution, actualres=resolution)
        show_facets(facet_image_name, directions, directions2=directions, r=rad)


    sizes = newsizes
    positions = directions



    for p,source in zip(poly,sourcelist):
        write_ds9_region('peel_facet_{source}.reg'.format(source=source), p, source)
        write_casapy_region('peel_facet_{source}.rgn'.format(source=source), p, source)
    write_ds9_allregions('peel_facets.reg', poly, sourcelist)

    #sys.exit()





    # make empty images - maxsize #
    NDPPPcmds = []
    Imcmds = []
    mslist = []
    imagelist = []
    sizelist = []
    for source_id,source in enumerate(sourcelist):

        tmpn  = ms.split('.')[0] + '.' + source + '.ms.tmp'
        mslist.append(tmpn)
        imagelist.append('templatemask_' + source)
        sizelist.append(str(sizes[source_id]))
        parset = create_phaseshift_parset_formasks(tmpms, tmpn, source, positions[source_id])
        NDPPPcmds.append('NDPPP ' + parset)

    run_parallel(NDPPPcmds, nthreads=8,logs='auto')

    mslistfile = 'mslist.npy'
    imagelistfile = 'magelist.npy'
    sizelistfile = 'sizelist.npy'
    numpy.save(mslistfile, numpy.array(mslist))
    numpy.save(imagelistfile, numpy.array(imagelist))
    numpy.save(sizelistfile, numpy.array(sizelist))
    Imcmds = ['casapy --nogui --nologfile  -c '+SCRIPTPATH+'/make_empty_image_many.py '+ mslistfile + ' ' + imagelistfile + ' ' + sizelistfile + ' ' +sresolution ]
    run_parallel(Imcmds, nthreads=1,logs='auto')

    ## clean up
    for source_id,source in enumerate(sourcelist):
        tmpn  = ms.split('.')[0] + '.' + source + '.ms.tmp'
        os.system('rm -rf '+tmpn)
        tmpn  = ms.split('.')[0] +'_ndppp_avgphaseshift.'+source+'.parset'
        os.system('rm -rf '+tmpn)

    #mask_cmds = ["casapy --nogui --nologfile -c ~/para/scripts/dde_weeren/bootes_hba/make_mask_from_region.py templatemask0_{source}  peel_facet_{source}.rgn templatemask0_{source}.masktmp False".format(source=source) for source in sourcelist]
    #run_parallel(mask_cmds, nthreads=1,logs='auto')

    for source_id,source in enumerate(sourcelist):
        make_facet_mask("templatemask_{source}".format(source=source), "templatemask_{source}.masktmp".format(source=source), poly[source_id], pad=True, edge=edges[source_id])


    # make a mosaic of the templates - at full res - this serves as a check of how the facets will be mosaiced
    if make_mosaic:
        print "making mask mosaic, this may take some time"
        cmd = 'python '+SCRIPTPATH+'/mos_masks.py'
        os.system(cmd)
