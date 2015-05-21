import matplotlib
#matplotlib.use('GTK')
import numpy
import os
import sys
from scipy import interpolate
import time
from subprocess import Popen, PIPE
import pyrap.tables as pt
import pyrap.images
import lofar.parmdb
from coordinates_mode import *
import pwd
from facet_utilities import run, bg
import glob
import logging

SCRIPTPATH = os.path.dirname(os.path.abspath(__file__))

pi = numpy.pi

def create_phaseshift_parset_field(msin, msout, numchanperms=20):
    ndppp_parset = msin.split('.')[0] +'ndppp_avgphaseshift_check.parset'
    os.system('rm -f ' + ndppp_parset)
    os.system('rm -rf ' + msout)
    f=open(ndppp_parset, 'w')
    f.write('msin ="%s"\n' % msin)
    f.write('msin.datacolumn = SUBTRACTED_DATA_ALL\n')
    f.write('msin.autoweight = false\n')
    f.write('msout ="%s"\n' % msout)
    f.write('msout.writefullresflag=False\n')
    f.write('steps = [uv,avg1]\n')
    #f.write('shift.type        = phaseshift\n')
    #f.write('shift.phasecenter = [%s]\n' % direction)

    f.write('uv.type=uvwflagger\n')
    f.write('uv.uvmmax=2500.0\n')
    #f.write('uv.uvmmin=20.0\n')
    f.write('avg1.type = squash\n')
    f.write('avg1.freqstep = %i\n'%numchanperms)
    f.write('avg1.timestep = 6\n')
    f.close()
    return ndppp_parset

def do_verify_subtract(mslist, res_val, source, numchanperms=20):
    ''' main function for verify_subtract '''

    msavglist = []
    for ms_id, ms in enumerate(mslist):
        msavglist.append(ms.split('.')[0] + '.' + source + '.ms.avgcheck')


    username = pwd.getpwuid(os.getuid())[0]

    ###########################################################################
    # NDPPP phase shift, less averaging (NEW: run 2 in parallel)
    b=bg(maxp=2)
    for ms_id, ms in enumerate(mslist):
        parset = create_phaseshift_parset_field(ms, msavglist[ms_id], numchanperms=numchanperms)

        #ncmd='NDPPP ' + parset+' &>'+ms.split('.')[0] +'.ndppp_avgphaseshift_check.log'
        ncmd='NDPPP ' + parset+' >'+ms.split('.')[0] +'.ndppp_avgphaseshift_check.log' + ' 2>&1'
        print 'Running',ncmd
        b.run(ncmd)

    # Check if all NDPPP processes are finished
    b.wait()
    ###########################################################################


    imsize = 2048


    ###########################################################################
    # IMAGE IN PARALLEL
    b=bg(maxp=8)
    for ms in msavglist:

        imout = 'im'+ '_residual_' + source + '_' + ms.split('.')[0]
        b.run('casapy --nogui -c '+SCRIPTPATH+'/casapy_cleanv4_checksubtract.py ' +\
                   ms + ' ' + imout + ' ' + str(imsize))
        time.sleep(20)


    # Check if all NDPPP processes are finished
    b.wait()

    #conver the images to FITS format
    for ms in msavglist:
        imout = 'im'+ '_residual_' + source + '_' + ms.split('.')[0]
        run('image2fits in=' + imout +'.image'     + ' ' + 'out='+ imout + '.fits')


    stopcal = False
    for ms in msavglist:


        # find the source which was done before the current one
        g=sorted(glob.glob('im_residual_*_' + ms.split('.')[0] + '.image'),key=os.path.getmtime,reverse=True)
        if len(g) > 1:
            pre_sourcename = g[1]
            #print 'Previous image was', pre_sourcename
        else:
            pre_sourcename=None

        image = 'im'+ '_residual_' + source + '_' + ms.split('.')[0] + '.image'

        img    = pyrap.images.image(image)
        pixels = numpy.copy(img.getdata())
        maxval = numpy.max(pixels)
        minval = numpy.min(pixels)

        maxvalpre = 1e9
        minvalpre = -maxvalpre
        if pre_sourcename is not None:
            imgpre    = pyrap.images.image(pre_sourcename)
            pixelspre = numpy.copy(imgpre.getdata())
            maxvalpre = numpy.max(pixelspre)
            minvalpre = numpy.min(pixelspre)

        logging.info('verify_subtract: For %s max, min residual were %f, %f' % (image,maxval,minval))
        print maxval, minval, ' ' + image
        if pre_sourcename is not None:
            print maxvalpre, minvalpre, ' ' + pre_sourcename
        if  (maxval > res_val) or ((maxval*0.92) > maxvalpre) :
            stopcal = True
            print 'WARNING RESIDUAL TOO LARGE, STOPPING', maxval, res_val
            print 'WARNING RESIDUAL TOO LARGE, STOPPING, previous max in image', maxvalpre

    return stopcal

if __name__=='__main__':

    mslist    = sys.argv[1:-2]
    res_val   = numpy.float(str(sys.argv[-2]))
    source    = str(sys.argv[-1])

    stop=do_verify_subtract(mslist,res_val,source)

    # reproduce old behaviour on the command line                 
    while(stop):
        time.sleep(100)
