#!/usr/bin/env python

import matplotlib
#matplotlib.use('GTK')
import numpy
import os
import sys
import lofar.parmdb
from scipy import interpolate
import time
#import subprocess
from subprocess import Popen, PIPE, STDOUT
import pyrap.tables as pt
import pyrap.images
import pwd
import logging
from numpy import pi

def image_size_from_mask(mask):
    im = pyrap.images.image(mask)
    sh = im.shape()
    if sh[-1] != sh[-2]:
        print "image is not square!"
        print sh[-1], sh[-2]
    npix = sh[-1]
    return npix

def rad2deg(rad):
    return rad*180./numpy.pi

def image_centre_from_mask(mask):
    im = pyrap.images.image(mask)
    ic = im.coordinates()
    lon,lat = ic.get_referencevalue()[2]

    ra = rad2deg(lat)
    dec = rad2deg(lon)

    centre = '{x:.4f},{y:.4f}'.format(x=ra,y=dec)
    facet_centre = '{x:s},{y:s}'.format(x=ra_to_str(ra,delim='h'),y=dec_to_str(dec,delim='d'))
    print centre
    print facet_centre
    return facet_centre

def make_image_wsclean_nomask(mslist, cluster, callnumber, threshpix, threshisl, nterms, atrous_do, imsize, inputmask, mscale, region,cellsize,uvrange,wsclean,WSCleanRobust):

    if imsize is None:
        imsize = image_size_from_mask(inputmask)

    niter   = numpy.int(20000 * (numpy.sqrt(numpy.float(len(mslist)))))
    cellsizeim = str(cellsize) +'arcsec'

    depth =  1e-3*0.7 / (numpy.sqrt(numpy.float(len(mslist))))
    print 'Cleaning to a noise level of',depth,'Jy: niter is',niter

    cleandepth2 = str(depth)     #+ 'mJy'

    wideband = False
    if len(mslist) > 5:
        wideband = True

    ms = ''
    for m in mslist:
        ms = ms + ' ' + m

    imout = 'im'+ callnumber +'_cluster'+cluster

    os.system('rm -rf ' + imout + '-*')


    outms      = 'field-'+cluster+'.ms'
    parsetname = 'concatforwsclean-'+cluster+'.parset'

    msinstr = ""

    for ms_id, ms in enumerate(mslist):
        msinstr = msinstr + "'" + ms + "'"
        if ms_id < len(mslist)-1:
            msinstr = msinstr + ", "
    os.system('rm -rf ' + parsetname)
    f=open(parsetname, 'w')
    f.write('msin = [%s]\n' % msinstr)
    f.write('msin.datacolumn = DATA\n')
    f.write('msin.missingdata=True\n')
    f.write('msin.orderms=False\n')
    f.write('msout=%s\n' % outms)
    f.write('steps=[]\n')
    f.close()
    os.system('rm -rf ' + outms)
    run('NDPPP ' + parsetname)

    if wideband:
        channelsout =  1 # there is a factor of 5 averaging
        cmd1 = wsclean + ' -reorder -name ' + imout + ' -size ' + str(imsize) + ' ' + str(imsize) + ' '
        cmd2 = '-scale ' + cellsizeim + ' -weight briggs '+str(WSCleanRobust)+' -niter ' + str(niter) + '-cleanborder 0 -threshold '+ cleandepth2 + ' '
        cmd3 = '-minuv-l '+ str(uvrange) + ' -casamask ' +  inputmask + ' '\
               +' -mgain 0.6 -fitbeam -datacolumn DATA -no-update-model-required -joinchannels -channelsout ' +\
               str(channelsout) + ' '  + outms
    else:
        cmd1 = wsclean + ' -reorder -name ' + imout + ' -size ' + str(imsize) + ' ' + str(imsize) + ' '
        cmd2 = '-scale ' + cellsizeim + ' -weight briggs '+str(WSCleanRobust)+' -niter ' + str(niter) + ' -cleanborder 0 -threshold '+ cleandepth2 + ' '
        cmd3 = '-minuv-l '+ str(uvrange) +' -casamask ' +  inputmask + ' -mgain 0.6 -fitbeam -datacolumn DATA -no-update-model-required ' + outms

    print cmd1+cmd2+cmd3
    run(cmd1+cmd2+cmd3)

    finalim=pyrap.images.image(imout+'-image.fits')
    finalim.saveas(imout +'.image')

    return imout, None, imsize


def create_phaseshift_parset_field_avg(msin, msout, source, direction):
    ndppp_parset = msin.split('.')[0] +'ndppp_avgphaseshift_field1.parset'
    os.system('rm -f ' + ndppp_parset)

    f=open(ndppp_parset, 'w')
    f.write('msin ="%s"\n' % msin)
    f.write('msin.datacolumn = CORRECTED_DATA\n')
    f.write('msin.autoweight = false\n')
    f.write('msout ="%s"\n' % msout)
    f.write('msout.writefullresflag=False\n')
    f.write('steps = [shift,avg1]\n')
    f.write('shift.type        = phaseshift\n')
    f.write('shift.phasecenter = [%s]\n' % direction)
    f.write('avg1.type = squash\n')
    f.write('avg1.freqstep = 5\n')
    f.write('avg1.timestep = 3\n')
    f.close()
    return ndppp_parset

def create_add_parset_field(outputcolumn, TEC, clock):
    bbs_parset = 'add.parset'
    os.system('rm -f ' + bbs_parset)
    f=open(bbs_parset, 'w')

    chunksize = 175

    f.write('Strategy.InputColumn = SUBTRACTED_DATA_ALL\n')
    f.write('Strategy.ChunkSize   = %s\n' % chunksize)
    f.write('Strategy.UseSolver   = F\n')
    f.write('Strategy.Steps       = [add]\n\n\n')
    f.write('Step.add.Model.Sources                   = [@ADDED_DATA_SOURCE]\n')  ## comes from ft_allbands, then copy over columns
    f.write('Step.add.Model.Cache.Enable              = T\n')
    f.write('Step.add.Model.Phasors.Enable            = F\n')
    f.write('Step.add.Model.DirectionalGain.Enable    = F\n')
    f.write('Step.add.Model.Gain.Enable               = T\n')
    f.write('Step.add.Model.Rotation.Enable           = F\n')
    f.write('Step.add.Model.CommonScalarPhase.Enable  = T\n')
    if TEC   == "True":
        f.write('Step.add.Model.TEC.Enable  = T\n')
    if clock == "True":
        f.write('Step.add.Model.Clock.Enable  = T\n')
    f.write('Step.add.Model.CommonRotation.Enable     = F\n')
    f.write('Step.add.Operation                       = ADD\n')
    f.write('Step.add.Model.Beam.Enable               = F\n')
    f.write('Step.add.Output.WriteCovariance          = F\n')
    f.write('Step.add.Output.Column                   = %s\n' % outputcolumn)
    f.close()
    return bbs_parset

#############
# MAIN CODE #
#############

if __name__=='__main__':

    if len(sys.argv)<2:
        raise Exception('Give the path to the setup code for the facet')

    username = pwd.getpwuid(os.getuid())[0]
    print 'Using',sys.argv[1],'as the setup code'
    execfile(sys.argv[1])
    print 'script path is',SCRIPTPATH
    try:
        StartAtStep
    except NameError:
        print 'No starting step specified, begin at the beginning'
        StartAtStep='preSC'

    if (len(do_sources) > 1) and (StartAtStep != 'preSC'):
       print 'For StartAtStep "' + StartAtStep + '" can only do a single source direction'
       raise Exception('do_sources not compatible with StartAtStep')
        
    try:
        WSCleanRobust
    except NameError:
        WSCleanRobust=-0.25 # default preserves old value
        print 'No WSClean robust set, defaulting to',WSCleanRobust

    try:
        BlankField
    except NameError:
        BlankField=False
        print 'BlankField not set, defaulting to',BlankField

    try:
        NoMask
    except NameError:
        NoMask=False
        print 'NoMask not set, defaulting to',NoMask

    try:
        WScleanWBgroup
    except NameError:
        print 'WScleanWBgroup is not set, not using wideband clean algorithm'
        # only print message here, because wideband is not used when len(mslist) <= WScleanWBgroup:
        WScleanWBgroup = 1000 # large number so wideband is never used

    print 'importing local modules....'

    if SCRIPTPATH not in sys.path:
        sys.path.append(SCRIPTPATH)
    from facet_utilities import run, bg, angsep, dec_to_str, ra_to_str
    from coordinates_mode import *
    from verify_subtract_v5 import do_verify_subtract
    from doDDE_v21_a2256 import create_phaseshift_parset_full,do_fieldFFT,runbbs,make_image_wsclean

    if os.path.exists("logging.conf"):
        logging.config.fileConfig('logging.conf')
        logger = logging.getLogger()
    else:
        # Start
        logger = logging.getLogger()
        logger.setLevel(logging.DEBUG)   
        formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
        # Log to STDOUT
        ch = logging.StreamHandler()
        ch.setFormatter(formatter)
        logger.addHandler(ch)
        # Log to file
        file_name = "dde.log"
        fh = logging.FileHandler(file_name) 
        fh.setFormatter(formatter)
        logger.addHandler(fh)
        logging.info('\n')


    source_info_rec = numpy.genfromtxt(peelsourceinfo, dtype="S10,S25,S5,S5,i8,i8,i8,i8,S2,S10,S10", names=["sourcelist","directions","atrous_do","mscale_field","imsizes","cellsizetime_p","cellsizetime_a","fieldsize","dynamicrange","regionselfc","regionfield"])

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

    sourcelist = sourcelist.tolist()



    mslistorig = ["{name:s}_SB{b1:03d}-{b2:03d}.{res:s}.ms".format(name=NAME,res=RES,b1=b,b2=b+9) for b in BANDS]
    mslistorigstr = ' '.join(mslistorig)

    mslist= [ms for ms in mslistorig if  os.path.isdir(ms)]  # filter out datasets that do not exist (takes also care of freq.gaps in field subtract))

    msliststr = ' '.join(mslist)



    tt = pt.table(mslist[0] + '/FIELD')
    pointingcenter = tt.getcol('REFERENCE_DIR')[0][0]
    pointingcenter = str(pointingcenter[0]) +'rad,' + str(pointingcenter[1])+'rad'
    print pointingcenter
    tt.close()

    freq_tab     = pt.table(mslist[0]  + '/SPECTRAL_WINDOW')
    numchanperms = freq_tab.getcol('NUM_CHAN')[0]
    logging.info('Number of channels per ms is {:d}'.format(numchanperms))
    freq_tab.close()

    ##################################

    logging.info('\n')
    logging.info('#######################################################################\n')
    logging.info('Running DDE-reimage')
    logging.info('Doing sources: '+','.join(do_sources))
    logging.info('Using MSlist: '+msliststr)

    for source in do_sources:
        ## I1

        source_id = sourcelist.index(source)

        print 'Re-imaging facet:', source
        logging.info('')
        logging.info('Re-imaging facet: '+ source )

        logging.info("removing any existing facet field1 images")
        os.system("rm -rf imfield1*_cluster"+source+"*")
        if StartAtStep in ['preSC','preFACET']:
            logging.info("removing any existing facet field1 imaging average MS")
            os.system("rm -rf *."+source+".ms.avgfield1*")

        #check if allbands.concat.shifted_'+source+'.ms is present
        if os.path.isdir('allbands.concat.shifted_'+source+'.ms'):
            logging.info('allbands.concat.shifted_'+source+'.ms already exists')
            logging.info('...but continuing because we are re-imaging')

        if not os.path.isdir('allbands.concat.ms'):
            print 'allbands.concat.ms does not exist'
            raise Exception('make measurement set and then restart')


        dummyskymodel   = SCRIPTPATH+'/dummy.skymodel' ## update every time again with new source, not used, just a dummy for correct

        msfieldavgfacetlist1 = []
        for ms in mslist:
            msfieldavgfacetlist1.append(ms.split('.')[0] + '.' + source + '.ms.avgfield1.facetdir')

        # combined SC and DDE solutions (basename - for use in runbbs - uses ms/basename)
        parmdb_master_out  = "instrument_master_" + source   # reset because runbbs uses basename of ms

        # set image mask region
        output_template_im = 'templatemask_' + source +'.masktmp'
        if not os.path.exists(output_template_im):
            raise  Exception('facet mask missing: '+output_template_im)

        # set directions #
        #selfcaldir = directions[source_id]
        facetdir = image_centre_from_mask(output_template_im)
        facetsize = image_size_from_mask(output_template_im)

        #logging.info("Selfcal direction: "+selfcaldir)
        logging.info("Facet direction: "+facetdir)
        logging.info("facetmask: "+str(facetsize))

        #continue

        ## NOTE: addback now needs to be done as subtract in reverse... with FT and allbands concat


        ## STEP 1: check  ##
        if StartAtStep in ['preSC']:
        # if we didn't keep the allbands.concat.shifted
            if not os.path.exists('allbands.concat.shifted_'+source+'.ms'):
                parset = create_phaseshift_parset_full('allbands.concat.ms', 'allbands.concat.shifted_'+source+'.ms', facetdir,'DATA')

                ndppplog = parset.replace('.parset','.log')
                #ndpppcmd = 'NDPPP ' + parset + ' > '+ ndppplog + ' 2>&1'
                #ndppprc = os.system(ndpppcmd)
                ndpppcmd = 'NDPPP ' + parset + ' > '+ ndppplog + ' 2>&1'
                run(ndpppcmd)

        ### STEP 3: prep for facet ##

        if StartAtStep in ['preSC','preFACET']:
            # image model to add back
            imsizef = image_size_from_mask(output_template_im)
            imout = 'im'+ 'field0' +'_cluster'+source

            if not os.path.exists(imout+'-model.fits'):
                raise Exception(imout+'-model.fits is missing')


            logging.info('running ft: '+imout)


            # DO THE FFT
            do_fieldFFT('allbands.concat.shifted_'+source+'.ms',imout, imsizef, cellsize, wsclean, mslist, WSCleanRobust,WScleanWBgroup, numchanperms)
            logging.info('FFTed model of DDE facet: ' + source)

            # SHIFT PHASE CENTER BACK TO ORIGINAL
            parset = create_phaseshift_parset_full('allbands.concat.shifted_'+source+'.ms', 'allbands.concat.shiftedback_'+source+'.ms',\
                                        pointingcenter,'MODEL_DATA')

            ndppplog = parset.replace('.parset','.log')
            ndpppcmd = 'NDPPP ' + parset  + ' > '+ ndppplog + ' 2>&1 '
            run(ndpppcmd)
            os.system('rm -rf allbands.concat.shifted_'+source+'.ms') # clean up

            # Add MODEL_DATA (allbands.concat.shiftedback_'+source+'.ms) into ADDED_DATA_SOURCE from mslist

            freq_tab1= pt.table('allbands.concat.ms' + '/SPECTRAL_WINDOW')
            numchan1    = freq_tab1.getcol('NUM_CHAN')
            freq_tab2= pt.table(mslist[0] + '/SPECTRAL_WINDOW')
            numchan2    = freq_tab2.getcol('NUM_CHAN')
            freq_tab1.close()
            freq_tab2.close()

            if (numchan1[0]) == (numchan2[0]*len(mslist)):
                run('python '+SCRIPTPATH+'/copy_over_columns.py '+ msliststr +\
                        ' ' +'allbands.concat.shiftedback_'+source+'.ms'+' ' + 'ADDED_DATA_SOURCE')
            else:
                run('python '+SCRIPTPATH+'/copy_over_columns.py '+ mslistorigstr +\
                        ' ' +'allbands.concat.shiftedback_'+source+'.ms'+' ' + 'ADDED_DATA_SOURCE')


            os.system('rm -rf allbands.concat.shiftedback_'+source+'.ms') # clean up

            addfieldparset   = create_add_parset_field('MODEL_DATA', TEC, clock) # SUBTRACTED_DATA_ALL + ADDED_DATA_SOURCE = MODEL_DATA


            logging.info('adding sources with solutions')
            runbbs(mslist, dummyskymodel, addfieldparset, parmdb_master_out+'_norm', False) # replace-sourcedb not needed since we use "@column"



            logging.info('correcting with solutions')
            ### STEP 3: prep for facet ##
            # apply master solutions, put in CORRECTED_DATA
            if TEC == "True":
                if clock == "True":
                    runbbs(mslist, dummyskymodel, SCRIPTPATH+'/correctfield2+TEC+clock.parset',parmdb_master_out+'_norm', False)
                else:
                    runbbs(mslist, dummyskymodel, SCRIPTPATH+'/correctfield2+TEC.parset',parmdb_master_out+'_norm', False)
            else:
                runbbs(mslist, dummyskymodel, SCRIPTPATH+'/correctfield2.parset',parmdb_master_out+'_norm', False)


            ###########################################################################
            # NDPPP phase shift, less averaging (NEW: run 2 in parallel)
            # CHANGE (v7) phaseshift to sc dir and avg
            # then phase shift to facetdir
            b=bg(maxp=2)
            for ms_id, ms in enumerate(mslist):
                phaseshiftfieldparset = create_phaseshift_parset_field_avg(ms, msfieldavgfacetlist1[ms_id], source, facetdir)
                ndppplog = phaseshiftfieldparset.replace('.parset','.log')
                ndpppcmd = 'NDPPP ' + phaseshiftfieldparset + ' > '+ ndppplog + ' 2>&1'
                b.run(ndpppcmd)
            b.wait()

        if StartAtStep in ['preSC','preFACET','doFACET']:
            logging.info('making image')
            ### STEP 4a: do facet ##
            # make large field image
            if NoMask:
                make_image_wsclean_nomask(msfieldavgfacetlist1, source,
                                          'field1', 5, 3, nterms, 'True',
                                          None, output_template_im,
                                          mscale_field[source_id],regionfield[source_id],cellsize,
                                          uvrange,wsclean,WSCleanRobust)
            else:
                make_image_wsclean(msfieldavgfacetlist1, source, 'field1',
                                   5, 3, nterms, 'True', None,
                                   output_template_im,
                                   mscale_field[source_id],
                                   regionfield[source_id],cellsize,
                                   uvrange, wsclean, WSCleanRobust,
                                   BlankField, WScleanWBgroup,
                                   numchanperms,path=SCRIPTPATH)
