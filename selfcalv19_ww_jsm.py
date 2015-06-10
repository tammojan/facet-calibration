import matplotlib
#matplotlib.use('GTK')
import numpy
import os
import sys
from scipy import interpolate
import time
from subprocess import Popen, PIPE, STDOUT
import pyrap.tables as pt
import pyrap.images as pim
import uuid
from numpy import pi
import pwd
import logging
from facet_utilities import run, bg

# Import functions from the main implementation
from selfcalv19_ww_cep3 import create_merged_parmdb_spline, create_merged_parmdb,\
    runbbs, create_scalarphase_parset, create_scalarphase_parset_p, create_amponly_parset,\
    create_scalarphase_parset_p2, get_group, make_image
                                
# v19
# - change ndppp job submit for better control over jobs running
# - add hardcoded groups (get_group)


def find_imagenoise(imagename):
    """
    Finds the noise level of an image
    """
    im    = pim.image(imagename)
    image = numpy.copy(im.getdata())
    mean, rms =  meanclip(image)
    #im.close()
    return rms,  numpy.abs(numpy.max(image)/numpy.min(image))


def meanclip(indata, clipsig=4.0, maxiter=10, converge_num=0.001, verbose=0):
   """
   Computes an iteratively sigma-clipped mean on a
   data set. Clipping is done about median, but mean
   is returned.

   .. note:: MYMEANCLIP routine from ACS library.

   :History:
       * 21/10/1998 Written by RSH, RITSS
       * 20/01/1999 Added SUBS, fixed misplaced paren on float call, improved doc. RSH
       * 24/11/2009 Converted to Python. PLL.

   Examples
   --------
   >>> mean, sigma = meanclip(indata)

   Parameters
   ----------
   indata: array_like
       Input data.

   clipsig: float
       Number of sigma at which to clip.

   maxiter: int
       Ceiling on number of clipping iterations.

   converge_num: float
       If the proportion of rejected pixels is less than
       this fraction, the iterations stop.

   verbose: {0, 1}
       Print messages to screen?

   Returns
   -------
   mean: float
       N-sigma clipped mean.

   sigma: float
       Standard deviation of remaining pixels.

   """
   # Flatten array
   skpix = indata.reshape( indata.size, )
 
   ct = indata.size
   iiter = 0
   c1 = 1.0
   c2 = 0.0
 
   while (c1 >= c2) and (iiter < maxiter):
       lastct = ct
       medval = numpy.median(skpix)
       sig = numpy.std(skpix)
       wsm = numpy.where( abs(skpix-medval) < clipsig*sig )
       ct = len(wsm[0])
       if ct > 0:
           skpix = skpix[wsm]
 
       c1 = abs(ct - lastct)
       c2 = converge_num * lastct
       iiter += 1
   # End of while loop
 
   mean  = numpy.mean( skpix )
   sigma = robust_sigma( skpix )
 
   #if verbose:
   prf = 'MEANCLIP:'
   logging.debug('%s %.1f-sigma clipped mean' % (prf, clipsig))
   logging.debug('%s Mean computed in %i iterations' % (prf, iiter))
   logging.debug('%s Mean = %.6f, sigma = %.6f' % (prf, mean, sigma))
 
   return mean, sigma


def robust_sigma(in_y, zero=0):
    """
    Calculate a resistant estimate of the dispersion of
    a distribution. For an uncontaminated distribution,
    this is identical to the standard deviation.
    
    Use the median absolute deviation as the initial
    estimate, then weight points using Tukey Biweight.
    See, for example, Understanding Robust and
    Exploratory Data Analysis, by Hoaglin, Mosteller
    and Tukey, John Wiley and Sons, 1983.
    
    .. note:: ROBUST_SIGMA routine from IDL ASTROLIB.
    
    :History:
        * H Freudenreich, STX, 8/90
        * Replace MED call with MEDIAN(/EVEN), W. Landsman, December 2001
        * Converted to Python by P. L. Lim, 11/2009
    
    Examples
    --------
    >>> result = robust_sigma(in_y, zero=1)
    
    Parameters
    ----------
    in_y: array_like
        Vector of quantity for which the dispersion is
        to be calculated
    
    zero: int
        If set, the dispersion is calculated w.r.t. 0.0
        rather than the central value of the vector. If
        Y is a vector of residuals, this should be set.
    
    Returns
    -------
    out_val: float
        Dispersion value. If failed, returns -1.
    
    """
    # Flatten array
    y = in_y.reshape(in_y.size, )
    
    eps = 1.0E-20
    c1 = 0.6745
    c2 = 0.80
    c3 = 6.0
    c4 = 5.0
    c_err = -1.0
    min_points = 3
 
    if zero:
        y0 = 0.0
    else:
        y0 = numpy.median(y)
 
    dy    = y - y0
    del_y = abs( dy )
    
    # First, the median absolute deviation MAD about the median:
    
    mad = numpy.median( del_y ) / c1
    
    # If the MAD=0, try the MEAN absolute deviation:
    if mad < eps:
        mad = numpy.mean( del_y ) / c2
    if mad < eps:
        return 0.0
 
    # Now the biweighted value:
    u  = dy / (c3 * mad)
    uu = u*u
    q  = numpy.where(uu <= 1.0)
    count = len(q[0])
    if count < min_points:
        logging.warning('ROBUST_SIGMA: This distribution is TOO WEIRD! Returning {}'.format(c_err))
        return c_err
 
    numerator = numpy.sum( (y[q]-y0)**2.0 * (1.0-uu[q])**4.0 )
    n    = y.size
    den1 = numpy.sum( (1.0-uu[q]) * (1.0-c4*uu[q]) )
    siggma = n * numerator / ( den1 * (den1 - 1.0) )
    
    if siggma > 0:
        out_val = numpy.sqrt( siggma )
    else:
        out_val = 0.0
 
    return out_val


def do_selfcal(mslist, cluster, atrous_do, imsize, nterms, cellsizetime_a, cellsizetime_p,
               TECi, clocki, HRi, region, clusterdesc, dbserver, dbuser, dbname, SCRIPTPATH):

    TEC  = False
    FFT  = False
    clock= False
    HR   = False # high dynamic range

    if TECi == "True":
        TEC = True

    if clocki == "True":
        clock = True

    if HRi == "HD":
        HR = True

    if imsize <=2048:
        wplanes = 1
        FFT = True  # FFT image into MODEL_DATA
        if imsize >= 512:
            wplanes = 64
        if imsize > 799:
            wplanes = 96
        if imsize > 1023:
            wplanes = 128
        if imsize > 1599:
            wplanes = 180
        if imsize > 1800:
            wplanes = 196
        if imsize > 2049:
            wplanes = 256
        #if imsize > 3000:
        #   wplanes = 448
        #if imsize > 4095:
        #   wplanes = 512


    logging.info('mslist {}'.format(mslist))
    logging.info('source {}'.format(cluster))
    logging.info('atrous_do {}'.format(atrous_do))
    logging.info('imsize {} '.format(imsize))
    logging.info('TEC is {} and clock is {}'.format(TEC, clock))

    msinputlist = ''
    for m in mslist:
        msinputlist = msinputlist + ' ' + m


    #if len(mslist) == 29:
        #group = "9,10,10"
    #elif len(mslist) == 28:
        #group = "7,7,7,7"
    #elif len(mslist) == 20:
        ##group = "10,10"
        #group = "7,7,6"
    #elif len(mslist) == 16:
        #group = "8,8"
    #else:
        #group = str(len(mslist))


    group = get_group(mslist)
    logging.info('GROUP {}'.format(group))

    uvrange = '80'
    #uvrange = '400'

    merge_parmdb = True
    phasors      = False   # if true only solve for amps on long timescales
    smooth       = False # seems that smooth does not help the selfcal (various reasons for that)
                         # 1. boundaries of flagged vs non-flagged data are sharp (should not be smoothed)
                         # 2. there sre some strong ampl various at low elevations
    smooth       = True # sometimes almost 0.0 amplitude, causes ripples
    phasezero    = True # reset phases from ap calibration
    
    # Loop parameters
    number_forced_selfcalcycles = 8
    rms_old          = 1.e9 # bad values to start with
    dynamicrange_old = 1. # low value to start with so we get into the while loop
    factor           = 1.0125 # demand 1.25% improvement
    im_count         = 4
    max_selfcalcycles = 16


    #####################
    #####################
    
    
    #### MAKE IMAGE 0 ###
    logging.info('Make image 0')
    imout,mask = make_image(mslist, cluster, '0', 10, 6, nterms, atrous_do, imsize, region, SCRIPTPATH)

    #####################

    ### CALIBRATE WITH BBS PHASE ONLY 1 ###
    # create skymodel for BBS
    run(SCRIPTPATH+'/casapy2bbs.py -m '+ mask + ' ' +'-t ' + str(nterms)+ ' ' + imout+'.model ' +  imout+'.skymodel')
    if FFT:
        os.system('casapy --nogui -c '+SCRIPTPATH+'/ft_v2.py ' + msinputlist + ' ' + imout+'.model' \
                  + ' ' + str(nterms) + ' '+ str(wplanes))

    # phase only calibrate
    skymodel = imout+'.skymodel'
    parset   = create_scalarphase_parset(cellsizetime_p, TEC, clock, group, FFT, uvrange)

    runbbs(mslist, skymodel, parset, 'instrument', False, TEC, clusterdesc, dbserver, dbuser, dbname)
    #NOTE WORK FROM MODEL_DATA (contains correct phase data from 10SB calibration)
    ######################################


    ### MAKE IMAGE 1 ###
    logging.info('Make image 1')
    imout,mask = make_image(mslist, cluster, '1', 15, 15, nterms, atrous_do, imsize, region, SCRIPTPATH)
    ####################


    ### CALIBRATE WITH BBS PHASE ONLY 2 ###
    # create skymodel for BBS
    run(SCRIPTPATH+'/casapy2bbs.py -m '+ mask + ' ' +'-t ' + str(nterms)+ ' ' + imout+'.model ' +  imout+'.skymodel')
    if FFT:
        run('casapy --nogui -c '+SCRIPTPATH+'/ft_v2.py ' + msinputlist + ' ' + imout+'.model' \
                  + ' ' + str(nterms) + ' '+ str(wplanes))


    # phase only calibrate
    skymodel = imout+'.skymodel'
    parset   = create_scalarphase_parset(cellsizetime_p, TEC, clock, group, FFT, uvrange)

    runbbs(mslist, skymodel, parset, 'instrument', False, TEC,clusterdesc, dbserver, dbuser, dbname) #NOTE WORK FROM MODEL_DATA (contains correct phase data from 10SB calibration)
    ######################################


    ### MAKE IMAGE 2 ###
    logging.info('Make image 2')
    imout,mask = make_image(mslist, cluster, '2', 15, 15, nterms, atrous_do, imsize, region, SCRIPTPATH)
    ####################
    
    ### CALIBRATE WITH BBS PHASE+AMP 1 ###
    run(SCRIPTPATH+'/casapy2bbs.py -m '+ mask + ' ' +'-t ' + str(nterms)+ ' ' + imout+'.model ' +  imout+'.skymodel')
    if FFT:
        run('casapy --nogui -c '+SCRIPTPATH+'/ft_v2.py ' + msinputlist + ' ' + imout+'.model' \
                  + ' ' + str(nterms) + ' '+ str(wplanes))

    skymodel = imout+'.skymodel'
    parset   = create_scalarphase_parset_p(cellsizetime_p, TEC, clock, group, FFT, uvrange)
    # solve +apply phases
    runbbs(mslist, skymodel, parset, 'instrument_phase0', False, TEC, clusterdesc, dbserver, dbuser, dbname)

    # solve amps
    parmdb = 'instrument_amps0'
    parset = create_amponly_parset(cellsizetime_a, FFT, uvrange)
    runbbs(mslist, skymodel, parset, parmdb, False, False, clusterdesc, dbserver, dbuser, dbname)

    for ms in mslist:
        # remove outliers from the solutions
        if phasors:
            run('python '+SCRIPTPATH+'/smoothcal_rx42.py ' + ms + ' ' + ms+'/'+parmdb + ' ' + ms+'/'+parmdb+'_smoothed'+' > '+ms+'_'+parmdb+'_smoothed.log')
        else:
            run('python '+SCRIPTPATH+'/smoothcal_a2256_nophasors.py ' + ms + ' ' + ms+'/'+parmdb + ' ' + ms+'/'+parmdb+'_smoothed'+' > '+ms+'_'+parmdb+'_smoothed.log')

    # apply amps
    if smooth:
        runbbs(mslist, skymodel,SCRIPTPATH+'/apply_amplitudeonly.parset', parmdb+'_smoothed', True, False, clusterdesc, dbserver, dbuser, dbname)
    else:
        runbbs(mslist, skymodel,SCRIPTPATH+'/apply_amplitudeonly.parset', parmdb, True, False, clusterdesc, dbserver, dbuser, dbname)

    ### MAKE IMAGE 3 ###
    logging.info('Make image 3')
    imout,mask = make_image(mslist, cluster, '3', 10, 10, nterms, atrous_do, imsize, region, SCRIPTPATH)

    
    ####################
    # LOOP 
    ####################
    
    rms, dynamicrange =  find_imagenoise(imout + '.image')
    
    while (((dynamicrange/factor) > dynamicrange_old) or ((rms*factor) < rms_old)):
        logging.info('Starting selfcal loop')
        #### CALIBRATE  BBS PHASE+AMP 2 (LOOP) ###
        # make model
        run(SCRIPTPATH+'/casapy2bbs.py -m '+ mask + ' ' +'-t ' + str(nterms)+ ' ' + imout+'.model ' +  imout+'.skymodel')
        if FFT:
            run('casapy --nogui -c '+SCRIPTPATH+'/ft_v2.py ' + msinputlist + ' ' + imout+'.model' \
                    + ' ' + str(nterms) + ' '+ str(wplanes))

        #parmdb keep from previous step
        skymodel = imout+'.skymodel'

        # phase only cal
        skymodel = imout+'.skymodel'
        parset   = create_scalarphase_parset_p(cellsizetime_p, TEC, clock, group, FFT, uvrange)
        runbbs(mslist, skymodel, parset, 'instrument_phase1', False, TEC, clusterdesc, dbserver, dbuser, dbname)

        # solve amps
        parmdb   = 'instrument_amps1'
        parset = create_amponly_parset(cellsizetime_a, FFT, uvrange)
        runbbs(mslist, skymodel, parset,parmdb, False, False, clusterdesc, dbserver, dbuser, dbname)

        for ms in mslist:
            # remove outliers from the solutions
            if phasors:
                run('python '+SCRIPTPATH+'/smoothcal_rx42.py ' + ms + ' ' + ms+'/'+parmdb + ' ' + ms+'/'+parmdb+'_smoothed'+' > '+ms+'_'+parmdb+'_smoothed.log')
            else:
                run('python '+SCRIPTPATH+'/smoothcal_a2256_nophasors.py ' + ms + ' ' + ms+'/'+parmdb + ' ' + ms+'/'+parmdb+'_smoothed'+' > '+ms+'_'+parmdb+'_smoothed.log')

        # apply amps
        if smooth:
            runbbs(mslist, skymodel,SCRIPTPATH+'/apply_amplitudeonly.parset',parmdb+'_smoothed', True, False, clusterdesc, dbserver, dbuser, dbname)
        else:
            runbbs(mslist, skymodel,SCRIPTPATH+'/apply_amplitudeonly.parset',parmdb, True, False, clusterdesc, dbserver, dbuser, dbname)

        ### MAKE IMAGE #N ###
        logging.info('Make image {}'.format(im_count))
        imout,mask = make_image(mslist, cluster, str(im_count), 10, 10, nterms, atrous_do, imsize, region, SCRIPTPATH)

        
        im_count += 1
        
        # save previous values to compare with
        rms_old          = rms
        dynamicrange_old = dynamicrange 
        if nterms < 2:    
            rms, dynamicrange =  find_imagenoise(imout + '.image')
        else:
            rms, dynamicrange =  find_imagenoise(imout + '.image.tt0')
        logging.info('IMAGE STATISTICS {}, {}'.format(rms, dynamicrange))

        if im_count < number_forced_selfcalcycles:
            rms_old          = 1.e9 # bad values to start with
            dynamicrange_old = 1.
            logging.debug("Count below the number of forced selfcal cycles. Force new loop.")
        
        if im_count >= max_selfcalcycles:
            logging.info("Maximum number of cycles ({}) reached. Leaving selfcal loop.".format(max_selfcalcycles))
            break


    ### CREATE FINAL MODEL ###
    logging.info('Create final model')

    skymodelf= 'im_cluster'+cluster+ '.final.skymodel'
    run(SCRIPTPATH+'/casapy2bbs.py -m '+ mask + ' ' +'-t ' + str(nterms)+ ' ' + imout+'.model ' +  skymodelf)
    if FFT:
        run('casapy --nogui -c '+SCRIPTPATH+'/ft_v2.py ' + msinputlist + ' ' + imout+'.model' \
                  + ' ' + str(nterms) + ' '+ str(wplanes))

    ### CREATED MERGED PARMDB SCALARPHASE+AMPS ###
    ### INCLUDES SPLINE INTERPOLARION OF AMPS ###
    if merge_parmdb:
        logging.info('Merge parmdb')
        if phasors:
            dummyparset = SCRIPTPATH+'/scalarphase+amp.parset'
        else:
            if TEC:
                dummyparset = SCRIPTPATH+'/scalarphase+ap+TEC.parset'
            else:
                dummyparset = SCRIPTPATH+'/scalarphase+ap.parset'

        if TEC:
            if clock:
                dummyparmdb = 'instrument_template_TECclock'
            else:
                dummyparmdb = 'instrument_template_Gain_TEC_CSphase'

        #if not os.path.isdir(dummyparmdb):
            #runbbs([mslist[0]], skymodel,dummyparset, dummyparmdb, True, False)
        # TO SPEED THINGS UP, hard coded for BOOTES - i.e. the above has already been run
        for ms in mslist:
            os.system('rm -rf ' + ms +'/' + dummyparmdb)
            os.system('cp -r ' + dummyparmdb + ' ' +  ms + '/instrument_template')

        if smooth:
            parmdb_a    = 'instrument_amps1_smoothed'  # last/best ampplitude(+phase) parmdb
        else:
            parmdb_a    = 'instrument_amps1'  # last/best ampplitude(+phase) parmdb
        parmdb_p    = 'instrument_phase1'          # last/best CommonScalarPhase parmdb
        parmdbout   = 'instrument_merged'

        #reset in case of instrument_template_TECclock
        dummyparmdb = 'instrument_template'

        for ms in mslist:
            create_merged_parmdb(ms, ms+'/'+parmdb_a, ms+'/'+parmdb_p, ms+'/'+dummyparmdb,ms+'/'+parmdbout,cellsizetime_a,cellsizetime_p)


if __name__=="__main__":

    # arguments are mslist, cluster, atrous_do, imsize, nterms, cellsizetime_a, cellsizetime_p, TECi, clocki, HRi, region, clusterdesc, dbserver, dbuser, dbname

    mslist    = sys.argv[1:-14]
    cluster   = str(sys.argv[-14])
    atrous_do = str(sys.argv[-13])
    imsize    = numpy.int(sys.argv[-12])

    nterms                  = numpy.int(sys.argv[-11])  # only 1 to 3 is supported !!
    cellsizetime_a          = numpy.int(sys.argv[-10])
    cellsizetime_p          = numpy.int(sys.argv[-9])
    TECi                    = str(sys.argv[-8])
    clocki                  = str(sys.argv[-7])
    HRi                     = str(sys.argv[-6])
    region                  = str(sys.argv[-5])
    clusterdesc = sys.argv[-4]
    dbserver=sys.argv[-3]
    dbuser=sys.argv[-2]
    dbname=sys.argv[-1]

    SCRIPTPATH = os.path.dirname(os.path.abspath(__file__))

    do_selfcal(mslist,cluster,atrous_do,imsize,nterms,cellsizetime_a,cellsizetime_p,TECi,clocki,HRi,region,clusterdesc,dbserver,dbuser,dbname,SCRIPTPATH)
