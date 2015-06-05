import matplotlib
#matplotlib.use('GTK')
import numpy
import os
import sys
from scipy import interpolate
import time
from subprocess import Popen, PIPE, STDOUT
import pyrap.tables as pt
import uuid
from numpy import pi
import pwd
import logging
from facet_utilities import run, bg

# Import functions from the main implementation
from selfcalv19_ww_cep3 import create_merged_parmdb_spline, create_merged_parmdb,\
    runbbs, create_scalarphase_parset, create_scalarphase_parset_p, create_amponly_parset,\
    create_scalarphase_parset_p2, get_group
                                
# v19
# - change ndppp job submit for better control over jobs running
# - add hardcoded groups (get_group)


def make_image(mslist, cluster, callnumber, threshpix, threshisl, nterms, atrous_do, imsize, region, SCRIPTPATH):
    """
    Test modifications to make_image
    """
    do_mask = True
    niter   = 1000 # 7500 causes nasty clean artifacts
    mscale  = 'False'
    if atrous_do == 'True':
        mscale = 'True'

    average = True  # average the data a lot to speed up the imaging,
    # ONLY average for small FOV, otherwise timesmearing is a problem
    #average data
    if average:
        b=bg(maxp=2)
        for ms in (mslist):
            ndppp_parset = ms + '_NDPPP.parset'
            ndppplog = ndppp_parset.replace('.parset','.log')
            os.system('rm -f ' + ndppp_parset)
            output = ms + '.tmpavg'
            os.system('rm -rf ' + output)
            f=open(ndppp_parset, 'w')
            f.write('msin = %s\n' % ms)
            f.write('msin.datacolumn = CORRECTED_DATA\n')
            f.write('msout = %s\n' % output)
            f.write('msout.writefullresflag=False\n')
            f.write('steps=[avg]\n')
            f.write('rficonsole.type=aoflagger\n')
            f.write('avg.type = squash\n')
            f.write('avg.freqstep = 1\n')
            f.write('avg.timestep = 12\n')      # is the default
            f.close()


            ndpppcmd = 'NDPPP ' + ndppp_parset+ ' >'+ndppplog+' 2>&1'
            b.run(ndpppcmd)
            #os.system (ndpppcmd)

        b.wait()

    ms = ''
    for m in mslist:
        if average:
            ms = ms + ' ' + m + '.tmpavg'
        else:
            ms = ms + ' ' + m

    imout = 'im'+ callnumber +'_cluster'+cluster+'nm'
    logging.debug(ms + ' ' + imout + ' ' + 'None' + ' ' + '1mJy' + ' ' + str(niter) + ' ' + str(nterms))

    if do_mask: ####
        if cluster == 'a2256': ## special case for a2256
            niter = niter*15 # clean very deep here

        run('casapy --nogui --logfile casapy-'+imout+'.log -c '+SCRIPTPATH+'/casapy_cleanv4.py ' + ms + ' ' + imout + ' ' + 'None' +\
                ' ' + '1mJy' + ' ' + str(niter) + ' ' + str(nterms) + ' ' + str(imsize) + ' ' + mscale)
        # make mask
        if nterms > 1:
            run('python '+SCRIPTPATH+'/makecleanmask.py --threshpix '+str(threshpix)+\
                      ' --threshisl '+str(threshisl) +' --atrous_do '+ str(atrous_do) +' '   +imout +'.image.tt0')
        else:
            run('python '+SCRIPTPATH+'/makecleanmask.py --threshpix '+str(threshpix)+\
                    ' --threshisl '+str(threshisl) +' --atrous_do '+ str(atrous_do) + ' '  + imout +'.image')

        # clean image with manual mask
        mask = imout+'.cleanmask'


    niter = 1000
    imout = 'im'+ callnumber +'_cluster'+cluster

    if region != 'empty' : ## special cases
        if region.startswith("region"):
            ####
            # Test to use only the region mask
            niter = niter*3
            run('casapy --nogui --logfile casapy-'+imout+'.log -c '+SCRIPTPATH+'/casapy_cleanv4.py '+ ms + ' ' + imout + ' ' +region + \
                    ' ' + '1mJy' + ' ' + str(niter) + ' ' + str(nterms) + ' ' + str(imsize) + ' ' + mscale)

        else:
            niter = niter*3
            run('casapy --nogui --logfile casapy-'+imout+'.log -c '+SCRIPTPATH+'/casapy_cleanv4.py '+ ms + ' ' + imout + ' ' + mask+','+region + \
                    ' ' + '1mJy' + ' ' + str(niter) + ' ' + str(nterms) + ' ' + str(imsize) + ' ' + mscale)
    else:
        run('casapy --nogui --logfile casapy-'+imout+'.log -c '+SCRIPTPATH+'/casapy_cleanv4.py '+ ms + ' ' + imout + ' ' + mask + \
                  ' ' + '1mJy' + ' ' + str(niter) + ' ' + str(nterms) + ' ' + str(imsize) + ' ' + mscale)

    # convert to FITS
    if nterms > 1:
        run('image2fits in=' + imout +'.image.tt0' + ' ' + 'out='+ imout + '.fits')
    else:
        run('image2fits in=' + imout +'.image'     + ' ' + 'out='+ imout + '.fits')

    if not os.path.exists(imout+'.fits'):
        raise Exception('Imaging Error: '+imout+'.fits does not exist' )


    if average:
        logging.debug('rm -rf {}'.format(ms))
        os.system('rm -rf ' + ms)

    return imout,mask


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



    #### CALIBRATE  BBS PHASE+AMP 2 ###
    # make model
    run(SCRIPTPATH+'/casapy2bbs.py -m '+ mask + ' ' +'-t ' + str(nterms)+ ' ' + imout+'.model ' +  imout+'.skymodel')
    if FFT:
        run('casapy --nogui -c '+SCRIPTPATH+'/ft_v2.py ' + msinputlist + ' ' + imout+'.model' \
                  + ' ' + str(nterms) + ' '+ str(wplanes))

    #parmdb keep from previous step
    skymodel = imout+'.skymodel'


    # reset the phases from instrument_amps0 to zero to prevent large phase corrections from incorrect AP solve
    if phasezero:
        inputparmdb  = parmdb +'_smoothed'
        outputparmdb = parmdb +'_smoothed_phasezero'
        for ms in mslist:
            run('python '+SCRIPTPATH+'/setphasezero.py ' + ms + ' ' + ms+'/'+inputparmdb +' ' + ms+'/'+outputparmdb)
    else:
        outputparmdb = parmdb +'_smoothed'


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

    ### MAKE IMAGE 4 ###
    logging.info('Make image 4')
    imout,mask = make_image(mslist, cluster, '4', 10, 10, nterms, atrous_do, imsize, region, SCRIPTPATH)


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
