import os
import os.path
import glob
import sys
import numpy
import pwd
from subprocess import Popen, PIPE
import time
import pyrap.images
import logging
import logging.config
from facet_utilities import run

username = pwd.getpwuid(os.getuid())[0]

def create_ndppp_parset(msin, msout):
    ndppp_parset = msin.split('.')[0] +'ndppp_lowresavg.parset'
    os.system('rm -f ' + ndppp_parset)

    f=open(ndppp_parset, 'w')
    f.write('msin ="%s"\n' % msin)
    f.write('msin.datacolumn = MODEL_DATA\n')
    f.write('msin.autoweight = false\n')
    f.write('msout ="%s"\n' % msout)
    f.write('msout.writefullresflag=False\n')
    f.write('steps = [avg1]\n')
    f.write('avg1.type = squash\n')
    f.write('avg1.freqstep = 5\n')
    f.write('avg1.timestep = 2\n')
    f.close()
    return ndppp_parset

def wsclean_wait():
    warned=False
      ######### check that no wsclean is running
    cmd = "ps -f -u " + username + " | grep wsclean | grep -v _wsclean_v2.py | grep -v wsclean.parset |grep -v grep |wc -l"
    while True: # "grep -v grep" to prevent counting the grep command
        output=numpy.int(Popen(cmd, shell=True, stdout=PIPE).communicate()[0])
        if output==0:
            break
        if not(warned):
            print 'A wsclean is already running, waiting for it to finish'
            warned=True
        time.sleep(2)

############ MAIN CODE ###########

if __name__=='__main__':

    if len(sys.argv)<2:
        raise Exception('Give the path to the setup code for the facet')

    # setup code must set SCRIPTPATH, wsclean, mslist, casaregion, and may
    # do anything else needed

    print 'Using',sys.argv[1],'as the setup code'
    execfile(sys.argv[1])
    print 'script path is',SCRIPTPATH

    print 'working on MS list',mslist
    try:
        tempdir
    except NameError:
        tempdir=None

    try:
        cleanup
    except NameError:
        cleanup=False

    try:
        logfile
    except NameError:
        logfile='subtract.log'

    ## Logger configuration
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
        fh = logging.FileHandler(logfile) 
        fh.setFormatter(formatter)
        logger.addHandler(fh)

    logging.info('Subtract starts')
    logging.debug('MS list is '+str(mslist))
        
    # parameters below should not need to be changed

    niterh = 40000
    niterl = 20000
    nterms = 1
    imsizeh= 6144
    imsizel= 4800
    cellh  = '7.5arcsec'
    celll  = '25arcsec'

    threshisl = 2.5
    threshpix = 5
    atrous_do = "True"

    for ms in mslist:
        # names of files that will be checked for if cleanup is False, or
        # deleted if cleanup is True
        imhigh = ms.split('.')[0] + 'highres'
        imlow  = ms.split('.')[0] + 'lowres'
        finalsky =  ms.split('.')[0] + '.skymodel'
        lr_bbslog=ms + '.lowressubbbslog'
        hr_bbslog=ms + '.highressubbbslog'

        if cleanup:
            os.system('rm -r '+imhigh+'*')
            os.system('rm -r '+imlow+'*')
            os.system('rm '+lr_bbslog+' '+hr_bbslog+' '+finalsky)

        if os.path.isfile(imhigh+'-image.fits'):
            logging.warning('High-resolution image exists, NOT remaking it')
        else:
    # ---------------------
    # image without mask
            cmd = wsclean + ' -reorder -name ' + imhigh + ' -size ' + str(imsizeh) + ' ' + str(imsizeh) + ' '
            if tempdir is not None:
                cmd+='-tempdir '+tempdir+' '
            cmd+= '-scale ' + cellh + ' -weight briggs 0.0 -niter ' + str(niterh) + ' '
            cmd+= '-maxuv-l 7e3 -mgain 0.65 -fitbeam -datacolumn CORRECTED_DATA -no-update-model-required ' + ms

            wsclean_wait()
            run(cmd)


    # create the mask
        mask_name  = imhigh + '.fitsmask'
        casa_mask  = imhigh + '.casamask'

        if os.path.isfile(mask_name):
            logging.warning('Mask exists, NOT remaking it')
        else:

            cmd='python '+SCRIPTPATH+'/makecleanmask_10sb_wsclean.py --threshpix '+str(threshpix)+\
                    ' --threshisl '+str(threshisl) +' --atrous_do '+ str(atrous_do)+' '
            if casaregion!='':
                cmd+='--casaregion  '+ casaregion +' '
            cmd+=imhigh + '-image.fits'

            run(cmd)

        if os.path.isdir(casa_mask):
            logging.warning('CASA mask exists, NOT remaking it')
        else:
            maskim=pyrap.images.image(mask_name)
            maskim.saveas(casa_mask)

        # include region file
            if casaregion != '':
                run('casapy --nogui -c '+SCRIPTPATH+'/fitsandregion2image.py '\
                     + mask_name + ' ' + casa_mask + ' ' + casaregion)
            else:
                run('casapy --nogui -c '+SCRIPTPATH+'/fitsandregion2image.py '\
                     + mask_name + ' ' + casa_mask + ' ' + 'None')

        imhigh = imhigh + 'withmask'

        if os.path.isfile(imhigh+'-image.fits'):
            logging.warning('Masked high-resolution image exists, NOT remaking it')
        else:
            cmd = wsclean + ' -reorder -name ' + imhigh + ' -size ' + str(imsizeh) + ' ' + str(imsizeh) + ' '
            if tempdir is not None:
                cmd+='-tempdir '+tempdir+' '
            cmd+= '-scale ' + cellh + ' -weight briggs 0.0 -niter ' + str(niterh) + ' '
            cmd+= '-maxuv-l 7e3 -no-update-model-required -mgain 0.65 -fitbeam -datacolumn CORRECTED_DATA -casamask ' + casa_mask  + ' ' + ms

            wsclean_wait()
            run(cmd)

        fits_model = imhigh + '-model.fits'
        casa_model = imhigh + '.model'
        if os.path.isdir(casa_model):
            logging.warning('CASA model exists, NOT remaking it')
        else:

        # convert model fits image to casa .model format
            logging.debug('Converting model to casapy format: '+fits_model+' ==> '+casa_model)

            modelim=pyrap.images.image(fits_model)
            modelim.saveas(casa_model)

        #run('casapy --nogui -c '+SCRIPTPATH+'/fits2image.py '\
        #            + fits_model + ' ' + casa_model)

        # ---------------------
        # make the skymodel
        skymodel = imhigh  + '.skymodel'
        if os.path.isfile(skymodel):
            logging.warning('Skymodel exists, NOT remaking it')
        else:
            run(SCRIPTPATH+'/casapy2bbs_one_patch_per_cc.py '  + casa_model + ' ' +  skymodel)

        if os.path.isfile(hr_bbslog):
            logging.warning('BBS log for subtraction exists, NOT re-doing subtraction. Be careful!')
        else:
        # ---------------------
        # subtract the cc
            parset = SCRIPTPATH+'/subtractall_highres_wsclean.parset'
            cmd = 'calibrate-stand-alone --replace-sourcedb --parmdb-name instrument_ap_smoothed '
            cmd+= ms + ' ' + parset + ' ' + skymodel + ' >' + hr_bbslog

            run(cmd)

        # ---------------------
        # average down for lowres image

        msout = ms.split('.')[0] + '.lowres.ms'
        if os.path.isdir(msout):
            logging.warning('Averaged lowres MS already exists, NOT remaking it')
        else:

            ndppp_parset = create_ndppp_parset(ms, msout)

            os.system('rm -rf ' + msout)
            run('NDPPP ' + ndppp_parset +' >'+ ms+'.NDPPPavelog')

        if os.path.isfile(imlow+'-image.fits'):
            logging.warning('Low-resolution image exists, NOT remaking it')
        else:

        # make the lowres image (no mask)
            cmd = wsclean + ' -reorder -name ' + imlow + ' -size ' + str(imsizel) + ' ' + str(imsizel) + ' '
            if tempdir is not None:
                cmd+='-tempdir '+tempdir+' '
            cmd+= '-scale ' + celll + ' -weight briggs 0.0 -niter ' + str(niterl) + ' '
            cmd+= '-maxuv-l 2e3 -no-update-model-required -mgain 0.65 -fitbeam -datacolumn DATA ' + msout

        #########
            run(cmd)

        # ---------------------
        # create the lowres mask
        mask_name  = imlow + '.fitsmask'
        casa_mask  = imlow + '.casamask'

        if os.path.isdir(casa_mask):
            logging.warning('Low-res CASA mask exists, NOT re-making it')
        else:
            cmd='python '+SCRIPTPATH+'/makecleanmask_10sb_wsclean.py --threshpix '+str(5.0)+\
                   ' --threshisl '+str(4.0) +' --atrous_do '+ str(atrous_do) + ' '
            if casaregion!='':
                cmd+='--casaregion  '+ casaregion +' '
            cmd+=imlow + '-image.fits'
            run(cmd)

        # convert to CASA format

            maskim=pyrap.images.image(mask_name)
            maskim.saveas(casa_mask)

        # include region file
            if casaregion != '':
                run('casapy --nogui -c '+SCRIPTPATH+'/fitsandregion2image.py '\
                     + mask_name + ' ' + casa_mask + ' ' + casaregion)
            else:
                run('casapy --nogui -c '+SCRIPTPATH+'/fitsandregion2image.py '\
                     + mask_name + ' ' + casa_mask + ' ' + 'None')

        imlow = imlow + 'withmask'
        if os.path.isfile(imlow+'-image.fits'):
            logging.warning('Low-res masked image exists, NOT re-making it')
        else:
            cmd = wsclean + ' -reorder -name ' + imlow + ' -size ' + str(imsizel) + ' ' + str(imsizel) + ' '
            if tempdir is not None:
                cmd+='-tempdir '+tempdir+' '
            cmd+= '-scale ' + celll + ' -weight briggs 0.0 -niter ' + str(niterl) + ' '
            cmd+= '-maxuv-l 2e3 -no-update-model-required -mgain 0.65 -fitbeam -datacolumn DATA -casamask ' + casa_mask  + ' ' + msout

        ######### check that no wsclean is running
            wsclean_wait()
        #########
            run(cmd)




        fits_model = imlow + '-model.fits'
        casa_model = imlow + '.model'

        if os.path.isdir(casa_model):
            logging.warning('Low-res CASA model exists, NOT re-making it')
        else:
        # convert model fits image to casa .model format
            logging.debug('Converting model to casapy format: '+fits_model+' ==> '+casa_model)

            modelim=pyrap.images.image(fits_model)
            modelim.saveas(casa_model)

        #run('casapy --nogui -c '+SCRIPTPATH+'/fits2image.py '\
        #            + fits_model + ' ' + casa_model)
        # ---------------------
        # create the low-res skymodel
        skymodel = imlow  + '.skymodel'

        if os.path.isfile(skymodel):
            logging.warning('Skymodel exists, NOT remaking it')
        else:
            run(SCRIPTPATH+'/casapy2bbs_one_patch_per_cc.py '  + casa_model  + ' ' +  skymodel)

        # ---------------------
        # subtract thelowres cc
        if os.path.isfile(lr_bbslog):
            logging.warning('BBS log exists, NOT doing low-res subtract (CAREFUL!)')
        else:
            parset = SCRIPTPATH+'/subtractall_lowres_wsclean.parset'
            cmd = 'calibrate-stand-alone --replace-sourcedb --parmdb-name instrument_ap_smoothed '
            cmd = cmd + ms + ' ' + parset + ' ' + skymodel + ' >' + lr_bbslog
            run(cmd)


        # ---------------------
        # merge the skymodels
        finalsky =  ms.split('.')[0] + '.skymodel'
        if os.path.isfile(finalsky):
            logging.warning('Final sky model exists, NOT replacing it (why did you run this script, exactly?)')
        else:
            run('cp ' + imhigh+ '.skymodel ' + finalsky)
            run("grep -v '#' "+  imlow+ ".skymodel > tmp.sky")
            run("cat tmp.sky >>" + finalsky)

            cmd = "sed -i 's/, , , /, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.5e8, [0.0]/g'"
            run(cmd + ' ' + finalsky)

            cmd = 'sed -i \"s/# (Name, Type, Patch, Ra, Dec, I, Q, U, V) = format/format = Name, Type, Patch, Ra, Dec, I, Q, U, V, MajorAxis, MinorAxis, Orientation, ReferenceFrequency=\'1.5e+08\', SpectralIndex=\'[]\'/g\"  '
            run(cmd + ' ' + finalsky)

        logging.info('MS '+ms+' all done!')
