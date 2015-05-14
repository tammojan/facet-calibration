import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pylab as pylab
import pyrap.tables as pt
import pyrap.images as pim
import os,sys
import time
scriptdir = os.path.dirname(os.path.realpath(__file__)+'/aux-lib')
sys.path.append(scriptdir)
scriptdir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(scriptdir)
from lofar_general import *
from pythoncodes.astronomy import *
from pythoncodes.inputs import *
from pythoncodes.maths import *
import lofar.bdsm as bdsm
import pyfits
import multiprocessing as mp

print 'Input .MS should have amplitude calibrated data in DATA and phase calibrated (off e.g. the GSM) in CORRECTED_DATA'

orig_lofar_msfile = str(raw_input('Give LOFAR msfile '))
workingdir = str(raw_input('Enter working directory '))
new_lofar_msfile = 'datafile.ms'
new_lofar_imagemsfile = 'datafile.ms'
doselfcal = yes_no('Do self calibration ')
imagingsize = float(raw_input('Enter imaging size initial (degrees) '))
highres_lowres = str(raw_input('Calibrate to lowist res (25'' -- big field) or highish res (8'' -- small field) [low/high]'))
solsmooth=int(raw_input('Enter amplitude solution smoothing '))
ncpus = int(raw_input('Enter number of cpus '))
deeper_clean = str(raw_input('Do a deeper clean on the faint sources [y/n]  '))
doampcal = str(raw_input('Do amplitude calibration [y/n] '))
brightestsource = float(raw_input('Enter brightest source '))
robust = float(raw_input('Enter imaging robust '))

startdir = os.getcwd()

if os.path.exists('%s'%workingdir):
    print '%s exists -- exiting '%workingdir
    sys.exit(0)

os.system('mkdir %s'%workingdir)

# Copy the dataset to the desired place
os.system('cp -r %s %s/%s'%(orig_lofar_msfile,workingdir,new_lofar_msfile))
os.chdir(workingdir)

# Set up logging file
loggingfilename = 'lofar-image-logging.txt'
loggingfile = open('lofar-image-logging.txt','w')
starttime = time.time()
loggingfile.write('Start time %s \n'%(starttime))
loggingfile.close()


infodict = find_fileinfo(new_lofar_msfile,'file-info.txt')

uvlimits = 10**(np.histogram(np.log10(infodict['Baseline distances']),8)[1])
wavelength = 3E8/infodict['Ref Freq']

#uvlimits = filter(lambda x: x > 4000.0*wavelength, uvlimits)
# What we want to do here is to find antennas not within the current uvlimits for the shortest ones and then gradually increase the number of antennas.
#print 4000.0*wavelength
#for i in range(0,len(uvlimits)):
#    uvlimits[i] = uvlimits[i] + 100.0

if highres_lowres == 'high':
    uvlimits = np.array([1.0,1.5,2.0,3.0,5.0])*6.0*(1000.0*wavelength) # Start at 6klambda (~35arcsec) and go to 32klambda 
if highres_lowres == 'low':
    uvlimits = np.array([1.0,1.5,2.0,4.0])*6.0*(1000.0*wavelength) # Start at 6klambda (~35arcsec) and go to 24klambda (~10arcsec)
    uvlimits = np.array([1.0,1.5,2.0])*6.0*(1000.0*wavelength) # Start at 6klambda (~35arcsec) and go to 12klambda (~20arcsec)
if highres_lowres == 'vlow':
    uvlimits = np.array([1.0])*6.0*(1000.0*wavelength) # Start at 6klambda (~35arcsec) and go to 12klambda (~20arcsec)

uvlimits.sort()




print 'UVlimits'
print uvlimits
os.system("echo 'Setup time %s \n' >> %s"%(time.time()-starttime,loggingfilename))
previoustime = time.time()


for i in range(0,len(uvlimits)):
    uvmin = 0.15
    uvmax = uvlimits[i]/(1000.0*wavelength)
    resolution = calc_resolution(infodict['Ref Freq'],uvlimits[i])*rad2arcsec
    cellsize = str(round(resolution/5.0,2)) + 'arcsec'
 
    #### Enough pixels so that bdsm will not say the image has uniform noise when using adaptive box sizes
    npix = int(imagingsize*deg2arcsec/(resolution/5.0) + 10.0) 
    if npix < 400:
        npix = 400
    npix = close_higher_bsmooth_number(int(npix),5)

    workingimage= 'IMAGE_%s_%s'%(round(resolution,2),i)

    if i == 0:
	iterations,threshold = 300000,'%smJy'%(np.max([infodict['Est Noise']*10.0,brightestsource*1000.0/400.0]))
        prevthreshold = np.max([infodict['Est Noise']*10.0,brightestsource*1000.0/400.0])
    else:
	iterations,threshold = 300000,'%smJy'%(np.max([rms*8.0*1000,brightestsource*1000.0/400.0])) # Do not clean so deeply without the clean boxes
        prevthreshold = np.max([rms*8.0*1000,brightestsource*1000.0/400.0])
        
    flag_measurementset(new_lofar_imagemsfile,new_lofar_imagemsfile,'CORRECTED_DATA','ndppp_flag_%s.log'%i)
        
    aw_imagedata('%s'%new_lofar_imagemsfile,'aw_imager_%s.inp'%i,'%s'%workingimage,iterations,threshold,cellsize,uvmin,uvmax,npix,'mfclark','',robust)
    os.system("echo 'Iter%s initial imaging time %s \n' >> %s"%(i,time.time()-previoustime,loggingfilename))
    previoustime = time.time()

    img = bdsm.process_image('%s.restored'%workingimage,advanced_opts='True',detection_image='%s.restored'%workingimage,thresh_isl=3,thresh_pix=5,blank_limit=1E-4,adaptive_rms_box='True',adaptive_thresh=200)#,adaptive_rms_box='True',atrous_do='True')
    img.export_image(outfile="mask_%s"%workingimage,img_type='island_mask',img_format='casa')
    img.export_image(outfile="rms_%s"%workingimage,img_type='rms',img_format='casa')
    os.system('image2fits in=%s.restored out=%s.restored.fits'%(workingimage,workingimage))
    f = pyfits.open('%s.restored.fits'%workingimage)
    noisearray = f[0].data.flatten()
    maxpixel = np.max(noisearray)
    noisearray = np.random.permutation(noisearray)[:10000]
    noisepix = np.array(filter(lambda x: abs(x) > 10E-8,noisearray))
    noisepix = np.array(filter(lambda x: abs(x)<infodict['Est Noise']*50.0/1000.0,noisepix))
    rms = fit_gaussian_histogram(noisepix,'n')
    print 'rms %s, maxpixel %s'%(rms,maxpixel)
    f.close()
    minthreshold = rms

    image = pim.image("rms_%s"%workingimage)
    imshape = image.shape()
    dataarray = image.getdata(blc=[0,0,0,0], trc=[imshape[0]-1,imshape[1]-1,imshape[2]-1,imshape[3]-1])
    dataarray[np.where(np.isnan(dataarray))] = 0
    dataarray = dataarray.flatten()
    threshold1 = str(np.min([np.max(dataarray)*1E3,prevthreshold]))+'mJy'
    dataarray = np.random.permutation(dataarray)[:10000]
    dataarray = np.array(filter(lambda x: abs(x) > minthreshold,dataarray))
    dataarray = np.random.permutation(dataarray)[:10000]
    dataarray.sort()

    threshold2 = dataarray[int(len(dataarray)*0.99)]
    threshold3 = dataarray[int(len(dataarray)*0.85)]
    threshold4 = dataarray[int(len(dataarray)*0.60)]
    thresholds = [str(threshold2*1E3)+'mJy',str(threshold3*1E3)+'mJy',str(threshold4*1E3)+'mJy']
    threshold2 = thresholds[0]
    threshold3 = thresholds[1]
    threshold4 = thresholds[2]
    print 'Thresholds', thresholds,threshold1

    os.system('mkdir orig-images')
    os.system('mv %s* orig-images/'%workingimage)
    os.system('mv "rms_%s" orig-images/'%workingimage)

    aw_imagedata('%s'%new_lofar_imagemsfile,'aw_imager_%s.inp'%i,'%s'%workingimage,iterations,threshold1,cellsize,uvmin,uvmax,npix,'mfclark','mask_%s'%workingimage,robust)
    os.system("echo 'Iter%s masked imaging time %s \n' >> %s"%(i,time.time()-previoustime,loggingfilename))
    previoustime = time.time()

    os.system('mv "rms_%s" orig-images/'%workingimage)
    img = bdsm.process_image('%s.restored'%workingimage,advanced_opts='True',detection_image='%s.restored'%workingimage,thresh_isl=3,thresh_pix=5,blank_limit=1E-4,adaptive_rms_box='True',adaptive_thresh=200)#,adaptive_rms_box='True',atrous_do='True')
    img.export_image(outfile="mask_%s"%workingimage,img_type='island_mask',img_format='casa')
    img.export_image(outfile="rms_%s"%workingimage,img_type='rms',img_format='casa')
    
    os.system('mkdir pre-deeper-images')
    os.system('cp -r %s* pre-deeper-images/'%workingimage)
    os.system('mv rms_%s pre-deeper-images/'%workingimage)
    os.system('cp -r mask_%s pre-deeper-images/'%workingimage)
    deeper_2nd_clean(workingimage,workingimage,new_lofar_imagemsfile,iterations,threshold2,cellsize,uvmin,uvmax,npix,'mfclark',threshold2,robust,200)
    os.system("echo 'Iter%s deeper masked imaging time %s \n' >> %s"%(i,time.time()-previoustime,loggingfilename))
    previoustime = time.time()

    os.system('mkdir pre-deeper-images2')
    os.system('cp -r %s* pre-deeper-images2/'%workingimage)
    os.system('mv rms_%s pre-deeper-images2/'%workingimage)
    os.system('cp -r mask_%s pre-deeper-images2/'%workingimage)
    deeper_2nd_clean(workingimage,workingimage,new_lofar_imagemsfile,iterations,threshold3,cellsize,uvmin,uvmax,npix,'mfclark',threshold3,robust,50)
    os.system("echo 'Iter%s deeper masked imaging 2 time %s \n' >> %s"%(i,time.time()-previoustime,loggingfilename))
    previoustime = time.time()
       
    os.system('mkdir pre-deeper-images3')
    os.system('cp -r %s* pre-deeper-images3/'%workingimage)
    os.system('mv rms_%s pre-deeper-images3/'%workingimage)
    os.system('cp -r mask_%s pre-deeper-images3/'%workingimage)
    deeper_2nd_clean(workingimage,workingimage,new_lofar_imagemsfile,iterations,threshold3,cellsize,uvmin,uvmax,npix,'mfclark',threshold4,robust,50)
    os.system("echo 'Iter%s deeper masked imaging 3 time %s \n' >> %s"%(i,time.time()-previoustime,loggingfilename))
    previoustime = time.time()
    #img = bdsm.process_image('%s.restored'%workingimage,adaptive_rms_box='True',advanced_opts='True',detection_image='%s.restored'%workingimage,thresh_isl=6,thresh_pix=8,blank_limit=1E-4)#,atrous_do='True')
    #img.export_image(outfile="mask_%s"%workingimage,img_type='island_mask',img_format='casa',clobber=True)
    #img.export_image(outfile="rms_%s"%workingimage,img_type='rms',img_format='casa',clobber=True)

    #os.system("echo 'Iter%s second source finding took %s \n' >> %s"%(i,time.time()-previoustime,loggingfilename))
    #previoustime = time.time()

    if doselfcal == 'y':
        os.system('casapy2bbs.py --mask=mask_%s %s.model.corr SELF%s.bbs'%(workingimage,workingimage,i))

        if i != len(uvlimits) -1:
            if i == 0:
                newfiles = split_equal_chunks(new_lofar_msfile,'./split',infodict['inttime'],ncpus,infodict['Obsdate_start'],infodict['Obsdate_end'])
            for newfile in newfiles:
                process=mp.Process(target=phase_calibratedata,args=(newfile,'SELF%s.bbs'%i,'BBSrun%s'%i,uvmin*1000.0,2*uvmax*1000.0,4,1))
                process.start()
                while len(mp.active_children()) >= ncpus:
                    print len(mp.active_children()),'still running phase calibration for iter',i,time.asctime(time.gmtime())
                    time.sleep(10)
            time.sleep(60) #Just to give a little break incase they are taking a min to start
            while len(mp.active_children()) >= 1:
                print len(mp.active_children()),'still running phase calibration for iter',i,time.asctime(time.gmtime())
                time.sleep(100)
	    print 'Finished phase calibration'
            new_lofar_imagemsfile = 'new_comb_%s.ms'%i
            combinelist = ' '.join(newfiles)
            os.system('python %s/aux-scripts/concat.py %s %s'%(scriptdir,new_lofar_imagemsfile,combinelist))
            concat_parmdbs_copy(newfiles,new_lofar_imagemsfile,'instrument','instrument')
            os.system("echo 'Iter%s phase calibration time %s \n' >> %s"%(i,time.time()-previoustime,loggingfilename))
            print 'Finished data combining', time.asctime(time.gmtime())
            previoustime = time.time()


        if i == len(uvlimits)-1 and doampcal == 'y':
            if i == 0:
                newfiles = split_equal_chunks(new_lofar_msfile,'./split',infodict['inttime'],ncpus,infodict['Obsdate_start'],infodict['Obsdate_end'])
            for newfile in newfiles:
                process=mp.Process(target=ampphase_calibratedata,args=(newfile,'SELF%s.bbs'%i,'BBSrun%s'%i,uvmin*1000.0,2*uvmax*1000.0,4,1))
                process.start()
                while len(mp.active_children()) >= ncpus:
                    print len(mp.active_children()),'still running amp phase calibration',time.asctime(time.gmtime())
                    time.sleep(100)
	    time.sleep(60) #Just a little break incase they take a while to start
            while len(mp.active_children()) >= 1:
                print len(mp.active_children()),'still running amp phase calibration',time.asctime(time.gmtime())
                time.sleep(10)
	    print 'Finished amp phase calibration'
            for newfile in newfiles:
                smooth_solutions(solsmooth,newfile,'T')
	    print 'Finished smoothing solution'
            new_lofar_imagemsfile = 'new_comb_amp_%s.ms'%i
            combinelist = ' '.join(newfiles)
            os.system('python %s/aux-scripts/concat.py %s %s'%(scriptdir,new_lofar_imagemsfile,combinelist))
            concat_parmdbs_copy(newfiles,new_lofar_imagemsfile,'instrument','instrument')
            os.system("echo 'Iter%s amplitude and phase calibration time %s \n' >> %s"%(i,time.time()-previoustime,loggingfilename))
            print 'Finished data combining', time.asctime(time.gmtime())
            previoustime = time.time()

            workingimage = 'FINAL_%s'%workingimage
            aw_imagedata('%s'%new_lofar_imagemsfile,'final_aw_imager_%s.inp'%i,'%s'%workingimage,iterations,threshold1,cellsize,uvmin,uvmax,npix,'mfclark','mask_%s'%workingimage,robust)
            os.system("echo 'Iter%s masked ampcal imaging time %s \n' >> %s"%(i,time.time()-previoustime,loggingfilename))
            previoustime = time.time()
    
            img = bdsm.process_image('%s.restored'%workingimage,advanced_opts='True',detection_image='%s.restored'%workingimage,thresh_isl=3,thresh_pix=5,blank_limit=1E-4,adaptive_rms_box='True',adaptive_thresh=200)#,adaptive_rms_box='True',atrous_do='True')
            img.export_image(outfile="mask_%s"%workingimage,img_type='island_mask',img_format='casa')
            img.export_image(outfile="rms_%s"%workingimage,img_type='rms',img_format='casa')

            os.system('cp -r %s* pre-deeper-images/'%workingimage)
            os.system('mv rms_%s pre-deeper-images/'%workingimage)
            os.system('cp -r mask_%s pre-deeper-images/'%workingimage)
            deeper_2nd_clean(workingimage,workingimage,new_lofar_imagemsfile,iterations,threshold2,cellsize,uvmin,uvmax,npix,'mfclark',threshold2,robust,200)
            os.system("echo 'Iter%s deeper masked ampcal imaging time %s \n' >> %s"%(i,time.time()-previoustime,loggingfilename))
            previoustime = time.time()

            os.system('cp -r %s* pre-deeper-images2/'%workingimage)
            os.system('mv rms_%s pre-deeper-images2/'%workingimage)
            os.system('cp -r mask_%s pre-deeper-images2/'%workingimage)
            deeper_2nd_clean(workingimage,workingimage,new_lofar_imagemsfile,iterations,threshold3,cellsize,uvmin,uvmax,npix,'mfclark',threshold3,robust,50)
            os.system("echo 'Iter%s deeper masked ampcal imaging 2 time %s \n' >> %s"%(i,time.time()-previoustime,loggingfilename))
            previoustime = time.time()

            os.system('cp -r %s* pre-deeper-images3/'%workingimage)
            os.system('mv rms_%s pre-deeper-images3/'%workingimage)
            os.system('cp -r mask_%s pre-deeper-images3/'%workingimage)
            deeper_2nd_clean(workingimage,workingimage,new_lofar_imagemsfile,iterations,threshold3,cellsize,uvmin,uvmax,npix,'mfclark',threshold4,robust,50)
            os.system("echo 'Iter%s deeper masked ampcal imaging 3 time %s \n' >> %s"%(i,time.time()-previoustime,loggingfilename))
            previoustime = time.time()



os.system('rm -r split')
loggingfile.write('End time %s \n'%(time.time()))
loggingfile.close()
