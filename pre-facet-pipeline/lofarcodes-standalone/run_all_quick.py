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
import glob
import multiprocessing as mp
import random





starttime = time.time()
randomid = random.random()
startdir = os.getcwd()
logfilename = '%s/run_all_%s.log'%(startdir,randomid)
loggingfile = open(logfilename,'w')
loggingfile.write('Start time %s \n'%starttime)
loggingfile.close()


downloadsourcelist = str(raw_input('Enter download list '))
downloadcallist = str(raw_input('Enter download list '))
calmodel = str(raw_input('Enter calibrator model '))
gsmmodel = str(raw_input('Enter the GSM model '))
processors = int(raw_input('Enter the number of processors '))
flagants = str(raw_input('Enter flagants '))
username = str(raw_input('Enter username '))
password = str(raw_input('Enter password '))


# Start monitoring in one of the cpus

processes1 = mp.Process(target=monitor_computer,args=('test',0.1,30))
processes1.start()


# Setting up the data structure
if os.path.exists('startid_%s'%randomid):
        print 'Directory startid_%s already exists -- exiting'%(randomid)
        sys.exit(0)
os.system('mkdir startid_%s'%randomid)
os.system('cp %s startid_%s'%(downloadsourcelist,randomid))
os.system('cp %s startid_%s'%(gsmmodel,randomid))
os.system('cp %s startid_%s'%(calmodel,randomid))
os.system('cp %s startid_%s'%(downloadcallist,randomid))
os.chdir('startid_%s'%randomid)
startdir = os.getcwd()
os.system('mkdir raw_files')

# Download the data
os.system('cp %s raw_files/'%(downloadsourcelist))
os.system('cp %s raw_files/'%(downloadcallist))
os.chdir('raw_files')
downloaddata(downloadcallist,'downloadcal.log',username,password)
os.system("echo 'Downloaded calibration files\n' >> %s"%logfilename)
os.system('rm -rf SRM*.tar')
calfiles = glob.glob('*.MS')
downloaddata(downloadsourcelist,'downloadtarget.log',username,password)
targetfiles = list(set(glob.glob('*.MS'))-set(calfiles))
os.system("echo 'Downloaded target files\n' >> %s"%logfilename)
os.system('rm -rf SRM*.tar')
downloadtime = time.time()-starttime
os.system("echo '@@@@@@@Download took %s seconds\n' >> %s"%(downloadtime,logfilename))


# Sort the data
fileids = {}
for calfile in calfiles:
    os.system("echo 'Downloaded CAL %s\n' >> %s"%(calfile,logfilename))
    fieldid = calfile.split('_') [0]
    os.system('mkdir ../%s'%fieldid)
    fileinfo = find_fileinfo(calfile,'%s.info'%calfile)
    os.system('mv %s ../%s'%(calfile,fieldid))
    os.system('mv %s.info ../%s'%(calfile,fieldid))
    try:
        fileids['Cal'].append([fileinfo['Ref Freq'],calfile,fieldid])
    except KeyError:
        fileids['Cal'] = [[fileinfo['Ref Freq'],calfile,fieldid]]

for targetfile in targetfiles:
    os.system("echo 'Downloaded TARGET %s\n' >> %s"%(targetfile,logfilename))
    fieldid = targetfile.split('_') [0]
    os.system('mkdir ../%s'%fieldid)
    fileinfo = find_fileinfo(targetfile,'%s.info'%targetfile)
    os.system('mv %s ../%s'%(targetfile,fieldid))
    os.system('mv %s.info ../%s'%(targetfile,fieldid))
    try:
        fileids['Target'].append([fileinfo['Ref Freq'],targetfile,fieldid])
    except KeyError:
        fileids['Target'] = [[fileinfo['Ref Freq'],targetfile,fieldid]]

# Match the calibrators and fields
for i in range(0,len(fileids['Target'])):
    for j in range (0,len(fileids['Cal'])):
        print fileids['Target'][i][0], fileids['Cal'][j][0]
        if fileids['Target'][i][0] == fileids['Cal'][j][0]:
            os.system("echo 'Matched %s and %s\n' >> %s"%(fileids['Target'][i],fileids['Cal'][j],logfilename))
            fileids['Target'][i].append(j)
sortingtime = time.time()-starttime-downloadtime
os.system("echo '@@@@@@@Sorting took %s seconds\n' >> %s"%(sortingtime,logfilename))



# Prepare the calibrators for calibration
for i in range(0,len(fileids['Cal'])):
    os.chdir(startdir)
    os.chdir(fileids['Cal'][i][2])
    calfile = fileids['Cal'][i][1]
    process = mp.Process(target=ndppp_prepare,args=(calfile,calfile.replace('dppp','avg_dppp'),2,2,flagants,'ndppp_prepare_%s'%calfile,'ndppp_prepare_%s.log'%calfile))
    process.start()
    while len(mp.active_children()) >= processors:
        print len(mp.active_children()),' still running preparing the calibration of the calibrators'
        time.sleep(10)

while len(mp.active_children()) > 0:
    print len(mp.active_children()),' still running preparing the calibration of the calibrators'
    time.sleep(10)

# Calibrate the calibrators
for i in range(0,len(fileids['Cal'])):
    os.chdir(startdir)
    os.system('cp %s %s/calmodel'%(calmodel,fileids['Cal'][i][2]))
    os.chdir(fileids['Cal'][i][2])
    calfile = fileids['Cal'][i][1].replace('dppp','avg_dppp')

    process = mp.Process(target=ampphase_calibratedata,args=(calfile,'calmodel','cal_%s.parset'%calfile,150,999999,1,1))
    process.start()
    while len(mp.active_children()) >= processors:
        print len(mp.active_children()),' still running calibrating the calibrators'
        time.sleep(10)

while len(mp.active_children()) > 0:
    print len(mp.active_children()),' still running calibrating the calibrators'
    time.sleep(10)

#Check calibration happened
for i in range(0,len(fileids['Cal'])):
    os.chdir(startdir)
    os.chdir(fileids['Cal'][i][2])
    calfile = fileids['Cal'][i][1].replace('dppp','avg_dppp')
    if not os.path.exists('%s/instrument'%calfile):
            print 'Calibration for calibrator %s failed -- exiting'%calfile
            sys.exit(0)
    
calibratortime = time.time()-starttime-sortingtime
os.system("echo '@@@@@@@Calibrating the calibrator took %s seconds\n' >> %s"%(calibratortime,logfilename))


# Calibrate the field with the calibrator
for i in range(0,len(fileids['Target'])):
    os.chdir(startdir)
    try:
	print fileids['Target'][i]
    except IndexError:
        print 'No match for element %s'%i
	print fileids['Target']         
	continue
    try:
        os.chdir(fileids['Target'][i][2])
    except IndexError:
        print fileids['Target'][i],' does not have a calibrator'
        del fileids['Target'][i]
        continue
    fieldfile = fileids['Target'][i][1]
    try:
        calid = fileids['Target'][i][3]
    except IndexError:
	print fileids['Target'][i], 'does not have a calibrator2'
        del fileids['Target'][i]
        continue

# Calibrate the field with the calibrator
for i in range(0,len(fileids['Target'])):
    os.chdir(startdir)
    os.chdir(fileids['Target'][i][2])
    calid = fileids['Target'][i][3]
    fieldfile = fileids['Target'][i][1]
    calfile = '../%s/%s'%(fileids['Cal'][calid][2],fileids['Cal'][calid][1])
    calfile = calfile.replace('dppp','avg_dppp')
    process = mp.Process(target=ndppp_prepare,args=(fieldfile,fieldfile.replace('dppp','avg_dppp'),2,2,flagants,'ndppp_prepare_%s'%fieldfile,'ndppp_prepare_%s.log'%fieldfile))
    process.start()
    while len(mp.active_children()) >= processors:
        print len(mp.active_children()),' still running preparing the field for calibration'
        time.sleep(10)
while len(mp.active_children()) > 0:
    print len(mp.active_children()),' still running  preparing the field for calibration'
    time.sleep(10)

for i in range(0,len(fileids['Target'])):
    os.chdir(startdir)
    os.chdir(fileids['Target'][i][2])
    fieldfile = fileids['Target'][i][1]
    calid = fileids['Target'][i][3]
    calfile = '../%s/%s'%(fileids['Cal'][calid][2],fileids['Cal'][calid][1])
    calfile = calfile.replace('dppp','avg_dppp')
    process = mp.Process(target=ateam_clip,args=('%s_ateam.parset'%fieldfile,fieldfile.replace('dppp','avg_dppp')))
    process.start()
    while len(mp.active_children()) >= processors:
        print len(mp.active_children()),' still running ateam clipping of the field'
        time.sleep(10)
while len(mp.active_children()) > 0:
    print len(mp.active_children()),' still running ateam clipping of the field'
    time.sleep(10)

#Check Ateam demixed
for i in range(0,len(fileids['Target'])):
    os.chdir(startdir)
    os.chdir(fileids['Target'][i][2])
    fieldfile = fileids['Target'][i][1]
    if not os.path.exists('%s/sky'%(fieldfile.replace('dppp','avg_dppp'))):
            print 'Ateam subtraction for field %s failed -- exiting'%fieldfile
            sys.exit(0)

    
for i in range(0,len(fileids['Target'])):
    os.chdir(startdir)
    os.chdir(fileids['Target'][i][2])
    fieldfile = fileids['Target'][i][1]
    calid = fileids['Target'][i][3]
    calfile = '../%s/%s'%(fileids['Cal'][calid][2],fileids['Cal'][calid][1])
    calfile = calfile.replace('dppp','avg_dppp')
    os.system('python %s/aux-scripts/fixparmdb_median_phaseandamp.py %s/instrument cal_instrument_%s'%(scriptdir,calfile,fileids['Cal'][calid][1]))
    process = mp.Process(target=bbs_correct,args=('%s_correct.parset'%fieldfile,fieldfile.replace('dppp','avg_dppp'),'cal_instrument_%s'%fileids['Cal'][calid][1]))
    process.start()
    while len(mp.active_children()) >= processors:
        print len(mp.active_children()),' still running flux calibrating the field'
        time.sleep(10)   
while len(mp.active_children()) > 0:
    print len(mp.active_children()),' still running flux calibrating the field'
    time.sleep(10)


for i in range(0,len(fileids['Target'])):
    os.chdir(startdir)
    os.chdir(fileids['Target'][i][2])
    fieldfile = fileids['Target'][i][1]
    calid = fileids['Target'][i][3]
    calfile = '../%s/%s'%(fileids['Cal'][calid][2],fileids['Cal'][calid][1])
    calfile = calfile.replace('dppp','avg_dppp')
    process = mp.Process(target=copy_colums,args=(fieldfile.replace('dppp','avg_dppp'),fieldfile.replace('dppp','avg_dppp_cols'),'ndppp_copycols_%s'%fieldfile.replace('dppp','avg_dppp'),'CORRECTED_DATA','DATA'))
    process.start()
    while len(mp.active_children()) >= processors:
        print len(mp.active_children()),' still running averaging the field'
        time.sleep(10)   
while len(mp.active_children()) > 0:
    print len(mp.active_children()),' still running averaging the field'
    time.sleep(10)

fluxcalibration = time.time()-starttime-calibratortime
os.system("echo '@@@@@@@Flux calibrating the field and running ateamclipper took %s seconds\n' >> %s"%(fluxcalibration,logfilename))

# Calibrate the field off the GSM and make simple images
os.chdir(startdir)
os.system('mkdir Imaging')
os.system('cp %s %s/gsmmodel'%(gsmmodel,'Imaging'))
os.chdir('Imaging')
os.system('mkdir tobe_comb')
for i in range(0,len(fileids['Target'])):
    fieldfile = '../%s/%s'%(fileids['Target'][i][2],fileids['Target'][i][1].replace('dppp','avg_dppp_cols'))
    os.system('mv %s tobe_comb'%fieldfile)
    delfile =  '../%s/%s'%(fileids['Target'][i][2],fileids['Target'][i][1])
    os.system('rm -rf %s'%delfile)
    delfile = '../%s/%s'%(fileids['Target'][i][2],fileids['Target'][i][1].replace('dppp','avg_dppp'))
    os.system('rm -rf %s'%delfile)
gsmcalfile=open('run_gsmcal.sh','w')
gsmcalfile.write('python %s/lofar-gsmcal.py %s %s %s %s\n'%(scriptdir,processors,'tobe_comb',2,2))
gsmcalfile.close()
os.system('chmod +x run_gsmcal.sh')
os.system('./run_gsmcal.sh')
gsmcaltime = time.time()-starttime-fluxcalibration
os.system("echo '@@@@@@@GSM calibrating the field took %s seconds\n' >> %s"%(gsmcaltime,logfilename))
os.system('python %s/aux-scripts/plot_solutions_all_stations.py mergedfile.ms/instrument -p sols'%scriptdir)

# Selfcal the data.
infodict = find_fileinfo('mergedfile.ms','infodict')
imagesize,brightestsource = find_imagesize('gsmmodel',infodict['RA(rad)'],infodict['DEC(rad)'])
imagesize = imagesize*2.0
selfcalfile=open('run_selfcal.sh','w')
selfcalfile.write('python %s/lofar-selfcal-image.py <<EOF \n'%scriptdir)
selfcalfile.write('mergedfile.ms\n')
selfcalfile.write('%s/selfcal\n'%os.getcwd())
selfcalfile.write('y\n')
selfcalfile.write('%s\n'%imagesize)
selfcalfile.write('vlow\n')
selfcalfile.write('30\n')	
selfcalfile.write('%s\n'%int(processors))
selfcalfile.write('y\n')
selfcalfile.write('n\n') # Dont bother doing the amp selfcal
selfcalfile.write('%s\n'%brightestsource)
selfcalfile.write('-0.5\n')#robust value
selfcalfile.write('<<EOF\n')
selfcalfile.close()
os.system('chmod +x run_selfcal.sh')
os.system('./run_selfcal.sh')
selfcaltime = time.time()-starttime-gsmcaltime
os.system("echo '@@@@@@@Selfcalibration of the field took %s seconds\n' >> %s"%(selfcaltime,logfilename))



# Finished the calibration
os.system("echo 'End time %s\n' >> %s"%(time.time(),logfilename))
os.system("echo '@@@@@@@Total time taken %s\n' >> %s"%(time.time()-starttime,logfilename))

