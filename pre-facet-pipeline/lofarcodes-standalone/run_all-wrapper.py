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
import glob

htmltarget = str(raw_input('Enter html target '))
htmlcal = str(raw_input('Enter html cal '))
gsmfile = str(raw_input('Enter gsmfile '))
calibratemodel = str(raw_input('Enter calibrator model '))
chunksize = int(raw_input('Number of SBs to process together '))
cpus = int(raw_input('Number of CPUs '))
numberofsubbbands = int(raw_input('Number of subbands in html files '))
flagants = str(raw_input('Enter flagants '))
username = str(raw_input('Enter username '))
password = str(raw_input('Enter password '))

os.system('mkdir htmlfiles')
startdir = os.getcwd()
subrange = np.arange(0,numberofsubbbands,chunksize)
subrange = np.append(subrange,numberofsubbbands)
for i in range(0,len(subrange)-1):
    os.chdir(startdir)

    newhtmltarget = 'htmlfiles/TARGET_html_%s_to_%s.txt'%(subrange[i],subrange[i+1])
    newhtmlcal = 'htmlfiles/CAL_html_%s_to_%s.txt'%(subrange[i],subrange[i+1])
    split_htmlfile(htmltarget,subrange[i],subrange[i+1],newhtmltarget,numberofsubbbands)
    split_htmlfile(htmlcal,subrange[i],subrange[i+1],newhtmlcal,numberofsubbbands)
    if not os.path.exists(newhtmltarget) or not os.path.exists(newhtmlcal):
        continue
    targetdir = 'run_%s_to_%s'%(subrange[i],subrange[i+1])
    if os.path.exists(targetdir):
        continue
    os.system('mkdir run_%s_to_%s'%(subrange[i],subrange[i+1]))
    runallfile = open('runall.sh','w')
    runallfile.write('python %s/run_all.py <<EOF\n'%scriptdir)
    runallfile.write('%s/%s\n'%(startdir,newhtmltarget))
    runallfile.write('%s/%s\n'%(startdir,newhtmlcal))
    runallfile.write('%s\n'%calibratemodel)
    runallfile.write('%s\n'%gsmfile)
    runallfile.write('%s\n'%cpus)
    runallfile.write('%s\n'%flagants)
    runallfile.write('%s\n'%username)
    runallfile.write('%s\n'%password)
    runallfile.write('<<EOF\n')
    runallfile.close()
    os.system('chmod +x runall.sh')
    os.system('mv runall.sh run_%s_to_%s'%(subrange[i],subrange[i+1]))
    os.system('cp image-torque-lofar-selfcal-new.sh  run_%s_to_%s'%(subrange[i],subrange[i+1]))
    os.system('cp %s  run_%s_to_%s'%(gsmfile,subrange[i],subrange[i+1]))

    #os.chdir('run_%s_to_%s'%(subrange[i],subrange[i+1]))
    #os.system('qsub image-torque-lofar-selfcal-new.sh')
    
    # Check download complete
    #currentdirs = glob.glob('start*/*')
    #while len(currentdirs)
