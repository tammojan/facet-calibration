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
import lofar.bdsm as bdsm
import pyfits
import multiprocessing as mp
import argparse

print 'Combines all data in the data directory'
print 'Calibrates the data off the gsmmodel'

parser = argparse.ArgumentParser() 
parser = argparse.ArgumentParser() 
parser.add_argument('ppn', help="Number of processors")
parser.add_argument('datadirectory',help="Directory that contains the data (and only the data)")
parser.add_argument('timeavg', help="Time averaging")
parser.add_argument('freqavg', help="Frequency averaging")


args = parser.parse_args()
ppn = int(args.ppn)
datadirectory = args.datadirectory
timeavg = int(args.timeavg)
freqavg = int(args.freqavg)

startdir = os.getcwd()
scriptdir = os.path.dirname(os.path.realpath(__file__))

combine_measurementsets('%s/*'%datadirectory,'mergedfile.ms','ndppp_comb.parset.log',timeavg,freqavg,'DATA')

flag_measurementset('mergedfile.ms','mergedfile.ms','ndppp_flag.parset.log')

infodict = find_fileinfo('mergedfile.ms','fileinfo.out')

newfiles = split_equal_chunks('mergedfile.ms','./split',infodict['inttime'],ppn,infodict['Obsdate_start'],infodict['Obsdate_end'])
for newfile in newfiles:
    process=mp.Process(target=phase_calibratedata,args=(newfile,'gsmmodel','gsmcal.parset',150,999999,4,1))
    process.start()
    while len(mp.active_children()) >= ppn:
        print len(mp.active_children()),'still running phase calibration of the field'
        time.sleep(10)
while len(mp.active_children()) >= 1:
    print len(mp.active_children()),'still running phase calibration of the field'
    time.sleep(10)
print 'NO processes running'
os.system('rm -rf mergedfile.ms')
new_lofar_imagemsfile = 'mergedfile.ms'
combinelist = ' '.join(newfiles)
os.system('python %s/aux-scripts/concat.py %s %s'%(scriptdir,new_lofar_imagemsfile,combinelist))
concat_parmdbs_copy(newfiles,new_lofar_imagemsfile,'instrument','instrument')
os.system('rm -rf split')
print 'FINISHED'
