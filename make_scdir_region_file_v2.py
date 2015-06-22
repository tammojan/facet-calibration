#!/usr/bin/env python

import matplotlib
#matplotlib.use('GTK')
import numpy
import os
import sys
#import lofar.parmdb
#from scipy import interpolate
#import time
#import subprocess
from subprocess import Popen, PIPE, STDOUT
#import pyrap.tables as pt
#import pyrap.images
import pwd
#import logging
#logging.basicConfig(filename='dde.log',level=logging.DEBUG, format='%(asctime)s -  %(message)s', datefmt='%Y-%d-%m %H:%M:%S')



pi       = numpy.pi
username = pwd.getpwuid(os.getuid())[0]

if len(sys.argv)<2:
    raise Exception('Give the path to the setup code')

print 'Using',sys.argv[1],'as the setup code'
execfile(sys.argv[1])
print 'script path is',SCRIPTPATH

def write_ds9_allregions(regfile, source_info_rec, col='yellow'):
    s= '''global color={col} dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1
fk5
'''.format(col=col)

    with open(regfile,'w') as f:
        f.write(s)
        for i in range(len(source_info_rec)):

            ddir = source_info_rec['directions'][i].split(',')
            size = source_info_rec['imsizes'][i]
            name = source_info_rec['sourcelist'][i]

            sra = ddir[0].replace('h',':').replace('m',':')
            sdec = ddir[1].replace('d',':').replace('m',':')


            ssize = size*1.5   ## boxsize in arcsec
            ssize2 = 0.8*ssize   ## boxsize in arcsec


            #box(14:26:29.673,+35:46:03.02,1200",600",0.00143176)

            s = 'box({ra},{dec},{s}",{s}",0)'.format(ra=sra, dec=sdec, s=ssize)
            s += r'# text={'
            s+= '{name}'.format(name=name)
            s+= r'} '+'\n'
            f.write(s)
            s = 'box({ra},{dec},{s}",{s}",0)\n'.format(ra=sra, dec=sdec, s=ssize2)
            f.write(s)
    return


try:
    region_file_name
except NameError:
    region_file_name='scdir_region_out.reg'


########## END USER INPUT ########

source_info_rec = numpy.genfromtxt(peelsourceinfo, \
                                   dtype="S50,S25,S5,S5,i8,i8,i8,i8,S2,S255,S255,S255,S5", \
                                   names=["sourcelist","directions","atrous_do","mscale_field","imsizes",\
                                   "cellsizetime_p","cellsizetime_a","fieldsize","dynamicrange",\
                                   "regionselfc","regionfield","peelskymodel","outliersource"],comments='#')
write_ds9_allregions(region_file_name, source_info_rec)
#write_ds9_allregions(SCRIPTPATH+"/dde_weeren/bootes_hba/facet_descrip/peel_source_info_v4_disp.reg", source_info_rec,col='black')
