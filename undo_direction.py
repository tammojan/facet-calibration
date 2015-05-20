#!/usr/bin/env python

import os
import time
from subprocess import Popen, PIPE
import pwd
import numpy
import sys
import glob

username = pwd.getpwuid(os.getuid())[0]

#mslist = glob.glob('BOOTES24_SB*.2ch8s.ms')

if len(sys.argv)<2:
    raise Exception('Give the path to the setup code for the facet')

print 'Using',sys.argv[1],'as the setup code'
execfile(sys.argv[1])

mslist = ["{name:s}_SB{b1:03d}-{b2:03d}.{res:s}.ms".format(name=NAME,res=RES,b1=b,b2=b+9) for b in BANDS]


# copy of last is stored in corrected_data
for ms in mslist:
    print "undoing "+ms
    if os.path.exists('dde.log'):
        with open('dde.log','a') as f:
            f.write("undoing "+ms+"\n")
    os.system("taql 'update " + ms + " set SUBTRACTED_DATA_ALL=CORRECTED_DATA' &")
    time.sleep(10)
    cmd = "ps -u " + username + " | grep taql | wc -l"
    output=numpy.int(Popen(cmd, shell=True, stdout=PIPE).communicate()[0])
    while output > 5 :   # dont run more than 5 at a time!! taql is disk intensive
        time.sleep(10)
        output=numpy.int(Popen(cmd, shell=True, stdout=PIPE).communicate()[0])
