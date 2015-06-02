#!/usr/bin/env python

import os
import pwd
import numpy
import sys
from facet_utilities import bg

username = pwd.getpwuid(os.getuid())[0]

if len(sys.argv)<2:
    raise Exception('Give the path to the setup code for the facet')

print 'Using',sys.argv[1],'as the setup code'
execfile(sys.argv[1])

mslist = ["{name:s}_SB{b1:03d}-{b2:03d}.{res:s}.ms".format(name=NAME,res=RES,b1=b,b2=b+9) for b in BANDS]

b=bg(maxp=5)
# copy of last is stored in corrected_data
for ms in mslist:
    print "undoing "+ms
    if os.path.exists('dde.log'):
        with open('dde.log','a') as f:
            f.write("undoing "+ms+"\n")
    b.run("taql 'update " + ms + " set SUBTRACTED_DATA_ALL=CORRECTED_DATA'")

b.wait()
