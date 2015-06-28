#!/usr/bin/python


# TO DO, need to fix input datasets with more than 1 freq. channels per MS

import numpy
import pyrap.tables as pt
import sys
import scipy.signal

msname         = str(sys.argv[1])
ionnormfactor  = numpy.float(sys.argv[2])
ionscalefactor = numpy.float(sys.argv[3])
nblockconcat   = numpy.int(sys.argv[4])
colname        = str(sys.argv[5])

t = pt.table(msname, readonly=False)


freq_tab   = pt.table(msname + '/SPECTRAL_WINDOW')
freq       = freq_tab.getcol('REF_FREQUENCY')
chanfreq   = freq_tab.getcol('CHAN_FREQ')[0]
wav        = 3e8/freq
anttab     = pt.table(msname + '/ANTENNA')
antlist    = anttab.getcol('NAME')

centerfreq = numpy.mean(chanfreq)
freq_res   = numpy.abs(chanfreq[0]-chanfreq[1])

if ((nblockconcat*2)+1) != len(chanfreq):
    raise Exception('More than 1 channel per block, needs to fixed.....')


for t2 in t.iter(["ANTENNA1","ANTENNA2"]):     
  if (t2.getcell('ANTENNA1',0)) < (t2.getcell('ANTENNA2',0)): 
    #print antlist[t2.getcell('ANTENNA1',0)],antlist[t2.getcell('ANTENNA2',0)] 
    weightscol = t2.getcol(colname)
    uvw        = t2.getcol('UVW')
    uvw_dist   = numpy.sqrt(uvw[:,0]**2 + uvw[:,1]**2 + uvw[:,2]**2)
    weightscol_modified = numpy.copy(weightscol)
    
    #timepersample    = t2[1]['TIME']-t2[0]['TIME']
    dist = numpy.mean(uvw_dist)/1e3
    stddev = ionnormfactor*((1./dist)**ionscalefactor)*(centerfreq/150e6)

    #print 'STEDDV', stddev

    fwhm = 2.3548*stddev/freq_res
    gauss = scipy.signal.gaussian(len(weightscol[0,:,0]),(stddev))
    #print gauss
    
    #sys.exit()
    for chan in range(0,len(weightscol[0,:,0])):
      weightscol_modified[:,chan,:] = weightscol[:,chan,:]*gauss[chan]
    
      #print numpy.mean(weightscol[:,chan,:]),numpy.mean(weightscol_modified[:,chan,:])
    

    t2.putcol(colname, weightscol_modified)

t.close()
freq_tab.close()
