#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys
import pyrap.tables as pt
import os
import lofar.parmdb
import numpy as numpy
import math
#import lofar.expion.parmdbmain
import scipy
import scipy.signal
#import matplotlib.pyplot as plt


def median_smooth(ampl, half_window):

    ampl_tot_copy = numpy.copy(ampl)
    ndata = len(ampl)
    flags = numpy.zeros(ndata, dtype=bool)
    sol = numpy.zeros(ndata + 2 * half_window)
    sol[half_window:half_window + ndata] = ampl

    for i in range(0, half_window):

        # Mirror at left edge.

        idx = min(ndata - 1, half_window - i)
        sol[i] = ampl[idx]

        # Mirror at right edge

        idx = max(0, ndata - 2 - i)
        sol[ndata + half_window + i] = ampl[idx]

    # fix oct 2012

    median_array = scipy.signal.medfilt(sol, half_window * 2. - 1)
    ampl_tot_copy = median_array[half_window:ndata + half_window]
    return ampl_tot_copy


def median_window_filter(ampl, half_window, threshold):

    ampl_tot_copy = numpy.copy(ampl)
    ndata = len(ampl)
    flags = numpy.zeros(ndata, dtype=bool)
    sol = numpy.zeros(ndata+2*half_window)
    sol[half_window:half_window+ndata] = ampl

    for i in range(0, half_window):
        # Mirror at left edge.
        idx = min(ndata-1, half_window-i)
        sol[i] = ampl[idx]

        # Mirror at right edge
        idx = max(0, ndata-2-i)
        sol[ndata+half_window+i] = ampl[idx]


    #fix oct 2012
    median_array  = scipy.signal.medfilt(sol,half_window*2-1)

    sol_flag = numpy.zeros(ndata+2*half_window, dtype=bool)
    sol_flag_val = numpy.zeros(ndata+2*half_window, dtype=bool)

    for i in range(half_window, half_window + ndata):
        # Compute median of the absolute distance to the median.
        window = sol[i-half_window:i+half_window+1]
        window_flag = sol_flag[i-half_window:i+half_window+1]
        window_masked = window[~window_flag]

        if len(window_masked) < math.sqrt(len(window)):
            # Not enough data to get accurate statistics.
            continue

        median = numpy.median(window_masked)
        q = 1.4826 * numpy.median(numpy.abs(window_masked - median))

        # Flag sample if it is more than 1.4826 * threshold * the
        # median distance away from the median.
        if abs(sol[i] - median) > (threshold * q):
            sol_flag[i] = True

        idx = numpy.where(sol == 0.0) # to remove 1.0 amplitudes
        #print idx
        #print 'sol', sol
        sol[idx] = True

    mask = sol_flag[half_window:half_window + ndata]

    for i in range(len(mask)):
        if mask[i]:
           ampl_tot_copy[i] = median_array[half_window+i] # fixed 2012
    return ampl_tot_copy


msname                   = str(sys.argv[1])                     
instrument_name          = str(sys.argv[2])
instrument_name_smoothed = str(sys.argv[3])  # msname +'.instrument_smoothed'

##### EDIT THESE PARAMETERS BELOW #####


pol_list = ['0:0','1:1']
gain     = 'Gain'


# I suggest to leave the above untouched !!

output_median = True  # if True the will be replaced with the median, using half_window2

             # if False than only the ampltudes will be filtered for outliers

output_phasezero = True  # if True the phases will be set to zero (for amplitude solutions transfer

             # if False the phases will be left untouched
#######################################

pdb = lofar.parmdb.parmdb(instrument_name)
parms = pdb.getValuesGrid('*')

key_names = parms.keys()
print key_names

anttab     = pt.table(msname + '/ANTENNA')
antenna_list    = anttab.getcol('NAME')
anttab.close()

print 'Stations available:', antenna_list
window = 4

for pol in pol_list:
    for antenna in antenna_list:
        print 'smoothing [antenna, polarization]:', antenna, pol

        #amp = numpy.copy(parms[gain + ':' + pol + ':Ampl:'+ antenna]['values'][:, 0])
   
        real = numpy.copy(parms[gain + ':' + pol + ':Real:'+ antenna]['values'][:, 0])
        imag = numpy.copy(parms[gain + ':' + pol + ':Imag:'+ antenna]['values'][:, 0])
     
        phase = numpy.arctan2(imag,real)
        amp   = numpy.sqrt(imag**2 + real**2)

        
        window_window = numpy.int(len(amp)/3.)

        amp = numpy.log10(amp)
       
        amp = median_window_filter(amp,window,6)  
        amp = median_window_filter(amp,window,6)         
        amp = median_window_filter(amp,7,6) # window of 7
        amp = median_window_filter(amp,4,6) # window of 4
        amp = median_window_filter(amp,3,6) # window of 3
        
        amp = 10**amp

        parms[gain + ':' + pol + ':Real:'+ antenna]['values'][:, 0] = amp*numpy.cos(phase)
        parms[gain + ':' + pol + ':Imag:'+ antenna]['values'][:, 0] = amp*numpy.sin(phase)

print 'writing the new database:', instrument_name_smoothed
print 'check your results with: parmdbplot.py', instrument_name_smoothed
print 'compare with: parmdbplot.py', instrument_name

pdbnew = lofar.parmdb.parmdb(instrument_name_smoothed, create=True)
pdbnew.addValues(parms)
pdbnew.flush()
