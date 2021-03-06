#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys
import os
import numpy as numpy
import math
import scipy
import scipy.signal
import pyrap.tables as pt
import lofar.parmdb


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


def my_solflag(ampl, half_window, threshold):

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

    sol_flag = numpy.zeros(ndata + 2 * half_window, dtype=bool)
    sol_flag_val = numpy.zeros(ndata + 2 * half_window, dtype=bool)

    for i in range(half_window, half_window + ndata):

        # Compute median of the absolute distance to the median.

        window = sol[i - half_window:i + half_window + 1]
        window_flag = sol_flag[i - half_window:i + half_window + 1]
        window_masked = window[~window_flag]

        if len(window_masked) < math.sqrt(len(window)):

            # Not enough data to get accurate statistics.

            continue

        median = numpy.median(window_masked)
        q = 1.4826 * numpy.median(numpy.abs(window_masked - median))

        # Flag sample if it is more than 1.4826 * threshold * the
        # median distance away from the median.

        if abs(sol[i] - median) > threshold * q:
            sol_flag[i] = True

    mask = sol_flag[half_window:half_window + ndata]

    for i in range(len(mask)):
        if mask[i]:
            ampl_tot_copy[i] = median_array[half_window + i]  # fixed 2012
    return ampl_tot_copy


def my_solflag_inv(ampl, half_window, threshold):

    ampl_tot_copy = 1. / numpy.copy(ampl)
    ndata = len(ampl)
    flags = numpy.zeros(ndata, dtype=bool)
    sol = numpy.zeros(ndata + 2 * half_window)
    sol[half_window:half_window + ndata] = 1. / ampl  # fixed 2012

    for i in range(0, half_window):

        # Mirror at left edge.

        idx = min(ndata - 1, half_window - i)
        sol[i] = 1. / ampl[idx]  # fixed 2012

        # Mirror at right edge

        idx = max(0, ndata - 2 - i)
        sol[ndata + half_window + i] = 1. / ampl[idx]  # fixed 2012

    # fix 2012 oct
    median_array = scipy.signal.medfilt(sol, half_window * 2. - 1)

    sol_flag = numpy.zeros(ndata + 2 * half_window, dtype=bool)
    sol_flag_val = numpy.zeros(ndata + 2 * half_window, dtype=bool)

    for i in range(half_window, half_window + ndata):

        # Compute median of the absolute distance to the median.

        window = sol[i - half_window:i + half_window + 1]
        window_flag = sol_flag[i - half_window:i + half_window + 1]
        window_masked = window[~window_flag]

        if len(window_masked) < math.sqrt(len(window)):

            # Not enough data to get accurate statistics.
            continue

        median = numpy.median(window_masked)
        q = 1.4826 * numpy.median(numpy.abs(window_masked - median))

        # Flag sample if it is more than 1.4826 * threshold * the
        # median distance away from the median.

        if abs(sol[i] - median) > threshold * q:
            sol_flag[i] = True

    mask = sol_flag[half_window:half_window + ndata]

    for i in range(len(mask)):
        if mask[i]:
            ampl_tot_copy[i] = median_array[half_window + i]  # fixed 2012
    return 1. / ampl_tot_copy

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
    median_array  = scipy.signal.medfilt(sol,half_window*2.-1)

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

    mask = sol_flag[half_window:half_window + ndata]

    for i in range(len(mask)):
        if mask[i]:
            ampl_tot_copy[i] = median_array[half_window+i] # fixed 2012
    return ampl_tot_copy




msname                   = str(sys.argv[1])
instrument_name          = str(sys.argv[2])
instrument_name_smoothed = str(sys.argv[3])  # msname +'.instrument_smoothed'

##### EDIT THESE PARAMETERS BELOW #####

gain             = 'Gain'

output_phasezero = True   # if True the phases will be set to zero (for amplitude solutions transfer
                          # if False the phases will be left untouched


#######################################

pdb = lofar.parmdb.parmdb(instrument_name)
parms = pdb.getValuesGrid('*')

key_names = parms.keys()
antenna_list = []
pol_list = []
sol_par = []
dir_list = []

# create the antenna+polarizations list
for ii in range(len(key_names)):

    string_a = str(key_names[ii])
    split_a = string_a.split(':')

    if gain == split_a[0]:
        print split_a
        antenna_list.append(split_a[4])
        pol_list.append(split_a[1] + ':' + split_a[2])
        sol_par.append(split_a[3])
        if gain == 'DirectionalGain':
            antenna_list[-1] = antenna_list[-1]+':'+split_a[5]

anttab          = pt.table(msname + '/ANTENNA')
antenna_list    = anttab.getcol('NAME')
anttab.close()

pol_list = numpy.unique(pol_list)
sol_par = numpy.unique(sol_par)
print 'Stations available:', antenna_list
print 'Polarizations:', pol_list, sol_par, gain

for pol in pol_list:
    for antenna in antenna_list:
        print 'setting to zero [antenna, polarization]:', antenna, pol
        real_val = parms[gain + ':' + pol + ':Real:'
                           + antenna]['values'][:, :]
        imag_val = parms[gain + ':' + pol + ':Imag:'
                           + antenna]['values'][:, :]

      # Check is we have freq dependent solutions

        print 'Found ', numpy.shape(real_val[0, :])[0], \
            'solution(s) along the frequency axis'

      # loop over channels

        for chan in range(numpy.shape(real_val[0, :])[0]):

            amp = numpy.sqrt((real_val[:,chan])**2 + (imag_val[:,chan])**2)
            phase = numpy.arctan2(imag_val[:,chan],real_val[:,chan])


        # ---------------------------------------------------------------
            if output_phasezero:
                parms[gain + ':' + pol + ':Imag:' + antenna]['values'][:,chan] = numpy.copy(amp*numpy.sin(0.0))
                parms[gain + ':' + pol + ':Real:' + antenna]['values'][:,chan] = numpy.copy(amp*numpy.cos(0.0))

            else:
                parms[gain + ':' + pol + ':Imag:' + antenna]['values'][:,chan] = numpy.copy(amp*numpy.sin(phase))
                parms[gain + ':' + pol + ':Real:' + antenna]['values'][:,chan] = numpy.copy(amp*numpy.cos(phase))



print 'writing the new database:', instrument_name_smoothed
print 'check your results with: parmdbplot.py', instrument_name_smoothed
print 'compare with: parmdbplot.py', instrument_name

pdbnew = lofar.parmdb.parmdb(instrument_name_smoothed, create=True)
pdbnew.addValues(parms)
pdbnew.flush()
