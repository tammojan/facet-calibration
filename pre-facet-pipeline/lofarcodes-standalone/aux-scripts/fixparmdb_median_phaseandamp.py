#################################################################################
#                                                                               #
# Written by Wendy Williams - 16 October 2014                                   #
#                                                                               #
# This script is intended to get rid of time dependence in gain solutions       #
# so they can be applied to another field                                       #
#                                                                               #
#       NOTE: THIS MEDIANS THE PHASES AND AMPLITUDES                            #
#                                                                               #
#################################################################################

# modified by W. Williams, 15 Jan 2014
# take median value instead of average and print out some values
# and make a copy of the parmdb

import lofar.parmdb as lp
import numpy
import sys
import os
import scipy.signal

if len(sys.argv) < 2:
    print 'Please give a parmdb name.'
    print 'Usage: python fixparmdb_zero_median.py <parmdbname>'

filename = sys.argv[1]
outfilename = sys.argv[2]

if os.path.isdir(outfilename):
    print 'output file exists'
    #sys.exit()
    os.system("mv {s} {s}.bak".format(s=outfilename))
os.system('cp -r %s %s' %(filename, outfilename))
filename = outfilename

parmdbmtable = lp.parmdb(filename)

dictionary = parmdbmtable.getValuesGrid('*')

real_names = parmdbmtable.getNames('*Real*')
imaginary_names = parmdbmtable.getNames('*Imag*')

names = parmdbmtable.getNames()

#import matplotlib.pyplot as pl

import numpy

def smooth(x,window_len=11,window='hanning'):
    """smooth the data using a window with requested size.
    
    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.
    
    input:
        x: the input signal 
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal
        
    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)
    
    see also: 
    
    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter
 
    TODO: the window parameter could be the window itself if an array instead of a string
    NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
    """

    if x.ndim != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."

    if x.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."


    if window_len<3:
        return x


    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"


    s=numpy.r_[x[window_len-1:0:-1],x,x[-1:-window_len:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=numpy.ones(window_len,'d')
    else:
        w=eval('numpy.'+window+'(window_len)')

    y=numpy.convolve(w/w.sum(),s,mode='valid')
    return y

def smooth_array(x, axis=0,window_len=11,window='hanning'):
    
    nx, ny = x.shape
    y = numpy.zeros((nx,ny))
    d = window_len/2
    for i in range(ny):
        xi = x[:,i]
        yi = smooth(xi,window_len=11,window='hanning')
        yi = yi[d:-d]
        y[:,i] = yi
    
    return y

#import matplotlib.pyplot as plt
for name in real_names:
    imaginary_name = name.replace('Real','Imag')
    real_values = dictionary[name]['values']
    imaginary_values = dictionary[imaginary_name]['values']
    #plt.figure()
    
    ntimes, nfreqs = real_values.shape
    vals = real_values +1.j*imaginary_values

    phases = numpy.angle(vals)
    amps = numpy.abs(vals)
    
    medamp = numpy.median(amps, axis=0)  ## for each frequency
    #smoothphases = smooth_array(phases)
    #smoothphase = scipy.signal.medfilt(phases.flatten(), 5 * 2 - 1)
    smedphase = phases[-1,:]  # take the last phase value... avoid edge
    
    #plt.plot(phases)
    #plt.plot(smoothphases)
    
    new_real_values = medamp*numpy.cos(smedphase)
    new_imaginary_values = medamp*numpy.sin(smedphase)
    parmdbmtable.addDefValues(name,new_real_values)
    parmdbmtable.addDefValues(imaginary_name,new_imaginary_values)
    print name, ": ", new_real_values, new_imaginary_values
    
    #if 'RS' in name:
        #pl.figure()
        #pl.plot(phases)
        #pl.plot(smoothphase)
        #pl.plot(len(smoothphase)-5, phase,'ko')
        
    #pl.show()
parmdbmtable.deleteValues('*')
parmdbmtable.flush()

