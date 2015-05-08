import matplotlib
matplotlib.use('GTK')
import lofar.expion.ionosphere
import lofar.parmdb
import sys
import os
import scipy
import time
import numpy
import math
import pyrap.tables
import scipy.signal
import scipy.ndimage
import scipy.interpolate
from pylab import *
#import lofar.expion.fitting as fitting
pi = numpy.pi
from scipy.interpolate import splprep, splev
from scipy.interpolate import Rbf, InterpolatedUnivariateSpline

def rebin( a, newshape ):
 	'''Rebin an array to a new shape.
 	'''
        assert len(a.shape) == len(newshape)

        slices = [ slice(0,old, float(old)/new) for old,new in zip(a.shape,newshape) ]
        coordinates = mgrid[slices]
        indices = coordinates.astype('i')   #choose the biggest smaller integer index
        return a[tuple(indices)]

	
 
	 
	 
def rebin2d(a, new_shape):
    """
    Resizes a 2d array by averaging or repeating elements, 
    new dimensions must be integral factors of original dimensions

    Parameters
    ----------
    a : array_like
        Input array.
    new_shape : tuple of int
        Shape of the output array

    Returns
    -------
    rebinned_array : ndarray
        If the new shape is smaller of the input array, the data are averaged, 
        if the new shape is bigger array elements are repeated

    See Also
    --------
    resize : Return a new array with the specified shape.

    Examples
    --------
    >>> a = np.array([[0, 1], [2, 3]])
    >>> b = rebin(a, (4, 6)) #upsize
    >>> b
    array([[0, 0, 0, 1, 1, 1],
           [0, 0, 0, 1, 1, 1],
           [2, 2, 2, 3, 3, 3],
           [2, 2, 2, 3, 3, 3]])

    >>> c = rebin(b, (2, 3)) #downsize
    >>> c
    array([[ 0. ,  0.5,  1. ],
           [ 2. ,  2.5,  3. ]])

    """
    M, N = a.shape
    m, n = new_shape
    if m<M:
        return a.reshape((m,M/m,n,N/n)).mean(3).mean(1)
    else:
        return numpy.repeat(numpy.repeat(a, m/M, axis=0), n/N, axis=1)


def savitzky_golay(y, window_size, order, deriv=0, rate=1):
    r"""Smooth (and optionally differentiate) data with a Savitzky-Golay filter.
    The Savitzky-Golay filter removes high frequency noise from data.
    It has the advantage of preserving the original shape and
    features of the signal better than other types of filtering
    approaches, such as moving averages techniques.
    Parameters
    ----------
    y : array_like, shape (N,)
        the values of the time history of the signal.
    window_size : int
        the length of the window. Must be an odd integer number.
    order : int
        the order of the polynomial used in the filtering.
        Must be less then `window_size` - 1.
    deriv: int
        the order of the derivative to compute (default = 0 means only smoothing)
    Returns
    -------
    ys : ndarray, shape (N)
        the smoothed signal (or it's n-th derivative).
    Notes
    -----
    The Savitzky-Golay is a type of low-pass filter, particularly
    suited for smoothing noisy data. The main idea behind this
    approach is to make for each point a least-square fit with a
    polynomial of high order over a odd-sized window centered at
    the point.
    Examples
    --------
    t = np.linspace(-4, 4, 500)
    y = np.exp( -t**2 ) + np.random.normal(0, 0.05, t.shape)
    ysg = savitzky_golay(y, window_size=31, order=4)
    import matplotlib.pyplot as plt
    plt.plot(t, y, label='Noisy signal')
    plt.plot(t, np.exp(-t**2), 'k', lw=1.5, label='Original signal')
    plt.plot(t, ysg, 'r', label='Filtered signal')
    plt.legend()
    plt.show()
    References
    ----------
    .. [1] A. Savitzky, M. J. E. Golay, Smoothing and Differentiation of
       Data by Simplified Least Squares Procedures. Analytical
       Chemistry, 1964, 36 (8), pp 1627-1639.
    .. [2] Numerical Recipes 3rd Edition: The Art of Scientific Computing
       W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery
       Cambridge University Press ISBN-13: 9780521880688
    """
    import numpy as np
    from math import factorial
    
    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
    except ValueError, msg:
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order+1)
    half_window = (window_size -1) // 2
    # precompute coefficients
    b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve( m[::-1], y, mode='valid')
    
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
    
    
def running_median(ampl,half_window) :

    ampl_tot_copy = numpy.copy(ampl)

    ndata = len(ampl)
    flags = numpy.zeros(ndata, dtype=bool)
    sol = numpy.zeros(ndata+2*half_window)
    sol[half_window:half_window+ndata] = ampl
    std = numpy.zeros(len(ampl))

    for i in range(0, half_window):
        # Mirror at left edge.
        idx = min(ndata-1, half_window-i)
        sol[i] = ampl[idx]

        # Mirror at right edge
        idx = max(0, ndata-2-i)
        sol[ndata+half_window+i] = ampl[idx]
    
    for i in range(len(ampl)):
      #print i, i+half_window
      std[i] =  numpy.median(sol[i:i+(2*half_window)])  

    return std


ionmodel = lofar.expion.ionosphere.IonosphericModel(['globaldb']) # 87 seems better

# save for quicker access
#numpy.save('amplitude_array.npy',ionmodel.amplitudes)

amplitude_arraytmp = numpy.load('amplitude_array.npy')


print numpy.shape(amplitude_arraytmp)

amplitude_array = rebin(amplitude_arraytmp, (75, 2912/8, 62, 1, 2))

freqsarray = rebin(numpy.copy(ionmodel.freqs),(2912/8,))

print numpy.shape(amplitude_array)
print numpy.shape(freqsarray)


source_id  = 0
show_plot  = False
freqs_tmp = freqsarray

# indices of SB in the ionmodel.freqs[:] list
goodfreq_el =[0,   1,   2,   3,   4,   5,   6,   7,   8,   9,  10,  11,  12,\
        13,  14,  15,  16,  17,  18,  19,  20,  21,  22,  23,  24,  25,\
        26,  27,  28,  29,  30,  31,  32,  33,  34,  35,  36,  37,  38,\
        39,  40,  41,  42,  43,  44,  45,  46,  47,  48,  49,  50,  51,\
        52,  53,  54,  55,  56,  57,  58,  59,  60,  61,  62,  63,  64,\
        65,  66,  67,  68,  69,  70,  71,  72,  73,  74,  75,  76,  77,\
        78,  79,  80,  81,  82,  83,  84,  85,  86,  87,  88,  89,  90,\
        91,  92,  93,  94,  95,  96,  97,  98,  99, 100, 101, 102, 103,\
       104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116,\
       117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129,\
       130, 131, 132, 133, 134, 136, 137, 138, 139, 140, 141, 142,\
       143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155,\
       156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168,\
       169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181,\
       182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194,\
       195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207,\
       208, 209, 210, 211, 212, 213, 214, 217, 218, 219, 220,\
       221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233,\
       234, 235, 236, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246,\
       247, 248, 249, 250, 251, 252, 253, 254, 255, 256, 257, 258, 259,\
       260, 261, 262, 263, 264, 265, 266, 267, 268, 269, 270, 271, 272,\
       273, 274, 275, 276, 277, 278, 279, 280, 281, 282, 283, 284, 285,\
       286, 287, 288, 289, 290, 291, 292, 293, 294, 295, 296, 297, \
       301, 302, 303, 304, 305, 306, 307, 308, 309, 310, 311,\
       312, 313, 314, 315, 316, 317, 318, 319, 320, 321, 322, 323, 324,\
       325, 326, 327, 328, 329, 330, 331, 332, 333, 334, 335, 336, 337,\
       338, 339, 349, 350,\
       351, 352, 353, 354, 355, 356, 357, 358, 359, 360]

# bad 302, 300
# bad 343,344,345,346,347
# bad 361,362,363

#goodfreq_el = numpy.arange(0,len(freqs_tmp))

freqs = numpy.copy(freqs_tmp[goodfreq_el])
times = numpy.arange(0,len(ionmodel.times[:]),1)
freqs_new  = numpy.arange(numpy.min(freqs)+1e3,numpy.max(freqs)-1e3, 195.3125e3)
amps_array = numpy.zeros( (len(ionmodel.stations[:]),len(ionmodel.times[:]),len(freqs_new),2), dtype='float')
numpy.save('freqs_for_amplitude_array.npy',freqs_new)



print numpy.shape(freqs)
print numpy.shape(times)

for antenna_id in range(0,len(ionmodel.stations[:])):
 amp_xx_tmp = numpy.copy(amplitude_array[:,:,antenna_id,source_id,0]) # array(time, freq)
 amp_yy_tmp = numpy.copy(amplitude_array[:,:,antenna_id,source_id,1]) # array(time, freq)

 amp_xx = numpy.copy(amp_xx_tmp[:,goodfreq_el])
 amp_yy = numpy.copy(amp_xx_tmp[:,goodfreq_el])
 
 print 'Doing', ionmodel.stations[antenna_id]
   
 if show_plot:
   minscale = numpy.median(amp_xx)*0.3
   maxscale = numpy.median(amp_xx)*2.0
   subplots_adjust(wspace = 0.6)

   matplotlib.pyplot.subplot(121)
   matplotlib.pyplot.imshow(amp_xx, vmax=maxscale, vmin=minscale, aspect='auto')
   matplotlib.pyplot.xlabel('calibrator SB (incrasing freq)')
   matplotlib.pyplot.ylabel('time')
   matplotlib.pyplot.title(ionmodel.stations[antenna_id]+ ' XX ampl')
   matplotlib.pyplot.colorbar()#(orientation='horizontal')
   #matplotlib.pyplot.tight_layout()
   matplotlib.pyplot.subplot(122)
   matplotlib.pyplot.imshow(amp_yy, vmax=maxscale, vmin=minscale, aspect='auto')
   matplotlib.pyplot.xlabel('calibrator SB (incrasing freq)')
   matplotlib.pyplot.ylabel('time')
   matplotlib.pyplot.title(ionmodel.stations[antenna_id]+ ' YY ampl')
   matplotlib.pyplot.colorbar()#(orientation='horizontal')
   #matplotlib.pyplot.tight_layout()
   matplotlib.pyplot.show()


 amp_xx = scipy.ndimage.filters.median_filter(amp_xx, (15,3))
 amp_xx = scipy.ndimage.filters.median_filter(amp_xx, (150,1))
 amp_yy = scipy.ndimage.filters.median_filter(amp_yy, (15,3))
 amp_yy = scipy.ndimage.filters.median_filter(amp_yy, (150,1))
 print numpy.shape(amp_xx)
 
 if show_plot:
   subplots_adjust(wspace = 0.6)

   matplotlib.pyplot.subplot(121)
   matplotlib.pyplot.imshow(amp_xx, vmax=maxscale, vmin=minscale, aspect='auto')
   matplotlib.pyplot.xlabel('calibrator SB (incrasing freq)')
   matplotlib.pyplot.ylabel('time')
   matplotlib.pyplot.title(ionmodel.stations[antenna_id]+' Smoothed XX ampl')
   matplotlib.pyplot.colorbar()#(orientation='horizontal')
   matplotlib.pyplot.subplot(122)
   matplotlib.pyplot.imshow(amp_yy, vmax=maxscale, vmin=minscale, aspect='auto')
   matplotlib.pyplot.xlabel('calibrator SB (incrasing freq)')
   matplotlib.pyplot.ylabel('time')
   matplotlib.pyplot.title(ionmodel.stations[antenna_id]+' Smoothed YY ampl')
   matplotlib.pyplot.colorbar()#(orientation='horizontal')
   matplotlib.pyplot.show()

  #2D interpol XX
 for time in range(0,len(ionmodel.times[:])):
  
   amps      = numpy.copy(amp_xx[time,:])
   fl = scipy.interpolate.interp1d(freqs,amps,kind='slinear',bounds_error=False)
   ynew= fl(freqs_new)
   amps_array[antenna_id,time,:,0] = numpy.copy(savitzky_golay(ynew, 17, 2))
   if show_plot:
      #print numpy.shape(freqs_new), numpy.shape(amps_array)
      matplotlib.pyplot.plot(freqs_new,amps_array[antenna_id,time,:,0])
 if show_plot:
   matplotlib.pyplot.show()

  #2D interpol YY
 for time in range(0,len(ionmodel.times[:])):
   amps      = numpy.copy(amp_yy[time,:])
   fl = scipy.interpolate.interp1d(freqs,amps,kind='slinear',bounds_error=False)
   ynew= fl(freqs_new)
   amps_array[antenna_id,time,:,1] = numpy.copy(savitzky_golay(ynew, 17, 2))
   if show_plot:
      matplotlib.pyplot.plot(freqs_new,amps_array[antenna_id,time,:,1])
 if show_plot:
   matplotlib.pyplot.show()





   #sys.exit()
 #for SB in range(len(ionmodel.freqs[:])):
 # print 'Doing SB', SB
 # amp_xx[:,SB] = median_window_filter(amp_xx[:,SB], 200,5.0) 
 # amp_xx[:,SB] = median_window_filter (amp_xx[:,SB], 50, 3.0)
 # amp_xx[:,SB] = median_window_filter(amp_xx[:,SB], 50, 2.5)

 #matplotlib.pyplot.imshow(amp_xx, vmax=numpy.median(amp_xx)*2.0, vmin=numpy.median(amp_xx)*0.3, aspect='auto')
 #matplotlib.pyplot.xlabel('calibrator SB (incrasing freq)')
 #matplotlib.pyplot.ylabel('time')
 #matplotlib.pyplot.colorbar()
 #matplotlib.pyplot.show()
numpy.save('3C196_amplitude_array.npy',amps_array)
print numpy.shape(amps_array)
