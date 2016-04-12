"""
Script to smooth and normalize amplitude solutions
"""
import argparse
from argparse import RawTextHelpFormatter
import pyrap.tables as pt
import numpy
import os
import lofar.parmdb
import math
import scipy.signal
import shutil
import multiprocessing
import itertools
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline
from scipy.interpolate import LSQUnivariateSpline
from scipy.interpolate import LSQBivariateSpline
import sys
from scipy.signal import medfilt
import scipy.ndimage #.filters.median_filter
import astropy.convolution
import matplotlib as mpl
import matplotlib.image as mpimg

plotting = True

if plotting:
  mpl.rc('font',size =8 )
  mpl.rc('figure.subplot',left=0.05, bottom=0.05, right=0.95, top=0.95 )


def std(inputData, Zero=False, axis=None, dtype=None):
	"""
	Robust estimator of the standard deviation of a data set.  
	
	Based on the robust_sigma function from the AstroIDL User's Library.
	
	.. versionchanged:: 1.0.3
		Added the 'axis' and 'dtype' keywords to make this function more
		compatible with numpy.std()
	"""
	epsilon = 1.0e-20
	if axis is not None:
		fnc = lambda x: std(x, dtype=dtype)
		sigma = numpy.apply_along_axis(fnc, axis, inputData)
	else:
		data = inputData.ravel()
		if type(data).__name__ == "MaskedArray":
			data = data.compressed()
		if dtype is not None:
			data = data.astype(dtype)
			
		if Zero:
			data0 = 0.0
		else:
			data0 = numpy.median(data)
		maxAbsDev = numpy.median(numpy.abs(data-data0)) / 0.6745
		if maxAbsDev < epsilon:
			maxAbsDev = (numpy.abs(data-data0)).mean() / 0.8000
		if maxAbsDev < epsilon:
			sigma = 0.0
			return sigma
			
		u = (data-data0) / 6.0 / maxAbsDev
		u2 = u**2.0
		good = numpy.where( u2 <= 1.0 )
		good = good[0]
		if len(good) < 3:
			print "WARNING:  Distribution is too strange to compute standard deviation"
			sigma = -1.0
			return sigma
			
		numerator = ((data[good]-data0)**2.0 * (1.0-u2[good])**2.0).sum()
		nElements = (data.ravel()).shape[0]
		denominator = ((1.0-u2[good])*(1.0-5.0*u2[good])).sum()
		sigma = nElements*numerator / (denominator*(denominator-1.0))
		if sigma > 0:
			sigma = math.sqrt(sigma)
		else:
			sigma = 0.0
			
	return sigma

def findscatter(datavector):

  shifted_vec = numpy.roll(datavector, 1)
  #scatter = sum(abs(shifted_vec - datavector))/numpy.float(len(datavector))
  scatter = numpy.median(abs(shifted_vec - datavector))
  return scatter

def findnoisevec(datavector):
  shifted_vec = numpy.roll(datavector, 1)
  scatter_vec = (abs(shifted_vec - datavector))
  #scatter_vec = medfilt(scatter_vec,21)
  scatter_vec = scipy.ndimage.filters.median_filter(scatter_vec,9, mode='mirror')
  
  # now smooth
  gauss = astropy.convolution.Gaussian1DKernel(stddev=4.0)
  scatter_vec = astropy.convolution.convolve(scatter_vec,gauss , boundary='extend')
  
  # normalize scatter_vec
  scatter_vec = scatter_vec/numpy.mean(scatter_vec)
  
  return scatter_vec


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
        sol[idx] = True

    mask = sol_flag[half_window:half_window + ndata]

    for i in range(len(mask)):
        if mask[i]:
           ampl_tot_copy[i] = median_array[half_window+i] # fixed 2012
    return ampl_tot_copy


def spline1D(amp_orig):
  # to compute knot points
  f = lambda m, n: [i*n//m + n//(2*m) for i in range(m)]
  
  # expand array and mirror full array around edges
  ndata = len(amp_orig)
  amp = numpy.zeros(ndata+2*ndata)
  amp[ndata:ndata+ndata] = amp_orig

  for i in range(0, ndata):
    # Mirror at left edge.
    idx = min(ndata-1, ndata-i)
    amp[i] = amp_orig[idx]
    # Mirror at right edge
    idx = max(0, ndata-2-i)
    amp[ndata+ndata+i] = amp_orig[idx]

  # work in log-sapce
  amp = numpy.log10(amp)
  weights = (0.*numpy.copy(amp)) + 1 # initialize weights to 1           
             
  # filter bad data and determine average scatter of amplitudes
  idx = numpy.where(amp != 0.0) # log10(1.0) = 0.0
  if len(idx) != 0:
    scatter = findscatter(amp[idx])
    # remove some really bad stuff, by putting weights to zero.
    idxbadi1 = numpy.where(amp > (numpy.median(amp) + (8.*std(amp)))) 
    weights[idxbadi1] = 1e-10 # small value, zero generates NaN in spline
    idxbadi2 = numpy.where(amp < (numpy.median(amp) - (8.*std(amp))))
    weights[idxbadi2] = 1e-10  # small value, zero generates NaN in spline
    #print 'scatter', scatter,  numpy.median(amp), idxbadi
    #sys.exit()
  else:
    scatter = 0.02 # just that we have a value to prevent crashes in case all amplitudes are 1.0
    print 'No valid data for found for this anntenna'


  # make the noisevec
  if len(idx) !=0 : # at least 1 good data point             
    if len(amp[idx]) > 30:  # so at least 30/3 = 10 good data points           
      # create noise vector           
      noisevec = findnoisevec(amp)
    else:
      noisevec = (numpy.copy(amp) * 0.) + 1.0 # just make constant noise, if we have too little datapoints
  else:
    noisevec = (numpy.copy(amp) * 0.) + 1.0 # just make constant noise, if we have too little datapoints

  
  if scatter < 0.005:
  #Interior knots t must satisfy Schoenberg-Whitney conditions
     scatter = 0.005 # otherwise we fit more parameters than we have data points
  knotfactor = 0.4e3*scatter  # normalize based on trial and error
             
  
  timevec = numpy.arange(0,len(amp))
  knotvec = f(numpy.int(len(amp)/knotfactor),len(amp))
  #print antenna, 'knots', knotvec, noisevec[knotvec]
             
  # simple optimization knot selection for vectors that have at least 30 data points
  # based on the noisevector
  # removes even numbered knots if the noise is high
  knotvec_copy = numpy.copy(knotvec) # otherwise tcopy is updated as well
  if len(timevec) > 30 and len(knotvec) > 2:
    for counter, knot in enumerate(knotvec_copy):
      #print counter, knot, noisevec[knot]
      if (counter % 2 == 0) and noisevec[knot] > 1.5: # even index and large noise
        knotvec.remove(knot)
        #print 'Removing knot because of local increase in noise'
             
  #print antenna, 'cleaned knots', knotvec, noisevec[knotvec]
  
  # asign midpoint if not enough data points/20
  if len (knotvec) < 3: # because we are working with a 3x larger mirrored array
    knotvec = [numpy.int(len(timevec)*0.25),numpy.int(len(timevec)/2),numpy.int(len(timevec)*0.75)]
    #print 'extending to', knotvec         
             
  splineorder =  5 #  default
  if len(knotvec) == 3 and scatter > 0.1:
    splineorder = 3 # reduce order, data is  bad
    if scatter > 0.2:
      splineorder = 1 # very bad data
             
  #print 'knots', knotvec
  spl2 = LSQUnivariateSpline(timevec, amp, knotvec, w=weights, k=splineorder)

  # now find bad data devatiating from the fit 30 x scatter 
  residual = numpy.abs(spl2(timevec)-amp)
  idx      = numpy.where(residual > 15.*scatter)

  # second iteration
  if len(idx) != 0:
    ampcopy = numpy.copy(amp)
    ampcopy[idx] = spl2(timevec[idx]) # replace bad amplitudes by model
    spl2 = LSQUnivariateSpline(timevec, ampcopy, knotvec,  w=weights, k=splineorder)

  residual = numpy.abs(spl2(timevec)-amp)
  idx      = numpy.where(residual > 8.*scatter)
  
  # third iteration
  if len(idx) != 0:
    ampcopy = numpy.copy(amp)
    ampcopy[idx] = spl2(timevec[idx]) # replace bad amplitudes by model
    spl2 = LSQUnivariateSpline(timevec, ampcopy, knotvec,  w=weights, k=splineorder)

  # again look at residual, go back to original amp again, find deviating data > 3x scatter
  residual = numpy.abs(spl2(timevec)-amp)
  idx      = numpy.where(residual > 3.*scatter)
  # replace the bad data with model
  model    =spl2(timevec)
  if len(idx) != 0:
    amp[idx] = model[idx]

  # go out of log-space
  amp = 10**amp
 
  amp_clean = amp[ndata:ndata + ndata]
  
  
  # remake for plot 
  #xforplot = numpy.arange(0,len(amp_orig))
  # find the replaced points for plotting
  
  idxbad = numpy.where(amp_clean != amp_orig)      
  n_knots = numpy.int(numpy.ceil(numpy.float(len(knotvec))/3.)) # approxmiate, just for plot

  # return cleaned amplitudes, model, scatter, number of knots, indices of replaced outliers
  return amp_clean, 10**(model[ndata:ndata + ndata]), noisevec[ndata:ndata + ndata], scatter, n_knots, idxbad, weights[ndata:ndata + ndata]


instrument_name_smoothed = 'test_parmdb_1Dspline'
instrument_name = 'L343226_SBgr025-10_uv.dppp.pre-cal_chunk9_12656E813t_0g.merge_amp_parmdbs2'
#instrument_name = 'testparmdb'
#instrument_name = '../facet_patch_185/L343226_SBgr025-10_uv.dppp.pre-cal_chunk9_12656E813t_0g.merge_amp_parmdbs2'
instrument_name = '../facet_patch_639/L343226_SBgr025-10_uv.dppp.pre-cal_chunk9_12656E813t_0g.merge_amp_parmdbs2'

instrument_name          = str(sys.argv[1]) # input
instrument_name_smoothed = str(sys.argv[2]) # ouput




pol_list = ['0:0','1:1']  # need to be updated to work with all four polarizations
gain = 'Gain'

pdb = lofar.parmdb.parmdb(instrument_name)
parms = pdb.getValuesGrid('*')

key_names = parms.keys()
nchans = len(parms[key_names[0]]['freqs'])

#print parms[key_names[0]]['freqs']

# Get station names
antenna_list = set([s.split(':')[-1] for s in pdb.getNames()])


# for plotting
Nr = int(numpy.ceil(numpy.sqrt(len(antenna_list))))
Nc = int(numpy.ceil(numpy.float(len(antenna_list))/Nr))



if nchans < 7: # do 1D spline in this case
 if plotting:
    fa, axa = plt.subplots(Nr, Nc, sharex=True, sharey=True, figsize=(16,12))
    axsa = axa.reshape((Nr*Nc,1)) 
  
 for pol in pol_list:
        for istat,antenna in enumerate(sorted(antenna_list)[::-1]):
            channel_parms_real = [parms[gain + ':' + pol + ':Real:'+ antenna]['values'][:, chan] for chan in range(nchans)]
            channel_parms_imag = [parms[gain + ':' + pol + ':Imag:'+ antenna]['values'][:, chan] for chan in range(nchans)]


            if len(channel_parms_real[0]) > 500:
              fmt = ','
            else:
               fmt = '.'
            ls='none'

            #loop of over channel
            for chan in range(nchans):
                #print 'doing channel', chan
		amp_orig = numpy.sqrt(channel_parms_real[chan]**2 + channel_parms_imag[chan]**2)
                phase = numpy.arctan2(channel_parms_imag[chan], channel_parms_real[chan]**2)

		amp_cleaned, model, noisevec, scatter, n_knots, idxbad, weights = spline1D(amp_orig)
		
		# put back the results
		parms[gain + ':' + pol + ':Real:' + antenna]['values'][:, chan] = numpy.copy(amp_cleaned*numpy.cos(phase))
                parms[gain + ':' + pol + ':Imag:' + antenna]['values'][:, chan] = numpy.copy(amp_cleaned*numpy.sin(phase))
		
		if pol in pol_list[0]:
		  cc = 'blue'
		  ccf = 'orange'
		else:
		  cc = 'green'
		  ccf= 'red'
		
		timevec = numpy.arange(0,len(amp_orig))
		
		# only plot channel 0, just to verify code
		if plotting and chan == 0:
		  axsa[istat][0].plot(timevec,amp_cleaned,  marker=fmt, ls=ls, markersize=200/len(amp_cleaned), c=cc)
		  axsa[istat][0].plot(timevec,noisevec, c=cc, lw=0.75, ls='--')
		  
		  if pol in pol_list[0]:
		    axsa[istat][0].annotate('scatter=' +'{:.2g}'.format(scatter), xy=(0.5,0.15), color=cc,textcoords='axes fraction')
		    axsa[istat][0].annotate('#knots=' +'{:d}'.format(n_knots), xy=(0.01,0.15), color=cc,textcoords='axes fraction') # we divded by three beucase we mirrored the array
		  else:
		    axsa[istat][0].annotate('scatter=' +'{:.2g}'.format(scatter), xy=(0.5,0.02), color=cc, textcoords='axes fraction')
		    axsa[istat][0].annotate('#knots=' +'{:d}'.format(n_knots), xy=(0.01,0.02), color=cc,textcoords='axes fraction')
		  
		  if len(idxbad) > 0:
		    axsa[istat][0].plot(timevec[idxbad],amp_orig[idxbad],  marker='o', c=ccf,ls=ls, markersize=4)
                  
                  idxbadi = numpy.where(weights < 1.0)
                  if len(idxbadi) > 0:
                    axsa[istat][0].plot(timevec[idxbadi],amp_orig[idxbadi],  marker='o', c='black',ls=ls, markersize=4)

		  axsa[istat][0].plot(timevec, model, c=ccf, lw=1.0)
		  axsa[istat][0].set_title(antenna)
		  axsa[istat][0].set_ylim(-0.3,2)
		  axsa[istat][0].set_xlim(0,max(timevec))
		  
		  # to make sure we really removed the bad stuff
		  #axsa[istat][0].plot(xforplot,amp_clean,  marker='o', c='purple',ls=ls, markersize=4)
		  
		  #axsa[istat][0].set_xlabel('time [sample]')
		  #axsa[istat][0].set_ylabel('amplitude')

 if plotting:
   plt.show()
   fa.savefig('test.png', dpi=100)

 if os.path.exists(instrument_name_smoothed):
    shutil.rmtree(instrument_name_smoothed)
 pdbnew = lofar.parmdb.parmdb(instrument_name_smoothed, create=True)
 pdbnew.addValues(parms)
 pdbnew.flush()

else:
 print '2D spline'




 #xsize  = 6
 #ysize  = 63
 
 #tx = numpy.linspace(-5,5,xsize)
 #ty = numpy.linspace(-5,5,ysize)
 #domain = numpy.meshgrid(tx,ty)
 ##print numpy.shape(domain)

 #X, Y = domain
 

 #Xp, Yp = domain
 #Z = myfunc (*domain)
 #Zp = Z
 #print numpy.shape(Z)

 #X = X.ravel()
 #Y = Y.ravel()
 #Z = Z.ravel()
 
 
 #kx = numpy.linspace(-5,5,3)[1:-1]
 #print kx
 ##sys.exit()
 #ky = kx.copy()
 #spl2 = LSQBivariateSpline(X, Y, Z, kx, ky, kx=3, ky=3)
 #print spl2
 #print numpy.shape(Z)
 ##print Z

 #img1D = spl2.ev(X,Y)
 #img2D = numpy.reshape(img1D,(ysize,xsize))
 

 #plt.subplot(121)
 #plt.imshow(img2D,interpolation='none')
 #plt.subplot(122)
 #plt.imshow(Zp,interpolation='none')
 #plt.show()
 ##sys.exit()

 for pol in pol_list:
        for antenna in sorted(antenna_list)[::-1]:
            channel_parms_real = [parms[gain + ':' + pol + ':Real:'+ antenna]['values'][:, chan] for chan in range(nchans)]
            channel_parms_imag = [parms[gain + ':' + pol + ':Imag:'+ antenna]['values'][:, chan] for chan in range(nchans)]


            channel_parms_real =  numpy.asarray(channel_parms_real)
            channel_parms_imag =  numpy.asarray(channel_parms_imag)
            amp_orig = numpy.sqrt(channel_parms_real[:]**2 + channel_parms_imag[:]**2)
            print numpy.shape(amp_orig)
            orinal_size = numpy.shape(amp_orig)
            
            
            # padd image by reflection around axis
            amp_orig = numpy.pad(amp_orig, ((numpy.shape(amp_orig)[0],numpy.shape(amp_orig)[0]),\
	               (numpy.shape(amp_orig)[1],numpy.shape(amp_orig)[1])),mode='reflect')
            
            print numpy.shape(amp_orig)
            
            
            # ------
            amp = numpy.log10(amp_orig)

            xsize  = numpy.shape(amp)[0]
            ysize  = numpy.shape(amp)[1]
            print '(xsize, ysize)', xsize, ysize
            amp_transpose = numpy.transpose(amp)
            
            tx = numpy.linspace(0,numpy.shape(amp)[0]-1,numpy.shape(amp)[0])
            ty = numpy.linspace(0,numpy.shape(amp)[1]-1,numpy.shape(amp)[1])

            domain = numpy.meshgrid(ty,tx)
            X, Y = domain
            
    
            tx = numpy.linspace(0,xsize,5)[1:-1]  # this is the freq axis
            ty = numpy.linspace(0,ysize,19)[1:-1]  # this is the time axis
            korderx = 3 # freq direction
            kordery = 5 # time direction
            

            print tx, ty
   
            X = X.ravel()
            Y = Y.ravel()
            Z = amp_transpose.ravel()
            #print numpy.shape(Z), numpy.shape(amp_transpose.ravel())
            
            # have to reverse ty, tx order here, not clear to me why, I confused the order X, Y and tranpose?
            spl2 = LSQBivariateSpline(X, Y, Z, ty, tx, kx=korderx, ky=kordery)
            #domain = numpy.meshgrid(xaxis,yaxis)
            #print numpy.shape(domain)
        
            img1D = spl2.ev(X,Y)
            img2D = numpy.reshape(img1D,(ysize,xsize))
            
            padidx = [orinal_size[0],2*orinal_size[0],orinal_size[1],2*orinal_size[1]]
            print padidx
            
            plt.subplot(1,3,1)
            amp_orig = numpy.transpose(amp_orig)
            
            print  orinal_size
            print numpy.shape(amp_orig)
            print numpy.shape(amp_orig[orinal_size[1]:2*orinal_size[1],orinal_size[0]:2*orinal_size[0]])
            
            plt.imshow(amp_orig[orinal_size[1]:2*orinal_size[1],orinal_size[0]:2*orinal_size[0]],interpolation='none',origin='lower',clim=(0.5, 1.5),aspect='auto')
            plt.xlabel('freq')
            plt.ylabel('time')
            plt.title('Original')
         
            
            plt.subplot(1,3,2)
            plt.imshow(10**(img2D[orinal_size[1]:2*orinal_size[1],orinal_size[0]:2*orinal_size[0]]),interpolation='none',origin='lower',aspect='auto', clim=(0.5,1.5))            
            plt.xlabel('freq')
            plt.ylabel('time')
            plt.title('2D spline fit')
            
            
            plt.subplot(1,3,3)
            plt.imshow(abs((10**(img2D[orinal_size[1]:2*orinal_size[1],orinal_size[0]:2*orinal_size[0]]))-\
	              (amp_orig[orinal_size[1]:2*orinal_size[1],orinal_size[0]:2*orinal_size[0]])),interpolation='none',origin='lower',clim=(0.0, 0.3),aspect='auto')
            
            plt.xlabel('freq')
            plt.ylabel('time')
            plt.title('2D spline fit')
            
            plt.show()
            
            sys.exit()
            amp = numpy.copy(numpy.log10(amp))
            
          
            
