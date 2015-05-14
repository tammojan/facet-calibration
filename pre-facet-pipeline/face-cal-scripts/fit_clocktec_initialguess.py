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
import pp
import scipy.signal
from pylab import *
#import lofar.expion.fitting as fitting

pi = numpy.pi
showplot = False

c  = 2.99792458e8

ionmodel = lofar.expion.ionosphere.IonosphericModel(['../startid_0.626704232754/L239640/L239640.globaldb'])

source_id     = 0  # source ID in global_db (usually 0)
ncpus         = 2 # number of CPU for parallel fitting

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

def fit_dTEC_dclock_dFR(phases_rr, phases_ll, amp_rr, amp_ll, freq, distance_station):
      c = 2.99792458e8
      freq_old = numpy.copy(freq)
      # NOTE THE *2 to use rr+ll instead of 0.5*(rr+ll)
      par3complex     = lambda p, freq, y: abs(numpy.cos((4.*pi*p[0]*freq) - (2.*8.44797245e9*p[1]/freq) + p[2]) - numpy.cos(y)) + abs(numpy.sin((4.*pi*p[0]*freq) - (2.*8.44797245e9*p[1]/freq) + p[2]) - numpy.sin(y))
      par2complex     = lambda p, freq, y: abs(numpy.cos((4.*pi*p[0]*freq) - (2.*8.44797245e9*p[1]/freq)) - numpy.cos(y)) + abs(numpy.sin((4.*pi*p[0]*freq) - (2.*8.44797245e9*p[1]/freq)) - numpy.sin(y))
      rmwavcomplex    = lambda RM, wav, y: abs(numpy.cos(2.*RM[0]*wav*wav)  - numpy.cos(y)) + abs(numpy.sin(2.*RM[0]*wav*wav)  - numpy.sin(y))

      # NOTE THE *2 to use rr+ll instead of 0.5*(rr+ll) 
      par3complex_w   = lambda p, freq, y: abs(numpy.cos((4.*pi*p[0]*freq) - (2.*8.44797245e9*p[1]/freq) + p[2]) - numpy.cos(y))*(freq/1e5) + abs(numpy.sin((4.*pi*p[0]*freq) - (2.*8.44797245e9*p[1]/freq) + p[2]) - numpy.sin(y))*(freq/1e5)
      par2complex_w   = lambda p, freq, y: abs(numpy.cos((4.*pi*p[0]*freq) - (2.*8.44797245e9*p[1]/freq)) - numpy.cos(y))*(freq/1e5) + abs(numpy.sin((4.*pi*p[0]*freq) - (2.*8.44797245e9*p[1]/freq)) - numpy.sin(y))*(freq/1e5)
      rmwavcomplex_w  = lambda RM, wav, y: abs(numpy.cos(2.*RM[0]*wav*wav)  - numpy.cos(y))/wav + abs(numpy.sin(2.*RM[0]*wav*wav)  - numpy.sin(y))/wav

      plotrm          = lambda RM, wav: numpy.mod( (2.*RM*wav*wav) +1.0*pi, 2.*pi)  -1.0*pi # notice the factor of 2
      fitfuncfastplot    = lambda p, freq: numpy.mod((4.*pi*p[0]*freq) - (2.*8.44797245e9*p[1]/freq) + (p[2])+ 1.0*pi, 2.*pi) -1.0*pi 


      idx_rr    = numpy.where(amp_rr != 1.0)
      if len(idx_rr) != 0:
       freq       = freq[:][idx_rr]
       phases_rr  = phases_rr[:][idx_rr]
       phases_ll  = phases_ll[:][idx_rr]
       amp_ll     = amp_ll[:][idx_rr] 

      idx_ll     = numpy.where(amp_ll != 1.0)
      if len(idx_ll) != 0:
       freq       = freq[:][idx_ll]
       phases_rr  = phases_rr[:][idx_ll]
       phases_ll  = phases_ll[:][idx_ll]
   
      if len(freq) < len(freq_old):
        print 'Number of filtered out data points:',  len(freq_old)-len(freq)
      # ---------- end filter bad data -----------

      if len(freq) != 0: # prepare and make arrays if there is valid data 
        freq      = (freq[0:len(freq_old)])
        phases_ll = (phases_ll[0: len(freq_old)])
        phases_rr = (phases_rr[0: len(freq_old)])
      
        phase       = (phases_rr + phases_ll)      # not divide by 2, then later fix this
        #phase       = (phases_ll + phases_ll)  # temp, just fit RR for now 
        phase_diff  = (phases_rr - phases_ll)      # not divide by 2, then later fix this
      
        wav      = c/freq
        pi = numpy.pi
        chi_old=1e9 
      else:
        phase = (phases_rr[:] + phases_ll[:])
	phase_diff = (phases_rr[:] - phases_ll[:])
      
      if len(freq) > 10: 
      
	# FIND INTIAL GUESS
	for dTEC in numpy.arange(-1.0,1.0, 0.01):
	 for dclock in numpy.arange(-200e-9,200e-9,5e-9):
            phase_model = numpy.mod ( (4.*pi*dclock*freq) - (2.*8.44797245e9*dTEC/freq), 2*pi)  # NOTE THE *2 to use rr+ll instead of 0.5*(rr+ll)
            phase_data  = numpy.mod (phase, 2*pi)
	    angle       = pi - numpy.abs(numpy.abs(phase_model - phase_data) - pi)
	    chi = numpy.sum(angle)

	    if chi < chi_old:
	      chi_old = chi
	      fitguess     = [dclock,dTEC]
	      #print 'Better fit', dclock, dTEC

	fitguess_1 = numpy.copy(fitguess)
	#print 'iter 1', fitguess

	for dTEC in numpy.arange(fitguess_1[1]-0.02,fitguess_1[1]+0.02, 0.002):
	 for dclock in numpy.arange(fitguess_1[0]-8e-9,fitguess_1[0]+ 8e-9,1e-9):
            phase_model = numpy.mod ( (4.*pi*dclock*freq) - (2.*8.44797245e9*dTEC/freq), 2*pi)  # NOTE THE *2 to use rr+ll instead of 0.5*(rr+ll)
            phase_data  = numpy.mod (phase, 2*pi)
	    angle       = pi - numpy.abs(numpy.abs(phase_model - phase_data) - pi)
	    chi = numpy.sum(angle)

	    if chi < chi_old:
	      chi_old = chi
	      fitguess     = [dclock,dTEC]
	      #print 'Better fit', dclock, dTEC


	#print 'iter 2', fitguess

	chi_old = 1e9
	for dFR in numpy.arange(-0.1,0.1,2e-4):
          phase_model = numpy.mod (2.*dFR*wav*wav, 2*pi) # notice the factor of 2
	  phase_data  = numpy.mod (phase_diff, 2*pi)
          angle       = pi - numpy.abs(numpy.abs(phase_model - phase_data) - pi)
	  chi = numpy.sum(angle)

	  if chi < chi_old:
	    chi_old = chi
	    fitrmguess     = dFR
	    #print 'Better fit', fitrmguess

	fitrmguess_1 = numpy.copy(fitrmguess)
	for dFR in numpy.arange(fitrmguess_1-5e-4,fitrmguess_1+5e-4,0.5e-5):
          phase_model = numpy.mod (2.*dFR*wav*wav, 2*pi) # notice the factor of 2
	  phase_data  = numpy.mod (phase_diff, 2*pi)
          angle       = pi - numpy.abs(numpy.abs(phase_model - phase_data) - pi)
	  chi = numpy.sum(angle)

	  if chi < chi_old:
	    chi_old = chi
	    fitrmguess     = dFR
	    #print 'Better fit', fitrmguess

        # DO THE FITTING 
        # SOLVE Clock-TEC anticorrelation problem on short baselines             
        freq = freq.astype(numpy.float64)
	phase= phase.astype(numpy.float64) #epsfcn=1e-7
      
        if distance_station < 0. :   #15.0*1e3: DOES NOT WORK, NEED 3 par FIT
          fitresult, success  = scipy.optimize.leastsq(par2complex, fitguess,    args=(freq, phase))
	  #fitresult = fitguess
        else:
          fitresult, success   = scipy.optimize.leastsq(par3complex, [fitguess[0], fitguess[1], 0.0], args=(freq, phase),maxfev=10000)
	  #fitresult = [fitguess[0], fitguess[1], 0.0]
          #print fitresult, success
        fitresultrm_wav, success = scipy.optimize.leastsq(rmwavcomplex, [fitrmguess], args=(wav, phase_diff))

      else:
         print 'No valid data found'
         fitresult = [0.0,0.0, 0.0]
	 fitresultrm_wav= 0.0
            
      
      
      showplot = False

      if showplot:  
        
	#if len(fitresult ==2): 
	#  fitresult =[fitresult[0], fitresult[1],0]
	#  print 'Here'
	#  #fitresult = fitresult #   [fitguess[0],  fitguess[1], 0]
	
	#fitresult = [[fitresultguess[0],fitresultguess[0],0.0]]
	
        matplotlib.pyplot.plot(freq, numpy.mod(phase + 1.0*pi, 2.*pi) -1.0*pi, 'or' )
        matplotlib.pyplot.plot(freq, numpy.mod(phase_diff + 1.0*pi, 2.*pi) -1.0*pi , '.', color='purple' )			   

        TEC   = numpy.mod((-8.44797245e9*(2.*fitresult[1])/freq)+numpy.pi, 2*pi) - pi   # notice factor of 2 because rr+ll
        Clock = numpy.mod((2.*numpy.pi*   2.*fitresult[0]*freq )+numpy.pi, 2*pi) - pi   # notice factor of 2 because rr+ll
    
        phase_total = (2.*numpy.pi*2.*fitresult[0]*freq)+(-8.44797245e9*(2.*fitresult[1])/freq)+fitresult[2]
        residual    = numpy.mod(phase-phase_total+pi,2.*pi)-pi
        matplotlib.pyplot.plot(freq, residual, '.', color='yellow')
    
        idxl = int(min(freq_old)/1e4) 
	idxh = int(max(freq_old)/1e4) 
        bigfreqaxistmp = range(idxl, idxh)
	bigfreqaxis    =  numpy.array([float(i) for i in bigfreqaxistmp])
	bigfreqaxis    = bigfreqaxis*1e4
	
	matplotlib.pyplot.plot (bigfreqaxis, fitfuncfastplot(fitresult, bigfreqaxis[:]), "r-")	
	matplotlib.pyplot.plot (bigfreqaxis, plotrm(fitresultrm_wav, c/bigfreqaxis[:]), "-", color='purple')
	#matplotlib.pyplot.plot (freq, fitfuncfastplot(fitresult, bigfreqaxis[:]), "r-")	
	
	matplotlib.pyplot.plot(freq, Clock,  ',g') 
	matplotlib.pyplot.plot(freq, TEC,  ',b') 
	matplotlib.pyplot.xlabel('freq')
	matplotlib.pyplot.ylabel('phase')

	matplotlib.pyplot.show()
      
      return [fitresult[0], fitresult[1],fitresult[2], fitresultrm_wav]
   





time_id       = 100
pol_id        = 0
#antenna_id    = 29
refantenna_id = 0
source_id     = 0
goodstartidx  = 0

CScorrect     = False

if CScorrect:
  csclockvals = numpy.load('../CS_clocks.npy')

A = numpy.zeros((len(ionmodel.freqs[:]), 2), dtype = float)
A[:,0] = ionmodel.freqs[:]*2*pi
A[:,1] = -8.44797245e9/ionmodel.freqs[:]
sol = numpy.zeros((len(ionmodel.stations), 2))

clockfit = 0.*numpy.copy(ionmodel.times)
TECfit   = 0.*numpy.copy(ionmodel.times)
RMfit   = 0.*numpy.copy(ionmodel.times)
phaseoffset = 0.*numpy.copy(ionmodel.times)

clockarray = numpy.zeros([len(ionmodel.times),len(ionmodel.stations)])
tecarray   = numpy.zeros([len(ionmodel.times),len(ionmodel.stations)])
rmarray    = numpy.zeros([len(ionmodel.times),len(ionmodel.stations)])
phaseoffsetarray = numpy.zeros([len(ionmodel.times),len(ionmodel.stations)])

print 'REF STATION:', ionmodel.stations[refantenna_id]
print '# TIMESLOTS  ', len(ionmodel.times)
print '# FREQUENCIES', len(ionmodel.freqs)

ppservers = ()
# Creates jobserver with ncpus workers
job_server = pp.Server(ncpus, ppservers=ppservers)
print "Starting pp with", job_server.get_ncpus(), "workers"


phases_all = numpy.copy(ionmodel.phases)
#phases_all = numpy.load('../phases_3C196.npy')


start_time_id = 0
stop_time_id  = len(ionmodel.phases[:,0,0,0])

if CScorrect:
 # CORRECT CS CLOCKS
 for time_id in range(start_time_id,stop_time_id):
   print time_id
   phases_rr = numpy.copy(phases_all[time_id,:,:,source_id,0])     
   phases_ll = numpy.copy(phases_all[time_id,:,:,source_id,1])

   #RR  correct
   for ss in range(0,len(ionmodel.stations[:])):
         sol[ss,0] = (csclockvals[ss,0]) # clock RR
         sol[ss,1] = 0.0  # TEC
	 #print 'RR', sol[ss,0]
   phases_rr = phases_rr - dot(A, sol.T)

   #LL  correct
   for ss in range(0,len(ionmodel.stations[:])):
       sol[ss,0] = (csclockvals[ss,1]) # clock LL
       sol[ss,1] = 0.0  # TEC 
       #print 'LL', sol[ss,0]  
   phases_ll = phases_ll - dot(A, sol.T) 

   phases_all[time_id,:,:,source_id,0] = numpy.copy(phases_rr)
   phases_all[time_id,:,:,source_id,1] = numpy.copy(phases_ll)


# count number of RS stations
N_RS = 0
for ss in range(0,len(ionmodel.stations[:])):
  sname = ionmodel.stations[ss]
  if sname[0:2] == 'RS':
    N_RS=N_RS+1

#for antenna_id in range(3,  4):
#for antenna_id in range(len(ionmodel.stations[:])-N_RS,  len(ionmodel.stations[:])) :
for antenna_id in range(1,  len(ionmodel.stations[:])) :


 if antenna_id != refantenna_id:
   stationspos      = ionmodel.station_positions[refantenna_id] - ionmodel.station_positions[antenna_id]
   distance_station = numpy.sqrt(stationspos[0]**2 + stationspos[1]**2 + stationspos[2]**2)
   print 'Distance stations to reference station', distance_station/1e3, ' km'

   jobs = []

   for time_id in range(start_time_id,stop_time_id):
   #for time_id in range(3050,3051):
     freq       = numpy.copy(ionmodel.freqs)
          
     
     phases_rr = numpy.copy((phases_all[time_id,:,antenna_id,source_id,0] \
	  -phases_all[time_id,:,refantenna_id,source_id,0])) 
     phases_ll = numpy.copy((phases_all[time_id,:,antenna_id,source_id,1] \
	  -phases_all[time_id,:,refantenna_id,source_id,1])) 

     # -------- filter bad data: where amplitudes equal 1.0 -----------
     amp_rr = numpy.copy(ionmodel.amplitudes[time_id,:,antenna_id,source_id,0]) 
     amp_ll = numpy.copy(ionmodel.amplitudes[time_id,:,antenna_id,source_id,1]) 

     jobs.append(job_server.submit(fit_dTEC_dclock_dFR,(phases_rr, phases_ll, amp_rr, amp_ll, freq, distance_station),(), ("numpy","scipy.optimize","matplotlib.pyplot",)))
     print 'Submitting job #', time_id


   i= 0
   for job in jobs:
     fitresult = job()

     clockfit[start_time_id+i]  = fitresult[0]
     TECfit[start_time_id+i]    = fitresult[1]
     RMfit[start_time_id+i]     = fitresult[3]
     phaseoffset[start_time_id+i]  =  fitresult[2]

     clockarray[start_time_id+i,antenna_id] = fitresult[0]
     tecarray[start_time_id+i,antenna_id]   = fitresult[1]
     rmarray[start_time_id+i,antenna_id]    = fitresult[3]
     phaseoffsetarray[start_time_id+i,antenna_id] = fitresult[2] 
     print 'TIME_ID', start_time_id+i, 'FIT (dclock, dTEC, offset, dRM)', fitresult
     i = i + 1

#os.system('rm -rf fitted_data_dclock_3C196.npy fitted_data_offset_3C196.npy \
#           fitted_data_dTEC_3C196.npy fitted_data_dFR_3C196.npy')

#numpy.save('fitted_data_offset_3C196_ll.npy', phaseoffsetarray)
numpy.save('fitted_data_dclock_3C196_1st.npy', clockarray)
numpy.save('fitted_data_dTEC_3C196_1st.npy', tecarray)
#numpy.save('fitted_data_dFR_3C196_ll.npy', rmarray)



clockarray= numpy.load('fitted_data_dclock_3C196_1st.npy')
tecarray  = numpy.load('fitted_data_dTEC_3C196_1st.npy')
if True:
 for antenna_id in range(0, len(ionmodel.stations[:])):
   print 'Cleaning up Clock and TEC values for: ', ionmodel.stations[antenna_id]
   clockfit = clockarray[:,antenna_id]
   TECfit   = tecarray[:,antenna_id]
 
   clockfit = median_window_filter(clockfit, 150, 5)
   TECfit   = median_window_filter(TECfit, 50, 5)

   
   clockfit = median_window_filter(clockfit, 50, 3)
   TECfit   = median_window_filter(TECfit, 10, 3)

   clockfit = running_median(clockfit, 40)
   TECfit = running_median(TECfit, 10)

   
   clockarray[:,antenna_id] = clockfit
   tecarray[:,antenna_id]   = TECfit
   matplotlib.pyplot.plot(clockarray[:,antenna_id])
   
      
 numpy.save('fitted_data_dclock_3C196_1st.sm.npy',clockarray)
 numpy.save('fitted_data_dTEC_3C196_1st.sm.npy',tecarray)

matplotlib.pyplot.show()

