import losoto._version
import losoto._logging
from losoto.h5parm import h5parm, solWriter, solFetcher
import matplotlib
#matplotlib.use('GTK')
#import lofar.parmdb
import sys
import os
import scipy.ndimage
import time
import numpy
import math

import scipy.signal
from pylab import *
import numpy as np
import shutil
import progressbar
import logging
import pp



def fit_dclock_doffset(phases, freq):
      c = 2.99792458e8
      freq_old = numpy.copy(freq)
      # NOTE THE *2 to use rr+ll instead of 0.5*(rr+ll)
      par3complex     = lambda p, freq, y: abs(numpy.cos((2.*pi*p[0]*freq) - (1.*8.44797245e9*p[1]/freq) + p[2]) - numpy.cos(y)) + abs(numpy.sin((2.*pi*p[0]*freq) - (1.*8.44797245e9*p[1]/freq) + p[2]) - numpy.sin(y))
      par2complex     = lambda p, freq, y: abs(numpy.cos((4.*pi*p[0]*freq) - (2.*8.44797245e9*p[1]/freq)) - numpy.cos(y)) + abs(numpy.sin((4.*pi*p[0]*freq) - (2.*8.44797245e9*p[1]/freq)) - numpy.sin(y))
      rmwavcomplex    = lambda RM, wav, y: abs(numpy.cos(2.*RM[0]*wav*wav)  - numpy.cos(y)) + abs(numpy.sin(2.*RM[0]*wav*wav)  - numpy.sin(y))
      par2slope       = lambda p, freq, y: abs(numpy.cos((2.*pi*p[0]*freq) + p[1]) - numpy.cos(y)) + abs(numpy.sin((2.*pi*p[0]*freq) + p[1]) - numpy.sin(y))
      
      plotrm          = lambda RM, wav: numpy.mod( (1.*RM*wav*wav) +1.0*pi, 2.*pi)  -1.0*pi # notice the factor of 2
      fitfuncfastplot    = lambda p, freq: numpy.mod((2.*pi*p[0]*freq) - (1.*8.44797245e9*p[1]/freq) + (p[2])+ 1.0*pi, 2.*pi) -1.0*pi 


      
      phase       = numpy.copy(phases)      
       
      pi = numpy.pi
      chi_old=1e9 
      
      if len(freq) > 10: 
      
        # FIND INTIAL GUESS
        for doffset in numpy.arange(-pi-0.1,pi+0.1,0.1):
         for dclock in numpy.arange(-80e-9,80e-9,2e-9):
            phase_model = numpy.mod ( (2.*pi*dclock*freq) + doffset, 2*pi)  # NOTE THE *2 to use rr+ll instead of 0.5*(rr+ll)
            phase_data  = numpy.mod (phase, 2*pi)
            angle       = pi - numpy.abs(numpy.abs(phase_model - phase_data) - pi)
            chi = numpy.sum(angle)

            if chi < chi_old:
              chi_old = chi
              fitguess     = [dclock,doffset]
              #print 'Better fit 1', dclock, doffset

        fitguess_1 = numpy.copy(fitguess)
        #print 'iter 1', fitguess

        for doffset in numpy.arange(fitguess_1[1]-0.2,fitguess_1[1]+ 0.2,0.01):
         for dclock in numpy.arange(fitguess_1[0]-4e-9,fitguess_1[0]+ 4e-9,0.1e-10):
            phase_model = numpy.mod ( (2.*pi*dclock*freq) + doffset, 2*pi)  # NOTE THE *2 to use rr+ll instead of 0.5*(rr+ll)
            phase_data  = numpy.mod (phase, 2*pi)
            angle       = pi - numpy.abs(numpy.abs(phase_model - phase_data) - pi)
            chi = numpy.sum(angle)

            if chi < chi_old:
              chi_old = chi
              fitguess     = [dclock,doffset]
              #print 'Better fit 2', dclock, doffset
          
        freq = freq.astype(numpy.float64)
        phase= phase.astype(numpy.float64) #epsfcn=1e-7
      
      	#print 'hello', fitguess
        

        fitresult, success   = scipy.optimize.leastsq(par2slope, [fitguess[0], fitguess[1]], args=(freq, phase),maxfev=10000)
        #print 'ok',fitresult


      else:
         print 'No valid data found'
         fitresult = [0.0,0.0]

            
      
      
   
      show_plot = False
      if show_plot:  
        
        #if len(fitresult ==2): 
        #  fitresult =[fitresult[0], fitresult[1],0]
        #  print 'Here'
        #  #fitresult = fitresult #   [fitguess[0],  fitguess[1], 0]

        #fitresult = [[fitresultguess[0],fitresultguess[0],0.0]]

        matplotlib.pyplot.plot(freq, numpy.mod(phase + 1.0*pi, 2.*pi) -1.0*pi, 'or' )
        matplotlib.pyplot.plot(freq, numpy.mod(phase_diff + 1.0*pi, 2.*pi) -1.0*pi , '.', color='purple' )                         

        TEC   = numpy.mod((-8.44797245e9*(fitresult[1])/freq)+numpy.pi, 2*pi) - pi   # notice factor of 2 because rr+ll
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
      
      return [fitresult[0], fitresult[1]]
   
def fit_dTEC_doffset(phases, freq):
      c = 2.99792458e8
      freq_old = numpy.copy(freq)
      # NOTE THE *2 to use rr+ll instead of 0.5*(rr+ll)
      par3complex     = lambda p, freq, y: abs(numpy.cos((2.*pi*p[0]*freq) - (1.*8.44797245e9*p[1]/freq) + p[2]) - numpy.cos(y)) + abs(numpy.sin((2.*pi*p[0]*freq) - (1.*8.44797245e9*p[1]/freq) + p[2]) - numpy.sin(y))
      par2complex     = lambda p, freq, y: abs(numpy.cos((4.*pi*p[0]*freq) - (2.*8.44797245e9*p[1]/freq)) - numpy.cos(y)) + abs(numpy.sin((4.*pi*p[0]*freq) - (2.*8.44797245e9*p[1]/freq)) - numpy.sin(y))
      rmwavcomplex    = lambda RM, wav, y: abs(numpy.cos(2.*RM[0]*wav*wav)  - numpy.cos(y)) + abs(numpy.sin(2.*RM[0]*wav*wav)  - numpy.sin(y))
      par2slope       = lambda p, freq, y: abs(numpy.cos((-8.44797245e9*p[0]/freq) + p[1]) - numpy.cos(y)) + abs(numpy.sin((-8.44797245e9*p[0]/freq) + p[1]) - numpy.sin(y))
      
      plotrm          = lambda RM, wav: numpy.mod( (1.*RM*wav*wav) +1.0*pi, 2.*pi)  -1.0*pi # notice the factor of 2
      fitfuncfastplot    = lambda p, freq: numpy.mod((2.*pi*p[0]*freq) - (1.*8.44797245e9*p[1]/freq) + (p[2])+ 1.0*pi, 2.*pi) -1.0*pi 


      
      phase       = numpy.copy(phases)      
       
      pi = numpy.pi
      chi_old=1e9 
      
      if len(freq) > 10: 
      
        # FIND INTIAL GUESS
        for doffset in numpy.arange(-pi-0.1,pi+0.1,0.1):
         for dclock in numpy.arange(-8.0,8.0,0.01):
            phase_model = numpy.mod ( (2.*pi*dclock*freq) + doffset, 2*pi)  # NOTE THE *2 to use rr+ll instead of 0.5*(rr+ll)
            phase_data  = numpy.mod (phase, 2*pi)
            angle       = pi - numpy.abs(numpy.abs(phase_model - phase_data) - pi)
            chi = numpy.sum(angle)

            if chi < chi_old:
              chi_old = chi
              fitguess     = [dclock,doffset]
              #print 'Better fit 1', dclock, doffset

        fitguess_1 = numpy.copy(fitguess)
        #print 'iter 1', fitguess

        for doffset in numpy.arange(fitguess_1[1]-0.2,fitguess_1[1]+ 0.2,0.01):
         for dclock in numpy.arange(fitguess_1[0]-0.1,fitguess_1[0]+ 0.1,0.001):
            phase_model = numpy.mod ( (2.*pi*dclock*freq) + doffset, 2*pi)  # NOTE THE *2 to use rr+ll instead of 0.5*(rr+ll)
            phase_data  = numpy.mod (phase, 2*pi)
            angle       = pi - numpy.abs(numpy.abs(phase_model - phase_data) - pi)
            chi = numpy.sum(angle)

            if chi < chi_old:
              chi_old = chi
              fitguess     = [dclock,doffset]
              #print 'Better fit 2', dclock, doffset
          
        freq = freq.astype(numpy.float64)
        phase= phase.astype(numpy.float64) #epsfcn=1e-7
      
      	#print 'hello', fitguess
        

        fitresult, success   = scipy.optimize.leastsq(par2slope, [fitguess[0], fitguess[1]], args=(freq, phase),maxfev=10000)
        
        #print 'ok',fitresult


      else:
         print 'No valid data found'
         fitresult = [0.0,0.0]

            
      
      
   
      show_plot = False
      if show_plot:  
        
        #if len(fitresult ==2): 
        #  fitresult =[fitresult[0], fitresult[1],0]
        #  print 'Here'
        #  #fitresult = fitresult #   [fitguess[0],  fitguess[1], 0]

        #fitresult = [[fitresultguess[0],fitresultguess[0],0.0]]

        matplotlib.pyplot.plot(freq, numpy.mod(phase + 1.0*pi, 2.*pi) -1.0*pi, 'or' )
        matplotlib.pyplot.plot(freq, numpy.mod(phase_diff + 1.0*pi, 2.*pi) -1.0*pi , '.', color='purple' )                         

        #TEC   = numpy.mod((-8.44797245e9*(2.*fitresult[1])/freq)+numpy.pi, 2*pi) - pi   # notice factor of 2 because rr+ll
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
      
      return [fitresult[0], fitresult[1]]


def fit_dTEC_dclock_doffset(phases, freq):
      c = 2.99792458e8
      freq_old = numpy.copy(freq)
      # NOTE THE *2 to use rr+ll instead of 0.5*(rr+ll)
      par3complex     = lambda p, freq, y: abs(numpy.cos((2.*pi*p[0]*freq) - (1.*8.44797245e9*p[1]/freq) + p[2]) - numpy.cos(y)) + abs(numpy.sin((2.*pi*p[0]*freq) - (1.*8.44797245e9*p[1]/freq) + p[2]) - numpy.sin(y))
      par2complex     = lambda p, freq, y: abs(numpy.cos((4.*pi*p[0]*freq) - (2.*8.44797245e9*p[1]/freq)) - numpy.cos(y)) + abs(numpy.sin((4.*pi*p[0]*freq) - (2.*8.44797245e9*p[1]/freq)) - numpy.sin(y))
      rmwavcomplex    = lambda RM, wav, y: abs(numpy.cos(2.*RM[0]*wav*wav)  - numpy.cos(y)) + abs(numpy.sin(2.*RM[0]*wav*wav)  - numpy.sin(y))
      par2slope       = lambda p, freq, y: abs(numpy.cos((-8.44797245e9*p[0]/freq) + p[1]) - numpy.cos(y)) + abs(numpy.sin((-8.44797245e9*p[0]/freq) + p[1]) - numpy.sin(y))
      
      plotrm          = lambda RM, wav: numpy.mod( (1.*RM*wav*wav) +1.0*pi, 2.*pi)  -1.0*pi # notice the factor of 2
      fitfuncfastplot    = lambda p, freq: numpy.mod((2.*pi*p[0]*freq) - (1.*8.44797245e9*p[1]/freq) + (p[2])+ 1.0*pi, 2.*pi) -1.0*pi 


      
      phase       = numpy.copy(phases)      
       
      pi = numpy.pi
      chi_old=1e9 
      
      if len(freq) > 10: 
        doffset = 0.0
        # FIND INTIAL GUESS
	# FIND INTIAL GUESS
	for dTEC in numpy.arange(-1.0,1.0, 0.01):
	 for dclock in numpy.arange(-150e-9,150e-9,5e-9):
	  for doffset in numpy.arange(-pi-0.2,pi+0.2,0.2):
            phase_model = numpy.mod ( (2.*pi*dclock*freq) - (8.44797245e9*dTEC/freq) + doffset, 2*pi)  # NOTE THE *2 to use rr+ll instead of 0.5*(rr+ll)
            phase_data  = numpy.mod (phase, 2*pi)
	    angle       = pi - numpy.abs(numpy.abs(phase_model - phase_data) - pi)
	    chi = numpy.sum(angle)

	    if chi < chi_old:
	      chi_old = chi
	      fitguess     = [dclock,dTEC,doffset]
	      #print 'Better fit', dclock, dTEC

	fitguess_1 = numpy.copy(fitguess)
	#print 'iter 1', fitguess

	for dTEC in numpy.arange(fitguess_1[1]-0.02,fitguess_1[1]+0.02, 0.004):
	 for dclock in numpy.arange(fitguess_1[0]-8e-9,fitguess_1[0]+8e-9,2e-9):
	  for doffset in numpy.arange(doffset+0.3,doffset-0.3,0.1):
            phase_model = numpy.mod ( (2.*pi*dclock*freq) - (8.44797245e9*dTEC/freq)+doffset, 2*pi)  # NOTE THE *2 to use rr+ll instead of 0.5*(rr+ll)
            phase_data  = numpy.mod (phase, 2*pi)
	    angle       = pi - numpy.abs(numpy.abs(phase_model - phase_data) - pi)
	    chi = numpy.sum(angle)

	    if chi < chi_old:
	      chi_old = chi
	      fitguess     = [dclock,dTEC,doffset]
	      #print 'Better fit', dclock, dTEC

          
        freq = freq.astype(numpy.float64)
        phase= phase.astype(numpy.float64) #epsfcn=1e-7
      
      	#print 'hello', fitguess
        

        #fitresult, success   = scipy.optimize.leastsq(par3complex, [fitguess[0], fitguess[1],0], args=(freq, phase),maxfev=10000)
        #fitresult, success   = scipy.optimize.leastsq(par3complex, [0,0,0], args=(freq, phase),maxfev=10000)
        fitresult, success   = scipy.optimize.leastsq(par3complex, [fitguess[0], fitguess[1],fitguess[2]], args=(freq, phase),maxfev=10000)
        
        #print 'ok',fitresult


      else:
         print 'No valid data found'
         fitresult = [0.0,0.0]

            
      
      
   
      show_plot = False
      if show_plot:  
        
        #if len(fitresult ==2): 
        #  fitresult =[fitresult[0], fitresult[1],0]
        #  print 'Here'
        #  #fitresult = fitresult #   [fitguess[0],  fitguess[1], 0]

        #fitresult = [[fitresultguess[0],fitresultguess[0],0.0]]

        matplotlib.pyplot.plot(freq, numpy.mod(phase + 1.0*pi, 2.*pi) -1.0*pi, 'or' )
        #matplotlib.pyplot.plot(freq, numpy.mod(phase_diff + 1.0*pi, 2.*pi) -1.0*pi , '.', color='purple' )                         

        TEC   = numpy.mod((-8.44797245e9*(fitresult[1])/freq)+numpy.pi, 2*pi) - pi   # notice factor of 2 because rr+ll
        Clock = numpy.mod((2.*numpy.pi*fitresult[0]*freq )+numpy.pi, 2*pi) - pi   # notice factor of 2 because rr+ll
    
        phase_total = (2.*numpy.pi*fitresult[0]*freq)+(-8.44797245e9*(fitresult[1])/freq)+fitresult[2]
        residual    = numpy.mod(phase-phase_total+pi,2.*pi)-pi
        matplotlib.pyplot.plot(freq, residual, '.', color='purple')
    
        idxl = int(min(freq_old)/1e4) 
        idxh = int(max(freq_old)/1e4) 
        bigfreqaxistmp = range(idxl, idxh)
        bigfreqaxis    =  numpy.array([float(i) for i in bigfreqaxistmp])
        bigfreqaxis    = bigfreqaxis*1e4

        matplotlib.pyplot.plot (bigfreqaxis, fitfuncfastplot(fitresult, bigfreqaxis[:]), "r-")
        #matplotlib.pyplot.plot (bigfreqaxis, plotrm(fitresultrm_wav, c/bigfreqaxis[:]), "-", color='purple')
        #matplotlib.pyplot.plot (freq, fitfuncfastplot(fitresult, bigfreqaxis[:]), "r-")

        matplotlib.pyplot.plot(freq, Clock,  ',g') 
        matplotlib.pyplot.plot(freq, TEC,  ',b') 
        matplotlib.pyplot.xlabel('freq')
        matplotlib.pyplot.ylabel('phase')

        matplotlib.pyplot.show()
      
      phase_total = (2.*numpy.pi*fitresult[0]*freq)+(-8.44797245e9*(fitresult[1])/freq)+fitresult[2]
      
      #angle       = pi - numpy.abs(numpy.abs(phase - phase_total) - pi)
      #chi = numpy.sum(angle)
	    
      residual    = numpy.mod(phase-phase_total+pi,2.*pi)-pi
      #print residual
      chi         = numpy.sum(numpy.abs(residual))/len(residual)
      return [fitresult[0], fitresult[1],fitresult[2], chi]




#import lofar.expion.fitting as fitting
pi = numpy.pi
c  = 2.99792458e8


ionmodel = h5parm('scalarphase.h',readonly=True)
solsetNames = ionmodel.getSolsets()
for solsetName in solsetNames:
    print solsetName
solset = ionmodel.getSolset(solsetName)

soltabs = ionmodel.getSoltabs('sol000')
scalarphasetab = ionmodel.getSoltab('sol000','scalarphase000')
anttab = ionmodel.getAnt('sol000')



print scalarphasetab.freq[:]
#print scalarphasetab.time[:]
#print scalarphasetab.val[:]
#print scalarphasetab.ant[:]

sphase = scalarphasetab.val[:]
print np.shape(sphase)


station = 52
ref_station = 0
time_id   = 710

print scalarphasetab.ant[station]

phase = sphase[0,station,:,:] - sphase[0,ref_station,:,:]

matplotlib.pyplot.plot(scalarphasetab.freq[:]/1e6, np.mod(phase[:,710:715]+pi, 2.*pi) -pi, 'o')
matplotlib.pyplot.xlabel('freq [MHz]')
matplotlib.pyplot.ylabel('phase [rad]')
matplotlib.pyplot.ylim(-1.*pi,pi)
matplotlib.pyplot.show()


#sys.exit()

phase = numpy.copy(np.rot90(phase))

#                                                               time, freq
#phase = numpy.copy(scipy.ndimage.filters.minimum_filter(phase, (0,0)))


matplotlib.pyplot.subplot(1,2,0)
matplotlib.pyplot.imshow(np.cos(phase), vmin=-1.0, vmax=1.0, aspect='auto')
matplotlib.pyplot.subplot(1,2,1)
matplotlib.pyplot.imshow(np.sin(phase), vmin=-1.0, vmax=1.0, aspect='auto')
matplotlib.pyplot.xlabel('block')
matplotlib.pyplot.ylabel('time')
#matplotlib.pyplot.show()

########################




refantenna_id = 0
source_id     = 0  # source ID in global_db (usually 0)
ncpus         = 24*2

clockfit = 0.*numpy.copy(scalarphasetab.time[:])
phaseoffset = 0.*numpy.copy(scalarphasetab.time[:])
TECfit = 0.*numpy.copy(scalarphasetab.time[:])
chifit = 0.*numpy.copy(scalarphasetab.time[:])

clockarray = numpy.zeros([len(scalarphasetab.time),len(scalarphasetab.ant)])
phaseoffsetarray = numpy.zeros([len(scalarphasetab.time),len(scalarphasetab.ant)])
TECarray = numpy.zeros([len(scalarphasetab.time),len(scalarphasetab.ant)])
chiarray = numpy.zeros([len(scalarphasetab.time),len(scalarphasetab.ant)])

print 'REF STATION:', scalarphasetab.ant[refantenna_id]
print '# TIMESLOTS  ', len(scalarphasetab.time)
print '# FREQUENCIES', len(scalarphasetab.freq)

ppservers = ()
# Creates jobserver with ncpus workers
job_server = pp.Server(ncpus, ppservers=ppservers)
print "Starting pp with", job_server.get_ncpus(), "workers"

phases_all = numpy.copy(scalarphasetab.val)
start_time_id = 0
stop_time_id  = len(scalarphasetab.time)

print numpy.shape(phases_all)

phases_freq = numpy.copy(scalarphasetab.freq[:])

for antenna_id in range(1,  len(scalarphasetab.ant[:])) :
#for antenna_id in range(1,  3) :


 if antenna_id != refantenna_id:
   stationspos      =  anttab[scalarphasetab.ant[refantenna_id]] - anttab[scalarphasetab.ant[antenna_id]]
   distance_station = numpy.sqrt(stationspos[0]**2 + stationspos[1]**2 + stationspos[2]**2)
   print 'Distance stations to reference station', distance_station/1e3, ' km'

   jobs = []

   for time_id in range(start_time_id,stop_time_id):
   #for time_id in range(0,10):
     
          
     
     phases = numpy.copy((phases_all[source_id,antenna_id,:,time_id] -phases_all[source_id,refantenna_id,:,time_id])) 
      

     #jobs.append(job_server.submit(fit_dclock_doffset,(phases,phases_freq),(), ("numpy","scipy.optimize","matplotlib.pyplot",)))
     jobs.append(job_server.submit(fit_dTEC_dclock_doffset,(phases,phases_freq),(), ("numpy","scipy.optimize","matplotlib.pyplot",)))

     print 'Submitting job #', time_id
     #print phases, phases_freq
 

   i= 0
   for job in jobs:
     fitresult = job()
     #print job()
     #print numpy.shape(clockfit)
     clockfit[start_time_id+i]  = fitresult[0]
     TECfit[start_time_id+i]  =  fitresult[1]
     phaseoffset[start_time_id+i]  =  fitresult[2]
     chifit[start_time_id+i]  =  fitresult[3]

     clockarray[start_time_id+i,antenna_id] = fitresult[0]
     TECarray[start_time_id+i,antenna_id] = fitresult[1]
     phaseoffsetarray[start_time_id+i,antenna_id] = fitresult[2] 
     chiarray[start_time_id+i,antenna_id] = fitresult[3]
     print 'TIME_ID', start_time_id+i, 'FIT (dclock,offset)', fitresult
     i = i + 1


numpy.save('ScPclock.npy', clockarray)
numpy.save('ScPoffset.npy', phaseoffsetarray)
numpy.save('ScPTEC.npy', TECarray)
numpy.save('ScPchi.npy', chiarray)


