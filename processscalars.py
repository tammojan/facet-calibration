import matplotlib.pyplot
import numpy
import sys, math, os
import scipy.signal
from losoto.h5parm import h5parm
import glob
import pyrap.tables as pt
import lofar.parmdb

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
    median_array  = scipy.signal.medfilt(sol,int(half_window*2.-1))

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



ionmodel = h5parm('scalarphase.h',readonly=True)

scalarphasetab = ionmodel.getSoltab('sol000','scalarphase000')


pi  = numpy.pi

clock = numpy.load('ScPclock.npy')
offset= numpy.load('ScPoffset.npy')
TEC   = numpy.load('ScPTEC.npy')
chi   = numpy.load('ScPchi.npy')

print clock[:,0]
print offset[:,0]
print TEC[:,0]


print numpy.shape(clock)

cut        = 3.0
freq       = 180e6
antenna_id = 56


# filter bad datapoints
phase_model = (2.*pi*clock[:,antenna_id]*freq)+(-8.44797245e9*TEC[:,antenna_id]/freq)+offset[:,antenna_id]
chi_vec     = chi[:,antenna_id]
idx         = numpy.where(chi_vec < (numpy.median(chi_vec)+cut*numpy.std(chi_vec)))

chi_vec     = chi_vec[idx]
phase_model = phase_model[idx]

# chi too good to be true
idx         = numpy.where((1./chi_vec) < (numpy.median(1./chi_vec)+1.5*cut*numpy.std(1./chi_vec)))
chi_vec     = chi_vec[idx]
phase_model = phase_model[idx]




matplotlib.pyplot.plot(numpy.mod(phase_model+pi, 2*pi)-pi, 'o')
matplotlib.pyplot.plot(chi_vec - 3.75, '.',)

matplotlib.pyplot.ylim(-pi-1.0,pi) # -1 to allow to show the chi_vec
matplotlib.pyplot.xlabel('time')
matplotlib.pyplot.title(scalarphasetab.ant[antenna_id])
#matplotlib.pyplot.show()




 ##################
 ##################

mslist          = sorted(glob.glob('A2256_SB3?0-3?9.2ch10s.ms'))
anttab          = pt.table(mslist[0] + '/ANTENNA')
antenna_list    = anttab.getcol('NAME')
anttab.close()


# weed out bad point and replace with nans so they get flagged
for antenna_id in range(1, len(clock[0,:])):
  # find the bad values, start with antenna 1 because of reference antenna 0
  chi_vec      = chi[:,antenna_id] 
  idx1         = numpy.where(chi_vec > (numpy.median(chi_vec)+cut*numpy.std(chi_vec)))
  idx2         = numpy.where((1./chi_vec) > (numpy.median(1./chi_vec)+1.5*cut*numpy.std(1./chi_vec)))

  clock[idx1,antenna_id] = numpy.nan
  clock[idx2,antenna_id] = numpy.nan
  
  TEC[idx1,antenna_id] = numpy.nan
  TEC[idx2,antenna_id] = numpy.nan  

  offset[idx1,antenna_id] = numpy.nan
  offset[idx2,antenna_id] = numpy.nan  
  



# fill the parmdb

for ms in mslist:
    
    newparmdb = ms + '/instrument_10blockphases'
    os.system('rm -rf ' + newparmdb)
    print 'Creating, ', newparmdb
    pdb   = lofar.parmdb.parmdb('/home/rvweeren/scripts/a2256_hba/instrumenttemplate_10blockphases')
    parms = pdb.getValuesGrid("*")

    for antenna_id, antenna in enumerate(scalarphasetab.ant):
      parms['Clock:' + antenna]['values'][:,0]             = clock[:,antenna_id]
      parms['TEC:' + antenna]['values'][:,0]               = TEC[:,antenna_id]
      parms['CommonScalarPhase:' + antenna]['values'][:,0] = offset[:,antenna_id]
      #print numpy.max(parms['CommonScalarPhase:' + antenna]['values'][:,0])
      print 'Filling values ', antenna
    
    print 'Writing the parmdb', newparmdb
    pdbnew = lofar.parmdb.parmdb(newparmdb, create=True)
    pdbnew.addValues(parms)
    pdbnew.flush()
    os.system("taql 'update " + newparmdb + " set ENDX=1.e12'")
    os.system("taql 'update " + newparmdb + " set STARTX=1.0'")
    
    # apply parmdb
    parset = '/home/rvweeren/scripts/a2256_hba/correct_10blockphases.parset'
    cmd = 'calibrate-stand-alone --parmdb-name instrument_10blockphases '
    os.system(cmd + ms + ' ' + parset + '&')
	      
	      
	      
sys.exit()


#matplotlib.pyplot.plot(numpy.mod(offset[:,52]+pi, 2*pi)-pi, '.')


for antenna_id in range(0, len(offset[0,:])):
#for antenna_id in range(10,11):
  
   print 'Cleaning up for antenna: ', antenna_id
   
  
   #slope[:,antenna_id] = median_window_filter(slope[:,antenna_id], 5, 3)
   #slope[:,antenna_id] = median_window_filter(slope[:,antenna_id], 5, 3)
   #slope[:,antenna_id] = running_median(slope[:,antenna_id] ,3) 
   
   real =   numpy.cos(offset[:,antenna_id])
   imag =   numpy.sin(offset[:,antenna_id])
   
   real = median_window_filter(real, 3, 3)
   imag = median_window_filter(imag, 3, 3)

   real = running_median(real, 3)
   imag = running_median(imag, 3)

   #print numpy.arctan2(imag[1], real[1]), offset[1,antenna_id]
   offset[:,antenna_id] = numpy.arctan2(imag, real)
 

#matplotlib.pyplot.plot(numpy.mod(offset[:,52]+pi, 2*pi)-pi, '.')

#matplotlib.pyplot.xlabel('time')
#matplotlib.pyplot.show()