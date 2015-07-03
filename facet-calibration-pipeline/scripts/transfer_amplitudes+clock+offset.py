#!/usr/bin/env python
import os
import numpy
import sys
import glob
import pyrap.tables as pt
import lofar.parmdb

def give_ampl(amps_in, freqs_amps, freqms):
  #given nearest element in freq_amps for freqms
  el  = min(range(len(freqs_amps)), key=lambda i: abs(freqs_amps[i]-freqms))
  amp = numpy.copy(amps_in[:,:,el,:]) # --> (antenna,time,pol) , el slices out freq axis 
  return amp

def give_phase(phases_in, freqs_phases, freqms):
  #given nearest element in freq_amps for freqms
  el  = min(range(len(freqs_phases)), key=lambda i: abs(freqs_phases[i]-freqms))
  phase = numpy.copy(phases_in[el,:]) # --> (antenna,time,freq) , el slices out freq axis 
  return phase

#mslist = sorted(glob.glob('L10000_SB026_uv_demix.MS'))


mslist = sys.argv[1]
#newparmdb = sys.argv[2]
storepath = sys.argv[3]
templatepath = sys.argv[4]
newparmdb = sys.argv[1] +'/' +sys.argv[2]
calibrator = sys.argv[5]
print mslist, newparmdb, storepath

mslist = [mslist]

anttab          = pt.table(mslist[0] + '/ANTENNA')
antenna_list    = anttab.getcol('NAME')
anttab.close()

freqs_ampl = numpy.load('{path}/freqs_for_amplitude_array.npy'.format(path=storepath))
amps_array  = numpy.load('{path}/{calname}_amplitude_array.npy'.format(path=storepath,calname=calibrator))
clock_array = numpy.load('{path}/fitted_data_dclock_{calname}_1st.sm.npy'.format(path=storepath,calname=calibrator))

freqs_phase = numpy.load('{path}/freqs_for_phase_array.npy'.format(path=storepath))
phases_array  = numpy.load('{path}/{calname}_phase_array.npy'.format(path=storepath,calname=calibrator))   #freq, stat

print numpy.shape(phases_array)
print numpy.shape(amps_array)
print numpy.shape(clock_array)



for ms in mslist:

 
 t          = pt.table(ms, readonly=False)
 freq_tab   = pt.table(ms + '/SPECTRAL_WINDOW')
 freq       = freq_tab.getcol('REF_FREQUENCY')
 t.close()
 print 'Frequency of MS', freq
 amp = give_ampl(amps_array, freqs_ampl, freq)
 
 phase = give_phase(phases_array, freqs_phase, freq)

 
 #SB = ms.split('_')[3]
 #newparmdb = 'Bootes_HBA_caltransfer_'+ SB + '.parmdb' 
 #os.system('rm -rf ' + newparmdb)
 
 pdb   = lofar.parmdb.parmdb(templatepath)
 parms = pdb.getValuesGrid("*")
 
 for antenna_id, antenna in enumerate(antenna_list):
   dummyvectmp = numpy.copy(parms['Gain:0:0:Real:' + antenna]['values'][:,0])
   dummyvec = 0.0*numpy.copy(dummyvectmp)
    

   amp_cal_00   = numpy.median(amp[antenna_id,:,0])
   amp_cal_11   = numpy.median(amp[antenna_id,:,1])
   
   # phase offset is stored as a single corr - to be applied to 11 (subtracted) or 00 (added)
   phase_cal_00   = 0.
   phase_cal_11   = -1*phase[antenna_id]

   print 'Amp value 00,11, Clock',amp_cal_00, amp_cal_11,  numpy.median(clock_array[:,antenna_id])*1e9
  
   
   real_00 = ((dummyvec) + (amp_cal_00))*numpy.cos(phase_cal_00)
   imag_00 = ((dummyvec) + (amp_cal_00))*numpy.sin(phase_cal_00)
    
   real_11 = ((dummyvec) + (amp_cal_11))*numpy.cos(phase_cal_11)
   imag_11 = ((dummyvec) + (amp_cal_11))*numpy.sin(phase_cal_11)

   
   parms['Gain:0:0:Real:' + antenna]['values'][:,0] = numpy.copy(real_00)
   parms['Gain:0:0:Imag:' + antenna]['values'][:,0] = numpy.copy(imag_00)        
   parms['Gain:1:1:Real:' + antenna]['values'][:,0] = numpy.copy(real_11)
   parms['Gain:1:1:Imag:' + antenna]['values'][:,0] = numpy.copy(imag_11)   
 
   parms['Clock:' + antenna]['values'][:,0] = dummyvec + numpy.median(clock_array[:,antenna_id]) 
   
 
 print parms['Gain:0:0:Real:CS001HBA0']['values'][:,0]
 print parms['Gain:0:0:Imag:CS001HBA0']['values'][:,0]
 print parms['Clock:RS208HBA']['values'][:,0]

 pdbnew = lofar.parmdb.parmdb(newparmdb, create=True)
 pdbnew.addValues(parms)
 pdbnew.flush()
 os.system("taql 'update " + newparmdb + " set ENDX=1.e12'")
 os.system("taql 'update " + newparmdb + " set STARTX=1.0'")
