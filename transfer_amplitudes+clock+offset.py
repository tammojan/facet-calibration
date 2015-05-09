import os
import numpy
import sys
import glob
import pyrap.tables as pt
import lofar.parmdb
#import pyrap.tables as pt
import argparse

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


def get_stored_info(storepath, calibrator = "3C295"):
    """
    Get the numpy data from the storage area
    """
    freqs_ampl = numpy.load('{path}/freqs_for_amplitude_array.npy'.format(path=storepath))
    amps_array  = numpy.load('{path}/{cal}_amplitude_array.npy'.format(path=storepath, cal=calibrator))
    clock_array = numpy.load('{path}/fitted_data_dclock_{cal}_1st.sm.npy'.format(path=storepath, cal=calibrator))

    freqs_phase = numpy.load('{path}/freqs_for_phase_array.npy'.format(path=storepath))
    phases_array  = numpy.load('{path}/{cal}_phase_array.npy'.format(path=storepath, cal=calibrator))   #freq, stat

    print "phases_array", numpy.shape(phases_array)
    print "amps_array", numpy.shape(amps_array)
    print "clock_array", numpy.shape(clock_array)
    print "freqs_phase", numpy.shape(freqs_phase)
    print "phases_array", numpy.shape(phases_array)

    return freqs_ampl, amps_array, clock_array, freqs_phase, phases_array


def transfer(mslist, storepath, output_parmdb=None, calibrator = "3C295"):
    """
    Transfer the clock, amplitudes and phase offsets to a list of MS entered as
    mslist.
    The data used for the transfer is loaded from storepath.
    The name of the output instrument is output_parmdb (optional with a default
    of 'instrument').
    The name of the calibrator can be entered as well.
    Note that the name of the input data files is hardcoded at the moment.
    """
    # Get the antenna information
    anttab          = pt.table(mslist[0] + '/ANTENNA')
    antenna_list    = anttab.getcol('NAME')
    anttab.close()

    # Load data
    freqs_ampl, amps_array, clock_array, freqs_phase, phases_array = get_stored_info(storepath, calibrator=calibrator)

    for ms in mslist:

        if output_parmdb is None:
            output_parmdb = os.path.join(ms, instrument)

        # Get the frequency information
        t          = pt.table(ms, readonly=False)
        freq_tab   = pt.table(ms + '/SPECTRAL_WINDOW')
        freq       = freq_tab.getcol('REF_FREQUENCY')
        t.close()
        print 'Frequency of MS', freq

        amp = give_ampl(amps_array, freqs_ampl, freq)
        phase = give_phase(phases_array, freqs_phase, freq)


        #SB = ms.split('_')[3]
        #newparmdb = 'Bootes_HBA_caltransfer_'+ SB + '.parmdb'
        os.system('mv %s %s_bck'%(output_parmdb, output_parmdb))

        pdb   = lofar.parmdb.parmdb('{path}/instrument'.format(path=storepath))
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

        pdbnew = lofar.parmdb.parmdb(output_parmdb, create=True)
        pdbnew.addValues(parms)
        pdbnew.flush()
        os.system("taql 'update " + output_parmdb + " set ENDX=1.e12'")
        os.system("taql 'update " + output_parmdb + " set STARTX=1.0'")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Transfer the clock, amplitudes and phase offsets')
    #parser.add_argument('ms', metavar='N', nargs='+', required=True, help='MS')
    parser.add_argument('ms', metavar='MS', help='MS')
    parser.add_argument('--store', required=True, help='Store path')
    parser.add_argument('--pdb', default="instrument", help='Name of the output parmdb')
    parser.add_argument('--calibrator', default="3C295", help='Name of the output parmdb')

    args = parser.parse_args()

    # Parse the input
    # Note that only one MS is acepted by now
    mslist = [os.path.abspath(args.ms)]
    storepath = os.path.abspath(args.store)
    output_parmdb = os.path.abspath(os.path.join(mslist[0], args.pdb))
    calibrator = args.calibrator

    print "mslist", mslist
    print "storepath", storepath
    print "calibrator", calibrator
    print "output_parmdb", output_parmdb

    transfer(mslist, storepath, output_parmdb=output_parmdb, calibrator=calibrator)
