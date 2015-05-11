import pyrap.tables as pt
import lofar.parmdb
import os
import numpy
from scipy import interpolate
import argparse

def create_merged_parmdb(ms, parmdb_a,parmdb_p, parmdb_t, parmdb_out,cellsizetime_a, cellsizetime_b):

    #pdb_pre = lofar.parmdb.parmdb(pre_apply_parmdb)
    pdb_a   = lofar.parmdb.parmdb(parmdb_a)
    pdb_p   = lofar.parmdb.parmdb(parmdb_p)
    pdb_t   = lofar.parmdb.parmdb(parmdb_t)

    #parms_pre = pdb_pre.getValuesGrid("*")
    parms_a   = pdb_a.getValuesGrid("*")
    parms_p   = pdb_p.getValuesGrid("*")
    parms_t   = pdb_t.getValuesGrid("*")


    os.system('rm -rf ' + parmdb_out)

    keynames_p = parms_p.keys()
    keynames_a = parms_a.keys()

    length = numpy.int(cellsizetime_b)
    for key in keynames_p:
        # copy over the CommonScalar phases, TEC, clock

        # just to get axis length
        tmp1  = numpy.copy(parms_t[key]['values'][:, 0])
        for idx in range(len(tmp1)):
            el =  numpy.float(idx)/numpy.float(length)
            el =  numpy.int(numpy.floor(el))

            parms_t[key]['values'][idx,0] = numpy.copy(parms_p[key]['values'][el,0])


    pol_list = ['0:0','1:1']
    gain     = 'Gain'
    anttab     = pt.table(ms + '/ANTENNA')
    antenna_list    = anttab.getcol('NAME')
    anttab.close()

    #length = numpy.int(cellsizetime_a/cellsizetime_b)
    length = numpy.int(cellsizetime_a)

    for pol in pol_list:
        for antenna in antenna_list:
            #real1 = numpy.copy(parms_pre[gain + ':' + pol + ':Real:'+ antenna]['values'][:, 0])
            real2 = numpy.copy(parms_a[gain + ':' + pol + ':Real:'+ antenna]['values'][:, 0])
            #imag1 = numpy.copy(parms_pre[gain + ':' + pol + ':Imag:'+ antenna]['values'][:, 0])
            imag2 = numpy.copy(parms_a[gain + ':' + pol + ':Imag:'+ antenna]['values'][:, 0])

            # just to get axis length
            tmp1  = numpy.copy(parms_t[gain + ':' + pol + ':Imag:'+ antenna]['values'][:, 0])

            #G1 = real1 + 1j*imag1
            G2 = real2 + 1j*imag2
            Gnew = G2

            for idx in range(len(tmp1)):
                el =  numpy.float(idx)/numpy.float(length)
                el =  numpy.int(numpy.floor(el))

                parms_t[gain + ':' + pol + ':Imag:'+ antenna]['values'][idx,0] = numpy.copy(numpy.imag(Gnew[el]))
                parms_t[gain + ':' + pol + ':Real:'+ antenna]['values'][idx,0] = numpy.copy(numpy.real(Gnew[el]))

    pdbnew = lofar.parmdb.parmdb(parmdb_out, create=True)
    pdbnew.addValues(parms_t)
    pdbnew.flush()


def default_path(instrument, ms, default):
    """
    Sets the default path to the parmdb if the argument is not set.
    """
    if instrument is None:
        out = os.path.join(ms, default)
    else:
        out = instrument
    return out


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Merge solutions from phase and amp parmdb') 
    parser.add_argument('ms', help="MS file")
    parser.add_argument('--instrument-amps', '-a', help="Amplitudes parmdb")
    parser.add_argument('--instrument-phase', '-p', help="Phases parmdb")
    parser.add_argument('--instrument-merged', '-m', help="Output parmdb")
    parser.add_argument('--instrument-template', '-t', help="Template parmdb")
    parser.add_argument('--asize', default=60, type=int, help="Amplitudes time size")
    parser.add_argument('--psize', default=1, type=int, help="Phases time size")
    
    args = parser.parse_args()
    
    # TODO: Check if the files exist and are directories
    
    # Get the full paths when needed
    instrument_amps = default_path(args.instrument_amps, args.ms, "instrument_amps1_smoothed")
    instrument_phase = default_path(args.instrument_phase, args.ms, "instrument_phase1")
    instrument_merged = default_path(args.instrument_merged, args.ms, "instrument_merged")
    instrument_template = default_path(args.instrument_template, args.ms, "instrument_template")

    
    debug = False
    if debug:
        print args.ms
        print instrument_amps
        print instrument_phase
        print instrument_template
        print instrument_merged
        print args.asize, type(args.asize)
        print args.psize, type(args.psize)
    else:
        create_merged_parmdb(args.ms, instrument_amps, instrument_phase,
                             instrument_template, instrument_merged,
                             args.asize, args.psize)
