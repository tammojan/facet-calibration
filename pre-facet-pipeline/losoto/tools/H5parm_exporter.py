#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This tool is used to convert an H5parm file to parmdb format by writing to
# existing parmdb instrument table(s).
#
# It handles Gain/DirectionalGain/RotationAngle/CommonRotationAngle/CommonScalarPhase solution types.
_author = "Francesco de Gasperin (fdg@hs.uni-hamburg.de), David Rafferty (drafferty@hs.uni-hamburg.de)"

import sys, os, glob, re
import numpy as np
import shutil
import progressbar
import logging
import pyrap.tables as pt
import lofar.parmdb
import losoto._version
import losoto._logging
from losoto.h5parm import h5parm, solWriter, solFetcher


def parmdbToAxes(solEntry):
    """
    Extract the information written as a string in the parmdb format
    """
    pol = None; pol1 = None; pol2 = None;
    dir = None; ant = None; parm = None

    thisSolType = solEntry.split(':')[0]

    # For CommonRotationAngle assuming [CommonRotationAngle:ant]
    if thisSolType == 'CommonRotationAngle':
        thisSolType, ant = solEntry.split(':')
        dir = 'pointing'

    # For RotationAngle assuming [RotationAngle:ant:sou]
    elif thisSolType == 'RotationAngle':
        thisSolType, ant, dir = solEntry.split(':')

    # For CommonScalarPhase assuming [CommonScalarPhase:ant]
    elif thisSolType == 'CommonScalarPhase':
        thisSolType, ant = solEntry.split(':')

    # For ScalarPhase assuming [ScalarPhase:ant:sou]
    elif thisSolType == 'ScalarPhase':
        thisSolType, ant, dir = solEntry.split(':')

    # For Gain assuming [Gain:pol1:pol2:parm:ant]
    elif thisSolType == 'Gain':
        thisSolType, pol1, pol2, parm, ant = solEntry.split(':')
        dir = 'pointing'

    # For DirectionalGain assuming [DirecitonalGain:pol1:pol2:parm:ant:sou]
    elif thisSolType == 'DirectionalGain':
        thisSolType, pol1, pol2, parm, ant, dir = solEntry.split(':')

    else:
        logging.error('Unknown solution type "'+thisSolType+'". Ignored.')

    if pol1 != None and pol2 != None:
        if pol1 == '0' and pol2 == '0': pol = 'XX'
        if pol1 == '1' and pol2 == '0': pol = 'YX'
        if pol1 == '0' and pol2 == '1': pol = 'XY'
        if pol1 == '1' and pol2 == '1': pol = 'YY'

    if pol != None:
        pol = re.escape(pol)
    if dir != None:
        dir = re.escape(dir)
    if ant != None:
        ant = re.escape(ant)
    return pol, dir, ant, parm


def getSoltabFromSolType(solType, solTabs, parm='ampl'):
    """Return a list of solution tables that corresponds to the type given

    solType - string defining solution type. E.g., "DirectionalGain"
    solTabs - solution tables returned by h5parm.getSoltabs()
    parm - parm to return if solType is "Gain": 'ampl'/'real' or
           'phase'/'imag'

    The soltab parmdb_type attribute is used as an additional filter to
    to distinguish multiple possible matches. If it is not available, all
    matches are returned.
    """
    solTabList = []
    if parm != None:
        parm = parm.lower()

    for name, st in solTabs.iteritems():
        # Handle gain table separately, as we need to distinguish ampl and phase
        if solType == 'DirectionalGain' or solType == 'Gain':
            if (parm == 'ampl' or parm == 'real') and st._v_title == 'amplitude':
                if hasattr(st._v_attrs, 'parmdb_type'):
                    if st._v_attrs['parmdb_type'] is not None:
                        if solType in st._v_attrs['parmdb_type'].split(', '):
                            solTabList.append(st)
                    else:
                        solTabList.append(st)
                else:
                    solTabList.append(st)
            elif (parm == 'phase' or parm == 'imag') and st._v_title == 'phase':
                if hasattr(st._v_attrs, 'parmdb_type'):
                    if st._v_attrs['parmdb_type'] is not None:
                        if solType in st._v_attrs['parmdb_type'].split(', '):
                            solTabList.append(st)
                    else:
                        solTabList.append(st)
                else:
                    solTabList.append(st)
        else:
            if hasattr(st._v_attrs, 'parmdb_type'):
                if st._v_attrs['parmdb_type'] is not None:
                    if solType in st._v_attrs['parmdb_type'].split(', '):
                        solTabList.append(st)
                else:
                    if (solType == 'RotationAngle' or solType == 'CommonRotationAngle') and st._v_title == 'rotation':
                        solTabList.append(st)
                    elif (solType == 'CommonScalarPhase' or solType == 'ScalarPhase') and st._v_title == 'scalarphase':
                        solTabList.append(st)
            else:
                if (solType == 'RotationAngle' or solType == 'CommonRotationAngle') and st._v_title == 'rotation':
                    solTabList.append(st)
                elif (solType == 'CommonScalarPhase' or solType == 'ScalarPhase') and st._v_title == 'scalarphase':
                    solTabList.append(st)

    if len(solTabList) == 0:
        return None
    else:
        return solTabList


def makeTECparmdb(H, solset, TECsolTab, timewidths, freq, freqwidth):
    """Returns TEC screen parmdb parameters

    H - H5parm object
    solset - solution set with TEC screen parameters
    TECsolTab = solution table with tecscreen values
    timewidths - time widths of output parmdb
    freq - frequency of output parmdb
    freqwidth - frequency width of output parmdb
    """
    from pylab import pinv
    global ipbar, pbar

    station_dict = H.getAnt(solset)
    station_names = station_dict.keys()
    station_positions = station_dict.values()
    source_dict = H.getSou(solset)
    source_names = source_dict.keys()
    source_positions = source_dict.values()

    tec_sf = solFetcher(TECsolTab)
    tec_screen, axis_vals = tec_sf.getValues()
    times = axis_vals['time']
    beta = TECsolTab._v_attrs['beta']
    r_0 = TECsolTab._v_attrs['r_0']
    height = TECsolTab._v_attrs['height']
    order = TECsolTab._v_attrs['order']
    pp = tec_sf.t.piercepoint

    N_sources = len(source_names)
    N_times = len(times)
    N_freqs = 1
    N_stations = len(station_names)
    N_piercepoints = N_sources * N_stations

    freqs = freq
    freqwidths = freqwidth
    parms = {}
    v = {}
    v['times'] = times
    v['timewidths'] = timewidths
    v['freqs'] = freqs
    v['freqwidths'] = freqwidths

    for station_name in station_names:
        for source_name in source_names:

            v['values'] = np.zeros((N_times, N_freqs), dtype=np.double)
            parmname = 'Piercepoint:X:%s:%s' % (station_name, source_name)
            parms[parmname] = v.copy()

            v['values'] = np.zeros((N_times, N_freqs), dtype=np.double)
            parmname = 'Piercepoint:Y:%s:%s' % (station_name, source_name)
            parms[parmname] = v.copy()

            v['values'] = np.zeros((N_times, N_freqs), dtype=np.double)
            parmname = 'Piercepoint:Z:%s:%s' % (station_name, source_name)
            parms[parmname] = v.copy()

            v['values'] = np.zeros((N_times, N_freqs), dtype=np.double)
            parmname = 'TECfit_white:%s:%s' % (station_name, source_name)
            parms[parmname] = v.copy()

            v['values'] = np.zeros((N_times, N_freqs), dtype=np.double)
            parmname = 'TECfit_white:0:%s:%s' % (station_name, source_name)
            parms[parmname] = v.copy()

            v['values'] = np.zeros((N_times, N_freqs), dtype=np.double)
            parmname = 'TECfit_white:1:%s:%s' % (station_name, source_name)
            parms[parmname] = v.copy()

    for k in range(N_times):
        D = np.resize(pp[k, :, :], (N_piercepoints, N_piercepoints, 3))
        D = np.transpose(D, ( 1, 0, 2 )) - D
        D2 = np.sum(D**2, axis=2)
        C = -(D2 / (r_0**2))**(beta / 2.0) / 2.0
        tec_fit_white = np.dot(pinv(C),
            tec_screen[:, k, :].reshape(N_piercepoints))
        pp_idx = 0
        for src, source_name in enumerate(source_names):
            for sta, station_name in enumerate(station_names):

                parmname = 'Piercepoint:X:%s:%s' % (station_name, source_name)
                parms[parmname]['values'][k, 0] = pp[k, pp_idx, 0]

                parmname = 'Piercepoint:Y:%s:%s' % (station_name, source_name)
                parms[parmname]['values'][k, 0] = pp[k, pp_idx, 1]

                parmname = 'Piercepoint:Z:%s:%s' % (station_name, source_name)
                parms[parmname]['values'][k, 0] = pp[k, pp_idx, 2]

                parmname = 'TECfit_white:%s:%s' % (station_name, source_name)
                parms[parmname]['values'][k, 0] = tec_fit_white[pp_idx]

                parmname = 'TECfit_white:0:%s:%s' % (station_name, source_name)
                parms[parmname]['values'][k, 0] = tec_fit_white[pp_idx]

                parmname = 'TECfit_white:1:%s:%s' % (station_name, source_name)
                parms[parmname]['values'][k, 0] = tec_fit_white[pp_idx]

                pp_idx += 1
        pbar.update(ipbar)
        ipbar += 1

    time_start = times[0] - timewidths[0]/2
    time_end = times[-1] + timewidths[-1]/2

    v['times'] = np.array([(time_start + time_end) / 2])
    v['timewidths'] = np.array([time_end - time_start])

    v_r0 = v.copy()
    v_r0['values'] = np.array(r_0, dtype=np.double, ndmin=2)
    parms['r_0'] = v_r0

    v_beta = v.copy()
    v_beta['values'] = np.array(beta, dtype=np.double, ndmin=2)
    parms['beta'] = v_beta

    v_height = v.copy()
    v_height['values'] = np.array(height, dtype=np.double, ndmin=2)
    parms['height'] = v_height

    return parms


if __name__=='__main__':
    # Options
    import optparse
    opt = optparse.OptionParser(usage='%prog <H5parm filename> <input globaldb/SB filename>\n'+
        _author, version='%prog '+losoto._version.__version__)
    opt.add_option('-v', '--verbose', help='Go VeRbOsE!',
        action='store_true', default=False)
    opt.add_option('-s', '--solset', help='Name of solution set to export '
        '(default=sol000)', type='string', default='sol000')
    opt.add_option('-o', '--outfile', help='Filename of globaldb/SB to export parmdb to '
        '(default=input globaldb/SB filename)', type='string', default=None)
    opt.add_option('-r', '--root', help='Root string to prepend to input parmdb '
        'instrument directories to make the output parmdb directories '
        '(default=solution-set name)', type='string', default=None)
    opt.add_option('-t', '--soltab', help='Solution tables to export; e.g., '
        '"amplitude000, phase000" (default=all)', type='string', default='all')
    opt.add_option('-i', '--instrument', help='Name of the instrument table '
        '(default=instrument*)', type='string', default='instrument*')
    opt.add_option('-c', '--clobber', help='Clobber exising files '
        '(default=False)', action='store_true', default=False)
    (options, args) = opt.parse_args()
    global ipbar, pbar

    # Check options
    if len(args) != 2:
        opt.print_help()
        sys.exit()
    if options.verbose: losoto._logging.setLevel("debug")

    # Check input H5parm file
    h5parmFile = args[0]
    if not os.path.exists(h5parmFile):
        logging.critical('Input H5parm file not found.')
        sys.exit(1)
    logging.info("Input H5parm filename = "+h5parmFile)

    # Open the h5parm file and get solution set names
    h5parm_in = h5parm(h5parmFile, readonly = True)
    solsetNames = h5parm_in.getSolsets()

    # Check input parmdb file
    globaldbFile = args[1]
    if not os.path.exists(globaldbFile):
        logging.critical('Input globaldb/SB file not found.')
        sys.exit(1)
    logging.info("Input globaldb/SB filename = "+globaldbFile)

    # Check input solution set name
    solsetName = options.solset
    if solsetName not in solsetNames:
        logging.critical('The solution set "'+solsetName+'" was not found in input H5parm file.')
        sys.exit(1)
    logging.info("Solution set name = "+solsetName)
    solset = h5parm_in.getSolset(solsetName)

    # Make output parmdb directory if needed
    out_globaldbFile = options.outfile
    if out_globaldbFile is None:
        out_globaldbFile = globaldbFile
    if not os.path.exists(out_globaldbFile):
        os.mkdir(out_globaldbFile)
    logging.info("Output globaldb/SB filename = "+out_globaldbFile)

    # Check output root
    outroot = options.root
    if outroot is None:
        outroot = solsetName

    # Make a list of all available instrument tables (only 1 for a standard MS)
    instrumentdbFiles = [ instrumentdbFile for instrumentdbFile in \
        glob.glob(os.path.join(globaldbFile, options.instrument)) \
        if os.path.isdir(instrumentdbFile) ]
    if len(instrumentdbFiles) == 0:
        logging.critical('No parmdb table(s) found in input globaldb/SB file.')
        sys.exit(1)
    instrumentdbFiles.sort()

    # Find solution table types using the first instrumentdb
    # TODO: is there a better solution which check all the instrumentdbs?
    pdb = lofar.parmdb.parmdb(instrumentdbFiles[0])
    solTypes = list(set(x[0] for x in  (x.split(":") for x in pdb.getNames())))
    solTabs = h5parm_in.getSoltabs(solset)
    if options.soltab != 'all':
        soltabs_to_use = [s.strip() for s in options.soltab.split(',')]
        logging.info('Using solution tables: {0}'.format(soltabs_to_use))
        solTabs_filt = {}
        for s, v in solTabs.iteritems():
            if s in soltabs_to_use:
                solTabs_filt[s] = v
        for s in soltabs_to_use:
            if s not in solTabs_filt.keys():
                logging.warning('Solution table {0} not found in input H5parm file.'.format(s))
        solTabs = solTabs_filt
    if len(solTabs) == 0:
        logging.critical('No solution tables found in input H5parm file')
        sys.exit(1)
    pdbSolTypes = solTypes[:]
    for solType in pdbSolTypes:
        solTabList = getSoltabFromSolType(solType, solTabs)
        if solTabList is None:
            logging.warning("Solution type {0} not found in solution set {1}. Skipping.".format(solType, solsetName))
            solTypes.remove(solType)

    # Look for tecscreen solution table in the solset. If
    # found, add to solTypes
    st_tec = None
    for name, st in solTabs.iteritems():
        if st._v_title == 'tecscreen':
            st_tec = st
    if st_tec is not None:
        solTypes.append('TECScreen')
    solTypes = list(set(solTypes))
    logging.info('Found solution types in input parmdb and H5parm: '+', '.join(solTypes))

    # For each solType, select appropriate solutions and construct
    # the dictionary to pass to pdb.addValues()
    len_sol = {}
    for solType in solTypes:
        if solType != 'TECScreen':
            len_sol[solType] = len(pdb.getNames(solType+':*'))
        else:
            tec_sf = solFetcher(st_tec)
            N_times = tec_sf.getAxisLen(axis='time')
            len_sol[solType] = N_times

    for instrumentdbFile in instrumentdbFiles:
        out_instrumentdbFile = out_globaldbFile + '/' + outroot + '_' + instrumentdbFile.split('/')[-1]
        logging.info('Filling '+out_instrumentdbFile+':')

        # Remove existing instrumentdb (if clobber) and create new one
        if os.path.exists(out_instrumentdbFile):
            if options.clobber:
                shutil.rmtree(out_instrumentdbFile)
            else:
                logging.critical('Output instrumentdb file exists and '
                    'clobber = False.')
                sys.exit(1)
        pdb_out = lofar.parmdb.parmdb(out_instrumentdbFile+'/', create=True)

        pbar = progressbar.ProgressBar(maxval=sum(len_sol.values())).start()
        ipbar = 0

        pdb_in = lofar.parmdb.parmdb(instrumentdbFile)

        # Add default values and steps
        DefValues = pdb_in.getDefValues()
        for k, v in DefValues.iteritems():
            pdb_out.addDefValues({k: pdb.makeDefValue(v.item(0))})
        pdb_out.setDefaultSteps(pdb_in.getDefaultSteps())

        for solType in solTypes:
            if len_sol[solType] == 0: continue

            if solType != 'TECScreen':
                solEntries = pdb_in.getNames(solType+':*')
                data = pdb_in.getValuesGrid(solType+':*')
                data_out = data.copy()
                for solEntry in solEntries:

                    pol, dir, ant, parm = parmdbToAxes(solEntry)
                    solTabList = getSoltabFromSolType(solType, solTabs, parm=parm)
                    if solTabList is None:
                        continue
                    if len(solTabList) > 1:
                        logging.warning('More than one solution table found in H5parm '
                            'matching parmdb entry "'+solType+'". Taking the first match.')
                    solTab = solTabList[0]
                    sf = solFetcher(solTab)

                    if pol == None and dir == None:
                        sf.setSelection(ant=ant)
                    elif pol == None and dir != None:
                        sf.setSelection(ant=ant, dir=dir)
                    elif pol != None and dir == None:
                        sf.setSelection(ant=ant, pol=pol)
                    else:
                        sf.setSelection(ant=ant, pol=pol, dir=dir)

                    # If needed, convert Amp and Phase to Real and Imag respectively
                    if parm == 'Real':
                        SolTabList = getSoltabFromSolType(solType, solTabs, parm='phase')
                        soltab_phase = SolTabList[0]
                        sf_phase = solFetcher(soltab_phase, ant=ant, pol=pol, dir=dir)
                        val_amp = sf.getValues()[0]
                        val_phase = sf_phase.getValues()[0]
                        val = val_amp * np.cos(val_phase)
                    elif parm == 'Imag':
                        SolTabList = getSoltabFromSolType(solType, solTabs, parm='ampl')
                        soltab_amp = SolTabList[0]
                        sf_amp = solFetcher(soltab_amp, ant=ant, pol=pol, dir=dir)
                        val_phase = sf.getValues()[0]
                        val_amp = sf_amp.getValues()[0]
                        val = val_amp * np.sin(val_phase)
                    else:
                        val = sf.getValues()[0]

                    # Match the frequency or frequencies of instrumentdb under
                    # consideration
                    sffreqs = sf.freq
                    freqs = data[solEntry]['freqs']
                    freq_list = [freq for freq in freqs if freq in sffreqs]
                    if len(freq_list) == 0:
                        for i in range(len(val.shape)):
                            freq_ind_list.append(slice(None))
                        freq_ind = tuple(freq_ind_list)
                    else:
                        freqAxisIdx = sf.getAxesNames().index('freq')
                        freq_ind = []
                        for i in range(len(val.shape)):
                            freq_ind.append(slice(None))
                        freq_ind[freqAxisIdx] = np.where(sffreqs == freq_list)
                        freq_ind = tuple(freq_ind)

                    shape = data_out[solEntry]['values'].shape
                    try:
                        data_out[solEntry]['values'] = val[freq_ind].T.reshape(shape)
                    except ValueError, err:
                        logging.critical('Mismatch between parmdb table and H5parm '
                        'solution table: Differing number of frequencies and/or times')
                        sys.exit(1)
                pbar.update(ipbar)
                ipbar += 1
            else:
                # Handle TECScreen parmdb
                #
                # Get timewidths, freqwidth and freq from first (non-TEC, phase)
                # solentry
                for nonTECsolType in pdbSolTypes:
                    if nonTECsolType != 'TECScreen' and 'Phase' in nonTECsolType:
                        break
                parmname = pdb_in.getNames(nonTECsolType+':*')[0]
                timewidths = pdb_in.getValuesGrid(parmname)[parmname]['timewidths']
                freqwidth = pdb.getValuesGrid(parmname)[parmname]['freqwidths'][0]
                freq = pdb.getValuesGrid(parmname)[parmname]['freqs'][0]
                data_out = makeTECparmdb(h5parm_in, solset, st_tec, timewidths, freq, freqwidth)

            pdb_out.addValues(data_out)

        pbar.finish()

    logging.info('Done.')






