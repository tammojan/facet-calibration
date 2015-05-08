#!/usr/bin/env python

"""
Plot solutions in a parmdb

Written by George Heald, modified by Oscar Martinez

For quality control plots, try:
solplot.py -q -a 0,0.06 -s RS\* -m -o gains-rs <MS>/instrument
solplot.py -q -a 0,0.06 -s CS\* -m -o gains-cs <MS>/instrument

version 0.2 11/02/2012: Initial version
version 0.3 22/02/2012: Added possibility to plot all correlations
                        Added handling of input data where paramdb is Ampl/Phase
                        Added handling of multi-channel solutions (they are averaged)
                        Added stats file
                        Phase is ignored if TEC data
                        Allow to plot even if only Phase is available
version 0.4 28/03/2012: Correct possible negative amplitudes and offset phase
version 0.5 15/08/2012: Try to eliminate failure when used in GNU screen
"""

import os
import lofar.parmdb as pdb
# pylab import statement is moved to main()
import optparse
import numpy

version_string = 'v0.5, 15 August 2012\nWritten by George Heald, modified by Oscar Martinez'
print 'solplot.py',version_string
print ''

# Some constants used during the script
PARAM_NAME_PATTERN = '%s:%s:%s:%s'
OFILE_PATTERN='%s\t%s\t%s\t%s\t%s\t%s\n'
AMPL_COORD = 'Ampl'
PHASE_COORD = 'Phase'
REAL_COORD = 'Real'
IMAG_COORD = 'Imag'
GAIN_TYPE = 'Gain'
TEC_TYPE = 'TEC'
AUTO_CORR = ['0:0', '1:1']
CROSS_CORR =['0:1','1:0']
CORR_COLOR_CODE = {'0:0':'r', '0:1':'g', '1:0':'y', '1:1':'b', }
CORR_COLOR_NAME = {'0:0':'Red', '0:1':'Green', '1:0':'Yellow', '1:1':'Blue', }

# This method gets the values of a paramdb for the specified key
def getParamDBData(p, pKey, cKey = 'values'):
    return p.getValuesGrid(pKey)[pKey][cKey]

# This method gets the (ampl,phase,averaged) of a paramdb for a specified station, ptype, corr, polar, last
# averaged indicates if the data has been averaged to a single channel 
def getAmplPhaseCorr(p, station, ptype, corr , polar, cords, last):
    valsamp = None
    valsphase = None
    if polar:
        if AMPL_COORD in cords:
            valsamp = getParamDBData(p, PARAM_NAME_PATTERN % (ptype,corr,AMPL_COORD,station))
        if ptype != TEC_TYPE and PHASE_COORD in cords:
            # In TEC type we do not get phase (we know they are real number)
            valsphase = getParamDBData(p, PARAM_NAME_PATTERN % (ptype,corr,PHASE_COORD,station))
            print 'TEST2'
    else:
        if ptype == TEC_TYPE:
            # In TEC type, since we do not have phase, we can get the real value as amplitude
            if REAL_COORD in cords:
                valsamp = getParamDBData(p, PARAM_NAME_PATTERN % (ptype,corr,REAL_COORD,station))
        else:
            if REAL_COORD in cords and IMAG_COORD in cords:
                valscorr = getParamDBData(p, PARAM_NAME_PATTERN%(ptype,corr,REAL_COORD,station))+1.j*getParamDBData(p, PARAM_NAME_PATTERN%(ptype,corr,IMAG_COORD,station))
                valsamp = numpy.abs(valscorr)
                valsphase = numpy.angle(valscorr)
    # We average the data if we detect the data have more than single channel
    averaged = False
    if (valsamp != None and valsamp.shape[1] > 1) or (valsphase != None and valsphase.shape[1] > 1):
        averaged = True
        if valsamp != None:
            valsamp = valsamp.mean(axis=1)
        if valsphase != None:
            valsphase = valsphase.mean(axis=1)
        
    if valsamp != None and not last:
        valsamp = valsamp[:-1]
    if valsphase != None and not last:
        valsphase = valsphase[:-1]

    if valsamp != None and valsphase != None: 
        valsphase[valsamp < 0.0] += numpy.pi
        valsamp = numpy.abs(valsamp)

    return (valsamp,valsphase,averaged)

# Get a dictionary which keys are the corrs with its values of ampl and phase
def getAmplPhase(p, station, ptype, corrs , polar, cords, last):
    valsamp = {}
    valsphase = {}
    for corr in corrs:
        (valsamp[corr], valsphase[corr], averaged) = getAmplPhaseCorr(p, station, ptype, corr , polar, cords, last)
    # We assume the times are the same for all the correlations and complex coord (they must be!)
    # We use first correlation and first complex coordinate
    times = getParamDBData(p, PARAM_NAME_PATTERN % (ptype, corrs[0], cords[0], station) , 'times')
    times -= times[0]
    if not last:
        times = times[:-1]
    return (valsamp, valsphase, times, averaged)

# Return True if at least one of the correlations contains valid data (Not None)
def isValid(vals):
    for corr in vals:
        if vals[corr] != None:
            return True
    return False

# MAIN script method
def main(opts,args):
    # Import pylab with a noninteractive backend if plotting to figure files
    if opts.out == '':
        import pylab
    else:
        import matplotlib as mpl
        mpl.use('Agg')
        import matplotlib.pyplot as pylab

    # Get the min and max y axis values
    l = opts.amplim.split(',')
    if l[0] == '': ampmin = None
    else: ampmin = float(l[0])
    if l[1] == '': ampmax = None
    else: ampmax = float(l[1])
    
    # Open the parmdb
    p = pdb.parmdb(args[0])
    
    if opts.title != '':
        title = opts.title
    else:
        title = args[0]
        
    ofile = None
    if opts.stats != '':
        if os.path.isfile(opts.stats):
            print 'Error: ' + opts.stats + ' already exists'
            exit()
        ofile = open(opts.stats, 'w')
        ofile.write(OFILE_PATTERN % ('#station', 'type:corr', 'num', 'mean', 'median', 'std'))
        
    # Get names (of the params) that are related to the selectes stations
    names = p.getNames('*:%s'%(opts.stations))
    if not len(names): print '\nNo stations matched the filter'
    ptype = GAIN_TYPE
    corrSet = set([])
    stationSet = set([])
    cordSet = set([])
    # We assume that at least
    for name in names:
        sname = name.split(':')
        # Add the station and correlation to its sets
        corrSet.add(sname[1] + ':' + sname[2])
        stationSet.add(sname[4])
        cordSet.add(sname[3])
        # Check that we are handling only Gain or TEC cases
        if sname[0] == TEC_TYPE: ptype=TEC_TYPE 
        elif sname[0] != GAIN_TYPE:
            print 'Error: parmdb is not '+GAIN_TYPE+' or '+TEC_TYPE
            exit()
    # We get if data is in polar format (preferable)
    if AMPL_COORD in cordSet or PHASE_COORD in cordSet:
        polar = True
    elif REAL_COORD in cordSet or IMAG_COORD in cordSet:
        polar = False
    else:
        print 'Error: paramdb contains ' + ','.join(cordSet)
        exit()

    # Check the found correlations: autocorrelations (0:0 and 1:1) are plotted always,
    # crosscorrelations are plotted only if the users asks for itand
    for corr in AUTO_CORR:
        if corr not in corrSet:
            print 'Warning: ' + corr + ' is not found in the paramdb'
    for corr in CROSS_CORR:
        if corr in corrSet:
            if not opts.allpol:
                print 'Warning: ' + corr + ' is found in paramdb but will be ignored'
                corrSet.remove(corr)
        elif opts.allpol:
            print 'Warning: ' + corr + ' was asked but is not found in paramdb'
    
    if not len(corrSet):
        print 'Error: no gain correlations to plot'
        exit()
    
    # Make lists of the set objects
    stations = sorted(list(stationSet))
    cords = list(cordSet)
    corrs = sorted(list(corrSet))
            
    # Get the reference stations
    if opts.reference == '':
        refstation = stations[0]
    else:
        refstation = opts.reference
    print 'Using %s for reference station'%(refstation)
    
    # Get ampltude and phase data for the reference station
    (refamp,refphase,reftimes,averaged) = getAmplPhase(p, refstation, ptype, corrs , polar, cords, opts.last)
    if averaged:
        print 'Warning: data is multi-channel. It is averaged to single channel for plotting purposes (this may also corrupt the statistics)'
        
    # This color code will be used in titles    
    colorCode = ''
    for corr in corrs:
        colorCode += CORR_COLOR_NAME[corr] + '=' + ptype + corr + ',' 
    colorCode = colorCode[:-1]
        
    if opts.quicklook:
        params = {'axes.labelsize': 6,
            'text.fontsize': 6,
            'legend.fontsize': 6,
            'xtick.labelsize': 6,
            'ytick.labelsize': 6,
            'figure.figsize': (8.27,11.69)}
        pylab.rcParams.update(params)
        ny = int(float(len(stations))/2.+0.5)
        
        pltamp = None
        pltphs = None
        deffigsize = (10.,6.)
        
        # Create figures and set the titles of the plots (they depend on the plotted correlations)
        if isValid(refamp):
            pltamp = pylab.figure(figsize=deffigsize)
            pylab.suptitle((title + ' Amplitude (' + colorCode + '), Ref = ' + refstation),fontsize=8)
        if isValid(refphase):
            pltphs = pylab.figure(figsize=deffigsize)
            pylab.suptitle((title + ' Phase (' + colorCode + '), Ref = ' + refstation),fontsize=8)
    
    c = 0 # subplot counter
    for station in stations:
        c += 1
        # Get the data for this station
        (valsamp,valsphase,times,averaged) = getAmplPhase(p, station, ptype, corrs , polar, cords, opts.last)
        
        fig = None
        
        if isValid(valsamp):
        
            # Create the sub-plot for the amplitude
            if opts.quicklook:
                pylab.figure(num=pltamp.number)
                pylab.subplot(ny,2,c)
            else:
                pylab.figure()
                pylab.subplot(211)
                pylab.suptitle('%s, Ref = %s'%(title,refstation))
            
            # Set the arguments for the pylab plot
            # And compute statistics if necessary
            plotArgs = []
            medians = {}
            for corr in corrs:
                valscorramp = valsamp[corr]
                plotArgs.extend([times,valscorramp,CORR_COLOR_CODE[corr]+'.'])
                if ofile != None:
                    numSamples = len(valscorramp)
                    mean = numpy.mean(valscorramp)
                    median = numpy.median(valscorramp)
                    medians[corr] = median
                    std = numpy.std(valscorramp)
                    ofile.write(OFILE_PATTERN % (station, (ptype+':'+corr), ('%d'%numSamples),('%.3f'% mean), ('%.3f'% median), ('%.3f'% std)))
            
            #Flush the stats file just in case the script is finished unexpectacly
            if ofile != None:
                ofile.flush()
            
            # Launch the pylab plotter por the amplitude for the selected correlations
            pylab.plot(*tuple(plotArgs))
            
            if opts.median:
                # Plot the median 
                for corr in corrs:
                    if len(medians):
                        median = medians[corr]
                    else:
                        median = numpy.median(valsamp[corr])
                    pylab.plot([times[0],times[-1]],[median,median],CORR_COLOR_CODE[corr]+'-')
            
            # Set labels 
            if opts.quicklook:
                pylab.ylabel(station,rotation='horizontal')
                if c<len(stations)-1: pylab.gca().set_xticklabels([])
                else: pylab.xlabel('Time [sec]')
            else:
                pylab.ylabel('Amplitude')
                pylab.gca().set_xticklabels([])
                pylab.title(station + colorCode)
            
            pylab.subplots_adjust(hspace=0.)
            
            # Set min and maximum values of amplitude plot (if provided)
            if ampmin is not None: pylab.ylim(ymin=ampmin)
            if ampmax is not None: pylab.ylim(ymax=ampmax)
            
        if isValid(valsphase):
            
            # Create the sub-plot for the phase (not for TEC, there is no phase)
            if opts.quicklook:
                pylab.figure(num=pltphs.number)
                pylab.subplot(ny,2,c)
            else:
                if fig==None:
                    # No figure was created in the AMPL
                    pylab.figure()
                    pylab.subplot(111)
                    pylab.suptitle('%s, Ref = %s'%(title,refstation))
                else:
                    pylab.subplot(212)
            
            # Set the arguments for the pylab plot
            plotArgs = []
            for corr in corrs:
                plotArgs.extend([times,((valsphase[corr]-refphase[corr]+numpy.pi)%(2.*numpy.pi)-numpy.pi),CORR_COLOR_CODE[corr]+'.'])
            pylab.plot(*tuple(plotArgs))
            
            # Set labels
            if opts.quicklook:
                pylab.ylabel(station,rotation='horizontal')
                if c<len(stations)-1: pylab.gca().set_xticklabels([])
                else: pylab.xlabel('Time [sec]')
            else:
                pylab.ylabel('Phase [rad]')
            pylab.xlabel('Time [sec]')
            pylab.ylim(ymin=-numpy.pi,ymax=numpy.pi)
            pylab.subplots_adjust(hspace=0.)
        
        # Show the plot (if not quicklook or displayall) or save it in a file        
        print 'Plotting %s'%(station)
        if opts.out == '' and not opts.quicklook and not opts.displayall:
            print 'Close figure to proceed'
            pylab.show()
        elif not opts.quicklook and opts.out != '':
            print '-> %s'%(opts.out+'-%s.%s'%(station,opts.ext))
            pylab.savefig(opts.out+'-%s.%s'%(station,opts.ext))
            
    # If we are in quicklook mode and no files have been generated we show all the plots or save tgem now
    if (opts.quicklook or opts.displayall) and opts.out == '':
        pylab.show()
    elif opts.quicklook:
        if pltamp != None:
            pylab.figure(num=pltamp.number)
            print '-> '+opts.out+'-amp.%s'%(opts.ext)
            pylab.savefig(opts.out+'-amp.%s'%(opts.ext))
        
        if pltphs != None:
            pylab.figure(num=pltphs.number)
            print '-> '+opts.out+'-phs.%s'%(opts.ext)
            pylab.savefig(opts.out+'-phs.%s'%(opts.ext))
    p = 0 # close the parmdb
    
    # Close the stats file if it is present
    if ofile != None:
        print '-> %s'%(opts.stats)
        ofile.close()

usage = 'Usage: %prog [options] parmdb'
description = "A program to plot solutions from a parmdb. For now the options are few (restricted to " +GAIN_TYPE+" or "+TEC_TYPE+") and this is really intended for a very specific type of parmdb (generated in MSSS processing) rather than a generic one that could contain anything. In particular DirectionalGain is currently NOT supported. For quality control plots, try the following two commands: solplot.py -q -a 0,0.06 -s RS\* -m -o gains-rs <MS>/instrument && solplot.py -q -a 0,0.06 -s CS\* -m -o gains-cs <MS>/instrument"
op = optparse.OptionParser(usage=usage, description=description)
op.add_option('-s','--stations',default='*',help='Filter for stations (e.g. CS*) [default *] Note that you have to escape asterisks, e.g. \"*\" or \\*',type='string')
op.add_option('-r','--reference',default='',help='Reference station [default \'\' means the first station in the list that matches the filter]',type='string')
op.add_option('-q','--quicklook',default=False,help='Quick look? [default False] This will plot all amplitudes on one page, and all phases on another',action='store_true')
op.add_option('-d','--displayall',default=False,help='Display all plots at once when not in quicklook mode? [default False] This is always True when in quicklook mode',action='store_true')
op.add_option('-o','--out',default='',help='Output plot basename [default \'\' means to show the plot instead of save it] Filenames will be <basename>-<station>.<extension> (or <basename>-<amp,phs>.<ext> in quicklook mode) if a basename is given',type='string')
op.add_option('-x','--ext',default='pdf',help='Image filename extension [default pdf]',type='string')
op.add_option('-a','--amplim',default=',',help='Comma-separated min,max for amplitude yaxis limits [default autorange]',type='string')
op.add_option('-p','--allpol',default=False,help='Plot all corr? [default False, which means only autocorrelations are plotted]',action='store_true')
op.add_option('-t','--title',default='',help='Title of the plots [default is to use the paramdb provided path]',type='string')
op.add_option('-l','--last',default=False,help='Include last value? [default False]',action='store_true')
op.add_option('-m','--median',default=False,help='Show median amplitude in the plots? [default False]',action='store_true')
op.add_option('-f','--stats',default='',help='Output amplitude statistics file name. If not specified no statistics file will be generated',type='string')
(opts, args) = op.parse_args()
if len(args) != 1:
    op.error('Incorrect number of arguments, please see help (-h)')
main(opts,args)

