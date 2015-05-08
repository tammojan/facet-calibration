#################################################################################
#                                                                               #
# Written by Leah Morabito, 3 May 2013                                          #
#                                                                               #
# This script is intended to get rid of time dependence in gain solutions       #
# so they can be applied to another field by calculating the                    #
# new Real Amplitude and setting the new Imaginary Amplitude to zero            #
# (Effectively getting rid of the phases for an amplitude-only gain transfer)   #
#                                                                               #
#       NOTE: THIS ZEROES OUT PHASES AND GETS RID OF TIME DEPENDENCE            #
#                                                                               #
#################################################################################

# modified by W. Williams, 15 Jan 2014
# take median value instead of average and print out some values
# and make a copy of the parmdb

import lofar.parmdb as lp
import numpy
import sys
import os

if len(sys.argv) < 2:
    print 'Please give a parmdb name.'
    print 'Usage: python fixparmdb_zero_median.py <parmdbname>'

filename = sys.argv[1]
outfilename = sys.argv[2]

if os.path.isdir(outfilename):
    print 'output file exists'
    sys.exit()
os.system('cp -r %s %s' %(filename, outfilename))
filename = outfilename

parmdbmtable = lp.parmdb(filename)

dictionary = parmdbmtable.getValuesGrid('*')

real_names = parmdbmtable.getNames('*Real*')
imaginary_names = parmdbmtable.getNames('*Imag*')

names = parmdbmtable.getNames()

for name in real_names:
    imaginary_name = name.replace('Real','Imag')
    real_values = dictionary[name]['values']
    imaginary_values = dictionary[imaginary_name]['values']
    new_real_values = numpy.sqrt(real_values**2. + imaginary_values**2.)
    #new_imaginary_values = numpy.zeros(dictionary[imaginary_name]['values'].shape)
    new_real_value = numpy.median(new_real_values)
    #new_imaginary_value = numpy.median(new_imaginary_values)
    new_imaginary_value = 0.
    parmdbmtable.addDefValues(name,new_real_value)
    parmdbmtable.addDefValues(imaginary_name,new_imaginary_value)
    print '%s %.5f %.5f' %(name, new_real_value, new_imaginary_value)

parmdbmtable.deleteValues('*')
parmdbmtable.flush()

# The output values can be accessed through getDefValues
