from scipy.integrate import quad
import pylab
from math import *
import numpy as np
import os,sys
sys.path.append('/Users/shi11n/Desktop/Total-Checked-Out-svn/my_library/trunk/pythoncodes/')
from maths import *
from astronomy import *
from lists import *
from constants import *
from inputs import *
from pyslalib import slalib, sladoc

# This program is to demonstrate how to integrate a function with python
# Here I present 3 examples, the source counts at 1.4GHz, 15GHz and a simple cos function


def intergrand1(x):
    return cos(x)

def source_counts_1400MHz(S):
    # Takes S as an np array or single value, input values in Jy. Outputs dn_ds according to Hopkins 2002.
    dn_ds = (S)**(-2.5)*10**(0.859*(np.log10(S/1E-3))**0.0 + 0.508*(np.log10(S/1E-3))**1.0 + 0.376*(np.log10(S/1E-3))**2.0 - 0.049*(np.log10(S/1E-3))**3.0 - 0.121*(np.log10(S/1E-3))**4.0 + 0.057*(np.log10(S/1E-3))**5.0 - 0.008*(np.log10(S/1E-3))**6.0)

    return dn_ds

def source_counts_1400MHz_plot(S):
    # Takes S as an np array or single value, input values in Jy. Outputs (dn_ds)/S^(-2.5) according to Hopkins 2002.
    dn_ds_S_25 = 10**(0.859*(np.log10(S/1E-3))**0.0 + 0.508*(np.log10(S/1E-3))**1.0 + 0.376*(np.log10(S/1E-3))**2.0 - 0.049*(np.log10(S/1E-3))**3.0 - 0.121*(np.log10(S/1E-3))**4.0 + 0.057*(np.log10(S/1E-3))**5.0 - 0.008*(np.log10(S/1E-3))**6.0)

    return dn_ds_S_25


def source_counts_15GHz(S):
    # Takes S as an np array or single value, inputs in Jy. Outputs dn_ds according to AMI Consortium Davis 1010.
    if S < 2.2E-3:
        dn_ds = 340.0*(S)**(-1.81)
    else:
        dn_ds = 48.0*(S)**(-2.13)
    return dn_ds


lowerlim = 0.0
upperlim = 1.0
result,error = quad(intergrand1,lowerlim,upperlim)

print 'Integrate cos(x) between %s and %s ' %(lowerlim,upperlim)
print 'Integral = %s, Error = %s ' %(result,error)


print 'Plotting an example of a function that could be integrated'

numbers = np.arange(0.01,1000.0,0.01) # X axis in mJy
pylab.plot(numbers,source_counts_1400MHz_plot(numbers/1000.0)) # Note that source_counts_1400MHz_plot takes input in Jy
pylab.semilogy()
pylab.semilogx()
pylab.ylim(ymin=1.0,ymax=1000.0)
pylab.xlim(xmin=0.05,xmax=1000.0)
pylab.xlabel('S(mJy)')
pylab.ylabel('dN/dS (/S^{-2/5})(Jy^(1.5)sr^(-1))')
pylab.title('Source counts at 1.4GHz according to Hopkins 2002')
pylab.show()


lowerlim = 4*4.0E-6 # in Jy
upperlim = 1.0 #in Jy
result,error = quad(source_counts_1400MHz,lowerlim,upperlim)

print 'Number of 1.4GHz sources from %sJy to %sJy per steradian is %s with error %s ' %(lowerlim,upperlim,result,error)
# Area of 6 arcmin
area = (((6.0/60.0)*(6.0/60.0))/2.0)*(1.0/steradians2degsquared)
print 'Number of 1.4GHz sources within area %s steradians is %s ' %(area,result*area)

lowerlim = 0.04E-3 #in Jy
upperlim = 1.0 #in Jy
result,error = quad(source_counts_1400MHz,lowerlim,upperlim)

print 'Number of 1.4GHz sources from %sJy to %sJy per steradian is %s with error %s ' %(lowerlim,upperlim,result,error)
# Area of 6 arcmin
area = (((6.0/60.0)*(6.0/60.0))/2.0)*(1.0/steradians2degsquared)
print 'Number of 1.4GHz sources within area %s steradians is %s ' %(area,result*area)


lowerlim = 9.0E-3 #in Jy
upperlim = 10.0E-3 #in Jy
result,error = quad(source_counts_15GHz,lowerlim,upperlim)
print 'Number of 15.0GHz sources from %sJy to %sJy per steradian is %s with error %s ' %(lowerlim,upperlim,result,error)
area = 47.83*(1.0/steradians2degsquared)
print 'Number of 15.0GHz sources within area %s steradians is %s ' %(area,result*area)

