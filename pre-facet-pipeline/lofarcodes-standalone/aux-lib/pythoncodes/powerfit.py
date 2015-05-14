from scipy.optimize import leastsq
import scipy
import scipy.stats
import pylab
import numpy as np
import random
import os,sys
sys.path.append('/Users/shi11n/Desktop/work/General/Total-Checked-Out-svn/my_library/trunk/pythoncodes/')
from astronomy import *
from constants import *
from fits import *
from maths import *
import re

# Program taken from http://wiki.scipy.org/Cookbook/FittingData#head-5eba0779a34c07f5a596bbcf99dbc7886eac18e5

class LogFormatterTeXExponent(pylab.LogFormatter, object):
    """Extends pylab.LogFormatter to use 
    tex notation for tick labels."""
    
    def __init__(self, *args, **kwargs):
        super(LogFormatterTeXExponent, 
              self).__init__(*args, **kwargs)
        
    def __call__(self, *args, **kwargs):
        """Wrap call to parent class with 
        change to tex notation."""
        label = super(LogFormatterTeXExponent, 
                      self).__call__(*args, **kwargs)
        label = re.sub(r'e(\S)0?(\d+)', 
                       r'\\times 10^{\1\2}', 
                       str(label))
        label = "$" + label + "$"
        return label



def model(t, coeffs):
    return coeffs[0] + t*coeffs[1]
def residuals(coeffs, y, t):
    return y - model(t, coeffs)

def powervals(freqrange,flux0,freq0,alpha):
    freqrange = np.array(freqrange)
    fluxvals = flux0*(freqrange/freq0)**alpha
    return fluxvals

def generate_testdata():
    sourceflux = 10.0
    sourcespec = 0.0
    freq0 = 1.4
    freqrange = np.arange(1.1,3.1,0.1)
    errorvals = np.zeros(len(freqrange))
    noerrorfluxvals = np.zeros(len(freqrange))
    channoise = 0.6
    specnoise = 0.0

    fittedpowerlaws  = []
    fittedfluxes = []

    fluxvals = np.zeros(len(freqrange))
    origfluxvals = np.zeros(len(freqrange))
    for j in range(0,len(freqrange)):
        noerrorfluxvals[j] = sourceflux*(freqrange[j]/freq0)**sourcespec
        fluxvals[j] = sourceflux*(freqrange[j]/freq0)**(random.gauss(sourcespec,specnoise))
        fluxvals[j] = random.gauss(fluxvals[j],channoise)
        errorvals[j] = abs(noerrorfluxvals[j]-fluxvals[j])
    return freqrange,fluxvals,errorvals,freq0

def fit_powerlaw(xdata,ydata,yerr,y0):
    removevals = []
    for i in range(0,len(xdata)):
        if xdata[i] < 0.0 or ydata[i] < 0.0:
            print 'Cant have negative values'
            removevals.append(i)
    xdata = np.delete(xdata,removevals)
    ydata = np.delete(ydata,removevals)
    yerr  = np.delete(yerr,removevals)
    logx = np.log10(xdata)
    logy = np.log10(ydata)
    #logyerr = yerr / ydata -- this is very similar to the logyerr below
    logyerr = (yerr)/(ydata*np.log(10)) # see e.g. http://en.wikipedia.org/wiki/Propagation_of_uncertainty for propagation of uncertainty
    # define our (line) fitting function
    fitfunc = lambda p, x: p[0] + p[1] * x
    errfunc = lambda p, x, y, err: (y - fitfunc(p, x)) / err

    pinit = [1.0, -1.0]
    out = scipy.optimize.leastsq(errfunc, pinit,
                    args=(logx, logy, logyerr), full_output=1)
    pfinal = out[0]
    covar = out[1]
    print pfinal
    print covar
    print logyerr

    indexErr = sqrt( covar[1][1] )
    ampErr = sqrt( covar[0][0] )
    
    index = pfinal[1]
    amps = 10**(index*np.log10(xdata)+pfinal[0])
    ampsmax = 10**((index+indexErr)*np.log10(xdata)+pfinal[0] + ampErr)
    ampsmin = 10**((index-indexErr)*np.log10(xdata)+pfinal[0] - ampErr)
    amp = 10**(index*np.log10(y0)+pfinal[0])
    ampmax = 10**((index+indexErr)*np.log10(y0)+pfinal[0] + ampErr)
    print ((index+indexErr)*np.log10(y0)+pfinal[0] + ampErr), ((index+indexErr)*np.log10(y0))
    ampmin = 10**((index-indexErr)*np.log10(y0)+pfinal[0] - ampErr)

    print 'Errors', amp, ampmax,ampmin,amp-ampmax,amp-ampmin
    ampErr = np.mean([abs(amp-ampmax),abs(amp-ampmin)])
    predicted = amp*(xdata/y0)**index
    chisquared = find_chi_squared_werrors(ydata,predicted,yerr)
    p_value = 1-scipy.stats.chi2.cdf(abs(chisquared),len(ydata))

    return amp,ampErr,index,indexErr,p_value

def plot_powerlaw(amp,ampErr,index,indexErr,xdata,ydata,yerr,y0,outputfilename):
    maxpower = index + indexErr
    minpower = index - indexErr
    meanpower = index
    maxflux = amp + ampErr
    minflux = amp - ampErr
    meanflux = amp
    print 'Power %s Power error %s'%(index,indexErr)
    print 'Flux %s Flux error %s'%(amp,ampErr)
    fig = pylab.figure()

    for j in range(0,len(xdata)):
        pylab.errorbar(xdata[j],ydata[j],yerr=yerr[j],xerr=None,color='r')
    pylab.plot(xdata,ydata,'ro')

    pylab.plot(xdata,powervals(xdata,meanflux,y0,meanpower),'k-')
    pylab.plot(xdata,powervals(xdata,maxflux,y0,maxpower),'k--')
    pylab.plot(xdata,powervals(xdata,minflux,y0,minpower),'k--')
    pylab.ylabel('Integrated flux (mJy)')
    pylab.xlabel('Frequency (GHz)')
    pylab.semilogx()
    pylab.semilogy()
    ax = fig.gca()
    pylab.xticks(np.arange(1.1,3.1,0.5))
    labels = np.arange(1.1,3.1,0.2)
    pylab.xticks(labels)
    ax.xaxis.set_minor_formatter(
        LogFormatterTeXExponent(base=10, 
                                labelOnlyBase=False))
    pylab.xlim(xmin=1.1,xmax=3.1)
    pylab.savefig(outputfilename)

    pylab.close()
    pylab.cla()

    return



#plotdata = True

#xdata,ydata,yerr,y0 = generate_testdata()
#amp,ampErr,index,indexErr = fit_powerlaw(xdata,ydata,yerr,y0)
#plot_powerlaw(amp,ampErr,index,indexErr,xdata,ydata,yerr,y0)



