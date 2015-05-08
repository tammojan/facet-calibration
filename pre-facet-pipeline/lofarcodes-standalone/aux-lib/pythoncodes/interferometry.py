from math import *
import numpy as np
import scipy.special

def time_smearing(resolution,delta_Theta):
    #http://www.cv.nrao.edu/course/astr534/Interferometers1.html
    # output and input units are equal
    # Condition that delta_Theta * delta_T << 1.37E4 * Resolution
    # where delta_Theta is the angular radius of the image and delta_T is the time smoothing.

    # Given resolution and field of view (delta_theta) this returns a condition on the time smearing at the correlator
    delta_T = 1.37E4*resolution/delta_Theta

    print 'Time averaging should be less than %s'%delta_T


    # Time average smearing amplitude loss in http://adsabs.harvard.edu/full/1989ASPC....6..247B
    Reduction = 1-1.22E-9*(delta_Theta/resolution)**2.0 * delta_T**2.0

    print 'At radius %s a source will only have %s percent of its flux if data smoothed in time to %s'%(delta_Theta,Reduction,delta_T)
    
    return delta_T
    
def bandwidth_smearing(resolution,freq,delta_Theta):
    #http://www.cv.nrao.edu/course/astr534/Interferometers1.html
    # Output and input units are equal 
    # Condition that delta_Theta * delta_freq << Freq * resolution
    # where delta_Theta is the offset from the pointing centre and delta_freq is the bandwidth smoothing.

    # Given resolution, freq and offset this gives the condition for the delta_freq

    delta_freq = freq*resolution/delta_Theta

    print 'Bandwidth averaging should be much less than %s'%delta_freq

    # Bandwidth smearing amplitude loss in http://adsabs.harvard.edu/full/1989ASPC....6..247B
    
    beta = (delta_freq/freq) * (delta_Theta/resolution)
    gamma = 2*(np.log(2)**0.5)
    Reduction = ((np.pi**0.5)/(gamma * beta)) * (scipy.special.erf(beta*gamma/2.0))

    print 'At radius %s a source will only have %s percent of its flux if data smoothed in freq to %s'%(delta_Theta,Reduction,delta_freq)

    return delta_freq
