from math import *
from constants import *
from pyslalib import slalib, sladoc
import numpy as np

#-----------------------------------------------------------


def sepn(r1,d1,r2,d2):

    """
    Calculate the separation between 2 sources, RA and Dec must be
    given in radians. Returns the separation in radians [TWS]
    """

    # NB slalib sla_dsep does this
    # www.starlink.rl.ac.uk/star/docs/sun67.htx/node72.html
    
    cos_sepn=np.sin(d1)*np.sin(d2) + np.cos(d1)*np.cos(d2)*np.cos(r1-r2)
    sepn = np.arccos(cos_sepn)
    
    return sepn    

#-----------------------------------------------------------

def RA_J2000_to_radians(hours,minutes,seconds):

    """
    Convert right asscension from J2000 to radians
    """
    #(rarad,status)  = slalib.sla_dtf2r(hours,minutes,seconds)
    #print rarad
    hours = float(hours)
    minutes= float(minutes)
    seconds = float(seconds)
    rarad = (hours*(360.0/24.0) + minutes*(360.0/(24.0*60.0)) + seconds*(360.0/(24.0*60.0*60.0)))*deg2rad

    return rarad

#-----------------------------------------------------------

def DEC_J2000_to_radians(degrees,minutes,seconds):

    """
    Convert declination from J2000 to radians

    """
    # Slalib library is not great with negative declination
    degrees = float(degrees)
    minutes = float(minutes)
    seconds = float(seconds)
    if degrees < 0.0:
        degrees = -degrees
        switchsign = True
    else:
        switchsign = False
    #(decrad,status) = slalib.sla_daf2r(degrees,minutes,seconds)
    #print decrad
    decrad = (degrees + minutes/60.0 + seconds/(60.0*60.0))*deg2rad

    if switchsign:
        decrad = -decrad

    return decrad

#-----------------------------------------------------------


def radians_to_RA_J2000(radians):

    """
    Convert from radians to J2000 for RA
    """

    RA_J2000 = slalib.sla_cr2tf(2,radians)[1]

    return RA_J2000


  #-----------------------------------------------------------


def radians_to_DEC_J2000(radians):

    """
    Convert from radians to J2000 for DEC
    """

    DEC_J2000 = slalib.sla_cr2af(2,radians)[1]

    if radians < 0.0:
        DEC_J2000[0] = DEC_J2000[0]

    return DEC_J2000

#-----------------------------------------------------------

def sep_2pos(r1,d1,del_r,d2):

    """
    Given a seperation between two sources in RA and the RA of one source the RA
 of the source source is given. All in radians.
    """

    cos_r1_minus_r2 = (cos(del_r) - (sin(d1)*sin(d2)))/(cos(d1)*cos(d2))

    if del_r < 0.0:
        r2 = r1 - (acos(cos_r1_minus_r2))
    else:
        r2 =  (acos(cos_r1_minus_r2)) + r1
    
    return r2

#-----------------------------------------------------------

def calc_frequency(wavelength):

    """
    Convert a wavelength into a frequency
    """
    freq = c/wavelength

    return freq

#-----------------------------------------------------------

def calc_wavelength(frequency):
    """
    Convert a frequency into a wavelength.
    """

    wavelength = c/frequency

    return wavelength

#-----------------------------------------------------------

def calc_resolution(frequency,diameter):

    """
    Given the frequency and the diameter this returns the resolution in radians
    """

    wavelength = calc_wavelength(frequency)

    resolution = wavelength/diameter

    print "Resolution {0:3f} radians".format(resolution)
    print "Resolution {0:3f} degrees".format(resolution*rad2deg)
    print "Resolution {0:3f} arcmin".format(resolution*rad2arcmin)
    print "Resolution {0:3f} arcsec".format(resolution*rad2arcsec)

    return resolution

#-----------------------------------------------------------

def calc_maxrecoverable(frequency,diameter_min):

    """
    Given the frequency (Hz) and the shortest baseline (m) of an interferometer this returns the maximum recoverable scale

    Equation taken from section 7.4 of the alma technical handbook -- available on https://almascience.nrao.edu/documents-and-tools
    """

    wavelength = calc_wavelength(frequency)
    
    max_coverable = 0.6*wavelength/(diameter_min)

    print "Maximum recoverable scale {0:3f} radians".format(max_coverable)
    print "Maximum recoverable scale {0:3f} degrees".format(max_coverable*rad2deg)
    print "Maximum recoverable scale {0:3f} arcmin".format(max_coverable*rad2arcmin)
    print "Maximum recoverable scale {0:3f} arcsec".format(max_coverable*rad2arcsec)

    return max_coverable

#-----------------------------------------------------

def calc_solid_angle(fwhm):

    """
    Given the FWHM in arcsec this calculates the solid angle in steradians
    """

    fwhm = fwhm*arcsec2rad

    solid_angle = (pi*fwhm**2.0)/(4*log(2))

    return solid_angle



#--------------------------

def brightness_sensitivity(frequency,beamfwhm,flux_sensitivity):

    """ 
    Given the frequency in Hz, synthesized beamf whm in arcsec and flux xssensitivity in mJy/beam this calculates the
    brightness sensitivity in mK
    """

    wavelength = calc_wavelength(frequency)

    solid_angle = calc_solid_angle(beamfwhm)

    flux_sensitivity = flux_sensitivity*jy2si
    
    brightness_sensitivity =(flux_sensitivity/solid_angle)*((wavelength**2)/(2*k_b))

    return brightness_sensitivity
    
#----------------------------

def removenoise(noisearray):

    noisearray = noisearray.flatten()

    noiseconverge = 0
    orignoise = np.std(noisearray)
    
    while noiseconverge == 0:
        thresholdupp = np.mean(noisearray) + 3*np.std(noisearray)
        thresholdlow = np.mean(noisearray) - 3*np.std(noisearray)
        noisepix = np.array(filter(lambda x: x<thresholdupp,noisearray))
        noisepix = np.array(filter(lambda x: x>thresholdlow,noisearray))
        rms = np.std(noisearray)
        if rms >= orignoise*0.99:
            noiseconverge = 1
        else:
            orignoise = rms

    return noisepix

#------------------------
def faraday_depth_max(lowestfreq,deltafreq):
    # Calculates the maximum observable faraday depth (see http://arxiv.org/pdf/1201.3161v2.pdf)

    deltalambda_sq = (3E8/(lowestfreq))**2.0 - (3E8/(deltafreq+lowestfreq))**2.0
    maxdepth = 3.0**0.5/(deltalambda_sq)
    return maxdepth
#------------------------
def faraday_largest_scale(lowestfreq):
    # Calculates the largest detectable scale in faraday depth (see http://arxiv.org/pdf/1201.3161v2.pdf)

    maxscale = pi/((3E8/lowestfreq)**2.0)
    return maxscale


#------------------------
def positional_accuracy(frequency,diameter,snr):
    # Given a signal to noise ratio, a frequency and diamter this returns the expected positional accuracy
    # See https://www.iram.fr/IRAMFR/IS/IS2002/html_2/node130.html or Reid et al. 1998 for formula
    resolution = calc_resolution(frequency,diameter)
    pos_accuracy = 0.5*resolution/snr #Radians

    return pos_accuracy
