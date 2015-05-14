import pyfits
from math import *
from astronomy import *
def read_fitsheader_keywords(filename,keywords):
    """
    Read the keywords from a fitheader
    """
    
    file = pyfits.open(filename)

    wordsarray = {}
    
    for word in keywords:
        wordsarray[word] = file[0].header[word]

    file.close()
    return wordsarray

def read_fitsimage(filename,decx,ray):
    """
    Read a fits image into a 2D array [RA][DEC]
    """
    file = pyfits.open(filename)

    data = decx*[0.0]
    for dec in range(0,decx):
        data[dec] = ray*[0.0]
        for ra in range(0,ray):
            data[dec][ra] = file[0].data[0][0][dec][ra]
            
    return data
    file.close()

def calc_beamarea(fitsfilename):
    # Given a fitsfile this calculates the beamarea in pixels
    
    f = pyfits.open(fitsfilename)
    
    beammaj = f[0].header['BMAJ'] # These are FWHM
    beammin = f[0].header['BMIN'] # These are FWHM
    beammaj = beammaj/(2.0*(2*log(2))**0.5) # Convert to sigma
    beammin = beammin/(2.0*(2*log(2))**0.5) # Convert to sigma
    pixarea = abs(f[0].header['CDELT1']*f[0].header['CDELT2'])

    beamarea = 2*pi*1.0*beammaj*beammin # Note that the volume of a two dimensional gaussian is simply 2*pi*A*sigma_x*sigma_y
    beamarea_pix = beamarea/pixarea

    return beamarea_pix

def find_pixel_coords(filename,decpix,rapix):
    # Given a filename and a declination and rightascension pixel number this will return the coordinates in J2000 degrees
    tmpf = pyfits.open(filename)
    naxis1 = tmpf[0].header['NAXIS1']
    naxis2 = tmpf[0].header['NAXIS2']
    racent = tmpf[0].header['CRVAL1']
    try:
        radelta = tmpf[0].header['CDELT1']
    except KeyError:
        radelta = tmpf[0].header['CD1_1']
    racentpix = tmpf[0].header['CRPIX1']
    deccent = tmpf[0].header['CRVAL2']
    try:
        decdelta = tmpf[0].header['CDELT2']
    except KeyError:
        decdelta = tmpf[0].header['CD2_2']
    deccentpix = tmpf[0].header['CRPIX2']

    rapixdif = rapix-racentpix
    decpixdif = decpix-deccentpix 

    radegdif = rapixdif*radelta
    decdegdif = decpixdif*decdelta

    decdegpix = deccent + decdegdif

    if rapixdif != 0:
        radegpix = sep_2pos(racent*deg2rad,decdegpix*deg2rad,radegdif*deg2rad,decdegpix*deg2rad)*rad2deg

    else:
        radegpix = racent

    tmpf.close()
    return decdegpix,radegpix


def find_pixel(filename,decdeg,radeg):
    # Given a filename and declination and right ascension this will return the corresponding pixels
    tmpf = pyfits.open(filename)
    naxis1 = tmpf[0].header['NAXIS1']
    naxis2 = tmpf[0].header['NAXIS2']
    racent = tmpf[0].header['CRVAL1']
    try:
        radelta = tmpf[0].header['CDELT1']
    except KeyError:
        radelta = tmpf[0].header['CD1_1']
    racentpix = tmpf[0].header['CRPIX1']
    deccent = tmpf[0].header['CRVAL2']
    try:
        decdelta = tmpf[0].header['CDELT2']
    except KeyError:
        decdelta = tmpf[0].header['CD2_2']
    deccentpix = tmpf[0].header['CRPIX2']        

    if radeg != racent:
        radiff = sepn(radeg*deg2rad,decdeg*deg2rad,racent*deg2rad,decdeg*deg2rad)*rad2deg
    else:
        radiff = 0.0
    if decdeg != deccent:
        decdiff = sepn(radeg*deg2rad,decdeg*deg2rad,radeg*deg2rad,deccent*deg2rad)*rad2deg
    else:
        decdiff =0.0
    
    if decdeg < deccent:
        decdiff = -decdiff
    if radeg < racent:
        radiff = -radiff
        
    offpixra = radiff/radelta
    offpixdec = decdiff/decdelta

    rapix = racentpix + offpixra
    decpix = deccentpix + offpixdec

    tmpf.close()

    return int(round(decpix,0)),int(round(rapix,0))
