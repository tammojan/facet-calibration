# Module for globally useful utilities

import os
from subprocess import Popen
import logging
import time
from math import floor

####
# Run parallel commands

def run(c,proceed=False,quiet=False):
    '''
    Run command c, throwing an exception if the return value is
    non-zero (which means that the command failed). By default
    (proceed==False) catches bad return values and raises an
    exception.
    '''

    if not(quiet):
        logging.debug('Running: '+c)
    retval=os.system(c)
    if retval!=0:
        report='FAILED to run '+c+' -- return value was '+str(retval)
        logging.error(report)
        if not(proceed):
            raise Exception(report)
    return retval

class bg:
    '''
    Simple wrapper class around Popen that lets you run a number of
    processes in the background and then wait for them to end
    successfully. The list of processes is maintained internally so
    you don't have to manage it yourself: just start by creating a
    class instance. By default (proceed==False) catches bad return
    values, kills all other background processes and raises an exception.

    Optionally you can set a value maxp on creating the instance. If
    you do this then run will block if you attempt to have more than
    this number of processes running concurrently.
    '''
    def __init__(self,quiet=False,pollint=1,proceed=False,maxp=None):
        self.pl=[]
        self.quiet=quiet
        self.pollint=pollint
        self.proceed=proceed
        self.maxp=maxp

    def run(self,c):
        if self.maxp:
            if len(self.pl)>=self.maxp:
                # too many processes running already. Wait till one finishes
                self.wait(queuelen=self.maxp-1)
        p=Popen('exec '+c,shell=True)
        if not(self.quiet):
            logging.debug('Process '+str(p.pid)+' started: '+c)
        self.pl.append(p)
        return p.pid
    def wait(self,queuelen=0):
        pl=self.pl
        while len(pl)>queuelen:
            for p in pl:
                retval=p.poll()
                if retval is not None:
                    pl.remove(p)
                    if retval==0:
                        if not(self.quiet):
                            logging.debug('Process '+str(p.pid)+' ended OK')
                    else:
                        report='Process '+str(p.pid)+' ended with return value '+str(retval)
                        logging.error(report)
                        if not(self.proceed):
                            for p2 in pl:
                                p2.kill()
                            raise Exception(report)
                        else:
                            print 'WARNING:',report
            time.sleep(self.pollint)
        self.pl=pl 
        return

####
# Coordinates

def angsep(ra1deg, dec1deg, ra2deg, dec2deg):
    """
    Returns angular separation between two coordinates (all in degrees)
    Input:
      * ra1deg - RA of the first position
      * dec1deg - dec of the first position
      * ra2deg - RA of the second position
      * dec2deg - dec of the second position
    Output:
      * Angular separation in degrees.
    """

    ra1rad = ra1deg*numpy.pi/180.0
    dec1rad = dec1deg*numpy.pi/180.0
    ra2rad = ra2deg*numpy.pi/180.0
    dec2rad = dec2deg*numpy.pi/180.0

    # calculate scalar product for determination
    # of angular separation
    x = numpy.cos(ra1rad)*numpy.cos(dec1rad)*numpy.cos(ra2rad)*numpy.cos(dec2rad)
    y = numpy.sin(ra1rad)*numpy.cos(dec1rad)*numpy.sin(ra2rad)*numpy.cos(dec2rad)
    z = numpy.sin(dec1rad)*numpy.sin(dec2rad)

    if x+y+z >= 1:
        rad = 0
    else:
        rad=numpy.acos(x+y+z)

    # Angular separation
    deg = rad*180/numpy.pi
    return deg


def ra_to_str(dra, ndec=2,delim=':'):
    '''
    converts a single decimal degrees ra to hh:mm:ss.s
    '''
    if delim == 'h':
        delim1 = 'h'
        delim2 = 'm'
    else:
        delim1 = delim
        delim2 = delim

    dra = dra/15.
    dd = floor(dra)
    dfrac = dra - dd
    dmins = dfrac*60.
    dm = floor(dmins)
    dsec = (dmins-dm)*60.
    if round(dsec, ndec) == 60.00:
        dsec = 0.
        dm += 1
    if dm == 60.:
        dm = 0.
        dd += 1
    sra = '%02d%s%02d%s%05.2f' %(dd,delim1,dm,delim2,dsec)
    return sra


def dec_to_str(ddec,ndec=1,delim=':'):
    '''
    converts a single decimal degrees dec to dd:mm:ss.s
    '''
    if delim == 'd':
        delim1 = 'd'
        delim2 = 'm'
    else:
        delim1 = delim
        delim2 = delim

    dd = floor(ddec)
    dfrac = ddec - dd
    dmins = dfrac*60.
    dm = floor(dmins)
    dsec = (dmins-dm)*60.
    if round(dsec, ndec) == 60.0:
        dsec = 0.
        dm += 1
    if dm == 60.:
        dm = 0.
        dd += 1
    sdec = '%02d%s%02d%s%04.1f' %(dd,delim1,dm,delim2,dsec)
    return sdec
