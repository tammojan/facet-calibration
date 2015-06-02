# Module for globally useful utilities

import os
from subprocess import Popen
import logging
import time

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
