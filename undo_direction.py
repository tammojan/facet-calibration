#!/usr/bin/env python

import os
import pwd
import numpy
import sys
import logging
from facet_utilities import bg

username = pwd.getpwuid(os.getuid())[0]


def undo_direction(mslist):
    """
    Undo the last direction in a list of MSs.
    """
    b = bg(maxp=5)
    # copy of last is stored in corrected_data
    for ms in mslist:
        logging.info("undoing "+ms)
        b.run("taql 'update " + ms + " set SUBTRACTED_DATA_ALL=CORRECTED_DATA'")
    b.wait()


if __name__ == "__main__":
    # Logging configuration
    if os.path.exists("logging.conf"):
        logging.config.fileConfig('logging.conf')
        logger = logging.getLogger()
    else:
        # Start
        logger = logging.getLogger()
        logger.setLevel(logging.DEBUG)   
        formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
        # Log to STDOUT
        ch = logging.StreamHandler()
        ch.setFormatter(formatter)
        logger.addHandler(ch)
        # Log to file
        file_name = "dde.log"
        fh = logging.FileHandler(file_name) 
        fh.setFormatter(formatter)
        logger.addHandler(fh)
        logging.info('\n')

    # Read parameters
    if len(sys.argv)<2:
        raise Exception('Give the path to the setup code for the facet')

    logging.info('Using {} as the setup code'.format(sys.argv[1]))
    config = {}
    execfile(sys.argv[1])

    mslistorig = ["{name:s}_SB{b1:03d}-{b2:03d}.{res:s}.ms".format(name=NAME,res=RES,b1=b,b2=b+9) for b in BANDS]
    
    mslist = [ms for ms in mslistorig if os.path.exists(ms)]

    undo_direction(mslist)

