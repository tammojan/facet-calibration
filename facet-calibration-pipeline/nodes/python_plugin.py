#                                                          LOFAR PIPELINE SCRIPT
#
#                                           running an executable with arguments
#                                                         Stefan Froehlich, 2014
#                                                      s.froehlich@fz-juelich.de
# ------------------------------------------------------------------------------

from __future__ import with_statement
from subprocess import CalledProcessError
import os
import shutil
import sys
import errno
import imp

from lofarpipe.support.pipelinelogging import CatchLog4CPlus
from lofarpipe.support.pipelinelogging import log_time
from lofarpipe.support.utilities import create_directory
from lofarpipe.support.utilities import catch_segfaults
from lofarpipe.support.lofarnode import LOFARnodeTCP
from lofarpipe.support.parset import Parset


class python_plugin(LOFARnodeTCP):
    """
    Basic script for running an executable with arguments.
    """

    def run(self, infile, executable, args, kwargs, work_dir='/tmp', parsetasfile=False, args_format='', environment=''):
        """
        This method contains all the needed functionality
        """
        # Debugging info
        self.logger.debug("infile            = %s" % infile)
        self.logger.debug("executable        = %s" % executable)
        self.logger.debug("working directory = %s" % work_dir)
        self.logger.debug("arguments         = %s" % args)
        self.logger.debug("arg dictionary    = %s" % kwargs)
        self.logger.debug("environment       = %s" % environment)

        self.environment.update(environment)

        # Time execution of this job
        with log_time(self.logger):
            #if os.path.exists(infile):
            self.logger.info("Processing %s" % infile)

            # Check if script is present
            if not os.path.isfile(executable):
                self.logger.error("Script %s not found" % executable)
                return 1

            # hurray! race condition when running with more than one process on one filesystem
            if not os.path.isdir(work_dir):
                try:
                    os.mkdir(work_dir, )
                except OSError as exc:  # Python >2.5
                    if exc.errno == errno.EEXIST and os.path.isdir(work_dir):
                        pass
                    else:
                        raise

            if parsetasfile:
                nodeparset = Parset()
                parsetname = os.path.join(work_dir, os.path.basename(infile) + '.parset')
                for k, v in kwargs.items():
                    nodeparset.add(k, v)
                nodeparset.writeFile(parsetname)
                args.insert(0, parsetname)

            try:
            # ****************************************************************
            # Run
                outdict = {}
                plugin = imp.load_source('main', executable)
                outdict = plugin.main(*args, **kwargs)

            except CalledProcessError, err:
                # CalledProcessError isn't properly propagated by IPython
                self.logger.error(str(err))
                return 1
            except Exception, err:
                self.logger.error(str(err))
                return 1

        if outdict:
            for k, v in outdict.items():
                self.outputs[k] = v
        # We need some signal to the master script that the script ran ok.
        self.outputs['ok'] = True
        return 0



if __name__ == "__main__":
    #   If invoked directly, parse command line arguments for logger information
    #                        and pass the rest to the run() method defined above
    # --------------------------------------------------------------------------
    jobid, jobhost, jobport = sys.argv[1:4]
    sys.exit(python_plugin(jobid, jobhost, jobport).run_with_stored_arguments())
