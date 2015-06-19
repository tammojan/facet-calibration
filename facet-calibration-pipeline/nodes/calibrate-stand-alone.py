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
import subprocess

from lofarpipe.support.pipelinelogging import CatchLog4CPlus
from lofarpipe.support.pipelinelogging import log_time
from lofarpipe.support.utilities import catch_segfaults
from lofarpipe.support.lofarnode import LOFARnodeTCP
from lofarpipe.support.parset import Parset


class calibrate_stand_alone(LOFARnodeTCP):
    """
    Basic script for running bbs-reducer.
    Conversion of the calibrate stand alone shell script.
    """
    def __init__(self, jobid, jobhost, jobport):
        super(calibrate_stand_alone, self).__init__(job_id=jobid,host=jobhost,port=jobport)
        self.sourcedb_path = ''
        self.replace_sourcedb = False
        self.sourcedb = None
        self.sourcedb_basename = "sky"
        self.parmdb_path = ''
        self.replace_parmdb = False
        self.parmdb = None
        self.parmdb_basename = "instrument"
        self.observation = None
        self.parset = None
        self.catalog = ''
        self.dry_run = False
        self.work_dir = None
        self.infile = None
        self.executable = None

    def run(self, infile, executable, args, kwargs, work_dir='/tmp', parsetasfile=True, args_format='', environment=''):
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

        self.work_dir = work_dir
        self.infile = infile
        self.executable = executable

        if 'replace-sourcedb' in kwargs:
            self.replace_sourcedb = kwargs['replace-sourcedb']
            kwargs.pop('replace-sourcedb')
        if 'replace-parmdb' in kwargs:
            self.replace_parmdb = kwargs['replace-parmdb']
            kwargs.pop('replace-parmdb')
        if 'dry-run' in kwargs:
            self.dry_run = kwargs['dry-run']
            kwargs.pop('dry-run')
        if 'sourcedb' in kwargs:
            self.sourcedb = kwargs['sourcedb']
            kwargs.pop('sourcedb')
        if 'parmdb' in kwargs:
            self.parmdb = kwargs['parmdb']
            kwargs.pop('parmdb')
        if 'sourcedb-name' in kwargs:
            self.sourcedb_basename = kwargs['sourcedb-name']
            self.replace_sourcedb = True
            kwargs.pop('sourcedb-name')
        if 'parmdb-name' in kwargs:
            self.parmdb_basename = kwargs['parmdb-name']
            self.replace_parmdb = True
            kwargs.pop('parmdb-name')
        if 'force' in kwargs:
            self.replace_parmdb = True
            self.replace_sourcedb = True
            kwargs.pop('force')
        numthreads = 1
        if 'numthreads' in kwargs:
            numthreads = kwargs['numthreads']
            kwargs.pop('numthreads')
        args.append('--numthreads='+str(numthreads))
        if 'observation' in kwargs:
            self.observation = kwargs.pop('observation')
        if 'catalog' in kwargs:
            self.catalog = kwargs.pop('catalog')

        self.createsourcedb()
        self.createparmdb()
        if not 'no-columns' in kwargs:
            #if not kwargs['no-columns']:
            self.addcolumns()
        else:
            kwargs.pop('no-columns')

        args.append('--sourcedb=' + self.sourcedb_path)
        args.append('--parmdb=' + self.parmdb_path)

        args.append(self.observation)
        #catalog = None


        # Time execution of this job
        with log_time(self.logger):
            #if os.path.exists(infile):
            self.logger.info("Processing %s" % infile)

            # Check if script is present
            if not os.path.isfile(executable):
                self.logger.error("Executable %s not found" % executable)
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
                #args.insert(0, parsetname)
                args.append(parsetname)

            #if catalog is not None:
            #    args.append(catalog)

            try:
            # ****************************************************************
            #Run
                cmd = [executable] + args
                with CatchLog4CPlus(
                    work_dir,
                    self.logger.name + "." + os.path.basename(infile),
                    os.path.basename(executable),
                ) as logger:
                    # Catch segfaults and retry
                    catch_segfaults(
                        cmd, work_dir, self.environment, logger
                    )
            except CalledProcessError, err:
                # CalledProcessError isn't properly propagated by IPython
                self.logger.error(str(err))
                return 1
            except Exception, err:
                self.logger.error(str(err))
                return 1
        # We need some signal to the master script that the script ran ok.
        self.outputs['ok'] = True
        return 0

    def createsourcedb(self):
        self.logger.debug("create sourcedb")
        #global sourcedb_path
        self.sourcedb_path = os.path.join(self.observation, self.sourcedb_basename)
        if os.path.exists(self.sourcedb_path) and not self.replace_sourcedb:
             self.logger.debug('warning: sourcedb exists and will not be replaced (specify --replace-sourcedb to force replacement)')
        elif self.sourcedb:
            self.logger.debug("removing default sourcedb and copying new one in its place")
            if not self.dry_run:
                shutil.rmtree(self.sourcedb_path)
                shutil.copytree(self.sourcedb,self.sourcedb_path)
        else:
            self.logger.debug("makesourcedb in="+self.catalog+" out="+self.sourcedb_path+" format=<")
            if not self.dry_run:
                if os.path.exists(self.sourcedb_path):
                    shutil.rmtree(self.sourcedb_path)
                #subprocess.call(["makesourcedb","in="+self.catalog,"out="+sourcedb_path,"format=<"])
                args = ["in="+self.catalog,"out="+self.sourcedb_path,"format=<"]
                prog = os.path.join(os.path.dirname(self.executable), 'makesourcedb')
                self.execute(prog, args)

    def createparmdb(self):
        print "create parmdb"
        #global parmdb_path
        self.parmdb_path = os.path.join(self.observation, self.parmdb_basename)
        if os.path.exists(self.parmdb_path) and not self.replace_parmdb:
            self.logger.debug("warning: parmdb exists and will not be replaced (specify --replace-parmdb to force replacement)")
        elif self.parmdb:
            self.logger.debug("some message about parmdb")
            if not self.dry_run:
                shutil.rmtree(self.parmdb_path)
                shutil.copytree(self.parmdb,self.parmdb_path)
        else:
            self.logger.debug("creating parmdb")
            if not self.dry_run:
                if os.path.exists(self.parmdb_path):
                    shutil.rmtree(self.parmdb_path)

                prog = os.path.join(os.path.dirname(self.executable), 'parmdbm')
                parmdbdefaults="create tablename=\""+self.parmdb_path+""""
                  adddef Gain:0:0:Ampl  values=1.0
                  adddef Gain:1:1:Ampl  values=1.0
                  adddef Gain:0:0:Real  values=1.0
                  adddef Gain:1:1:Real  values=1.0
                  adddef DirectionalGain:0:0:Ampl  values=1.0
                  adddef DirectionalGain:1:1:Ampl  values=1.0
                  adddef DirectionalGain:0:0:Real  values=1.0
                  adddef DirectionalGain:1:1:Real  values=1.0
                  adddef Clock values=0.0, pert=1e-15
                  quit"""
                pipe = subprocess.Popen([prog], stdin=subprocess.PIPE)
                pipe.communicate(input=parmdbdefaults)
                pipe.stdin.close()
                if pipe.wait() != 0:
                    self.logger.error("Error while executing parmdbm")

    def addcolumns(self):
        self.logger.debug("add CASA imaging columns to MS")
        if not self.dry_run:
            #subprocess.call(['addImagingColumns.py',observation])
            args = [self.observation]
            prog = os.path.join(os.path.dirname(self.executable), 'addImagingColumns.py')
            self.execute(prog, args)

    def execute(self, executable, args):
        try:
        # ****************************************************************
        # Run
            cmd = [executable] + args
            with CatchLog4CPlus(
                self.work_dir,
                self.logger.name + "." + os.path.basename(self.infile),
                os.path.basename(self.executable),
            ) as logger:
                # Catch segfaults and retry
                catch_segfaults(
                    cmd, self.work_dir, self.environment, logger
                )
        except CalledProcessError, err:
            # CalledProcessError isn't properly propagated by IPython
            self.logger.error(str(err))
            return 1
        except Exception, err:
            self.logger.error(str(err))
            return 1

if __name__ == "__main__":
    #   If invoked directly, parse command line arguments for logger information
    #                        and pass the rest to the run() method defined above
    # --------------------------------------------------------------------------
    jobid, jobhost, jobport = sys.argv[1:4]
    sys.exit(calibrate_stand_alone(jobid, jobhost, jobport).run_with_stored_arguments())
