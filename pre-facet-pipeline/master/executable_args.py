#                                                          LOFAR PIPELINE SCRIPT
#
#                                           running an executable with arguments
#                                                         Stefan Froehlich, 2014
#                                                      s.froehlich@fz-juelich.de
# ------------------------------------------------------------------------------

import copy
import sys
import os

import lofarpipe.support.lofaringredient as ingredient

from lofarpipe.support.baserecipe import BaseRecipe
from lofarpipe.support.remotecommand import RemoteCommandRecipeMixIn
from lofarpipe.support.remotecommand import ComputeJob
from lofarpipe.support.data_map import DataMap, validate_data_maps, align_data_maps
from lofarpipe.support.parset import Parset


class executable_args(BaseRecipe, RemoteCommandRecipeMixIn):
    """
    Basic script for running an executable with arguments.
    Passing a mapfile along so the executable can process MS.
    """
    inputs = {
        'executable': ingredient.ExecField(
            '--executable',
            help="The full path to the relevant executable",
            optional=True
        ),
        'arguments': ingredient.ListField(
            '-a', '--arguments',
            help="List of arguments for the executable. Will be added as ./exe arg0 arg1...",
            default='',
            optional=True
        ),
        'nodescript': ingredient.StringField(
            '--nodescript',
            help="Name of the node script to execute",
            default='executable_args',
            optional=True
        ),
        'parset': ingredient.FileField(
            '-p', '--parset',
            help="Path to the arguments for this executable. Will be converted to --key=value",
            optional=True
        ),
        'inputkey': ingredient.StringField(
            '-i', '--inputkey',
            help="Parset key that the executable will recognize as key for inputfile",
            default='',
            optional=True
        ),
        'outputkey': ingredient.StringField(
            '-0', '--outputkey',
            help="Parset key that the executable will recognize as key for outputfile",
            default='',
            optional=True
        ),
        'inputkeys': ingredient.ListField(
            '--inputkeys',
            help="List of parset keys that the executable will recognize as key for inputfile",
            default=[],
            optional=True
        ),
        'outputkeys': ingredient.ListField(
            '--outputkeys',
            help="List of parset keys that the executable will recognize as key for outputfile",
            default=[],
            optional=True
        ),
        'mapfiles_in': ingredient.ListField(
            '--mapfiles-in',
            help="List of the input mapfiles containing the names of the "
                 "data to run the recipe on",
            default=[],
            optional=True
        ),
        'mapfiles_out': ingredient.ListField(
            '--mapfiles-out',
            help="List of the output mapfiles containing the names of the "
                 "data produced by the recipe",
            default=[],
            optional=True
        ),
        'mapfile_in': ingredient.StringField(
            '--mapfile-in',
            help="Name of the input mapfile containing the names of the "
                 "MS-files to run the recipe",
            default='',
            optional=True
        ),
        'mapfile_out': ingredient.StringField(
            '--mapfile-out',
            help="Name of the output mapfile containing the names of the "
                 "MS-files produced by the recipe",
            default='',
            optional=True
        ),
        'skip_infile': ingredient.BoolField(
            '--skip-infile',
            help="Dont give the input file to the executable.",
            default=False,
            optional=True
        ),
        'skip_outfile': ingredient.BoolField(
            '--skip-outfile',
            help="Dont produce an output file",
            default=False,
            optional=True
        ),
        'inplace': ingredient.BoolField(
            '--inplace',
            help="Manipulate input files inplace",
            default=False,
            optional=True
        ),
        'outputsuffixes': ingredient.ListField(
            '--outputsuffixes',
            help="Suffixes for the outputfiles",
            default=[]
        ),
        'parsetasfile': ingredient.BoolField(
            '--parsetasfile',
            help="Will the argument be a parsetfile or --opt=var",
            default=False
        ),
        'args_format': ingredient.StringField(
            '--args_format',
            help="Will change the format of the arguments",
            default='gnu'
        ),
        'max_per_node': ingredient.IntField(
            '--max_per_node',
            help="Sets the number of jobs per node",
            default=0
        )
    }

    outputs = {
        'mapfile': ingredient.FileField(
            help="The full path to a mapfile describing the processed data"
        )
    }

    def go(self):
        if 'executable' in self.inputs:
            executable = self.inputs['executable']

        self.logger.info("Starting %s run" % executable)
        super(executable_args, self).go()

        # *********************************************************************
        # try loading input/output data file, validate output vs the input location if
        #    output locations are provided
        try:
            indatas = []
            if self.inputs['mapfile_in']:
                indata = DataMap.load(self.inputs['mapfile_in'])
                indatas.append(indata)

            if self.inputs['mapfiles_in']:
                for item in self.inputs['mapfiles_in']:
                    indatas.append(DataMap.load(item))
                self.inputs['mapfile_in'] = self.inputs['mapfiles_in'][0]
            #else:
            #    indatas.append(indata)
        except Exception:
            self.logger.error('Could not load input Mapfile %s' % self.inputs['mapfile_in'])
            return 1
        if self.inputs['mapfile_out']:
            try:
                outdata = DataMap.load(self.inputs['mapfile_out'])
            except Exception:
                self.logger.error('Could not load output Mapfile %s' % self.inputs['mapfile_out'])
                return 1
            # sync skip fields in the mapfiles
            align_data_maps(indatas[0], outdata)
        else:
            # ouput will be directed in the working directory if no output mapfile is specified
            outdata = copy.deepcopy(indatas[0])
            if not self.inputs['inplace']:
                for item in outdata:
                    item.file = os.path.join(
                        self.inputs['working_directory'],
                        self.inputs['job_name'],
                        os.path.basename(item.file) + '.' + os.path.split(str(executable))[1]
                    )
                self.inputs['mapfile_out'] = os.path.join(os.path.dirname(self.inputs['mapfile_in']), os.path.basename(executable) + '.' + 'mapfile')
            else:
                self.inputs['mapfile_out'] = self.inputs['mapfile_in']

        if not validate_data_maps(indatas[0], outdata):
            self.logger.error(
                "Validation of data mapfiles failed!"
            )
            return 1

        # Handle multiple outputfiles
        outputsuffix = self.inputs['outputsuffixes']
        outputmapfiles = []
        prefix = os.path.join(self.inputs['working_directory'], self.inputs['job_name'])
        for name in outputsuffix:
            outputmapfiles.append(copy.deepcopy(indatas[0]))
            for item in outputmapfiles[-1]:
                item.file = os.path.join(
                    prefix,
                    os.path.basename(item.file) + '.' + os.path.split(str(executable))[1] + '.' + name
                )

        # prepare arguments
        arglist = self.inputs['arguments']
        parsetdict = {}
        if 'parset' in self.inputs:
            parset = Parset()
            parset.adoptFile(self.inputs['parset'])
            for k in parset.keys:
                parsetdict[k] = str(parset[k])

        #for k in parset.keys:
        #    arglist.append('--' + k + '=' + parset.getString(k))
        #if not self.inputs['inputkey'] and not self.inputs['skip_infile']:
        #    arglist.insert(0, None)

        # construct multiple input data
        inputlist = []
        if not self.inputs['inputkeys'] and self.inputs['inputkey']:
            self.inputs['inputkeys'].append(self.inputs['inputkey'])

        if indatas:
            for item in indatas:
                item.iterator = DataMap.SkipIterator
            for mfile in indatas:
                inputlist.append([])
                for inp in mfile:
                    inputlist[-1].append(inp.file)

        # ********************************************************************
        # Call the node side of the recipe
        # Create and schedule the compute jobs
        command = "python %s" % (self.__file__.replace('master', 'nodes')).replace('executable_args', self.inputs['nodescript'])
        indatas[0].iterator = outdata.iterator = DataMap.SkipIterator
        jobs = []
        for i, (outp, inp,) in enumerate(zip(
            outdata, indatas[0])
        ):
            arglist_copy = copy.deepcopy(arglist)
            parsetdict_copy = copy.deepcopy(parsetdict)

            if self.inputs['inputkeys'] and not self.inputs['skip_infile']:
                for name, value in zip(self.inputs['inputkeys'], inputlist):
                    if arglist_copy and name in arglist_copy:
                        ind = arglist_copy.index(name)
                        arglist_copy[ind] = value[i]
                    else:
                        parsetdict_copy[name] = value[i]

            if self.inputs['outputkey'] and not self.inputs['skip_infile']:
                if arglist_copy and self.inputs['outputkey'] in arglist_copy:
                    ind = arglist_copy.index(self.inputs['outputkey'])
                    arglist_copy[ind] = outp.file
                else:
                    parsetdict_copy[self.inputs['outputkey']] = outp.file

            jobs.append(
                ComputeJob(
                    inp.host, command,
                    arguments=[
                        inp.file,
                        executable,
                        arglist_copy,
                        parsetdict_copy,
                        prefix,
                        self.inputs['parsetasfile'],
                        self.inputs['args_format'],
                        #self.inputs['working_directory'],
                        self.environment
                    ]
                )
            )
        max_per_node = self.inputs['max_per_node']
        self._schedule_jobs(jobs, max_per_node)
        for job, outp in zip(jobs, outdata):
            if job.results['returncode'] != 0:
                outp.skip = True
            #print 'JOBRESULTS: ', job.results
        # *********************************************************************
        # Check job results, and create output data map file
        if self.error.isSet():
            # Abort if all jobs failed
            if all(job.results['returncode'] != 0 for job in jobs):
                self.logger.error("All jobs failed. Bailing out!")
                return 1
            else:
                self.logger.warn(
                    "Some jobs failed, continuing with succeeded runs"
                )
        self.logger.debug("Writing data map file: %s" % self.inputs['mapfile_out'])
        #outdata.save(self.inputs['mapfile_out'])
        #self.outputs['mapfile'] = self.inputs['mapfile_out']
        mapdict = {}
        for item, name in zip(outputmapfiles, outputsuffix):
            item.save(os.path.join(prefix, name + '.' + 'mapfile'))
            mapdict[name] = os.path.join(prefix, name + '.' + 'mapfile')
            #self.outputs[name] = name + '.' + 'mapfile'
        if not outputsuffix:
            outdata.save(self.inputs['mapfile_out'])
            self.outputs['mapfile'] = self.inputs['mapfile_out']
        else:
            self.outputs.update(mapdict)
            self.outputs['mapfile'] = os.path.join(prefix, outputsuffix[0] + '.' + 'mapfile')
        return 0

if __name__ == '__main__':
    sys.exit(executable_args().main())
