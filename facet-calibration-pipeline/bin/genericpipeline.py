import os
import sys
import copy
from lofarpipe.support.parset import Parset
from lofarpipe.support.control import control

from lofarpipe.support.loggingdecorators import duration
from lofarpipe.support.data_map import DataMap, DataProduct, validate_data_maps
from lofarpipe.support.lofarexceptions import PipelineException
from lofarpipe.support.utilities import create_directory
import logging
from lofarpipe.support.pipelinelogging import getSearchingLogger
import lofarpipe.support.lofaringredient as ingredient
import loader
import lofarpipe.support.utilities as utilities
from ConfigParser import NoOptionError, NoSectionError
from ConfigParser import SafeConfigParser as ConfigParser
overwrite = False

class GenericPipeline(control):

    inputs = {
        'loglevel': ingredient.StringField(
            '--loglevel',
            help="loglevel",
            default='INFO',
            optional=True
        )
    }

    def __init__(self):
        control.__init__(self)
        self.parset = Parset()
        self.input_data = {}
        self.output_data = {}
        self.parset_feedback_file = None
        #self.logger = None#logging.RootLogger('DEBUG')
        self.name = ''

        #if not overwrite:
        #    self.inputs['job_name'] = 'generic-pipeline'
        # if not self.inputs.has_key("start_time"):
        #     import datetime
        #     self.inputs["start_time"] = datetime.datetime.utcnow().replace(microsecond=0).isoformat()
        # if not hasattr(self, "config"):
        #     self.config = self._read_config()
        # #self._read_config()
        # # ...and task files, if applicable
        # if not self.inputs.has_key("task_files"):
        #     try:
        #         self.inputs["task_files"] = utilities.string_to_list(
        #             self.config.get('DEFAULT', "task_files")
        #         )
        #     except NoOptionError:
        #         self.inputs["task_files"] = []
        # self.task_definitions = ConfigParser(self.config.defaults())
        # print >> sys.stderr, "Reading task definition file(s): %s" % \
        #                      ",".join(self.inputs["task_files"])
        # self.task_definitions.read(self.inputs["task_files"])
        #    self.go()

    def usage(self):
        """
        Display usage
        """
        print >> sys.stderr, "Usage: %s [options] <parset-file>" % sys.argv[0]
        print >> sys.stderr, "Parset structure should look like:\n" \
                             "NYI"
        #return 1

    def go(self):
        #"""
        #Read the parset-file that was given as input argument, and set the
        #jobname before calling the base-class's `go()` method.
        #"""
        try:
            parset_file = os.path.abspath(self.inputs['args'][0])
        except IndexError:
            #return self.usage()
            self.usage()

        # Set job-name to basename of parset-file w/o extension, if it's not
        # set on the command-line with '-j' or '--job-name'
        if not 'job_name' in self.inputs:
            self.inputs['job_name'] = (
                os.path.splitext(os.path.basename(parset_file))[0])
            self.name = self.inputs['job_name']
        try:
            self.logger
        except:
            self.logger = getSearchingLogger(self.name)
            self.logger.setLevel(self.inputs['loglevel'])
        # Call the base-class's `go()` method.
        return super(GenericPipeline, self).go()

#    def pipeline_logic(self):
#        print 'Dummy because of stupid wrapping inside the framework'
#        if overwrite:
#            self.execute_pipeline()

    #def execute_pipeline(self):
    def pipeline_logic(self):
        try:
            parset_file = os.path.abspath(self.inputs['args'][0])
        except IndexError:
            return self.usage()
        try:
            self.parset.adoptFile(parset_file)
            self.parset_feedback_file = parset_file + "_feedback"
        except RuntimeError:
            print >> sys.stderr, "Error: Parset file not found!"
            return self.usage()
        self._replace_values()
        # just a reminder that this has to be implemented
        validator = GenericPipelineParsetValidation(self.parset)
        if not validator.validate_pipeline():
            self.usage()
            exit(1)
        if not validator.validate_steps():
            self.usage()
            exit(1)

        #set up directories
        job_dir = self.config.get("layout", "job_directory")
        parset_dir = os.path.join(job_dir, "parsets")
        mapfile_dir = os.path.join(job_dir, "mapfiles")
        # Create directories for temporary parset- and map files
        create_directory(parset_dir)
        create_directory(mapfile_dir)

        # *********************************************************************
        # maybe we dont need a subset but just a steplist
        # at the moment only a list with stepnames is given for the pipeline.steps parameter
        # pipeline.steps=[vdsmaker,vdsreader,setupparmdb1,setupsourcedb1,ndppp1,....]
        # the names will be the prefix for parset subsets
        pipeline_args = self.parset.makeSubset(
            self.parset.fullModuleName('pipeline') + '.')

        # *********************************************************************
        # forward declaration of things. just for better overview and understanding whats in here.
        # some of this might be removed in upcoming iterations, or stuff gets added.
        step_name_list = pipeline_args.getStringVector('steps')
        step_control_dict = {}
        step_parset_files = {}
        step_parset_obj = {}
        activeloop = ['']
        # construct the list of step names and controls
        self._construct_steps(step_name_list, step_control_dict, step_parset_files, step_parset_obj, parset_dir)
        # initial parameters to be saved in resultsdict so that recipes have access to this step0
        # double init values. 'input' should be considered deprecated
        # self.name would be consistent to use in subpipelines
        resultdicts = {'input': {
            'parset': parset_file,
            'parsetobj': self.parset,
            'job_dir': job_dir,
            'parset_dir': parset_dir,
            'mapfile_dir': mapfile_dir}}

        resultdicts.update({self.name: {
            'parset': parset_file,
            'parsetobj': self.parset,
            'job_dir': job_dir,
            'parset_dir': parset_dir,
            'mapfile_dir': mapfile_dir}})

        if 'pipeline.mapfile' in self.parset.keys:
            resultdicts['input']['mapfile'] = str(self.parset['pipeline.mapfile'])
            resultdicts[self.name]['mapfile'] = str(self.parset['pipeline.mapfile'])

        # *********************************************************************
        # main loop
        # there is a distinction between recipes and plugins for user scripts.
        # plugins are not used at the moment and might better be replaced with master recipes
        while step_name_list:
            stepname = step_name_list.pop(0)
            step = step_control_dict[stepname]
            #step_parset = step_parset_obj[stepname]
            inputdict = {}
            inputargs = []
            resultdict = {}
            # default kind_of_step to recipe.
            try:
                kind_of_step = step.getString('kind')
            except:
                kind_of_step = 'recipe'
            try:
                typeval = step.getString('type')
            except:
                typeval = ''
            #self._construct_cmdline(inputargs, step, resultdicts)

            additional_input = {}

            if stepname in step_parset_obj:
                additional_input = self._construct_step_parset(step_parset_obj[stepname],
                                                               resultdicts,
                                                               step_parset_files[stepname],
                                                               stepname)

            # stepname not a valid input for old recipes
            if kind_of_step == 'recipe':
                if self.task_definitions.get(typeval, 'recipe') == 'executable_args':
                    inputdict = {'stepname': stepname}
                    inputdict.update(additional_input)

            self._construct_cmdline(inputargs, step, resultdicts)

            if stepname in step_parset_files:
                inputdict['parset'] = step_parset_files[stepname]


            self._construct_input(inputdict, step, resultdicts)

            # hack, popping 'type' is necessary, why? because you deleted kind already in parsets
            try:
                inputdict.pop('type')
            except:
                pass
            try:
                inputdict.pop('kind')
            except:
                pass
            # \hack
            # more hacks. Frameworks DictField not properly implemented. Construct your own dict from input.
            # python buildin functions cant handle the string returned from parset class.
            if 'environment' in inputdict.keys():
                val = inputdict['environment'].rstrip('}').lstrip('{').replace(' ', '')
                splitval = str(val).split(',')
                valdict = {}
                for item in splitval:
                    valdict[item.split(':')[0]] = item.split(':')[1]
                inputdict['environment'] = valdict

            # subpipeline. goal is to specify a pipeline within a pipeline.
            # load other existing pipeline parset and add them to your own.
            if kind_of_step == 'pipeline':
                subpipeline_parset = Parset()
                subpipeline_parset.adoptFile(typeval)
                submapfile = ''
                subpipeline_steplist = subpipeline_parset.getStringVector('pipeline.steps')

                if 'pipeline.mapfile' in subpipeline_parset.keys:
                    submapfile = subpipeline_parset['pipeline.mapfile']
                    subpipeline_parset.remove('pipeline.mapfile')
                if 'mapfile' in inputdict.keys():
                    submapfile = inputdict.pop('mapfile')
                resultdicts.update({os.path.splitext(os.path.basename(typeval))[0]: {
                    'parset': typeval,
                    'mapfile': submapfile,
                }})
                #todo: take care of pluginpathes and everything other then individual steps
                # make a pipeline parse methods that returns everything needed.
                # maybe as dicts to combine them to one

                subpipeline_parset.remove('pipeline.steps')
                if 'pipeline.pluginpath' in subpipeline_parset.keys:
                    subpipeline_parset.remove('pipeline.pluginpath')
                checklist = copy.deepcopy(subpipeline_steplist)
                for k in subpipeline_parset.keys:
                    if 'loopsteps' in k:
                        for item in subpipeline_parset.getStringVector(k):
                            checklist.append(item)
                # *********************************************************************
                # master parset did not handle formatting and comments in the parset.
                # proper format only after use of parset.makesubset. then it is a different object
                # from a different super class :(. this also explains use of parset.keys and parset.keys()
                # take the parset from subpipeline and add it to the master parset.
                # *********************************************************************
                # replace names of steps with the subpipeline stepname to create a unique identifier.
                # replacement values starting with ! will be taken from the master parset and overwrite
                # the ones in the subpipeline. only works if the ! value is already in the subpipeline
                for k in subpipeline_parset.keys:
                    if not str(k).startswith('#'):
                        val = subpipeline_parset[k]
                        if not str(k).startswith('!'):
                            for item in checklist:
                                if item in str(val):
                                    val = str(val).replace(item, stepname + '-' + item)

                            self.parset.add(stepname + '-' + k, str(val))
                        else:
                            self.parset.add(k, str(val))
                for i, item in enumerate(subpipeline_steplist):
                    subpipeline_steplist[i] = stepname + '-' + item
                for item in step_parset_obj[stepname].keys():
                    for k in self.parset.keys:
                        if str(k).startswith('!') and item in k:
                            self.parset.remove(k)
                            self.parset.add('! ' + item, str(step_parset_obj[stepname][item]))

                self._replace_values()

                self._construct_steps(subpipeline_steplist, step_control_dict, step_parset_files, step_parset_obj, parset_dir)
                for j in reversed(subpipeline_steplist):
                    name = j
                    step_control_dict[name] = step_control_dict[j]
                    step_name_list.insert(0, name)

            # loop
            if kind_of_step == 'loop':
                # remember what loop is running to stop it from a conditional step
                if activeloop[0] is not stepname:
                    activeloop.insert(0, stepname)
                # prepare
                counter = 0
                breakloop = False
                if stepname in resultdicts:
                    counter = int(resultdicts[stepname]['counter']) + 1
                    breakloop = resultdicts[stepname]['break']
                loopsteps = step.getStringVector('loopsteps')

                # break at max iteration or when other step sets break variable
                if counter is step.getInt('loopcount'):
                    breakloop = True
                if not breakloop:
                    # add loop steps to the pipeline including the loop itself
                    step_name_list.insert(0, stepname)
                    self._construct_steps(loopsteps, step_control_dict, step_parset_files, step_parset_obj, parset_dir)
                    for j in reversed(loopsteps):
                        name = j
                        step_control_dict[name] = step_control_dict[j]
                        step_name_list.insert(0, name)
                    # results for other steps to check and write states
                    resultdict = {'counter': counter, 'break': breakloop}
                else:
                    # reset values for second use of the loop (but why would you do that?)
                    resultdict = {'counter': -1, 'break': False}
                    activeloop.pop(0)

            # recipes
            if kind_of_step == 'recipe':
                with duration(self, stepname):
                    resultdict = self.run_task(
                        typeval,
                        inputargs,
                        **inputdict
                    )

            # plugins
            if kind_of_step == 'plugin':
                with duration(self, stepname):
                    resultdict = loader.call_plugin(typeval, pipeline_args.getString('pluginpath'),
                                                    inputargs,
                                                    **inputdict)
            resultdicts[stepname] = resultdict

            # breaking the loopstep
            # if the step has the keyword for loopbreaks assign the value
            if resultdict is not None and 'break' in resultdict:
                if resultdict['break']:
                    resultdicts[activeloop[0]]['break'] = resultdict['break']

    # *********************************************************************
    # build the inputs for the master recipes.
    def _construct_input(self, inoutdict, controlparset, resdicts):
        # intermediate backward compatibility for opts subparset
        if controlparset.fullModuleName('opts'):
            argsparset = controlparset.makeSubset(controlparset.fullModuleName('opts') + '.')
        # hack
        elif 'loopcount' not in controlparset.keys():
            argsparset = controlparset
        else:
            argsparset = controlparset.makeSubset(controlparset.fullModuleName('imaginary') + '.')
        # \hack

        self._replace_output_keyword(inoutdict, argsparset, resdicts)

    def _construct_cmdline(self, inoutargs, controlparset, resdicts):
        argsparset = controlparset.makeSubset(controlparset.fullModuleName('cmdline') + '.')
        for k in argsparset.keys():
            if argsparset.getString(k).__contains__('.output.'):
                step, outvar = argsparset.getString(k).split('.output.')
                inoutargs.append(resdicts[step][outvar])
            else:
                inoutargs.append(argsparset.getString(k))
        try:
            controlparset.remove('cmdline.inmap')
        except:
            pass

    def _construct_steps(self, step_name_list, step_control_dict, step_parset_files, step_parset_obj, parset_dir):
        step_list_copy = (copy.deepcopy(step_name_list))
        counter = 0
        while step_list_copy:
            counter -= 1
            stepname = step_list_copy.pop(-1)
            fullparset = self.parset.makeSubset(self.parset.fullModuleName(str(stepname)) + '.')
            subparset = fullparset.makeSubset(fullparset.fullModuleName('control') + '.')
            number = 0
            for item in step_list_copy:
                if item == stepname:
                    number += 1
            if number != 0:
                stepname += str(number)
            step_name_list[counter] = stepname
            step_control_dict[stepname] = subparset
            # double implementation for intermediate backward compatibility
            if fullparset.fullModuleName('parsetarg') or fullparset.fullModuleName('argument'):
                if fullparset.fullModuleName('parsetarg'):
                    stepparset = fullparset.makeSubset(fullparset.fullModuleName('parsetarg') + '.')
                if fullparset.fullModuleName('argument'):
                    stepparset = fullparset.makeSubset(fullparset.fullModuleName('argument') + '.')
                # *********************************************************************
                # save parsets
                # either a filename is given in the main parset
                # or files will be created from subsets with stepnames.parset as filenames
                # for name, parset in step_parset_dict.iteritems():
                try:
                    file_parset = Parset(stepparset.getString('parset'))
                    for k in file_parset.keys:
                        if not k in stepparset.keys():
                            stepparset.add(k, str(file_parset[k]))
                    stepparset.remove('parset')
                except:
                    pass
                # parset from task.cfg
                try:
                    file_parset = Parset(self.task_definitions.get(str(subparset['type']), 'parset'))
                    for k in file_parset.keys:
                        if not k in stepparset.keys():
                            stepparset.add(k, str(file_parset[k]))
                except:
                    pass
                # for parset in control section
                try:
                    file_parset = Parset(subparset.getString('parset'))
                    for k in file_parset.keys:
                        if not k in stepparset.keys():
                            stepparset.add(k, str(file_parset[k]))
                    subparset.remove('parset')
                except:
                    pass
                step_parset = os.path.join(parset_dir, stepname + '.parset')
                stepparset.writeFile(step_parset)
                step_parset_files[stepname] = step_parset
                step_parset_obj[stepname] = stepparset

    def _replace_output_keyword(self, inoutdict, argsparset, resdicts):
        for k in argsparset.keys():
            keystring = argsparset.getString(k)
            if keystring.__contains__('.output.'):
                if keystring.__contains__(','):
                    keystring = keystring.rstrip(']')
                    keystring = keystring.lstrip('[')
                    vec = []
                    for item in keystring.split(','):
                        if item.__contains__('.output.'):
                            step, outvar = item.split('.output.')
                            vec.append(resdicts[step][outvar])
                        else:
                            vec.append(item)
                    inoutdict[k] = vec
                else:
                    step, outvar = argsparset.getString(k).split('.output.')
                    inoutdict[k] = resdicts[step][outvar]
            else:
                inoutdict[k] = argsparset.getString(k)

    def _construct_step_parset(self, argsparset, resdicts, filename, stepname):
        addvals = {'inputkeys': [], 'mapfiles_in': [], 'arguments': []}
        # hack for original order of args
        tmp_keys = argsparset.keys()
        ordered_keys = []
        for orig in self.parset.keys:
            for item in tmp_keys:
                if (stepname + '.') in orig and ('argument.'+item in orig and not 'argument.'+item+'.' in orig):
                    ordered_keys.append(item)
                    continue
        # \hack
        for k in ordered_keys:
            keystring = argsparset.getString(k)
            if keystring.__contains__('.output.'):
                if keystring.__contains__(','):
                    keystring = keystring.rstrip(']')
                    keystring = keystring.lstrip('[')
                    vec = []
                    for item in keystring.split(','):
                        if item.__contains__('.output.'):
                            step, outvar = item.split('.output.')
                            vec.append(resdicts[step][outvar])
                            addvals['inputkeys'].append(resdicts[step][outvar])
                            addvals['mapfiles_in'].append(resdicts[step][outvar])
                        else:
                            vec.append((item))
                    argsparset.replace(k,str(vec))
                    if k == 'flags':
                        addvals['arguments'] = (vec)
                        argsparset.remove(k)
                else:
                    step, outvar = argsparset.getString(k).split('.output.')
                    argsparset.replace(k,str(resdicts[step][outvar]))
                    if isinstance(resdicts[step][outvar], str):
                        addvals['inputkeys'].append(resdicts[step][outvar])
                        addvals['mapfiles_in'].append(resdicts[step][outvar])
                    if k == 'flags':
                        addvals['arguments'] = str(argsparset[k])
                        argsparset.remove(k)
            else:
                if k == 'flags':
                    addvals['arguments'] = str(argsparset[k])
                    argsparset.remove(k)
        argsparset.writeFile(filename)
        return addvals

    def _get_parset_dicts(self):
        return {}

    def show_tasks(self):
        tasklist = []
        tasklist = self.task_definitions.sections()
        for item in tasklist:
            print item
        #return tasklist

    def show_task(self, task):
        task_parset = Parset()
        if self.task_definitions.has_option(task,'parset'):
            task_parset.adoptFile(self.task_definitions.get(task,'parset'))
            print 'possible arguments: key    =    value'
            for k in task_parset.keys:
                print '                   ',k,'    ','=','    ',task_parset[k]

    def _add_step(self):
        steplist = []

    def _replace_values(self):
        replacedict = {}
        try:
            import imp
            plugin = imp.load_source('main', str(self.parset['prepare']))
            replacedict = plugin.main()
        except:
            pass
        for check in self.parset.keys:
            if str(check).startswith('!'):
                replacedict[str(check).lstrip('!').lstrip(' ')] = str(self.parset[check])
        #print 'REPLACEDICT: ',replacedict
        for check in self.parset.keys:
            if not str(check).startswith('#'):
                for k, v in replacedict.iteritems():
                    if '{{ '+k+' }}' in str(self.parset[check]):
                        replacestring = str(self.parset[check]).replace('{{ '+k+' }}',v)
                        self.parset.replace(check,replacestring)


class GenericPipelineParsetValidation():

    def __init__(self, parset):
        self.parset = parset
        #self.validate_pipeline()
        #self.validate_steps()

    def validate_pipeline(self):
        try:
            self.parset.getStringVector('pipeline.steps')
            return True
        except:
            print "Error: No pipeline steps defined"
            return None

    def validate_steps(self):
        try:
            print 'NYI: validate_steps'
            return True
        except:
            print "Error: Steps validation failed"
            return None


if __name__ == '__main__':
    overwrite = True
    sys.exit(GenericPipeline().main())