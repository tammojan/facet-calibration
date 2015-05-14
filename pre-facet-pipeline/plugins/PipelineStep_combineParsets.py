from lofarpipe.support.parset import Parset


def plugin_main(*args, **kwargs):
    parset = Parset(kwargs['first_parset'])
    parset.adoptFile(kwargs['second_parset'])
    parset.writeFile(kwargs['result_parset'] + '_feedback_file')

