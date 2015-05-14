def plugin_main(*args, **kwargs):
    print 'From example1:'
    print args,kwargs
    resdict = {}
    resdict['outtext'] = 'Output from ex1'
    return resdict
