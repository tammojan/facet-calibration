def plugin_main(*args, **kwargs):
    print 'From example2:'
    print args,kwargs
    resdict = {}
    if 'counter' in kwargs:
        print 'LOOPCOUNT: ', kwargs['counter']
        if kwargs['counter'] == 2:
            print 'REACHED BREAK'
            resdict['break'] = True
    return resdict
