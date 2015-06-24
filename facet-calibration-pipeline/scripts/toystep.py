def main(positionalargument, optional=0.0, threshpix=0):
    outdict = {}
    print 'FILENAME FOR TOYSTEP: ', str(positionalargument)
    # cast because we only pass strings in pipelines. stuff for later...
    derivedval = float(optional) / 2.0
    # you could as well do something else here depending on your inputdata
    # in positionalargument. that is different for every call of this script.

    # this one just gets a pass
    passthrough = threshpix

    # names in dict get saved as 'optionalhalf.mapfile' and 'threspix.mapfile'
    # naming scheme not fixed yet. and they get saved in the wrong directory...
    outdict['optionalhalf'] = derivedval
    outdict['threspix'] = passthrough
    
    return outdict
