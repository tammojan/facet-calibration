import os
from lofarpipe.support.data_map import DataMap
from lofarpipe.support.data_map import DataProduct


# mandatory arguments:
# cmdline for type of mapfile creation
# options: mapfile-dir, filename, identifier(name in parsetparset)
def plugin_main(args, **kwargs):
    #print 'PLUGIN KWARG: ', kwargs
    result = {}
    datamap = None
    if args[0] == 'mapfile_split_list':
        tmpmap = _create_mapfile_ato(kwargs['mapfile_in'])
        datamap = _split_listmap(tmpmap, int(kwargs['listsize']))
    if args[0] == 'mapfile_all_to_one':
        datamap = _create_mapfile_ato(kwargs['mapfile_in'])
    if args[0] == 'mapfile_list_ms':
        datamap = _create_mapfile_list(kwargs['folder'])
    if args[0] == 'mapfile_pythonlist_ms':
        datamap = _create_mapfile_pythonlist(kwargs['folder'])
    if args[0] == 'mapfile_from_folder':
        datamap = _create_mapfile_from_folder(kwargs['folder'])
    if args[0] == 'mapfile_from_parset':
        datamap = _create_mapfile_from_parset(kwargs['parset'], kwargs['identifier'])
    if args[0] == 'changed_mapfile_from_parset':
        datamap = _create_mapfile_from_parset(kwargs['parset'], kwargs['identifier'])
        folder = '/private/regression_test_runner_workdir/msss_imager_genericpipeline'
        for item in datamap:
            item.file = os.path.join(folder, 'concat.ms')
    if args[0] == 'mapfile_empty':
        fileid = os.path.join(kwargs['mapfile_dir'], kwargs['filename'])
        DataMap().save(fileid)
        result['mapfile'] = fileid
    if datamap:
        fileid = os.path.join(kwargs['mapfile_dir'], kwargs['filename'])
        datamap.save(fileid)
        result['mapfile'] = fileid
    return result


# helper function
def _create_mapfile_ato(inmap):
    maps = DataMap([])
    mapsin = DataMap.load(inmap)
    mapsin.iterator = DataMap.SkipIterator
    newlist = ''
    for i, item in enumerate(mapsin):
        newlist = newlist + item.file + ','
    newlist = newlist.rstrip(',')
    newlist = '[' + newlist + ']'
    maps.data.append(DataProduct('localhost', newlist, False))
    return maps

def _combine_local_map(inmap):
    map_out = DataMap([])
    map_in = DataMap.load(inmap)
    map_in.iterator = DataMap.SkipIterator
    local_files = {}
    for item in map_in:
        if item.host in local_files:
            local_files[item.host] += item.file + ','
        else:
            local_files[item.host] = item.file + ','
    for k, v in local_files.iteritems():
        v = v.rstrip(',')
        v = '[' + v + ']'
        map_out.data.append(DataProduct(k, v, False))
    return map_out

def _split_listmap(map_in, number):
    print 'MAP_IN: ', map_in
    map_out = DataMap([])
    for item in map_in:
        filelist = ((item.file.rstrip(']')).lstrip('[')).split(',')
        chunks = [filelist[i:i+number] for i in xrange(0, len(filelist), number)]
        print 'FILELIST: ', filelist
        print 'CHUNKS: ', chunks
        for slist in chunks:
            for i, name in enumerate(slist):
                #print 'NAMEB: ', name
                slist[i] = '"' + name + '"'
                #print 'NAMEA: ', name
            print 'SLIST: ', slist
            map_out.data.append(DataProduct(item.host, slist, False))
    return map_out


def _create_mapfile_from_folder(folder):
    maps = DataMap([])
    measurements = os.listdir(folder)
    measurements.sort()
    for ms in measurements:
        maps.data.append(DataProduct('localhost', folder + '/' + ms, False))
    return maps

def _create_mapfile_list(folder):
    maps = DataMap([])
    measurements = os.listdir(folder)
    measurements.sort()
    msfulll = []
    msfull = ''
    for ms in measurements:
        msfulll.append(os.path.join(folder , ms))
        #msfull.append(os.path.join(folder,ms).replace("'",""))
        msfull = msfull +os.path.join(folder , ms)+' '
    #msfull = msfull.rstrip(',')
    #msfull = '[' + msfull + ']'
    maps.data.append(DataProduct('localhost', msfull, False))
    #maps.file = msfulll
    return maps

def _create_mapfile_pythonlist(folder):
    maps = DataMap([])
    measurements = os.listdir(folder)
    measurements.sort()
    msfull = ''
    for ms in measurements:
        msfull = msfull + os.path.join(folder, ms)+','
    msfull = msfull.rstrip(',')
    msfull = '[' + msfull + ']'
    maps.data.append(DataProduct('localhost', msfull, False))
    return maps

def _create_mapfile_from_parset(parset, identifier):
    dps = parset.makeSubset(
        parset.fullModuleName('DataProducts') + '.'
    )
    datamap = DataMap([
        tuple(os.path.join(location, filename).split(':')) + (skip,)
        for location, filename, skip in zip(
            dps.getStringVector(identifier + '.locations'),
            dps.getStringVector(identifier + '.filenames'),
            dps.getBoolVector(identifier + '.skip'))
    ])
    return datamap
