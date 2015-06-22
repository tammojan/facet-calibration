import os
from lofarpipe.support.data_map import DataMap
from lofarpipe.support.data_map import DataProduct


# mandatory arguments:
# cmdline for type of mapfile creation
# options: mapfile-dir, filename, identifier(name in parsetparset)
def plugin_main(args, **kwargs):
    print 'PLUGIN KWARG: ', kwargs
    result = {}
    datamap = None
    fileid = kwargs['mapfile_in']
    datamap = DataMap.load(fileid)
    #if kwargs['change_files']:
    #    for item in datamap:
    #        item.file = kwargs['change_files']
    if kwargs['join_files']:
        for item in datamap:
            item.file = os.path.join(item.file,kwargs['join_files'])
    if kwargs['newname']:
        fileid = os.path.join(os.path.dirname(fileid), kwargs['newname'])
    if datamap:
        print 'Wrinting mapfile: ',fileid
        datamap.save(fileid)
        result['mapfile'] = fileid
    return result
