import os
from lofarpipe.support.data_map import DataMap


def plugin_main(args, **kwargs):
    result = {}
    folder = '/home/sfroehli/scratch'
    datamap = kwargs[map_to_copy]
    for item in datamap:
        item.file = os.path.join(folder, 'concat.ms')
    fileid = kwargs['folder'] + kwargs['filename']
    datamap.save(fileid)
    result['mapfile'] = fileid
    return result