import os
from lofarpipe.support.data_map import DataMap
from lofarpipe.support.data_map import DataProduct
import argparse
from argparse import RawTextHelpFormatter


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
    if args[0] == 'mapfile_from_folder':
        datamap = _create_mapfile_from_folder(kwargs['folder'])
    if args[0] == 'mapfile_from_parset':
        datamap = _create_mapfile_from_parset(kwargs['parset'], kwargs['identifier'])
    if args[0] == 'add_suffix_to_file':
        datamap == _add_name(kwargs['mapfile_in'], kwargs['add_suffix_to_file'])
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
    return MultiDataMap(DataMap.load(inmap))


def _split_listmap(dmap, number):
    dmap.split_list(number)
    return dmap


def _create_mapfile_from_folder(path):
    dmap = MapfileManager(folder=path)
    return dmap


def _create_mapfile_from_parset(parset, identifier):
    dmap = MapfileManager()
    dmap.set_data_from_parset(parset, identifier)
    return dmap

def _add_name(inmap, suffix):
    dmap = DataMap.load(inmap)
    for item in dmap:
        item.file += suffix
    return dmap


class MapfileManager(DataMap):

    def __init__(self, folder=None, pattern=None):
        super(MapfileManager, self).__init__()
        if folder:
            self.from_folder(folder, pattern)
        #self.map = DataMap([])

    def expand(self, number, hostlist=None, filelist=None):
        if hostlist:
            if len(hostlist) != number:
                print 'Error: length of hostlist should correspond to number of expansions'
                exit(1)
        else:
            print 'Info: no hostlist given. Will use "localhost" instead'
            hostlist = []
            for item in range(number):
                hostlist.append('localhost')

        if filelist:
            if len(filelist) != number:
                print 'Error: length of hostlist should correspond to number of expansions'
                exit(1)
        else:
            print 'Info: no filelist given. Will use "dummy" instead'
            filelist = []
            for item in range(number):
                filelist.append('dummy')

        prodlist = []
        for h, f in zip(hostlist, filelist):
            prodlist.append(DataProduct(h, f, False))

        self._set_data(prodlist)

    def insert(self, place, data):
        self._insert(place, data)

    def _insert(self, place, data, dtype=DataProduct):
        try:
            if isinstance(data, dtype):
                self._data.insert(place, data)
            elif isinstance(data, dict):
                self._data.insert(place, dtype(**data))
            elif all(isinstance(item, tuple) for item in data):
                self._data.insert(place, dtype(*data))
            else:
                raise TypeError
        except TypeError:
            raise DataMapError("Failed to validate data map: %s" % repr(data))

    def append(self, data):
        self._append(data)

    def _append(self, data, dtype=DataProduct):
        try:
            if isinstance(data, dtype):
                self._data.append(data)
            elif isinstance(data, dict):
                self._data.append(dtype(**data))
            elif all(isinstance(item, tuple) for item in data):
                self._data.append(dtype(*data))
            else:
                raise TypeError
        except TypeError:
            raise DataMapError("Failed to validate data map: %s" % repr(data))

    def delete(self, host=None, data=None, skip=None, pattern=None):
        for i, item in enumerate(self._data):
            if item.host == host or item.file == data or item.skip == skip:
                del self._data[i]
            if pattern:
                if pattern in item.file:
                    del self._data[i]

    def from_folder(self, folder, pattern=None, exclude_pattern=False):
        measurements = os.listdir(folder)
        measurements.sort()
        for ms in measurements:
            if pattern:
                if pattern in ms and not exclude_pattern:
                    self._append(DataProduct('localhost', folder + '/' + ms, False))
                elif pattern in ms and exclude_pattern:
                    pass
                elif pattern not in ms and exclude_pattern:
                    self._append(DataProduct('localhost', folder + '/' + ms, False))
            else:
                self._append(DataProduct('localhost', folder + '/' + ms, False))

    def set_data_from_parset(self, parset, identifier):
        dps = parset.makeSubset(
            parset.fullModuleName('DataProducts') + '.'
        )
        self._set_data([
            tuple(os.path.join(location, filename).split(':')) + (skip,)
            for location, filename, skip in zip(
                dps.getStringVector(identifier + '.locations'),
                dps.getStringVector(identifier + '.filenames'),
                dps.getBoolVector(identifier + '.skip'))
        ])

    def from_parts(self, host='localhost', data='dummy', skip=False, ntimes=1):
        hostlist = self._input_to_list(host)
        datalist = self._input_to_list(data)
        skiplist = self._input_to_list(skip)
        if len(hostlist) is not len(datalist) or len(hostlist) is not len(skiplist) or len(hostlist) is not ntimes:
            print 'Length of parts is not equal. Will expand to max length given.'
            maxval = max(len(hostlist), len(datalist), len(skiplist), ntimes)
            lastval = hostlist[-1]
            if len(hostlist) is not maxval:
                for x in range(len(hostlist), maxval):
                    hostlist.append(lastval)

            lastval = datalist[-1]
            if len(datalist) is not maxval:
                for x in range(len(datalist), maxval):
                    datalist.append(lastval)

            lastval = skiplist[-1]
            if len(skiplist) is not maxval:
                for x in range(len(skiplist), maxval):
                    skiplist.append(lastval)
        prodlist = []
        for h, f, z in zip(hostlist, datalist, skiplist):
            prodlist.append(DataProduct(h, f, z))
        self._set_data(prodlist)

    def _input_to_list(self, data):
        datalist = []
        if isinstance(data, list):
            datalist = data
        else:
            datalist.append(data)
        return datalist


class MultiDataProduct(DataProduct):
    """
    Class representing multiple files in a DataProduct.
    """
    def __init__(self, host=None, file=None, skip=True):
        super(MultiDataProduct, self).__init__(host, file, skip)
        if not file:
            self.file = list()
        else:
            self._set_file(file)
        print 'FILE: ', self.file

    def __repr__(self):
        """Represent an instance as a Python dict"""
        return (
            "{{'host': '{0}', 'file': {1}, 'skip': {2}}}".format(self.host, self.file, str(self.skip))
        )

    def __str__(self):
        """Print an instance as 'host:[filelist]'"""
        return ':'.join((self.host, str(self.file)))

    def _set_file(self, data):
        try:
            # Try parsing as a list
            if isinstance(data, list):
                self.file = data
            if isinstance(data, DataProduct):
                self._from_dataproduct(data)
            if isinstance(data, DataMap):
                self._from_datamap(data)

        except TypeError:
            raise DataProduct("No known method to set a filelist from %s" % str(file))

    def _from_dataproduct(self, prod):
        print 'setting filelist from DataProduct'
        self.host = prod.host
        self.file = prod.file
        self.skip = prod.skip

    def _from_datamap(self, inmap):
        print 'setting filelist from DataMap'
        filelist = {}
        for item in inmap:
            if not item.host in filelist:
                filelist[item.host] = []
            filelist[item.host].append(item.file)
        self.file = filelist['i am']

    def append(self, item):
        self.file.append(item)


class MultiDataMap(DataMap):
    """
    Class representing a specialization of data-map, a collection of data
    products located on the same node, skippable as a set and individually
    """
    @DataMap.data.setter
    def data(self, data):
        if isinstance(data, DataMap):
            mdpdict = {}
            data.iterator = DataMap.SkipIterator
            for item in data:
                if not item.host in mdpdict:
                    mdpdict[item.host] = []
                mdpdict[item.host].append(item.file)
            mdplist = []
            for k, v in mdpdict.iteritems():
                mdplist.append(MultiDataProduct(k, v, False))
            self._set_data(mdplist, dtype=MultiDataProduct)
        else:
            #print 'HELP: ', data
            self._set_data(data, dtype=MultiDataProduct)

    def split_list(self, number):
        mdplist = []
        for item in self.data:
            for i in xrange(0, len(item.file), number):
                chunk = item.file[i:i+number]
                mdplist.append(MultiDataProduct(item.host, chunk, item.skip))
        self._set_data(mdplist)

# all test stuff down here
if __name__ == '__main__':
    descriptiontext = "This script lets you create Mapfiles for the Lofar Pipeline Framework"

    parser = argparse.ArgumentParser(description=descriptiontext, formatter_class=RawTextHelpFormatter)
    parser.add_argument('name', help='Give the name of the file to output')
    parser.add_argument('-n', '--number',help='Number of times to expand', type=int, default=0)
    #parser.add_argument('-v','--verbose',help='More detailed information',action='store_true')
    #parser.add_argument('-f','--faillog',help='Name of a file which will contain a list of failed commands from the list.',default=None)
    #parser.add_argument('-N','--NumberOfTasks',help='Number of concurrent commands.',type=int,default=0)
    #parser.add_argument('-l','--Logs',help='Individual log files for each process.',action='store_true')
    #parser.add_argument('-R','--retry',help='Number of times failed commands should be retried after all commands ran through',type=int,default=-1)
    #parser.add_argument('-L','--low',help='Low index of the commandlist. Start from here.',type=int,default=0)
    #parser.add_argument('-H','--high',help='High index of the commandlist. End execution at this index',type=int,default=None)
    args = parser.parse_args()

    mm = MapfileManager()
    #print 'MAP: ', mm.map
    #mm.expand(args.number)
    #mm.from_parts(ntimes=args.number)
    mm.from_parts(data=['d1','d2','d3'],ntimes=args.number)
    dp = DataProduct('i am', 'last', False)
    dmtest = DataMap([dp])
    mm.insert(2, {'host': 'i am', 'file': 'number two', 'skip': False})
    mm.append(dp)
    print 'MAP: ', mm.data
    mm.save(args.name)
    dm = DataMap.load(args.name)
    print 'LOADED: ', dm
    md = MultiDataProduct('localhost', dm, False)
    md2 = MultiDataProduct('foreignhost', dm, False)
    print 'MULTIprod', md
    mm.append(md)
    print 'BLA: ', mm.data
    mdm = MultiDataMap([md])
    print 'MULTIMAP: ', mdm
    mdm.split_list(1)
    print 'MULTIMAP SPLIT: ', mdm
    mm.delete(pattern=' two')
    print 'MM: ', mm
    mm2 = MapfileManager('/home/zam/sfroehli/PipelineExample/data')
    print 'FOLDERMAP: ', mm2