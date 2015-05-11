import os
import sys
import numpy
import time

#casapy --nologger -c ft.py ms1 ms2 modimage nterms wplanes

execfile('mytasks.py')

el=len(sys.argv)

print sys.argv[:]


mslistfull = sys.argv[4:el-3]
modimage   = sys.argv[el-3]
ntermsi    = numpy.int(sys.argv[el-2])
wplanes    = numpy.int(sys.argv[el-1])


async     = True

if ntermsi == 1:
    mod = [modimage]
if ntermsi == 2:
    mod = [modimage+'.tt0',modimage+'.tt1']
if ntermsi == 3:
    mod = [modimage+'.tt0',modimage+'.tt1',modimage+'.tt2']


os.system('rm -rf *.tmpod')


if wplanes <= 127:
    nparallel = len(mslistfull)
if (wplanes > 127) and (wplanes <= 179):
    nparallel = 15
if (wplanes > 179) and (wplanes <= 195):
    nparallel = 14
if (wplanes > 195):
    nparallel = 6

groups = range(0,numpy.int(numpy.ceil(numpy.float(len(mslistfull))/numpy.float(nparallel))))
print numpy.float(len(mslistfull)),numpy.float(nparallel)




for group_id,group in enumerate(groups):
    print 'GROUP', group_id, group
    mslist = mslistfull[group_id*nparallel:nparallel*(group_id+1)]
    print mslist
    returnval = numpy.zeros(len(mslist))
    for ms_id,ms in enumerate(mslist):

        modimtmp = []
        for mod_id,modim in enumerate(mod): # to avoid file lock on model
            tmpname = modim +'.'+ ms +'.tmpod'
            os.system('cp -r ' + modim + ' ' + tmpname)
            modimtmp.append(tmpname)

        if wplanes > 1:
            returnval[ms_id] = ftw(vis=ms,field="",spw="",model=modimtmp,nterms=ntermsi,reffreq="",\
                                   wprojplanes=wplanes,complist="",incremental=False,usescratch=True, async=True)
        else:
            returnval[ms_id] = ft(vis=ms,field="",spw="",model=modimtmp,nterms=ntermsi,reffreq="",\
                                  complist="",incremental=False,usescratch=True, async=True)
        time.sleep(3)
        del modimtmp

    running = True
    while(running):
        running = False # set so that we go out here unless the part below sets it to True
        for ms_id,ms in enumerate(mslist):
            result = tm.retrieve(returnval[ms_id])
            print ms, result
            if result['status'] != 'done':
                running = True
        time.sleep(4)

os.system('rm -rf *.tmpod')

# 15% MEM 2048 --> max 6 bands
# 6.5% MEM 1600 --> max 14 bands
# 2.3% MEM 1200 --> no issue
# if wplanes > 179 but < 195 dan 1600
# if wplanes > 195 but < 255 dan 2048
