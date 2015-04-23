import sys
import pyrap.tables as pt
import os
import numpy

el=len(sys.argv)


mslist     = sys.argv[1:el-2]
inputms    = str(sys.argv[el-2])
outputcol  = str(sys.argv[el-1])


if os.path.isdir(inputms):
   datain     = pt.table(inputms)
   data       = datain.getcol('DATA', nrow=1)
   numberofchans = numpy.int(numpy.shape(data)[1])
   chanperms     = numberofchans/numpy.int(len(mslist))

else:
 print 'Nothing to do,', inputms, 'does not exist'

print 'TOTAL #chans & chans per ms', numberofchans, chanperms


for ms_id,ms in enumerate(mslist):

 if os.path.isdir(ms):

   data    = datain.getcolslice('DATA', [chanperms*ms_id,0], [(chanperms*(ms_id+1))-1,3])
   dataout = pt.table(ms, readonly=False)

   print '[blc],[trc]', [chanperms*ms_id,0], [(chanperms*(ms_id+1))-1,3]
   print ms, outputcol, chanperms*ms_id, chanperms*(ms_id+1), numpy.shape(data)
   ###dataout.putcol(outputcol, data[:,chanperms*ms_id:chanperms*(ms_id+1),:])
   dataout.putcol(outputcol, data)
   dataout.flush()
   dataout.close()



