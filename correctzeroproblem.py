import numpy
import pyrap.tables as pt
import sys


msname     = str(sys.argv[1])
ref_msname = 'L99083_SB200_uv.dppp.MS'

t         = pt.table(msname, readonly=False)
t_ref     = pt.table(ref_msname, readonly=True)
time      = t.getcol('TIME')
time_ref  = t_ref.getcol('TIME')

print len(numpy.unique(time))
print len(numpy.unique(time_ref))


timelist = numpy.unique(time)

count = 0
for timeval in timelist:
  print count
  t1 = t.query('TIME > ' + str(timeval-1e-1)+'&&'+ 'TIME < ' + str(timeval+1e-1), columns='DATA')
  data = t1.getcol('DATA')
  print numpy.shape(data)
  idx    = numpy.where(numpy.abs(data) == 0.0)
  if len(idx) != 0:
    t1_ref  = t_ref.query('TIME > ' + str(timeval-1e-1)+'&&'+ 'TIME < ' + str(timeval+1e-1), columns='DATA')
    data_ref= t1_ref.getcol('DATA')
    if numpy.shape(data) != numpy.shape(data_ref):
      print 'Problem detected with data shape'
      sys.exit()
    data    = numpy.copy(data_ref)
    t1.putcol('DATA',data)
    t1_ref.close()
  t1.close() 
  count = count+1


t_ref.close()
t.close() 
