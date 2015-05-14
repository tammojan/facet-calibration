import pyrap.tables as pt

tablename  = str(raw_input('Enter input msfile '))
outputname = str(raw_input('Enter output file name: '))
start_out = float(raw_input('Enter start time (hours) ')) #hours
end_out= float(raw_input('Enter end time (hours) '))

t = pt.table(tablename)
starttime = t[0]['TIME']
endtime   = t[t.nrows()-1]['TIME']
print '===================='
print 'Input Measurement Set is '+tablename
print 'Start time (sec) = '+str(starttime)
print 'End time   (sec) = '+str(endtime)
print 'Total time duration (hrs)  = '+str((endtime-starttime)/3600)
print '====================='
print 'Output Measurement Set is '+outputname
print 'Start time (relative to input ms start) = '+str(start_out)
print 'End time   (relative to input ms start) = '+str(end_out)
print 'Total time duration (hrs)  ='+str(end_out-start_out)
print '====================='
print 'Now going to do the Querry to select the required time range'

t1 = t.query('TIME > ' +str(starttime+start_out*3600) + ' && TIME < '  + str(starttime+end_out*3600), sortlist='TIME,ANTENNA1,ANTENNA2')

print 'Total rows in Input MS = '+str(t.nrows())
print 'Total rows in Output MS = '+str(t1.nrows())

print 'Now Writing the output MS'
t1.copy(outputname,True)
t1.close()
t.close()
print 'Copy completed'
