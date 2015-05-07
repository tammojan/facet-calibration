import os
import numpy
import sys
import glob
import pyrap.tables as pt
import lofar.parmdb


### STAR USER INPUT ###
ms                 = 'A2256_SB050-059.2ch10s.ms'
phase_only_parmdb  = 'A2256_SB050-059.2ch10s.ms/instrument_phaseonly' # contains

ap_parmdb_template = 'A2256_SB050-059.2ch10s.ms/instrument_ap_template'
                     # you need to make this one
                     # same(!) timegrid as phase_only_parmdb, needs keys Real:0:0, Imag:0:0
output_parmdb      = 'A2256_SB050-059.2ch10s.ms/instrument_ap_smoothed'
### END USER INPUT ###





# remove the output in case it exists
os.system('rm -rf ' + output_parmdb)


# make antenna list
anttab          = pt.table(ms + '/ANTENNA')
antenna_list    = anttab.getcol('NAME')
anttab.close()

pdb_phase   = lofar.parmdb.parmdb(phase_only_parmdb)
parms_phase = pdb_phase.getValuesGrid("*")

pdb_ap   = lofar.parmdb.parmdb(ap_parmdb_template)
parms_ap = pdb_ap.getValuesGrid("*")

# fill in the parmdb keys
for antenna_id, antenna in enumerate(antenna_list):
            
            dummyvectmp = numpy.copy(parms_ap['Gain:0:0:Real:' + antenna]['values'][:,0])
            dummyvec = 0.0*numpy.copy(dummyvectmp)
                
            phase_cal_00   = numpy.copy(parms_phase['Gain:0:0:Phase:' + antenna]['values'][:,0])
            phase_cal_11   = numpy.copy(parms_phase['Gain:1:1:Phase:' + antenna]['values'][:,0])

            real_00 = ((dummyvec) + (1.0))*numpy.cos(phase_cal_00)
            imag_00 = ((dummyvec) + (1.0))*numpy.sin(phase_cal_00)
                
            real_11 = ((dummyvec) + (1.0))*numpy.cos(phase_cal_11)
            imag_11 = ((dummyvec) + (1.0))*numpy.sin(phase_cal_11)

            parms_ap['Gain:0:0:Real:' + antenna]['values'][:,0] = numpy.copy(real_00)
            parms_ap['Gain:0:0:Imag:' + antenna]['values'][:,0] = numpy.copy(imag_00)        
            parms_ap['Gain:1:1:Real:' + antenna]['values'][:,0] = numpy.copy(real_11)
            parms_ap['Gain:1:1:Imag:' + antenna]['values'][:,0] = numpy.copy(imag_11)  


# write the parmdb
pdbnew = lofar.parmdb.parmdb(output_parmdb, create=True)
pdbnew.addValues(parms_ap)
pdbnew.flush()

# to make sure there are no issues with templates made at a different frequency
os.system("taql 'update " + output_parmdb + " set ENDX=1.e12'")
os.system("taql 'update " + output_parmdb + " set STARTX=1.0'")


### PARSET TO MAKE TEMPLATE "instrument_ap_template"
#Strategy.InputColumn = DATA
#Strategy.ChunkSize   = 800
#Strategy.Steps       = [solve]

#Step.solve.Model.Sources                = [] # all in skymodel 
#Step.solve.Model.Cache.Enable           = T
#Step.solve.Model.Phasors.Enable         = F
#Step.solve.Model.Gain.Enable            = T
#Step.solve.Operation                    = SOLVE
#Step.solve.Solve.Parms                  = ["Gain:0:0:*","Gain:1:1:*"] 
#Step.solve.Solve.CellSize.Freq          = 0
#Step.solve.Solve.CellSize.Time          = 1 
#Step.solve.Solve.CellChunkSize          = 100
#Step.solve.Solve.PropagateSolutions     = T
#Step.solve.Solve.Options.MaxIter        = 1
#Step.solve.Solve.Options.EpsValue       = 1e-9
#Step.solve.Solve.Options.EpsDerivative  = 1e-9
#Step.solve.Solve.Options.ColFactor      = 1e-9
#Step.solve.Solve.Options.LMFactor       = 1.0
#Step.solve.Solve.Options.BalancedEqs    = F
#Step.solve.Solve.Options.UseSVD         = T
#Step.solve.Solve.Mode                   = COMPLEX
#Step.solve.Model.Beam.Enable            = F
