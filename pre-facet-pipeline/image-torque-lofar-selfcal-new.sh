#! /bin/csh

# send email alerts to this address
#PBS -M twshimwell@gmail.com
# mail alerts for end of job or aborted job
#PBS -m ea

# set name of job
#PBS -N P23

# set the error log file location
#PBS -j oe

# set max wallclock time and memory
#PBS -l walltime=40:00:00

# set the number of nodes and processes per node
#PBS -l nodes=1:ppn=12
#PBS -k oe

######################## END PBS ##############################################


# start job from the directory it was submitted
cd $PBS_O_WORKDIR

# source the lofar software
module load casa
module load lofar
source /soft/lofar-091114/lofarinit.csh
setenv PYTHONPATH /soft/pyrap:$PYTHONPATH
setenv PYTHONPATH /home/shimwell/lofarcodes-standalone:$PYTHONPATH
setenv PYTHONPATH /home/shimwell/lofarcodes-standalone/aux-lib:$PYTHONPATH
setenv PYTHONPATH /home/shimwell/survey_pipeline:$PYTHONPATH
setenv PYTHONPATH /home/shimwell/survey_pipeline/losoto/tools:$PYTHONPATH

./run_P1.sh > run_P1.out
