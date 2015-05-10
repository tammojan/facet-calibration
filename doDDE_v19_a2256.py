import matplotlib
#matplotlib.use('GTK')
import numpy
import os
import sys
from scipy import interpolate
import time
#import subprocess
from subprocess import Popen, PIPE
import pyrap.tables as pt
import pyrap.images
import pwd
import logging
pi       = numpy.pi
import sys


if len(sys.argv)<2:
    raise Exception('Give the path to the setup code for the facet')

print 'Using',sys.argv[1],'as the setup code'
execfile(sys.argv[1])
print 'script path is',SCRIPTPATH

try:
    StartAtStep
except NameError:
    print 'No starting step specified, begin at the beginning'
    StartAtStep='preSC'

try:
    WSCleanRobust
except NameError:
    WSCleanRobust=-0.25 # default preserves old value
    print 'No WSClean robust set, defaulting to',WSCleanRobust

try:
    BlankField
except NameError:
    BlankField=False
    print 'BlankField not set, defaulting to',BlankField

print 'StartAtStep is',StartAtStep

logging.basicConfig(filename='dde.log',level=logging.DEBUG, format='%(asctime)s -  %(message)s', datefmt='%Y-%d-%m %H:%M:%S')
logging.info('\n')

os.system('cp ' + SCRIPTPATH + '/coordinates_mode.py .')
os.system('cp ' + SCRIPTPATH + '/blank.py .')
os.system('cp ' + SCRIPTPATH + '/ftw.xml .')
os.system('cp ' + SCRIPTPATH + '/task_ftw.py .')

os.system(buildmytasks) # make casapy tasks


from coordinates_mode import *
import blank

# check high-DR



# TO DO
# - HIGH-DYNAMIC RANGE (need to adjust merger/join parmdb, solution smoohting is ok)
# - make freq averaging for selfcal fieldsize dependent (in do_dde (selfcal is not needed))
# - increase SNR in slow A&P by adding nearby blocks?
# - parallel instrument_merged in selfcal (with pp?)


# - smoothcal in parallel
# - selfcal switch do wsclean?, but can do predict in parallel?
# - nterms/multiscale > 1 in Facet imaging WSclean
# - convolve images all to the same resolution (use casapy for that, easy to do)
# - What about second imaging and calibration cycle?
#    - 1. ADD back skymodel using "master solutions", then CORRECT using "master solutions"
#    - 2. image that again using the same setting
#    - 3. redo the subtract (will be slightly better....but solutions remain the same, just better noise) or just proceed to the next field?


def find_newsize(mask):

    img    = pyrap.images.image(mask)
    pixels = numpy.copy(img.getdata())
    sh     = numpy.shape(pixels)[3:4]
    newsize = numpy.copy(sh[0])
    sh      = sh[0]
    print newsize

    trysizes = numpy.copy(sorted([6400,6144,5600,5400,5184,4800,4608,4320,4096,3840,3600,3200,3072,2880,2560,2304,2048, 1600, 1536, 1200, 1024, 800, 512]))
    idx = numpy.where(trysizes < sh)
    print idx
    trysizes = numpy.copy(trysizes[idx]) # remove sizes larger than image
    trysizes = numpy.copy(trysizes[::-1]) # reverse sorted
    print trysizes

    for size in trysizes:
        print 'Trying', size
        cutedge = numpy.int((sh - size)/2.)
        print cutedge

        idx1 = numpy.size(numpy.where(pixels[0,0,  0:cutedge,0:sh]      != 0))
        idx2 = numpy.size(numpy.where(pixels[0,0,  sh-cutedge:sh,0:sh]  != 0))
        idx3 = numpy.size(numpy.where(pixels[0,0,  0:sh,0:cutedge]      != 0))
        idx4 = numpy.size(numpy.where(pixels[0,0,  0:sh,sh-cutedge:sh]  != 0))

        print idx1, idx2, idx3, idx4

        if ((idx1) == 0) and ((idx2)  == 0) and ((idx3)  == 0) and ((idx4) == 0):
            # UPDATE THE IMAGE SIZE
            newsize  = numpy.copy(size)
            print 'Found new size', newsize, ' fits within the mask'

    return newsize



def runbbs(mslist, skymodel, parset, parmdb, replacesource):
    #NOTE WORK FROM MODEL_DATA (contains correct phase data from 10SB calibration)
    for ms in mslist:
        log      =  ms + '.bbslog'
        if replacesource:
            cmd = 'calibrate-stand-alone --replace-sourcedb --parmdb-name ' + parmdb + ' ' + ms + ' ' + parset + ' ' + skymodel + '>' + log + ' 2>&1 &'
        else:
            cmd = 'calibrate-stand-alone --parmdb-name ' + parmdb + ' ' + ms + ' ' + parset + ' ' + skymodel + '>' + log + ' 2>&1 &'
        print cmd
        os.system(cmd)
    time.sleep(10)

    done = 0
    while(done < len(mslist)):
        done = 0
        for ms in mslist:
            cmd = "grep 'bbs-reducer terminated successfully.' " + ms + ".bbslog"
            output=Popen(cmd, shell=True, stdout=PIPE).communicate()[0]
            if 'INFO' in output:
                done = done + 1
                print ms, 'is done'
        time.sleep(5)
    return

def runbbs16(mslist, skymodel, parset, parmdb, replacesource):
    #NOTE WORK FROM MODEL_DATA (contains correct phase data from 10SB calibration)

    mslist1 = mslist[0:15]
    mslist2 = mslist[15:len(mslist)]

    # PART1
    for ms in mslist1:
        log      =  ms + '.bbslog'
        if replacesource:
            cmd = 'calibrate-stand-alone --replace-sourcedb --parmdb-name ' + parmdb + ' ' + ms + ' ' + parset + ' ' + skymodel + '>' + log + ' 2>&1 &'
        else:
            cmd = 'calibrate-stand-alone --parmdb-name ' + parmdb + ' ' + ms + ' ' + parset + ' ' + skymodel + '>' + log + ' 2>&1 &'
        print cmd
        os.system(cmd)
    time.sleep(10)

    done = 0
    while(done < len(mslist1)):
        done = 0
        for ms in mslist1:
            cmd = "grep 'bbs-reducer terminated successfully.' " + ms + ".bbslog"
            output=Popen(cmd, shell=True, stdout=PIPE).communicate()[0]
            if 'INFO' in output:
                done = done + 1
                print ms, 'is done'
        time.sleep(5)


    # PART2
    for ms in mslist2:
        log      =  ms + '.bbslog'
        if replacesource:
            cmd = 'calibrate-stand-alone --replace-sourcedb --parmdb-name ' + parmdb + ' ' + ms + ' ' + parset + ' ' + skymodel + '>' + log + ' 2>&1 &'
        else:
            cmd = 'calibrate-stand-alone --parmdb-name ' + parmdb + ' ' + ms + ' ' + parset + ' ' + skymodel + '>' + log + ' 2>&1 &'
        print cmd
        os.system(cmd)
    time.sleep(10)

    done = 0
    while(done < len(mslist2)):
        done = 0
        for ms in mslist2:
            cmd = "grep 'bbs-reducer terminated successfully.' " + ms + ".bbslog"
            output=Popen(cmd, shell=True, stdout=PIPE).communicate()[0]
            if 'INFO' in output:
                done = done + 1
                print ms, 'is done'
        time.sleep(5)
    return

def create_subtract_parset_field_outlier(outputcolumn,TEC):
    bbs_parset = 'sub.parset'
    os.system('rm -f ' + bbs_parset)
    f=open(bbs_parset, 'w')

    chunksize = 200

    f.write('Strategy.InputColumn = MODEL_DATA\n')
    f.write('Strategy.ChunkSize   = %s\n' % chunksize)
    f.write('Strategy.UseSolver   = F\n')
    f.write('Strategy.Steps       = [subtract]\n\n\n')
    f.write('Step.subtract.Model.Sources                   = []\n')
    f.write('Step.subtract.Model.Cache.Enable              = T\n')
    f.write('Step.subtract.Model.Phasors.Enable            = F\n')
    f.write('Step.subtract.Model.DirectionalGain.Enable    = F\n')
    f.write('Step.subtract.Model.Gain.Enable               = T\n')
    f.write('Step.subtract.Model.Rotation.Enable           = F\n')
    f.write('Step.subtract.Model.CommonScalarPhase.Enable  = T\n')
    if TEC   == "True":
        f.write('Step.subtract.Model.TEC.Enable    = T\n')
    #if clock == "True":
    #   f.write('Step.subtract.Model.Clock.Enable  = T\n')
    f.write('Step.subtract.Model.CommonRotation.Enable     = F\n')
    f.write('Step.subtract.Operation                       = SUBTRACT\n')
    f.write('Step.subtract.Model.Beam.Enable               = F\n')
    f.write('Step.subtract.Output.WriteCovariance          = F\n')
    f.write('Step.subtract.Output.Column                   = %s\n' % outputcolumn)
    f.close()
    return bbs_parset



def runbbs_diffskymodel_addback(mslist, parmdb, replacesource, direction, imsize, output_template_im, do_ap):

    for ms in mslist:
        log      =  ms + '.bbslog'

        #set skymodel # ~weeren does not work in numpy.load
        skymodel =  ms.split('.')[0] + '.skymodel'

        # find sources to add back, make parset

        #callist, callistarraysources = cal_return_slist(output_template_im +'.masktmp',skymodel, direction, imsize)
        cmd = 'python ' + SCRIPTPATH + '/cal_return_slist.py '+ output_template_im +'.masktmp ' +skymodel +' "'+str(direction) +'" ' + str(imsize)
        output = Popen(cmd, shell=True, stdout=PIPE).communicate()[0]
        callist = output.strip()
        callistarraysources = callist.split(',')

        print 'Adding back for calibration:', callist

        if len(callist) != 0: # otherwise do not have to add
            parset = create_add_parset_ms(callist, ms, do_ap)
            if replacesource:
                cmd = 'calibrate-stand-alone --replace-sourcedb --parmdb-name ' + parmdb + ' ' + ms + ' ' + parset + ' ' + skymodel + '>' + log + ' 2>&1 &'
            else:
                cmd = 'calibrate-stand-alone --parmdb-name ' + parmdb + ' ' + ms + ' ' + parset + ' ' + skymodel + '>' + log + ' 2>&1 &'
            print cmd
            os.system(cmd)
            time.sleep(10)  # otherwise add.parset is deleted (takes time for BBS to start up)
        else:
            print 'No source to add back, are you sure the DDE position is correct?'
            os.system("taql 'update " + ms + " set ADDED_DATA_SOURCE=SUBTRACTED_DATA_ALL'")


    time.sleep(10)
    done = 0
    while(done < len(mslist)):
        done = 0
        for ms in mslist:
            cmd = "grep 'bbs-reducer terminated successfully.' " + ms + ".bbslog"
            output=Popen(cmd, shell=True, stdout=PIPE).communicate()[0]
            if 'INFO' in output:
                done = done + 1
                print ms, 'is done'
        time.sleep(5)
    return

def runbbs_diffskymodel_addback16(mslist, parmdb, replacesource, direction, imsize, output_template_im, do_ap):
    mslist1 = mslist[0:16]
    mslist2 = mslist[16:len(mslist)]
    #part1
    for ms in mslist1:
        log      =  ms + '.bbslog'

        #set skymodel # ~weeren does not work in numpy.load
        #skymodel =  ms.split('.')[0] + '.skymodel'
        skymodel = ms.split('.')[0] + '.skymodel'


        # find sources to add back, make parset
        #callist, callistarraysources = cal_return_slist(output_template_im +'.masktmp',skymodel, direction, imsize)
        cmd = 'python ' + SCRIPTPATH + '/cal_return_slist.py '+ output_template_im +'.masktmp ' +skymodel +' "'+str(direction) +'" ' + str(imsize)
        output = Popen(cmd, shell=True, stdout=PIPE).communicate()[0]
        callist = output.strip()
        callistarraysources = callist.split(',')

        if len(callist) != 0: # otherwise do not have to add
            print 'Adding back for calibration:', callist
            parset = create_add_parset_ms(callist, ms, do_ap)

            if replacesource:
                cmd = 'calibrate-stand-alone --replace-sourcedb --parmdb-name ' + parmdb + ' ' + ms + ' ' + parset + ' ' + skymodel + '>' + log + ' 2>&1 &'
            else:
                cmd = 'calibrate-stand-alone --parmdb-name ' + parmdb + ' ' + ms + ' ' + parset + ' ' + skymodel + '>' + log + ' 2>&1 &'
            print cmd
            os.system(cmd)
            time.sleep(10)  # otherwise add.parset is deleted (takes time for BBS to start up)
        else:
            if do_ap:
                print 'do_ap not supported'
                sys.exit()
            print 'No source to add back, are you sure the DDE position is correct?'
            os.system("taql 'update " + ms + " set ADDED_DATA_SOURCE=SUBTRACTED_DATA_ALL'")

    time.sleep(10)
    done = 0
    while(done < len(mslist1)):
        done = 0
        for ms in mslist1:
            cmd = "grep 'bbs-reducer terminated successfully.' " + ms + ".bbslog"
            output=Popen(cmd, shell=True, stdout=PIPE).communicate()[0]
            if 'INFO' in output:
                done = done + 1
                print ms, 'is done'
        time.sleep(5)

    #part2
    for ms in mslist2:
        log      =  ms + '.bbslog'

        #set skymodel # ~weeren does not work in numpy.load
        skymodel =  ms.split('.')[0] + '.skymodel'

        # find sources to add back, make parset
        #callist, callistarraysources = cal_return_slist(output_template_im +'.masktmp',skymodel, direction, imsize)
        cmd = 'python '+ SCRIPTPATH +'/cal_return_slist.py '+ output_template_im +'.masktmp ' +skymodel +' "'+str(direction) +'" ' + str(imsize)
        output = Popen(cmd, shell=True, stdout=PIPE).communicate()[0]
        callist = output.strip()
        callistarraysources = callist.split(',')

        print 'Adding back for calibration:', callist

        if len(callist) != 0: # otherwise do not have to add
            parset = create_add_parset_ms(callist, ms, do_ap)
            if replacesource:
                cmd = 'calibrate-stand-alone --replace-sourcedb --parmdb-name ' + parmdb + ' ' + ms + ' ' + parset + ' ' + skymodel + '>' + log + ' 2>&1 &'
            else:
                cmd = 'calibrate-stand-alone --parmdb-name ' + parmdb + ' ' + ms + ' ' + parset + ' ' + skymodel + '>' + log + ' 2>&1 &'
            print cmd
            os.system(cmd)
            time.sleep(10)  # otherwise add.parset is deleted (takes time for BBS to start up)
        else:
            if do_ap:
                print 'do_ap not supported'
                sys.exit()
            print 'No source to add back, are you sure the DDE position is correct?'
            os.system("taql 'update " + ms + " set ADDED_DATA_SOURCE=SUBTRACTED_DATA_ALL'")

    time.sleep(10)
    done = 0
    while(done < len(mslist2)):
        done = 0
        for ms in mslist2:
            cmd = "grep 'bbs-reducer terminated successfully.' " + ms + ".bbslog"
            output=Popen(cmd, shell=True, stdout=PIPE).communicate()[0]
            if 'INFO' in output:
                done = done + 1
                print ms, 'is done'
        time.sleep(5)
    return


def runbbs_diffskymodel_addbackfield(mslist, parmdb, replacesource, direction, imsize, output_template_im, do_ap):

    for ms in mslist:
        log      =  ms + '.bbslog'

        #set skymodel
        skymodel =  ms.split('.')[0] + '.skymodel'

        # find peeling sources (from previous step)
        #callist, callistarraysources = cal_return_slist(output_template_im +'.masktmp',skymodel, direction, imsize)
        cmd = 'python '+ SCRIPTPATH + '/cal_return_slist.py '+ output_template_im +'.masktmp ' +skymodel +' "'+str(direction) +'" ' + str(imsize)
        output = Popen(cmd, shell=True, stdout=PIPE).communicate()[0]
        callist = output.strip()
        callistarraysources = callist.split(',')

        print 'Add field back step 1'

        # return the source list from the source to be added back sourrinding the peeling source and which fall within the mask boundaries
        # put in MODEL_DATA


        #addback_sourcelist = return_slist(output_template_im +'.masktmp', skymodel, callistarraysources)
        cmd = 'python '+ SCRIPTPATH +'/return_slist.py '+ output_template_im +'.masktmp ' +skymodel +' "'+str(callist)+'"'
        output = Popen(cmd, shell=True, stdout=PIPE).communicate()[0]
        addback_sourcelist = output.strip()




        print 'Field source added back are: ', addback_sourcelist

        if len(addback_sourcelist) != 0: # otherwise do not have to add
            parset = create_add_parset_field_ms(addback_sourcelist, ms, do_ap)
            if replacesource:
                cmd = 'calibrate-stand-alone --replace-sourcedb --parmdb-name ' + parmdb + ' ' + ms + ' ' + parset + ' ' + skymodel + '>' + log + ' 2>&1 &'
            else:
                cmd = 'calibrate-stand-alone --parmdb-name ' + parmdb + ' ' + ms + ' ' + parset + ' ' + skymodel + '>' + log + ' 2>&1 &'
            print cmd
            os.system(cmd)
        else:
            os.system("taql 'update " + ms + " set MODEL_DATA=ADDED_DATA_SOURCE'") # in case no sources are put back

        time.sleep(10) # otherwise addfield.parset is deleted (takes time for BBS to start up)


    time.sleep(10)
    done = 0
    while(done < len(mslist)):
        done = 0
        for ms in mslist:
            cmd = "grep 'bbs-reducer terminated successfully.' " + ms + ".bbslog"
            output=Popen(cmd, shell=True, stdout=PIPE).communicate()[0]
            if 'INFO' in output:
                done = done + 1
                print ms, 'is done'
        time.sleep(5)
    return


def runbbs_2(mslist, msparmdb, skymodel, parset, parmdb):
    for ms_id, ms in enumerate(mslist):
        log      =  ms + '.bbslog'
        cmd = 'calibrate-stand-alone --parmdb ' + msparmdb[ms_id]+'/'+parmdb + ' ' + ms + ' ' + parset + ' ' + skymodel + '>' + log + ' 2>&1 &'
        print cmd
        os.system(cmd)
    time.sleep(10)

    done = 0
    while(done < len(mslist)):
        done = 0
        for ms in mslist:
            cmd = "grep 'bbs-reducer terminated successfully.' " + ms + ".bbslog"
            output=Popen(cmd, shell=True, stdout=PIPE).communicate()[0]
            if 'INFO' in output:
                done = done + 1
                print ms, 'is done'
        time.sleep(5)
    return


def create_phaseshift_parset_full(msin, msout, direction, column):
    ndppp_parset = 'ndppp_phaseshiftfull.parset'
    os.system('rm -f ' + ndppp_parset)
    f=open(ndppp_parset, 'w')
    f.write('msin ="%s"\n' % msin)
    f.write('msin.datacolumn = "%s"\n' % column)
    f.write('msin.autoweight = false\n')
    f.write('msout ="%s"\n' % msout)
    f.write('msout.writefullresflag=False\n')
    f.write('steps = [shift]\n')
    f.write('shift.type        = phaseshift\n')
    f.write('shift.phasecenter = [%s]\n' % direction)
    f.close()
    return ndppp_parset


def create_phaseshift_parset(msin, msout, source, direction, imsize, dynamicrange):
    ndppp_parset = (msin.split('.')[0]) +'_ndppp_avgphaseshift.parset'
    os.system('rm -f ' + ndppp_parset)

    f=open(ndppp_parset, 'w')
    f.write('msin ="%s"\n' % msin)
    f.write('msin.datacolumn = ADDED_DATA_SOURCE\n')
    f.write('msin.autoweight = false\n')
    f.write('msout ="%s"\n' % msout)
    f.write('msout.writefullresflag=False\n')
    f.write('steps = [shift,avg1]\n')
    f.write('shift.type        = phaseshift\n')
    f.write('shift.phasecenter = [%s]\n' % direction)
    f.write('avg1.type = squash\n')
    f.write('avg1.freqstep = 20\n')
    f.write('avg1.timestep = 1\n')
    f.close()
    return ndppp_parset



def create_phaseshift_parset2(msin, msout, source, direction, imsize, dynamicrange):
    ndppp_parset = (msin.split('.')[0]) +'_ndppp_avgphaseshift.parset'
    os.system('rm -f ' + ndppp_parset)

    f=open(ndppp_parset, 'w')
    f.write('msin ="%s"\n' % msin)
    f.write('msin.datacolumn = ADDED_DATA_SOURCE\n')
    f.write('msin.autoweight = false\n')
    f.write('msout ="%s"\n' % msout)
    f.write('msout.writefullresflag=False\n')
    f.write('steps = [shift,avg1]\n')
    f.write('shift.type        = phaseshift\n')
    f.write('shift.phasecenter = [%s]\n' % direction)
    f.write('avg1.type = squash\n')

    if dynamicrange == 'LD':
        if imsize <= 800:
            f.write('avg1.freqstep = 20\n')
        else:
            if imsize <= 1200:
                f.write('avg1.freqstep = 10\n')
            else:
                f.write('avg1.freqstep = 5\n')  # we have a large image  2048 is more or less max expected
    else:
        if dynamicrange != 'HD':
            print 'dynamicrange ', dynamicrange
            raise Exception('Wrong dynamicrange code, use "LD" or "HD"')
        print 'High dynamic range DDE cycle, eveything with be slow...'
        f.write('avg1.freqstep = 2\n') # one channel per SB



    f.write('avg1.timestep = 1\n')
    f.close()
    return ndppp_parset

def create_phaseshift_parset_formasks(msin, msout, source, direction):
    ndppp_parset = (msin.split('.')[0]) +'_ndppp_avgphaseshift.parset'
    os.system('rm -f ' + ndppp_parset)

    f=open(ndppp_parset, 'w')
    f.write('msin ="%s"\n' % msin)
    f.write('msin.datacolumn = DATA\n')
    f.write('msin.autoweight = False\n')
    f.write('msin.baseline   = 0&1\n') # only one baseline
    f.write('msout ="%s"\n' % msout)
    f.write('msout.writefullresflag=False\n')
    f.write('steps = [shift,avg1]\n')
    f.write('shift.type        = phaseshift\n')
    f.write('shift.phasecenter = [%s]\n' % direction)
    f.write('avg1.type = squash\n')
    f.write('avg1.freqstep = 1\n')
    f.write('avg1.timestep = 1\n')
    f.close()
    return ndppp_parset


def create_phaseshift_parset_field(msin, msout, source, direction):
    ndppp_parset = msin.split('.')[0] +'ndppp_avgphaseshift_field.parset'
    os.system('rm -f ' + ndppp_parset)

    f=open(ndppp_parset, 'w')
    f.write('msin ="%s"\n' % msin)
    f.write('msin.datacolumn = CORRECTED_DATA\n')
    f.write('msin.autoweight = false\n')
    f.write('msout ="%s"\n' % msout)
    f.write('msout.writefullresflag=False\n')
    f.write('steps = [shift,avg1]\n')
    f.write('shift.type        = phaseshift\n')
    f.write('shift.phasecenter = [%s]\n' % direction)
    f.write('avg1.type = squash\n')
    f.write('avg1.freqstep = 5\n')
    f.write('avg1.timestep = 3\n')
    f.close()
    return ndppp_parset




def create_add_parset_ms(source, ms, do_ap):
    bbs_parset = ms + '_add.parset'
    os.system('rm -f ' + bbs_parset)
    f=open(bbs_parset, 'w')
    f.write('Strategy.InputColumn = SUBTRACTED_DATA_ALL\n')
    f.write('Strategy.ChunkSize   = 200\n')
    f.write('Strategy.UseSolver   = F\n')
    f.write('Strategy.Steps       = [add]\n\n\n')
    f.write('Step.add.Model.Sources                   = [%s]\n' % source)
    f.write('Step.add.Model.Cache.Enable              = T\n')
    f.write('Step.add.Model.Phasors.Enable            = F\n')
    f.write('Step.add.Model.DirectionalGain.Enable    = F\n')
    if do_ap:
        f.write('Step.add.Model.Gain.Enable               = T\n')
    else:
        f.write('Step.add.Model.Gain.Enable               = F\n')
    f.write('Step.add.Model.Rotation.Enable           = F\n')
    f.write('Step.add.Model.CommonScalarPhase.Enable  = F\n')
    f.write('Step.add.Model.CommonRotation.Enable     = F\n')
    f.write('Step.add.Operation                       = ADD\n')
    f.write('Step.add.Model.Beam.Enable               = F\n')
    f.write('Step.add.Output.WriteCovariance          = F\n')
    f.write('Step.add.Output.Column                   = ADDED_DATA_SOURCE\n')
    f.close()
    return bbs_parset


def create_add_parset_field(source):
    bbs_parset = 'addfield.parset'
    os.system('rm -f ' + bbs_parset)
    f=open(bbs_parset, 'w')
    f.write('Strategy.InputColumn = ADDED_DATA_SOURCE\n') # already contains peeling source
    f.write('Strategy.ChunkSize   = 200\n')
    f.write('Strategy.UseSolver   = F\n')
    f.write('Strategy.Steps       = [add]\n\n\n')
    f.write('Step.add.Model.Sources                   = [%s]\n' % source)
    f.write('Step.add.Model.Cache.Enable              = T\n')
    f.write('Step.add.Model.Phasors.Enable            = F\n')
    f.write('Step.add.Model.DirectionalGain.Enable    = F\n')
    f.write('Step.add.Model.Gain.Enable               = T\n')
    f.write('Step.add.Model.Rotation.Enable           = F\n')
    f.write('Step.add.Model.CommonScalarPhase.Enable  = F\n')
    f.write('Step.add.Model.CommonRotation.Enable     = F\n')
    f.write('Step.add.Operation                       = ADD\n')
    f.write('Step.add.Model.Beam.Enable               = F\n')
    f.write('Step.add.Output.WriteCovariance          = F\n')
    f.write('Step.add.Output.Column                   = MODEL_DATA\n') # use use to save disk space
    f.close()
    return bbs_parset

def create_add_parset_field_ms(source, ms, do_ap):
    bbs_parset = ms + '_addfield.parset'
    os.system('rm -f ' + bbs_parset)
    f=open(bbs_parset, 'w')
    f.write('Strategy.InputColumn = ADDED_DATA_SOURCE\n') # already contains peeling source
    f.write('Strategy.ChunkSize   = 200\n')
    f.write('Strategy.UseSolver   = F\n')
    f.write('Strategy.Steps       = [add]\n\n\n')
    f.write('Step.add.Model.Sources                   = [%s]\n' % source)
    f.write('Step.add.Model.Cache.Enable              = T\n')
    f.write('Step.add.Model.Phasors.Enable            = F\n')
    f.write('Step.add.Model.DirectionalGain.Enable    = F\n')
    if do_ap:
        f.write('Step.add.Model.Gain.Enable               = T\n')
    else:
        f.write('Step.add.Model.Gain.Enable               = F\n')
    f.write('Step.add.Model.Rotation.Enable           = F\n')
    f.write('Step.add.Model.CommonScalarPhase.Enable  = F\n')
    f.write('Step.add.Model.CommonRotation.Enable     = F\n')
    f.write('Step.add.Operation                       = ADD\n')
    f.write('Step.add.Model.Beam.Enable               = F\n')
    f.write('Step.add.Output.WriteCovariance          = F\n')
    f.write('Step.add.Output.Column                   = MODEL_DATA\n') # use use to save disk space
    f.close()
    return bbs_parset

def create_subtract_parset(outputcolumn):
    bbs_parset = 'sub.parset'
    os.system('rm -f ' + bbs_parset)
    f=open(bbs_parset, 'w')
    f.write('Strategy.InputColumn = ADDED_DATA_SOURCE\n')
    f.write('Strategy.ChunkSize   = 100\n')
    f.write('Strategy.UseSolver   = F\n')
    f.write('Strategy.Steps       = [subtract]\n\n\n')
    f.write('Step.subtract.Model.Sources                   = []\n')
    f.write('Step.subtract.Model.Cache.Enable              = T\n')
    f.write('Step.subtract.Model.Phasors.Enable            = F\n')
    f.write('Step.subtract.Model.DirectionalGain.Enable    = F\n')
    f.write('Step.subtract.Model.Gain.Enable               = T\n')
    f.write('Step.subtract.Model.Rotation.Enable           = F\n')
    f.write('Step.subtract.Model.CommonScalarPhase.Enable  = T\n')
    f.write('Step.subtract.Model.CommonRotation.Enable     = F\n')
    f.write('Step.subtract.Operation                       = SUBTRACT\n')
    f.write('Step.subtract.Model.Beam.Enable               = F\n')
    f.write('Step.subtract.Output.WriteCovariance          = F\n')
    f.write('Step.subtract.Output.Column                   = %s\n' % outputcolumn)
    f.close()
    return bbs_parset

def create_subtract_parset_field(outputcolumn,TEC):
    bbs_parset = 'sub.parset'
    os.system('rm -f ' + bbs_parset)
    f=open(bbs_parset, 'w')

    chunksize = 175

    f.write('Strategy.InputColumn = MODEL_DATA\n')
    f.write('Strategy.ChunkSize   = %s\n' % chunksize)
    f.write('Strategy.UseSolver   = F\n')
    f.write('Strategy.Steps       = [subtract]\n\n\n')
    f.write('Step.subtract.Model.Sources                   = [@ADDED_DATA_SOURCE]\n')
    f.write('Step.subtract.Model.Cache.Enable              = T\n')
    f.write('Step.subtract.Model.Phasors.Enable            = F\n')
    f.write('Step.subtract.Model.DirectionalGain.Enable    = F\n')
    f.write('Step.subtract.Model.Gain.Enable               = T\n')
    f.write('Step.subtract.Model.Rotation.Enable           = F\n')
    f.write('Step.subtract.Model.CommonScalarPhase.Enable  = T\n')
    if TEC   == "True":
        f.write('Step.subtract.Model.TEC.Enable  = T\n')
    #if clock == "True":
    #   f.write('Step.subtract.Model.Clock.Enable  = T\n')
    f.write('Step.subtract.Model.CommonRotation.Enable     = F\n')
    f.write('Step.subtract.Operation                       = SUBTRACT\n')
    f.write('Step.subtract.Model.Beam.Enable               = F\n')
    f.write('Step.subtract.Output.WriteCovariance          = F\n')
    f.write('Step.subtract.Output.Column                   = %s\n' % outputcolumn)
    f.close()
    return bbs_parset

def join_parmdb(ms, parmdb_selfcal,parmdb_nondde, parmdb_template, parmdb_out, TEC, clock):
    import lofar.parmdb
    pdb_s = lofar.parmdb.parmdb(parmdb_selfcal)
    pdb_p = lofar.parmdb.parmdb(parmdb_nondde)
    pdb_t = lofar.parmdb.parmdb(parmdb_template)

    parms_s = pdb_s.getValuesGrid("*")
    parms_p = pdb_p.getValuesGrid("*")
    parms_t = pdb_t.getValuesGrid("*")

    keynames = parms_s.keys()
    os.system('rm -rf ' + parmdb_out)

    for key in keynames:
                # copy over the selfcal solutions, can copy all (Real, Imag, CommonScalarPhase, TEC, clock)
        parms_t[key]['values'][:,0] = numpy.copy(parms_s[key]['values'][:,0])

    pol_list = ['0:0','1:1']
    gain     = 'Gain'
    anttab     = pt.table(ms + '/ANTENNA')
    antenna_list    = anttab.getcol('NAME')
    anttab.close()

    for pol in pol_list:
        for antenna in antenna_list:

            real1 = numpy.copy(parms_p[gain + ':' + pol + ':Real:'+ antenna]['values'][:, 0])
            imag1 = numpy.copy(parms_p[gain + ':' + pol + ':Imag:'+ antenna]['values'][:, 0])

            real2 = numpy.copy(parms_s[gain + ':' + pol + ':Real:'+ antenna]['values'][:, 0])
            imag2 = numpy.copy(parms_s[gain + ':' + pol + ':Imag:'+ antenna]['values'][:, 0])


            G1 = real1 + 1j*imag1
            G2 = real2 + 1j*imag2

            #G_new = G_nondde*G_selfcal

            if TEC == "True":
                Gnew = numpy.copy(G2)
            else:
                Gnew = numpy.copy(G1*G2)


            parms_t[gain + ':' + pol + ':Imag:'+ antenna]['values'][:, 0] = numpy.copy(numpy.imag(Gnew))
            parms_t[gain + ':' + pol + ':Real:'+ antenna]['values'][:, 0] = numpy.copy(numpy.real(Gnew))


    #lofar.expion.parmdbmain.store_parms(parmdb_out, parms_t, create_new = True)
    pdbnew = lofar.parmdb.parmdb(parmdb_out, create=True)
    pdbnew.addValues(parms_t)
    pdbnew.flush()
    return

def normalize_parmdbs(mslist, parmdbname, parmdboutname):

    import lofar.parmdb
    amplist = []

    # create antenna list
    pol_list        = ['0:0','1:1']
    gain            = 'Gain'
    anttab          = pt.table(mslist[0] + '/ANTENNA')
    antenna_list    = anttab.getcol('NAME')
    anttab.close()

    for ms in mslist:

        pdb = lofar.parmdb.parmdb(ms + '/' + parmdbname)
        parms = pdb.getValuesGrid("*")

        for pol in pol_list:
            for antenna in antenna_list:
                real = numpy.copy(parms[gain + ':' + pol + ':Real:'+ antenna]['values'][:, 0])
                imag = numpy.copy(parms[gain + ':' + pol + ':Imag:'+ antenna]['values'][:, 0])
                amp  = numpy.copy(numpy.sqrt(real**2 + imag**2))
                amplist.append(amp)

    norm_factor = 1./(numpy.mean(amplist))
    print 'Normalizing gains: average gain value is', 1./norm_factor
    print 'Multiplying gains by:', norm_factor

    if (norm_factor > 1.5) or (norm_factor < (1./1.5)):
        print 'Check normalization'
        raise Exception('Wrong normalization')

    # now normalize the parmdbs
    for ms in mslist:
        pdb = lofar.parmdb.parmdb(ms + '/' + parmdbname)
        parms = pdb.getValuesGrid("*")

        os.system('rm -rf ' + ms + '/' + parmdboutname)

        for pol in pol_list:
            for antenna in antenna_list:
                real = numpy.copy(parms[gain + ':' + pol + ':Real:'+ antenna]['values'][:, 0])
                imag = numpy.copy(parms[gain + ':' + pol + ':Imag:'+ antenna]['values'][:, 0])
                parms[gain + ':' + pol + ':Imag:'+ antenna]['values'][:, 0] = numpy.copy(imag*norm_factor)
                parms[gain + ':' + pol + ':Real:'+ antenna]['values'][:, 0] = numpy.copy(real*norm_factor)

        pdbnew = lofar.parmdb.parmdb(ms + '/' +parmdboutname, create=True)
        pdbnew.addValues(parms)
        pdbnew.flush()

    return numpy.mean(amplist)



def return_slist(imagename, skymodel, ref_source):

    fluxweight = False

    data = load_bbs_skymodel(skymodel)

    if len(numpy.shape(data)) == 1:  # in this case not issue and we do not use Pythonlibs
        patchest,ra_patches,dec_patches, flux_patches =  compute_patch_center(data,fluxweight)
        print 'option 1'
    if len(numpy.shape(data)) == 2:
        patchest,ra_patches,dec_patches, flux_patches = compute_patch_center_libsproblem(data,fluxweight)
        print 'option 2'

    # remove sources already in the field and convert to radians

    if len(ref_source) == 1:
        idx = numpy.where(patchest != ref_source)
        ralist  = pi*(ra_patches[idx])/180.
        declist = pi*(dec_patches[idx])/180.
        patches = patchest[idx]
    else:
        idx = numpy.asarray([numpy.where(patchest == y)[0][0] for y in ref_source])
        accept_idx = sorted(set(range(patchest.size)) - set(idx))

        ralist  = pi*(ra_patches[accept_idx])/180.
        declist = pi*(dec_patches[accept_idx])/180.
        patches = patchest[accept_idx]


    img    = pyrap.images.image(imagename)
    pixels = numpy.copy(img.getdata())
    plist = []
    sh    = numpy.shape(pixels)[2:4]


    for patch_id,patch in enumerate(patches):
        coor = [0,1,declist[patch_id],ralist[patch_id]]
        pix  = img.topixel(coor)[2:4]

        if (pix[0] >= 0) and (pix[0] <= (sh[0]-1)) and \
           (pix[1] >= 0) and (pix[1] <= (sh[1]-1)):
            if pixels[0,0,pix[0],pix[1]] != 0.0:  # only include if withtin the clean mask (==1)
                plist.append(patches[patch_id])

    sourcess = ''
    if len(plist) == 1:
        sourcess = str(plist[0])
    else:
        for patch in plist:
            sourcess = sourcess+patch+','
        sourcess = sourcess[:-1]

    return sourcess

def angsep(ra1deg, dec1deg, ra2deg, dec2deg):
    """Returns angular separation between two coordinates (all in degrees)"""
    import math

    ra1rad=ra1deg*math.pi/180.0
    dec1rad=dec1deg*math.pi/180.0
    ra2rad=ra2deg*math.pi/180.0
    dec2rad=dec2deg*math.pi/180.0

    # calculate scalar product for determination
    # of angular separation
    x=math.cos(ra1rad)*math.cos(dec1rad)*math.cos(ra2rad)*math.cos(dec2rad)
    y=math.sin(ra1rad)*math.cos(dec1rad)*math.sin(ra2rad)*math.cos(dec2rad)
    z=math.sin(dec1rad)*math.sin(dec2rad)

    if x+y+z >= 1: rad = 0
    else: rad=math.acos(x+y+z)

    # Angular separation
    deg=rad*180/math.pi
    return deg


def cal_return_slist(imagename,skymodel, direction, imsize):

    factor = 0.8 # only add back in the center 80%
    cut = 1.5*(imsize/2.)*factor/3600.

    ra  = direction.split(',')[0]
    dec = direction.split(',')[1]

    ra1   = float(ra.split('h')[0])*15.
    ratmp = (ra.split('h')[1])
    ra2   = float(ratmp.split('m')[0])*15./60
    ra3   = float(ratmp.split('m')[1])*15./3600.
    ref_ra= ra1 + ra2 +ra3

    dec1   = float(dec.split('d')[0])
    dectmp = (dec.split('d')[1])
    dec2   = float(dectmp.split('m')[0])/60
    dec3   = float(dectmp.split('m')[1])/3600.
    ref_dec= dec1 + dec2 +dec3

    fluxweight = False


    data = load_bbs_skymodel(skymodel)

    if len(numpy.shape(data)) == 1:  # in this case not issue and we do not use Pythonlibs
        patches,ra_patches,dec_patches, flux_patches =  compute_patch_center(data,fluxweight)
        print 'option 1'
    if len(numpy.shape(data)) == 2:
        patches,ra_patches,dec_patches, flux_patches = compute_patch_center_libsproblem(data,fluxweight)
        print 'option 2'

    ralist  = pi*(ra_patches)/180.
    declist = pi*(dec_patches)/180.


    plist = []

    print ref_ra, ref_dec



    # load image to check if source within boundaries
    img    = pyrap.images.image(imagename)
    pixels = numpy.copy(img.getdata())
    plist = []
    sh    = numpy.shape(pixels)[2:4]


    # CHECK TWO THINGS
    #  - sources fall within the image size
    #  - sources fall within the mask from the tessellation
    for patch_id,patch in enumerate(patches):
        coor = [0,1,declist[patch_id],ralist[patch_id]]
        pix  = img.topixel(coor)[2:4]

        # compute radial distance to image center
        dis = angsep(ra_patches[patch_id],dec_patches[patch_id], ref_ra, ref_dec)

        if dis < cut: # ok sources is within image
            # check if the sources is within the mask region (because mask can be smaller than image)
            if (pix[0] >= 0) and (pix[0] <= (sh[0]-1)) and \
               (pix[1] >= 0) and (pix[1] <= (sh[1]-1)):
                if pixels[0,0,pix[0],pix[1]] != 0.0:  # only include if withtin the clean mask (==1)
                    plist.append(patches[patch_id])

    # make the string type source list
    sourcess = ''
    if len(plist) == 1:
        sourcess = str(plist[0])
    else:
        for patch in plist:
            sourcess = sourcess+patch+','
        sourcess = sourcess[:-1]

    return sourcess, plist



def make_image(mslist, cluster, callnumber, threshpix, threshisl, nterms, atrous_do, imsize, inputmask, mscale, region):

    niter   = numpy.int(2000 * (numpy.sqrt(numpy.float(len(mslist)))))

    depth =  0.7 / (numpy.sqrt(numpy.float(len(mslist))))
    cleandepth1 = str(depth*1.5) + 'mJy'
    cleandepth2 = str(depth)    + 'mJy'

    # speed up the imaging if possible by reducing image size within the mask region
    newsize = find_newsize(inputmask)
    if newsize < imsize: # ok so we can use a smaller image size then
        #make a new template
        os.system('casapy --nologger -c '+ SCRIPTPATH +'/make_empty_image.py '+ str(mslist[0]) + ' ' + inputmask+'2' + ' ' + str(newsize) + ' ' +'1.5arcsec')
        os.system('casapy --nologger -c '+ SCRIPTPATH +'/regrid_image.py '    + inputmask      + ' ' + inputmask+'2' + ' ' + inputmask+'3')

        # reset the imsize and the mask
        imsize    = newsize
        inputmask = inputmask+'3'


    ms = ''
    for m in mslist:
        ms = ms + ' ' + m

    imout = 'im'+ callnumber +'_cluster'+cluster+'nm'

    os.system('casapy --nologger -c ' + SCRIPTPATH + '/casapy_cleanv4.py ' + ms + ' ' + imout + ' ' + 'None' +\
               ' ' + cleandepth1 + ' ' + str(niter) + ' ' + str(nterms) + ' ' + str(imsize) + ' ' + mscale)


    # make mask
    if nterms > 1:
        os.system('python ' + SCRIPTPATH +'/makecleanmask_field.py --threshpix '+str(threshpix)+\
                    ' --threshisl '+str(threshisl) +' --atrous_do '+ str(atrous_do) +' '   +imout +'.image.tt0')
    else:
        os.system('python ' + SCRIPTPATH +'/makecleanmask_field.py --threshpix '+str(threshpix)+\
                  ' --threshisl '+str(threshisl) +' --atrous_do '+ str(atrous_do) + ' '  + imout +'.image')

    mask_sources = imout+'.cleanmask'
    os.system('rm -rf ' + mask_sources + 'field')
    os.system('cp -r ' + mask_sources + ' ' + mask_sources + 'field')

    #Merge the two masks
    img    = pyrap.images.image(mask_sources+'field')
    pixels = numpy.copy(img.getdata())

    img2    = pyrap.images.image(inputmask)
    pixels2 = numpy.copy(img2.getdata())

    idx = numpy.where(pixels2 == 0.0)
    pixels[idx] = 0.0
    img.putdata(pixels)
    img.unlock()
    del img
    del img2

    niter = numpy.int(niter*1.2) # clean a bit deeper (will actually be quite a bit deeper because of mask)
    imout = 'im'+ callnumber +'_cluster'+cluster


    if region != 'empty': # in that case we have a extra region file for the clean mask
        niter = niter*3 # increase niter, tune manually if needed
        os.system('casapy --nologger -c ' + SCRIPTPATH +'/casapy_cleanv4.py '+ ms + ' ' + imout + ' ' + mask_sources+'field,'+region + \
                  ' ' + cleandepth2 + ' ' + str(niter) + ' ' + str(nterms) + ' ' + str(imsize) + ' ' + mscale)

    else:
        os.system('casapy --nologger -c '+ SCRIPTPATH + '/casapy_cleanv4.py '+ ms + ' ' + imout + ' ' + mask_sources+'field' + \
                   ' ' + cleandepth2 + ' ' + str(niter) + ' ' + str(nterms) + ' ' + str(imsize) + ' ' + mscale)

    # convert to FITS
    if nterms > 1:
        os.system('image2fits in=' + imout +'.image.tt0' + ' ' + 'out='+ imout + '.fits')
    else:
        os.system('image2fits in=' + imout +'.image'     + ' ' + 'out='+ imout + '.fits')

    return imout, mask_sources+'field', imsize

def make_image_wsclean(mslist, cluster, callnumber, threshpix, threshisl, nterms, atrous_do, imsize, inputmask, mscale, region,cellsize,uvrange,wsclean,WSCleanRobust,BlankField):

    if imsize is None:
        imsize = image_size_from_mask(inputmask)

    niter   = numpy.int(5000 * (numpy.sqrt(numpy.float(len(mslist)))))
    cellsizeim = str(cellsize) +'arcsec'
    #wsclean =  '/home/rvweeren/software/WSClean/wsclean-1.6+MORESANE+MASKS4/build/wsclean'
    #wsclean = '/home/rvweeren/software/WSClean/wsclean-1.7/build/wsclean'

    depth =  1e-3*0.7 / (numpy.sqrt(numpy.float(len(mslist))))
    cleandepth1 = str(depth*1.5) #+ 'mJy'
    cleandepth2 = str(depth)     #+ 'mJy'

    wideband = False
    if len(mslist) > 5:
        wideband = True

    # speed up the imaging if possible by reducing image size within the mask region
    newsize = find_newsize(inputmask)
    if newsize < imsize: # ok so we can use a smaller image size then
        #make a new template
        os.system('casapy --nologger -c ' + SCRIPTPATH + '/make_empty_image.py '+ str(mslist[0]) + ' ' + inputmask+'2' + ' ' + str(newsize) + ' ' +'1.5arcsec')
        os.system('casapy --nologger -c ' + SCRIPTPATH + '/regrid_image.py '    + inputmask      + ' ' + inputmask+'2' + ' ' + inputmask+'3')

        # reset the imsize and the mask
        imsize    = newsize
        inputmask = inputmask+'3'


    ms = ''
    for m in mslist:
        ms = ms + ' ' + m

    imout = 'im'+ callnumber +'_cluster'+cluster+'nm'

    os.system('rm -rf ' + imout + '-*')

    # NDPPP concat
    outms      = 'field.ms'
    parsetname = 'concatforwsclean.parset'

    msinstr = ""

    for ms_id, ms in enumerate(mslist):
        msinstr = msinstr + "'" + ms + "'"
        if ms_id < len(mslist)-1:
            msinstr = msinstr + ", "
    os.system('rm -rf ' + parsetname)
    f=open(parsetname, 'w')
    f.write('msin = [%s]\n' % msinstr)
    f.write('msin.datacolumn = DATA\n')
    f.write('msin.missingdata=True\n')
    f.write('msin.orderms=False\n')
    f.write('msout=%s\n' % outms)
    f.write('steps=[]\n')
    f.close()
    os.system('rm -rf ' + outms)
    os.system('NDPPP ' + parsetname)

    if wideband:
        channelsout =  1 # there is a factor of 5 averaging
        cmd1 = wsclean + ' -reorder -name ' + imout + ' -size ' + str(imsize) + ' ' + str(imsize) + ' '
        cmd2 = '-scale ' + cellsizeim + ' -weight briggs '+str(WSCleanRobust)+' -niter ' + str(niter) + '-cleanborder 0 -threshold '+ cleandepth1 + ' '
        cmd3 = '-minuv-l '+ str(uvrange) \
               +' -mgain 0.75 -fitbeam -datacolumn DATA -no-update-model-required -joinchannels -channelsout ' +\
               str(channelsout) + ' '  + outms
    else:
        cmd1 = wsclean + ' -reorder -name ' + imout + ' -size ' + str(imsize) + ' ' + str(imsize) + ' '
        cmd2 = '-scale ' + cellsizeim + ' -weight briggs '+str(WSCleanRobust)+' -niter ' + str(niter) + ' -cleanborder 0 -threshold '+ cleandepth1 + ' '
        cmd3 = '-minuv-l '+ str(uvrange) +' -mgain 0.75 -fitbeam -datacolumn DATA -no-update-model-required ' + outms

    print cmd1+cmd2+cmd3
    os.system(cmd1+cmd2+cmd3)

    if BlankField:
        mask_image=blank.blank_facet(imout+'-image.fits',inputmask)
    else:
        mask_image=imout+'-image.fits'

    # create the mask
    os.system('python ' + SCRIPTPATH + '/makecleanmask_field_wsclean.py --threshpix '+str(threshpix)+\
              ' --threshisl '+str(threshisl) +' --atrous_do '+ str(atrous_do) + \
              ' --casaregion  '+ region + ' '  + mask_image)

    mask_name  = mask_image + '.fitsmask'
    casa_mask  = imout + '.casamask'

    maskim=pyrap.images.image(mask_name)
    maskim.saveas(casa_mask)

    # Convert to casapy format and includ region file
    if region != 'empty':
        os.system('casapy --nologger -c ' + SCRIPTPATH+'/fitsandregion2image.py '\
                  + mask_name + ' ' + casa_mask + ' ' + region)
    else:
        os.system('casapy --nologger -c ' + SCRIPTPATH+'/fitsandregion2image.py '\
                  + mask_name + ' ' + casa_mask + ' ' + 'None')

    mask_sources = imout+'.casamask'
    os.system('rm -rf ' + mask_sources + 'field')
    os.system('cp -r ' + mask_sources + ' ' + mask_sources + 'field')

    # Merge the two masks and blank outside template field
    img    = pyrap.images.image(mask_sources+'field')
    pixels = numpy.copy(img.getdata())

    img2    = pyrap.images.image(inputmask)
    pixels2 = numpy.copy(img2.getdata())

    idx = numpy.where(pixels2 == 0.0)
    pixels[idx] = 0.0
    img.putdata(pixels)
    img.unlock()
    del img
    del img2


    imout = 'im'+ callnumber +'_cluster'+cluster
    os.system('rm -rf ' + imout + '-*')
    niter = niter*5 # increase niter, tune manually if needed, try to reach threshold

    if wideband:
        cmd1 = wsclean + ' -reorder -name ' + imout + ' -size ' + str(imsize) + ' ' + str(imsize) + ' '
        cmd2 = '-scale ' + cellsizeim + ' -weight briggs '+str(WSCleanRobust)+' -niter ' + str(niter) + ' -cleanborder 0 -threshold '+ cleandepth2 + ' '
        cmd3 = '-minuv-l '+ str(uvrange) +' -mgain 0.6 -fitbeam -datacolumn DATA -no-update-model-required -casamask ' + \
               mask_sources+'field' + ' -joinchannels -channelsout ' + str(channelsout) + ' ' + outms
    else:
        cmd1 = wsclean + ' -reorder -name ' + imout + ' -size ' + str(imsize) + ' ' + str(imsize) + ' '
        cmd2 = '-scale ' + cellsizeim + ' -weight briggs '+str(WSCleanRobust)+' -niter ' + str(niter) + ' -cleanborder 0 -threshold '+ cleandepth2 + ' '
        cmd3 = '-minuv-l '+ str(uvrange) +' -mgain 0.6 -fitbeam -datacolumn DATA -no-update-model-required -casamask ' + \
               mask_sources+'field' + ' '+ outms

    print cmd1+cmd2+cmd3
    os.system(cmd1+cmd2+cmd3)

    # convert from FITS to casapy format
    # os.system('casapy --nologger -c ' + SCRIPTPATH +'/fits2image.py ' + \
    #           imout + '-image.fits' + ' ' + imout +'.image')
    finalim=pyrap.images.image(imout+'-image.fits')
    finalim.saveas(imout +'.image')

    return imout, mask_sources+'field', imsize


def do_fieldFFT(ms,image,imsize,cellsize,wsclean,mslist,WSCleanRobust):
    niter   = 1
    cellsizeim = str(cellsize)+ 'arcsec'

    # note no uvrange here!
    # also no re-order

    wideband = False
    if len(mslist) > 5:
        wideband = True
    if wideband:
        channelsout =  5 # DO NOT CHANGE !! (in make_image_wsclean_wideband there is averaging)
        cmd1 = wsclean + ' -predict -name ' + image + ' -size ' + str(imsize) + ' ' + str(imsize) + ' '
        cmd2 = '-scale ' + cellsizeim + ' -weight briggs '+str(WSCleanRobust)+' -niter ' + str(niter) + ' '
        cmd3 = '-cleanborder 0 -mgain 0.85 -fitbeam -datacolumn DATA '+ '-joinchannels -channelsout ' + str(channelsout) + ' ' + ms
    else:
        cmd1 = wsclean + ' -predict -name ' + image + ' -size ' + str(imsize) + ' ' + str(imsize) + ' '
        cmd2 = '-scale ' + cellsizeim + ' -weight briggs '+str(WSCleanRobust)+' -niter ' + str(niter) + ' '
        cmd3 = '-cleanborder 0 -mgain 0.85 -fitbeam -datacolumn DATA '+ ' ' + ms

    print cmd1+cmd2+cmd3
    os.system(cmd1+cmd2+cmd3)
    return

def image_size_from_mask(mask):
    im = pyrap.images.image(mask)
    sh = im.shape()
    if sh[-1] != sh[-2]:
        print "image is not square!"
        print sh[-1], sh[-2]
    npix = sh[-1]
    return npix


### END of FUNCTION DEFS, MAIN SCRIPT STARTS HERE#


source_info_rec = numpy.genfromtxt(peelsourceinfo, \
                                   dtype="S50,S25,S5,S5,i8,i8,i8,i8,S2,S255,S255,S255,S5", \
                                   names=["sourcelist","directions","atrous_do","mscale_field","imsizes",\
                                   "cellsizetime_p","cellsizetime_a","fieldsize","dynamicrange",\
                                   "regionselfc","regionfield","peelskymodel","outliersource"],comments='#')

sourcelist = source_info_rec["sourcelist"]
directions = source_info_rec["directions"]
atrous_do = source_info_rec["atrous_do"]
mscale_field = source_info_rec["mscale_field"]
imsizes = source_info_rec["imsizes"]
cellsizetime_p = source_info_rec["cellsizetime_p"]
cellsizetime_a = source_info_rec["cellsizetime_a"]
fieldsize = source_info_rec["fieldsize"]
dynamicrange = source_info_rec["dynamicrange"]
regionselfc = source_info_rec["regionselfc"]
regionfield = source_info_rec["regionfield"]
peelskymodel = source_info_rec["peelskymodel"]
outliersource = source_info_rec["outliersource"]

sourcelist = sourcelist.tolist()


# order is important (from low to high freq)
# frequencies should be contigeous, so also add names of bands that do not exist

#mslistorig = ['A2256_SB200-209.2ch10s.ms',\
#             'A2256_SB210-219.2ch10s.ms', 'A2256_SB220-229.2ch10s.ms', 'A2256_SB230-239.2ch10s.ms']



#mslistorig = ['A2256_SB050-059.2ch10s.ms',\
#             'A2256_SB060-069.2ch10s.ms', 'A2256_SB070-079.2ch10s.ms', 'A2256_SB080-089.2ch10s.ms',\
#             'A2256_SB090-099.2ch10s.ms', 'A2256_SB100-109.2ch10s.ms', 'A2256_SB110-119.2ch10s.ms',\
#             'A2256_SB120-129.2ch10s.ms', 'A2256_SB130-139.2ch10s.ms', 'A2256_SB140-149.2ch10s.ms',\
#             'A2256_SB150-159.2ch10s.ms', 'A2256_SB160-169.2ch10s.ms', 'A2256_SB170-179.2ch10s.ms',\
#             'A2256_SB180-189.2ch10s.ms', 'A2256_SB190-199.2ch10s.ms', 'A2256_SB200-209.2ch10s.ms',\
#             'A2256_SB210-219.2ch10s.ms', 'A2256_SB220-229.2ch10s.ms', 'A2256_SB230-239.2ch10s.ms',\
#             'A2256_SB240-249.2ch10s.ms', 'A2256_SB250-259.2ch10s.ms', 'A2256_SB260-269.2ch10s.ms',\
#             'A2256_SB270-279.2ch10s.ms', 'A2256_SB280-289.2ch10s.ms', 'A2256_SB290-299.2ch10s.ms',\
#             'A2256_SB300-309.2ch10s.ms',\
#             'A2256_SB310-319.2ch10s.ms', 'A2256_SB320-329.2ch10s.ms', 'A2256_SB330-339.2ch10s.ms',\
#             'A2256_SB340-349.2ch10s.ms',\
#             'A2256_SB350-359.2ch10s.ms']


mslistorig = ["{name:s}_SB{b1:03d}-{b2:03d}.{res:s}.ms".format(name=NAME,res=RES,b1=b,b2=b+9) for b in BANDS]
mslistorigstr = ' '.join(mslistorig)

# filter out datasets that do not exist (takes also care of freq.gaps in field subtract))
mslist= [ms for ms in mslistorig if  os.path.isdir(ms)]
msliststr = ' '.join(mslist)


#mslistorigstr = ''
#for ms in mslistorig:
#  mslistorigstr =  mslistorigstr + ' ' + ms

#mslist= []
#for ms in mslistorig:
#  if os.path.isdir(ms): # filter out datasets that do not exist (takes also care of freq.gaps in field subtract))
#    mslist.append(ms)

msliststr = ''
for ms in mslist:
    msliststr = msliststr + ' ' + ms



tt = pt.table(mslist[0] + '/FIELD')
pointingcenter = tt.getcol('REFERENCE_DIR')[0][0]
pointingcenter = str(pointingcenter[0]) +'rad,' + str(pointingcenter[1])+'rad'
print pointingcenter
tt.close()

#if len(mslist) == 1:
#  TEC    = "False" # no TEC fitting for one (channel) dataset
#  nterms = 1

if len(mslist) > 8:
    nterms = 2
if len(mslist) > 40:
    nterms = 3
if len(mslist) < 20:
    clock = "False"



numpy.save('directions.npy', directions)


### MAKE ALL THE MASKS AT ONCE ###
if makemasks:
    for source_id,source in enumerate(sourcelist):
        msavglist = []
        for ms_id, ms in enumerate(mslist):
            msavglist.append(ms.split('.')[0] + '.' + source + '.ms')

        tmpn  = str(msavglist[0])
        parset = create_phaseshift_parset_formasks(mslist[0], tmpn, source, directions[source_id])
        os.system('NDPPP ' + parset)
        output_template_im = 'templatemask_' + source
        print output_template_im
        os.system('casapy --nologger -c ' + SCRIPTPATH + '/make_empty_image.py '+ tmpn + ' ' + output_template_im + ' ' + str(fieldsize[source_id]) + ' ' +'1.5arcsec')
        os.system('rm -rf ' + tmpn)
        # now generate the mask
        os.system(SCRIPTPATH + '/make_facet_mask.py ' + output_template_im +' ' + 'directions.npy' + ' ' + str(source_id) + ' ' + '1.5arcsec'  +' ' + '&')
    sys.exit()
##################################

logging.info('\n')
logging.info('#######################################################################\n')
logging.info('Running DDE')
logging.info('Doing sources: '+','.join(do_sources))
logging.info('Using MSlist: '+msliststr)
logging.info('Uvrange: '+ str(uvrange) + ' lambda')
logging.info('Cellsize: '+str(cellsize)+ ' arcsec')


for source in do_sources:

    source_id = sourcelist.index(source)

    print 'DOING DDE patch:', source
    logging.info('')
    logging.info('DOING DDE patch: '+ source)

# tidying-up code from Wendy's version

    if StartAtStep in ['preSC']:
            # remove selfcal images #
        logging.info("removing any existing sc ms")
        os.system("rm -rf *."+source+".ms")
    if StartAtStep in ['preSC', 'doSC']:
        # remove selfcal images #
        logging.info("removing any existing selfcal images")
        os.system("rm -rf im*_cluster"+source+"*")
        #os.system("rm -rf allbands.concat.shifted_'+source+'.ms")
    if StartAtStep in ['preSC', 'doSC', 'postSC','preFACET','doFACET']:
        # remove selfcal images #
        logging.info("removing any existing facet images")
        os.system("rm -rf imfield*_cluster"+source+"*")
    if StartAtStep in ['preSC', 'doSC', 'postSC','preFACET']:
        logging.info("removing any existing facet imaging average MS")
        os.system("rm -rf *."+source+".ms.avgfield*")

    #check if allbands.concat.shifted.ms is present
    if os.path.isdir('allbands.concat.shifted.ms'):
        print 'allbands.concat.shifted.ms already exists'
        if StartAtStep in ['preSC']:
            #raise Exception('delete measurement set and then restart')
            os.system('rm -rf allbands.concat.shifted_'+source+'.ms')
            logging.info('removing')
        else:
            logging.info('...but continuing because we are redoing selfcal')

    if not os.path.isdir('allbands.concat.ms'):
        print 'allbands.concat.ms does not exist'
        raise Exception('make measurement set and then restart')

    dummyskymodel   = SCRIPTPATH + '/dummy.skymodel' ## update every time again with new source, not used, just a dummy for correct

    msavglist = []
    for ms_id, ms in enumerate(mslist):
        msavglist.append(ms.split('.')[0] + '.' + source + '.ms')

    # set image region
    tmpn  = str(msavglist[0])

    output_template_im = 'templatemask_' + source
    if not os.path.exists(output_template_im+'.masktmp'):
        raise  Exception('facet mask missing: '+output_template_im+'.masktmp')

    selfcaldir = directions[source_id]
    facetsize = image_size_from_mask(output_template_im+'.masktmp')

    logging.info("Selfcal direction: "+selfcaldir)
    logging.info("facetmask: "+str(facetsize))


    ## STEP 1: prep for SC ##
    if StartAtStep in ['preSC']:
## FIXME -- hard-wired CPU limit in what follows
        if len(mslist) > 22:
            runbbs_diffskymodel_addback16(mslist, 'instrument_ap_smoothed', True, directions[source_id],imsizes[source_id],output_template_im, do_ap)
        else:
            runbbs_diffskymodel_addback(mslist, 'instrument_ap_smoothed', True, directions[source_id],imsizes[source_id],output_template_im, do_ap)

        # average and phaseshift with NDPPP
        for ms_id, ms in enumerate(mslist):
            parset = create_phaseshift_parset(ms, msavglist[ms_id], source, directions[source_id],imsizes[source_id], dynamicrange[source_id])
            os.system('rm -rf ' + msavglist[ms_id])
            os.system('NDPPP ' + parset)


    ### PHASESHIFT the FULL resolution dataset, for MODEL_DATA FFT subtract
        if outliersource[source_id] == 'False':
            parset = create_phaseshift_parset_full('allbands.concat.ms', 'allbands.concat.shifted.ms', directions[source_id],'DATA')
            os.system('NDPPP ' + parset + '&') # run in background


    ## END STEP 1
    ## STEP 2a: SC ##
    if StartAtStep in ['preSC', 'doSC']:

        # correct with amps and phases from selfcal, needs to be done here because CORRECTED_DATA needs to be reset for the selfcal
        if do_ap:
            runbbs_2(msavglist, mslist, dummyskymodel ,SCRIPTPATH+'/correct.parset','instrument_ap_smoothed')

        logging.info('selfcal started '+source)
        # do the selfcal loop
        # make string ms list for input selfcal
        inputmslist = ''
        for ms in msavglist:
            inputmslist = inputmslist + ' ' + ms

        logging.info('Start selfcal DDE patch: '+ source)
        logging.info('Solint CommonScalarPhase: '+ str(cellsizetime_p[source_id]))
        logging.info('Solint A&P: '+ str(cellsizetime_a[source_id]))
        logging.info('Region file: '+ str(regionselfc[source_id]))

        #os.system('python ' + SCRIPTPATH + '/selfcalv19.py ' + inputmslist + ' ' + source + ' ' + atrous_do[source_id] + ' ' + str(imsizes[source_id]) + ' ' + \
        #             str(nterms) + ' ' + str(cellsizetime_a[source_id]) + ' ' + str(cellsizetime_p[source_id]) + ' ' + TEC + ' ' + clock + ' ' + \
        #         str(dynamicrange[source_id]) + ' ' + regionselfc[source_id] + ' ' + str(uvrange) + ' ' + str(peelskymodel[source_id]) + ' ' +\
        #         str(cellsize))


        os.system('python '+SCRIPTPATH+'/selfcalv19_ww_cep3.py ' + inputmslist + ' ' + source + ' ' + atrous_do[source_id] + ' ' + str(imsizes[source_id]) + ' ' + \
                         str(nterms) + ' ' + str(cellsizetime_a[source_id]) + ' ' + str(cellsizetime_p[source_id]) + ' ' + TEC + ' ' + clock + ' ' + \
                         str(dynamicrange[source_id]) + ' ' + regionselfc[source_id])



        logging.info('Finished selfcal DDE patch: '+ source)

    ## STEP 2b:  SC wrap up ##
    if StartAtStep in ['preSC', 'doSC', 'postSC']:
        # combine selfcal solutions with non-DDE phases
        for ms_id, ms in enumerate(mslist):
            parmdb_selfcal     = msavglist[ms_id]+"/"+"instrument_merged"
                #parmdb_nondde      = ms+"/"+"instrument_ap_smoothed_ampto1" # contains only phase corrections
            parmdb_nondde      = ms+"/"+"instrument_ap_smoothed" # contains only phase corrections
            this_parmdb_master_out  = ms+"/"+"instrument_master_" + source
            parmdb_template    = msavglist[ms_id]+"/"+"instrument_template"
            join_parmdb(ms, parmdb_selfcal,parmdb_nondde, parmdb_template, this_parmdb_master_out, TEC, clock)
            print 'joined SC and DDE parmdbs for {ms}'.format(ms=ms)
            parmdb_master_out  = "instrument_master_" + source   # reset because runbbs uses basename of ms

# stefcal code, doesn't work at the moment
#        parmdb_selfcal     = msavglist[ms_id]+"/"+"instrument_merged"
#        parmdb_master_out  = ms+"/"+"instrument_master_" + source
#        parmdb_template    = msavglist[ms_id]+"/"+"instrument_template"
#        join_parmdb(ms, parmdb_selfcal,parmdb_template, parmdb_master_out)
#

        # maybe there are some issues with the frequency boundaries if you solve on averaged data
        for ms_id, ms in enumerate(mslist):
            parmdb_master_outtmp  = ms+"/"+"instrument_master_" + source
            os.system("taql 'update " + parmdb_master_outtmp + " set ENDX=1.e12'")
            os.system("taql 'update " + parmdb_master_outtmp + " set STARTX=1.0'")

        if outliersource[source_id] == 'False':
            # normalize the solutions to 1.0, for outlier sources do not(!) normalize
            # this is a global normalization (one --single-- factor for all SB combined)
            nvalue = normalize_parmdbs(mslist,parmdb_master_out, parmdb_master_out+'_norm')
            logging.info('Normalized amps to 1.0, found mean amplitude value of ' + str(nvalue))

        # plot the solutions
        for ms in mslist:
            parmdb_master_plot  = ms+"/"+"instrument_master_" + source
            plotim_base         = ms.split('.')[0] + "_instrument_master_" + source
            if TEC=='True':
                os.system(SCRIPTPATH + '/plot_solutions_all_stations_v2.py \
                        -t -a -p --freq 150 ' + parmdb_master_plot + ' ' + plotim_base +'&')
            else:
                os.system(SCRIPTPATH + '/plot_solutions_all_stations_v2.py \
                        -s -a -p --freq 150 ' + parmdb_master_plot + ' ' + plotim_base +'&')

        print 'Updated frequency boundaries parmdb and normalized amps to 1.0'
        time.sleep(5)
        # make new mslist for field averaged data
        msavglist = []
        for ms_id, ms in enumerate(mslist):
            msavglist.append(ms.split('.')[0] + '.' + source + '.ms.avgfield')


        ######### check if all plot_solutions_all_stations_v2.py processes is finished
        cmd = "ps -f -u " + username + " | grep plot_solutions_all_stations_v2.py | grep -v grep |wc -l"
        output=numpy.int(Popen(cmd, shell=True, stdout=PIPE).communicate()[0])
        while output > 0 : # "grep -v grep" to prevent counting the grep command
            time.sleep(2)
            output=numpy.int(Popen(cmd, shell=True, stdout=PIPE).communicate()[0])
        #########


    parmdb_master_out="instrument_master_" + source
    if outliersource[source_id] == 'False':
        if StartAtStep in ['preSC', 'doSC', 'postSC','preFACET']:
            ## STEP 3: prep for facet ##
            parmdb_master_out="instrument_master_" + source
            runbbs_diffskymodel_addbackfield(mslist, 'instrument_ap_smoothed', True,  directions[source_id],imsizes[source_id], output_template_im, do_ap)
            logging.info('Adding back rest of the field for DDE facet ' + source)

            if TEC=='True':
                runbbs(mslist, dummyskymodel, SCRIPTPATH + '/correctfield2+TEC.parset',parmdb_master_out+'_norm', False)
            else:
                runbbs(mslist, dummyskymodel, SCRIPTPATH + '/correctfield2.parset',parmdb_master_out+'_norm', False)
            ###########################################################################
            # NDPPP phase shift, less averaging (NEW: run 2 in parallel)
            msavglist = []
            for ms_id, ms in enumerate(mslistorig): # remake msavglist from mslistorig
                msavglist.append(ms.split('.')[0] + '.' + source + '.ms.avgfield')

            for ms_id, ms in enumerate(mslist):
                parset = create_phaseshift_parset_field(ms, msavglist[ms_id], source, directions[source_id])

                cmd = "ps -u " + username + " | grep NDPPP | wc -l"
                output=numpy.int(Popen(cmd, shell=True, stdout=PIPE).communicate()[0])

                while output > 1 : # max 2 processes (max 2, instead of 3, because we also phaseshift)
                    time.sleep(10)
                    output=numpy.int(Popen(cmd, shell=True, stdout=PIPE).communicate()[0])
                    pid = (Popen('pidof NDPPP', shell=True, stdout=PIPE).communicate()[0])
                    pid_list = pid.split(' ')
                # START NDPPP BECAUSE LESS/EQ 2 PROCESSES ARE RUNNING
                os.system('rm -rf ' + msavglist[ms_id])
                os.system('NDPPP ' + parset + '&')

            # Check if all NDPPP processes are finished
            output=numpy.int(Popen(cmd, shell=True, stdout=PIPE).communicate()[0])
            while output > 0 :
                time.sleep(10)
                output=numpy.int(Popen(cmd, shell=True, stdout=PIPE).communicate()[0])
                pid = (Popen('pidof NDPPP', shell=True, stdout=PIPE).communicate()[0])
                pid_list = pid.split(' ')
            ###########################################################################
        ## STEP 4a -- do facet ##

        ###### MAKE FACET IMAGE #####

        # imsize None forces the code to work out the image size from the mask size
        if StartAtStep in ['preSC', 'doSC', 'postSC','preFACET','doFACET']:
            msavglist = []
            for ms_id, ms in enumerate(mslistorig): # remake msavglist from mslistorig
                msavglist.append(ms.split('.')[0] + '.' + source + '.ms.avgfield')

            imout,mask_out, imsizef = make_image_wsclean(msavglist, source, 'field0', 5, 3, nterms, 'True', None, output_template_im +'.masktmp', mscale_field[source_id],regionfield[source_id],cellsize, uvrange,wsclean,WSCleanRobust,BlankField)
            logging.info('Imaged full DDE facet: ' + source)

        ## STEP 4b -- post facet ##

        if StartAtStep in ['preSC', 'doSC', 'postSC','preFACET','doFACET','postFACET']:

# if we are restarting, it's possible that allbands.concat.shifted.ms may have been deleted earlier. So re-create it if it doesn't exist

            if not(os.path.isdir('allbands.concat.shifted.ms')):
                print 'allbands.concat.shifted.ms does not exist, re-creating'
                parset = create_phaseshift_parset_full('allbands.concat.ms', 'allbands.concat.shifted.ms', directions[source_id],'DATA')
                os.system('NDPPP ' + parset)
            if StartAtStep=='postFACET':
                # imout won't be set, so guess it
                imout='imfield0_cluster'+source
                imsizef=image_size_from_mask(output_template_im +'.masktmp')

            # BACKUP SUBTRACTED DATA IN CASE OF CRASH
            ###########################################################################
            # (NEW: run in parallel)
            for ms_id, ms in enumerate(mslist):
                cmd = "ps -u " + username + " | grep taql | wc -l"
                output=numpy.int(Popen(cmd, shell=True, stdout=PIPE).communicate()[0])

                while output > 3 : # max 4 processes
                    time.sleep(10)
                    output=numpy.int(Popen(cmd, shell=True, stdout=PIPE).communicate()[0])

                # START taql BECAUSE LESS/EQ 2 PROCESSES ARE RUNNING
                os.system("taql 'update " + ms + " set CORRECTED_DATA=SUBTRACTED_DATA_ALL' &")

            # Check if all taql processes are finished
            output=numpy.int(Popen(cmd, shell=True, stdout=PIPE).communicate()[0])
            while output > 0 :
                time.sleep(10)
                output=numpy.int(Popen(cmd, shell=True, stdout=PIPE).communicate()[0])
            logging.info('Backup SUBTRACTED_DATA_ALL')
            ###########################################################################


            # DO THE FFT
            do_fieldFFT('allbands.concat.shifted.ms',imout, imsizef, cellsize, wsclean, msavglist, WSCleanRobust)
            logging.info('FFTed model of DDE facet: ' + source)

            # SHIFT PHASE CENTER BACK TO ORIGINAL
            logging.info('Shift model back to pointing centre')
            parset = create_phaseshift_parset_full('allbands.concat.shifted.ms', 'allbands.concat.shiftedback.ms',\
                                          pointingcenter,'MODEL_DATA')

            os.system('NDPPP ' + parset)
            os.system('rm -rf allbands.concat.shifted.ms') # clean up

            # Add MODEL_DATA (allbands.concat.shiftedback.ms) into ADDED_DATA_SOURCE from mslist

            freq_tab1= pt.table('allbands.concat.ms' + '/SPECTRAL_WINDOW')
            numchan1    = freq_tab1.getcol('NUM_CHAN')
            freq_tab2= pt.table(mslist[0] + '/SPECTRAL_WINDOW')
            numchan2    = freq_tab2.getcol('NUM_CHAN')
            freq_tab1.close()
            freq_tab2.close()

            if (numchan1[0]) == (numchan2[0]*len(mslist)):
                os.system('python ' + SCRIPTPATH + '/copy_over_columns.py '+ msliststr +\
                          ' ' +'allbands.concat.shiftedback.ms'+' ' + 'ADDED_DATA_SOURCE')
            else:
                os.system('python ' + SCRIPTPATH + '/copy_over_columns.py '+ mslistorigstr +\
                          ' ' +'allbands.concat.shiftedback.ms'+' ' + 'ADDED_DATA_SOURCE')

            os.system('rm -rf allbands.concat.shiftedback.ms') # clean up

    #### OUTLIER CASE ####
    else:  # do this because we are not going to add back field sources
        logging.info('Do not add field back for outlier source')
        for ms in mslist:
            os.system("taql 'update " + ms + " set MODEL_DATA=ADDED_DATA_SOURCE'")

        # BACKUP SUBTRACTED DATA IN CASE OF CRASH
        for ms in mslist:
            os.system("taql 'update " + ms + " set CORRECTED_DATA=SUBTRACTED_DATA_ALL'")
        logging.info('Backup SUBTRACTED_DATA_ALL for outliersource')

    if StartAtStep in ['preSC', 'doSC', 'postSC','preFACET','doFACET','postFACET']:
    #### DO THE SUBTRACT ####
        if peelskymodel[source_id] != 'empty': # should also cover "outliersource"
            print 'Subtracting source with a user defined skymodel', peelskymodel[source_id]
            parset   = create_subtract_parset_field_outlier('SUBTRACTED_DATA_ALL',TEC)
            runbbs(mslist, peelskymodel[source_id], parset, parmdb_master_out, True) # NOTE: no 'normalization' and replace sourcedb
            logging.info('Subtracted outlier source from data for DDE : ' + source)
        else:
            parset   = create_subtract_parset_field('SUBTRACTED_DATA_ALL',TEC)
            runbbs(mslist, ' ', parset, parmdb_master_out+'_norm', False) # replace-sourcedb not needed since we use "@column"
            logging.info('Subtracted facet model from data for DDE : ' + source)

        # CHECK IF THE SUBTRACT WORKED OK by making dirty low-res images
        inputmslist = ''
        for ms in mslist:
            inputmslist = inputmslist + ' ' + ms
        #os.system('python ' + SCRIPTPATH + '/verify_subtract_v3.py ' + inputmslist + ' 0.3 ' + source)
        os.system('python '+ SCRIPTPATH+'/verify_subtract_v5.py ' + inputmslist + ' 0.15 ' + source)

    logging.info('finished '+source)
