import matplotlib
#matplotlib.use('GTK')
import numpy
import os
import sys
import lofar.parmdb
from scipy import interpolate
import time
#import subprocess
from subprocess import Popen, PIPE, STDOUT
import pyrap.tables as pt
import pyrap.images
import pwd
import logging
import blank
from coordinates_mode import *
pi = numpy.pi




# check high-DR
# check old selfcalv15

# check point source model


# TO DO
# - problem with e1 mask (sources not in clean mask)
# - problem with s3, diffuse source not in clean mask causes neg. bowl --> FIXED
# - check for mask gaps --> DONE
# - fit for clock in selfcal --> DONE
# - clean component overlap check -> avoid by giving up patches
# - reset phases to zero in selfcal part --> DONE
# - gain normalization --> DONE
# - HIGH-DYNAMIC RANGE (need to adjust merger/join parmdb, solution smoohting as well)
# - move m10 direction --> DONE
# - make normal template names files --> DONE
# - remove source list loading to prevent segfaults --> DONE

# - run the addback sources in seperate python script to avoid frequent Bus Errors DONE
# - run the addback sources FIELD in seperate python script to avoid frequent Bus Errors DONE
# - make sure parrallel FFT does not run out of memory (limit max. processes)/facets? DONE
# - automatic w-planes computation in imaging
# - how to deal with outlier sources low-freqs? # peel them first using low-freq bands only?
# - combine images into big image --> DONE
# - make fake avg pb images from mask to work with msss mosaic --> DONE
# - convolve images all to the same resolution (use casapy for that, easy to do)
# - What about second imaging and calibration cycle?
#    - 1. ADD back skymodel using "master solutions", then CORRECT using "master solutions"
#    - 2. image that again using the same setting
#    - 3. redo the subtract (will be slightly better....but solutions remained the same, just better noise) or just proceed to the next field?
# - how to deal with Toothbrush, how to subtract it from the data properly (make low-res images for that?)
# - memory control in full resolution subtract step --> DONE


# v6
# - added SCRIPTPATH - makes it easier to switch between cep and para
# - added check of return values on some system calls (esp NDPPP which returns 0 on success)

# v7
# increase niter (for deeper clean)
# change scheme 2 to 3
#2. your new scheme
#- correct with solutions
#- phaseshift to updated position +avg
#- facet imaging

#3. proposed scheme
#- correct with solutions
#- phaseshift to selfcal center + avg
#- phaseshift to updated position (from selfcal center)
#- facet imaging

# v8
# go back... to scheme 2


# v9
# for difficult to low snr facet, but not so far from another facet try just using the other facet solutions

def ra_to_str(dra, ndec=2,delim=':'):
    '''
    converts a single decimal degrees ra to hh:mm:ss.s
    '''
    if delim == 'h':
        delim1 = 'h'
        delim2 = 'm'
    else:
        delim1 = delim
        delim2 = delim

    dra = dra/15.
    dd = math.floor(dra)
    dfrac = dra - dd
    dmins = dfrac*60.
    dm = math.floor(dmins)
    dsec = (dmins-dm)*60.
    if round(dsec, ndec) == 60.00:
        dsec = 0.
        dm += 1
    if dm == 60.:
        dm = 0.
        dd += 1
    sra = '%02d%s%02d%s%05.2f' %(dd,delim1,dm,delim2,dsec)
    return sra
def dec_to_str(ddec,ndec=1,delim=':'):
    '''
    converts a single decimal degrees dec to dd:mm:ss.s
    '''
    if delim == 'd':
        delim1 = 'd'
        delim2 = 'm'
    else:
        delim1 = delim
        delim2 = delim

    dd = math.floor(ddec)
    dfrac = ddec - dd
    dmins = dfrac*60.
    dm = math.floor(dmins)
    dsec = (dmins-dm)*60.
    if round(dsec, ndec) == 60.0:
        dsec = 0.
        dm += 1
    if dm == 60.:
        dm = 0.
        dd += 1
    sdec = '%02d%s%02d%s%04.1f' %(dd,delim1,dm,delim2,dsec)
    return sdec


def find_newsize(mask):

    img    = pyrap.images.image(mask)
    pixels = numpy.copy(img.getdata())
    sh     = numpy.shape(pixels)[3:4]
    newsize = numpy.copy(sh[0])
    sh      = sh[0]
    #print mask, newsize

    trysizes = numpy.copy(sorted([6400,6144,5600,5400,5184,5000,4800,4608,4320,4096,3840,3600,3200,3072,2880,2560,2304,2048, 1600, 1536, 1200, 1024, 800, 512]))

    trysizes1 = numpy.copy(trysizes[::-1]) # reverse sorted

    idx = numpy.where(pixels  != 0)
    idx_x = idx[2]
    idx_y = idx[3]
    span_x1 = 2*abs(max(idx_x) - sh/2)
    span_x2 = 2*abs(sh/2 - min(idx_x))
    span_y1 = 2*abs(max(idx_y)- sh/2)
    span_y2 = 2*abs(sh/2 - min(idx_y))
    max_span = max(span_x1, span_y1, span_x2, span_y2)
    this_size = trysizes1[numpy.sum(trysizes1>= max_span)-1]
    newsize = this_size

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

def create_subtract_parset_field_outlier(outputcolumn, TEC, clock):
    bbs_parset = 'sub.parset'
    os.system('rm -f ' + bbs_parset)
    f=open(bbs_parset, 'w')

    chunksize = 100

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
    if clock == "True":
        f.write('Step.subtract.Model.Clock.Enable  = T\n')
    f.write('Step.subtract.Model.CommonRotation.Enable     = F\n')
    f.write('Step.subtract.Operation                       = SUBTRACT\n')
    f.write('Step.subtract.Model.Beam.Enable               = F\n')
    f.write('Step.subtract.Output.WriteCovariance          = F\n')
    f.write('Step.subtract.Output.Column                   = %s\n' % outputcolumn)
    f.close()
    return bbs_parset





def runbbs_diffskymodel_addback(mslist, parmdb, replacesource, direction, imsize, output_template_im):

    for ms in mslist:
        log      =  ms + '.bbslog'

        #set skymodel # ~weeren does not work in numpy.load
        skymodel = ms.split('.')[0] + '.skymodel'

        # find sources to add back, make parset

        #callist, callistarraysources = cal_return_slist(output_template_im +'.masktmp',skymodel, direction, imsize)
        cmd = 'python '+SCRIPTPATH+'/cal_return_slist.py '+ output_template_im+' '+skymodel +' "'+str(direction) +'" ' + str(imsize)
        output = Popen(cmd, shell=True, stdout=PIPE).communicate()[0]
        callist = output.strip()
        callistarraysources = callist.split(',')

        print 'Adding back for calibration:', callist


        parset = create_add_parset_ms(callist, ms)


        if replacesource:
            cmd = 'calibrate-stand-alone --replace-sourcedb --parmdb-name ' + parmdb + ' ' + ms + ' ' + parset + ' ' + skymodel + '>' + log + ' 2>&1 &'
        else:
            cmd = 'calibrate-stand-alone --parmdb-name ' + parmdb + ' ' + ms + ' ' + parset + ' ' + skymodel + '>' + log + ' 2>&1 &'
        print cmd
        os.system(cmd)
        time.sleep(10)  # otherwise add.parset is deleted (takes time for BBS to start up)

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

def runbbs_diffskymodel_addback16(mslist, parmdb, replacesource, direction, imsize, output_template_im):
    mslist1 = mslist[0:16]
    mslist2 = mslist[16:len(mslist)]
    #part1
    for ms in mslist1:
        log      =  ms + '.bbslog'

        #set skymodel # ~weeren does not work in numpy.load
        skymodel = ms.split('.')[0] + '.skymodel'

        # find sources to add back, make parset
        #callist, callistarraysources = cal_return_slist(output_template_im +'.masktmp',skymodel, direction, imsize)
        cmd = 'python '+SCRIPTPATH+'/cal_return_slist.py '+ output_template_im+' '+skymodel +' "'+str(direction) +'" ' + str(imsize)
        output = Popen(cmd, shell=True, stdout=PIPE).communicate()[0]
        callist = output.strip()
        callistarraysources = callist.split(',')


        print 'Adding back for calibration:', callist
        parset = create_add_parset_ms(callist, ms)

        if replacesource:
            cmd = 'calibrate-stand-alone --replace-sourcedb --parmdb-name ' + parmdb + ' ' + ms + ' ' + parset + ' ' + skymodel + '>' + log + ' 2>&1 &'
        else:
            cmd = 'calibrate-stand-alone --parmdb-name ' + parmdb + ' ' + ms + ' ' + parset + ' ' + skymodel + '>' + log + ' 2>&1 &'
        print cmd
        os.system(cmd)
        time.sleep(10)  # otherwise add.parset is deleted (takes time for BBS to start up)

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
        skymodel = ms.split('.')[0] + '.skymodel'

        # find sources to add back, make parset
        #callist, callistarraysources = cal_return_slist(output_template_im +'.masktmp',skymodel, direction, imsize)
        cmd = 'python '+SCRIPTPATH+'/cal_return_slist.py '+ output_template_im+' '+skymodel +' "'+str(direction) +'" ' + str(imsize)
        output = Popen(cmd, shell=True, stdout=PIPE).communicate()[0]
        callist = output.strip()
        callistarraysources = callist.split(',')

        print 'Adding back for calibration:', callist
        parset = create_add_parset_ms(callist, ms)

        if replacesource:
            cmd = 'calibrate-stand-alone --replace-sourcedb --parmdb-name ' + parmdb + ' ' + ms + ' ' + parset + ' ' + skymodel + '>' + log + ' 2>&1 &'
        else:
            cmd = 'calibrate-stand-alone --parmdb-name ' + parmdb + ' ' + ms + ' ' + parset + ' ' + skymodel + '>' + log + ' 2>&1 &'
        print cmd
        os.system(cmd)
        time.sleep(10)  # otherwise add.parset is deleted (takes time for BBS to start up)

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


def runbbs_diffskymodel_addbackfield(mslist, parmdb, replacesource, direction, imsize, output_template_im):

    for ms in mslist:
        log      =  ms + '.bbslog'

        #set skymodel
        skymodel = ms.split('.')[0] + '.skymodel'

        # find peeling sources (from previous step)
        #callist, callistarraysources = cal_return_slist(output_template_im +'.masktmp',skymodel, direction, imsize)
        cmd = 'python '+SCRIPTPATH+'/cal_return_slist.py '+ output_template_im+' '+skymodel +' "'+str(direction) +'" ' + str(imsize)
        output = Popen(cmd, shell=True, stdout=PIPE).communicate()[0]
        callist = output.strip()
        callistarraysources = callist.split(',')

        print 'Add field back step 1'

        # return the source list from the source to be added back sourrinding the peeling source and which fall within the mask boundaries
        # put in MODEL_DATA


        #addback_sourcelist = return_slist(output_template_im +'.masktmp', skymodel, callistarraysources)
        cmd = 'python '+SCRIPTPATH+'/return_slist.py '+ output_template_im+' '+skymodel +' "'+str(callist)+'"'
        output = Popen(cmd, shell=True, stdout=PIPE).communicate()[0]
        addback_sourcelist = output.strip()




        print 'Field source added back are: ', addback_sourcelist

        if len(addback_sourcelist) != 0: # otherwise do not have to add
            parset = create_add_parset_field_ms(addback_sourcelist, ms)
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


def create_phaseshift_parset(msin, msout, source, direction, freqstep=30):
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

    # just to test

    f.write('avg1.freqstep = %i\n' % freqstep)
    f.write('avg1.timestep = 1\n')
    f.close()
    return ndppp_parset


def create_phaseshift_parset_formasks(msin, msout, source, direction, freqstep=30):
    ndppp_parset = (msin.split('.')[0]) +'_ndppp_avgphaseshift.parset'
    os.system('rm -f ' + ndppp_parset)

    f=open(ndppp_parset, 'w')
    f.write('msin ="%s"\n' % msin)
    f.write('msin.datacolumn = DATA\n')
    f.write('msin.autoweight = false\n')
    f.write('msout ="%s"\n' % msout)
    f.write('msout.writefullresflag=False\n')
    f.write('steps = [shift,avg1]\n')
    f.write('shift.type        = phaseshift\n')
    f.write('shift.phasecenter = [%s]\n' % direction)
    f.write('avg1.type = squash\n')
    f.write('avg1.freqstep = %i\n' % freqstep)
    f.write('avg1.timestep = 1\n')
    f.close()
    return ndppp_parset


def create_phaseshift_parset_field(msin, msout, source, direction):
    ndppp_parset = msin.split('.')[0] +'ndppp_avgphaseshift_field2.parset'
    os.system('rm -f ' + ndppp_parset)

    f=open(ndppp_parset, 'w')
    f.write('msin ="%s"\n' % msin)
    f.write('msin.datacolumn = DATA\n')  # this is the second phaseshift
    f.write('msin.autoweight = false\n')
    f.write('msout ="%s"\n' % msout)
    f.write('msout.writefullresflag=False\n')
    f.write('steps = [shift]\n')
    f.write('shift.type        = phaseshift\n')
    f.write('shift.phasecenter = [%s]\n' % direction)
    f.close()
    return ndppp_parset

def create_phaseshift_parset_field_avg(msin, msout, source, direction):
    ndppp_parset = msin.split('.')[0] +'ndppp_avgphaseshift_field1.parset'
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



def create_add_parset_ms(source, ms):
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
    f.write('Step.add.Model.Gain.Enable               = T\n')
    f.write('Step.add.Model.Rotation.Enable           = F\n')
    f.write('Step.add.Model.CommonScalarPhase.Enable  = F\n')
    f.write('Step.add.Model.CommonRotation.Enable     = F\n')
    f.write('Step.add.Operation                       = ADD\n')
    f.write('Step.add.Model.Beam.Enable               = F\n')
    f.write('Step.add.Output.WriteCovariance          = F\n')
    f.write('Step.add.Output.Column                   = ADDED_DATA_SOURCE\n')
    f.close()
    return bbs_parset


#def create_add_parset_field(source):
    #bbs_parset = 'addfield.parset'
    #os.system('rm -f ' + bbs_parset)
    #f=open(bbs_parset, 'w')
    #f.write('Strategy.InputColumn = ADDED_DATA_SOURCE\n') # already contains peeling source
    #f.write('Strategy.ChunkSize   = 200\n')
    #f.write('Strategy.UseSolver   = F\n')
    #f.write('Strategy.Steps       = [add]\n\n\n')
    #f.write('Step.add.Model.Sources                   = [%s]\n' % source)
    #f.write('Step.add.Model.Cache.Enable              = T\n')
    #f.write('Step.add.Model.Phasors.Enable            = F\n')
    #f.write('Step.add.Model.DirectionalGain.Enable    = F\n')
    #f.write('Step.add.Model.Gain.Enable               = T\n')
    #f.write('Step.add.Model.Rotation.Enable           = F\n')
    #f.write('Step.add.Model.CommonScalarPhase.Enable  = F\n')
    #f.write('Step.add.Model.CommonRotation.Enable     = F\n')
    #f.write('Step.add.Operation                       = ADD\n')
    #f.write('Step.add.Model.Beam.Enable               = F\n')
    #f.write('Step.add.Output.WriteCovariance          = F\n')
    #f.write('Step.add.Output.Column                   = MODEL_DATA\n') # use use to save disk space
    #f.close()
    #return bbs_parset

def create_add_parset_field_ms(source, ms):
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

def create_subtract_parset_field(outputcolumn, TEC, clock):
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
    if clock == "True":
        f.write('Step.subtract.Model.Clock.Enable  = T\n')
    f.write('Step.subtract.Model.CommonRotation.Enable     = F\n')
    f.write('Step.subtract.Operation                       = SUBTRACT\n')
    f.write('Step.subtract.Model.Beam.Enable               = F\n')
    f.write('Step.subtract.Output.WriteCovariance          = F\n')
    f.write('Step.subtract.Output.Column                   = %s\n' % outputcolumn)
    f.close()
    return bbs_parset


def create_add_parset_field(outputcolumn, TEC, clock):
    bbs_parset = 'add.parset'
    os.system('rm -f ' + bbs_parset)
    f=open(bbs_parset, 'w')

    chunksize = 175

    f.write('Strategy.InputColumn = SUBTRACTED_DATA_ALL\n')
    f.write('Strategy.ChunkSize   = %s\n' % chunksize)
    f.write('Strategy.UseSolver   = F\n')
    f.write('Strategy.Steps       = [add]\n\n\n')
    f.write('Step.add.Model.Sources                   = [@ADDED_DATA_SOURCE]\n')  ## comes from ft_allbands, then copy over columns
    f.write('Step.add.Model.Cache.Enable              = T\n')
    f.write('Step.add.Model.Phasors.Enable            = F\n')
    f.write('Step.add.Model.DirectionalGain.Enable    = F\n')
    f.write('Step.add.Model.Gain.Enable               = T\n')
    f.write('Step.add.Model.Rotation.Enable           = F\n')
    f.write('Step.add.Model.CommonScalarPhase.Enable  = T\n')
    if TEC   == "True":
        f.write('Step.add.Model.TEC.Enable  = T\n')
    if clock == "True":
        f.write('Step.add.Model.Clock.Enable  = T\n')
    f.write('Step.add.Model.CommonRotation.Enable     = F\n')
    f.write('Step.add.Operation                       = ADD\n')
    f.write('Step.add.Model.Beam.Enable               = F\n')
    f.write('Step.add.Output.WriteCovariance          = F\n')
    f.write('Step.add.Output.Column                   = %s\n' % outputcolumn)
    f.close()
    return bbs_parset


def join_parmdb(ms, parmdb_selfcal,parmdb_nondde, parmdb_template, parmdb_out, TEC, clock):
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

    return



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

    #depth =  0.7 / (numpy.sqrt(numpy.float(len(mslist))))
    depth =  0.5 / (numpy.sqrt(numpy.float(len(mslist))))
    cleandepth1 = str(depth*1.5) + 'mJy'
    cleandepth2 = str(depth)    + 'mJy'

    # get size from mask image
    if imsize is None:
        imsize = image_size_from_mask(inputmask)

    # already done in facetting step
    ## speed up the imaging if possible by reducing image size within the mask region
    #newsize = find_newsize(inputmask)
    #if newsize < imsize: # ok so we can use a smaller image size then
        ##make a new template
        #os.system('casapy --nologger -c '+SCRIPTPATH+'/make_empty_image.py '+ str(mslist[0]) + ' ' + inputmask+'2' + ' ' + str(newsize) + ' ' +'1.5arcsec')
        #os.system('casapy --nologger -c '+SCRIPTPATH+'/regrid_image.py '    + inputmask      + ' ' + inputmask+'2' + ' ' + inputmask+'3')

        ## reset the imsize and the mask
        #imsize    = newsize
        #inputmask = inputmask+'3'


    ms = ''
    for m in mslist:
        ms = ms + ' ' + m

    imout = 'im'+ callnumber +'_cluster'+cluster+'nm'

    os.system('casapy --nologger --logfile casapy-'+imout+'.log -c '+SCRIPTPATH+'/dde_weeren/a2256_hba/casapy_cleanv4.py ' + ms + ' ' + imout + ' ' + 'None' +\
               ' ' + cleandepth1 + ' ' + str(niter) + ' ' + str(nterms) + ' ' + str(imsize) + ' ' + mscale)


    # make mask
    if nterms > 1:
        os.system('python '+SCRIPTPATH+'/makecleanmask_field.py --threshpix '+str(threshpix)+\
                    ' --threshisl '+str(threshisl) +' --atrous_do '+ str(atrous_do) +' '   +imout +'.image.tt0')
    else:
        os.system('python '+SCRIPTPATH+'/makecleanmask_field.py --threshpix '+str(threshpix)+\
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

    #niter = numpy.int(niter*1.2) # clean a bit deeper (will actually be quite a bit deeper because of mask)
    niter = numpy.int(niter*2.0) # clean a bit deeper (will actually be quite a bit deeper because of mask)
    imout = 'im'+ callnumber +'_cluster'+cluster


    if region != 'empty': # in that case we have a extra region file for the clean mask
        if 'fs' in region:
            niter = numpy.int(niter*3.0) # increase niter, tune manually if needed (change back)
        else:
            niter = numpy.int(niter*3.0) # increase niter, tune manually if needed
        os.system('casapy --nologger --logfile casapy-'+imout+'.log -c '+SCRIPTPATH+'/dde_weeren/a2256_hba/casapy_cleanv4.py '+ ms + ' ' + imout + ' ' + mask_sources+'field,'+region + \
                  ' ' + cleandepth2 + ' ' + str(niter) + ' ' + str(nterms) + ' ' + str(imsize) + ' ' + mscale)

    else:
        os.system('casapy --nologger --logfile casapy-'+imout+'.log -c '+SCRIPTPATH+'/dde_weeren/a2256_hba/casapy_cleanv4.py '+ ms + ' ' + imout + ' ' + mask_sources+'field' + \
                   ' ' + cleandepth2 + ' ' + str(niter) + ' ' + str(nterms) + ' ' + str(imsize) + ' ' + mscale)

    # convert to FITS
    if nterms > 1:
        os.system('image2fits in=' + imout +'.image.tt0' + ' ' + 'out='+ imout + '.fits')
    else:
        os.system('image2fits in=' + imout +'.image'     + ' ' + 'out='+ imout + '.fits')

    # if the output image file (fits) does not exist - something went wrong - DIE
    if not os.path.exists(imout+'.fits'):
        raise Exception('Imaging Error: '+imout+'.fits does not exist' )


    return imout, mask_sources+'field', imsize

def image_size_from_mask(mask):
    im = pyrap.images.image(mask)
    sh = im.shape()
    if sh[-1] != sh[-2]:
        print "image is not square!"
        print sh[-1], sh[-2]
    npix = sh[-1]
    return npix

def rad2deg(rad):
    return rad*180./numpy.pi

def image_centre_from_mask(mask):
    im = pyrap.images.image(mask)
    ic = im.coordinates()
    lon,lat = ic.get_referencevalue()[2]

    ra = rad2deg(lat)
    dec = rad2deg(lon)

    centre = '{x:.4f},{y:.4f}'.format(x=ra,y=dec)
    facet_centre = '{x:s},{y:s}'.format(x=ra_to_str(ra,delim='h'),y=dec_to_str(dec,delim='d'))
    print centre
    print facet_centre
    return facet_centre

def make_image_wsclean_nomask(mslist, cluster, callnumber, threshpix, threshisl, nterms, atrous_do, imsize, inputmask, mscale, region,cellsize,uvrange,wsclean,WSCleanRobust):

    if imsize is None:
        imsize = image_size_from_mask(inputmask)

    niter   = numpy.int(20000 * (numpy.sqrt(numpy.float(len(mslist)))))
    cellsizeim = str(cellsize) +'arcsec'

    depth =  1e-3*0.7 / (numpy.sqrt(numpy.float(len(mslist))))
    print 'Cleaning to a noise level of',depth,'Jy: niter is',niter

    cleandepth2 = str(depth)     #+ 'mJy'

    wideband = False
    if len(mslist) > 5:
        wideband = True

    ms = ''
    for m in mslist:
        ms = ms + ' ' + m

    imout = 'im'+ callnumber +'_cluster'+cluster

    os.system('rm -rf ' + imout + '-*')


    outms      = 'field-'+cluster+'.ms'
    parsetname = 'concatforwsclean-'+cluster+'.parset'

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
        cmd2 = '-scale ' + cellsizeim + ' -weight briggs '+str(WSCleanRobust)+' -niter ' + str(niter) + '-cleanborder 0 -threshold '+ cleandepth2 + ' '
        cmd3 = '-minuv-l '+ str(uvrange) \
               +' -mgain 0.75 -fitbeam -datacolumn DATA -no-update-model-required -joinchannels -channelsout ' +\
               str(channelsout) + ' '  + outms
    else:
        cmd1 = wsclean + ' -reorder -name ' + imout + ' -size ' + str(imsize) + ' ' + str(imsize) + ' '
        cmd2 = '-scale ' + cellsizeim + ' -weight briggs '+str(WSCleanRobust)+' -niter ' + str(niter) + ' -cleanborder 0 -threshold '+ cleandepth2 + ' '
        cmd3 = '-minuv-l '+ str(uvrange) +' -mgain 0.75 -fitbeam -datacolumn DATA -no-update-model-required ' + outms

    print cmd1+cmd2+cmd3
    os.system(cmd1+cmd2+cmd3)

    finalim=pyrap.images.image(imout+'-image.fits')
    finalim.saveas(imout +'.image')

    return imout, None, imsize


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

def make_image_wsclean(mslist, cluster, callnumber, threshpix, threshisl, nterms, atrous_do, imsize, inputmask, mscale, region,cellsize,uvrange,wsclean,WSCleanRobust,BlankField):

    if imsize is None:
        imsize = image_size_from_mask(inputmask)

    niter   = numpy.int(5000 * (numpy.sqrt(numpy.float(len(mslist)))))
    cellsizeim = str(cellsize) +'arcsec'

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

#############
# MAIN CODE #
#############

if __name__=='__main__':

    logging.basicConfig(filename='dde_reimage.log',level=logging.DEBUG, format='%(asctime)s -  %(message)s', datefmt='%Y-%d-%m %H:%M:%S')

    if len(sys.argv)<2:
        raise Exception('Give the path to the setup code for the facet')

    username = pwd.getpwuid(os.getuid())[0]
    print 'Using',sys.argv[1],'as the setup code'
    execfile(sys.argv[1])
    print 'script path is',SCRIPTPATH
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

    try:
        NoMask
    except NameError:
        NoMask=False
        print 'NoMask not set, defaulting to',NoMask

    source_info_rec = numpy.genfromtxt(peelsourceinfo, dtype="S10,S25,S5,S5,i8,i8,i8,i8,S2,S10,S10", names=["sourcelist","directions","atrous_do","mscale_field","imsizes","cellsizetime_p","cellsizetime_a","fieldsize","dynamicrange","regionselfc","regionfield"])

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

    sourcelist = sourcelist.tolist()



    mslistorig = ["{name:s}_SB{b1:03d}-{b2:03d}.{res:s}.ms".format(name=NAME,res=RES,b1=b,b2=b+9) for b in BANDS]
    mslistorigstr = ' '.join(mslistorig)

    mslist= [ms for ms in mslistorig if  os.path.isdir(ms)]  # filter out datasets that do not exist (takes also care of freq.gaps in field subtract))

    msliststr = ' '.join(mslist)



    tt = pt.table(mslist[0] + '/FIELD')
    pointingcenter = tt.getcol('REFERENCE_DIR')[0][0]
    pointingcenter = str(pointingcenter[0]) +'rad,' + str(pointingcenter[1])+'rad'
    print pointingcenter
    tt.close()

    if len(mslist) == 1:
        TEC    = "False" # no TEC fitting for one (channel) dataset
        nterms = 1
    if len(mslist) > 16:
        nterms = 2
    if len(mslist) > 40:
        nterms = 3
    if len(mslist) < 20:
        clock = "False"


    ##################################

    logging.info('\n')
    logging.info('#######################################################################\n')
    logging.info('Running DDE-reimage')
    logging.info('Doing sources: '+','.join(do_sources))
    logging.info('Using MSlist: '+msliststr)

    for source in do_sources:
        ## I1

        source_id = sourcelist.index(source)

        print 'Re-imaging facet:', source
        logging.info('')
        logging.info('Re-imaging facet: '+ source )

        logging.info("removing any existing facet field1 images")
        os.system("rm -rf imfield1*_cluster"+source+"*")
        logging.info("removing any existing facet field1 imaging average MS")
        os.system("rm -rf *."+source+".ms.avgfield1*")

        #check if allbands.concat.shifted_'+source+'.ms is present
        if os.path.isdir('allbands.concat.shifted_'+source+'.ms'):
            logging.info('allbands.concat.shifted_'+source+'.ms already exists')
            logging.info('...but continuing because we are re-imaging')

        if not os.path.isdir('allbands.concat.ms'):
            print 'allbands.concat.ms does not exist'
            raise Exception('make measurement set and then restart')


        dummyskymodel   = SCRIPTPATH+'/dummy.skymodel' ## update every time again with new source, not used, just a dummy for correct

        msfieldavgfacetlist1 = []
        for ms_id, ms in enumerate(mslist):
            msfieldavgfacetlist1.append(ms.split('.')[0] + '.' + source + '.ms.avgfield1.facetdir')

        # combined SC and DDE solutions (basename - for use in runbbs - uses ms/basename)
        parmdb_master_out  = "instrument_master_" + source   # reset because runbbs uses basename of ms

        # set image mask region
        output_template_im = 'templatemask_' + source +'.masktmp'
        if not os.path.exists(output_template_im):
            raise  Exception('facet mask missing: '+output_template_im)

        # set directions #
        #selfcaldir = directions[source_id]
        facetdir = image_centre_from_mask(output_template_im)
        facetsize = image_size_from_mask(output_template_im)

        #logging.info("Selfcal direction: "+selfcaldir)
        logging.info("Facet direction: "+facetdir)
        logging.info("facetmask: "+str(facetsize))

        #continue

        ## NOTE: addback now needs to be done as subtract in reverse... with FT and allbands concat


        ## STEP 1: check  ##
        # if we didn't keep the allbands.concat.shifted
        if not os.path.exists('allbands.concat.shifted_'+source+'.ms'):
            parset = create_phaseshift_parset_full('allbands.concat.ms', 'allbands.concat.shifted_'+source+'.ms', facetdir,'DATA')

            ndppplog = parset.replace('.parset','.log')
            #ndpppcmd = 'NDPPP ' + parset + ' > '+ ndppplog + ' 2>&1'
            #ndppprc = os.system(ndpppcmd)
            ndpppcmd = 'NDPPP ' + parset + ' > '+ ndppplog + ' 2>&1'
            print ndpppcmd
            with open(ndppplog,'w') as f:
                ndppp_proc = Popen(ndpppcmd.split(), stdout=f, stderr=STDOUT ).wait()   # DONT run in background
            if ndppp_proc != 0:
                raise Exception("NDPPP failed: "+ndpppcmd+' > '+ndppplog)


        ### STEP 3: prep for facet ##

        # image model to add back
        imsizef = image_size_from_mask(output_template_im)
        imout = 'im'+ 'field0' +'_cluster'+source

        if not os.path.exists(imout+'-model.fits'):
            raise Exception(imout+'-model.fits is missing')


        logging.info('running ft: '+imout)


        # DO THE FFT
        do_fieldFFT('allbands.concat.shifted_'+source+'.ms',imout, imsizef, cellsize, wsclean, mslist, WSCleanRobust)
        logging.info('FFTed model of DDE facet: ' + source)

        # SHIFT PHASE CENTER BACK TO ORIGINAL
        parset = create_phaseshift_parset_full('allbands.concat.shifted_'+source+'.ms', 'allbands.concat.shiftedback_'+source+'.ms',\
                                    pointingcenter,'MODEL_DATA')

        ndppplog = parset.replace('.parset','.log')
        ndpppcmd = 'NDPPP ' + parset  + ' > '+ ndppplog + ' 2>&1 '
        print ndpppcmd
        os.system(ndpppcmd)
        os.system('rm -rf allbands.concat.shifted_'+source+'.ms') # clean up

        # Add MODEL_DATA (allbands.concat.shiftedback_'+source+'.ms) into ADDED_DATA_SOURCE from mslist

        freq_tab1= pt.table('allbands.concat.ms' + '/SPECTRAL_WINDOW')
        numchan1    = freq_tab1.getcol('NUM_CHAN')
        freq_tab2= pt.table(mslist[0] + '/SPECTRAL_WINDOW')
        numchan2    = freq_tab2.getcol('NUM_CHAN')
        freq_tab1.close()
        freq_tab2.close()

        if (numchan1[0]) == (numchan2[0]*len(mslist)):
            os.system('python '+SCRIPTPATH+'/copy_over_columns.py '+ msliststr +\
                    ' ' +'allbands.concat.shiftedback_'+source+'.ms'+' ' + 'ADDED_DATA_SOURCE')
        else:
            os.system('python '+SCRIPTPATH+'/copy_over_columns.py '+ mslistorigstr +\
                    ' ' +'allbands.concat.shiftedback_'+source+'.ms'+' ' + 'ADDED_DATA_SOURCE')


        os.system('rm -rf allbands.concat.shiftedback_'+source+'.ms') # clean up

        addfieldparset   = create_add_parset_field('MODEL_DATA', TEC, clock) # SUBTRACTED_DATA_ALL + ADDED_DATA_SOURCE = MODEL_DATA


        logging.info('adding sources with solutions')
        runbbs(mslist, dummyskymodel, addfieldparset, parmdb_master_out+'_norm', False) # replace-sourcedb not needed since we use "@column"



        logging.info('correcting with solutions')
        ### STEP 3: prep for facet ##
        # apply master solutions, put in CORRECTED_DATA
        if TEC == "True":
            if clock == "True":
                runbbs(mslist, dummyskymodel, SCRIPTPATH+'/correctfield2+TEC+clock.parset',parmdb_master_out+'_norm', False)
            else:
                runbbs(mslist, dummyskymodel, SCRIPTPATH+'/correctfield2+TEC.parset',parmdb_master_out+'_norm', False)
        else:
            runbbs(mslist, dummyskymodel, SCRIPTPATH+'/correctfield2.parset',parmdb_master_out+'_norm', False)


        ###########################################################################
        # NDPPP phase shift, less averaging (NEW: run 2 in parallel)
        # CHANGE (v7) phaseshift to sc dir and avg
        # then phase shift to facetdir
        for ms_id, ms in enumerate(mslist):
            phaseshiftfieldparset = create_phaseshift_parset_field_avg(ms, msfieldavgfacetlist1[ms_id], source, facetdir)

            cmd = "ps -u " + username + " | grep NDPPP | wc -l"
            output=numpy.int(Popen(cmd, shell=True, stdout=PIPE).communicate()[0])

            while output > 1 : # max 2 processes (max 2, instead of 3, because we also phaseshift)
                time.sleep(10)
                output=numpy.int(Popen(cmd, shell=True, stdout=PIPE).communicate()[0])
                #pid = (Popen('pidof NDPPP', shell=True, stdout=PIPE).communicate()[0])
                #pid_list = pid.split(' ')
            # START NDPPP BECAUSE LESS/EQ 2 PROCESSES ARE RUNNING
            ndppplog = phaseshiftfieldparset.replace('.parset','.log')
            ndpppcmd = 'NDPPP ' + phaseshiftfieldparset + ' > '+ ndppplog + ' 2>&1 &'
            print ndpppcmd
            os.system(ndpppcmd)

        # Check if all NDPPP processes are finished
        output=numpy.int(Popen(cmd, shell=True, stdout=PIPE).communicate()[0])
        while output > 0 :
            time.sleep(10)
            output=numpy.int(Popen(cmd, shell=True, stdout=PIPE).communicate()[0])
            #pid = (Popen('pidof NDPPP', shell=True, stdout=PIPE).communicate()[0])
            #pid_list = pid.split(' ')


        logging.info('making image')
        ### STEP 4a: do facet ##
        # make large field image
        if NoMask:
            make_image_wsclean_nomask(msfieldavgfacetlist1, source, 'field1', 5, 3, nterms, 'True', None, output_template_im, mscale_field[source_id],regionfield[source_id],cellsize, uvrange,wsclean,WSCleanRobust)
        else:
            make_image_wsclean(msfieldavgfacetlist1, source, 'field1', 5, 3, nterms, 'True', None, output_template_im, mscale_field[source_id],regionfield[source_id],cellsize, uvrange,wsclean,WSCleanRobust,BlankField)
