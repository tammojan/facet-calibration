import matplotlib
matplotlib.use('GTK')
import numpy
import os
import sys
import lofar.parmdb
from scipy import interpolate
import time
from subprocess import Popen, PIPE, STDOUT
import pyrap.tables as pt
import uuid
pi = numpy.pi
import pwd
username = pwd.getpwuid(os.getuid())[0]

# location of this script - all scripts/parsets it uses are contained in subdirectory 'use'
SCRIPTPATH = os.path.dirname(sys.argv[0])

# to run casa when not logged in
# Step 0. Run in screen
# Step 1. Xvfb :5 &
# Step 2. setenv DISPLAY :5
# run the program in screen and to exit screen: screen -d

# v19
# - change ndppp job submit for better control over jobs running
# - add hardcoded groups (get_group)

def create_merged_parmdb_spline(parmdb_a,parmdb_p, parmdb_t, parmdb_out,cellsizetime_a, cellsizetime_b):


    pdb_a = lofar.parmdb.parmdb(parmdb_a)
    pdb_p = lofar.parmdb.parmdb(parmdb_p)
    pdb_t = lofar.parmdb.parmdb(parmdb_t)

    parms_a = pdb_a.getValuesGrid("*")
    parms_p = pdb_p.getValuesGrid("*")
    parms_t = pdb_t.getValuesGrid("*")

    os.system('rm -rf ' + parmdb_out)

    keynames_p = parms_p.keys()
    keynames_a = parms_a.keys()

    for key in keynames_p:
        # copy over the scalar phases
        scalarp                     = numpy.copy(parms_p[key]['values'][:,0])
        parms_t[key]['values'][:,0] = scalarp


    for key in keynames_a:
        print key
        tmp1 = numpy.copy(parms_t[key]['values'][:,0])
        tmp2 = numpy.copy(parms_a[key]['values'][:,0])

        time_t = 1.0*numpy.arange(len(tmp1))
        time_a = (numpy.float(cellsizetime_a/cellsizetime_p)*numpy.arange(len(tmp2))) + (0.5*cellsizetime_a/cellsizetime_p)

        el =   ((numpy.max(time_t) -  numpy.max(time_a[:-1]))/2 )+ numpy.max(time_a[:-1])
        time_a[len(time_a)-1] = el

        # to avoid edge effects in spline add values edges to be the last values from BBS
        time_a = numpy.copy(numpy.append(0,time_a))
        time_a = numpy.append(time_a, numpy.max(time_t))

        val_start = tmp2[0]
        val_end   = tmp2[len(tmp2)-1]

        tmp2 = numpy.copy(numpy.append(val_start,tmp2))
        tmp2 = numpy.append(tmp2,val_end)

        x = numpy.copy(time_a)
        y = numpy.copy(tmp2)


        # SPLINE INTERPOL of AMPS
        tck = interpolate.splrep(x,y,s=0)
        xnew = time_t
        ynew = interpolate.splev(xnew,tck,der=0)


        parms_t[key]['values'][:,0] = ynew

    pdbnew = lofar.parmdb.parmdb(parmdb_out, create=True)
    pdbnew.addValues(parms_t)
    pdbnew.flush()
    #lofar.expion.parmdbmain.store_parms(parmdb_out, parms_t, create_new = True)
    return

def create_merged_parmdb(ms, parmdb_a,parmdb_p, parmdb_t, parmdb_out,cellsizetime_a, cellsizetime_b):

    #pdb_pre = lofar.parmdb.parmdb(pre_apply_parmdb)
    pdb_a   = lofar.parmdb.parmdb(parmdb_a)
    pdb_p   = lofar.parmdb.parmdb(parmdb_p)
    pdb_t   = lofar.parmdb.parmdb(parmdb_t)

    #parms_pre = pdb_pre.getValuesGrid("*")
    parms_a   = pdb_a.getValuesGrid("*")
    parms_p   = pdb_p.getValuesGrid("*")
    parms_t   = pdb_t.getValuesGrid("*")


    os.system('rm -rf ' + parmdb_out)

    keynames_p = parms_p.keys()
    keynames_a = parms_a.keys()

    length = numpy.int(cellsizetime_b)
    for key in keynames_p:
        # copy over the CommonScalar phases, TEC, clock

        # just to get axis length
        tmp1  = numpy.copy(parms_t[key]['values'][:, 0])
        for idx in range(len(tmp1)):
            el =  numpy.float(idx)/numpy.float(length)
            el =  numpy.int(numpy.floor(el))

            parms_t[key]['values'][idx,0] = numpy.copy(parms_p[key]['values'][el,0])


    pol_list = ['0:0','1:1']
    gain     = 'Gain'
    anttab     = pt.table(ms + '/ANTENNA')
    antenna_list    = anttab.getcol('NAME')
    anttab.close()



    #length = numpy.int(cellsizetime_a/cellsizetime_b)
    length = numpy.int(cellsizetime_a)

    for pol in pol_list:
        for antenna in antenna_list:
            #real1 = numpy.copy(parms_pre[gain + ':' + pol + ':Real:'+ antenna]['values'][:, 0])
            real2 = numpy.copy(parms_a[gain + ':' + pol + ':Real:'+ antenna]['values'][:, 0])
            #imag1 = numpy.copy(parms_pre[gain + ':' + pol + ':Imag:'+ antenna]['values'][:, 0])
            imag2 = numpy.copy(parms_a[gain + ':' + pol + ':Imag:'+ antenna]['values'][:, 0])

            # just to get axis length
            tmp1  = numpy.copy(parms_t[gain + ':' + pol + ':Imag:'+ antenna]['values'][:, 0])

            #G1 = real1 + 1j*imag1
            G2 = real2 + 1j*imag2
            Gnew = G2

            for idx in range(len(tmp1)):
                el =  numpy.float(idx)/numpy.float(length)
                el =  numpy.int(numpy.floor(el))

                parms_t[gain + ':' + pol + ':Imag:'+ antenna]['values'][idx,0] = numpy.copy(numpy.imag(Gnew[el]))
                parms_t[gain + ':' + pol + ':Real:'+ antenna]['values'][idx,0] = numpy.copy(numpy.real(Gnew[el]))


    #for key in keynames_a:
    #  print key
    #  tmp1 = numpy.copy(parms_t[key]['values'][:,0])
    #  tmp2 = numpy.copy(parms_a[key]['values'][:,0])
    #  length = numpy.int(cellsizetime_a/cellsizetime_b)

    #  for idx in range(len(tmp1)):
    #
    #    el =  numpy.float(idx)/numpy.float(length)
    #    el =  numpy.int(numpy.floor(el))
    #    #print idx, el
    #    parms_t[key]['values'][idx,0] = tmp2[el]

    #lofar.expion.parmdbmain.store_parms(parmdb_out, parms_t, create_new = True)
    pdbnew = lofar.parmdb.parmdb(parmdb_out, create=True)
    pdbnew.addValues(parms_t)
    pdbnew.flush()

    return



def make_image(mslist, cluster, callnumber, threshpix, threshisl, nterms, atrous_do, imsize):

    do_mask = True
    niter   = 1000 # 7500 causes nasty clean artifacts
    mscale  = 'False'
    if atrous_do == 'True':
        mscale = 'True'

    average = True  # average the data a lot to speed up the imaging,
    # ONLY average for small FOV, otherwise timesmearing is a problem
    #average data
    if average:


        processes = set()
        max_processes = 2

        for ms in (mslist):
            ndppp_parset = ms + '_NDPPP.parset'
            ndppplog = ndppp_parset.replace('.parset','.log')
            os.system('rm -f ' + ndppp_parset)
            output = ms + '.tmpavg'
            f=open(ndppp_parset, 'w')
            f.write('msin = %s\n' % ms)
            f.write('msin.datacolumn = CORRECTED_DATA\n')
            f.write('msout = %s\n' % output)
            f.write('msout.writefullresflag=False\n')
            f.write('steps=[avg]\n')
            f.write('rficonsole.type=aoflagger\n')
            f.write('avg.type = squash\n')
            f.write('avg.freqstep = 1\n')
            f.write('avg.timestep = 12\n')      # is the default
            f.close()

            #cmd = "ps -u "+ username+ " | grep NDPPP | wc -l"
            #output=numpy.int(Popen(cmd, shell=True, stdout=PIPE).communicate()[0])

            #while output > 2 :
              #time.sleep(10)
              #output=numpy.int(Popen(cmd, shell=True, stdout=PIPE).communicate()[0])

              #pid = (Popen('pgrep NDPPP', shell=True, stdout=PIPE).communicate()[0])
              ##print pid
              #pid_list = pid.split()
              #for pidnr in pid_list:
                ##print 'PID', pidnr
                #os.system('kill -s CONT ' + str(pidnr))

            ndpppcmd = 'NDPPP ' + ndppp_parset
            #print ndpppcmd
            #os.system (ndpppcmd)

            with open(ndppplog,'w') as f:
                print 'exec: {cmd} > {log}'.format(cmd=ndpppcmd, log=ndppplog)
                p = Popen(ndpppcmd.split(), stdout=f, stderr=STDOUT)
                processes.add(p)
                pid=p.pid

            if len(processes) >= max_processes:
                os.wait()
                processes.difference_update([p for p in processes if p.poll() is not None])

        for p in processes:
            if p.poll() is None:
                p.wait()




    #output=numpy.int(Popen(cmd, shell=True, stdout=PIPE).communicate()[0])
    #while output > 0 :
            #time.sleep(10)
            #output=numpy.int(Popen(cmd, shell=True, stdout=PIPE).communicate()[0])
            #pid = (Popen('pgrep NDPPP', shell=True, stdout=PIPE).communicate()[0])
            ##print pid
            #pid_list = pid.split()
            #for pidnr in pid_list:
                ##print 'PID', numpy.int(pidnr)
                #os.system('kill -s CONT ' + str(pidnr))




    ms = ''
    for m in mslist:
        if average:
            ms = ms + ' ' + m + '.tmpavg'
        else:
            ms = ms + ' ' + m

    imout = 'im'+ callnumber +'_cluster'+cluster+'nm'
    print ms + ' ' + imout + ' ' + 'None' + ' ' + '1mJy' + ' ' + str(niter) + ' ' + str(nterms)

    if do_mask:
        if cluster == 'a2256': ## special case for a2256
            niter = niter*15 # clean very deep here

        os.system('casapy --nologger --logfile casapy-'+imout+'.log -c '+SCRIPTPATH+'/casapy_cleanv4.py ' + ms + ' ' + imout + ' ' + 'None' +\
                ' ' + '1mJy' + ' ' + str(niter) + ' ' + str(nterms) + ' ' + str(imsize) + ' ' + mscale)
        # make mask
        if nterms > 1:
            os.system('python '+SCRIPTPATH+'/makecleanmask.py --threshpix '+str(threshpix)+\
                      ' --threshisl '+str(threshisl) +' --atrous_do '+ str(atrous_do) +' '   +imout +'.image.tt0')
        else:
            os.system('python '+SCRIPTPATH+'/makecleanmask.py --threshpix '+str(threshpix)+\
                    ' --threshisl '+str(threshisl) +' --atrous_do '+ str(atrous_do) + ' '  + imout +'.image')

        # clean image with manual mask
        mask = imout+'.cleanmask'


    niter = 1000
    imout = 'im'+ callnumber +'_cluster'+cluster

    if region != 'empty' : ## special cases
        niter = niter*3
        os.system('casapy --nologger --logfile casapy-'+imout+'.log -c '+SCRIPTPATH+'/casapy_cleanv4.py '+ ms + ' ' + imout + ' ' + mask+','+region + \
                  ' ' + '1mJy' + ' ' + str(niter) + ' ' + str(nterms) + ' ' + str(imsize) + ' ' + mscale)
    else:
        os.system('casapy --nologger --logfile casapy-'+imout+'.log -c '+SCRIPTPATH+'/casapy_cleanv4.py '+ ms + ' ' + imout + ' ' + mask + \
                  ' ' + '1mJy' + ' ' + str(niter) + ' ' + str(nterms) + ' ' + str(imsize) + ' ' + mscale)

    # convert to FITS
    if nterms > 1:
        os.system('image2fits in=' + imout +'.image.tt0' + ' ' + 'out='+ imout + '.fits')
    else:
        os.system('image2fits in=' + imout +'.image'     + ' ' + 'out='+ imout + '.fits')

    if not os.path.exists(imout+'.fits'):
        raise Exception('Imaging Error: '+imout+'.fits does not exist' )


    if average:
        print 'rm -rf ' + ms
        os.system('rm -rf ' + ms)

    return imout,mask

def runbbs(mslist, skymodel, parset, parmdb, applycal, TEC):
    #NOTE WORK FROM MODEL_DATA (contains correct phase data from 10SB calibration)
    #NOTE WORK FROM DATA  if TEC

    if TEC:

        key = str(uuid.uuid4())
        key = key[0:12]

        # postgres db settings are different for paracluster
        if 'para' in os.uname()[1]:
            clusterdesc = '/home/wwilliams/para/paranew.clusterdesc'
            db = 'kolkje'
            dbuser = 'shimwell'
            dbname = 'shimwell'
        else:
            clusterdesc = '/home/williams/clusterdesc/cep3.clusterdesc'
            db = 'ldb002.offline.lofar'
            dbuser = 'postgres'
            dbname = username

        if len(mslist) == 10000:  # 34 does not work now!!!!!!
            # this is a very special case for a full run, manually here to run with 3 solver controls
            vdslist = ''
            os.system('makevds '+clusterdesc+' ' + ms)
            vdslist = vdslist + ms + '.vds '


        else:
            vdslist = ''
            for ms in mslist:
                os.system('makevds '+clusterdesc+' ' + ms)
                vdslist = vdslist + ms + '.vds '
            gds = 'tec.gds'
            print vdslist
            os.system('combinevds ' + gds + ' '+ vdslist)

            cmd1= 'calibrate -f --key ' + key + ' --cluster-desc '+clusterdesc +' --instrument-name ' + parmdb + ' '
            cmd2= '--db '+db + ' --db-user '+ dbuser + ' '+ ' --db-name '+ dbname + ' ' + gds + ' ' + parset + ' ' + skymodel + ' ' + '. > ' + gds + '.bbslog'
            bbscmd = cmd1 + cmd2

            ntries = 0
            done = 0
            while (ntries < 10)  and (done < 1):
                print 'calibrate try ', ntries
                print bbscmd
                os.system(bbscmd)

                #done = 0
                #while(done < 1):
                cmd = "grep '\[OK\] done!' " + gds + ".bbslog"
                output=Popen(cmd, shell=True, stdout=PIPE).communicate()[0]
                if 'OK' in output:
                    done += 1
                    print gds, 'is done'
                else:
                #cmd = "grep 'FAIL' " + gds + ".bbslog"
                #output=Popen(cmd, shell=True, stdout=PIPE).communicate()[0]
                #if 'FAIL' in output:
                    print 'calibrate failed, trying again...'
                    ntries += 1
                time.sleep(5)

            # after ntries exceeded is it still at FAIL

            cmd = "grep 'FAIL' " + gds + ".bbslog"
            output=Popen(cmd, shell=True, stdout=PIPE).communicate()[0]
            if 'FAIL' in output:
                raise Exception('calibrate has failed too many times' )








    else:
        for ms in mslist:
            log      =  ms + '.bbslog'
            if applycal:
                cmd = 'calibrate-stand-alone --parmdb-name ' + parmdb + ' ' + ms + ' ' + parset + ' ' + skymodel + '>' + log + '&'
                #cmd = 'calibrate-stand-alone -t 4 --parmdb-name ' + parmdb + ' ' + ms + ' ' + parset + ' ' + skymodel + '>' + log + '&'
            else:
                cmd = 'calibrate-stand-alone -f --parmdb-name ' + parmdb + ' ' + ms + ' ' + parset + ' ' + skymodel + '>' + log + '&'
                #cmd = 'calibrate-stand-alone -t 4 -f --parmdb-name ' + parmdb + ' ' + ms + ' ' + parset + ' ' + skymodel + '>' + log + '&'
            print cmd
            os.system(cmd)
        time.sleep(10)

        done = 0
        while(done < len(mslist)):
            done = 0
            for ms in mslist:
                cmd = "grep 'INFO - bbs-reducer terminated successfully.' " + ms + ".bbslog"
                output=Popen(cmd, shell=True, stdout=PIPE).communicate()[0]
                if output[0:4] == 'INFO':
                    done = done + 1
                    print ms, 'is done'
            time.sleep(5)

    return


def create_scalarphase_parset(timestep, TEC, clock, groups, FFT, uvrange):
    bbs_parset = 'scalarphase.parset'
    os.system('rm -f ' + bbs_parset)
    f=open(bbs_parset, 'w')

    if TEC:
        f.write('Strategy.InputColumn = DATA\n')
        f.write('Strategy.UseSolver   = T\n')

    else:
        f.write('Strategy.InputColumn = MODEL_DATA\n')
        f.write('Strategy.UseSolver   = F\n')

    f.write('Strategy.ChunkSize   = 960\n')
    f.write('Strategy.Steps       = [solve,correct]\n')

    if FFT:
        f.write('Step.solve.Model.Sources                = [@MODEL_DATA]\n')
    else:
        f.write('Step.solve.Model.Sources                = []\n')
    f.write('Step.solve.Model.Cache.Enable           = T\n')
    f.write('Step.solve.Model.Phasors.Enable         = F\n')
    f.write('Step.solve.Model.DirectionalGain.Enable = F\n')
    f.write('Step.solve.Model.Gain.Enable            = F\n')

    if TEC:
        f.write('Step.solve.Model.TEC.Enable            = T\n')
        f.write('Step.solve.Solve.CalibrationGroups     = [%s]\n' % groups)
    if clock:
        f.write('Step.solve.Model.Clock.Enable          = T\n')

    f.write('Step.solve.Model.Rotation.Enable        = F\n')
    f.write('Step.solve.Model.CommonScalarPhase.Enable = T \n')
    f.write('Step.solve.Operation                    = SOLVE\n')

    if TEC:
        if clock:
            f.write('Step.solve.Solve.Parms                  = ["CommonScalarPhase:*","TEC:*","Clock:*"]\n')
        else:
            f.write('Step.solve.Solve.Parms                  = ["CommonScalarPhase:*","TEC:*"]\n')
    else:
        f.write('Step.solve.Solve.Parms                  = ["CommonScalarPhase:*"]\n')

    f.write('Step.solve.Solve.CellSize.Freq          = 0\n')
    f.write('Step.solve.Solve.CellSize.Time          = %s\n' % str(timestep))
    f.write('Step.solve.Solve.CellChunkSize          = 120\n')
    f.write('Step.solve.Solve.PropagateSolutions     = T\n')
    f.write('Step.solve.Solve.Options.MaxIter        = 100\n')
    f.write('Step.solve.Solve.Options.LMFactor       = 1.0\n')
    f.write('Step.solve.Solve.Options.BalancedEqs    = F\n')
    f.write('Step.solve.Solve.Options.UseSVD         = T\n')
    f.write('Step.solve.Model.Beam.Enable            = F\n')
    f.write('Step.solve.Solve.UVRange                = [%s]\n' % uvrange)
    f.write('Step.solve.Solve.Mode                   = COMPLEX \n')

    f.write('Step.correct.Model.Sources                 = []\n')
    f.write('Step.correct.Model.CommonScalarPhase.Enable= T\n')
    f.write('Step.correct.Model.Cache.Enable            = T\n')

    if TEC:
        f.write('Step.correct.Model.TEC.Enable  = T\n')
    if clock:
        f.write('Step.correct.Model.Clock.Enable  = T\n')

    f.write('Step.correct.Model.DirectionalGain.Enable  = F\n')
    f.write('Step.correct.Model.Gain.Enable             = F\n')
    f.write('Step.correct.Model.Phasors.Enable          = F\n')
    f.write('Step.correct.Operation                     = CORRECT\n')
    f.write('Step.correct.Output.Column                 = CORRECTED_DATA\n')
    f.write('Step.correct.Model.Beam.Enable             = F\n')
    f.write('Step.correct.Output.WriteCovariance        = F\n')

    f.close()
    return bbs_parset

def create_scalarphase_parset_p(timestep, TEC, clock, groups, FFT, uvrange):
    bbs_parset = 'scalarphase_p.parset'
    os.system('rm -f ' + bbs_parset)
    f=open(bbs_parset, 'w')
    if TEC:
        f.write('Strategy.InputColumn = DATA\n')
        f.write('Strategy.UseSolver   = T\n')
    else:
        f.write('Strategy.InputColumn = MODEL_DATA\n')
        f.write('Strategy.UseSolver   = F\n')

    f.write('Strategy.ChunkSize   = 960\n')

    f.write('Strategy.Steps       = [solve,correct]\n')

    if FFT:
        f.write('Step.solve.Model.Sources                = [@MODEL_DATA]\n')
    else:
        f.write('Step.solve.Model.Sources                = []\n')
    f.write('Step.solve.Model.Cache.Enable           = T\n')
    f.write('Step.solve.Model.Phasors.Enable         = F\n')
    f.write('Step.solve.Model.DirectionalGain.Enable = F\n')
    f.write('Step.solve.Model.Gain.Enable            = F\n')

    if TEC:
        f.write('Step.solve.Model.TEC.Enable            = T\n')
        f.write('Step.solve.Solve.CalibrationGroups     = [%s]\n' % groups)
    if clock:
        f.write('Step.solve.Model.Clock.Enable          = T\n')

    f.write('Step.solve.Model.Rotation.Enable        = F\n')
    f.write('Step.solve.Model.CommonScalarPhase.Enable = T \n')
    f.write('Step.solve.Operation                    = SOLVE\n')

    if TEC:
        if clock:
            f.write('Step.solve.Solve.Parms                  = ["CommonScalarPhase:*","TEC:*","Clock:*"]\n')
        else:
            f.write('Step.solve.Solve.Parms                  = ["CommonScalarPhase:*","TEC:*"]\n')
    else:
        f.write('Step.solve.Solve.Parms                  = ["CommonScalarPhase:*"]\n')

    f.write('Step.solve.Solve.CellSize.Freq          = 0\n')
    f.write('Step.solve.Solve.CellSize.Time          = %s\n' % str(timestep))
    f.write('Step.solve.Solve.CellChunkSize          = 120\n')
    f.write('Step.solve.Solve.PropagateSolutions     = T\n')
    f.write('Step.solve.Solve.Options.MaxIter        = 100\n')
    f.write('Step.solve.Solve.Options.LMFactor       = 1.0\n')
    f.write('Step.solve.Solve.Options.BalancedEqs    = F\n')
    f.write('Step.solve.Solve.Options.UseSVD         = T\n')
    f.write('Step.solve.Model.Beam.Enable            = F\n')
    f.write('Step.solve.Solve.UVRange                = [%s]\n'  % uvrange)
    f.write('Step.solve.Solve.Mode                   = COMPLEX \n')

    f.write('Step.correct.Model.Sources                 = []\n')
    f.write('Step.correct.Model.CommonScalarPhase.Enable= T\n')
    f.write('Step.correct.Model.Cache.Enable            = T\n')
    f.write('Step.correct.Model.DirectionalGain.Enable  = F\n')

    if TEC:
        f.write('Step.correct.Model.TEC.Enable  = T\n')
    if clock:
        f.write('Step.correct.Model.Clock.Enable  = T\n')

    f.write('Step.correct.Model.Gain.Enable             = F\n')
    f.write('Step.correct.Model.Phasors.Enable          = F\n')
    f.write('Step.correct.Operation                     = CORRECT\n')
    f.write('Step.correct.Output.Column                 = CORRECTED_DATA_PHASE\n')
    f.write('Step.correct.Model.Beam.Enable             = F\n')
    f.write('Step.correct.Output.WriteCovariance        = F\n')

    f.close()
    return bbs_parset

def create_amponly_parset(timestep, FFT, uvrange):
    bbs_parset = 'amplitudeonly.parset'
    os.system('rm -f ' + bbs_parset)
    f=open(bbs_parset, 'w')

    f.write('Strategy.InputColumn = CORRECTED_DATA_PHASE\n')
    f.write('Strategy.ChunkSize   = 1440\n')
    f.write('Strategy.UseSolver   = F\n')
    f.write('Strategy.Steps       = [solve]\n')


    if FFT:
        f.write('Step.solve.Model.Sources                = [@MODEL_DATA]\n')
    else:
        f.write('Step.solve.Model.Sources                  = []\n')
    f.write('Step.solve.Model.Cache.Enable             = T\n')
    f.write('Step.solve.Model.Phasors.Enable           = F\n')   ## T
    f.write('Step.solve.Model.DirectionalGain.Enable   = F\n')
    f.write('Step.solve.Model.Gain.Enable              = T\n')
    f.write('Step.solve.Model.Rotation.Enable          = F\n')
    f.write('Step.solve.Model.CommonScalarPhase.Enable = F\n')
    f.write('Step.solve.Operation                      = SOLVE\n')
    f.write('Step.solve.Solve.Parms                    = ["Gain:0:0:*","Gain:1:1:*"]\n')   ##  ["Gain:0:0:Ampl:*","Gain:1:1:Ampl:*"]
    #  to test
    f.write('Step.solve.Solve.CellSize.Freq            = 0\n')
    f.write('Step.solve.Solve.CellSize.Time            = %s\n' %str(timestep))
    f.write('Step.solve.Solve.CellChunkSize            = 12\n')
    f.write('Step.solve.Solve.PropagateSolutions       = T\n')
    f.write('Step.solve.Solve.Options.MaxIter          = 100\n')
    f.write('Step.solve.Solve.Options.LMFactor         = 1.0\n')
    f.write('Step.solve.Solve.Options.BalancedEqs      = F\n')
    f.write('Step.solve.Solve.Options.UseSVD           = T\n')
    f.write('Step.solve.Model.Beam.Enable              = F\n')
    f.write('Step.solve.Solve.UVRange                  = [%s]\n'  % uvrange) ## [600]
    f.write('Step.solve.Solve.Mode                     = COMPLEX\n')

    f.close()
    return bbs_parset

def create_scalarphase_parset_p2(timestep, TEC, clock, groups, FFT, uvrange):
    bbs_parset = 'scalarphase_p2.parset'
    os.system('rm -f ' + bbs_parset)
    f=open(bbs_parset, 'w')
    f.write('Strategy.InputColumn = CORRECTED_DATA_AMP\n')
    f.write('Strategy.ChunkSize   = 960\n')

    if TEC:
        f.write('Strategy.UseSolver   = T\n')
    else:
        f.write('Strategy.UseSolver   = F\n')

    f.write('Strategy.Steps       = [solve,correct]\n')

    if FFT:
        f.write('Step.solve.Model.Sources                = [@MODEL_DATA]\n')
    else:
        f.write('Step.solve.Model.Sources                = []\n')
    f.write('Step.solve.Model.Cache.Enable           = T\n')
    f.write('Step.solve.Model.Phasors.Enable         = F\n')
    f.write('Step.solve.Model.DirectionalGain.Enable = F\n')
    f.write('Step.solve.Model.Gain.Enable            = F\n')
    if TEC:
        f.write('Step.solve.Model.TEC.Enable            = T\n')
        f.write('Step.solve.Solve.CalibrationGroups     = [%s]\n' % groups)
    if clock:
        f.write('Step.solve.Model.Clock.Enable          = T\n')

    f.write('Step.solve.Model.Rotation.Enable        = F\n')
    f.write('Step.solve.Model.CommonScalarPhase.Enable = T \n')
    f.write('Step.solve.Operation                    = SOLVE\n')

    if TEC:
        if clock:
            f.write('Step.solve.Solve.Parms                  = ["CommonScalarPhase:*","TEC:*","Clock:*"]\n')
        else:
            f.write('Step.solve.Solve.Parms                  = ["CommonScalarPhase:*","TEC:*"]\n')
    else:
        f.write('Step.solve.Solve.Parms                  = ["CommonScalarPhase:*"]\n')

    f.write('Step.solve.Solve.CellSize.Freq          = 0\n')
    f.write('Step.solve.Solve.CellSize.Time          = %s\n' % str(timestep))
    f.write('Step.solve.Solve.CellChunkSize          = 120\n')
    f.write('Step.solve.Solve.PropagateSolutions     = T\n')
    f.write('Step.solve.Solve.Options.MaxIter        = 100\n')
    f.write('Step.solve.Solve.Options.LMFactor       = 1.0\n')
    f.write('Step.solve.Solve.Options.BalancedEqs    = F\n')
    f.write('Step.solve.Solve.Options.UseSVD         = T\n')
    f.write('Step.solve.Model.Beam.Enable            = F\n')
    f.write('Step.solve.Solve.UVRange                = [%s]\n'  % uvrange)
    f.write('Step.solve.Solve.Mode                   = COMPLEX \n')

    f.write('Step.correct.Model.Sources                 = []\n')
    f.write('Step.correct.Model.CommonScalarPhase.Enable= T\n')
    f.write('Step.correct.Model.Cache.Enable            = T\n')
    f.write('Step.correct.Model.DirectionalGain.Enable  = F\n')
    f.write('Step.correct.Model.Gain.Enable             = F\n')

    if TEC:
        f.write('Step.correct.Model.TEC.Enable  = T\n')
    if clock:
        f.write('Step.correct.Model.Clock.Enable  = T\n')

    f.write('Step.correct.Model.Phasors.Enable          = F\n')
    f.write('Step.correct.Operation                     = CORRECT\n')
    f.write('Step.correct.Output.Column                 = CORRECTED_DATA_PHASE\n')
    f.write('Step.correct.Model.Beam.Enable             = F\n')
    f.write('Step.correct.Output.WriteCovariance        = F\n')

    f.close()
    return bbs_parset

def get_group(thismslist):
    nms = len(thismslist)

    if nms > 37:
        raise Exception('too many blocks....')

    if   nms == 37:    group = "5,5,5,5,5,6,6"
    elif nms == 36:    group = "6,6,6,6,6,6"
    elif nms == 35:    group = "5,6,6,6,6,6"
    elif nms == 34:    group = "5,5,6,6,6,6"
    elif nms == 33:    group = "5,5,5,6,6,6"
    elif nms == 32:    group = "5,5,5,5,6,6"
    elif nms == 31:    group = "5,5,5,5,5,6"
    elif nms == 30:    group = "6,6,6,6,6"

    elif nms == 29:    group = "5,6,6,6,6"
    elif nms == 28:    group = "5,5,6,6,6"
    elif nms == 27:    group = "5,5,5,6,6"
    elif nms == 26:    group = "5,5,5,5,6"
    elif nms == 25:    group = "5,5,5,5,5"
    elif nms == 24:    group = "6,6,6,6"
    elif nms == 23:    group = "5,6,6,6"
    elif nms == 22:    group = "5,5,6,6,"
    elif nms == 21:    group = "5,5,5,6"
    elif nms == 20:    group = "5,5,5,5"

    elif nms == 19:    group = "4,5,5,5"
    elif nms == 18:    group = "6,6,6"
    elif nms == 17:    group = "5,6,6"
    elif nms == 16:    group = "5,5,6"
    elif nms == 15:    group = "5,5,5"
    elif nms == 14:    group = "4,5,5"
    elif nms == 13:    group = "4,4,5"
    elif nms == 12:    group = "6,6"
    elif nms == 11:    group = "5,6"
    elif nms == 10:    group = "5,5"

    elif nms ==  9:    group = "4,5"
    elif nms ==  8:    group = "4,4"
    elif nms ==  7:    group = "3,4"

    else:              group = str(nms)

    return group


el=len(sys.argv)

mslist    = sys.argv[1:el-10]
cluster   = str(sys.argv[el-10])
atrous_do = str(sys.argv[el-9])
imsize    = numpy.int(sys.argv[el-8])

nterms                  = numpy.int(sys.argv[el-7])  # only 1 to 3 is supported !!
cellsizetime_a          = numpy.int(sys.argv[el-6])
cellsizetime_p          = numpy.int(sys.argv[el-5])
TECi                    = str(sys.argv[el-4])
clocki                  = str(sys.argv[el-3])
HRi                     = str(sys.argv[el-2])
region                  = str(sys.argv[el-1])

TEC  = False
FFT  = False
clock= False
HR   = False # high dynamic range

if TECi == "True":
    TEC = True

if clocki == "True":
    clock = True

if HRi == "HD":
    HR = True

if imsize <=2048:
    wplanes = 1
    FFT = True  # FFT image into MODEL_DATA
    if imsize >= 512:
        wplanes = 64
    if imsize > 799:
        wplanes = 96
    if imsize > 1023:
        wplanes = 128
    if imsize > 1599:
        wplanes = 180
    if imsize > 1800:
        wplanes = 196
    if imsize > 2049:
        wplanes = 256
    #if imsize > 3000:
    #   wplanes = 448
    #if imsize > 4095:
    #   wplanes = 512


print 'mslist', mslist
print 'source', cluster
print 'atrous_do', atrous_do
print 'imsize', imsize

msinputlist = ''
for m in mslist:
    msinputlist = msinputlist + ' ' + m


#if len(mslist) == 29:
    #group = "9,10,10"
#elif len(mslist) == 28:
    #group = "7,7,7,7"
#elif len(mslist) == 20:
    ##group = "10,10"
    #group = "7,7,6"
#elif len(mslist) == 16:
    #group = "8,8"
#else:
    #group = str(len(mslist))


group = get_group(mslist)
print 'GROUP', group


uvrange = '80'
#uvrange = '400'

merge_parmdb = True
phasors      = False   # if true only solve for amps on long timescales
smooth       = False # seems that smooth does not help the selfcal (various reasons for that)
                     # 1. boundaries of flagged vs non-flagged data are sharp (should not be smoothed)
                     # 2. there sre some strong ampl various at low elevations
smooth       = True # sometimes almost 0.0 amplitude, causes ripples
phasezero    = True # reset phases from ap calibration



#### MAKE IMAGE 0 ###
imout,mask = make_image(mslist, cluster, '0', 10, 6, nterms, atrous_do, imsize)

#####################

### CALIBRATE WITH BBS PHASE ONLY 1 ###
# create skymodel for BBS
os.system(SCRIPTPATH+'/casapy2bbs.py -m '+ mask + ' ' +'-t ' + str(nterms)+ ' ' + imout+'.model ' +  imout+'.skymodel')
if FFT:
    os.system('casapy --nologger -c '+SCRIPTPATH+'/ft_v2.py ' + msinputlist + ' ' + imout+'.model' \
              + ' ' + str(nterms) + ' '+ str(wplanes))

# phase only calibrate
skymodel = imout+'.skymodel'
parset   = create_scalarphase_parset(cellsizetime_p, TEC, clock, group, FFT, uvrange)

runbbs(mslist, skymodel, parset, 'instrument', False, TEC)
#NOTE WORK FROM MODEL_DATA (contains correct phase data from 10SB calibration)
######################################


### MAKE IMAGE 1 ###
imout,mask = make_image(mslist, cluster, '1', 15, 15, nterms, atrous_do, imsize)
####################


### CALIBRATE WITH BBS PHASE ONLY 2 ###
# create skymodel for BBS
os.system(SCRIPTPATH+'/casapy2bbs.py -m '+ mask + ' ' +'-t ' + str(nterms)+ ' ' + imout+'.model ' +  imout+'.skymodel')
if FFT:
    os.system('casapy --nologger -c '+SCRIPTPATH+'/ft_v2.py ' + msinputlist + ' ' + imout+'.model' \
              + ' ' + str(nterms) + ' '+ str(wplanes))


# phase only calibrate
skymodel = imout+'.skymodel'
parset   = create_scalarphase_parset(cellsizetime_p, TEC, clock, group, FFT, uvrange)

runbbs(mslist, skymodel, parset, 'instrument', False, TEC) #NOTE WORK FROM MODEL_DATA (contains correct phase data from 10SB calibration)
######################################


### MAKE IMAGE 2 ###
imout,mask = make_image(mslist, cluster, '2', 15, 15, nterms, atrous_do, imsize)
####################



### CALIBRATE WITH BBS PHASE+AMP 1 ###
os.system(SCRIPTPATH+'/casapy2bbs.py -m '+ mask + ' ' +'-t ' + str(nterms)+ ' ' + imout+'.model ' +  imout+'.skymodel')
if FFT:
    os.system('casapy --nologger -c '+SCRIPTPATH+'/ft_v2.py ' + msinputlist + ' ' + imout+'.model' \
              + ' ' + str(nterms) + ' '+ str(wplanes))

skymodel = imout+'.skymodel'
parset   = create_scalarphase_parset_p(cellsizetime_p, TEC, clock, group, FFT, uvrange)
# solve +apply phases
runbbs(mslist, skymodel, parset, 'instrument_phase0', False, TEC)

# solve amps
parmdb = 'instrument_amps0'
parset = create_amponly_parset(cellsizetime_a, FFT, uvrange)
runbbs(mslist, skymodel, parset, parmdb, False, False)

for ms in mslist:
    # remove outliers from the solutions
    if phasors:
        os.system('python '+SCRIPTPATH+'/smoothcal_rx42.py ' + ms + ' ' + ms+'/'+parmdb + ' ' + ms+'/'+parmdb+'_smoothed'+' > '+ms+'_'+parmdb+'_smoothed.log')
    else:
        os.system('python '+SCRIPTPATH+'/smoothcal_rx42_nophasors.py ' + ms + ' ' + ms+'/'+parmdb + ' ' + ms+'/'+parmdb+'_smoothed'+' > '+ms+'_'+parmdb+'_smoothed.log')

# apply amps
if smooth:
    runbbs(mslist, skymodel,SCRIPTPATH+'/apply_amplitudeonly.parset', parmdb+'_smoothed', True, False)
else:
    runbbs(mslist, skymodel,SCRIPTPATH+'/apply_amplitudeonly.parset', parmdb, True, False)

### MAKE IMAGE 3 ###
imout,mask = make_image(mslist, cluster, '3', 10, 10, nterms, atrous_do, imsize)



#### CALIBRATE  BBS PHASE+AMP 2 ###
# make model
os.system(SCRIPTPATH+'/casapy2bbs.py -m '+ mask + ' ' +'-t ' + str(nterms)+ ' ' + imout+'.model ' +  imout+'.skymodel')
if FFT:
    os.system('casapy --nologger -c '+SCRIPTPATH+'/ft_v2.py ' + msinputlist + ' ' + imout+'.model' \
              + ' ' + str(nterms) + ' '+ str(wplanes))

#parmdb keep from previous step
skymodel = imout+'.skymodel'


# reset the phases from instrument_amps0 to zero to prevent large phase corrections from incorrect AP solve
if phasezero:
    inputparmdb  = parmdb +'_smoothed'
    outputparmdb = parmdb +'_smoothed_phasezero'
    for ms in mslist:
        os.system('python '+SCRIPTPATH+'/setphasezero.py ' + ms + ' ' + ms+'/'+inputparmdb +' ' + ms+'/'+outputparmdb)
else:
    outputparmdb = parmdb +'_smoothed'


# phase only cal
skymodel = imout+'.skymodel'
parset   = create_scalarphase_parset_p(cellsizetime_p, TEC, clock, group, FFT, uvrange)
runbbs(mslist, skymodel, parset, 'instrument_phase1', False, TEC)

# solve amps
parmdb   = 'instrument_amps1'
parset = create_amponly_parset(cellsizetime_a, FFT, uvrange)
runbbs(mslist, skymodel, parset,parmdb, False, False)

for ms in mslist:
    # remove outliers from the solutions
    if phasors:
        os.system('python '+SCRIPTPATH+'/smoothcal_rx42.py ' + ms + ' ' + ms+'/'+parmdb + ' ' + ms+'/'+parmdb+'_smoothed'+' > '+ms+'_'+parmdb+'_smoothed.log')
    else:
        os.system('python '+SCRIPTPATH+'/smoothcal_rx42_nophasors.py ' + ms + ' ' + ms+'/'+parmdb + ' ' + ms+'/'+parmdb+'_smoothed'+' > '+ms+'_'+parmdb+'_smoothed.log')

# apply amps
if smooth:
    runbbs(mslist, skymodel,SCRIPTPATH+'/apply_amplitudeonly.parset',parmdb+'_smoothed', True, False)
else:
    runbbs(mslist, skymodel,SCRIPTPATH+'/apply_amplitudeonly.parset',parmdb, True, False)

### MAKE IMAGE 4 ###
imout,mask = make_image(mslist, cluster, '4', 10, 10, nterms, atrous_do, imsize)


### CREATE FINAL MODEL ###
skymodelf= 'im_cluster'+cluster+ '.final.skymodel'
os.system(SCRIPTPATH+'/casapy2bbs.py -m '+ mask + ' ' +'-t ' + str(nterms)+ ' ' + imout+'.model ' +  skymodelf)
if FFT:
    os.system('casapy --nologger -c '+SCRIPTPATH+'/ft_v2.py ' + msinputlist + ' ' + imout+'.model' \
              + ' ' + str(nterms) + ' '+ str(wplanes))

### CREATED MERGED PARMDB SCALARPHASE+AMPS ###
### INCLUDES SPLINE INTERPOLARION OF AMPS ###
if merge_parmdb:

    if phasors:
        dummyparset = SCRIPTPATH+'/scalarphase+amp.parset'
    else:
        if TEC:
            dummyparset = SCRIPTPATH+'/scalarphase+ap+TEC.parset'
        else:
            dummyparset = SCRIPTPATH+'/scalarphase+ap.parset'

    if TEC:
        if clock:
            dummyparmdb = 'instrument_template_TECclock'
        else:
            dummyparmdb = 'instrument_template_Gain_TEC_CSphase'

    #if not os.path.isdir(dummyparmdb):
        #runbbs([mslist[0]], skymodel,dummyparset, dummyparmdb, True, False)
    # TO SPEED THINGS UP, hard coded for BOOTES - i.e. the above has already been run
    for ms in mslist:
        os.system('rm -rf ' + ms +'/' + dummyparmdb)
        os.system('cp -r ' + dummyparmdb + ' ' +  ms + '/instrument_template')

    if smooth:
        parmdb_a    = 'instrument_amps1_smoothed'  # last/best ampplitude(+phase) parmdb
    else:
        parmdb_a    = 'instrument_amps1'  # last/best ampplitude(+phase) parmdb
    parmdb_p    = 'instrument_phase1'          # last/best CommonScalarPhase parmdb
    parmdbout   = 'instrument_merged'

    #reset in case of instrument_template_TECclock
    dummyparmdb = 'instrument_template'

    for ms in mslist:
        create_merged_parmdb(ms, ms+'/'+parmdb_a, ms+'/'+parmdb_p, ms+'/'+dummyparmdb,ms+'/'+parmdbout,cellsizetime_a,cellsizetime_p)
