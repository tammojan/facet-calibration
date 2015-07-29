import matplotlib
matplotlib.use('GTK')
import numpy
import os
import sys
from scipy import interpolate
import time
from subprocess import Popen, PIPE
import pyrap.tables as pt
import pyrap.images as pi
import uuid
pi = numpy.pi
import pwd
username = pwd.getpwuid(os.getuid())[0]
from facet_utilities import run, bg, angsep, getcpu, getmem

SCRIPTPATH = os.path.dirname(os.path.abspath(__file__))


# to run casa when not logged in
# Step 0. Run in screen
# Step 1. Xvfb :5 &
# Step 2. setenv DISPLAY :5
# run the program in screen and to exit screen: screen -d

def robust_sigma(in_y, zero=0):
    """
    Calculate a resistant estimate of the dispersion of
    a distribution. For an uncontaminated distribution,
    this is identical to the standard deviation.

    Use the median absolute deviation as the initial
    estimate, then weight points using Tukey Biweight.
    See, for example, Understanding Robust and
    Exploratory Data Analysis, by Hoaglin, Mosteller
    and Tukey, John Wiley and Sons, 1983.

    .. note:: ROBUST_SIGMA routine from IDL ASTROLIB.

    Examples
    --------
    >>> result = robust_sigma(in_y, zero=1)

    Parameters
    ----------
    in_y : array_like
        Vector of quantity for which the dispersion is
        to be calculated

    zero : int
        If set, the dispersion is calculated w.r.t. 0.0
        rather than the central value of the vector. If
        Y is a vector of residuals, this should be set.

    Returns
    -------
    out_val : float
        Dispersion value. If failed, returns -1.

    """
    # Flatten array
    y = in_y.ravel()

    eps = 1.0E-20
    c1 = 0.6745
    c2 = 0.80
    c3 = 6.0
    c4 = 5.0
    c_err = -1.0
    min_points = 3

    if zero:
        y0 = 0.0
    else:
        y0 = numpy.median(y)

    dy    = y - y0
    del_y = abs( dy )

    # First, the median absolute deviation MAD about the median:

    mad = numpy.median( del_y ) / c1

    # If the MAD=0, try the MEAN absolute deviation:
    if mad < eps:
        mad = del_y.mean() / c2
    if mad < eps:
        return 0.0

    # Now the biweighted value:
    u  = dy / (c3 * mad)
    uu = u * u
    q  = numpy.where(uu <= 1.0)
    count = len(q[0])
    if count < min_points:
        module_logger.warn('ROBUST_SIGMA: This distribution is TOO WEIRD! '
                           'Returning {}'.format(c_err))
        return c_err

    numerator = numpy.sum( (y[q] - y0)**2.0 * (1.0 - uu[q])**4.0 )
    n    = y.size
    den1 = numpy.sum( (1.0 - uu[q]) * (1.0 - c4 * uu[q]) )
    siggma = n * numerator / ( den1 * (den1 - 1.0) )

    if siggma > 0:
        out_val = numpy.sqrt( siggma )
    else:
        out_val = 0.0

    return out_val

def meanclip(indata, clipsig=4.0, maxiter=10, converge_num=0.001, verbose=0):
   """
   Computes an iteratively sigma-clipped mean on a
   data set. Clipping is done about median, but mean
   is returned.

   .. note:: MYMEANCLIP routine from ACS library.

   :History:
       * 21/10/1998 Written by RSH, RITSS
       * 20/01/1999 Added SUBS, fixed misplaced paren on float call, improved doc. RSH
       * 24/11/2009 Converted to Python. PLL.

   Examples
   --------
   >>> mean, sigma = meanclip(indata)

   Parameters
   ----------
   indata: array_like
       Input data.

   clipsig: float
       Number of sigma at which to clip.

   maxiter: int
       Ceiling on number of clipping iterations.

   converge_num: float
       If the proportion of rejected pixels is less than
       this fraction, the iterations stop.

   verbose: {0, 1}
       Print messages to screen?

   Returns
   -------
   mean: float
       N-sigma clipped mean.

   sigma: float
       Standard deviation of remaining pixels.

   """
   # Flatten array
   skpix = indata.reshape( indata.size, )
 
   ct = indata.size
   iter = 0; c1 = 1.0 ; c2 = 0.0
 
   while (c1 >= c2) and (iter < maxiter):
       lastct = ct
       medval = numpy.median(skpix)
       sig = numpy.std(skpix)
       wsm = numpy.where( abs(skpix-medval) < clipsig*sig )
       ct = len(wsm[0])
       if ct > 0:
           skpix = skpix[wsm]
 
       c1 = abs(ct - lastct)
       c2 = converge_num * lastct
       iter += 1
   # End of while loop
 
   mean  = numpy.mean( skpix )
   sigma = robust_sigma( skpix )
 
   if verbose:
       prf = 'MEANCLIP:'
       print '%s %.1f-sigma clipped mean' % (prf, clipsig)
       print '%s Mean computed in %i iterations' % (prf, iter)
       print '%s Mean = %.6f, sigma = %.6f' % (prf, mean, sigma)
 
   return mean, sigma


def create_merged_parmdb_spline(parmdb_a,parmdb_p, parmdb_t, parmdb_out,cellsizetime_a, cellsizetime_b):
 import lofar.parmdb

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
 import lofar.parmdb
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
  
  

 pdbnew = lofar.parmdb.parmdb(parmdb_out, create=True)
 pdbnew.addValues(parms_t)
 pdbnew.flush()
  
 return

def create_merged_parmdb_scalarphase(ms, parmdb_p, parmdb_t, parmdb_out,cellsizetime_b):
 import lofar.parmdb
 pdb_p   = lofar.parmdb.parmdb(parmdb_p)
 pdb_t   = lofar.parmdb.parmdb(parmdb_t)
 parms_p   = pdb_p.getValuesGrid("*")
 parms_t   = pdb_t.getValuesGrid("*")
 
 os.system('rm -rf ' + parmdb_out)
 
 keynames_p = parms_p.keys()
 
 length = numpy.int(cellsizetime_b)
 for key in keynames_p:
   # copy over the CommonScalar phases
   
   # just to get axis length   
   tmp1  = numpy.copy(parms_t[key]['values'][:, 0]) 
   for idx in range(len(tmp1)):
       el =  numpy.float(idx)/numpy.float(length)
       el =  numpy.int(numpy.floor(el))  
       parms_t[key]['values'][idx,0] = numpy.copy(parms_p[key]['values'][el,0])

 pdbnew = lofar.parmdb.parmdb(parmdb_out, create=True)
 pdbnew.addValues(parms_t)
 pdbnew.flush()
  
 return


def make_image(mslist, cluster, callnumber, threshpix, threshisl, nterms, atrous_do, imsize, \
               column, uvrange, increaseSNR_via_freqcoverage, numchanperms, nblockconcat):

 do_mask = True
 niter   = 1000 # 7500 causes nasty clean artifacts
 mscale  = 'False'
 if atrous_do == 'True':
   mscale = 'True'
  
 average = True  # average the data a lot to speed up the imaging, 
 # ONLY average for small FOV, otherwise timesmearing is a problem
 #average data

 if average:
        b=bg(maxp=10)
        for ms in (mslist):
            ndppp_parset = ms + '_NDPPP.parset'
            ndppplog = ndppp_parset.replace('.parset','.log')
            os.system('rm -f ' + ndppp_parset)
            output = ms + '.tmpavg'
            os.system('rm -rf ' + output)
            f=open(ndppp_parset, 'w')
            f.write('msin = %s\n' % ms)
            f.write('msin.datacolumn = %s\n' % column)
            
            if increaseSNR_via_freqcoverage == 'True':
                f.write('msin.startchan = %s\n' % str(numchanperms*nblockconcat))
                f.write('msin.nchan = %s\n' % str(numchanperms))                
                
            f.write('msout = %s\n' % output)
            f.write('msout.writefullresflag=False\n') 
            f.write('steps=[avg]\n')
            f.write('avg.type = squash\n')
            f.write('avg.freqstep = 1\n')
    
            if imsize <= 800:  
                f.write('avg.timestep = 12\n')      # is the default
            else:  
                if imsize <= 1200:
                    f.write('avg.timestep = 6\n')
                else:
                    f.write('avg.timestep = 3\n') 
            f.close()


            ndpppcmd = 'NDPPP ' + ndppp_parset+ ' >'+ndppplog+' 2>&1'
            b.run(ndpppcmd)
            #os.system (ndpppcmd)

        b.wait()


 ms = ''
 for m in mslist:
   if average:
     ms = ms + ' ' + m + '.tmpavg'
   else:
     ms = ms + ' ' + m

 imout = 'im'+ callnumber +'_cluster'+cluster+'nm'
 print ms + ' ' + imout + ' ' + 'None' + ' ' + '1mJy' + ' ' + str(niter) + ' ' + str(nterms)
 os.system('rm -rf ' + imout + '.*')
 
 if do_mask:
   if cluster == 'a2256': ## special case for a2256
      niter = niter*15 # clean very deep here
      
   os.system('casapy --nogui -c ' + SCRIPTPATH +'/casapy_cleanv4.py ' + ms + ' ' + imout + ' ' + 'None' +\
             ' ' + '1mJy' + ' ' + str(niter) + ' ' + str(nterms) + ' ' + str(imsize) + ' ' + mscale + ' ' + uvrange + ' ' +\
	     str(cellsize))
   # make mask
   if nterms > 1:
     os.system('python ' + SCRIPTPATH + '/makecleanmask.py --threshpix '+str(threshpix)+\
               ' --threshisl '+str(threshisl) +' --atrous_do '+ str(atrous_do) +' '   +imout +'.image.tt0')
   else:
     os.system('python '+ SCRIPTPATH + '/makecleanmask.py --threshpix '+str(threshpix)+\
             ' --threshisl '+str(threshisl) +' --atrous_do '+ str(atrous_do) + ' '  + imout +'.image')

   # clean image with manual mask
   mask = imout+'.cleanmask'


 niter = 1000
 imout = 'im'+ callnumber +'_cluster'+cluster
 os.system('rm -rf ' + imout + '.*')
 
 if region != 'empty' : ## special cases
    niter = niter*3 # deeper clean in region as they are usualy used for extended sources
    os.system('casapy --nogui -c '+ SCRIPTPATH +'/casapy_cleanv4.py '+ ms + ' ' + imout + ' ' + mask+','+region + \
              ' ' + '1mJy' + ' ' + str(niter) + ' ' + str(nterms) + ' ' + str(imsize) + ' ' + mscale + ' ' + uvrange + ' ' +\
	      str(cellsize))
 else:
    os.system('casapy --nogui -c '+ SCRIPTPATH +'/casapy_cleanv4.py '+ ms + ' ' + imout + ' ' + mask + \
              ' ' + '1mJy' + ' ' + str(niter) + ' ' + str(nterms) + ' ' + str(imsize) + ' ' + mscale + ' ' + uvrange + ' ' +\
	     str(cellsize))
 
 # convert to FITS
 if nterms > 1:
   os.system('image2fits in=' + imout +'.image.tt0' + ' ' + 'out='+ imout + '.fits')
 else:
   os.system('image2fits in=' + imout +'.image'     + ' ' + 'out='+ imout + '.fits')
 
 
 if average:
   print 'rm -rf ' + ms
   os.system('rm -rf ' + ms)
 mask = imout+'.mask'
 return imout,mask


def run_stefcal(mslist, skymodel, parmdb, solvetype, FFT, uvrange, timestep, inputcol, \
                increaseSNR_via_freqcoverage, numchanperms, nblockconcat):
 # solve type can be "diagonal", "scalarphase"
 # if FFT == True then solve out of MODEL_DATA
  import pyrap.tables as pt
  freq_tab = pt.table(mslist[len(mslist)/2] + '/SPECTRAL_WINDOW')
  freq     = freq_tab.getcol('REF_FREQUENCY')[0]
  freq_tab.close()
  wav      = 299792458./freq
  
  timestep = numpy.int(timestep)
 
  # make the parset
  os.system('rm -f stefcal.parset')

  f=open('stefcal.parset', 'w')
  f.write('msin="tmp.ms"\n')
  f.write('msin.datacolumn=%s\n' % inputcol)
  if FFT:
    f.write('msin.modelcolumn=MODEL_DATA\n')
  
  f.write('msout=.\n')
  f.write('steps=[f,c]\n')
  
  f.write('f.type="filter"\n')
  f.write('f.blrange=[%s,1e30]\n' % (numpy.float(uvrange)*wav)) # convert to meters
  
  f.write('c.type="gaincal"\n')
  if FFT:
    f.write('c.usemodelcolumn=True\n')
  else:  
    skydb = skymodel + '.skydb'
    #print "makesourcedb in=" + skymodel + " out=" +skydb + " format='<'"
    os.system('rm -rf ' + skydb)
    os.system("makesourcedb in=" + skymodel + " out=" +skydb + " format='<'")
    f.write('c.usemodelcolumn=False\n')
    f.write('c.sourcedb=%s\n' % skydb)
  f.write('c.parmdb="tmp_instrumentdb"\n')
  f.write('c.maxiter=100\n')
  f.write('c.tolerance=1e-6\n')
  f.write('c.usebeammodel=False\n')
  f.write('c.usechannelfreq=False\n')
  f.write('c.caltype=%s\n' % solvetype) # can be diagonal, fulljones, phaseonly or scalarphase
  f.write('c.solint=%s\n' % timestep)
  f.write('c.debuglevel=1\n')
  f.close()


  if increaseSNR_via_freqcoverage == 'True':
      print 'Updating weights for solvetype', solvetype
      b=bg(maxp=numpy.int(128/(((2*nblockconcat) +1)*numchanperms))) # to limit number of taql processes
      for ms in mslist:
          if solvetype == 'phaseonly' or solvetype == 'scalarphase': 
              b.run("taql 'update " + ms + " set WEIGHT_SPECTRUM=WEIGHT_SPECTRUM_SC'")
          if solvetype == 'diagonal' or solvetype == 'fulljones':   
              b.run("taql 'update " + ms + " set WEIGHT_SPECTRUM=WEIGHT_SPECTRUM_AP'")
      # Check if all taql processes are finished
      b.wait()

  ndppp_stefcal=bg(maxp=32)
  for ms in mslist:
     log  =  ms + '.stefcallog'
     instrument = ms + '/' + parmdb
     os.system('rm -rf ' + instrument)
     cmd = 'NDPPP stefcal.parset ' + 'msin=' + ms + ' c.parmdb=' + instrument #+ ' & ' #>' + log
     print cmd
     ndppp_stefcal.run(cmd)  
     #os.system(cmd)
     #time.sleep(2)
  #time.sleep(30)
  ndppp_stefcal.wait()

  # check if NDPPP has finished
  #cmd = "ps -f -u "+ username+ " | grep NDPPP | grep stefcal |grep -v grep | wc -l"
  #output=numpy.int(Popen(cmd, shell=True, stdout=PIPE).communicate()[0])
  #while output > 0 :
  #  time.sleep(10)
  #  output=numpy.int(Popen(cmd, shell=True, stdout=PIPE).communicate()[0])
    
  return


def runbbs(mslist, skymodel, parset, parmdb, applycal, TEC):
 #NOTE WORK FROM MODEL_DATA (contains correct phase data from 10SB calibration) 
 #NOTE WORK FROM DATA  if TEC

 if TEC:
  if len(mslist) == 10000:  # 34 does not work now!!!!!!
   # this is a very special case for a full run, manually here to run with 3 solver controls 
   vdslist = ''
   os.system('makevds /home/'+username+'/cep2.clusterdesc ' + ms)
   vdslist = vdslist + ms + '.vds '
  
  
  else:
   vdslist = ''
   for ms in mslist:
     os.system('makevds /home/'+username+'/cep2.clusterdesc ' + ms)
     vdslist = vdslist + ms + '.vds '
   gds = 'tec.gds'  
   print vdslist
   os.system('combinevds ' + gds + ' '+ vdslist)
 
   key = str(uuid.uuid4())
   key = key[0:12]
   cmd1= 'calibrate -f --key ' + key + ' --cluster-desc /home/'+username+'/cep2.clusterdesc --instrument-name ' + parmdb + ' '
   cmd2= '--db ldb001 --db-user postgres ' + gds + ' ' + parset + ' ' + skymodel + ' ' + '. > ' + gds + '.bbslog'  
   cmd = cmd1 + cmd2
   print cmd
   os.system(cmd)



 else:  
   for ms in mslist:
     log      =  ms + '.bbslog'
     if applycal:
       cmd = 'calibrate-stand-alone --parmdb-name ' + parmdb + ' ' + ms + ' ' + parset + ' ' + skymodel + '>' + log + ' 2>&1 &'
       #cmd = 'calibrate-stand-alone -t 4 --parmdb-name ' + parmdb + ' ' + ms + ' ' + parset + ' ' + skymodel + '>' + log + '&'
     else:
       cmd = 'calibrate-stand-alone -f --parmdb-name ' + parmdb + ' ' + ms + ' ' + parset + ' ' + skymodel + '>' + log + ' 2>&1 &'
       #cmd = 'calibrate-stand-alone -t 4 -f --parmdb-name ' + parmdb + ' ' + ms + ' ' + parset + ' ' + skymodel + '>' + log + '&'
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
  f.write('Step.solve.Solve.UVRange		   = [%s]\n' % uvrange)
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
  f.write('Step.solve.Solve.UVRange		   = [%s]\n'  % uvrange)
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
  f.write('Step.solve.Solve.UVRange		     = [%s]\n'  % uvrange) ## [600]
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
  f.write('Step.solve.Solve.UVRange		   = [%s]\n'  % uvrange)
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

def find_imagenoise(imagename):
  import pyrap.images as pi
  im    = pi.image(imagename)
  image = numpy.copy(im.getdata())
  mean, rms =  meanclip(image)
  #im.close()
  return rms,  numpy.abs(numpy.max(image)/numpy.min(image))

def plotds9image(mslist, cluster):
  import glob
  
  sclmax = 0.04*numpy.sqrt(8.)/numpy.sqrt(numpy.float(len(mslist)))
  sclmin = -0.5*sclmax
  
  listofimages = glob.glob('im*_cluster' + cluster + '.fits')
  
  if len(listofimages) > 9: # to get the image order right in ds9
    cmd = 'ds9 im?_cluster' + cluster + '.fits ' +'im??_cluster' + cluster + '.fits'+ ' -view buttons no -view info no ' + \
          '-view panner no -view magnifier no -geometry 800x800 ' + \
          '-scale limits ' + str(sclmin) + ' ' + str(sclmax) + ' ' + \
          '-scale match -cmap bb -cmap match -zoom to fit -match frame wcs -colorbar no ' + \
          '-saveimage jpg 90 im_cluster' + cluster +'.jpg' +' -quit&'
  else:
    cmd = 'ds9 im?_cluster' + cluster + '.fits' + ' -view buttons no -view info no ' + \
          '-view panner no -view magnifier no -geometry 800x800 ' + \
          '-scale limits ' + str(sclmin) + ' ' + str(sclmax) + ' ' + \
          '-scale match -cmap bb -cmap match -zoom to fit -match frame wcs -colorbar no ' + \
          '-saveimage jpg 90 im_cluster' + cluster +'.jpg' +' -quit&'    
  print cmd
  os.system(cmd)  
  return


def concat_neighbor_block(basename, resolution, source, mslist, nblockconcat):
    ndp=bg(maxp=32)

    for ms in mslist:
        number = numpy.int(ms.split(basename)[1].split('.')[0].split('SB')[1].split('-')[0])
        BANDS = range(number-(10*nblockconcat),number+(10*nblockconcat)+9,10)
        ndppplist = ["{name:s}_SB{b1:03d}-{b2:03d}.{sr:s}.ms".format(name=basename,sr=source,b1=b,b2=b+9) for b in BANDS]
        f=open(ms+'.concatparset','w')
        f.write('msin='+str(ndppplist)+'\nmsin.datacolumn = DATA\nmsin.missingdata=True\nmsin.orderms=False\nmsout='+ ms+'.concat'+'\nsteps=[]\n')
        f.close()
        os.system('rm -rf ' + ms + '.concat')
        cmd = 'NDPPP ' + ms + '.concatparset'
        print cmd
        ndp.run(cmd)
    ndp.wait()

    return


def create_correct_dde_scalarphase_ap_parset(parset, numchanperms):
    f=open(parset, 'w')
    f.write('Strategy.InputColumn = DATA\n')
    f.write('Strategy.ChunkSize   = %s\n' % numpy.int(7200/numchanperms))
    f.write('Strategy.Steps       = [correct]\n')
    f.write('Step.correct.Model.Sources                 = []\n')
    f.write('Step.correct.Model.CommonScalarPhase.Enable= T\n')
    f.write('Step.correct.Model.Cache.Enable            = T\n') 
    f.write('Step.correct.Model.DirectionalGain.Enable  = F\n') 
    f.write('Step.correct.Model.Gain.Enable             = T\n')
    f.write('Step.correct.Model.Phasors.Enable          = F\n')
    f.write('Step.correct.Operation                     = CORRECT\n')
    f.write('Step.correct.Output.Column                 = CORRECTED_DATA\n')
    f.write('Step.correct.Model.Beam.Enable             = F\n') 
    f.write('Step.correct.Output.WriteCovariance        = F\n') 
    f.close()
    return

def create_correct_scalarphasep2_parset(parset, numchanperms):
    f=open(parset, 'w')
    f.write('Strategy.InputColumn = DATA\n')
    f.write('Strategy.ChunkSize   = %s\n' % numpy.int(7200/numchanperms))
    f.write('Strategy.Steps       = [correct]\n')
    f.write('Step.correct.Model.Sources                 = []\n')
    f.write('Step.correct.Model.CommonScalarPhase.Enable= T\n')
    f.write('Step.correct.Model.Cache.Enable            = T\n') 
    f.write('Step.correct.Model.DirectionalGain.Enable  = F\n') 
    f.write('Step.correct.Model.Gain.Enable             = F\n')
    f.write('Step.correct.Model.Phasors.Enable          = F\n')
    f.write('Step.correct.Operation                     = CORRECT\n')
    f.write('Step.correct.Output.Column                 = CORRECTED_DATA_PHASE\n')
    f.write('Step.correct.Model.Beam.Enable             = F\n') 
    f.write('Step.correct.Output.WriteCovariance        = F\n') 
    f.close()
    return

def create_correct_scalarphase_parset(parset, numchanperms):
    f=open(parset, 'w')
    f.write('Strategy.InputColumn = DATA\n')
    f.write('Strategy.ChunkSize   = %s\n' % numpy.int(7200/numchanperms))
    f.write('Strategy.Steps       = [correct]\n')
    f.write('Step.correct.Model.Sources                 = []\n')
    f.write('Step.correct.Model.CommonScalarPhase.Enable= T\n')
    f.write('Step.correct.Model.Cache.Enable            = T\n') 
    f.write('Step.correct.Model.DirectionalGain.Enable  = F\n') 
    f.write('Step.correct.Model.Gain.Enable             = F\n')
    f.write('Step.correct.Model.Phasors.Enable          = F\n')
    f.write('Step.correct.Operation                     = CORRECT\n')
    f.write('Step.correct.Output.Column                 = CORRECTED_DATA\n')
    f.write('Step.correct.Model.Beam.Enable             = F\n') 
    f.write('Step.correct.Output.WriteCovariance        = F\n') 
    f.close()
    return  


#print sys.argv
el=len(sys.argv)

mslist    = sys.argv[1:el-17]
cluster   = str(sys.argv[el-17])
atrous_do = str(sys.argv[el-16])
imsize    = numpy.int(sys.argv[el-15])

nterms                  = numpy.int(sys.argv[el-14])  # only 1 to 3 is supported !!
cellsizetime_a          = numpy.int(sys.argv[el-13])
cellsizetime_p          = numpy.int(sys.argv[el-12])
TECi                    = str(sys.argv[el-11])
clocki                  = str(sys.argv[el-10])
HRi                     = str(sys.argv[el-9])
region                  = str(sys.argv[el-8])
uvrange                 = str(sys.argv[el-7]) # in wavelengths
peelskymodel            = str(sys.argv[el-6])
cellsize                = str(sys.argv[el-5])
basename                = str(sys.argv[el-4])
resolutionstr           = str(sys.argv[el-3])
increaseSNR_via_freqcoverage = str(sys.argv[el-2])
nblockconcat            = numpy.int(sys.argv[el-1])

 
ScaleMaskThreshold = (numpy.float(len(mslist))/5.)**(0.15)

#uvrange = str(120.)

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

casaclean = True


if imsize <=3072:
    wplanes = 1
    FFT = True  # FFT image into MODEL_DATA
    if imsize > 512:
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

if peelskymodel != 'empty':
  FFT = False # use peelskymodel in this case

print 'mslist', mslist
print 'source', cluster
print 'atrous_do', atrous_do 
print 'imsize', imsize








if len(mslist) > 60:
     raise Exception('too many blocks....')


freq_tab     = pt.table(mslist[0]  + '/SPECTRAL_WINDOW')
numchanperms = freq_tab.getcol('NUM_CHAN')[0]
freq_tab.close()

#merge_parmdb = True
phasors      = False   # if true only solve for amps on long timescales
smooth       = True # sometimes almost 0.0 amplitude, causes ripples
phasezero    = True # reset phases from ap calibration
imagenoise   = []




#nblockconcat = 1 # add this number of block on each side (so 1 means concat 3 ms)

if increaseSNR_via_freqcoverage == 'True':
    print 'Going to attach neighboring blocks'
    concat_neighbor_block(basename, resolutionstr, cluster, mslist, nblockconcat)

    print 'Create WEIGHT_SPECTRUM_SC'
    mw=bg(maxp=32)
    for ms in mslist:
        cmd = 'python ' + SCRIPTPATH + '/backup_weights_column.py -i ' + ms + '.concat -o WEIGHT_SPECTRUM_SC'
        mw.run(cmd)
    mw.wait()

    print 'Create WEIGHT_SPECTRUM_AP'
    mw=bg(maxp=32)
    for ms in mslist:
        cmd = 'python ' + SCRIPTPATH + '/backup_weights_column.py -i ' + ms + '.concat -o WEIGHT_SPECTRUM_AP'
        mw.run(cmd)
    mw.wait()


    # scalarphase (fast)
    ionnormfactor = 1.0
    ionscalefactor = 0.25
    print 'Modify WEIGHT_SPECTRUM_SC'
    mw=bg(maxp=numpy.int(16/numchanperms))
    for ms in mslist:
        cmd = 'python ' + SCRIPTPATH + '/modweights_freq.py ' + ms + '.concat ' + \
              str(ionnormfactor) + ' ' + str(ionscalefactor) + ' ' + str(nblockconcat) \
              + ' ' + 'WEIGHT_SPECTRUM_SC' + ' ' + str(numchanperms)
        mw.run(cmd)
    mw.wait()

    # a&p (slow)
    ionnormfactor = 1.1
    ionscalefactor = 0.0
    # "0.7  0.0"  gives [ 0.46893629  1.  0.46893629]
    # "1.0 0.0"  gives [ 0.18901805  0.65936489  1.          0.65936489  0.18901805]
    print 'Modify WEIGHT_SPECTRUM_AP'
    mw=bg(maxp=numpy.int(16/numchanperms))
    for ms in mslist:
        cmd = 'python ' + SCRIPTPATH + '/modweights_freq.py ' + ms + '.concat ' + \
              str(ionnormfactor) + ' ' + str(ionscalefactor) + ' ' + str(nblockconcat) \
              + ' ' + 'WEIGHT_SPECTRUM_AP' + ' ' + str(numchanperms)
        mw.run(cmd)
    mw.wait()  


if increaseSNR_via_freqcoverage == 'True':
    # reset mslist, do seperatly from above so we can restart more easily
    print 'RESET mslist, ".concat" is added'
    for ms_id,ms in enumerate(mslist):
      mslist[ms_id] = ms + '.concat'


msinputlist = ''
for m in mslist:
     msinputlist = msinputlist + ' ' + m



print 'Create the parsets'
if increaseSNR_via_freqcoverage == 'True':
    create_correct_dde_scalarphase_ap_parset('correct_dde_scalarphase_ap.parset', ((2*nblockconcat)+1)*numchanperms)
    create_correct_scalarphasep2_parset('correct_scalarphasep2.parset',((2*nblockconcat)+1)*numchanperms)
    create_correct_scalarphase_parset('correct_scalarphase.parset',((2*nblockconcat)+1)*numchanperms)
else:
    create_correct_dde_scalarphase_ap_parset('correct_dde_scalarphase_ap.parset', numchanperms)
    create_correct_scalarphasep2_parset('correct_scalarphasep2.parset',numchanperms)
    create_correct_scalarphase_parset('correct_scalarphase.parset',numchanperms)



print 'Copying over instrument template'

# add the instrument_template
dummyparmdb = 'instrument_template_Gain_CSphase'
dummyparmdb_csphase = 'instrument_template_CSphase'
for ms in mslist:
    os.system('rm -rf ' + ms +'/' + dummyparmdb)
    os.system('cp -r ' + dummyparmdb + ' ' +  ms + '/instrument_template')
    os.system('rm -rf ' + ms +'/' + dummyparmdb_csphase)
    os.system('cp -r ' + dummyparmdb_csphase + ' ' +  ms + '/instrument_template_csphase')


# reset name for rest of script
dummyparmdb         = 'instrument_template'
dummyparmdb_csphase = 'instrument_template_csphase'



# ------------- peelskymodel------------  
if peelskymodel != 'empty':
  print 'Special case, skymodel given'
  # RUN STEFCAL scalarphase

  parmdb_p = 'instrument_phase0'
  run_stefcal(mslist, peelskymodel,'instrument', 'scalarphase', FFT, \
              uvrange, cellsizetime_p, 'DATA', increaseSNR_via_freqcoverage, numchanperms, nblockconcat)

  # merge the parmdbs fixes NDPPP problem
  for ms in mslist:
    create_merged_parmdb_scalarphase(ms,  ms+'/'+'instrument', ms+'/'+dummyparmdb_csphase,ms+'/'+parmdb_p,cellsizetime_p)

# ottherwise data will be flagged
  for ms in mslist: 
    parmdb_master_outtmp  = ms+"/"+ parmdb_p
    os.system("taql 'update " + parmdb_master_outtmp + " set ENDX=1.e12'")
    os.system("taql 'update " + parmdb_master_outtmp + " set STARTX=1.0'")

  # APPLYCAL, put in 'CORRECTED_DATA_PHASE'
  runbbs(mslist, peelskymodel, 'correct_scalarphasep2.parset', parmdb_p, True, False) 

  # RUN STEFCAL gains
  parmdb_a = 'instrument_amps0'
  run_stefcal(mslist, peelskymodel,parmdb_a, 'diagonal', FFT, uvrange, cellsizetime_a, \
             'CORRECTED_DATA_PHASE', numchanperms, nblockconcat)

  # OUTLIER REMOVAL
  b=bg(maxp=32)
  for ms in mslist:
  # remove outliers from the solutions
    b.run('python ' + SCRIPTPATH + '/smoothcal_a2256_nophasors.py ' + ms + ' ' + ms+'/'+parmdb_a +\
              ' ' + ms+'/'+parmdb_a+'_smoothed')
  b.wait()
  

  parmdb_out = 'instrument_merged'
  # merge the parmdbs fixes NDPPP problem
  for ms in mslist:
    # do not use parmdb_p, but "instrument" due to create_merged_parmdb_scalarphase (important for solint_p ne 1)
    create_merged_parmdb(ms, ms+'/'+parmdb_a+'_smoothed', ms+'/'+'instrument', ms+'/'+dummyparmdb,ms+'/'+parmdb_out,cellsizetime_a,cellsizetime_p)

# ottherwise data will be flagged
  for ms in mslist: 
    parmdb_master_outtmp  = ms+"/"+ parmdb_out
    os.system("taql 'update " + parmdb_master_outtmp + " set ENDX=1.e12'")
    os.system("taql 'update " + parmdb_master_outtmp + " set STARTX=1.0'")

  # APPLYCAL amps
  runbbs(mslist, peelskymodel, 'correct_dde_scalarphase_ap.parset', parmdb_out, True, False)
 
  ### MAKE IMAGE 0 ###
  if casaclean:
    imout,mask = make_image(mslist, cluster, '0', 6.*ScaleMaskThreshold, 6.*ScaleMaskThreshold, nterms, \
                            atrous_do, imsize, 'CORRECTED_DATA', uvrange, increaseSNR_via_freqcoverage, \
                            numchanperms, nblockconcat)

    if nterms < 2:    
      print 'IMAGE STATISTICS ', find_imagenoise(imout + '.image')
      rms, dynamicrange =  find_imagenoise(imout + '.image')
    else:
      print 'IMAGE STATISTICS ', find_imagenoise(imout + '.image.tt0')
      rms, dynamicrange =  find_imagenoise(imout + '.image.tt0')
  else: # wsclean
    print 'Not ready'
    print 'IMAGE STATISTICS ', find_imagenoise(imout + '.image')

  plotds9image(mslist, cluster)
  # GET OUT OF SELFCAL script
  sys.exit()   
# ------------- peelskymodel------------  
   



#### MAKE IMAGE 0 ###
imout,mask = make_image(mslist, cluster, '0', 6, 5, nterms, atrous_do, \
                        imsize, 'DATA', uvrange, increaseSNR_via_freqcoverage, \
			numchanperms, nblockconcat)
#####################


#imout = 'im0_clusters1'
#mask  = 'im0_clusters1.mask'

if nterms < 2: 
   print 'IMAGE STATISTICS ', find_imagenoise(imout + '.image')
   rms, dynamicrange =  find_imagenoise(imout + '.image')
else:
   print 'IMAGE STATISTICS ', find_imagenoise(imout + '.image.tt0')
   rms, dynamicrange =  find_imagenoise(imout + '.image.tt0')

plotds9image(mslist, cluster)

### CALIBRATE WITH BBS PHASE ONLY 1 ###
# create skymodel for BBS
os.system(SCRIPTPATH + '/casapy2bbs.py -m '+ mask + ' ' +'-t ' + str(nterms)+ ' ' + imout+'.model ' +  imout+'.skymodel')
if FFT:
 os.system('/home/rvweeren/software/casapy/casapy-42.2.30986-1-64b/casapy --nogui -c ' + SCRIPTPATH + '/ft_v2.py ' + msinputlist + ' ' + imout+'.model' \
           + ' ' + str(nterms) + ' '+ str(wplanes))

# phase only calibrate
skymodel = imout+'.skymodel'

# RUN STEFCAL scalarphase
run_stefcal(mslist, skymodel, 'instrument', 'scalarphase', FFT, uvrange, cellsizetime_p, \
            'DATA', increaseSNR_via_freqcoverage, numchanperms, nblockconcat)

parmdb_out = 'instrument_phase'
# merge the parmdbs fixes NDPPP problem
for ms in mslist:
  create_merged_parmdb_scalarphase(ms,  ms+'/'+'instrument', ms+'/'+dummyparmdb_csphase,ms+'/'+parmdb_out,cellsizetime_p)

# ottherwise data will be flagged
for ms in mslist: 
    parmdb_master_outtmp  = ms+"/"+ parmdb_out
    os.system("taql 'update " + parmdb_master_outtmp + " set ENDX=1.e12'")
    os.system("taql 'update " + parmdb_master_outtmp + " set STARTX=1.0'")
 

# APPLYCAL
runbbs(mslist, skymodel, 'correct_scalarphase.parset', 'instrument_phase', True, False) 


### MAKE IMAGE 1 ###
imout,mask = make_image(mslist, cluster, '1', 8.*ScaleMaskThreshold, 7.*ScaleMaskThreshold, nterms, 
			atrous_do, imsize, 'CORRECTED_DATA', uvrange, increaseSNR_via_freqcoverage, \
			numchanperms, nblockconcat)
####################

if nterms < 2:    
   print 'IMAGE STATISTICS ', find_imagenoise(imout + '.image')
   rms, dynamicrange =  find_imagenoise(imout + '.image')
else:
   print 'IMAGE STATISTICS ', find_imagenoise(imout + '.image.tt0')
   rms, dynamicrange =  find_imagenoise(imout + '.image.tt0')

plotds9image(mslist, cluster) 

### CALIBRATE WITH BBS PHASE ONLY 2 ###
# create skymodel for BBS
os.system(SCRIPTPATH + '/casapy2bbs.py -m '+ mask + ' ' +'-t ' + str(nterms)+ ' ' + imout+'.model ' +  imout+'.skymodel')
if FFT:
  os.system('/home/rvweeren/software/casapy/casapy-42.2.30986-1-64b/casapy --nogui -c ' + SCRIPTPATH + '/ft_v2.py ' + msinputlist + ' ' + imout+'.model' \
            + ' ' + str(nterms) + ' '+ str(wplanes))
 
# phase only calibrate
skymodel = imout+'.skymodel'



# RUN STEFCAL scalarphase
run_stefcal(mslist, skymodel, 'instrument', 'scalarphase', FFT, uvrange, cellsizetime_p, \
            'DATA', increaseSNR_via_freqcoverage, numchanperms, nblockconcat)

parmdb_out = 'instrument_phase'
# merge the parmdbs fixes NDPPP problem
for ms in mslist:
  create_merged_parmdb_scalarphase(ms,  ms+'/'+'instrument', ms+'/'+dummyparmdb_csphase,ms+'/'+parmdb_out,cellsizetime_p)

# ottherwise data will be flagged
for ms in mslist: 
    parmdb_master_outtmp  = ms+"/"+ parmdb_out
    os.system("taql 'update " + parmdb_master_outtmp + " set ENDX=1.e12'")
    os.system("taql 'update " + parmdb_master_outtmp + " set STARTX=1.0'")

# APPLYCAL
runbbs(mslist, skymodel, 'correct_scalarphase.parset', 'instrument_phase', True, False) 

### MAKE IMAGE 2 ###
imout,mask = make_image(mslist, cluster, '2', 8.*ScaleMaskThreshold, 7.*ScaleMaskThreshold, nterms, \
                        atrous_do, imsize, 'CORRECTED_DATA', uvrange, increaseSNR_via_freqcoverage, \
                        numchanperms, nblockconcat)
####################

if nterms < 2:    
   print 'IMAGE STATISTICS ', find_imagenoise(imout + '.image')
   rms, dynamicrange =  find_imagenoise(imout + '.image')
else:
   print 'IMAGE STATISTICS ', find_imagenoise(imout + '.image.tt0')
   rms, dynamicrange =  find_imagenoise(imout + '.image.tt0')

plotds9image(mslist, cluster)

### CALIBRATE WITH BBS PHASE+AMP 1 ###
os.system(SCRIPTPATH +'/casapy2bbs.py -m '+ mask + ' ' +'-t ' + str(nterms)+ ' ' + imout+'.model ' +  imout+'.skymodel')
if FFT:
  os.system('/home/rvweeren/software/casapy/casapy-42.2.30986-1-64b/casapy --nogui -c ' + SCRIPTPATH +'/ft_v2.py ' + msinputlist + ' ' + imout+'.model' \
            + ' ' + str(nterms) + ' '+ str(wplanes)) 

skymodel = imout+'.skymodel'

# RUN STEFCAL scalarphase
parmdb_p = 'instrument_phase0'
run_stefcal(mslist, skymodel,'instrument', 'scalarphase', FFT, uvrange, cellsizetime_p, \
            'DATA', increaseSNR_via_freqcoverage, numchanperms, nblockconcat)

# merge the parmdbs fixes NDPPP problem
for ms in mslist:
  create_merged_parmdb_scalarphase(ms,  ms+'/'+'instrument', ms+'/'+dummyparmdb_csphase,ms+'/'+parmdb_p,cellsizetime_p)

# ottherwise data will be flagged
for ms in mslist: 
    parmdb_master_outtmp  = ms+"/"+ parmdb_p
    os.system("taql 'update " + parmdb_master_outtmp + " set ENDX=1.e12'")
    os.system("taql 'update " + parmdb_master_outtmp + " set STARTX=1.0'")

# APPLYCAL, put in 'CORRECTED_DATA_PHASE'
runbbs(mslist, skymodel, 'correct_scalarphasep2.parset', parmdb_p, True, False) 

# RUN STEFCAL gains
parmdb_a = 'instrument_amps0'
run_stefcal(mslist, skymodel,parmdb_a, 'diagonal', FFT, uvrange, cellsizetime_a, \
            'CORRECTED_DATA_PHASE', increaseSNR_via_freqcoverage, numchanperms, nblockconcat)

# OUTLIER REMOVAL
b=bg(maxp=32)
for ms in mslist:
    # remove outliers from the solutions
    b.run('python ' + SCRIPTPATH + '/smoothcal_a2256_nophasors.py ' + ms + ' ' + \
          ms+'/'+parmdb_a + ' ' + ms+'/'+parmdb_a+'_smoothed')
b.wait()


parmdb_out = 'instrument_merged'
# merge the parmdbs fixes NDPPP problem
for ms in mslist:
  # do not use parmdb_p, but "instrument" due to create_merged_parmdb_scalarphase (important for solint_p ne 1)
  create_merged_parmdb(ms, ms+'/'+parmdb_a+'_smoothed', ms+'/'+'instrument', ms+'/'+dummyparmdb,ms+'/'+parmdb_out,cellsizetime_a,cellsizetime_p)

# ottherwise data will be flagged
for ms in mslist: 
    parmdb_master_outtmp  = ms+"/"+ parmdb_out
    os.system("taql 'update " + parmdb_master_outtmp + " set ENDX=1.e12'")
    os.system("taql 'update " + parmdb_master_outtmp + " set STARTX=1.0'")


# APPLYCAL amps
runbbs(mslist, skymodel, 'correct_dde_scalarphase_ap.parset', parmdb_out, True, False)

   
### MAKE IMAGE 3 ###
# use a somewhat higher island threshold here to prevent artifacts entering the model 
imout,mask = make_image(mslist, cluster, '3', 10.*ScaleMaskThreshold, 10.*ScaleMaskThreshold, nterms, \
                        atrous_do, imsize, 'CORRECTED_DATA', uvrange, increaseSNR_via_freqcoverage, \
                        numchanperms, nblockconcat)

#imout = 'im3_clusters1'
#mask  = 'im3_clusters1.mask

if nterms < 2:    
   print 'IMAGE STATISTICS ', find_imagenoise(imout + '.image')
   rms, dynamicrange =  find_imagenoise(imout + '.image')
else:
   print 'IMAGE STATISTICS ', find_imagenoise(imout + '.image.tt0')
   rms, dynamicrange =  find_imagenoise(imout + '.image.tt0')

plotds9image(mslist, cluster)


number_forced_selfcalcycles = 7
rms_old          = 1.e9 # bad values to start with
dynamicrange_old = 1. # low value to start with so we get into the while loop
factor           = 1.0125 # demand 1.25% improvement
im_count         = 4

#### GO INTO THE WHILE LOOP
#while(((dynamicrange/factor) > dynamicrange_old) or ((rms*factor) < rms_old)):
while ((dynamicrange/factor) > dynamicrange_old) :

  #### CALIBRATE  BBS PHASE+AMP ###
  # make model
  os.system(SCRIPTPATH +'/casapy2bbs.py -m '+ mask + ' ' +'-t ' + str(nterms)+ ' ' + imout+'.model ' +  imout+'.skymodel')
  if FFT:
    os.system('/home/rvweeren/software/casapy/casapy-42.2.30986-1-64b/casapy --nogui -c '+ SCRIPTPATH +'/ft_v2.py ' + msinputlist + ' ' + imout+'.model' \
              + ' ' + str(nterms) + ' '+ str(wplanes))
	    
  #parmdb keep from previous step
  skymodel = imout+'.skymodel'


  # RUN STEFCAL scalarphase
  parmdb_p = 'instrument_phase1'
  run_stefcal(mslist, skymodel, 'instrument', 'scalarphase', FFT, uvrange, cellsizetime_p, 'DATA', \
              increaseSNR_via_freqcoverage, numchanperms, nblockconcat)


  # merge the parmdbs fixes NDPPP problem
  for ms in mslist:
    create_merged_parmdb_scalarphase(ms,  ms+'/'+'instrument', ms+'/'+dummyparmdb_csphase,ms+'/'+parmdb_p,cellsizetime_p)

  # ottherwise data will be flagged
  for ms in mslist: 
    parmdb_master_outtmp  = ms+"/"+ parmdb_p
    os.system("taql 'update " + parmdb_master_outtmp + " set ENDX=1.e12'")
    os.system("taql 'update " + parmdb_master_outtmp + " set STARTX=1.0'")

  # APPLYCAL, put in 'CORRECTED_DATA_PHASE'
  runbbs(mslist, skymodel, 'correct_scalarphasep2.parset', parmdb_p, True, False) 

  # RUN STEFCAL gains
  parmdb_a = 'instrument_amps1'
  run_stefcal(mslist, skymodel,parmdb_a, 'diagonal', FFT, uvrange, cellsizetime_a, \
              'CORRECTED_DATA_PHASE', increaseSNR_via_freqcoverage, numchanperms, nblockconcat)

  # OUTLIER REMOVAL
  b=bg(maxp=32)
  for ms in mslist:
      # remove outliers from the solutions
      b.run('python ' + SCRIPTPATH + '/smoothcal_a2256_nophasors.py ' + ms + ' ' + \
            ms+'/'+parmdb_a + ' ' + ms+'/'+parmdb_a+'_smoothed')
  b.wait()


  parmdb_out = 'instrument_merged'
  # merge the parmdbs fixes NDPPP problem
  for ms in mslist:
    create_merged_parmdb(ms, ms+'/'+parmdb_a+'_smoothed', ms+'/'+'instrument', ms+'/'+dummyparmdb,ms+'/'+parmdb_out,cellsizetime_a,cellsizetime_p)

# ottherwise data will be flagged
  for ms in mslist: 
    parmdb_master_outtmp  = ms+"/"+ parmdb_out
    os.system("taql 'update " + parmdb_master_outtmp + " set ENDX=1.e12'")
    os.system("taql 'update " + parmdb_master_outtmp + " set STARTX=1.0'")

  # APPLYCAL amps
  runbbs(mslist, skymodel, 'correct_dde_scalarphase_ap.parset', parmdb_out, True, False)
   
  ### MAKE IMAGE #N ###
  imout,mask = make_image(mslist, cluster, str(im_count), 8.*ScaleMaskThreshold, 8.*ScaleMaskThreshold, nterms, \
                          atrous_do, imsize, 'CORRECTED_DATA', uvrange, increaseSNR_via_freqcoverage,\
                          numchanperms, nblockconcat)
  im_count = im_count + 1

  # save previous values to compare with
  rms_old          = rms
  dynamicrange_old = dynamicrange 
  if nterms < 2:    
      print 'IMAGE STATISTICS ', find_imagenoise(imout + '.image')
      rms, dynamicrange =  find_imagenoise(imout + '.image')
  else:
     print 'IMAGE STATISTICS ', find_imagenoise(imout + '.image.tt0')
     rms, dynamicrange =  find_imagenoise(imout + '.image.tt0')

  if im_count < number_forced_selfcalcycles:
    rms_old          = 1.e9 # bad values to start with
    dynamicrange_old = 1.
    
  plotds9image(mslist, cluster) 


##################################################################
##################################################################
##################################################################
##################################################################

