#!/usr/bin/python

from pylab import *
import pyrap.tables as tab 
from pyrap.measures import measures
import pyrap.quanta as qa;
me=measures()
light_speed=299792458.
class stationBeam():

    def __init__(self,ms,times=[],direction=[]):
        myt=tab.table(ms+'/LOFAR_ANTENNA_FIELD')
        self.flags=myt.getcol("ELEMENT_FLAG")
        self.offsets=myt.getcol("ELEMENT_OFFSET")
        myt=tab.table(ms+'/ANTENNA')
        self.fieldcenter=myt.getcol("LOFAR_PHASE_REFERENCE")
        self.stationcenter=myt.getcol("POSITION")
        self.offsetshift=self.fieldcenter-self.stationcenter
        self.ms=ms
        self.refdir=tab.table(ms+'FIELD').getcol('PHASE_DIR')[0][0]
        itrfdir=[]
        self.times=times
        self.direction=direction
        self.initiated=False
        self.init_times()

    def init_times(self):
        if not list(self.direction):
            print "cannot init, set direction first"
        if not list(self.times):
            print "cannot init, set times first"
        self.delayxyz=[]
        self.phasexyz=[]
        self.k=[]
        for ist in range(self.stationcenter.shape[0]):
            print "getting station",ist,"from",self.stationcenter.shape[0]
            self.delayxyz.append([])
            self.phasexyz.append([])
            self.k.append([])
            pos=self.stationcenter[ist]
            p = me.position('ITRF',str(pos[0])+'m',str(pos[1])+'m',str(pos[2])+'m')
            me.do_frame(p);
            ra=str(self.refdir[0])+'rad'
            dec=str(self.refdir[1])+'rad'
            phasedir=me.direction('J2000',ra,dec)
            ra=str(self.direction[0])+'rad'
            dec=str(self.direction[1])+'rad'
            delaydir=me.direction('J2000',ra,dec)
            for itm,time in enumerate(self.times):
                t=me.epoch("UTC",qa.quantity(str(time)+'s'))
                me.do_frame(t)
                itrfdir = me.measure(phasedir,'ITRF')
                itrfdelaydir = me.measure(delaydir,'ITRF')
                coslon=np.cos(itrfdir['m0']['value'])
                coslat=np.cos(itrfdir['m1']['value'])
                sinlon=np.sin(itrfdir['m0']['value'])
                sinlat=np.sin(itrfdir['m1']['value'])
                phasexyz=np.array([coslat*coslon,coslat*sinlon,sinlat])
                coslon=np.cos(itrfdelaydir['m0']['value'])
                coslat=np.cos(itrfdelaydir['m1']['value'])
                sinlon=np.sin(itrfdelaydir['m0']['value'])
                sinlat=np.sin(itrfdelaydir['m1']['value'])
                delayxyz=np.array([coslat*coslon,coslat*sinlon,sinlat])
                k=delayxyz-phasexyz
                self.delayxyz[-1].append(delayxyz)
                self.phasexyz[-1].append(phasexyz)
                self.k[-1].append(k)
        self.delayxyz=np.array(self.delayxyz)
        self.phasexyz=np.array(self.phasexyz)
        self.k=np.array(self.k)
        self.initiated=True

    def settimes(self,times):
        self.times=times
        if list(self.direction):
            self.init_times()
    def setdirection(self,direction):
        self.direction=direction
        if list(self.times):
            self.init_times()

    def getAF(self,ist,freqs):
        AF=np.zeros((2,len(self.times),len(freqs)),dtype=np.complex64)
        if not self.initiated:
            print "initialize first with times and direction"
            return AF
        
        elements=self.offsets[ist]
        flags=self.flags[ist]
        phase=2*np.pi*freqs*np.sum(self.k[ist,:,np.newaxis,:]*(elements+self.offsetshift[ist])[np.newaxis,:,:],axis=2)[:,:,np.newaxis]/light_speed
        cdata=np.cos(phase)+1j*np.sin(phase)
        AF[0,:]=np.sum(cdata[:,flags[:,0]==0],axis=1)/np.sum(np.logical_not(flags[:,0]))
        AF[1,:]=np.sum(cdata[:,flags[:,1]==0],axis=1)/np.sum(np.logical_not(flags[:,1]))
            
        return AF


    def get_element_phases(self,ist,freqs):
        if not self.initiated:
            print "initialize first with times and direction"
            AF=np.zeros((len(self.times),self.flags.shape[0],len(freqs)),dtype=np.complex64)
            return AF
        
        elements=self.offsets[ist]
        
        phase=2*np.pi*freqs*np.sum(self.k[ist,:,np.newaxis,:]*(elements+self.offsetshift[ist])[np.newaxis,:,:],axis=2)[:,:,np.newaxis]/light_speed
        cdata=np.cos(phase)+1j*np.sin(phase)
        return cdata


class multistationBeam(stationBeam):
    #def __init__(self,ms,times=[],direction=[]):
    #    super(multistationBeam,self).__init__(ms,times,direction)



    def init_times(self):
        if not list(self.direction):
            print "cannot init, set direction first"
        if not list(self.times):
            print "cannot init, set times first"
        self.delayxyz=[[] for i in self.direction]
        self.phasexyz=[[] for i in self.direction]
        self.k=[[] for i in self.direction]
        for ist in range(self.stationcenter.shape[0]):
            print "getting station",ist,"from",self.stationcenter.shape[0]
            for i in range(len(self.direction)):
                self.delayxyz[i].append([])
                self.phasexyz[i].append([])
                self.k[i].append([])
            pos=self.stationcenter[ist]
            p = me.position('ITRF',str(pos[0])+'m',str(pos[1])+'m',str(pos[2])+'m')
            me.do_frame(p);
            ra=str(self.refdir[0])+'rad'
            dec=str(self.refdir[1])+'rad'
            phasedir=me.direction('J2000',ra,dec)
            delaydir=[]
            for idir,dir in enumerate(self.direction):
                ra=str(dir[0])+'rad'
                dec=str(dir[1])+'rad'
                delaydir.append(me.direction('J2000',ra,dec))
            for itm,time in enumerate(self.times):
                t=me.epoch("UTC",qa.quantity(str(time)+'s'))
                me.do_frame(t)
                itrfdir = me.measure(phasedir,'ITRF')
                coslon=np.cos(itrfdir['m0']['value'])
                coslat=np.cos(itrfdir['m1']['value'])
                sinlon=np.sin(itrfdir['m0']['value'])
                sinlat=np.sin(itrfdir['m1']['value'])
                phasexyz=np.array([coslat*coslon,coslat*sinlon,sinlat])
                for idir,dir in enumerate(self.direction):
                    itrfdelaydir = me.measure(delaydir[idir],'ITRF')
                    coslon=np.cos(itrfdelaydir['m0']['value'])
                    coslat=np.cos(itrfdelaydir['m1']['value'])
                    sinlon=np.sin(itrfdelaydir['m0']['value'])
                    sinlat=np.sin(itrfdelaydir['m1']['value'])
                    delayxyz=np.array([coslat*coslon,coslat*sinlon,sinlat])
                    k=delayxyz-phasexyz
                    self.delayxyz[idir][-1].append(delayxyz)
                    self.phasexyz[idir][-1].append(phasexyz)
                    self.k[idir][-1].append(k)
        self.delayxyz=np.array(self.delayxyz)
        self.phasexyz=np.array(self.phasexyz)
        self.k=np.array(self.k)
        self.initiated=True

  
    def getAF(self,ist,freqs,isrc):
        AF=np.zeros((2,len(self.times),len(freqs)),dtype=np.complex64)
        if not self.initiated:
            print "initialize first with times and direction"
            return AF
        
        elements=self.offsets[ist]
        flags=self.flags[ist]
        phase=2*np.pi*freqs*np.sum(self.k[isrc,ist,:,np.newaxis,:]*(elements+self.offsetshift[ist])[np.newaxis,:,:],axis=2)[:,:,np.newaxis]/light_speed
        cdata=np.cos(phase)+1j*np.sin(phase)
        AF[0,:]=np.sum(cdata[:,flags[:,0]==0],axis=1)/np.sum(np.logical_not(flags[:,0]))
        AF[1,:]=np.sum(cdata[:,flags[:,1]==0],axis=1)/np.sum(np.logical_not(flags[:,1]))
            
        return AF


    def get_element_phases(self,ist,freqs,isrc):
        if not self.initiated:
            print "initialize first with times and direction"
            AF=np.zeros((len(self.times),self.flags.shape[0],len(freqs)),dtype=np.complex64)
            return AF
        
        elements=self.offsets[ist]
        
        phase=2*np.pi*freqs*np.sum(self.k[isrc,ist,:,np.newaxis,:]*(elements+self.offsetshift[ist])[np.newaxis,:,:],axis=2)[:,:,np.newaxis]/light_speed
        cdata=np.cos(phase)+1j*np.sin(phase)
        return cdata
