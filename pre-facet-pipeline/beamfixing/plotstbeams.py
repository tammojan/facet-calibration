#!/usr/bin/python

"""
Written 6 Feb 2014 by G Heald
Looks into LOFAR MSs, gets missing tile info, plots predicted array
beams, and compares with the perfect (no missing tile) situation
"""

import argparse
import pylab as plt
import pyrap.tables as pt
import stationBeam2 as stbeam
import numpy

def main(args):
	plt.figure(figsize=(11.69,8.27))
	times = pt.table(args.MSname).sort('unique desc TIME').getcol('TIME')
	#print times
	time = times[len(times)/2]
	#print time
	freq = pt.table(args.MSname+'/SPECTRAL_WINDOW').getcol('REF_FREQUENCY')[0]
	direction = pt.table(args.MSname+'/FIELD').getcol('PHASE_DIR')[0][0]
	print direction
	nx = 60
	ny = 60
	fieldrad = 30./180.*numpy.pi
	dx = numpy.linspace(direction[0]-fieldrad,direction[0]+fieldrad,nx)
	dy = numpy.linspace(direction[1]-fieldrad,direction[1]+fieldrad,ny)
	XX,YY=numpy.meshgrid(dx,dy)
	XXf = XX.flatten()
	YYf = YY.flatten()
	thisbeamxy = numpy.zeros((nx,ny),dtype=numpy.complex)
	goodbeamxy = numpy.zeros((nx,ny),dtype=numpy.complex)
	laf = pt.table(args.MSname+'/LOFAR_ANTENNA_FIELD')
	ant = pt.table(args.MSname+'/ANTENNA')
	sn = ant.getcol('NAME')


	if args.stations == '__ALL__':
		sttoplot = range(len(sn))
	else:
		sttoplot = where(sn==args.stations)
	myst = stbeam.multistationBeam(ms=args.MSname+'/',times=times[:1],direction=zip(XXf,YYf))
	for st in sttoplot:
		if 'HBA0' in sn[st]:
			ind = range(24)
			goodflags = [False]*24+[True]*24
		elif 'HBA1' in sn[st]:
			ind = range(24,48)
			goodflags = [True]*24+[False]*24
		else:
			ind = range(48)
			goodflags = [ True, True, True ,True ,True,False,False, True ,True, True,False,False,False,False ,True ,True ,True,False,False,False,False,False,False ,True ,True,False,False,False,False,False,False ,True, True, True,False,False,False,False, True, True, True,False,False ,True ,True ,True, True, True] 
		offsets = laf.getcell('ELEMENT_OFFSET',st)
		flags = laf.getcell('ELEMENT_FLAG',st)
		print flags[ind,0].astype(int)
		#plt.subplot(2,2,1)
		plt.clf()
		plt.subplot(2,2,1)
		plt.scatter(offsets[ind,0],offsets[ind,1],marker='s',c=1-(flags[ind,0].astype(int)),cmap='RdYlGn',s=225)
		plt.title(sn[st]+', RED=FLAGGED')
		plt.subplot(2,2,2)
	
		count = 0
		for i in range(nx):
			for j in range(ny):
				elements=[]
				for ist in range(len(sn)):
					elements.append(myst.get_element_phases(ist,numpy.array([freq,]),count))
				fx=lambda weights,station_number,polarization: numpy.sum(weights[numpy.newaxis,:]*elements[station_number][:,:,polarization],axis=1)/numpy.sum(weights)
				stbeamflags = myst.flags[:]

				thisbeamxy[i,j] = fx(1-stbeamflags[st,:,0],st,0)[0]
				goodbeamxy[i,j] = fx(1-numpy.array(goodflags),st,0)[0]
				count += 1
		#print thisbeamxy
		# left right bottom top
		plt.imshow(numpy.abs(thisbeamxy),extent=[dx[0],dx[-1],dy[0],dy[-1]],vmin=0.,vmax=1.)
		plt.colorbar()
		plt.title('Beam with MS element flags')
		plt.subplot(2,2,3)
		plt.imshow(numpy.abs(goodbeamxy),extent=[dx[0],dx[-1],dy[0],dy[-1]],vmin=0.,vmax=1.)
		plt.colorbar()
		plt.title('Beam with no flags')
		plt.subplot(2,2,4)
		#plt.imshow((numpy.abs(thisbeamxy)-numpy.abs(goodbeamxy))/numpy.abs(goodbeamxy),extent=[dx[0],dx[-1],dy[0],dy[-1]],vmin=0.,vmax=1.)
		plt.imshow(numpy.abs(thisbeamxy)-numpy.abs(goodbeamxy),extent=[dx[0],dx[-1],dy[0],dy[-1]],vmin=-0.25,vmax=0.25)
		plt.colorbar()
		plt.title('Beam difference (MS-unflagged)')
		#plt.show()
		print sn[st]
		plt.savefig(args.outname+'-'+sn[st]+'.'+args.format)
	

p = argparse.ArgumentParser()
p.add_argument("MSname",help="Name of Measurement Set")
p.add_argument('-s','--stations',help='Comma separated list of stations to plot [default all]',default='__ALL__')
p.add_argument('-o','--outname',help='Output plot basename [default stplot]',default='stplot')
p.add_argument('-f','--format',help='File format for plots [default png]',default='png')
args = p.parse_args()
main(args)

