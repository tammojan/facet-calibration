import numpy
import os
import sys
from coordinates_mode import *
import pyrap.images

pi = numpy.pi


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



def load_bbs_skymodel(infilename):
        tmp_input = infilename + '.tmp'
	
	# remove empty lines sed '/^$/d'
	# remove format line grep -v 'format'
	# remove comment lines  grep -v '#'
        os.system("grep -v '00:00:00, +00.00.00' " + infilename+ " | grep -v '#' | grep -v '00:00:00, +90.00.00' | grep -v 'format' | sed '/^$/d'>" + tmp_input) # to remove patches headers from skymodel

        types = numpy.dtype({'names':['Name', 'Type','Patch','Ra', 'Dec', 'I', 'Q', 'U', 'V', 'Maj', 'Min', 'PA', 'RefFreq', 'Spidx'],\
	       'formats':['S100','S100','S100','S100','S100',numpy.float,numpy.float,numpy.float,numpy.float,numpy.float,numpy.float,numpy.float,numpy.float,'S100']})
        

        data  = numpy.loadtxt(tmp_input, comments='format', unpack=True, delimiter=', ', dtype=types)
	os.system('rm ' + tmp_input)
	
        return data	

def compute_patch_center(data,fluxweight):
   #print 'These are the input patches', numpy.unique(data['Patch'])
   patches = numpy.unique(data['Patch'])

   ra_patches   = numpy.zeros(len(patches))
   dec_patches  = numpy.zeros(len(patches))
   flux_patches = numpy.zeros(len(patches))

   for (patch_id,patch) in enumerate(patches):
    idx = numpy.where(data['Patch'] == patch)
    #print 'Patch', patch, 'has', len(idx[0]), 'components'
    #print numpy.shape( data['Patch'][idx])
  
    # set to zero
    ra_patch   = 0.
    dec_patch  = 0.
    ra_weights = 0.
    dec_weights= 0.
    flux_patch = 0.
    
    for component in idx[0]:
    
      # conver RA, DEC to degrees for component
      ra_comp  = data['Ra'][component]
      dec_comp = data['Dec'][component]
      ra_comp  = (ra_comp.split(':'))
      dec_comp = (dec_comp.split('.'))
      flux_comp= numpy.float(data['I'][component])
      ra_comp  = hmstora(numpy.float(ra_comp[0]),numpy.float(ra_comp[1]),numpy.float(ra_comp[2]))
 
      if '-' in dec_comp[0]: sign = '-'
      else: sign = '+'
      if len(dec_comp) == 4:  # decimal arcsec in Dec
        dec_comp = dmstodec(numpy.float(dec_comp[0]),numpy.float(dec_comp[1]),numpy.float(str(dec_comp[2]+"."+dec_comp[3])),sign=sign) 
      else:
        dec_comp = dmstodec(numpy.float(dec_comp[0]),numpy.float(dec_comp[1]),numpy.float(dec_comp[2]),sign=sign)
       
      # calculate the average weighted patch center, and patch flux 
      flux_patch = flux_patch + flux_comp
      #print ra_comp, dec_comp, flux_comp
      
      if fluxweight:
       ra_patch = ra_patch + (flux_comp*ra_comp)
       dec_patch= dec_patch+ (flux_comp*dec_comp)
       ra_weights = ra_weights + flux_comp
       dec_weights= dec_weights + flux_comp
      else:
       ra_patch = ra_patch + (1.*ra_comp)
       dec_patch= dec_patch+ (1.*dec_comp)
       ra_weights = ra_weights + 1.
       dec_weights= dec_weights + 1.
    
    #print 'Center RA, Center DEC, flux', ra_patch/ra_weights, dec_patch/dec_weights,  flux_patch
    
    
    ra_patches[patch_id]  = ra_patch/ra_weights
    dec_patches[patch_id] = dec_patch/dec_weights
    
    
    
    flux_patches[patch_id]= flux_patch
    
   return patches, ra_patches,dec_patches, flux_patches

def compute_patch_center_libsproblem(data,fluxweight):

   ### FIX FOR USE PYTHONLIBS
   fixid_patch = 2
   fixid_ra    = 3
   fixid_dec   = 4
   fixid_I     = 5   
   
   #print 'These are the input patches',numpy.unique(data[fixid_patch]) 
   patches = numpy.unique(data[fixid_patch])
   ra_patches   = numpy.zeros(len(patches))
   dec_patches  = numpy.zeros(len(patches))
   flux_patches = numpy.zeros(len(patches))

   for (patch_id,patch) in enumerate(patches):
    idx = numpy.where(data[fixid_patch] == patch)
    #print 'Patch', patch, 'has', len(idx[0]), 'components'
    #print numpy.shape( data['Patch'][idx])
  
    # set to zero
    ra_patch   = 0.
    dec_patch  = 0.
    ra_weights = 0.
    dec_weights= 0.
    flux_patch = 0.
    
    for component in idx[0]:
    
      # conver RA, DEC to degrees for component
      ra_comp  = data[fixid_ra][component]
      dec_comp = data[fixid_dec][component]
      ra_comp  = (ra_comp.split(':'))
      dec_comp = (dec_comp.split('.'))
      flux_comp= numpy.float(data[fixid_I][component])
      ra_comp  = hmstora(numpy.float(ra_comp[0]),numpy.float(ra_comp[1]),numpy.float(ra_comp[2]))

      if '-' in dec_comp[0]: sign = '-'
      else: sign = '+' 
      if len(dec_comp) == 4:  # decimal arcsec in Dec
        dec_comp = dmstodec(numpy.float(dec_comp[0]),numpy.float(dec_comp[1]),numpy.float(str(dec_comp[2]+"."+dec_comp[3])),sign=sign) 
      else:
        dec_comp = dmstodec(numpy.float(dec_comp[0]),numpy.float(dec_comp[1]),numpy.float(dec_comp[2]),sign=sign)
       
      # calculate the average weighted patch center, and patch flux 
      flux_patch = flux_patch + flux_comp
      #print ra_comp, dec_comp, flux_comp
      
      if fluxweight:
       ra_patch = ra_patch + (flux_comp*ra_comp)
       dec_patch= dec_patch+ (flux_comp*dec_comp)
       ra_weights = ra_weights + flux_comp
       dec_weights= dec_weights + flux_comp
      else:
       ra_patch = ra_patch + (1.*ra_comp)
       dec_patch= dec_patch+ (1.*dec_comp)
       ra_weights = ra_weights + 1.
       dec_weights= dec_weights + 1.
    
    #print 'Center RA, Center DEC, flux', ra_patch/ra_weights, dec_patch/dec_weights,  flux_patch
    
    
    ra_patches[patch_id]  = ra_patch/ra_weights
    dec_patches[patch_id] = dec_patch/dec_weights
    
    
    
    flux_patches[patch_id]= flux_patch
    
   return patches, ra_patches,dec_patches, flux_patches



#def cal_return_slist(imagename,skymodel, direction, imsize):

#print sys.argv


 
imagename   = str(sys.argv[1])
skymodel    = str(sys.argv[2])
ref_sourcet = str(sys.argv[3])
ref_source = ref_sourcet .split(',')


  
fluxweight = False

data = load_bbs_skymodel(skymodel) 

if len(numpy.shape(data)) == 1:  # in this case not issue and we do not use Pythonlibs
   patchest,ra_patches,dec_patches, flux_patches =  compute_patch_center(data,fluxweight)
   #print 'option 1'
if len(numpy.shape(data)) == 2: 
   patchest,ra_patches,dec_patches, flux_patches = compute_patch_center_libsproblem(data,fluxweight)
   #print 'option 2'

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

print sourcess




