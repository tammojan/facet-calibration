import numpy
import os
import sys
from coordinates_mode import *
import pyrap.images

pi = numpy.pi

def make_image(ra,dec,data,cdelt):

    # make a dummy image to get a starting set of co-ords
    # avoids using casacore co-ords directly

    im=pyrap.images.image('',shape=data.shape)
    cs=im.coordinates()
    cs.dict()['direction0']['units']=['rad','rad']
    cs.dict()['direction0']['crval']=[ra,dec]
    cs.dict()['direction0']['cdelt']=[-cdelt,cdelt]
    im=pyrap.images.image('',values=data,coordsys=cs)
    return im

def load_bbs_skymodel(infilename):
    tmp_input = infilename + '.tmp'

    # remove empty lines sed '/^$/d'
    # remove format line grep -v 'format'
    # remove comment lines  grep -v '#'
    os.system("grep -v '00:00:00, +00.00.00' " + infilename+ " | grep -v '#' |  grep -v '00:00:00, +90.00.00' | grep -v 'format' | sed '/^$/d'>" + tmp_input) # to remove patches headers from skymodel

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
                dec_comp = dmstodec(numpy.float(dec_comp[0]),numpy.float(dec_comp[1]),numpy.float(str(dec_comp[2]+"."+dec_comp[3])), sign=sign)
            else:
                dec_comp = dmstodec(numpy.float(dec_comp[0]),numpy.float(dec_comp[1]),numpy.float(dec_comp[2]), sign=sign)

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
                dec_comp = dmstodec(numpy.float(dec_comp[0]),numpy.float(dec_comp[1]),numpy.float(str(dec_comp[2]+"."+dec_comp[3])), sign=sign)
            else:
                dec_comp = dmstodec(numpy.float(dec_comp[0]),numpy.float(dec_comp[1]),numpy.float(dec_comp[2]), sign=sign)

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



def cal_return_slist(imagename,skymodel, direction, imsize):

    pixelsize=1.5 # arcsec
    factor = 0.8 # only add back in the center 80%
    cut = pixelsize*(imsize/2.)*factor/3600.

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

    if '-' in dec.split('d')[0]:
      ref_dec= dec1 + dec2 +dec3
    else:
      ref_dec= dec1 - dec2 -dec3

    fluxweight = False

    data = load_bbs_skymodel(skymodel)

    if len(numpy.shape(data)) == 1:  # in this case not issue and we do not use Pythonlibs
        patches,ra_patches,dec_patches, flux_patches =  compute_patch_center(data,fluxweight)
        #print 'option 1'
    if len(numpy.shape(data)) == 2:
        patches,ra_patches,dec_patches, flux_patches = compute_patch_center_libsproblem(data,fluxweight)
        #print 'option 2'

    ralist  = pi*(ra_patches)/180.
    declist = pi*(dec_patches)/180.

    ref_ra_rad=pi*ref_ra/180.0
    ref_dec_rad=pi*ref_dec/180.0
    maskimage=numpy.zeros((imsize,imsize))
    minpix=int(imsize*(1.0-factor)/2.0)
    maxpix=int(imsize*(factor+(1.0-factor)/2.0))
    maskimage[minpix:maxpix,minpix:maxpix]=1.0
    mask_img=make_image(ref_ra_rad,ref_dec_rad,maskimage,pixelsize*pi/(180.0*3600.0))

    plist = []

    #print ref_ra, ref_dec

         # load image to check if source within boundaries
    facet_img    = pyrap.images.image(imagename)
    pixels = facet_img.getdata()
    sh    = numpy.shape(pixels)[2:4]


     # CHECK TWO THINGS
     #  - sources fall within the calibration image
     #  - sources fall within the mask from the tessellation
     # both of these are done by considering the co-ordinates
    for patch_id,patch in enumerate(patches):
        # first the calibration image
        coor = [declist[patch_id],ralist[patch_id]]
        pix=mask_img.topixel(coor)
        if pix[0]<0 or pix[0]>=imsize or pix[1]<0 or pix[1]>=imsize:
            continue
        if maskimage[pix[0],pix[1]]<0.5:
            continue

        # now the facet image
        coor = [0,1,declist[patch_id],ralist[patch_id]]
        pix  = facet_img.topixel(coor)[2:4]

        if (pix[0] >= 0) and (pix[0] <= (sh[0]-1)) and \
           (pix[1] >= 0) and (pix[1] <= (sh[1]-1)):
            if pixels[0,0,pix[0],pix[1]]>0.5:  # only include if within the clean mask (==1)
                plist.append(patches[patch_id])
    
 # make the string type source list
    sourcess = ''
    if len(plist) == 1:
        sourcess = str(plist[0])
    else:
        for patch in plist:
            sourcess = sourcess+patch+','
        sourcess = sourcess[:-1]
    return sourcess,plist

def return_slist(imagename,skymodel,ref_source):
    """return a list of the sky model components in the mask defined by
    imagename, excluding all of those listed in ref_source, which have
    already been added.
    """

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
            if pixels[0,0,pix[0],pix[1]] > 0.5:  # only include if withtin the clean mask (==1)
                plist.append(patches[patch_id])

    sourcess= ''
    if len(plist) == 1:
        sourcess = str(plist[0])
    else:
        for patch in plist:
            sourcess = sourcess+patch+','
        sourcess = sourcess[:-1]

    return sourcess,plist

if __name__=='__main__':

    if len(sys.argv)==5:
        # does old cal_return_slist behaviour

        imagename = sys.argv[1]
        skymodel  = sys.argv[2]
        direction = sys.argv[3]
        imsize    = int(sys.argv[4])

        sourcess,plist=cal_return_slist(imagename,skymodel,direction,imsize)

        print sourcess
    elif len(sys.argv)==4:
        # does old return_slist behaviour

        imagename = sys.argv[1]
        skymodel  = sys.argv[2]
        ref_sourcet = sys.argv[3]
        ref_source = ref_sourcet.split(',')
        sourcess,plist=return_slist(imagename,skymodel,ref_source)
        
        print sourcess
    else:
        print 'Must have 3 or 4 arguments'
        sys.exit(-1)
