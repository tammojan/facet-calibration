SCRIPTPATH      = '/home/'+username+'/scripts/a2256_hba/'

peelsourceinfo  =  '/home/'+username+'/scripts/a2256_hba/peel_source_info_a2256.txt'

ms              =  '/data2/rvweeren/a2256_hba/SB050-359/A2256_SB050-059.2ch10s.ms'

resolution      = 1.5        # arcsec per pixel
lowresolution   = 15      # arcsec per pixel
max_fieldsize   = 5600   # don't allow any final masks to be bigger than this (at fullres)

maxsize_fieldimage = 12 ### degrees

# set max sizes of inner and outer facets (within and out fwhm)
maxoutliersize = 4800
maxcentralsize = 6400

FREQ = 150.     # MHz  - used for plotting approx FWHM sizes

make_mosaic = True      # end with a final (SLOW) step of mosaicing the templates - useful check
