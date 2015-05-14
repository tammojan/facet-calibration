#!/usr/bin/env python

import sys
import numpy
import scipy.ndimage
import pyrap.images
from optparse import OptionParser

# Return the fractional part of the floating point number x.
def fraction(x):
    return numpy.modf(x)[0]

# Convert an angle from degrees to radians.
def deg2rad(angle):
    return (angle * numpy.pi) / 180.0

# Convert an angle from radians to degrees.
def rad2deg(angle):
    return (angle * 180.0) / numpy.pi

# Compute hours, min, sec from an angle in radians.
def rad2ra(angle):
    deg = numpy.fmod(rad2deg(angle), 360.0)
    if deg < 0.0:
        deg += 360.0
    # Ensure positive output (deg could equal -0.0).
    deg = abs(deg)
    assert(deg < 360.0)

    return (int(deg / 15.0), int(fraction(deg / 15.0) * 60.0), fraction(deg * 4.0) * 60.0)

# Compute degrees, arcmin, arcsec from an angle in radians, in the range [-90.0, +90.0].
def rad2dec(angle):
    deg = numpy.fmod(rad2deg(angle), 360.0)
    if deg > 180.0:
        deg -= 360.0
    elif deg < -180.0:
        deg += 360.0

    sign = (deg < 0.0)
    deg = abs(deg)
    assert(deg <= 90.0)

    return (sign, int(deg), int(fraction(deg) * 60.0), fraction(deg * 60.0) * 60.0)

# Return string representation of the input right ascension (as returned by rad2ra).
def ra2str(ra):
    return "%02d:%02d:%07.4f" % ra

# Return string representation of the input declination (as returned by rad2dec).
def dec2str(dec):
    return "%s%02d.%02d.%07.4f" % ("-" if dec[0] else "+", dec[1], dec[2], dec[3])

# Recursively flatten a list.
def flatten_list(l):
    def flatten_rec(l, out):
        if type(l) == list:
            for el in l:
                flatten_rec(el, out)
        else:
            out.append(l)

    out = []
    flatten_rec(l, out)
    return out

# Compute an index that maps each element in requested to the index of the same
# element in available. If such an element cannot be found and mandatory is
# set to True an exception is raised. If mandatory is set to False, None is
# stored in the index instead.
def build_index(available, requested, mandatory=True):
    index = [None for i in range(len(requested))]
    for i in range(len(requested)):
        try:
            index[i] = available.index(requested[i])
        except ValueError:
            if mandatory:
                raise

    return index

# Return a sequence that contains the index of each element in the input
# sequence if the input sequence would be sorted (ascendingly).
#
# For example: [4,2,9,11] => [1,0,2,3]
#
def renumber(sequence):
    index = sorted(range(len(sequence)), key=lambda x: sequence[x])
    return sorted(range(len(index)), key=lambda x: index[x])

# Load pixel data from the input CASA image and create a slice that contains
# only the requested axes in order. Zero is used as the coordinate for each
# additional axis that the image contains when creating the slice. Both the
# slice and the index of each requested axis in the original image is returned.
def load_data_and_reformat(im, axes):
    # Get image coordinate system and find the index of the requested axes.
    csys = im.coordinates()
    try:
        index = build_index(flatten_list(csys.get_axes()), axes)
    except:
        raise Exception("error: incompatible image format: one or more required axes undefined.")

    # Load pixel data.
    data = im.getdata()

    # Create slice.
    slicer = [0 for i in range(len(data.shape))]
    for i in range(len(axes)):
        slicer[index[i]] = slice(data.shape[index[i]])

    return (numpy.transpose(data[slicer], renumber(index)), index)

def main(options, args):
    assert(len(args) >= 1)

    # Load CLEAN component map and reformat.


    if options.nterms > 1:
        component_map_im      = pyrap.images.image(args[0]+'.tt0')
        component_map_im_tt1  = pyrap.images.image(args[0]+'.tt1')
        if options.nterms > 2:
            component_map_im_tt2  = pyrap.images.image(args[0]+'.tt2')
    else:
        component_map_im = pyrap.images.image(args[0])

    ref_freq = component_map_im.info()['coordinates']['spectral2']['restfreqs'][0]


    if options.nterms > 1:
        component_map, axis_index = load_data_and_reformat(component_map_im, ["Stokes", "Declination", "Right Ascension"])
        component_map_tt1, axis_index_tt1 = load_data_and_reformat(component_map_im_tt1, ["Stokes", "Declination", "Right Ascension"])
        if options.nterms > 2:
            component_map_tt2, axis_index_tt2 = load_data_and_reformat(component_map_im_tt2, ["Stokes", "Declination", "Right Ascension"])
    else:
        component_map, axis_index = load_data_and_reformat(component_map_im, ["Stokes", "Declination", "Right Ascension"])

    # Ugly kludge to find out mapping from Stokes coordinates to indices. No
    # better alternative seems to be available using pyrap at this time.
    stokes = component_map_im.coordinates().get_coordinate("stokes")._coord["stokes"]
    stokes_index = build_index(stokes, ["I", "Q", "U", "V"], False)

    # Make sure Stokes I is available.
    if stokes_index[0] is None:
        print "error: incompatible CLEAN component image format: Stokes I unavailable."
        sys.exit(1)

    # Find locations of all CLEAN components (pixels with non-zero flux).
    component_flux = component_map[stokes_index[0], :, :]
    components = numpy.nonzero(component_flux)

    if options.nterms > 1:
        component_flux_tt1 = component_map_tt1[stokes_index[0], :, :]
        components_tt1 = numpy.nonzero(component_flux_tt1)
        if options.nterms > 2:
            component_flux_tt2 = component_map_tt2[stokes_index[0], :, :]
            components_tt2 = numpy.nonzero(component_flux_tt2)



    # Compute statistics.
    total_component_count = len(components[0])
    print "info: total number of CLEAN components: %d" % total_component_count
    if total_component_count == 0:
        print "error: no CLEAN components found in CASA image:", args[0]
        sys.exit(1)

    total_component_flux = numpy.sum(component_flux[components])
    print "info: total flux in CLEAN components: %.2f Jy" % total_component_flux

    have_mask = not options.mask is None
    if have_mask:
        # Load CASA mask image and reformat.
        mask_im = pyrap.images.image(options.mask)
        mask, mask_index = load_data_and_reformat(mask_im, ["Declination", "Right Ascension"])

        # Sanity check mask image dimensions.
        if mask.shape != component_flux.shape:
            print "error: CASA mask image should have the same dimensions as the CLEAN component image."
            sys.exit(1)

        # Find islands in the mask image.
        (island_label, island_count) = scipy.ndimage.measurements.label(mask)
        if island_count == 0:
            print "error: no islands found in CASA mask image:", options.mask
            sys.exit(1)

        # Assign CLEAN components to the corresponding (co-located) islands and
        # keep track of the total flux in each island. Also compute the total
        # flux present in all the islands.
        total_flux = 0.0
        islands = [[0.0, []] for i in range(island_count)]
        for component in zip(*components):
            # Get the label of the island at the position of the CLEAN
            # component.
            label = island_label[component]

            # Label 0 is used for the background of the mask image and as
            # such does not define an island, so skip it.
            if label == 0:
                continue

            # Get the absolute flux for the current CLEAN component.


            flux = component_flux[component]

            # Append CLEAN component to the right patch and update the patch
            # flux.
            islands[label - 1][0] += flux
            islands[label - 1][1].append(component)

            # Update total flux.
            total_flux += flux

        # Prune islands that do not contain any CLEAN components.
        islands = filter(lambda x: len(x[1]) > 0, islands)

        # Print total number of non-empty island.
        print "info: total number of non-empty islands: %d" % len(islands)
        print "info: total flux in non-empty islands: %.2f Jy" % total_flux

        # Determine clip level.
        clip_flux = options.clip_level / 100.0 * total_flux
        print "info: clipping at: %.2f Jy" % clip_flux

        # Sort islands by absolute flux (descending).
        islands = sorted(islands, key = lambda x: abs(x[0]), reverse = True)

        # Add islands to the patch list until the absolute flux clip level is
        # reached.
        patches = []
        component_count = 0
        sum_flux = 0.0
        for island in islands:
            if options.clip_level < 100.0 and sum_flux + island[0] > clip_flux:
                break

            sum_flux += island[0]
            component_count += len(island[1])
            patches.append(island[1])

        print "info: number of non-empty islands selected: %d" % len(patches),
        if len(islands) > 0:
            print "(%.2f%%)" % ((100.0 * len(patches)) / len(islands))
        else:
            print "(100.00%)"
        if len(patches) == 0:
            print "warning: no non-empty islands selected; you may need to raise the clip level (option -c)."

        print "info: number of CLEAN components in selected non-empty islands: %d (%.2f%%)" % (component_count, (100.0 * component_count) / total_component_count)
        print "info: flux in selected non-empty islands: %.2f Jy (%.2f%%)" % (sum_flux, (100.0 * sum_flux) / total_flux)
    else:
        # Without a mask we are selecting CLEAN components, so the total
        # absolute flux is equal to the sum of the absolute flux of all CLEAN
        # components.
        total_flux = total_component_flux

        # Determine clip level.
        clip_flux = options.clip_level / 100.0 * total_flux
        print "info: clipping at: %.2f Jy" % clip_flux

        # Sort CLEAN components on absolute flux (descending).
        components = sorted(zip(*components),
            key = lambda x: abs(component_flux[x]), reverse = True)

        # Create a single patch and add CLEAN components to it until the
        # absolute flux clip level is reached.
        patch = []
        sum_flux = 0.0
        for component in components:
            flux = component_flux[component]
            if options.clip_level < 100.0 and sum_flux + flux > clip_flux:
                break

            sum_flux += flux
            patch.append(component)

        patches = []
        if len(patch) > 0:
            patches.append(patch)

        print "info: number of CLEAN components selected: %d (%.2f%%)" % (len(patch), (100.0 * len(patch)) / total_component_count)
        if len(patch) == 0:
            print "warning: no CLEAN components selected; you may need to raise the clip level (option -c)."
        print "info: flux in selected CLEAN components: %.2f Jy (%.2f%%)" % (sum_flux, (100.0 * sum_flux) / total_flux)

    # Open output file.
    out = sys.stdout
    if len(args) == 1:
        out = file("%s.catalog" % args[0], 'w')
    elif args[1] != "-":
        out = file(args[1], 'w')

    # Write the catalog header.
    if options.use_patches:
        #print >>out, "# (Name, Type, Patch, Ra, Dec, I, Q, U, V) = format"
        print >>out, "format = Name, Type, Patch, Ra, Dec, I, Q, U, V, MajorAxis, MinorAxis, Orientation, ReferenceFrequency='%f', SpectralIndex='[]'" %ref_freq
    else:
        print >>out, "# (Name, Type, Ra, Dec, I, Q, U, V) = format"

    print >>out
    print >>out, "# CLEAN component list converted from:", args[0]
    if not options.mask is None:
        print >>out, "# Mask:", options.mask
    print >>out, "# Total flux in CLEAN components: %.2f Jy" % total_flux
    print >>out, "# Percentage of total flux kept: %.2f%%" % options.clip_level
    print >>out

    # Output all the patches in BBS catalog file format.
    patch_count = 0
    component_count = 0
    stokes_desc = ["" for i in range(len(stokes_index))]
    pixel_coord = [0 for i in range(len(component_map_im.shape()))]
    for patch in patches:
        # Output patch definition.
        if options.use_patches:
            print >>out, ", , patch-%d, 00:00:00, +90.00.00" % patch_count

        # When using patches, reset the component count at the start of each
        # patch. Thus, components within a patch are counted from zero.
        # Components are labelled both with a patch count and a component count
        # to avoid name clashes. If not using patches, just continue counting.
        if options.use_patches:
            component_count = 0

        # Output all the CLEAN components in the patch.
        for component in patch:
            if options.use_patches:
                print >>out, "patch-%d-%d, POINT, patch-%d," % (patch_count, component_count, patch_count),
            else:
                print >>out, "component-%d, POINT," % (component_count),

            pixel_coord[axis_index[1]] = component[0]
            pixel_coord[axis_index[2]] = component[1]
            world_coord = component_map_im.toworld(pixel_coord)
            print >>out, "%s," % ra2str(rad2ra(world_coord[axis_index[2]])),
            print >>out, "%s," % dec2str(rad2dec(world_coord[axis_index[1]])),

            for i in range(len(stokes_index)):
                if stokes_index[i] is None:
                    stokes_desc[i] = "0.0"
                else:
                    stokes_desc[i] = "%f" % component_map[(stokes_index[i], component[0], component[1])]
            print >>out, ", ".join(stokes_desc),
            print >>out, ",0.0,0.0,0.0,",
            print >>out, "%f," % ref_freq,
            ## ALPHA = TT1/TT0 from http://casa.nrao.edu/docs/userman/UserMansu247.html
            #print numpy.shape(component_map[(stokes_index[i], component[0], component[1])])
            alpha = 0.0 # in case we do not use spectral index (nterms=1)
            if options.nterms > 1:
                alpha = -0.7 # default spectral index
                alpha = numpy.copy( (component_map_tt1[(0,component[0], component[1])]) / (component_map[(0, component[0], component[1])]) )
            # to prevent crazy values
                if (numpy.abs(alpha) > 6.0): # 10, noticed unstable behaviour if 1000
                    alpha = -0.7
                if options.nterms > 2:
                    beta = numpy.copy( ((component_map_tt2[(0,component[0], component[1])])/(component_map[(0, component[0], component[1])])) -\
                                         (0.5*alpha*(alpha-1.0)))
                    if (numpy.abs(beta) > 6.0): # 10, noticed unstable behaviour if 1000
                        beta = 0.0
                        if alpha == -0.7: # because then the beta is flawed anyway by wrong alpha
                            beta = 0.0

            if (options.nterms) == 1 or (options.nterms == 2):
                print >>out, "[%f]" % alpha
            if (options.nterms) == 3:
                spix_str = str(alpha)+','+str(beta)
                print >>out, "[%s]" % spix_str

            component_count += 1

        print >>out
        patch_count += 1

parser = OptionParser(usage="%prog [options] <CLEAN component image> [output catalog file]")
parser.add_option("-n", "--no-patches", action="store_false", dest="use_patches", default=True, help="Do not group components into patches. All CLEAN components will be written as separate components.")
parser.add_option("-m", "--mask", dest="mask", help="CASA mask image; If provided, the islands in the mask image define separate patches; CLEAN components that do not fall on an island are discarded, as are islands that do not contain any CLEAN components.")
parser.add_option("-c", "--clip-level", type="float", dest="clip_level", default=100.0, help="Percentage [0.0, 100.0] of the total flux that should be kept (default: %default). If there are CLEAN components with negative flux in the CLEAN component image, you may not recover all CLEAN components even if you specify 100% here. This is due to the way the cumulative flux builds up in the presence of negative components.")
parser.add_option("-t", "--nterms", type="int", dest="nterms", default=1, help="Use in casa nterms>1, NOTE: currently only works with 2")
(options, args) = parser.parse_args()

if options.clip_level < 0.0 or options.clip_level > 100.0:
    parser.error("option -c: invalid percentage: %.2f" % options.clip_level)

if len(args) < 1 or len(args) > 2:
    parser.error("incorrect number of arguments")

main(options, args)
