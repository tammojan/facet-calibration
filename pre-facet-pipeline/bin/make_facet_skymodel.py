#! /usr/bin/env python
import argparse
from argparse import RawTextHelpFormatter
import lsmtool
from lsmtool.operations_lib import radec2xy
import matplotlib.path as mplPath
from numpy import array, zeros
import sys
import os


def main(fullskymodel, outmodel, cal_only=False, vertices=None, facet_ra=0.0, facet_dec=0.0, cal_radius_deg=0.0):
    print 'Running {0} on input data {1}'.format(__file__, str(fullskymodel))

    s = lsmtool.load(fullskymodel)

    if cal_only:
        # Get calibrator model
        dist = s.getDistance(facet_ra, facet_dec)
        s.select(dist < cal_radius_deg)
    else:
        # Get all facet sources
        x, y, midRA, midDec = s._getXY()
        xv, yv = radec2xy(vertices[0], vertices[1], midRA, midDec)
        xyvertices = array([[xp, yp] for xp, yp in zip(xv, yv)])
        bbPath = mplPath.Path(xyvertices)
        inside = zeros(len(s), dtype=bool)
        for i in range(len(s)):
            inside[i] = bbPath.contains_point((x[i], y[i]))
        s.select(inside, force=True)

    if len(s) == 0:
        print('No sources found for this facet')
    else:
        s.write(outmodel, clobber=True)


if __name__ == '__main__':
    descriptiontext = "Put a wall of text here.\n" + \
        "This will show when called with option     --help\n" + \
        "Describe what this script is for.\n" + \
        "And maybe what the parameters are for."

    parser = argparse.ArgumentParser(description=descriptiontext, formatter_class=RawTextHelpFormatter)
    # positional arguments
    parser.add_argument('fullskymodel', help='name of the full skymodel')
    parser.add_argument('outmodel', help='name for the output')
    # optional arguments with short and long names and defaults.
    # long names suitable for lofar parset files
    parser.add_argument('-c', '--cal_only', help='a bool for...', type=bool, default=False)
    parser.add_argument('-t', '--vertices', help='', default=None)
    parser.add_argument('-r', '--facet_ra', help='', type=float, default=0.0)
    parser.add_argument('-d', '--facet_dec', help='', type=float, default=0.0)
    parser.add_argument('-d', '--cal_radius_deg', help='', type=float, default=0.0)

    args = parser.parse_args()

    # do funny stuff to the parameters here if the formatting had to be
    # inacurate for historical/compatibility reasons and does not match your functions

    main(args.fullskymodel, args.outmodel, args.cal_only, args.vertices, args.facet_ra, args.facet_dec, args.cal_radius_deg)