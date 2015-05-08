#!/usr/bin/python

# msconcat.py
#
# concatenate N ms's in time
#
# usage: python msconcat.py output.ms input1.ms .. inputN.ms
# written by: Wendy Williams wwilliams@strw.leidenuniv.nl
# last modified: 25 feb 2015



import sys
from pyrap.tables import table

if __name__ == "__main__":
    print "Concatenate MeasurementSets\n"
    t = table(sys.argv[2:])
    t.sort('TIME').copy(sys.argv[1], deep = True)
