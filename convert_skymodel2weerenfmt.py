#!/usr/bin/env python
"""
Script to combine the low and high resolution skymodels
"""
import os
import sys
import argparse

SCRIPTPATH = os.path.dirname(os.path.abspath(__file__))


#bands=[17]

# this script assumes your highres and lowres images live in a directory BAND{b}/suball_nobeam/

NAME = "L138658_T2"
RES = "3ch10s"


#for b in bands:
def transform(b, name=NAME, res=RES):
    """
    Run in the directory with the data
    """
    ms = "{name}_SB{b:02d}0-{b:02d}9.{res}.ms".format(b=b,name=name,res=res)

    rootname = ms.split('.')[0]
    highresmodel = "{name}highresmask.skymodel".format(name=rootname)
    lowresmodel = "{name}lowresmask.skymodel".format(name=rootname)
    model = "{name}.skymodel".format(name=rootname)

    cmd=SCRIPTPATH+"/use/casapy2bbs_one_patch_per_cc.py {name}highresmask.model {name}highresmask.skymodel".format(name=rootname,b=b)
    print cmd
    os.system(cmd)
    cmd=SCRIPTPATH+"/use/casapy2bbs_one_patch_per_cc.py {name}lowresmask.model {name}lowresmask.skymodel".format(name=rootname,b=b)
    print cmd
    os.system(cmd)
    cmd="cp {name}highresmask.skymodel {name}.skymodel".format(name=rootname,b=b)
    os.system(cmd)
    cmd="grep -v '#' {name}lowresmask.skymodel >tmp${b}.sky ; cat tmp${b}.sky >> {name}.skymodel".format(name=rootname,b=b)
    print cmd
    os.system(cmd)
    cmd="sed -i 's/, , , /, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.5e8, [0.0]/g' {name}.skymodel".format(name=rootname,b=b)
    print cmd
    os.system(cmd)
    cmd="sed -i \"s/# (Name, Type, Patch, Ra, Dec, I, Q, U, V) = format/format = Name, Type, Patch, Ra, Dec, I, Q, U, V, MajorAxis, MinorAxis, Orientation, ReferenceFrequency=\'1.5e+08\', SpectralIndex=\'[]\'/g\" {name}.skymodel".format(name=rootname,b=b)
    print cmd
    os.system(cmd)
    cmd="sed -i 's/BAND{b:02d}\/suball_nobeam\///g' {name}.skymodel".format(name=rootname,b=b)
    print cmd
    os.system(cmd)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Transform the skymodels to the correct format')
    parser.add_argument('band', metavar='b', type=int, help='band')
    parser.add_argument('--name', default="L138658_T2", help='Name prefix')
    parser.add_argument('--res', default="3ch10s", help='resolution configuration of the ms')

    args = parser.parse_args()

    print args.band
    print args.name
    print args.res

    transform(int(args.band), name=args.name, res=args.res)
