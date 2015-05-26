#!/usr/bin/env python

# Convert ds9 region file to facet format

import numpy as np
import sys
from subprocess import check_output

#-----------------------------------------------------------

def prime_factors(n):
    """Find the prime factors of a number
    """
    factors=[]
    d=2
    while(d*d<n):
        while(n>1):            
            while n%d==0:
                factors.append(d)
                n=n/d
            d+=1
    return factors[-1]

def close_higher_bsmooth_number(n,B):
    """Finds the closest higher B-smooth number (good numbers for FFTs) http://en.wikipedia.org/wi$
        """
    fivesmooths = np.array([1])
    i = np.max([n-300,5])
    while np.max(fivesmooths) < n:
        if prime_factors(i) <=B and (i/2.0 == int(i)/2): #Image must have an even number of pixels
            fivesmooths = np.append(fivesmooths,i)
        i+=1
    fivesmooths=np.sort(fivesmooths)
    closest_smooth_index = np.where(abs(fivesmooths-n)==np.min(abs(fivesmooths-n)))
    if len(fivesmooths[closest_smooth_index[0]]) == 2:
	closest_smooth = fivesmooths[closest_smooth_index[0][1]]
        return closest_smooth
    if fivesmooths[closest_smooth_index[0]] < n:
        closest_smooth = fivesmooths[closest_smooth_index[0]+1]
    else:
        closest_smooth = fivesmooths[closest_smooth_index[0]]
    return closest_smooth[0]


if __name__=='__main__':

    infile=sys.argv[1]
    outfile=open(sys.argv[2],'w')
    cmd='grep box '+infile+' | sed -e '
    cmd+='''\'s/box(//;s/:/h/;s/:/m/;s/:/d/;s/:/m/;s/"//g\' | awk -F, 'BEGIN{i=1}{print "s"i,$1","$2,"False False",0.5*($3+$4)/1.5,"1 60 3600 LD empty empty empty False";i++}\''''
    print cmd
    output=check_output(cmd, shell=True)
    lines=output.split('\n')
    print lines
    for line in lines:
        if len(line)>0:
            origline = line
            line = line[:-1]
            line = line.split(' ')
            while '' in line:
                line.remove('')
            origsize = float(line[4])	
            newsize = close_higher_bsmooth_number(int(origsize),5)
            if newsize < 512:
                newsize = 512
            outfile.write(origline.replace(str(origsize),str(newsize))+'\n')


