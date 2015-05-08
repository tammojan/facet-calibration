import os,sys
import glob

filenames = glob.glob('./*tar')
for filename in filenames:
	os.system('tar -xf %s'%filename)
	os.system('rm -rf %s'%filename)
