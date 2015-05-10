import glob
import os
import sys


# example script how to restore a bad subtract (when verify subtract indicates problems)

mslist = sorted(glob.glob('A2256_SB*0-*9.2ch10s.ms'))

for ms in mslist:
  print "taql 'update " +ms  +" set SUBTRACTED_DATA_ALL=CORRECTED_DATA'"
  os.system("taql 'update " +ms  +" set SUBTRACTED_DATA_ALL=CORRECTED_DATA' &")


