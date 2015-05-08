import os
import glob
import sys


filelist = glob.glob('3C196_24_SB???.ms')


for ms in filelist:
  print ms
  os.system('makevds ~/cep2.clusterdesc ' + ms)
