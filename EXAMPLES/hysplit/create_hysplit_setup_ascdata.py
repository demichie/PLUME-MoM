import numpy as np
import subprocess
import datetime
import os,sys
import re
import shutil
from part_density import calc_density
from input_file import *

with open('SETUP.CFG','w') as setup:
    setup.write('&SETUP  \n')
    setup.write('kmsl='+str(kmsl)+'\n')
    setup.write('ninit='+str(ninit)+'\n')
    setup.write('ndump='+str(ndump)+'\n')
    setup.write('ncycl='+str(ncycl)+'\n')
    setup.write("efile = 'EMITIMES', \n")
    setup.write('numpar='+str(numpar)+'\n')
    setup.write('maxpar='+str(maxpar)+'\n')
    setup.write("pinpf = 'PARINIT', \n")
    setup.write("poutf = 'PARDUMP', \n")
    setup.write('/ \n')
setup.close()

with open('ASCDATA.CFG','w') as ascdata:
    ascdata.write('-90.0  -180.0 \n')	
    ascdata.write('1.0     1.0 \n')
    ascdata.write('180     360 \n')
    ascdata.write('2 \n')
    ascdata.write('0.2 \n')
    ascdata.write("'"+hysplit_dir+"/trunk/bdyfiles/'\n")
ascdata.close()


