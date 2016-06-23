import glob, os
from input_file import *

filelist = glob.glob(runname+'*')
for f in filelist:
    os.remove(f)


filelist = glob.glob('atm*')
for f in filelist:
    os.remove(f)


filelist = glob.glob('cdump*')
for f in filelist:
    os.remove(f)


filelist = glob.glob('*.CFG')
for f in filelist:
    os.remove(f)


filelist = glob.glob('*.ps')
for f in filelist:
    os.remove(f)


filelist = glob.glob('*.inp')
for f in filelist:
    os.remove(f)


filelist = glob.glob('*.pyc')
for f in filelist:
    os.remove(f)

filelist = glob.glob('*~')
for f in filelist:
    os.remove(f)

                              
os.remove('EMITIMES')
os.remove('MESSAGE')
os.remove('PARDUMP')
os.remove('WARNING')
os.remove('CONTROL')
os.remove('STARTUP')
os.remove('VMSDIST')                                                         
                       
                       
     

