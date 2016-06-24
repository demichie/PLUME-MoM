import numpy as np
import subprocess
import datetime
import os,sys
import re
import shutil
from extract_wind import write_atm

from input_file import *

diam1=np.ones(npart)*diam1
diam2=np.ones(npart)*diam2

cp_part=np.ones(npart)*cp_part
rho1 = np.ones(npart)*rho1
rho2 = np.ones(npart)*rho2

diam2=np.ones(npart)
diam1=.1*diam2

# create a second template with the parameters constant in time
f = open('plume_model.template','r')
filedata = f.read()
f.close()

filedata = filedata.replace("{vent_radius}", str(vent_radius) )

filedata = filedata.replace("{log10_mfr}", str(log10_mfr) )

filedata = filedata.replace("{npart}", str(npart) )

filedata = filedata.replace("{shapefactor}", str(shapefactor) )

filedata = filedata.replace("{deltaz_release}", str(deltaz_release) )


filedata = filedata.replace("{vent_height}", str(vent_height) )

filedata = filedata.replace("{diam1}", ",".join(np.char.mod('%f', diam1)) )
filedata = filedata.replace("{diam2}", ",".join(np.char.mod('%f', diam2)) )

filedata = filedata.replace("{rho1}", ",".join(np.char.mod('%f', rho1)) )
filedata = filedata.replace("{rho2}", ",".join(np.char.mod('%f', rho2)) )

filedata = filedata.replace("{cp_part}", ",".join(np.char.mod('%f', cp_part)) )

filedata = filedata.replace("{solid_partial_mass_fraction}", ",".join(np.char.mod('%f', partial_mass_fractions)) )
filedata = filedata.replace("{diam_constant_phi}", ",".join(np.char.mod('%f', diam_phi)) )

f = open('plume_model.temp1','w')
f.write(filedata)

f.close()


time_format = "%y %m %d %H %M"

# compute the total simulation time
runtime = datetime.datetime.strptime(endemittime,time_format) - datetime.datetime.strptime(starttime,time_format)

# compute the number of plumemom runs to do
n_runs = np.int(np.floor( runtime.total_seconds() / deltat_plumemom ) )

for i in range(n_runs):

    runnamenew = runname + '_{0:03}'.format(i+1)
    print 'runname',runnamenew

    f = open('plume_model.temp1','r')
    filedata = f.read()
    f.close()

    # create a third template with the parameters changing with time
    f = open('plume_model.temp2','w')

    filedata = filedata.replace("{runname}", '"'+str(runnamenew)+'"' )

    f.write(filedata)
    f.close()

    timei =  datetime.datetime.strptime(starttime,time_format)+datetime.timedelta(seconds=i*deltat_plumemom)
	
    timei_str = timei.strftime("%y %m %d %H %M")

    print 'Time',timei_str

    write_atm(timei_str)


    # append the atmospheric data to the input file
    filenames = ['plume_model.temp2', 'atm.txt']
    with open('plume_model.inp', 'w') as outfile:
        for fname in filenames:
            with open(fname) as infile:
                for line in infile:
                    outfile.write(line)


    subprocess.call(plumemom_dir+"/bin/PLUMEMoM", shell=True) 

subprocess.call("rm plume_model.temp1", shell=True) 
subprocess.call("rm plume_model.temp2", shell=True) 




