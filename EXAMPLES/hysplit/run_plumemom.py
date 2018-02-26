import numpy as np
import subprocess
import datetime
import os,sys
import re
import shutil
from extract_wind import write_atm

from input_file import *

def round_minutes(dt, direction, resolution):

    if ( dt.minute%resolution == 0 ):

        rounded_time = dt

    else: 

        new_minute = (dt.minute // resolution + (1 if direction == 'up' else 0)) * resolution

        rounded_time = dt + datetime.timedelta(minutes=new_minute - dt.minute)

    return rounded_time


# create a backup of the input file
src = 'input_file.py'
dst = runname+'.bak'

shutil.copyfile(src, dst)

diam1=np.ones(npart)*diam1
diam2=np.ones(npart)*diam2

cp_part=np.ones(npart)*cp_part
rho1 = np.ones(npart)*rho1
rho2 = np.ones(npart)*rho2

# create a second template with the parameters constant in time
f = open('plume_model.template','r')
filedata = f.read()
f.close()

filedata = filedata.replace("{vent_velocity}", str(vent_velocity) )

filedata = filedata.replace("{npart}", str(npart) )

filedata = filedata.replace("{shapefactor}", str(shapefactor) )

filedata = filedata.replace("{deltaz_release}", str(deltaz_release) )

filedata = filedata.replace("{ncloud}", str(ncloud) )

filedata = filedata.replace("{vent_height}", str(vent_height) )

filedata = filedata.replace("{diam1}", ",".join(np.char.mod('%4f', diam1)) )
filedata = filedata.replace("{diam2}", ",".join(np.char.mod('%4f', diam2)) )

filedata = filedata.replace("{rho1}", ",".join(np.char.mod('%4f', rho1)) )
filedata = filedata.replace("{rho2}", ",".join(np.char.mod('%4f', rho2)) )

filedata = filedata.replace("{cp_part}", ",".join(np.char.mod('%4f', cp_part)) )

filedata = filedata.replace("{solid_partial_mass_fraction}", ",".join(np.char.mod('%f', partial_mass_fractions)) )
filedata = filedata.replace("{diam_constant_phi}", ",".join(np.char.mod('%f', diam_phi)) )


filedata = filedata.replace("{ngas}", str(ngas) )

if ngas>0:

    filedata = filedata.replace("{rvolcgas}", ",".join(np.char.mod('%f', rvolcgas)) )
    filedata = filedata.replace("{cpvolcgas}", ",".join(np.char.mod('%f', cpvolcgas)) )
    filedata = filedata.replace("{volcgas_mol_wt}", ",".join(np.char.mod('%f', volcgas_mol_wt)) )
    filedata = filedata.replace("{volcgas_mass_fraction}", ",".join(np.char.mod('%f', volcgas_mass_fraction)) )

else:

    filedata = filedata.replace("{rvolcgas}", "" )
    filedata = filedata.replace("{cpvolcgas}", "" )
    filedata = filedata.replace("{volcgas_mol_wt}", "" )
    filedata = filedata.replace("{volcgas_mass_fraction}", "" )


filedata = filedata.replace("{water_mass_fraction0}", str(water_mass_fraction0) )

f = open('plume_model.temp1','w')
f.write(filedata)

f.close()


time_format = "%y %m %d %H %M"

starttime_hhmm = datetime.datetime.strptime(starttime,time_format)
starttime_round = round_minutes(starttime_hhmm, 'down', 60) # arrotonda per difetto starttime

endemittime_hhmm = datetime.datetime.strptime(endemittime,time_format)
endemittime_round = round_minutes(endemittime_hhmm, 'up', 60) # arrotonda per eccesso endemittime

print 'starttime',starttime_hhmm,starttime_round
print 'endemittime',endemittime_hhmm,endemittime_round


runtime=endemittime_round-starttime_round # numero ore arrotondate tra inizio e fine emissione 
n_runs = np.int(np.floor( runtime.total_seconds() / deltat_plumemom ) ) # numero run di PlumeMoM


if isinstance(log10_mfr, (np.ndarray) ):

    if ( len(log10_mfr) != n_runs ):

        print 'WARNING: check numbers of values of log10_mfr',len(log10_mfr),n_runs
        sys.exit()

else:

    log10_mfr = np.ones(n_runs)*log10_mfr


for i in range(n_runs):

    runnamenew = runname + '_{0:03}'.format(i+1)
    print 'runname',runnamenew

    f = open('plume_model.temp1','r')
    filedata = f.read()
    f.close()

    # create a third template with the parameters changing with time
    f = open('plume_model.temp2','w')

    filedata = filedata.replace("{runname}", '"'+str(runnamenew)+'"' )


    filedata = filedata.replace("{log10_mfr}", str(log10_mfr[i]) )

    
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




