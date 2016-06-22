import numpy as np
import subprocess
import datetime
import os,sys
import re
import shutil
from part_density import calc_density
from input_file import *

time_format = "%y %m %d %H %M"

# compute the total simulation time
runtime = datetime.datetime.strptime(endtime,time_format) - datetime.datetime.strptime(starttime,time_format)

d = datetime.datetime(2000,1,1) + runtime
runtime_hh = '{0:02}'.format(int(str(d.strftime("%H"))))

# compute the number of plumemom runs to do
n_runs = np.int(np.floor( runtime.total_seconds() / deltat_plumemom ) )

d = datetime.datetime(2000,1,1) + datetime.timedelta(seconds=deltat_plumemom)
duration_hhmm = str(d.strftime("%H%M"))

duration_hhhh = '{0:04}'.format(int(str(d.strftime("%H"))))

diam = 2**(-np.asarray(diam_phi))

# density in g/cc (calc density compute it in kg/m^3)
density = calc_density(diam)/1000

shapefactor = np.ones(npart)*shapefactor

with open('EMITIMES','w') as emitimes:    
	emitimes.write('YYYY MM DD HH    DURATION(hhhh) #RECORDS \nYYYY MM DD HH MM DURATION(hhmm) LAT LON HGT(m) RATE(/h) AREA(m2) HEAT(w) \n')
emitimes.close()

# search for the maximum number of lines in the .hy files
max_lines = 0

for i in range(n_runs):

    plume_hy = runname + '_{0:03}'.format(i+1)+'.hy'

    with open(plume_hy) as f:
        max_lines = max(max_lines,sum(1 for _ in f)-1)

# loop over the .hy files to write the blocks in EMITIMES
for i in range(n_runs):

    # name of the .hy file
    plume_hy = runname + '_{0:03}'.format(i+1)+'.hy'

    # time of the block
    timei =  datetime.datetime.strptime(starttime,time_format)+datetime.timedelta(seconds=i*deltat_plumemom)
	
    timei_str = timei.strftime("%Y %m %d %H")
    timei_str_mm = timei.strftime("%Y %m %d %H %M")

    print 'Time',timei_str,'File',plume_hy

    # read the whole plumemom .hy file
    with open(plume_hy, 'r') as fin:
	    data = fin.read().splitlines(True)
    fin.close()

    # delete the header line and save to temp.hy
    with open('temp.hy', 'w') as fout:
	    fout.writelines(data[1:])
    fout.close()

    # load the data from temp.hy
    with open('temp.hy', 'r') as fin:
        data = np.loadtxt(fin)
    fin.close()
    subprocess.call('rm temp.hy', shell=True)
    # put the data in a numpy array
    data=np.asarray(data)

    # data1: array containing data from .hy file, without x,z,h
    data1=np.delete(data, [0,1,2], 1)

    # array containing lat,lon and height for time i
    b=[]

    for i0 in range(len(data)):
        x=data[i0,0] #[m]
        y=data[i0,1] #[m] 
        height=data[i0,2] #[m] 	

        # convert from m to lat lon	  
        lon_col = lon + ((x*10**-3)/float(100))
        lat_col = lat - ((y*10**-3)/float(100))

        b.append([lat_col, lon_col, height])

    b = np.asarray(b)
    b = b.reshape((-1,3))	

    # add lines in order to have all the blocks with the same lenght
    for i in range(max_lines-len(b)):

        b = np.vstack((b,b[len(b)-1,:]))
        data1 = np.vstack((data1,np.zeros(npart)))

    # b1 is an array containing lat,lon and height for time i repeated npart times
    b1=[]

    for i0 in range(len(b)):    
        for i1 in range(npart):
	        b1.append([b[i0,0],b[i0,1],b[i0,2]])

    b1=np.asarray(b1)
    b1=b1.reshape((-1,3))	

    # data3 is the array to be written in EMITTIMES for every time interval
    data3 = np.zeros((max_lines*npart,4))

    for i0 in range(max_lines):

        for i1 in range(npart):

            data3[i1+i0*npart,0:3] = b1[i1+i0*npart,0:3]

            data3[i1+i0*npart,3] = data1[i0,i1]
	
	# mass released in one hour [kg]
    emission_rate = data3[:,3]*3600/1000

    with open('EMITIMES','a') as emitimes:	

	    emitimes.write(timei_str+' '+duration_hhhh+' '+str(len(data3))+'\n')	

	    for h in range(len(data3)):
	        emitimes.write(timei_str_mm+' '+duration_hhmm+' '+
                           str(data3[h,0]) + ' ' + str(data3[h,1]) + ' ' +
                           str(data3[h,2]) + ' ' + str(emission_rate[h]) +
                           ' 0.0 0.0\n')

emitimes.close()

# write CONTROL file
file_control=open('CONTROL','w')

file_control.writelines(starttime+'\n')
file_control.writelines('%d\n'%max_lines)
for i in range(max_lines):
    file_control.writelines("%f %f %f\n"%(lat,lon,vent_height))
file_control.writelines(str(runtime_hh)+'\n')
file_control.writelines('0\n')
file_control.writelines(str(model_top)+'\n')
file_control.writelines('1\n')
file_control.writelines('./\n')
file_control.writelines(meteo_file+'\n')
file_control.writelines('%d\n'%npart)
for i in range(npart):
    file_control.writelines('CL%02d\n'%i)
    file_control.writelines('0.0\n')
    file_control.writelines('0\n')
    file_control.writelines('00 00 00 00 00\n')
file_control.writelines('1\n')
file_control.writelines('0.0 0.0\n')
file_control.writelines(str(spacing_lat)+' '+str(spacing_lon)+'\n')
file_control.writelines(str(span_lat)+' '+str(span_lon)+'\n')
file_control.writelines('./\n')
file_control.writelines('cdump_'+runname+'\n')
file_control.writelines('1\n')
file_control.writelines('15000\n')
file_control.writelines(starttime+'\n')
file_control.writelines('00 00 00 00 00\n')
file_control.writelines('1 1 0\n')
file_control.writelines('%d\n'%npart)
for i in range(npart):
    file_control.writelines('%f %f %f \n'%(diam[i],density[i],shapefactor[i]))#50.0 6.0 1.0
    file_control.writelines('0 0.0 0.0 0.0 0.0 \n')#0 0.0 0.0 0.0 0.0
    file_control.writelines('0.0 0.0 0.0 \n')#0.0 1.0E+06 1.0E-06
    file_control.writelines('0\n')#0
    file_control.writelines('0.0\n')#0.0
file_control.close()







