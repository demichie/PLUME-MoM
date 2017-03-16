import numpy as np
import subprocess
import datetime
import os,sys,time
import shutil

def calc_atm(profile,fields_list):
    """create atmosperic profile for plumemom """
    nlev=len(profile[:,0])

    print 'nlev',nlev
    grav = 9.81
    Rd = 287.04
    Rw = 461.5

    P0 = 1011*100.0
    T0 = 25.8+273.15
    RH0 = 25.8
    z0 = 24.0

    Es0 = np.exp(21.4-(5351/T0))
    Q0 = RH0 / P0 * ( 0.622 * Es0 )
    Tv0 = T0 * ( 1 + 0.61*Q0)

    P = profile[:,0]*100.0

    T_idx = fields_list.index("TEMP") + 1

    T = profile[:,T_idx] + 273.15

    Es = np.exp(21.4-(5351/T))

    # dry air density
    rho_da = P / ( Rd * T )

    if 'SPHU' in fields_list:

        SPHU_idx = fields_list.index("SPHU") + 1

        Q = profile[:,SPHU_idx]/1000.0

        RH = (Q * P) / ( 0.622 * Es )

    else:

        RH_idx = fields_list.index("RELH") + 1

        RH = profile[:,RH_idx]

        Q = RH / P * ( 0.622 * Es )

    # Mixture density, Eq. (3) http://www.engineeringtoolbox.com/density-air-d_680.html
    rho = rho_da * ( 1 + Q ) / ( 1 + Q * Rw / Rd )

    if 'HGTS' in fields_list:

        z_idx = fields_list.index("HGTS") + 1

        z = profile[:,z_idx]

    else:

        # Tv = T .* ( 1 + 0.61*min(1,Q))
        Tv = T * ( 1 + 0.61*Q)

        Tv01 = 0.5*(Tv0+Tv[0])
        deltaz0 = np.log(P0/P[0]) *Rd *Tv01 / grav
        z=np.zeros(nlev)
        z[0] = z0 + deltaz0

        Tv12 = 0.5*(Tv[0:-1]+Tv[1:])
        deltaz = np.log(P[0:-1]/P[1:]) *Rd *Tv12 / grav

        z[1:] = z[0] + np.cumsum(deltaz)

	# W->E component of horizontal velocity (m/s)
    U_idx = fields_list.index("UWND") + 1
    U = profile[:,U_idx]

	# S->N component of horizontal velocity (m/s)
    V_idx = fields_list.index("VWND") + 1
    V = profile[:,V_idx]


#	print(len(z),len(rho),len(P),len(T),len(RH),len(profile[:,10]))
    return z/1000.0, rho, P/100, T, RH, U, V

def write_atm(time_input):
    """create atmosperic profile for plumemom """

    from input_file import hysplit_dir, meteo_file, lat, lon

    run_path = './'

    profile = os.path.join(hysplit_dir,'exec','profile')

    subprocess.call(profile+" -d"+run_path+" -f"+meteo_file+" -y"
                    +str(lat)+" -x"+str(lon)+" -o0 -p00", shell=True) # wind profile at vent position at the beginnig of the wind file

    with open("profile_00.txt","r") as file1:
        line=file1.readlines()
        time_start_line=line[1]
        time_start=time_start_line[19:35]

        # add 0 to year if it is a single digit (ex: 6 for 2006 becomes 06)
        time_split = time_start.split()
        year = time_split[0]
        if len(year)==1:
            time_start = '0'+time_start[1:14]

    file1.close()

    for i in range(len(line)):
        if ( '3D Field' in line[i] ):
            print 'success',i
            idx_fields_line = i+1
            fields_list = line[idx_fields_line].split()

    a= datetime.datetime.strptime(str(time_start[0:14]), "%y %m %d %H %M")
    b= datetime.datetime.strptime(str(time_input), "%y %m %d %H %M")

    delta = b-a
    time_difference_in_hours = int(delta.total_seconds() / 3600) # hours of difference between time and the start of the wind file

    subprocess.call(profile+ " -d"+run_path+" -f"+meteo_file
                    + " -y"+str(lat)+" -x"+str(lon)+" -o"
                    +str(time_difference_in_hours)+" -p01", shell=True)
    
    with open('atm_profile.txt','w') as file: # output file: wind profile at vent position at time
        with open("profile_01.txt") as fp:
            for i, line in enumerate(fp):
                if i > 13:
                    file.write(str(line[0:64])+'\n')
        fp.close()
    file.close()

             
    subprocess.call("rm " +run_path+"profile_00.txt", shell=True) 
    subprocess.call("rm " +run_path+"profile_01.txt", shell=True) 
    

    # read the file 'atm_profile.txt'
    arrays = [np.array(map(float, line.split())) for line in open('atm_profile.txt')]

    # compute the number of values in the first level
    n0 = len(arrays[0])
    b = np.zeros((len(arrays),n0))

    # search for the number of levels with the same number of values
    nfull = 0
    for i in range(len(arrays)):
        if (len(arrays[i]) == n0):
            b[i,:] = arrays[i]
            nfull = nfull+1

    # extract only the lines with the same number of values (n0)
    profile = np.asarray(b[0:nfull,:])

    # compute the missing fields
    a = calc_atm(profile,fields_list)
    a = np.asarray(a)


    nrows = a.shape[1]

    f=open('atm.txt','w')

    f.write(str(nrows)+'\n')
    f.close()
    
    f=open('atm.txt','a')

    np.savetxt(f,a.transpose(),fmt='%.4e')

    return







