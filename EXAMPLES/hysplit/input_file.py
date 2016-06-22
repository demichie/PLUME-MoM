hysplit_dir = "/home/demichie/Codes/hysplit"
plumemom_dir = "/home/demichie/Codes/PLUME_MoM"
runname = 'Eya2010'
starttime="10 04 14 00 00" # Year,month,day,hour,minute
endtime = "10 04 14 12 00"
deltat_plumemom = 3600  # seconds

lat = 63.63
lon = -19.62
model_top = 32000.0
meteo_file = 'apr1420.bin'
spacing_lat = 0.02 
spacing_lon = 0.02
span_lat = 30.0
span_lon = 30.0
vent_height = 3000.0
deltaz_release = 500.0

# setup.cfg parameters
kmsl=0
ninit=1
ndump=1
ncycl=1
numpar = 200
maxpar = 50000
# particles parameters
npart = 10
diam1 = 8.0E-6
rho1 = 2000
diam2 = 2.E-3
rho2 = 2000
cp_part = 1610
shapefactor = 1.0
partial_mass_fractions = [3.0E-004, 2.2E-003, 1.2E-002, 4.47E-002, 0.1152 , 
                          0.2031 , 0.2453 , 0.2031 , 0.1152 , 5.9E-002]
diam_phi = [-5.0 , -4.0 , -3.0 , -2.0 , -1.0 , 0.0 , 1.0 , 2.0 , 3.0 , 4.0]



