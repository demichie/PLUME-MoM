hysplit_dir = "/home/demichie/Codes/hysplit/trunk"
plumemom_dir = "/home/demichie/Codes/PLUME_MoM"
runname = 'Calbuco2015'
starttime="15 04 22 21 00" # Year,month,day,hour,minute
endemittime = "15 04 22 22 30"
endruntime = "15 04 24 00 00"
deltat_plumemom = 3600  # seconds

lat = -41.33
lon = -72.61
model_top = 32000.0
meteo_file = 'wind_GFSmodel.bin'
spacing_lat = 0.05 
spacing_lon = 0.05
span_lat = 30.0
span_lon = 30.0
vent_height = 2003.0
vent_radius = 20.0
log10_mfr = 7.0
deltaz_release = 200.0

# setup.cfg parameters
kmsl=0  #starting heights default to AGL=0 or MSL=1
ninit=1  #particle initialization(0-none; 1-once; 2-add; 3-replace)
ndump=1  #dump particles to/from file 0-none or nhrs-output intervall
ncycl=1 # pardump output cycle time
numpar = 10 # number of puffs or particles to released per cycle
maxpar = 50000 # maximum number of particles carried in simulation
initd = 3 # initial distribution, particle, puff, or combination
delt = 30 # hysplit integration step (minutes)

# particles parameters
npart = 10
diam1 = 0.0262
rho1 = 760
diam2 = 0.024
rho2 = 1670
cp_part = 1610
shapefactor = 1.0
partial_mass_fractions = [3.0E-004, 2.2E-003, 1.2E-002, 4.47E-002, 0.1152 , 
                          0.2031 , 0.2453 , 0.2031 , 0.1152 , 5.9E-002]
diam_phi = [-5.0 , -4.0 , -3.0 , -2.0 , -1.0 , 0.0 , 1.0 , 2.0 , 3.0 , 4.0]


# CONTROL parameters
#SAMPLING INTERVAL
SI_TYPE = 0
SI_HOUR = 6
SI_MINUTE = 0
#HEIGHT OF EACH CONCENTRATION LEVEL (m-msl)
H_LEVELS = '0 15000'

#SAMPLING POINTS
P01=[-41.00,-65.00]
P02=[-41.00,-67.50]
P03=[-41.00,-70.00]
