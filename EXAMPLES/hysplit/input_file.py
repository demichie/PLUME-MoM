hysplit_dir = "/home/demichie/Codes/hysplit/trunk"
plumemom_dir = "/home/demichie/Codes/PLUME_MoM"
runname = 'Calbuco2015_phase_one'
starttime="15 04 22 21 30" # Year,month,day,hour,minute
endemittime = "15 04 22 23 30"
endruntime = "15 04 24 12 00"
deltat_plumemom = 3600  # seconds

lat = -41.33
lon = -72.61
model_top = 32000.0
meteo_file = 'calbuco_wind_emcwf.bin'
spacing_lat = 0.05 
spacing_lon = 0.05
span_lat = 30.0
span_lon = 30.0
vent_height = 2003.0
vent_radius = 20.0
log10_mfr = 6.78
deltaz_release = 200.0

# setup.cfg parameters
kmsl=1  #starting heights default to AGL=0 or MSL=1
ninit=1  #particle initialization(0-none; 1-once; 2-add; 3-replace)
ndump=1  #dump particles to/from file 0-none or nhrs-output intervall
ncycl=1 # pardump output cycle time
numpar = 1000 # number of puffs or particles to released per cycle
maxpar = 5000 # maximum number of particles carried in simulation
initd = 3 # initial distribution, particle, puff, or combination
delt = 30 # hysplit integration step (minutes)

# particles parameters
npart = 10
diam1 = 0.000004
rho1 = 1500
diam2 = 0.002
rho2 = 900
cp_part = 1610
shapefactor = 0.8
partial_mass_fractions = [0.08, 0.15, 0.005, 0.015, 0.06, 0.12, 0.14, 0.24, 0.13, 0.04]
diam_phi = [ -1.0, 0  ,1.0 , 2.0 , 3.0 , 4.0, 5.0, 6.0, 7.0, 8.0]


# CONTROL parameters
#SAMPLING INTERVAL
SI_TYPE = 0 # Avg:0 Now:1 Max:2
SI_HOUR = 6 # hrs
SI_MINUTE = 0 # min
#HEIGHT OF EACH CONCENTRATION LEVEL (m-msl)
H_LEVELS = '0 30000'

#SAMPLING POINTS
P01=[-41.21,-72.53]
P02=[-41.13,-72.41]
