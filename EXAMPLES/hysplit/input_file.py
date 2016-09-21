hysplit_dir = "/Users/demichie/Codes/Hysplit4/"
plumemom_dir = "/Users/demichie/Codes/PLUME-MoM/"
runname = 'Etna_23Feb2013'
starttime="13 02 23 18 00" # Year,month,day,hour,minute
endemittime = "13 02 23 19 30"
endruntime = "13 02 25 00 00"
deltat_plumemom = 3600  # seconds

lat = 37.73
lon = 15.00
model_top = 32000.0
meteo_file = '23Feb2013.bin'
spacing_lat = 0.05 
spacing_lon = 0.05
span_lat = 30.0
span_lon = 30.0
vent_height = 3300.0
vent_radius = 20.0
log10_mfr = 6.01
deltaz_release = 200.0

# setup.cfg parameters
kmsl=0  #starting heights default to AGL=0 or MSL=1
ninit=1  #particle initialization(0-none; 1-once; 2-add; 3-replace)
ndump=1  #dump particles to/from file 0-none or nhrs-output intervall
ncycl=1 # pardump output cycle time
numpar = 1000 # number of puffs or particles to released per cycle
maxpar = 5000 # maximum number of particles carried in simulation
initd = 3 # initial distribution, particle, puff, or combination
delt = 5 # hysplit integration step (minutes)

# particles parameters
npart = 12
diam1 = 0.0000078
rho1 = 2500
diam2 = 0.001
rho2 = 1000
cp_part = 1610
shapefactor = 0.8
partial_mass_fractions = [0.03, 0.22, 0.28, 0.20, 0.12, 0.07, 0.04, 0.02, 0.01, 0.005, 0.003, 0.002]
diam_phi = [ -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6]


# CONTROL parameters
#SAMPLING INTERVAL
SI_TYPE = 0 # Avg:0 Now:1 Max:2
SI_HOUR = 1 # hrs
SI_MINUTE = 0 # min
#HEIGHT OF EACH CONCENTRATION LEVEL (m-msl)
H_LEVELS = '0 30000'

#SAMPLING POINTS
P01=[38.19, 15.55]
P02=[38.65, 16.38]
P03=[40.63, 17.93]
