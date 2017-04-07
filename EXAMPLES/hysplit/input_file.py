hysplit_dir = "/home/demichie/Codes/hysplit/trunk"
plumemom_dir = "/home/demichie/Codes/PLUME-MoM-master"
runname = 'Etna_23Nov2013'
starttime="13 11 23 09 00"      # Year,month,day,hour,minute
endemittime = "13 11 23 10 15"  # Year,month,day,hour,minute
endruntime = "13 11 23 16 00"   # Year,month,day,hour,minute
deltat_plumemom = 3600          # seconds

lat = 37.73
lon = 15.00
model_top = 32000.0             # meters
meteo_file = 'INGVC_2311.arl'
spacing_lat = 0.01              # degrees 
spacing_lon = 0.01              # degrees
span_lat = 10.0                 # degrees 
span_lon = 10.0                 # degrees
vent_height = 3300.0            # meters
vent_velocity = 100.0           # m/s
log10_mfr = [6.6, 6.7]          # log10 of kg/s
gas_mass_fraction = 0.03        # wt. fraction 
deltaz_release = 200.0          # meters

# setup.cfg parameters
kmsl=0                          # starting heights default to AGL=0 or MSL=1
ninit=1                         # particle initialization(0-none; 1-once; 2-add; 3-replace)
ndump=1                         # dump particles to/from file 0-none or nhrs-output intervall
ncycl=1                         # pardump output cycle time
numpar = 1000                   # number of puffs or particles to released per cycle
maxpar = 50000                  # maximum number of particles carried in simulation
initd = 3                       # initial distribution, particle, puff, or combination
delt = 2                        # hysplit integration step (minutes)

# particles parameters
npart = 8                       # particle classes
diam1 = 0.0000078               # diameter in meters for density function (Bonadonna and Phillips)
rho1 = 2500                     # density in kg/m3
diam2 = 0.001                   # diameter in meters
rho2 = 1000                     # density in kg/m3
cp_part = 1610                  # J / ( kg K )
shapefactor = 0.8               # for Ganser model
partial_mass_fractions = [0.02, 0.125, 0.32, 0.28, 0.18, 0.06, 0.01, 0.005]
diam_phi = [ -5, -4, -3, -2, -1, 0, 1, 2]


# CONTROL parameters
#SAMPLING INTERVAL
SI_TYPE = 0                     # Avg:0 Now:1 Max:2
SI_HOUR = 1                     # hrs
SI_MINUTE = 0                   # min
#HEIGHT OF EACH CONCENTRATION LEVEL (m-msl)
H_LEVELS = '0 30000'

#SAMPLING POINTS
P01=[37.76, 15.05]
P02=[37.81, 15.17]
P03=[37.83, 15.26]
