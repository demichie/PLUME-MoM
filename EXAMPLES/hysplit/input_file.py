hysplit_dir = "/home/demichie/Codes/hysplit/trunk"
plumemom_dir = "/home/demichie/Codes/PLUME-MoM-multigas"
runname = 'etna_test'
starttime="15 12 03 10 00" # Year,month,day,hour,minute
endemittime = "15 12 03 14 00"
endruntime = "15 12 03 17 00"
deltat_plumemom = 3600  # seconds

lat = 37.73   # center latitude of the grid
lon = 15.00  # center longitude of the grid
model_top = 32000.0
meteo_file = 'extract_22666.bin'

spacing_lat = 0.01 # degrees between nodes of the sampling grid
spacing_lon = 0.01 # degrees between nodes of the sampling grid
span_lat = 5.00   # the total span of the grid in x direction. For instance, a span of 10 degrees would cover 5 degrees on each side of the center grid location
span_lon = 5.00   # the total span of the grid in y direction. For instance, a span of 10 degrees would cover 5 degrees on each side of the center grid location


vent_lat = 37.73  	# vent latitude
vent_lon = 15.00       # vent longitude
vent_height = 3300.00    # vent height above sea level (it can be different from ground level of meteo data at vent lat,lon)
vent_velocity = 100.0
log10_mfr = 5.10

# volcanic gas parameters
ngas = 2   # in addition to H2O
rvolcgas = [189, 130 ] # CO2 and SO2 R constant [J/kgK]
cpvolcgas = [844, 640]
volcgas_mol_wt = [0.044, 0.064]
volcgas_mass_fraction = [0.01, 0.01]

#initial water mass fraction
water_mass_fraction0 = 0.03

# hysplit parameters
deltaz_release = 200.0
ncloud = 1

# setup.cfg parameters
kmsl=1  	# starting heights default to AGL=0 or MSL=1
ninit=1  	# particle initialization(0-none; 1-once; 2-add; 3-replace)
ndump=1  	# dump particles to/from file 0-none or nhrs-output intervall
ncycl=1 	# pardump output cycle time
numpar = 1000 	# number of puffs or particles to released per cycle
maxpar = 30000 # maximum number of particles carried in simulation
initd = 3 	# initial distribution, particle, puff, or combination
delt = 10 	# hysplit integration step (minutes)
pinpf = '' 

# CONTROL parameters
#SAMPLING INTERVAL
SI_TYPE = 1 # Avg:0 Now:1 Max:2
SI_HOUR = 1 # hrs
SI_MINUTE = 0 # min
#HEIGHT OF EACH CONCENTRATION LEVEL (m-msl)
H_LEVELS = '0 30000'



# particles parameters
npart = 9
diam1 = 0.000004
rho1 = 2200
diam2 = 0.002
rho2 = 1800
cp_part = 1610
shapefactor = 0.6




# Fuego (subplinian)
#partial_mass_fractions = [	0.07, 	0.12, 	0.27,	0.265, 	0.16, 	0.045,	0.02, 	0.025,	0.0125,	0.0075,	0.005 ]
#diam_phi = 		 [ 	-2 ,	-1 , 	0 ,	1 , 	2 , 	3 , 	4 , 	5 , 	6 , 	7 , 	8 ]	

# Ruapehu (subplinian)
#partial_mass_fractions = [	0.05, 	0.09, 	0.14 ,	0.13, 	0.155,	0.15,	0.10, 	0.08,	0.03,	0.02,	0.01,	0.005,	0.0025,	0.0013,	0.0012 	]
#diam_phi = 		 [ 	-4 ,	-3 , 	-2 ,	-1 , 	0 , 	1 , 	2 , 	3 , 	4 , 	5 , 	6 ,	7 ,	8,	9,	10	]	

# Stromboli (violent strombolian, 2003)
#partial_mass_fractions = [	0.05, 	0.05, 	0.03,	0.10, 	0.30,	0.255,	0.115, 	0.045,	0.055]
#diam_phi = 		 [ 	-4 ,	-3 , 	-2 ,	-1 , 	0 , 	1 , 	2 , 	3 , 	4 ]	

# Cordon Caulle (subplininan, 2011)
#partial_mass_fractions = [	0.13, 	0.155, 	0.164,	0.146, 	0.093,	0.093,	0.018, 	0.014,	0.031,	0.029,	0.022,	0.015,	0.012,	0.007,	0.005	]
#diam_phi = 		 [ 	-4 ,	-3 , 	-2 ,	-1 , 	0 , 	1 , 	2 , 	3 , 	4 , 	5 , 	6 ,	7 ,	8,	9,	10	]	

# Eyjafjallajokull (VEI4, 2010)
#partial_mass_fractions = [	0.055,	0.095, 	0.10,	0.16,	0.13, 	0.11,	0.09,	0.085,	0.065,	0.05,	0.03,	0.01,	0.02	]
#diam_phi = 		 [ 	-2 ,	-1 , 	0 , 	1 , 	2 , 	3 , 	4 , 	5 , 	6 ,	7 ,	8,	9,	10	]	

# Vesuvio (subplinian, 1631)
#partial_mass_fractions = [	0.04, 	0.07, 	0.06,	0.06, 	0.12,	0.17,	0.10, 	0.06,	0.07,	0.09,	0.09,	0.05,	0.01	]
#diam_phi = 		 [ 	-4 ,	-3 , 	-2 ,	-1 , 	0 , 	1 , 	2 , 	3 , 	4 , 	5 , 	6 ,	7 ,	8,	]	

# Etna (violent strmobolian, 2002)
partial_mass_fractions = [	0.020, 	0.090, 	0.140,	0.240, 	0.195,	0.170,	0.085, 	0.035,	0.025	]
diam_phi = 		 [ 	-3 , 	-2 ,	-1 , 	0 , 	1 , 	2 , 	3 , 	4 , 	5	]	

# Vesuvio (plinian, 79)
#partial_mass_fractions = [	0.055,	0.075, 	0.12, 	0.125,	0.13, 	0.15,	0.095,	0.08, 	0.045,	0.135	]
#diam_phi = 		 [ 	-4 ,	-3 , 	-2 ,	-1 , 	0 , 	1 , 	2 , 	3 , 	4 , 	5 	]	




#SAMPLING POINTS
#P01=[-38.97,-67.82]
#P02=[-38.70,-68.02]
#P03=[-38.96,-68.04]
#P04=[-38.96,-68.06]
#P05=[-38.94,-68.11]
#P06=[-38.96,-68.05]
#P07=[-38.95,-68.23]
#P08=[-39.02,-67.58]
#P09=[-39.04,-67.57]


