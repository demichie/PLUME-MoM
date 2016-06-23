!********************************************************************
!> \brief Input/Output module
!
!> This module contains all the input/output subroutine and the 
!> realted variables.
!> \date 28/10/2013
!> @author 
!> Mattia de' Michieli Vitturi
!********************************************************************

MODULE inpout

    USE variables

    USE plume_module, ONLY: vent_height, alpha_inp , beta_inp , particles_loss ,&
         r0 , w0 , z , mfr_exp0

    USE particles_module, ONLY: n_part , n_mom , mom0 , rhop_mom

    USE particles_module, ONLY : solid_partial_mass_fraction , diam1 , rho1 ,   &
         diam2 , rho2 , cp_part , settling_model , distribution ,               &
         distribution_variable , solid_mass_fraction , shape_factor

    USE meteo_module, ONLY: gt , gs , p0 , t0 , h1 , h2 , u_atm0 , duatm_dz0 ,  &
         visc_atm0 , rair , cpair , read_atm_profile , u_r , z_r , exp_wind ,   &
         wind_mult_coeff

    USE solver_module, ONLY: ds0

    USE mixture_module, ONLY: tp0 , gas_mass_fraction0 , rwvapour , cpwvapour , &
         initial_neutral_density

  IMPLICIT NONE

  !> Counter for unit files
  INTEGER :: n_unit

  !> Name of input file
  CHARACTER(LEN=30) :: inp_file

  !> Name of output file for backup of input parameters
  CHARACTER(LEN=30) :: bak_file   

  !> Name of the run (used for the output and backup files)
  CHARACTER(LEN=30) :: run_name            

  !> Name of output file for backup of input parameters
  CHARACTER(LEN=30) :: col_file

  !> Name of output file for hysplit
  CHARACTER(LEN=30) :: hy_file

  !> Name of output file for backup of input parameters
  CHARACTER(LEN=30) :: mom_file

  !> Name of output file for the variables used by dakota
  CHARACTER(LEN=30) :: dak_file

  !> Name of output file for the parameters of the beta distribution
  CHARACTER(LEN=30) :: mat_file

  !> Name of file for the parameters of the atmosphere
  CHARACTER(LEN=50) :: atm_file

  !> Atmosphere input unit
  INTEGER :: atm_unit


  !> Backup input unit
  INTEGER :: bak_unit

  !> Beta distribution parameters file unit
  INTEGER :: mat_unit

  !> Input data unit
  INTEGER :: inp_unit

  !> Output values along the column data unit
  INTEGER :: col_unit

  !> hysplit data unit
  INTEGER :: hy_unit

  INTEGER :: hy_lines

  !> hysplit scratch unit
  INTEGER :: temp_unit

  !> Moments values along the column data unit
  INTEGER :: mom_unit

  !> Dakota variables data unit
  INTEGER :: dak_unit

  REAL*8, ALLOCATABLE :: mu_lognormal(:) , sigma_lognormal(:)

  REAL*8 :: month
  REAL*8 :: lat

  REAL*8 :: hy_deltaz , hy_z , hy_z_old , hy_x , hy_y , hy_x_old , hy_y_old 

  REAL*8, ALLOCATABLE :: solid_mfr(:) , solid_mfr_old(:), solid_mfr_init(:) ,   &
        solid_mfr_oldold(:)

  NAMELIST / control_parameters / run_name , verbose_level , dakota_flag ,      &
        inversion_flag , hysplit_flag

  NAMELIST / inversion_parameters / height_weight , height_obj , mu_weight ,    &
       mu_obj , sigma_weight , sigma_obj , skew_weight , skew_obj
  
  NAMELIST / plume_parameters / alpha_inp , beta_inp , particles_loss
  
  NAMELIST / atm_parameters / visc_atm0 , rair , cpair , wind_mult_coeff ,      &
       read_atm_profile , settling_model , shape_factor
  
  NAMELIST / std_atm_parameters / gt , gs , p0 , t0 , h1 , h2 , u_atm0 ,        &
       duatm_dz0 , u_r , z_r , exp_wind
  
  NAMELIST / table_atm_parameters / month , lat , u_r , z_r , exp_wind

  NAMELIST / initial_values / r0 , w0 , mfr_exp0 , tp0 ,                        &
       initial_neutral_density , gas_mass_fraction0 , vent_height , ds0 ,       &
       n_part , distribution , distribution_variable , n_mom
 
  NAMELIST / hysplit_parameters / hy_deltaz , nbl_stop
 
  NAMELIST / mixture_parameters / diam1 , rho1 , diam2 , rho2 , cp_part ,       &
       rwvapour , cpwvapour
  
  NAMELIST / lognormal_parameters / solid_partial_mass_fraction ,               &
       mu_lognormal , sigma_lognormal
  

  SAVE

CONTAINS


  !*****************************************************************************
  !> \brief Initialize variables
  !
  !> This subroutine check if the input file exists and if it does not then it
  !> it initializes the input variables with default values and creates an input
  !> file.
  !> \date 28/10/2013
  !> @author 
  !> Mattia de' Michieli Vitturi
  !******************************************************************************

  SUBROUTINE initialize

    ! External procedures

    USE particles_module, ONLY: allocate_particles

    IMPLICIT NONE

    LOGICAL :: lexist

    gi = 9.81d0               ! Gravity acceleration
    pi_g = 4.D0 * ATAN(1.D0) 

    WIND_MULT_COEFF = 1.D0


    n_unit = 10

    inp_file = 'plume_model.inp'

    INQUIRE (FILE=inp_file,exist=lexist)

    IF (lexist .EQV. .FALSE.) THEN

       !
       !***  Initialization of variables readed in the input file (any version of the
       !***  input file)
       !

       !---------- parameters of the CONTROL_PARAMETERS namelist -------------------
       run_name = 'default_run'
       verbose_level = 0
       dakota_flag = .FALSE.
       hysplit_flag = .FALSE.
       inversion_flag = .FALSE.

       !---------- parameters of the INERSION_PARAMETERS namelist ------------------
       height_weight = 0.D0
       height_obj = 0.D0
       mu_weight = 0.D0
       mu_obj = 0.D0 
       sigma_weight = 0.D0 
       sigma_obj = 0.D0
       skew_weight = 0.D0
       skew_obj = 0.D0

       !---------- parameters of the PLUME_PARAMETERS namelist ---------------------
       alpha_inp = 9.0D-2
       beta_inp = 0.6D0
       particles_loss = .TRUE.

       !---------- parameters of the ATM_PARAMETERS namelist -----------------------
       VISC_ATM0 =  1.8D-5
       RAIR=  287.026
       CPAIR=  998.000000  
       WIND_MULT_COEFF = 1.D0
       READ_ATM_PROFILE = "standard" 
       SETTLING_MODEL = "textor"
       SHAPE_FACTOR = 0.43

       !---------- parameters of the STD_ATM_PARAMETERS namelist -------------------
       gt = -6.5D-3              ! Temp gradient Troposphere
       gs = 1.0D-3               ! Temp gradient Stratosphere
       p0 = 101325.D0            ! Pressure at sea level
       t0 = 288.15D0             ! Temperature at sea level
       h1 = 11.D3
       h2 = 20.D3
       u_atm0 = 0.D0
       duatm_dz0 = 0.D0

       !---------- parameters of the INITIAL_VALUES namelist -----------------------

       R0 = 0.D0 
       W0 = 0.D0
       MFR_exp0 = -1.0
       TP0 = 1273.D0
       INITIAL_NEUTRAL_DENSITY = .FALSE.
       GAS_MASS_FRACTION0 = 3.0D-2
       VENT_HEIGHT =  1500.D0
       DS0 =  5.D0
       N_PART = 1
       DISTRIBUTION = 'lognormal'
       DISTRIBUTION_VARIABLE = 'particles_number'
       N_MOM = 6

       CALL allocate_particles

       !---------- parameters of the MIXTURE_PARAMETERS namelist -------------------

       DIAM1 = 8.D-6
       RHO1 = 2000.D0
       DIAM2 = 2.D-3
       RHO2 = 2600.D0
       CP_PART = 1100.D0
       RWVAPOUR = 462.D0
       CPWVAPOUR = 1810.D0

       !---------- parameters of the LOGNORMAL_PARAMETERS namelist -----------------

       ALLOCATE( mu_lognormal(n_part) )
       ALLOCATE( sigma_lognormal(n_part) )

       SOLID_PARTIAL_MASS_FRACTION =  1.D0
       MU_LOGNORMAL=  2.D0
       SIGMA_LOGNORMAL=  1.6D0

       inp_unit = n_unit

       OPEN(inp_unit,FILE=inp_file,STATUS='NEW')

       WRITE(inp_unit, control_parameters )
       WRITE(inp_unit, plume_parameters )
       WRITE(inp_unit, atm_parameters )
       WRITE(inp_unit, std_atm_parameters )
       WRITE(inp_unit, initial_values )
       WRITE(inp_unit, mixture_parameters )
       WRITE(inp_unit, lognormal_parameters )

       CLOSE(inp_unit)

       WRITE(*,*) 'Input file plume_model.inp not found'
       WRITE(*,*) 'A new one with default values has been created'
       STOP

    END IF

  END SUBROUTINE initialize

  !******************************************************************************
  !> \brief Read Input data 
  !
  !> This subroutine reads input data from the file plume_model.inp and writes a
  !> backup file of the input data.
  !> \date 28/10/2013
  !> @author 
  !> Mattia de' Michieli Vitturi
  !******************************************************************************

  SUBROUTINE read_inp

    ! External variables

    USE meteo_module, ONLY: rho_atm , pa , atm_profile , n_atm_profile

    USE moments_module, ONLY : n_nodes

    USE mixture_module, ONLY: gas_volume_fraction0 , rgasmix

    ! External procedures

    USE meteo_module, ONLY : zmet

    USE meteo_module, ONLY : h_levels

    USE meteo_module, ONLY : rho_atm_month_lat , pres_atm_month_lat ,  temp_atm_month_lat

    USE moments_module, ONLY : beta_function , wheeler_algorithm , coefficient

    USE particles_module, ONLY: particles_density , allocate_particles

    IMPLICIT NONE

    LOGICAL :: tend1
    CHARACTER(LEN=80) :: card

    INTEGER :: i , k , j

    REAL*8, DIMENSION(max_n_part) :: solid_volume_fraction0
    REAL*8, ALLOCATABLE :: d_max(:) 
    REAL*8, ALLOCATABLE :: p_beta(:) , q_beta(:)

    REAL*8, ALLOCATABLE :: mu_bar(:) , sigma_bar(:)
    REAL*8, ALLOCATABLE :: diam_constant(:)
    REAL*8, ALLOCATABLE :: diam_constant_phi(:)

    REAL*8 :: solid_tot_volume_fraction0

    REAL*8 :: C0

    REAL*8, DIMENSION(max_n_part) :: rho_solid_avg

    REAL*8 :: rho_solid_tot_avg

    REAL*8 :: rho_gas
    REAL*8 :: rho_mix

    REAL*8 :: alfa_s

    REAL*8, ALLOCATABLE :: xi(:) , wi(:) 
    REAL*8, ALLOCATABLE :: part_dens_array(:)

    REAL*8, ALLOCATABLE :: atm_profile0(:,:)

    INTEGER :: i_part

    INTEGER*8 :: fact2

    INTEGER, ALLOCATABLE :: coeff(:,:)

    REAL*8, ALLOCATABLE :: rho_atm_month(:,:)

    REAL*8 :: rho_atm_jan(100,13)
    REAL*8 :: rho_atm_apr(100,13)
    REAL*8 :: rho_atm_jul(100,13)
    REAL*8 :: rho_atm_oct(100,13)

    REAL*8, ALLOCATABLE :: pres_atm_month(:,:)

    REAL*8 :: pres_atm_jan(100,13)
    REAL*8 :: pres_atm_apr(100,13)
    REAL*8 :: pres_atm_jul(100,13)
    REAL*8 :: pres_atm_oct(100,13)

    REAL*8, ALLOCATABLE :: temp_atm_month(:,:)

    REAL*8 :: temp_atm_jan(100,13)
    REAL*8 :: temp_atm_apr(100,13)
    REAL*8 :: temp_atm_jul(100,13)
    REAL*8 :: temp_atm_oct(100,13)

    INTEGER :: atm_level

    INTEGER :: n_atm_levels

    REAL*8 :: coeff_lat

    INTEGER :: io

    NAMELIST / beta_parameters / solid_partial_mass_fraction , p_beta , q_beta ,&
         d_max

    NAMELIST / constant_parameters / solid_partial_mass_fraction ,              &
         diam_constant_phi

    WRITE(*,*) 'PLUME_MODEL: *** Starting the run ***' 

    n_unit = n_unit + 1

    inp_unit = n_unit

    inp_file = 'plume_model.inp'

    OPEN(inp_unit,FILE=inp_file,STATUS='old')


    READ(inp_unit, control_parameters)

    n_unit = n_unit + 1
    bak_unit = n_unit
    bak_file = TRIM(run_name)//'.bak'

    OPEN(bak_unit,file=bak_file,status='unknown')


    WRITE(bak_unit, control_parameters)
    
    
    IF ( verbose_level .GE. 1 ) WRITE(*,*) 'read control_parameters: done'

    IF ( inversion_flag ) THEN

       READ(inp_unit, inversion_parameters)
       WRITE(bak_unit, inversion_parameters)

    END IF

    READ(inp_unit, plume_parameters)
    WRITE(bak_unit, plume_parameters)

    IF ( verbose_level .GE. 1 ) WRITE(*,*) 'read plume_parameters: done'

    READ(inp_unit, atm_parameters)
    WRITE(bak_unit, atm_parameters)

    IF ( read_atm_profile .EQ. 'table' ) THEN

       n_atm_levels = 0

       READ( inp_unit, table_atm_parameters )
       WRITE( bak_unit, table_atm_parameters )

       n_unit = n_unit + 1

       atm_unit = n_unit

       atm_file = '../AtmProfile_info/Density_April.txt'

       OPEN(atm_unit,FILE=atm_file,STATUS='old')

       READ(atm_unit,*) 

       atm_level = 0

       atm_read_levels_apr: DO

          atm_level = atm_level + 1
          
          READ(atm_unit,*,IOSTAT=io ) rho_atm_apr(atm_level,1:8)

          IF ( io > 0 ) EXIT

       END DO atm_read_levels_apr

       CLOSE(atm_unit)

       n_unit = n_unit + 1

       atm_unit = n_unit

       atm_file = '../AtmProfile_info/Density_Jan.txt'

       OPEN(atm_unit,FILE=atm_file,STATUS='old')

       READ(atm_unit,*) 

       atm_level = 0

       atm_read_levels_jan: DO

          atm_level = atm_level + 1

          READ(atm_unit,*,IOSTAT=io) rho_atm_jan(atm_level,1:8)

          IF ( io > 0 ) EXIT

       END DO atm_read_levels_jan

       CLOSE(atm_unit)

       n_unit = n_unit + 1

       atm_unit = n_unit

       atm_file = '../AtmProfile_info/Density_July.txt'

       OPEN(atm_unit,FILE=atm_file,STATUS='old')

       READ(atm_unit,*) 

       atm_level = 0

       atm_read_levels_jul: DO

          atm_level = atm_level + 1

          READ(atm_unit,*,IOSTAT=io) rho_atm_jul(atm_level,1:8)

          IF ( io > 0 ) EXIT

       END DO atm_read_levels_jul

       CLOSE(atm_unit)

       n_unit = n_unit + 1

       atm_unit = n_unit

       atm_file = '../AtmProfile_info/Density_Oct.txt'

       OPEN(atm_unit,FILE=atm_file,STATUS='old')

       READ(atm_unit,*) 

       atm_level = 0

       atm_read_levels_oct: DO

          atm_level = atm_level + 1

          READ(atm_unit,*,IOSTAT=io) rho_atm_oct(atm_level,1:8)

          IF ( io > 0 ) EXIT

          n_atm_levels = atm_level

       END DO atm_read_levels_oct

       CLOSE(atm_unit)

       ! ----- READ PRESSURES -------

       n_unit = n_unit + 1

       atm_unit = n_unit

       atm_file = '../AtmProfile_info/Pressure_April.txt'

       OPEN(atm_unit,FILE=atm_file,STATUS='old')

       READ(atm_unit,*) 

       atm_level = 0

       pres_read_levels_apr: DO

          atm_level = atm_level + 1

          READ(atm_unit,*,IOSTAT=io) pres_atm_apr(atm_level,1:8)

          IF ( io > 0 ) EXIT

       END DO pres_read_levels_apr

       CLOSE(atm_unit)

       n_unit = n_unit + 1

       atm_unit = n_unit

       atm_file = '../AtmProfile_info/Pressure_Jan.txt'

       OPEN(atm_unit,FILE=atm_file,STATUS='old')

       READ(atm_unit,*) 

       atm_level = 0

       pres_read_levels_jan: DO
          
          atm_level = atm_level + 1
          
          READ(atm_unit,*,IOSTAT=io) pres_atm_jan(atm_level,1:8)
          
          IF ( io > 0 ) EXIT
          
       END DO pres_read_levels_jan
       
       CLOSE(atm_unit)

       n_unit = n_unit + 1

       atm_unit = n_unit

       atm_file = '../AtmProfile_info/Pressure_July.txt'

       OPEN(atm_unit,FILE=atm_file,STATUS='old')

       READ(atm_unit,*) 

       atm_level = 0

       pres_read_levels_jul: DO

          atm_level = atm_level + 1

          READ(atm_unit,*,IOSTAT=io) pres_atm_jul(atm_level,1:8)

          IF ( io > 0 ) EXIT

       END DO pres_read_levels_jul

       CLOSE(atm_unit)


       n_unit = n_unit + 1

       atm_unit = n_unit

       atm_file = '../AtmProfile_info/Pressure_Oct.txt'

       OPEN(atm_unit,FILE=atm_file,STATUS='old')

       READ(atm_unit,*) 

       atm_level = 0

       pres_read_levels_oct: DO

          atm_level = atm_level + 1

          READ(atm_unit,*,IOSTAT=io) pres_atm_oct(atm_level,1:8)

          IF ( io > 0 ) EXIT

          n_atm_levels = atm_level

       END DO pres_read_levels_oct

       CLOSE(atm_unit)



       ! ----- READ TEMPERATURES -------

       n_unit = n_unit + 1

       atm_unit = n_unit

       atm_file = '../AtmProfile_info/Temp_April.txt'

       OPEN(atm_unit,FILE=atm_file,STATUS='old')

       READ(atm_unit,*) 

       atm_level = 0

       temp_read_levels_apr: DO

          atm_level = atm_level + 1

          READ(atm_unit,*,IOSTAT=io) temp_atm_apr(atm_level,1:8)

          IF ( io > 0 ) EXIT

       END DO temp_read_levels_apr

       CLOSE(atm_unit)

       n_unit = n_unit + 1

       atm_unit = n_unit

       atm_file = '../AtmProfile_info/Temp_Jan.txt'

       OPEN(atm_unit,FILE=atm_file,STATUS='old')

       READ(atm_unit,*) 

       atm_level = 0

       temp_read_levels_jan: DO

          atm_level = atm_level + 1

          READ(atm_unit,*,IOSTAT=io) temp_atm_jan(atm_level,1:8)

          IF ( io > 0 ) EXIT

       END DO temp_read_levels_jan

       CLOSE(atm_unit)

       n_unit = n_unit + 1

       atm_unit = n_unit

       atm_file = '../AtmProfile_info/Temp_July.txt'

       OPEN(atm_unit,FILE=atm_file,STATUS='old')

       READ(atm_unit,*) 

       atm_level = 0

       temp_read_levels_jul: DO

          atm_level = atm_level + 1

          READ(atm_unit,*,IOSTAT=io) temp_atm_jul(atm_level,1:8)

          IF ( io > 0 ) EXIT

       END DO temp_read_levels_jul

       CLOSE(atm_unit)


       n_unit = n_unit + 1

       atm_unit = n_unit

       atm_file = '../AtmProfile_info/Temp_Oct.txt'

       OPEN(atm_unit,FILE=atm_file,STATUS='old')

       READ(atm_unit,*) 

       atm_level = 0

       temp_read_levels_oct: DO

          atm_level = atm_level + 1

          READ(atm_unit,*,IOSTAT=io) temp_atm_oct(atm_level,1:8)

          IF ( io > 0 ) EXIT

          n_atm_levels = atm_level

       END DO temp_read_levels_oct

       CLOSE(atm_unit)

       ALLOCATE( h_levels(n_atm_levels) )

       ALLOCATE( rho_atm_month_lat(n_atm_levels) , rho_atm_month(n_atm_levels,8) )
       ALLOCATE( pres_atm_month_lat(n_atm_levels) , pres_atm_month(n_atm_levels,8) )
       ALLOCATE( temp_atm_month_lat(n_atm_levels) , temp_atm_month(n_atm_levels,8) )

       IF ((month .GE. 0.d0) .and. (month .LE. 1.d0)) THEN
          WRITE(*,*)  'winter'
          rho_atm_month(1:n_atm_levels,1:8) = rho_atm_jan(1:n_atm_levels,1:8)
          pres_atm_month(1:n_atm_levels,1:8) = pres_atm_jan(1:n_atm_levels,1:8)
          temp_atm_month(1:n_atm_levels,1:8) = temp_atm_jan(1:n_atm_levels,1:8)
          
       ELSEIF ((month .GT. 1.d0) .and. (month .LE. 2.d0)) THEN
          WRITE(*,*)  'spring'
          rho_atm_month(1:n_atm_levels,1:8) = rho_atm_apr(1:n_atm_levels,1:8)
          pres_atm_month(1:n_atm_levels,1:8) = pres_atm_apr(1:n_atm_levels,1:8)
          temp_atm_month(1:n_atm_levels,1:8) = temp_atm_apr(1:n_atm_levels,1:8)
          
       ELSEIF ((month .GT. 2.d0) .and. (month .LE. 3.d0)) THEN
          WRITE(*,*)  'summer'
          rho_atm_month(1:n_atm_levels,1:8) = rho_atm_jul(1:n_atm_levels,1:8)
          pres_atm_month(1:n_atm_levels,1:8) = pres_atm_jul(1:n_atm_levels,1:8)
          temp_atm_month(1:n_atm_levels,1:8) = temp_atm_jul(1:n_atm_levels,1:8)
          
       ELSEIF ((month .GT. 3.d0) .and. (month .LE. 4.d0)) THEN
          WRITE(*,*)  'autumn'
          rho_atm_month(1:n_atm_levels,1:8) = rho_atm_apr(1:n_atm_levels,1:8)
          pres_atm_month(1:n_atm_levels,1:8) = pres_atm_apr(1:n_atm_levels,1:8)
          temp_atm_month(1:n_atm_levels,1:8) = temp_atm_apr(1:n_atm_levels,1:8)
          
       END IF

       IF ( ( lat .GE. 0.d0 ) .AND. ( lat .LE. 15.d0 ) ) THEN
          
          coeff_lat = 1.d0 - ( lat - 0.d0 ) / ( 15.d0 - 0.d0 )
          
          rho_atm_month_lat(1:n_atm_levels) = coeff_lat * rho_atm_month(1:n_atm_levels,2) &
               + ( 1.d0 - coeff_lat ) * rho_atm_month(1:n_atm_levels,3)

          pres_atm_month_lat(1:n_atm_levels) = coeff_lat * pres_atm_month(1:n_atm_levels,2) &
               + ( 1.d0 - coeff_lat ) * pres_atm_month(1:n_atm_levels,3)

          temp_atm_month_lat(1:n_atm_levels) = coeff_lat * temp_atm_month(1:n_atm_levels,2) &
               + ( 1.d0 - coeff_lat ) * temp_atm_month(1:n_atm_levels,3)
          
       ELSEIF ( ( lat .GT. 15.d0 ) .AND. ( lat .LE. 30.d0 ) ) THEN
          
          coeff_lat = 1.d0 - ( lat - 15.d0 ) / ( 30.d0 - 15.d0 )
          
          rho_atm_month_lat(1:n_atm_levels) = coeff_lat * rho_atm_month(1:n_atm_levels,3) &
               + ( 1.d0 - coeff_lat ) * rho_atm_month(1:n_atm_levels,4)
          
          pres_atm_month_lat(1:n_atm_levels) = coeff_lat * pres_atm_month(1:n_atm_levels,3) &
               + ( 1.d0 - coeff_lat ) * pres_atm_month(1:n_atm_levels,5)
          
          temp_atm_month_lat(1:n_atm_levels) = coeff_lat * temp_atm_month(1:n_atm_levels,3) &
               + ( 1.d0 - coeff_lat ) * temp_atm_month(1:n_atm_levels,5)
          
       ELSEIF ( ( lat .GT. 30.d0 ) .AND. ( lat .LE. 45.d0 ) ) THEN
          
          coeff_lat = 1.d0 - ( lat - 30.d0 ) / ( 45.d0 - 30.d0 )
          
          rho_atm_month_lat(1:n_atm_levels) = coeff_lat * rho_atm_month(1:n_atm_levels,4) &
               + ( 1.d0 - coeff_lat ) * rho_atm_month(1:n_atm_levels,5)
          
          pres_atm_month_lat(1:n_atm_levels) = coeff_lat * pres_atm_month(1:n_atm_levels,4) &
               + ( 1.d0 - coeff_lat ) * pres_atm_month(1:n_atm_levels,5)
          
          temp_atm_month_lat(1:n_atm_levels) = coeff_lat * temp_atm_month(1:n_atm_levels,4) &
               + ( 1.d0 - coeff_lat ) * temp_atm_month(1:n_atm_levels,5)
          
       ELSEIF ( ( lat .GT. 45.d0 ) .AND. ( lat .LE. 60.d0 ) ) THEN
          
          coeff_lat = 1.d0 - ( lat - 45.d0 ) / ( 60.d0 - 45.d0 )
          
          rho_atm_month_lat(1:n_atm_levels) = coeff_lat * rho_atm_month(1:n_atm_levels,5) &
               + ( 1.d0 - coeff_lat ) * rho_atm_month(1:n_atm_levels,6)
          
          pres_atm_month_lat(1:n_atm_levels) = coeff_lat * pres_atm_month(1:n_atm_levels,5) &
               + ( 1.d0 - coeff_lat ) * pres_atm_month(1:n_atm_levels,6)
          
          temp_atm_month_lat(1:n_atm_levels) = coeff_lat * temp_atm_month(1:n_atm_levels,5) &
               + ( 1.d0 - coeff_lat ) * temp_atm_month(1:n_atm_levels,6)
          
       ELSEIF ( ( lat .GT. 60.d0 ) .AND. ( lat .LE. 75.d0 ) ) THEN
          
          coeff_lat = 1.d0 - ( lat - 60.d0 ) / ( 75.d0 - 60.d0 )
          
          rho_atm_month_lat(1:n_atm_levels) = coeff_lat * rho_atm_month(1:n_atm_levels,6) &
               + ( 1.d0 - coeff_lat ) * rho_atm_month(1:n_atm_levels,7)
          
          pres_atm_month_lat(1:n_atm_levels) = coeff_lat * pres_atm_month(1:n_atm_levels,6) &
               + ( 1.d0 - coeff_lat ) * pres_atm_month(1:n_atm_levels,7)
          
          temp_atm_month_lat(1:n_atm_levels) = coeff_lat * temp_atm_month(1:n_atm_levels,6) &
               + ( 1.d0 - coeff_lat ) * temp_atm_month(1:n_atm_levels,7)
          
       ELSEIF ( ( lat .GT. 75.d0 ) .AND. ( lat .LE. 90.d0 ) ) THEN
          
          coeff_lat = 1.d0 - ( lat - 75.d0 ) / ( 90.d0 - 75.d0 )
          
          rho_atm_month_lat(1:n_atm_levels) = coeff_lat * rho_atm_month(1:n_atm_levels,7) &
               + ( 1.d0 - coeff_lat ) * rho_atm_month(1:n_atm_levels,8)
          
          pres_atm_month_lat(1:n_atm_levels) = coeff_lat * pres_atm_month(1:n_atm_levels,7) &
               + ( 1.d0 - coeff_lat ) * pres_atm_month(1:n_atm_levels,8)
          
          temp_atm_month_lat(1:n_atm_levels) = coeff_lat * temp_atm_month(1:n_atm_levels,7) &
               + ( 1.d0 - coeff_lat ) * temp_atm_month(1:n_atm_levels,8)
          
       END IF
       
       pres_atm_month_lat(1:n_atm_levels) = 100.d0 * pres_atm_month_lat(1:n_atm_levels)

       h_levels(1:n_atm_levels) = 1000.d0 * temp_atm_month(1:n_atm_levels,1)

	!WRITE(*,*) 'rho_atm_month_lat(1:n_atm_levels)', rho_atm_month_lat(1:n_atm_levels)

    ELSEIF ( read_atm_profile .EQ. 'card' ) THEN

       tend1 = .FALSE.

       WRITE(*,*) 'search atm_profile'

       atm_profile_search: DO

          READ(inp_unit,*, END = 200 ) card

          IF( TRIM(card) == 'ATM_PROFILE' ) THEN

             EXIT atm_profile_search

          END IF

       END DO atm_profile_search

       READ(inp_unit,*) n_atm_profile

       IF ( verbose_level .GE. 1 ) WRITE(*,*) 'n_atm_profile',n_atm_profile

       ALLOCATE( atm_profile(7,n_atm_profile) )
       ALLOCATE( atm_profile0(7,n_atm_profile) )

       DO i = 1, n_atm_profile

          READ(inp_unit,*) atm_profile0(1:7,i)


          atm_profile(1:7,i) = atm_profile0(1:7,i)
          ! convert from km to meters
          atm_profile(1,i) = atm_profile(1,i) * 1000.D0

          ! convert from hPa to Pa
          atm_profile(3,i) = atm_profile(3,i) * 100.D0

          atm_profile(6,i) = atm_profile(6,i) * wind_mult_coeff
          atm_profile(7,i) = atm_profile(7,i) * wind_mult_coeff

          IF ( verbose_level .GE. 1 ) WRITE(*,*) i,atm_profile(1,i)

       END DO

       GOTO 210
200    tend1 = .TRUE.
210    CONTINUE

       REWIND(inp_unit)

    ELSEIF ( read_atm_profile .EQ. 'standard' ) THEN

       READ( inp_unit,std_atm_parameters )
       WRITE( bak_unit,std_atm_parameters )

    END IF

    IF ( verbose_level .GE. 1 ) WRITE(*,*) 'read atm_parameters: done'

    READ(inp_unit, initial_values)
    WRITE(bak_unit, initial_values)

    IF ( ( mfr_exp0 .LT. 0.d0 ) .AND. ( r0 .EQ. 0.d0 ) .AND. ( w0 .GT. 0.D0 ) ) THEN
       
       WRITE(*,*) 'WARNING: initial radius calculated from MER and velocity'

    END IF

    IF ( ( mfr_exp0 .LT. 0.d0 ) .AND. ( r0 .EQ. 0.d0 ) .AND. ( w0 .EQ. 0.d0 ) ) THEN
       
       WRITE(*,*) 'WARNING: initial radius and velocity calculated from MER and gas mass fraction'
       STOP

    END IF

    IF ( ( mfr_exp0 .GT. 0.d0 ) .AND. ( w0 .GT. 0.d0 )  .AND. ( r0 .GT. 0.d0 ) ) THEN

       WRITE(*,*) 'ERROR: too many unknown input parameters: input mfr_exp0 or w0 and r0'
       STOP

    END IF

    IF ( distribution .EQ. 'constant' ) THEN

       n_mom = 4
       n_nodes = 1

    ELSE

       n_nodes = NINT(0.5D0 * n_mom)

    END IF

    IF ( hysplit_flag ) THEN

       IF (  distribution .NE. 'constant' ) THEN

          WRITE(*,*) 'hysplit run requires constant distribution'
          STOP
          
       ELSE
          
          READ(inp_unit, hysplit_parameters)
          WRITE(bak_unit, hysplit_parameters)
          
          ALLOCATE( solid_mfr(n_part) , solid_mfr_old(n_part) )
          
          hy_z = vent_height + hy_deltaz
          hy_z_old = vent_height
          hy_x_old = 0.D0
          hy_y_old = 0.D0
          
       END IF
       
    END IF

    z = vent_height

    CALL allocate_particles

    IF ( verbose_level .GE. 1 ) WRITE(*,*) 'read initial_parameters: done'

    READ(inp_unit, mixture_parameters) 
    WRITE(bak_unit, mixture_parameters) 

    IF ( verbose_level .GE. 1 ) WRITE(*,*) 'read mixture_parameters: done'


    IF ( distribution .EQ. 'beta' ) THEN

       ALLOCATE( p_beta(n_part) )
       ALLOCATE( q_beta(n_part) )
       ALLOCATE( d_max(n_part) )

       READ(inp_unit, beta_parameters)
       WRITE(bak_unit, beta_parameters)

    ELSEIF ( distribution .EQ. 'lognormal' ) THEN

       ALLOCATE( mu_lognormal(n_part) )
       ALLOCATE( sigma_lognormal(n_part) )

       ALLOCATE( mu_bar(n_part) )
       ALLOCATE( sigma_bar(n_part) )

       READ(inp_unit, lognormal_parameters)
       WRITE(bak_unit, lognormal_parameters)

       mu_bar = -log( 2.D0 ) * mu_lognormal
       sigma_bar = log( 2.D0 ) * sigma_lognormal

       WRITE(*,*) 'mu_bar',mu_bar
       WRITE(*,*) 'sigma_bar',sigma_bar


    ELSEIF ( distribution .EQ. 'constant' ) THEN

       ALLOCATE( diam_constant(n_part) )
       ALLOCATE( diam_constant_phi(n_part) )

       READ(inp_unit, constant_parameters)
       WRITE(bak_unit, constant_parameters)

       diam_constant = 1.D-3 * 2.D0**(-diam_constant_phi)

    END IF

    IF ( SUM( solid_partial_mass_fraction(1:n_part) ) .NE. 1.D0 ) THEN

       WRITE(*,*) 'WARNING: Sum of solid mass fractions :',                     &
            SUM( solid_partial_mass_fraction(1:n_part) )

    END IF

    solid_partial_mass_fraction(1:n_part)=solid_partial_mass_fraction(1:n_part) &
         / SUM( solid_partial_mass_fraction(1:n_part) )

    solid_mass_fraction(1:n_part) = ( 1.d0 - gas_mass_fraction0 ) *             &
          solid_partial_mass_fraction(1:n_part)

    ALLOCATE( xi(n_nodes) )
    ALLOCATE( wi(n_nodes) )
    ALLOCATE( part_dens_array(n_nodes) )

    ! evaluate the moments from the parameters of the beta distribution
    ! these moments have to be corrected for the mass fractions give in
    ! the input file

    IF ( distribution_variable .EQ. 'mass_fraction' ) THEN
       
       ALLOCATE( coeff(0:n_mom,0:n_mom) )
       CALL coefficient(n_mom,coeff)

    END IF

    DO i_part = 1,n_part

       IF ( verbose_level .GE. 1 ) WRITE(*,*) 'i_part',i_part

       DO i = 0, n_mom-1

          IF ( distribution_variable .EQ. 'particles_number' ) THEN

             IF ( distribution .EQ. 'beta' ) THEN

                mom0(i_part,i) = d_max(i_part)**i *beta_function(p_beta(i_part) &
                     +i,q_beta(i_part)) / beta_function(p_beta(i_part) ,        &
                     q_beta(i_part))

             ELSEIF ( distribution .EQ. 'lognormal' ) THEN

                mom0(i_part,i) = 6.d0 / pi_g * 10.D0**(-3*(i-3)) *              &
                     EXP( ( i-3.D0 ) * mu_bar(i_part) + ( i - 3.D0 )**2 / 2.D0 *&
                     sigma_bar(i_part)**2 ) 

                IF ( verbose_level .GE. 1 ) WRITE(*,*) 'before correction',     &
                     mom0(i_part,i) 

             ELSEIF ( distribution .EQ. 'constant' ) THEN

                mom0(i_part,i) = diam_constant(i_part)**i

             END IF

          ELSEIF ( distribution_variable .EQ. 'mass_fraction' ) THEN

             IF ( distribution .EQ. 'lognormal' ) THEN

                IF ( mu_lognormal(i_part) .EQ. 0.D0) mu_lognormal(i_part) = 1.D-5

                mom0(i_part,i) = 0.d0
                
                DO k=0,floor(0.5D0 * i)
                
                   fact2 = product ((/(j, j = 2*k-1,1,-2)/))

                   mom0(i_part,i) = mom0(i_part,i) + coeff(i,2*k) * fact2       &
                        * sigma_lognormal(i_part)**(2*k) *                      &
                        mu_lognormal(i_part)**(i-2*k)

                END DO

             ELSE

             END IF

          END IF

       END DO

       CALL zmet

       IF ( initial_neutral_density ) THEN
          
          rgasmix = rair
          
       ELSE
          
          rgasmix = rwvapour
          
       END IF
       
       rho_gas = pa / ( rgasmix * tp0 )

       IF ( distribution .EQ. 'constant' ) THEN

          CALL wheeler_algorithm( mom0(i_part,0:1) , distribution , xi , wi )

       ELSE

          CALL wheeler_algorithm( mom0(i_part,0:n_mom-1) , distribution , xi ,  &
               wi )
 
       END IF


       DO i=1,n_nodes

          part_dens_array(i) = particles_density( i_part , xi(i) )

       END DO

       IF ( verbose_level .GE. 1 ) THEN

          WRITE(*,*) 'part_dens_array'
          WRITE(*,*) 'xi',xi
          WRITE(*,*) 'wi',wi
          WRITE(*,*) 'rho',part_dens_array

       END IF

       ! the density of the particles phases are evaluated here. It is 
       ! independent from the mass fraction of the particles phases, so
       ! it is possible to evaluate them with the "uncorrected" moments

       IF ( distribution_variable .EQ. 'particles_number' ) THEN

          rho_solid_avg(i_part) = SUM( part_dens_array * wi * xi**3 ) /         &
               mom0(i_part,3)

       ELSEIF ( distribution_variable .EQ. 'mass_fraction' ) THEN

          rho_solid_avg(i_part) = 1.D0 / ( SUM( wi / part_dens_array ) /        &
               mom0(i_part,0) )

       ELSE

          WRITE(*,*) 'input_file: distribution_variable',distribution_variable
          STOP

       END IF

       IF ( verbose_level .GE. 1 ) THEN

          WRITE(*,*) 'rho avg',rho_solid_avg(i_part)

       END IF

    END DO

    ! the average solid density is evaluated through the mass fractions and 
    ! the densities of the particles phases

    rho_solid_tot_avg = 1.D0 / SUM( solid_partial_mass_fraction(1:n_part) /     &
         rho_solid_avg(1:n_part) )

    IF ( verbose_level .GE. 1 ) THEN

       WRITE(*,*) 
       WRITE(*,*) '******* CHECK ON MASS AND VOLUME FRACTIONS *******'
       WRITE(*,*) 'rho solid avg', rho_solid_tot_avg

    END IF

    IF ( initial_neutral_density ) THEN

       rho_mix = rho_atm

       solid_tot_volume_fraction0 = ( rho_mix - rho_gas ) /                     &
            ( rho_solid_tot_avg - rho_gas )

       gas_volume_fraction0 = 1.D0 - solid_tot_volume_fraction0

       gas_mass_fraction0 =  gas_volume_fraction0 * rho_gas / rho_mix

    ELSE

       gas_volume_fraction0 = rho_solid_tot_avg / ( rho_gas * ( 1.D0 /          &
            gas_mass_fraction0 - 1.D0 ) + rho_solid_tot_avg )

       solid_tot_volume_fraction0 = 1.D0 - gas_volume_fraction0

       rho_mix = gas_volume_fraction0 * rho_gas + solid_tot_volume_fraction0    &
            * rho_solid_tot_avg

    END IF

    IF ( verbose_level .GE. 1 ) THEN
       
       WRITE(*,*) 'gas_volume_fraction0',gas_volume_fraction0
       WRITE(*,*) 'solid_tot_volume_fraction0',solid_tot_volume_fraction0
       WRITE(*,*) 'rho_gas',rho_gas
       WRITE(*,*) 'rho_mix',rho_mix
       WRITE(*,*) 
       
    END IF
    
    DO i_part = 1,n_part

       ! the volume fraction of the particle phases ( with respect to the
       ! solid phase only) is evaluated

       alfa_s = solid_partial_mass_fraction(i_part) * rho_solid_tot_avg /       &
            rho_solid_avg(i_part)

       ! this is the volume fraction of the particles phases in the mixture

       solid_volume_fraction0(i_part) = solid_tot_volume_fraction0 * alfa_s

       ! the coefficient C0 (=mom0) for the particles size distribution is
       ! evaluated in order to have the corrected volume or mass fractions

       ! initialization only
       C0 = 1.D0

       IF ( distribution_variable .EQ. 'particles_number' ) THEN

          C0 = 6.D0 / pi_g * solid_volume_fraction0(i_part) / mom0(i_part,3)
          
       ELSEIF ( distribution_variable .EQ. 'mass_fraction' ) THEN
          
          C0 = ( 1.D0 - gas_mass_fraction0 ) / mom0(i_part,0) *                 &
               solid_partial_mass_fraction(i_part)
      
       END IF
       
       ! the moments are corrected with the factor C0

       DO i = 0, n_mom-1

          mom0(i_part,i) = C0 * mom0(i_part,i)

          IF ( verbose_level .GE. 1 ) WRITE(*,*) 'mom',i,mom0(i_part,i)

       END DO

       IF ( verbose_level .GE. 1 ) THEN

          WRITE(*,*) 'i_part =',i_part 
          WRITE(*,*) 'alfa_s',i_part,alfa_s
          WRITE(*,*) 'solid_volume_fraction0',solid_volume_fraction0(i_part)
          WRITE(*,*) 'solid_partial_mass_fract',                                &
               solid_partial_mass_fraction(i_part)
          WRITE(*,*) 'solid_mass_fract', solid_mass_fraction(i_part)
          WRITE(*,*) 

       END IF

    END DO

    gas_volume_fraction0 = 1.D0 - solid_tot_volume_fraction0

    gas_mass_fraction0 = gas_volume_fraction0 * rho_gas / rho_mix

    WRITE(*,*) 'solid volume fraction',solid_tot_volume_fraction0
    WRITE(*,*) 'solid total mass_fraction', solid_tot_volume_fraction0 *        &
         rho_solid_tot_avg / rho_mix

    IF ( distribution_variable .EQ. 'particles_number' ) THEN

       WRITE(*,*) 'solid_mass_fractions', mom0(1:n_part,3) *                    &
            rho_solid_avg(1:n_part) / ( SUM( mom0(1:n_part,3) *                 &
            rho_solid_avg(1:n_part)) )
          
    ELSEIF ( distribution_variable .EQ. 'mass_fraction' ) THEN

       WRITE(*,*) 'solid_mass_fractions', mom0(1:n_part,0)
       
    END IF

    WRITE(*,*) 'gas volume fraction', gas_volume_fraction0
    WRITE(*,*) 'gas mass fraction', gas_mass_fraction0

    ! the parameters of the particles phases distributions are saved in a file 
    ! readable by Matlab

    IF ( .NOT.dakota_flag ) THEN

       n_unit = n_unit + 1
       
       mat_unit = n_unit
       
       mat_file = TRIM(run_name)//'.m'
       
       OPEN(mat_unit,file=mat_file,status='unknown')
       
       WRITE(mat_unit,*) 'n_part = ',n_part,';'
       WRITE(mat_unit,*) 'gas_volume_fraction = ',gas_volume_fraction0,';'
       
       IF ( distribution .EQ. 'beta' ) THEN
          
          WRITE(mat_unit,*) 'p = [',p_beta(1:n_part),'];'
          WRITE(mat_unit,*) 'q = [',q_beta(1:n_part),'];' 
          WRITE(mat_unit,*) 'd_max = [',d_max(1:n_part),'];'
          
          
       ELSEIF ( distribution .EQ. 'lognormal' ) THEN
          
          WRITE(mat_unit,*) 'mu = [',mu_lognormal(1:n_part),'];'
          WRITE(mat_unit,*) 'sigma = [',sigma_lognormal(1:n_part),'];' 
          
       ELSEIF ( distribution .EQ. 'constant' ) THEN
          
          WRITE(mat_unit,*) 'diam = [',diam_constant(1:n_part),'];'
          
       END IF
       
       WRITE(mat_unit,*) 'solid_mass_fractions = [',                               &
            solid_partial_mass_fraction(1:n_part),'];'
       
       WRITE(mat_unit,*) 'd1 = [',diam1(1:n_part),'];'
       WRITE(mat_unit,*) 'd2 = [',diam2(1:n_part),'];'
       WRITE(mat_unit,*) 'rho1 = [',rho1(1:n_part),'];'
       WRITE(mat_unit,*) 'rho2 = [',rho2(1:n_part),'];'
       
       
       IF ( verbose_level .GE. 0 ) WRITE(*,*) 'write matlab file: done' 
       
       CLOSE(mat_unit)
       
    END IF

    tend1 = .FALSE.
    
    IF ( distribution .EQ. 'moments' ) THEN

       moments_search: DO

          READ(inp_unit , *, END = 300 ) card

          IF( TRIM(card) == 'MOMENTS' ) THEN

             EXIT moments_search

          END IF

       END DO moments_search

       READ(inp_unit,*) n_mom

       WRITE(*,*) 'input_moments'

       READ(inp_unit,*) solid_partial_mass_fraction(1:n_part)

       DO i = 0, n_mom-1

          READ(inp_unit,*) mom0(1:n_part,i)

          WRITE(*,*) mom0(1:n_part,i)

       END DO

       GOTO 310
300    tend1 = .TRUE.
310    CONTINUE

    END IF

    ! Close input file

    CLOSE(inp_unit)

    ! Write a backup of the input file 

    IF ( distribution .EQ. 'moments' ) THEN
       
       IF (( tend1 ) .OR. ( n_mom .EQ. 0 )) THEN

          WRITE(*,*) 'WARNING: input ', ' SAMPLING POINTS not found '

       ELSE

          WRITE(bak_unit,*) '''MOMENTS'''
          WRITE(bak_unit,*) n_mom

           DO i = 0, n_mom-1

             WRITE(bak_unit,*) mom0(1:n_part,i)

          END DO

       END IF

    END IF

    IF ( read_atm_profile .EQ. 'card' ) THEN

       WRITE(bak_unit,*) '''ATM_PROFILE'''
       WRITE(bak_unit,*) n_atm_profile
       
       DO i = 1, n_atm_profile
          
          WRITE(bak_unit,107) atm_profile0(1:7,i)
  
107 FORMAT(7(1x,e14.7))

        
       END DO
       
    END IF

    CLOSE(bak_unit)

    IF ( verbose_level .GE. 1 ) WRITE(*,*) 'end subroutine reainp'

    RETURN

  END SUBROUTINE read_inp

  !*****************************************************************************
  !> \brief Initialize output units
  !
  !> This subroutine set the names of the output files and open the output units
  !> \date 28/10/2013
  !> @author 
  !> Mattia de' Michieli Vitturi
  !******************************************************************************

  SUBROUTINE open_file_units

    ! External variables
    USE particles_module, ONLY : n_part
    USE moments_module, ONLY : n_mom
    USE variables, ONLY : dakota_flag , hysplit_flag

    IMPLICIT NONE
    
    
    col_file = TRIM(run_name)//'.col'
    mom_file = TRIM(run_name)//'.mom'
    dak_file = TRIM(run_name)//'.dak' 
    hy_file = TRIM(run_name)//'.hy'

    IF ( .NOT.dakota_flag ) THEN

       n_unit = n_unit + 1
       col_unit = n_unit
       
       OPEN(col_unit,FILE=col_file)

       n_unit = n_unit + 1
       mom_unit = n_unit
       
       OPEN(mom_unit,FILE=mom_file)
       
       
       WRITE(mom_unit,*) n_part
       WRITE(mom_unit,*) n_mom

    END IF

    IF ( hysplit_flag ) THEN

       n_unit = n_unit + 1
       hy_unit = n_unit
       
       OPEN(hy_unit,FILE=hy_file)

    END IF

    n_unit = n_unit + 1
    dak_unit = n_unit

    OPEN(dak_unit,FILE=dak_file)

    RETURN
    
  END SUBROUTINE open_file_units

  !*****************************************************************************
  !> \brief Close output units
  !
  !> This subroutine close the output units
  !> \date 28/10/2013
  !> @author 
  !> Mattia de' Michieli Vitturi
  !******************************************************************************

  SUBROUTINE close_file_units

    USE variables, ONLY : dakota_flag , hysplit_flag

    IMPLICIT  NONE

    IF ( .not.dakota_flag ) THEN

       CLOSE(col_unit)
       CLOSE(mom_unit)

    END IF

    IF ( hysplit_flag ) CLOSE ( hy_unit )

    CLOSE(dak_unit)

    RETURN

  END SUBROUTINE close_file_units

  !*****************************************************************************
  !> \brief Write outputs
  !
  !> This subroutine writes the output values on the output files. The values
  !> are saved along the column.
  !> \date 28/10/2013
  !> @author 
  !> Mattia de' Michieli Vitturi
  !******************************************************************************

  SUBROUTINE write_column

    USE meteo_module, ONLY: rho_atm , ta, pa

    USE particles_module, ONLY: n_mom , n_part , solid_partial_mass_fraction , &
         mom , set_mom

    USE plume_module, ONLY: x , y , z , w , r , mag_u
    USE mixture_module, ONLY: rho_mix , tp , atm_mass_fraction ,               &
         wvapour_mass_fraction

    ! USE plume_model, ONLY : gas_mass_fraction


    USE variables, ONLY: verbose_level

    IMPLICIT NONE

    REAL*8 :: mfr

    INTEGER :: i_part

    mfr = 3.14 * r**2 * rho_mix * mag_u

    ! WRITE(*,*) 'INPOUT: atm_mass_fraction',atm_mass_fraction
    ! READ(*,*)

    IF ( z .EQ. vent_height ) THEN

       WRITE(col_unit,97,advance="no")
       
       DO i_part=1,n_part
          
          WRITE(col_unit,98,advance="no")
          
       END DO
       
       WRITE(col_unit,99)
       
97     FORMAT(1x,'     z (m)     ',1x,'       r (m)    ',1x,'      x (m)    ',    &
            1x,'     y (m)     ',1x,'mix.dens(kg/m3)',1x,'temperature(C)',         &
            1x,' vert vel (m/s)',1x,' mag vel (m/s) ',1x,' atm.mass fract',         &
            1x,' wv mass fract ',1x)
       
98     FORMAT(1x,'sol.mass fract ')
       
99     FORMAT(1x,'atm.rho(kg/m3)',1x,' MFR (kg/s)     ',1x,'atm.temp (K)  ',         &
            1x,' atm pres (Pa) ')
       

    END IF

    WRITE(col_unit,100) z , r , x , y , rho_mix , tp - 273.15D0 , w , mag_u,&
         atm_mass_fraction , wvapour_mass_fraction ,                        &
         solid_partial_mass_fraction(1:n_part) , rho_atm , mfr , ta, pa

!********* format for plume models comparison ********************************
!
!    WRITE(col_unit,100) z , r , x , y , rho_mix , tp - 273.15D0 , w ,       &
!         atm_mass_fraction , wvapour_mass_fraction ,                        &
!         solid_partial_mass_fraction(1:n_part) * ( 1.D0 - gas_mass_fraction)&
!         , rho_atm , mfr , ta

    WRITE(mom_unit,*) z , mom(1:n_part,0:n_mom-1),set_mom(1:n_part,0)

100 FORMAT(33(1x,e15.8))
    
    IF ( verbose_level .GE. 1 ) THEN
       
       WRITE(*,*) '******************'
       WRITE(*,*) 'z',z
       WRITE(*,*) 'x',x
       WRITE(*,*) 'y',y
       WRITE(*,*) 'r',r
       WRITE(*,*) 'w',w
       WRITE(*,*) '******************'
       READ(*,*)
       
    END IF
    
    RETURN

  END SUBROUTINE write_column

  !*****************************************************************************
  !> \brief Write outputs
  !
  !> This subroutine writes the output values on the output files. The values
  !> are saved along the column.
  !> \date 28/10/2013
  !> @author 
  !> Mattia de' Michieli Vitturi
  !******************************************************************************

  SUBROUTINE write_hysplit(x,y,z,last)

    USE particles_module, ONLY: n_part , solid_partial_mass_fraction

    USE plume_module, ONLY: r , mag_u
    USE mixture_module, ONLY: rho_mix , gas_mass_fraction



    IMPLICIT NONE
    
    REAL*8, INTENT(IN) :: x,y,z

    LOGICAL, INTENT(IN) :: LAST

    INTEGER :: i_part ,i

    CHARACTER(len=8) :: x1 ! format descriptor

    REAL*8 :: xtemp,ytemp,ztemp
    REAL*8 :: zold
    REAL*8 :: znew

    REAL*8, ALLOCATABLE :: solid_temp(:)

    INTEGER :: nbl_lines

    IF ( z .EQ. vent_height ) THEN

       hy_lines = 0

       hy_x_old = x
       hy_y_old = y
       hy_z_old = z

  
       solid_mfr(1:n_part) = solid_partial_mass_fraction(1:n_part) *  ( 1.D0 -   &
            gas_mass_fraction) * pi_g * mag_u * r**2.D0 * rho_mix 

       IF ( nbl_stop ) THEN

          ALLOCATE( solid_mfr_init(1:n_part) ) 
          ALLOCATE( solid_mfr_oldold(1:n_part) ) 

          solid_mfr_init(1:n_part) = solid_mfr(1:n_part)

          n_unit = n_unit + 1
          temp_unit = n_unit
          OPEN(temp_unit, status='SCRATCH' )

       END IF

       WRITE(*,*) 'Solid mass flow rate: ',solid_mfr(1:n_part)

       WRITE(hy_unit,107,advance="no")
       
       DO i_part=1,n_part
          
       WRITE(x1,'(I2.2)') i_part ! converting integer to string using a 'internal file'

          WRITE(hy_unit,108,advance="no") 'S mfr'//trim(x1)//' (kg/s)'
          
       END DO

       WRITE(hy_unit,*) ''
      
    ELSEIF ( z .GE. hy_z ) THEN
       
       hy_lines = hy_lines + 1

       solid_mfr_old(1:n_part) = solid_mfr(1:n_part)
       
       solid_mfr(1:n_part) = solid_partial_mass_fraction(1:n_part) * ( 1.D0 -   &
            gas_mass_fraction) * pi_g * mag_u * r**2.D0 * rho_mix 
       
       WRITE(hy_unit,110) 0.5D0 * ( x +  hy_x_old ) ,                           &
            0.5D0 * ( y +  hy_y_old ) , 0.5D0 * ( z +  hy_z_old ) ,             &
            solid_mfr_old(1:n_part) - solid_mfr(1:n_part)

       hy_x_old = hy_x
       hy_y_old = hy_y
       hy_z_old = hy_z
       
       hy_x = x
       hy_y = y
       hy_z = hy_z + hy_deltaz
       
    END IF
    
    IF ( last ) THEN

       IF ( ( nbl_stop ) .AND. ( z .LT. hy_z_old ) ) THEN

          ALLOCATE( solid_temp(1:n_part) )

          REWIND(hy_unit)
          READ(hy_unit,*) 

          solid_mfr_old(1:n_part) = solid_mfr_init(1:n_part)

          nbl_lines = -1
          
          READ(hy_unit,110), xtemp,ytemp,ztemp,solid_temp(1:n_part)

          nbl_lines = nbl_lines + 1
          
          znew = ztemp - 0.5D0 * hy_deltaz
      
          solid_mfr_oldold(1:n_part) = solid_mfr_old(1:n_part)
          solid_mfr_old(1:n_part) = solid_mfr_old(1:n_part) - solid_temp(1:n_part)
          
          nbl_loop:DO WHILE ( ( znew .LT. z ) .AND. ( nbl_lines + 1 < hy_lines ) )
             
             READ(hy_unit,110), xtemp,ytemp,ztemp,solid_temp(1:n_part)
             
             zold = znew
            
             IF  ( z .GT. ztemp - 0.5D0 * hy_deltaz ) THEN
                
                nbl_lines = nbl_lines + 1
    
                znew = ztemp - 0.5D0 * hy_deltaz
                
                solid_mfr_oldold(1:n_part) = solid_mfr_old(1:n_part)
                solid_mfr_old(1:n_part) = solid_mfr_old(1:n_part) - solid_temp(1:n_part)
             
             ELSE
                
                EXIT nbl_loop
                
             END IF
             
          END DO nbl_loop

          hy_z_old = znew

          solid_mfr_old(1:n_part) = solid_mfr_oldold(1:n_part)

          REWIND(hy_unit)
          READ(hy_unit,*) 

          DO i = 1,nbl_lines

             READ(hy_unit,110), xtemp,ytemp,ztemp,solid_temp(1:n_part)
             WRITE(temp_unit,110), xtemp,ytemp,ztemp,solid_temp(1:n_part)

          END DO

          REWIND(temp_unit)                     ! back to the beginning of SCRATCH

          CLOSE(hy_unit, STATUS = 'DELETE' )   ! delete original

          OPEN(hy_unit,FILE=hy_file)

          WRITE(hy_unit,107,advance="no")
          
          DO i_part=1,n_part
             
             WRITE(x1,'(I2.2)') i_part ! converting integer to string using a 'internal file'
             
             WRITE(hy_unit,108,advance="no") 'S mfr'//trim(x1)//' (kg/s)'
             
          END DO
          
          WRITE(hy_unit,*) ''
       
          DO i = 1,nbl_lines

             READ(temp_unit,110), xtemp,ytemp,ztemp,solid_temp(1:n_part)
             WRITE(hy_unit,110), xtemp,ytemp,ztemp,solid_temp(1:n_part)

          END DO
 
          CLOSE(temp_unit, STATUS = 'DELETE')                      ! delete


       ELSE

          solid_mfr_old(1:n_part) = solid_mfr(1:n_part)
          
       END IF

       solid_mfr(1:n_part) = solid_partial_mass_fraction(1:n_part) * ( 1.D0 -   &
            gas_mass_fraction) * pi_g * mag_u * r**2.D0 * rho_mix 
       
       hy_x = x
       hy_y = y
       hy_z = z
       
       WRITE(hy_unit,110) 0.5D0 * ( hy_x +  hy_x_old ) , 0.5D0 * ( hy_y +       &
            hy_y_old ) , 0.5D0 * ( hy_z +  hy_z_old ) , solid_mfr_old(1:n_part) &
            - solid_mfr(1:n_part)
       
       WRITE(hy_unit,110) hy_x , hy_y , hy_z , solid_mfr(1:n_part)
       
    END IF

107     FORMAT(1x,'     x (m)     ',1x,'      y (m)    ', 1x,'     z (m)     ')
       
108     FORMAT(2x,A)

110 FORMAT(33(1x,e15.8))
 
    RETURN

  END SUBROUTINE write_hysplit

  !*****************************************************************************
  !> \brief Dakota outputs
  !
  !> This subroutine writes the output values used for the sensitivity analysis
  !> by dakota.  
  !> \param[in]    description     descriptor of the output variable
  !> \param[in]    value           value of the output variable
  !> \date 28/10/2013
  !> @author 
  !> Mattia de' Michieli Vitturi
  !******************************************************************************

  SUBROUTINE write_dakota(description,value)
    
    USE variables, ONLY : verbose_level

    IMPLICIT NONE

    CHARACTER(20), INTENT(IN) :: description
    REAL*8, INTENT(IN) :: value

    WRITE(dak_unit,*) description,value
    
    IF ( verbose_level .GE. 2 ) THEN

       WRITE(*,*) description,value
       
    END IF

    RETURN

  END SUBROUTINE write_dakota
  


END MODULE inpout
