!********************************************************************************
!> \brief Particles module
!
!> This module contains the procedures and the variables related to the solid
!> particles. In particular, the statistical moments of the properties of the 
!> particles are defined and evaluated in this module.
!> \date 22/10/2013
!> @author 
!> Mattia de' Michieli Vitturi
!********************************************************************************
MODULE particles_module
  !
  USE moments_module, ONLY : n_mom , n_nodes

  IMPLICIT NONE

  INTEGER :: n_part

  !> mass fraction of the particle phases with respect to the total solid
  REAL*8, ALLOCATABLE, DIMENSION(:) :: solid_partial_mass_fraction

  !> volume fraction of the particle phases with respect to the total solid
  REAL*8, ALLOCATABLE, DIMENSION(:) :: solid_partial_volume_fraction
  
  !> mass fraction of the particle phases with respect to the mixture
  REAL*8, ALLOCATABLE, DIMENSION(:) :: solid_mass_fraction

  !> volume fraction of the particle phases with respect to the mixture
  REAL*8, ALLOCATABLE, DIMENSION(:) :: solid_volume_fraction
    
  !> Moments of the particles diameter
  REAL*8, ALLOCATABLE, DIMENSION(:,:) :: mom

  !> Moments of the settling velocities
  REAL*8, ALLOCATABLE, DIMENSION(:,:) :: set_mom

  !> Moments of the densities
  REAL*8, ALLOCATABLE, DIMENSION(:,:) :: rhop_mom

  !> Moments of the densities times the settling velocities
  REAL*8, ALLOCATABLE, DIMENSION(:,:) :: set_rhop_mom

  !> Moments of the heat capacities times the densities
  REAL*8, ALLOCATABLE, DIMENSION(:,:) :: cp_rhop_mom 

  !> Moments of the settling velocities times the heat cap times the densities
  REAL*8, ALLOCATABLE, DIMENSION(:,:) :: set_cp_rhop_mom

  !> Moments of the settling velocities times the heat capacity
  REAL*8, ALLOCATABLE, DIMENSION(:,:) :: set_cp_mom

  !> Term accounting for the birth of aggregates in the moments equations
  REAL*8, ALLOCATABLE, DIMENSION(:,:) :: birth_term

  !> Term accounting for the loss of particles because of aggregation 
  REAL*8, ALLOCATABLE, DIMENSION(:,:) :: death_term

  !> shape factor for settling velocity (Pfeiffer)
  REAL*8 :: shape_factor

  !> First diameter for the density function
  REAL*8, ALLOCATABLE :: diam1(:)

  !> Density at d=diam1
  REAL*8, ALLOCATABLE :: rho1(:)

  !> Second diameter for the density function
  REAL*8, ALLOCATABLE :: diam2(:)

  !> Density at d=diam1
  REAL*8, ALLOCATABLE :: rho2(:)

  REAL*8, ALLOCATABLE :: cp_part(:)

  !> Initial (at the base) moments of the particles diameter
  REAL*8, DIMENSION(1:50,0:100) :: mom0

  !> Settling model:\n
  !> - 'textor'    => Textor et al. 2006
  !> - 'pfeiffer'  => Pfeiffer et al. 2005
  !> .
  CHARACTER*10 :: settling_model

  !> Ditribution of the particles:\n
  !> - 'beta'      => beta distribution
  !> - 'lognormal' => lognormal distribution
  !> - 'constant'  => 
  !> .
  CHARACTER(LEN=20) :: distribution

  CHARACTER(LEN=20) :: distribution_variable

  !> Flag for the aggregation:\n
  !> - 'TRUE'   => aggregation enabled
  !> - 'FALSE'  => aggregation disabled
  LOGICAL :: aggregation

  SAVE

CONTAINS

  !******************************************************************************
  !> \brief Particles variables inizialization
  !
  !> This subroutine allocate and evaluate the variables defining the moments for  
  !> the particles. The moments are then corrected, if needed, and then the 
  !> abscissas and weights for the quadrature formulas are computed.
  !> \date 22/10/2013
  !> @author 
  !> Mattia de' Michieli Vitturi
  !******************************************************************************
 
  SUBROUTINE allocate_particles

    IMPLICIT NONE

    ALLOCATE ( solid_partial_mass_fraction(1:n_part) )
    ALLOCATE ( solid_partial_volume_fraction(1:n_part) )
    ALLOCATE ( solid_mass_fraction(1:n_part) )
    ALLOCATE ( solid_volume_fraction(1:n_part) )

    ! Allocation of the arrays for the moments
    ALLOCATE ( mom(1:n_part,0:n_mom-1) )
    ALLOCATE ( set_mom(1:n_part,0:n_mom-1) )
    ALLOCATE ( rhop_mom(1:n_part,0:n_mom-1) )
    ALLOCATE ( set_rhop_mom(1:n_part,0:n_mom-1) )
    ALLOCATE ( set_cp_rhop_mom(1:n_part,0:n_mom-1) )
    ALLOCATE ( set_cp_mom(1:n_part,0:n_mom-1) )
    ALLOCATE ( cp_rhop_mom(1:n_part,0:n_mom-1) )
    ALLOCATE ( birth_term(1:n_part,0:n_mom-1) )
    ALLOCATE ( death_term(1:n_part,0:n_mom-1) )

    ! Allocation of the parameters for the variable density
    ALLOCATE ( diam1(n_part) )
    ALLOCATE ( rho1(n_part) )
    ALLOCATE ( diam2(n_part) )
    ALLOCATE ( rho2(n_part) )
    
    ALLOCATE ( cp_part(n_part) )

  END SUBROUTINE allocate_particles

  !******************************************************************************
  !> \brief Particles variables inizialization
  !
  !> This subroutine allocate and evaluate the variables defining the moments for  
  !> the particles. The moments are then corrected, if needed, and then the 
  !> abscissas and weights for the quadrature formulas are computed.
  !> \date 22/10/2013
  !> @author 
  !> Mattia de' Michieli Vitturi
  !******************************************************************************
 
  SUBROUTINE initialize_particles

    USE moments_module, ONLY : wheeler_algorithm , moments_correction
    USE variables, ONLY : verbose_level

    IMPLICIT NONE

    REAL*8, DIMENSION(n_part,n_nodes) :: xi
    REAL*8, DIMENSION(n_part,n_nodes) :: w

    INTEGER :: i_part

    DO i_part=1,n_part
 
       mom(i_part,:) = mom0(i_part,0:n_mom-1)

       ! CALL moments_correction(mom(i_part,:),iter)

       IF ( distribution .EQ. 'constant' ) THEN

          CALL wheeler_algorithm( mom(i_part,0:1) , distribution , xi(i_part,:),&
               w(i_part,:) )

       ELSE

          CALL wheeler_algorithm( mom(i_part,:) , distribution , xi(i_part,:) , &
               w(i_part,:) )

       END IF

       IF ( verbose_level .GE. 1 ) THEN

          WRITE(*,*) 'part ',i_part
          WRITE(*,*) 'mom',mom(i_part,:)
          WRITE(*,*) 'abscissas',xi(i_part,:)
          WRITE(*,*) 'weights',w(i_part,:)
          READ(*,*)

       END IF

          
    END DO

    CALL eval_particles_moments( xi(1:n_part,:) , w(1:n_part,:) ) 


  END SUBROUTINE initialize_particles


  !******************************************************************************
  !> \brief Settling velocity
  !
  !> This function evaluates the settling velocity of a particle given the size
  !> (diameter), using the expression given in Textor et al. 2006 or in Pfeiffer 
  !> et al 2005, accordingly with the variable SETTLING MODEL specified in the
  !> input file.
  !> \param[in]   i_part   particle phase index 
  !> \param[in]   diam_in  particle diameter (m or phi)
  !> \date 22/10/2013
  !> @author 
  !> Mattia de' Michieli Vitturi
  !******************************************************************************

  FUNCTION particles_settling_velocity(i_part,diam_in)
    !
    USE meteo_module, ONLY : rho_atm , rho_atm0 , visc_atm

    USE variables, ONLY : gi , pi_g , verbose_level

    IMPLICIT NONE
    
    REAL*8 :: particles_settling_velocity

    INTEGER, INTENT(IN) :: i_part
    REAL*8, INTENT(IN) :: diam_in
    
    REAL*8 :: diam

    REAL*8 :: rhop
    
    REAL*8 :: k1 , k2 , k3
    REAL*8 :: CD

    !> cross sectional area
    REAL*8 :: A_cs

    !> Drag coefficients at Rey=100,1000
    REAL*8 :: Cd_100 , Cd_1000

    !> Drag coefficent for intermediate values of Re
    REAL*8 :: Cd_interp

    !> Settling velocity at Rey=100,1000
    REAL*8 :: Us_100 , Us_1000

    !> Mass of the particle
    REAL*8 :: mass

    !> Settling velocities
    REAL*8 :: Us , Us_1 ,Us_2

    !> Reynolds numbers for the two solutions of the settling equation
    REAL*8 :: Rey1 , Rey2

    !> Coefficients of the settling equation
    REAL*8 :: c0 , c1 , c2

    !> Square root of the discriminant of the settling equation
    REAL*8 :: sqrt_delta

    IF ( distribution_variable .EQ. 'mass_fraction' ) THEN

       diam = 1.D-3 * 2.D0 ** ( - diam_in )

    ELSE

       diam = diam_in
       
    END IF

    rhop = particles_density(i_part,diam_in)
    
    IF ( settling_model .EQ. 'textor' ) THEN

    ! Textor et al. 2006
 
    IF ( diam .LE. 1.D-4 ) THEN

       k1 = 1.19D5   ! (m^2 kg^-1 s^-1 )

       particles_settling_velocity = k1 * rhop * SQRT( rho_atm0 / rho_atm ) *   &
            ( 0.5D0 * diam )**2

    ELSEIF ( diam .LE. 1.D-3 ) THEN

       k2 = 8.D0    ! (m^3 kg^-1 s^-1 )

       particles_settling_velocity = k2 * rhop * SQRT( rho_atm0 / rho_atm ) *   &
            ( 0.5D0 * diam )

    ELSE 
    
       k3 = 4.833D0 ! (m^2 kg^-0.5 s^-1 )
       CD = 0.75D0

       particles_settling_velocity = k3 * SQRT( rhop / CD ) * SQRT(  rho_atm0   &
            / rho_atm ) * SQRT( 0.5D0 * diam )

    END IF

    ELSEIF ( settling_model .EQ. 'pfeiffer' ) THEN

       k1 = shape_factor**(-0.828)
       k2 = 2.D0 * DSQRT( 1.07 - shape_factor )

       mass = rhop * 4.D0/3.D0 * pi_g * ( 0.5*diam )**3

       A_cs = pi_g * ( 0.5*diam )**2

       c0 = -2.D0 * diam * mass * gi
       c1 = 24.D0 * visc_atm * k1 * A_cs
       c2 = rho_atm * diam * k2 * A_cs

       sqrt_delta = sqrt( c1**2 - 4 * c0*c2 )

       Us_1 = ( - c1 + sqrt_delta ) / ( 2 * c2 )
       Us_2 = ( - c1 - sqrt_delta ) / ( 2 * c2 )


       Cd_100 = 24.D0/100.D0 * k1 + k2
       Us_100 = sqrt( 2 * mass * gi / ( Cd_100*rho_atm * A_cs ) )

       Cd_1000 = 1.D0
       Us_1000 = sqrt( 2 * mass * gi / ( Cd_1000*rho_atm * A_cs ) )

       Rey1 = rho_atm * diam * Us_1 / visc_atm
       Rey2 = rho_atm * diam * Us_2 / visc_atm

       IF ( verbose_level .GE. 4 ) THEN

          WRITE(*,*) 'rho_atm , diam , Us_1 , visc_atm',rho_atm , diam , Us_1 , &
               visc_atm
          WRITE(*,*) 'Rey1,Rey2',Rey1,Rey2
          READ(*,*)

       END IF

       ! Initialization only
       Us = Us_1000

       IF ( ( Rey1 .GT. 0.D0 ) .AND. ( Rey1 .LE. 100.D0 ) ) THEN
    
          ! For small Reynolds numbers the drag coefficient is given by Eq.8
          ! of Pfeiffer et al. 2005 and the settling velocity is Us_1

          Us = Us_1  
             
       ELSEIF ( ( Rey1 .GT. 100.D0 ) .AND. ( Rey1 .LE. 1000.D0 ) ) THEN
    
          ! For intermediate Reyonlds numbers, 100<Re<1000, the drag coefficient 
          ! is linearly interpolated between Cd_100 and Cd_1000

          Cd_interp = Cd_100 + ( Rey1 - 100 ) / ( 1000 - 100 ) *                &
               ( Cd_1000 - Cd_100)
          Us = sqrt( 2 * mass * gi / ( Cd_interp * rho_atm * A_cs ) )

       ELSEIF ( Rey1 .GT. 1000.D0 ) THEN
    
          ! For large Reynolds numbers the drag coefficient is taken as Cd=1,
          ! as in Pfeiffer et al. 2005 with the settling velocity is Us_1000

          Us = Us_1000
       
       END IF

       IF ( ( Rey2 .GT. 0.D0 ) .AND. ( Rey2 .LE. 100.D0 ) ) THEN 
    
          Us = Us_2
    
       ELSEIF ( ( Rey2 .GT. 100.D0 ) .AND. ( Rey2 .LE. 1000.D0 ) ) THEN 
    
          Cd_interp = Cd_100 + ( Rey2 - 100 ) / ( 1000 - 100 ) * ( Cd_1000 - Cd_100)
          Us = sqrt( 2 * mass * gi / ( Cd_interp * rho_atm * A_cs ) )

       ELSEIF ( Rey2 .GT. 1000.D0 ) THEN
    
          Us = Us_1000
       
       END IF

       particles_settling_velocity = Us

    ELSE

       WRITE(*,*) 'wrong settling model'
       STOP

    END IF

  END FUNCTION particles_settling_velocity

  !******************************************************************************
  !> \brief Heat capacity
  !
  !> This function evaluates the heat capacity of the particles given the size
  !> (diameter). 
  !> \param[in]   i_part   particle phase index 
  !> \param[in]   diam_in  particle diameter (m or phi)
  !> \date 22/10/2013
  !> @author 
  !> Mattia de' Michieli Vitturi
  !******************************************************************************

  FUNCTION particles_heat_capacity(i_part,diam_in)
    !
    IMPLICIT NONE

    REAL*8 :: particles_heat_capacity
    INTEGER, INTENT(IN) :: i_part
    REAL*8, INTENT(IN) :: diam_in
    REAL*8 :: diam

    IF ( distribution_variable .EQ. 'mass_fraction' ) THEN

       diam = 1.D-3 * 2.D0 ** ( - diam_in )

    ELSE

       diam = diam_in

    END IF

    particles_heat_capacity = cp_part(i_part)

  END FUNCTION particles_heat_capacity

  !******************************************************************************
  !> \brief Particle density
  !
  !> This function evaluates the density of a particle given the size (diameter),
  !> using the expression given in Bonadonna and Phillips, 2003.
  !> \param[in]   i_part   particle phase index 
  !> \param[in]   diam_in  particle diameter (m or phi)
  !> \date 22/10/2013
  !> @author 
  !> Mattia de' Michieli Vitturi
  !******************************************************************************

  FUNCTION particles_density(i_part,diam_in)
    !
    IMPLICIT NONE

    REAL*8 :: particles_density

    INTEGER, INTENT(IN) :: i_part
    REAL*8, INTENT(IN) :: diam_in

    REAL*8 :: diam

    REAL*8 :: diam_phi , diam1_phi , diam2_phi


    IF ( distribution_variable .EQ. 'mass_fraction' ) THEN

       diam = 1.D-3 * 2.D0 ** ( - diam_in )

    ELSE

       diam = diam_in

    END IF

    IF ( diam .LE. diam1(i_part) ) THEN

       particles_density = rho1(i_part)

    ELSEIF ( diam .LE. diam2(i_part) ) THEN

       diam_phi = -log(diam*1000)/log(2.D0)
       diam1_phi = -log(diam1(i_part)*1000)/log(2.D0)
       diam2_phi = -log(diam2(i_part)*1000)/log(2.D0)

       particles_density = rho1(i_part) + ( diam_phi - diam1_phi ) /            &
            ( diam2_phi - diam1_phi ) * ( rho2(i_part) - rho1(i_part) )
       
    ELSE

       particles_density = rho2(i_part)
       
    END IF

    RETURN

  END FUNCTION particles_density

  !******************************************************************************
  !> \brief Brownian aggregation
  !
  !> This function evaluates the aggregation kernel using a Brownian formulation
  !> given in Marchisio et al., 2003.
  !> \param[in]   i_part   particle phase index 
  !> \param[in]   diam_i   first particle diameter (m or phi)
  !> \param[in]   diam_j   second particle diameter (m or phi) 
  !> \date 05/05/2015
  !> @author 
  !> Mattia de' Michieli Vitturi
  !******************************************************************************

  FUNCTION particles_beta(diam_i,diam_j)
    !
    IMPLICIT NONE

    REAL*8 :: particles_beta

    REAL*8, INTENT(IN) :: diam_i
    REAL*8, INTENT(IN) :: diam_j

    ! Diameters in meters
    REAL*8 :: diam_im , diam_jm


    IF ( distribution_variable .EQ. 'mass_fraction' ) THEN

       diam_im = 1.D-3 * 2.D0 ** ( - diam_i )
       diam_jm = 1.D-3 * 2.D0 ** ( - diam_j )

    ELSE

       diam_im = diam_i
       diam_jm = diam_j

    END IF

    particles_beta = ( diam_im + diam_jm ) ** 2 / ( diam_im + diam_jm ) 

    RETURN

  END FUNCTION particles_beta


  !******************************************************************************
  !> \brief Aggregation kernel 
  !
  !> This function evaluates the aggregation kernel, using the expression given 
  !> in Textor et al. 2006.
  !> \param[in]   i_part   particle phase index 
  !> \param[in]   diam1    particle diameter (m)
  !> \param[in]   diam2    particle diameter (m)
  !> \date 24/01/2014
  !> @author 
  !> Mattia de' Michieli Vitturi
  !******************************************************************************

  FUNCTION aggregation_kernel(i_part,diam1,diam2)

    IMPLICIT NONE
    
    REAL*8 :: aggregation_kernel

    INTEGER, INTENT(IN) :: i_part
    REAL*8, INTENT(IN) :: diam1
    REAL*8, INTENT(IN) :: diam2

    REAL*8 :: beta
    REAL*8 :: alfa

    beta = collision_kernel(i_part,diam1,diam2)
    
    alfa = coalescence_efficiency(i_part,diam1,diam2)

    aggregation_kernel = beta * alfa

  END FUNCTION aggregation_kernel

  !******************************************************************************
  !> \brief Collision kernel 
  !
  !> \param[in]   i_part   particle phase index 
  !> \param[in]   diam1    particle diameter (m)
  !> \param[in]   diam2    particle diameter (m)
  !> \date 24/01/2014
  !> @author 
  !> Mattia de' Michieli Vitturi
  !******************************************************************************

  FUNCTION collision_kernel(i_part,diam1,diam2)

    USE meteo_module, ONLY : visc_atm

    USE variables, ONLY: pi_g

    IMPLICIT NONE
    
    REAL*8 :: collision_kernel

    INTEGER, INTENT(IN) :: i_part
    REAL*8, INTENT(IN) :: diam1
    REAL*8, INTENT(IN) :: diam2

    !> Brownian motion collisions kernel
    REAL*8 :: beta_B   

    !> Laminar and turbulent fluid shear collisions kernel
    REAL*8 :: beta_S

    !> Differential sedimentation kernel
    REAL*8 :: beta_DS

    !> Boltzmann constant
    REAL*8 :: k_b

    !> Partciles settling velocities
    REAL*8 :: Vs_1 , Vs_2

    !> Gravitational collision efficiency
    REAL*8 :: E_coll

    !> Rate of dissipation of turbulent kinetic energy
    REAL*8 :: epsilon

    !> Fluid shear
    REAL*8 :: Gamma_s

    !> Air kinematic viscosity
    REAL*8 :: air_kin_viscosity

    REAL*8 :: tp

    !!! WARNING: uninitialized variable
    air_kin_viscosity = 1.5D-5
    tp = 1000.D0
    epsilon = 1.D0
    !!!

    k_b =1.3806488D-23 

    beta_B = 2.D0 / 3.D0 * k_b * tp / visc_atm * ( diam1 + diam2 )**2 / ( diam1*diam2 ) 

    Gamma_s = DSQRT( 1.3D0 * epsilon * air_kin_viscosity )

    E_coll = collision_efficiency(i_part,diam1,diam2)

    beta_S = 1.D0 / 6.D0 * Gamma_s * ( diam1 + diam2 )**3

    Vs_1 = particles_settling_velocity(i_part,diam1)

    Vs_2 = particles_settling_velocity(i_part,diam2)

    beta_DS = pi_g / 4.D0 * ( diam1 + diam2 )**2 * ABS( Vs_2 - Vs_1 )

    collision_kernel = beta_B + beta_S + beta_DS

  END FUNCTION collision_kernel

  !******************************************************************************
  !> \brief Collision efficiency 
  !
  !> \param[in]   i_part   particle phase index 
  !> \param[in]   diam1    particle diameter (m)
  !> \param[in]   diam2    particle diameter (m)
  !> \date 24/01/2014
  !> @author 
  !> Mattia de' Michieli Vitturi
  !******************************************************************************

  FUNCTION collision_efficiency(i_part,diam1,diam2)

    USE variables, ONLY : gi

    IMPLICIT NONE
    
    REAL*8 :: collision_efficiency

    INTEGER, INTENT(IN) :: i_part
    REAL*8, INTENT(IN) :: diam1
    REAL*8, INTENT(IN) :: diam2

    REAL*8 :: E_V , E_A

    REAL*8 :: Re

    REAL*8 :: Stokes

    REAL*8 :: kin_visc_air

    !> Partciles settling velocities
    REAL*8 :: Vs_1 , Vs_2


    !!! WARNING: uninitialized variable
    kin_visc_air = 1.5D-5
    E_A = 0.D0
    !!!

    Vs_1 = particles_settling_velocity(i_part,diam1)

    Vs_2 = particles_settling_velocity(i_part,diam2)

    IF ( diam1 .GT. diam2 ) THEN

       Re = diam1 * Vs_1 / kin_visc_air
       
       Stokes = 2.D0 * Vs_2 * ABS( Vs_1 - Vs_2 ) / diam1 * gi

    ELSE

       Re = diam2 * Vs_2 / kin_visc_air 
       
       Stokes = 2.D0 * Vs_1 * ABS( Vs_2 - Vs_1 ) / diam2 * gi

    END IF



    IF ( Stokes > 1.214 ) THEN

       E_V = ( 1.D0 + ( 0.75 * LOG( 2.D0 * Stokes ) / ( Stokes - 1.214 ) ) )** &
            ( -2.D0 )

    ELSE

       E_V = 0.D0

    END IF

    collision_efficiency = ( 60.D0 * E_V + E_A * Re ) / ( 60.D0 * Re )

  END FUNCTION collision_efficiency


  !******************************************************************************
  !> \brief Coalescence efficiency 
  !
  !> \param[in]   i_part   particle phase index 
  !> \param[in]   diam1    particle diameter (m)
  !> \param[in]   diam2    particle diameter (m)
  !> \date 24/01/2014
  !> @author 
  !> Mattia de' Michieli Vitturi
  !******************************************************************************

  FUNCTION coalescence_efficiency(i_part,diam1,diam2)

    USE variables, ONLY: gi

    IMPLICIT NONE
    
    REAL*8 :: coalescence_efficiency

    INTEGER :: i_part
    REAL*8, INTENT(IN) :: diam1
    REAL*8, INTENT(IN) :: diam2

    !> particle Stokes number
    REAL*8 :: Stokes

    !> Critical Stokes number
    REAL*8 :: Stokes_cr

    !> Efficiency exponent
    REAL*8 :: q

    !> Partciles settling velocities
    REAL*8 :: Vs_1 , Vs_2

    REAL*8 :: tp

    !!! UNINITIALIZED VARIALBES: CHECK
    tp = 1000
    !!!

    IF ( tp .LE. 273 ) THEN

       coalescence_efficiency = 0.09D0
       
    ELSE

       Vs_1 = particles_settling_velocity(i_part,diam1)

       Vs_2 = particles_settling_velocity(i_part,diam2)

       IF ( diam1 .GT. diam2 ) THEN

          Stokes = 2.D0 * Vs_2 * ABS( Vs_1 - Vs_2 ) / diam1 * gi
          
       ELSE
          
          Stokes = 2.D0 * Vs_1 * ABS( Vs_2 - Vs_1 ) / diam2 * gi
          
       END IF

       Stokes_cr = 1.3D0
     
       q = 0.8D0

       coalescence_efficiency = 1.D0 / ( 1.D0 + ( Stokes / Stokes_cr ) ) ** q 

    END IF

  END FUNCTION coalescence_efficiency



  !******************************************************************************
  !> \brief Particles moments computation
  !
  !> This subroutine compute the moments of the particles properties (density,
  !> heat capacity and settling velocity) using the quadrature formulas.
  !> \param[in]   xi     abscissas for the quadrature
  !> \param[out]  w      weights for the quadrature
  !> \date 22/10/2013
  !> @author 
  !> Mattia de' Michieli Vitturi
  !******************************************************************************

  SUBROUTINE eval_particles_moments( xi , w )

    ! external variables
    USE variables, ONLY : verbose_level

    ! external procedures
    USE meteo_module, ONLY : zmet

    IMPLICIT NONE

    REAL*8, DIMENSION(n_part,n_nodes), INTENT(IN) :: xi
    REAL*8, DIMENSION(n_part,n_nodes), INTENT(IN) :: w

    REAL*8, DIMENSION(n_part,n_nodes) :: part_dens_array
    REAL*8, DIMENSION(n_part,n_nodes) :: part_set_vel_array
    REAL*8, DIMENSION(n_part,n_nodes) :: part_cp_array
    REAL*8, DIMENSION(n_part,n_nodes,n_nodes) :: part_beta_array

    INTEGER :: i , j , j1 , j2
    INTEGER :: i_part

    CALL zmet

    DO i_part=1,n_part

       DO j=1,n_nodes
          
          part_dens_array(i_part,j) = particles_density( i_part , xi(i_part,j) )
          
          part_set_vel_array(i_part,j) = particles_settling_velocity( i_part ,  &
               xi(i_part,j) ) 

          part_cp_array(i_part,j) = particles_heat_capacity( i_part,xi(i_part,j))  

          IF ( aggregation) THEN
         
             DO j2=1,n_nodes

                part_beta_array(i_part,j,j2) = particles_beta( xi(i_part,j) ,   &
                     xi(i_part,j2) )

             END DO

          END IF
 
       END DO

       IF ( verbose_level .GE. 2 ) THEN

          WRITE(*,*) 'i_part',i_part
          WRITE(*,*) 'abscissas', xi(i_part,1:n_nodes)
          WRITE(*,*) 'weights', w(i_part,1:n_nodes)
          WRITE(*,*) 'part_dens_array',part_dens_array(i_part,:)
          WRITE(*,*) 'part_set_vel_array',part_set_vel_array(i_part,:)
          WRITE(*,*) 'part_cp_array',part_cp_array(i_part,:)

       END IF

    END DO


    DO i_part=1,n_part
       
       DO i=0,n_mom-1
          
          set_mom(i_part,i) = SUM( part_set_vel_array(i_part,:) * w(i_part,:)   &
               * xi(i_part,:)**i ) / mom(i_part,i)
          
          rhop_mom(i_part,i) = SUM( part_dens_array(i_part,:) * w(i_part,:)     &
               * xi(i_part,:)**i ) / mom(i_part,i)
        
          cp_rhop_mom(i_part,i) = SUM( part_cp_array(i_part,:)                  &
               * part_dens_array(i_part,:) * w(i_part,:) * xi(i_part,:)**i )    &
               / mom(i_part,i) 

          set_rhop_mom(i_part,i) = SUM( part_set_vel_array(i_part,:)            &
               *  part_dens_array(i_part,:) * w(i_part,:) * xi(i_part,:)**i )   &
               / mom(i_part,i) 
          
          set_cp_rhop_mom(i_part,i) = SUM( part_set_vel_array(i_part,:)         &
               * part_cp_array(i_part,:) * part_dens_array(i_part,:) *          &
               w(i_part,:) * xi(i_part,:)**i ) / mom(i_part,i) 

          set_cp_mom(i_part,i) = SUM( part_set_vel_array(i_part,:)              &
               * part_cp_array(i_part,:) * w(i_part,:) * xi(i_part,:)**i )      &
               / mom(i_part,i) 

          IF ( aggregation ) THEN

             birth_term(i_part,i) = 0.D0
             death_term(i_part,i) = 0.D0
             
             DO j1=1,n_nodes
                
                DO j2=1,n_nodes
                   
                   birth_term(i_part,i) = birth_term(i_part,i) + w(i_part,j1)   &
                        * w(i_part,j2) * part_beta_array(i_part,j1,j2)          &
                        *( xi(i_part,j1)**3 + xi(i_part,j2)**3 ) ** ( i / 3.D0 )
                   
                   death_term(i_part,i) = death_term(i_part,i) - w(i_part,j1)   &
                        * xi(i_part,j1) * part_beta_array(i_part,j1,j2)         &
                        * w(i_part,j1) * w(i_part,j2) 
                   
                END DO
                
             END DO

             birth_term(i_part,i) = 0.5D0 * birth_term(i_part,i)

          END IF

       END DO
              
    END DO

    RETURN

  END SUBROUTINE eval_particles_moments

END MODULE particles_module

