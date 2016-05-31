!********************************************************************************
!> \brief Gas/particles mixture module
!
!> This module contains all the variables and the procedures related to the 
!> gas-particles mixture.
!> \date 28/10/2013
!> @author 
!> Mattia de' Michieli Vitturi
!********************************************************************************

MODULE mixture_module

  USE variables, ONLY: gi , pi_g

  IMPLICIT NONE

  !> gas mass fraction in the mixture
  REAL*8 :: gas_mass_fraction
  
  !> gas phase density
  REAL*8 :: rho_gas   
  
  !> universal constant for the mixture
  REAL*8 :: rgasmix  
  
  !> specific heat of the mixture
  REAL*8 :: cpmix     
  
  !> mixture density
  REAL*8 :: rho_mix  

  !> logical defining if the plume has neutral density at the base
  LOGICAL :: initial_neutral_density

  !> mixture temperature
  REAL*8 :: tp       
  
  !> gas (initial gas+entrained air) volume fraction in the mixture
  REAL*8 :: gas_volume_fraction
  
  !> solid volume fraction in the mixture
  REAL*8 :: solid_tot_volume_fraction
  
  !> initial temperature 
  REAL*8 :: tp0      
  
  !> initial gas volume fraction
  REAL*8 :: gas_volume_fraction0
  
  !> initial gas mass fraction
  REAL*8 :: gas_mass_fraction0

  !> initial solid mass fraction in the mixture
  REAL*8 :: solid_tot_mass_fraction0
  
  ! mass flow rate
  REAL*8 :: mass_flow_rate

  !> perfect gas constant for water vapour ( J/(kg K) )
  REAL*8 :: rwvapour

  !> specific heat capacity for water vapour
  REAL*8 :: cpwvapour

  !> mass fraction of the entrained air in the mixture
  REAL*8 :: atm_mass_fraction

  !> mass fraction of the volcanic gas in the mixture
  REAL*8 :: wvapour_mass_fraction

  SAVE

CONTAINS

  !******************************************************************************
  !> \brief Mixture properties initialization
  !
  !> This subroutine initialize the properties of the gas-particles mixture.
  !
  !> \date 22/10/2013
  !> @author 
  !> Mattia de' Michieli Vitturi
  !******************************************************************************

  SUBROUTINE initialize_mixture

    ! external variables
    USE meteo_module, ONLY : ta , pa , rho_atm , rair

    USE moments_module, ONLY : n_nodes

    USE particles_module, ONLY: n_part , solid_partial_mass_fraction ,          &
         cp_rhop_mom , mom , rhop_mom , distribution

    USE particles_module, ONLY: distribution_variable

    USE plume_module, ONLY: w , r , u , mag_u , phi , mfr_exp0, r0

    USE variables, ONLY: verbose_level

    ! external procedures
    USE moments_module, ONLY: wheeler_algorithm
    USE particles_module, ONLY: eval_particles_moments 
    USE particles_module, ONLY: particles_density

    IMPLICIT NONE

    REAL*8 :: rho_solid_avg(n_part)

    REAL*8 :: rho_solid_tot_avg

    REAL*8 :: rho_wvapour

    REAL*8 :: rho_atm_tp

    REAL*8 :: alfa_s(n_part)

    REAL*8 :: alfa_g_atm , alfa_g_wvapour

    REAL*8 :: wvapour_volume_fraction , atm_volume_fraction 

    REAL*8 :: cp_solid0

    REAL*8, DIMENSION(n_part,n_nodes) :: xi , wi

    INTEGER :: i_part

    INTEGER :: i

    REAL*8 :: part_dens_array(n_nodes)

    IF ( verbose_level .GE. 1 ) WRITE(*,*) 'initialize_mixture'

    tp = tp0

    IF ( initial_neutral_density ) THEN

       rgasmix = rair

    ELSE

       rgasmix = rwvapour

    END IF

    rho_atm_tp = pa / ( rair * tp )

    rho_gas = pa / ( rgasmix * tp )

    rho_wvapour = pa / ( rwvapour * tp )

    alfa_g_wvapour = ( rho_gas - rho_atm_tp ) / ( rho_wvapour - rho_atm_tp )

    alfa_g_atm = 1.D0 - alfa_g_wvapour

    DO i_part=1,n_part
       
       IF ( distribution .EQ. 'constant' ) THEN

          CALL wheeler_algorithm( mom(i_part,0:1) , xi(i_part,:) ,              &
               wi(i_part,:) )

       ELSE

          CALL wheeler_algorithm( mom(i_part,:) , xi(i_part,:) , wi(i_part,:) )
          
       END IF

       IF ( verbose_level .GE. 2 ) THEN

          WRITE(*,*) 'i_part',i_part
          WRITE(*,*) 'xi',xi(i_part,:)
          WRITE(*,*) 'wi',wi(i_part,:)
       
       END IF

       DO i=1,n_nodes

          part_dens_array(i) = particles_density( i_part , xi(i_part,i) )

       END DO


       IF ( distribution_variable .EQ. 'particles_number' ) THEN
          
          rho_solid_avg(i_part) = SUM( part_dens_array * wi(i_part,:)           &
               * xi(i_part,:)**3 ) / mom(i_part,3)

       ELSEIF ( distribution_variable .EQ. 'mass_fraction' ) THEN
          
          rho_solid_avg(i_part) = 1.D0 / ( SUM( wi(i_part,:) / part_dens_array )&
               / mom(i_part,0) )
          
       END IF

    END DO

    rho_solid_tot_avg = 1.D0 / SUM( solid_partial_mass_fraction(1:n_part) /     &
         rho_solid_avg(1:n_part) )


    DO i_part = 1,n_part

       alfa_s(i_part) = solid_partial_mass_fraction(i_part) *                   &
            rho_solid_tot_avg / rho_solid_avg(i_part)

       IF ( verbose_level .GE. 1 ) THEN

          WRITE(*,*) 'i_part',i_part
          WRITE(*,*) 'rho_solid_avg',rho_solid_avg(i_part)
          WRITE(*,*) 'rho_solid_avg',rho_solid_avg(i_part)
          WRITE(*,*) 'alfa_s',i_part,alfa_s(i_part)
       
       END IF



    END DO


    rho_solid_tot_avg = 1.D0 / SUM( solid_partial_mass_fraction(1:n_part) /     &
         rho_solid_avg(1:n_part) )

    gas_volume_fraction = rho_solid_tot_avg / ( rho_gas * ( 1.D0 /              &
         gas_mass_fraction0 - 1.D0 ) + rho_solid_tot_avg )

    solid_tot_volume_fraction = 1.D0 - gas_volume_fraction0

    solid_tot_mass_fraction0 = 1.D0 - gas_mass_fraction0

    IF ( verbose_level .GE. 1 ) THEN

       WRITE(*,*) 'solid volume fractions',solid_tot_volume_fraction*alfa_s
       WRITE(*,*) 'solid_tot_volume_fraction',solid_tot_volume_fraction
       WRITE(*,*) 'gas_volume_fraction',gas_volume_fraction
       WRITE(*,*) 'solid_tot_mass_fraction0',solid_tot_mass_fraction0
       READ(*,*)

    END IF

    CALL eval_particles_moments( xi , wi ) 
       
    rho_solid_tot_avg = SUM( alfa_s * rho_solid_avg )

    ! ... Plume density
    rho_mix = gas_volume_fraction * rho_gas + solid_tot_volume_fraction *       &
         rho_solid_tot_avg

    atm_volume_fraction = gas_volume_fraction * alfa_g_atm

    wvapour_volume_fraction = gas_volume_fraction * alfa_g_wvapour

    atm_mass_fraction = atm_volume_fraction * rho_atm / rho_mix

    ! WRITE(*,*) 'gas_volume_fraction * alfa_g_atm',gas_volume_fraction , alfa_g_atm
    ! WRITE(*,*) 'atm_volume_fraction * rho_atm / rho_mix',atm_volume_fraction , rho_atm , rho_mix
    ! WRITE(*,*) 'MIXTURE: atm_mass_fraction',atm_mass_fraction
    ! READ(*,*)


    wvapour_mass_fraction = wvapour_volume_fraction * rho_wvapour / rho_mix


    ! ... At the beginning the plume is a mixture of water-vapour and particles
    ! ... only (without entrained air)

    WRITE(*,*) 'cp_rhop_mom',cp_rhop_mom(:,3)
    WRITE(*,*) 'rhop_mom',rhop_mom(:,3)

    cp_solid0 = SUM(solid_partial_mass_fraction * cp_rhop_mom(:,3) /            &
         rhop_mom(:,3) )

    cpmix = gas_mass_fraction0 * cpwvapour + ( 1.D0 - gas_mass_fraction0 ) *    &
         cp_solid0

    IF ( mfr_exp0 .GT. 0.d0 ) THEN

       mass_flow_rate = 10.0**mfr_exp0
  
       WRITE(*,*) 'WARNING: Fixed mfr =',mass_flow_rate

       IF ( r0 .EQ. 0.D0 ) THEN

          IF ( w .EQ. 0.D0 ) THEN
     
             ! Equation 4 from Carazzo et al. 2008
             w = 138 * DSQRT( gas_mass_fraction0 * 100.d0 )
             mag_u = DSQRT(u*u+w*w)
             phi = ATAN(w/u)

             WRITE(*,*) 'Calculated initial velocity =',w

          END IF
                         
          r = DSQRT( mass_flow_rate / ( pi_g * rho_mix * mag_u ) )
          r0=r

          WRITE(*,*) 'Calculated radius =',r

       END IF

    ELSE

       mass_flow_rate = pi_g * rho_mix * mag_u * (r**2)
       
    END IF

    IF ( verbose_level .GE. 1 ) THEN

       WRITE(*,*) 'cpsolid',cp_solid0
       WRITE(*,*) 'rho_atm',rho_atm
       WRITE(*,*) 'rho_gas',rho_gas
       WRITE(*,*) 'rho_mix',rho_mix
       WRITE(*,*) 'mass_flow_rate',mass_flow_rate
       
       READ(*,*)

    END IF

    gas_mass_fraction = gas_volume_fraction * rho_gas / rho_mix

    WRITE(*,*) 'initial mfr',mass_flow_rate

    RETURN
  END SUBROUTINE initialize_mixture

END MODULE mixture_module

