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
  
  !> gas vlume fraction in the mixture
  REAL*8 :: gas_volume_fraction

  !> gas phase density
  REAL*8 :: rho_gas   
  
  !> universal constant for the mixture
  REAL*8 :: rgasmix  
  
  !> mixture density
  REAL*8 :: rho_mix  

  !> logical defining if the plume has neutral density at the base
  LOGICAL :: initial_neutral_density

  !> mixture temperature
  REAL*8 :: tp       
  
  !> water volume fraction in the mixture
  REAL*8 :: water_volume_fraction
  
  !> solid volume fraction in the mixture
  REAL*8 :: solid_tot_volume_fraction
  
  !> initial temperature 
  REAL*8 :: tp0      
  
  !> initial water volume fraction
  REAL*8 :: water_volume_fraction0
  
  !> initial water mass fraction
  REAL*8 :: water_mass_fraction0

  !> solid mass fraction in the mixture
  REAL*8 :: solid_tot_mass_fraction
  
  ! mass flow rate
  REAL*8 :: mass_flow_rate

  !> volcanic gas species number
  INTEGER :: n_gas

  !> volcanic gases densities
  REAL*8, ALLOCATABLE, DIMENSION(:) :: rhovolcgas

  !> gas constants for volcanic gases
  REAL*8, ALLOCATABLE, DIMENSION(:) :: rvolcgas

  !> specific heat capacity for volcanic gases
  REAL*8, ALLOCATABLE, DIMENSION(:) :: cpvolcgas

  !> molecular weight of additional volcanic gases
  REAL*8, ALLOCATABLE, DIMENSION(:) :: volcgas_mol_wt 

  !> initial mass fractions of volcanic gases
  REAL*8, ALLOCATABLE, DIMENSION(:) :: volcgas_mass_fraction0

  !> mass fractions of volcanic gases
  REAL*8, ALLOCATABLE, DIMENSION(:) :: volcgas_mass_fraction

  !> volcanic gases mixture density
  REAL*8 :: rhovolcgas_mix

  !> gas constant of volcanic gases mixture ( J/(kg K) )
  REAL*8 :: rvolcgas_mix

  !> specific heat of volcanic gases mixture
  REAL*8 :: cpvolcgas_mix

  !> mass fraction of the entrained air in the mixture
  REAL*8 :: atm_mass_fraction

  !> mass fraction of the volcanic gas in the mixture
  REAL*8 :: volcgas_mix_mass_fraction

  !> mass fraction of dry air in the mixture
  REAL*8 :: dry_air_mass_fraction

  !> mass fraction of water in the mixture
  REAL*8 :: water_mass_fraction

  !> mass fraction of liquid water in the mixture
  REAL*8 :: liquid_water_mass_fraction

  !> mass fraction of water vapor in the mixture
  REAL*8 :: water_vapor_mass_fraction

  REAL*8 :: volcgas_mix_mol_fract

  REAL*8 :: volcgas_mix_mol_wt

  REAL*8 :: mixture_enthalpy

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
    USE meteo_module, ONLY : pa , rho_atm , rair , rwv , c_wv

    USE moments_module, ONLY : n_nodes

    USE particles_module, ONLY: n_part , solid_partial_mass_fraction ,          &
         cp_rhop_mom , mom , rhop_mom , distribution , cp_mom

    USE particles_module, ONLY: distribution_variable , cpsolid

    USE plume_module, ONLY: w , r , u , mag_u , phi , log10_mfr, r0

    USE variables, ONLY: verbose_level , write_flag

    ! external procedures
    USE moments_module, ONLY: wheeler_algorithm
    USE particles_module, ONLY: eval_particles_moments 
    USE particles_module, ONLY: particles_density
    

    IMPLICIT NONE

    REAL*8 :: rho_solid_avg(n_part)

    REAL*8 :: rho_solid_tot_avg

    REAL*8 :: alfa_s(n_part)

    REAL*8 :: atm_volume_fraction 

    REAL*8 :: volcgas_mix_volume_fraction 

    REAL*8 :: cp_solid0

    REAL*8, DIMENSION(n_part,n_nodes) :: xi , wi

    INTEGER :: i_part

    INTEGER :: i

    INTEGER :: i_gas

    REAL*8 :: part_dens_array(n_nodes)

    REAL*8 :: Rrhovolcgas_mix

    REAL*8 :: rhowv

    IF ( verbose_level .GE. 1 ) WRITE(*,*) 'initialize_mixture'

    tp = tp0

    ! Compute density of gas species and mixture of gas species

    Rrhovolcgas_mix = 0.D0
    rvolcgas_mix = 0.D0
    cpvolcgas_mix = 0.D0

    ! WRITE(*,*) 'n_gas other than H2O',n_gas

    volcgas_mass_fraction(1:n_gas) = volcgas_mass_fraction0(1:n_gas)

    IF ( n_gas .GT. 0 ) THEN
       
       DO i_gas = 1,n_gas

          rvolcgas_mix = rvolcgas_mix + volcgas_mass_fraction(i_gas)               &
               * rvolcgas(i_gas)
       
          cpvolcgas_mix = cpvolcgas_mix + volcgas_mass_fraction(i_gas)             &
               * cpvolcgas(i_gas)

          ! WRITE(*,*)   pa , rvolcgas(i_gas) , tp0

          Rrhovolcgas_mix = Rrhovolcgas_mix + volcgas_mass_fraction(i_gas)         &
               / (  pa / ( rvolcgas(i_gas) * tp0 ) )

       END DO
    
       rvolcgas_mix = rvolcgas_mix / SUM( volcgas_mass_fraction(1:n_gas) )
    
       cpvolcgas_mix = cpvolcgas_mix / SUM( volcgas_mass_fraction(1:n_gas) )

       rhovolcgas_mix =  SUM(volcgas_mass_fraction(1:n_gas)) / Rrhovolcgas_mix

       volcgas_mix_mass_fraction = SUM(volcgas_mass_fraction(1:n_gas))

    ELSE

       rvolcgas_mix = 0.D0
       
       cpvolcgas_mix = 0.D0
       
       rhovolcgas_mix =  0.D0
       
       volcgas_mix_mass_fraction = 0.D0
    
    END IF

    ! WRITE(*,*) 'rvolcgas_mix',rvolcgas_mix
    ! WRITE(*,*) 'cpvolcgas_mix',cpvolcgas_mix
    ! WRITE(*,*) 'rhovolcgas_mix',rhovolcgas_mix
    ! WRITE(*,*) 'volcgas_mix_mass_fraction',volcgas_mix_mass_fraction

    
    rhowv = pa / ( rwv * tp0 )

    water_mass_fraction = water_mass_fraction0

    ! ---- We assume all volcanic H2O at the vent is water vapor 
    water_vapor_mass_fraction = water_mass_fraction0
    liquid_water_mass_fraction = 0.D0

    ! ---- No air is entrained
    dry_air_mass_fraction = 0.D0

    gas_mass_fraction = water_vapor_mass_fraction + volcgas_mix_mass_fraction 

    !WRITE(*,*) '-------->gas_mass_fraction',gas_mass_fraction
    !WRITE(*,*) '-------->water_vapor_mass_fraction',water_vapor_mass_fraction
    !WRITE(*,*) '-------->rhowv ',rhowv 
    !WRITE(*,*) '--------> volcgas_mix_mass_fraction', volcgas_mix_mass_fraction
    !WRITE(*,*) '-------->rhovolcgas_mix ',rhovolcgas_mix 


    IF ( n_gas .GT. 0 ) THEN 

       rho_gas = gas_mass_fraction / (  water_vapor_mass_fraction / rhowv          &
            + volcgas_mix_mass_fraction / rhovolcgas_mix ) 
    ELSE

       rho_gas = gas_mass_fraction / (  water_vapor_mass_fraction / rhowv)

    END IF

    !WRITE(*,*) '-------->rho_gas',rho_gas
    !!READ(*,*)

    DO i_part=1,n_part
       
       IF ( distribution .EQ. 'constant' ) THEN

          CALL wheeler_algorithm( mom(i_part,0:1), distribution, xi(i_part,:),  &
               wi(i_part,:) )

       ELSE

          CALL wheeler_algorithm( mom(i_part,:) , distribution , xi(i_part,:) , &
               wi(i_part,:) )
          
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
         gas_mass_fraction - 1.D0 ) + rho_solid_tot_avg )

    solid_tot_volume_fraction = 1.D0 - gas_volume_fraction

    solid_tot_mass_fraction = 1.D0 - gas_mass_fraction



    IF ( verbose_level .GE. 1 ) THEN

       WRITE(*,*) 'solid volume fractions',solid_tot_volume_fraction*alfa_s
       WRITE(*,*) 'solid_tot_volume_fraction',solid_tot_volume_fraction
       WRITE(*,*) 'water_volume_fraction',water_volume_fraction
       WRITE(*,*) 'solid_tot_mass_fraction',solid_tot_mass_fraction
       !READ(*,*)

    END IF

    CALL eval_particles_moments( xi , wi ) 
 


    IF ( distribution_variable .EQ. "particles_number" ) THEN

       cpsolid = ( SUM( solid_partial_mass_fraction(1:n_part) * cp_rhop_mom(1:n_part,3)     &
            / rhop_mom(1:n_part,3) ) ) / ( SUM( solid_partial_mass_fraction(1:n_part) ) ) 

    ELSEIF ( distribution_variable .EQ. "mass_fraction" ) THEN

       cpsolid = ( SUM( solid_partial_mass_fraction(1:n_part) * cp_mom(1:n_part,0) ) )      &
            / ( SUM( solid_partial_mass_fraction(1:n_part) ) ) 

    END IF 

      
    rho_solid_tot_avg = SUM( alfa_s * rho_solid_avg )

    ! ... Plume density
    rho_mix = gas_volume_fraction * rho_gas + solid_tot_volume_fraction *       &
         rho_solid_tot_avg

    volcgas_mix_mass_fraction = SUM( volcgas_mass_fraction(1:n_gas) )
    !WRITE(*,*) 'volcgas_mix_mass_fraction',volcgas_mix_mass_fraction
    !READ(*,*)

    ! ... At the beginning the plume is a mixture of water-vapour , volcanic gas
    ! ... and particles only (without entrained air)

    IF ( verbose_level .GE. 1 ) THEN

       WRITE(*,*) 'cp_rhop_mom',cp_rhop_mom(:,3)
       WRITE(*,*) 'rhop_mom',rhop_mom(:,3)

    END IF

    IF ( distribution_variable .EQ. "particles_number" ) THEN

       cp_solid0 = SUM(solid_partial_mass_fraction * cp_rhop_mom(:,3) /         &
         rhop_mom(:,3) )

    ELSEIF ( distribution_variable .EQ. "mass_fraction" ) THEN

       cp_solid0 = SUM(solid_partial_mass_fraction * cp_mom(:,0) )

    END IF
    
    IF ( log10_mfr .GT. 0.d0 ) THEN

       mass_flow_rate = 10.0** log10_mfr
  
       WRITE(*,*) 'WARNING: Fixed MER =',mass_flow_rate

       IF ( r0 .EQ. 0.D0 ) THEN

          IF ( w .EQ. 0.D0 ) THEN
     
             ! Equation 4 from Carazzo et al. 2008
             w = 138 * DSQRT( water_mass_fraction0 * 100.d0 )
             mag_u = DSQRT(u*u+w*w)
             phi = ATAN(w/u)

             WRITE(*,*) 'WARNING: calculated initial velocity =',w

          END IF
                         
          r = DSQRT( mass_flow_rate / ( pi_g * rho_mix * mag_u ) )
          r0=r

          IF ( write_flag) WRITE(*,*) 'WARNING: Initial radius [m] computed from MER and w0 =',r

       END IF

    ELSE

       mass_flow_rate = pi_g * rho_mix * mag_u * (r**2)
       IF ( write_flag) WRITE(*,'(1x,A,1x,es15.8)') 'Initial MER [kgs-1] computed from r0 and w0 =',mass_flow_rate
       
    END IF

    IF ( verbose_level .GE. 1 ) THEN

       WRITE(*,*) 'cpsolid',cp_solid0
       WRITE(*,*) 'rho_atm',rho_atm
       WRITE(*,*) 'rho_gas',rho_gas
       WRITE(*,*) 'rho_mix',rho_mix
       WRITE(*,*) 'mass_flow_rate',mass_flow_rate
       WRITE(*,*) 'solid_mass_flow_rates',mass_flow_rate *                      &
            ( 1.D0 - gas_mass_fraction ) * solid_partial_mass_fraction(1:n_part)
       
       !READ(*,*)

    END IF

    RETURN
  END SUBROUTINE initialize_mixture


  SUBROUTINE eval_wv(pres,temp,w_mf,wv_mf)
    
    USE meteo_module, ONLY : wv_mol_wt
    
    IMPLICIT none
    
    REAL*8, INTENT(IN) :: pres
    REAL*8, INTENT(IN) :: temp
    REAL*8, INTENT(IN) :: w_mf

    ! Mass fraction of liquid water in the mixture
    REAL*8, INTENT(OUT) :: wv_mf
    
    !> saturation pressure (hPa)
    REAL*8 :: el

    !> water vapor molar fraction
    REAL*8 :: wv_mol_fract
    
    ! WRITE(*,*) 'temp',temp

    el = 6.112D0 * DEXP( 17.67D0 * ( temp - 273.16D0 ) / ( temp - 29.65D0 ) )

    wv_mol_fract = el / ( pres / 100.D0 )

    ! WRITE(*,*)  wv_mol_fract , el , pres

    volcgas_mix_mol_wt = SUM( volcgas_mass_fraction(1:n_gas) ) /                &
         SUM( volcgas_mass_fraction(1:n_gas) / volcgas_mol_wt(1:n_gas ) ) 

    wv_mf = wv_mol_fract / ( 1.D0 - wv_mol_fract ) * wv_mol_wt                  &
         / volcgas_mix_mol_wt * volcgas_mix_mass_fraction

    wv_mf = MAX( 0.D0 , wv_mf )
    wv_mf = MIN( wv_mf , w_mf )

    ! WRITE(*,*) 'water vapor mass fraction',wv_mf 
    ! !READ(*,*)

  END SUBROUTINE eval_wv


  SUBROUTINE eval_temp(enth,pres,cpsolid,temp,wv_mf)



   
    USE meteo_module, ONLY : cpair , T_ref , h_wv0 , c_wv , h_lw0 , c_lw ,       &
         da_mol_wt , wv_mol_wt

    ! USE meteo_module

    IMPLICIT none

    !> mixture enthalpy
    REAL*8, INTENT(IN) :: enth

    !> pressure in Pa
    REAL*8, INTENT(IN) :: pres
    
    REAL*8, INTENT(IN) :: cpsolid

    !> mixture temperature in K
    REAL*8, INTENT(OUT) :: temp

    ! Mass fraction of liquid water in the mixture
    REAL*8, INTENT(OUT) :: wv_mf

    !> water vapor molar fraction
    REAL*8 :: wv_mol_fract

    !> dry air molar fraction
    REAL*8 :: da_mol_fract

    !> saturation pressure (hPa)
    REAL*8 :: el

    !> pressure in hPa
    REAL*8 :: hPres

    REAL*8 :: lw_mf0 , lw_mf1 , lw_mf2

    REAL*8 :: f0,f1,f2
   
    REAL*8 :: temp0 , temp1 , temp2

    hPres = pres / 100.D0



    ! WRITE(*,*) 'water_mass_fraction', water_mass_fraction
    ! WRITE(*,*) 'volcgas_mix_mass_fraction', volcgas_mix_mass_fraction
    ! WRITE(*,*) 'dry_air_mass_fraction', dry_air_mass_fraction
    ! WRITE(*,*) 'water_vapor_mass_fraction', water_vapor_mass_fraction

    ! WRITE(*,*)
    !WRITE(*,*) '************** EVAL TEMP **************' 

    IF ( dry_air_mass_fraction .EQ. 0.D0 ) THEN


       liquid_water_mass_fraction = 0.D0
       water_vapor_mass_fraction = water_mass_fraction - liquid_water_mass_fraction 

       temp = ( enth - liquid_water_mass_fraction * ( h_lw0 - c_lw * T_ref )       &
            - water_vapor_mass_fraction * ( h_wv0 - c_wv * T_ref ) ) /             &
            ( dry_air_mass_fraction * cpair + solid_tot_mass_fraction * cpsolid    &
            + liquid_water_mass_fraction * c_lw + water_vapor_mass_fraction * c_wv   &
            +  volcgas_mix_mass_fraction * cpvolcgas_mix)
       
       wv_mf = water_mass_fraction - liquid_water_mass_fraction

        !WRITE(*,*)  dry_air_mass_fraction , cpair , solid_tot_mass_fraction ,       &
        !    cpsolid , water_vapor_mass_fraction , liquid_water_mass_fraction ,     &
        !     volcgas_mix_mass_fraction , cpvolcgas_mix
        !WRITE(*,*) 

       
        !WRITE(*,*) 'cpvolcgas_mix',cpvolcgas_mix      
        !WRITE(*,*) 'temp',temp
        !WRITE(*,*) 'wv_mf',wv_mf
        !!READ(*,*)

       RETURN

    END IF


    IF ( n_gas .GT. 0) THEN

       volcgas_mix_mol_wt = SUM( volcgas_mass_fraction(1:n_gas) ) /                &
            SUM( volcgas_mass_fraction(1:n_gas) / volcgas_mol_wt(1:n_gas ) ) 

    ELSE

       volcgas_mix_mol_wt=0

    END IF



    ! ---- All water is liquid 

    !WRITE(*,*) '! ---- All water is liquid'

    lw_mf0 = water_mass_fraction 

    liquid_water_mass_fraction = lw_mf0
    water_vapor_mass_fraction = water_mass_fraction - liquid_water_mass_fraction

    

    !volcgas_mix_mol_wt = SUM( volcgas_mass_fraction(1:n_gas) ) /                &
    !     SUM( volcgas_mass_fraction(1:n_gas) / volcgas_mol_wt(1:n_gas ) ) 

     !WRITE(*,*) '--->volcgas_mix_mass_fraction', volcgas_mix_mass_fraction
     !WRITE(*,*) '--->volcgas_mix_mol_wt',volcgas_mix_mol_wt
     !WRITE(*,*) '--->water_vapor_mass_fraction',water_vapor_mass_fraction
     !WRITE(*,*) '--->wv_mol_wt ', wv_mol_wt 
     !WRITE(*,*) '--->volcgas_mix_mol_wt',volcgas_mix_mol_wt    
     !WRITE(*,*) '--->dry_air_mass_fraction',dry_air_mass_fraction
     !!READ(*,*)

    IF ( n_gas .GT. 0) THEN

       wv_mol_fract = ( water_vapor_mass_fraction / wv_mol_wt ) /                  &
            ( water_vapor_mass_fraction / wv_mol_wt                                &
            + volcgas_mix_mass_fraction / volcgas_mix_mol_wt                       &
            + dry_air_mass_fraction / da_mol_wt )

       volcgas_mix_mol_fract = ( volcgas_mix_mass_fraction / volcgas_mix_mol_wt ) / &
            ( water_vapor_mass_fraction / wv_mol_wt                                 &
            + volcgas_mix_mass_fraction / volcgas_mix_mol_wt                        &
            + dry_air_mass_fraction / da_mol_wt )
    
       da_mol_fract = ( dry_air_mass_fraction / da_mol_wt ) /                       &
            ( water_vapor_mass_fraction / wv_mol_wt                                 &
            + volcgas_mix_mass_fraction / volcgas_mix_mol_wt                        &
            + dry_air_mass_fraction / da_mol_wt )
    ELSE

       wv_mol_fract = ( water_vapor_mass_fraction / wv_mol_wt ) /                  &
            ( water_vapor_mass_fraction / wv_mol_wt                                &
            + dry_air_mass_fraction / da_mol_wt )

       volcgas_mix_mol_fract = 0
    
       da_mol_fract = ( dry_air_mass_fraction / da_mol_wt ) /                       &
            ( water_vapor_mass_fraction / wv_mol_wt                                 &
            + dry_air_mass_fraction / da_mol_wt )

    END IF
    

    
  
     !WRITE(*,*)  dry_air_mass_fraction , water_vapor_mass_fraction

     !WRITE(*,*) 'wv_mol_fract',wv_mol_fract
     !WRITE(*,*) 'da_mol_fract',da_mol_fract
     !WRITE(*,*) 'volcgas_mix_mol_fract',volcgas_mix_mol_fract

  
    temp0 = ( enth - liquid_water_mass_fraction * ( h_lw0 - c_lw * T_ref )       &
         - water_vapor_mass_fraction * ( h_wv0 - c_wv * T_ref ) ) /             &
         ( dry_air_mass_fraction * cpair + solid_tot_mass_fraction * cpsolid    &
         + liquid_water_mass_fraction * c_lw + water_vapor_mass_fraction * c_wv   &
         +  volcgas_mix_mass_fraction * cpvolcgas_mix)
    
    !WRITE(*,*) 'temp0',temp0

    IF ( temp0 .GT. 29.65D0 ) THEN

       IF ( temp0 .LT. 273.16D0 ) THEN
          
          el = 6.112D0 * DEXP( 17.67D0 * ( temp0 - 273.16D0 ) / ( temp0 - 29.65D0 ) )

       ELSE
          
          el = 1.D0

       END IF
                    
    ELSE

       el = 0.D0
       
    END IF

    el = 6.112D0 * DEXP( 17.67D0 * ( temp0 - 273.16D0 ) / ( temp0 - 29.65D0 ) )

    !WRITE(*,*) 'el',el

    f0 = ( hPres - el ) * wv_mol_fract - el * da_mol_fract - el * volcgas_mix_mol_fract

    ! ---- All water is vapor

    !WRITE(*,*) '! ---- All water is vapor'

    lw_mf2 = 0.D0

    liquid_water_mass_fraction = lw_mf2
    water_vapor_mass_fraction = water_mass_fraction - liquid_water_mass_fraction 

     !WRITE(*,*) '--->volcgas_mix_mass_fraction', volcgas_mix_mass_fraction
     !WRITE(*,*) '--->volcgas_mix_mol_wt',volcgas_mix_mol_wt
     !WRITE(*,*) '--->water_vapor_mass_fraction',water_vapor_mass_fraction
     !WRITE(*,*) '--->wv_mol_wt ', wv_mol_wt 
     !WRITE(*,*) '--->volcgas_mix_mol_wt',volcgas_mix_mol_wt    
     !WRITE(*,*) '--->dry_air_mass_fraction',dry_air_mass_fraction

    IF ( n_gas .GT. 0) THEN

       wv_mol_fract = ( water_vapor_mass_fraction / wv_mol_wt ) /               &
            ( water_vapor_mass_fraction / wv_mol_wt                             &
            + volcgas_mix_mass_fraction / volcgas_mix_mol_wt                    &
            + dry_air_mass_fraction / da_mol_wt )

       volcgas_mix_mol_fract = ( volcgas_mix_mass_fraction /                    &
            volcgas_mix_mol_wt ) / ( water_vapor_mass_fraction / wv_mol_wt      &
            + volcgas_mix_mass_fraction / volcgas_mix_mol_wt                    &
            + dry_air_mass_fraction / da_mol_wt )
    
       da_mol_fract = ( dry_air_mass_fraction / da_mol_wt ) /                   &
            ( water_vapor_mass_fraction / wv_mol_wt                             &
            + volcgas_mix_mass_fraction / volcgas_mix_mol_wt                    &
            + dry_air_mass_fraction / da_mol_wt )

    ELSE

       wv_mol_fract = ( water_vapor_mass_fraction / wv_mol_wt ) /               &
            ( water_vapor_mass_fraction / wv_mol_wt                             &
            + dry_air_mass_fraction / da_mol_wt )

       volcgas_mix_mol_fract = 0.d0
    
       da_mol_fract = ( dry_air_mass_fraction / da_mol_wt ) /                   &
            ( water_vapor_mass_fraction / wv_mol_wt                             &
            + dry_air_mass_fraction / da_mol_wt )

    END IF

    temp2 = ( enth - liquid_water_mass_fraction * ( h_lw0 - c_lw * T_ref )      &
         - water_vapor_mass_fraction * ( h_wv0 - c_wv * T_ref ) ) /             &
         ( dry_air_mass_fraction * cpair + solid_tot_mass_fraction * cpsolid    &
         + liquid_water_mass_fraction * c_lw + water_vapor_mass_fraction * c_wv &
         +  volcgas_mix_mass_fraction * cpvolcgas_mix)
    
    !WRITE(*,*) 'water_vapor_mass_fraction',water_vapor_mass_fraction
    !WRITE(*,*) 'wv_mol_fract',wv_mol_fract
    !WRITE(*,*) 'da_mol_fract',da_mol_fract
    !WRITE(*,*) 'volcgas_mix_mol_fract',volcgas_mix_mol_fract
    !WRITE(*,*) 'temp2',temp2

    IF ( temp2 .GT. 29.65D0 ) THEN

       IF ( temp2 .LT. 273.16D0 ) THEN
          
          el = 6.112D0 * DEXP( 17.67D0 * ( temp2 - 273.16D0 ) / ( temp2 - 29.65D0 ) )

       ELSE
          
          el = 1.D0

       END IF
                    
    ELSE

       el = 0.D0
       
    END IF

    el = 6.112D0 * DEXP( 17.67D0 * ( temp2 - 273.16D0 ) / ( temp2 - 29.65D0 ) )

       
    f2 = ( hPres - el ) * wv_mol_fract - el * da_mol_fract - el * volcgas_mix_mol_fract

       
    !WRITE(*,*) 'f0,f2,el',f0,f2,el
    
    IF ( ( f0 .LT. 0.D0 ) .AND. ( f2 .LT. 0.D0 ) ) THEN
        
       liquid_water_mass_fraction = 0.D0
       water_vapor_mass_fraction = water_mass_fraction - liquid_water_mass_fraction 

       temp = temp2
       wv_mf = water_vapor_mass_fraction

       RETURN
       
    ELSEIF ( ( f0 .GT. 0.D0 ) .AND. ( f2 .GT. 0.D0 ) ) THEN

       liquid_water_mass_fraction = water_mass_fraction
       water_vapor_mass_fraction = water_mass_fraction - liquid_water_mass_fraction 

       temp = temp0
       wv_mf = water_vapor_mass_fraction

       RETURN

    ELSE

       
       find_temp:DO  

          lw_mf1 = 0.5D0 * ( lw_mf0 + lw_mf2 )
          
          liquid_water_mass_fraction = lw_mf1
          water_vapor_mass_fraction = water_mass_fraction - liquid_water_mass_fraction 

         IF ( n_gas .GT. 0) THEN
          
              wv_mol_fract = ( water_vapor_mass_fraction / wv_mol_wt ) /        &
                   ( water_vapor_mass_fraction / wv_mol_wt                      &
                   + volcgas_mix_mass_fraction / volcgas_mix_mol_wt             &
                   + dry_air_mass_fraction / da_mol_wt )
              
              volcgas_mix_mol_fract = ( volcgas_mix_mass_fraction /             &
                   volcgas_mix_mol_wt ) / ( water_vapor_mass_fraction /         &
                   wv_mol_wt + volcgas_mix_mass_fraction / volcgas_mix_mol_wt   &
                   + dry_air_mass_fraction / da_mol_wt )
              
              da_mol_fract = ( dry_air_mass_fraction / da_mol_wt ) /            &
                   ( water_vapor_mass_fraction / wv_mol_wt                      &
                   + volcgas_mix_mass_fraction / volcgas_mix_mol_wt             &
                   + dry_air_mass_fraction / da_mol_wt )

   
          ELSE

              wv_mol_fract = ( water_vapor_mass_fraction / wv_mol_wt ) /        &
                   ( water_vapor_mass_fraction / wv_mol_wt                      &
                   + dry_air_mass_fraction / da_mol_wt )
              
              volcgas_mix_mol_fract = 0
              
              da_mol_fract = ( dry_air_mass_fraction / da_mol_wt ) /            &
                   ( water_vapor_mass_fraction / wv_mol_wt                      &
                   + dry_air_mass_fraction / da_mol_wt )            
                  
          END IF

          !WRITE(*,*) 'water_vapor_mass_fraction',water_vapor_mass_fraction
          !WRITE(*,*) 'wv_mol_fract',wv_mol_fract
          !WRITE(*,*) 'da_mol_fract',da_mol_fract
          !WRITE(*,*) 'volcgas_mix_mol_fract',volcgas_mix_mol_fract


          temp1 = ( enth - liquid_water_mass_fraction * ( h_lw0 - c_lw*T_ref )  &
               - water_vapor_mass_fraction * ( h_wv0 - c_wv * T_ref ) ) /       &
               ( dry_air_mass_fraction * cpair + solid_tot_mass_fraction *      &
               cpsolid + liquid_water_mass_fraction * c_lw +                    &
               water_vapor_mass_fraction * c_wv +  volcgas_mix_mass_fraction *  &
               cpvolcgas_mix)

          IF ( temp1 .GT. 29.65D0 ) THEN
             
             IF ( temp1 .LT. 273.16D0 ) THEN
                
                el = 6.112D0 * DEXP( 17.67D0 * ( temp1 - 273.16D0 ) / ( temp1 - 29.65D0 ) )
                
             ELSE
                
                el = 1.D0
                
             END IF
             
          ELSE
             
             el = 0.D0
             
          END IF

          el = 6.112D0 * DEXP( 17.67D0 * ( temp1 - 273.16D0 ) / ( temp1 - 29.65D0 ) )

          
          f1 = ( hPres - el ) * wv_mol_fract - el * da_mol_fract - el *         &
               volcgas_mix_mol_fract


          !WRITE(*,*) 'volcgas_mix_mol_fract ',volcgas_mix_mol_fract
          !WRITE(*,*) wv_mol_fract+volcgas_mix_mol_fract+da_mol_fract

          !WRITE(*,*) 't0,t1,t2',temp0,temp1,temp2
          !WRITE(*,*) 'lw_mf0,lw_mf1,lw_mf2',lw_mf0,lw_mf1,lw_mf2
          !WRITE(*,*) 'f0,f1,f2',f0,f1,f2
          !READ(*,*)

          IF (  f1 * f0 .LT. 0.D0 ) THEN

             lw_mf2 = lw_mf1
             f2 = f1
             temp2 = temp1

          ELSE

             lw_mf0 = lw_mf1
             f0 = f1
             temp0 = temp1


          END IF



          IF ( DABS(temp2-temp0) .LT. 1.D-3 ) THEN

             temp = temp1
             
             wv_mf = water_mass_fraction - lw_mf1


             RETURN

          ELSEIF ( DABS(lw_mf2 - lw_mf0) .LT. 1.D-5 ) THEN

             temp = temp1
             
             wv_mf = water_mass_fraction - lw_mf1

             RETURN

          END IF
                                 
       END DO find_temp

    END IF



  END SUBROUTINE eval_temp


END MODULE mixture_module

