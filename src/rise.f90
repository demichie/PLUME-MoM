!********************************************************************************
!> \brief Predictor-corrector module
!
!> This module contains the main subroutine of the code, i.e. the solver for the
!> predictor-corrector integration scheme. 
!> \date 23/12/2013
!> @author 
!> Mattia de' Michieli Vitturi
!********************************************************************************
MODULE rise

  IMPLICIT NONE
  SAVE

CONTAINS

  !******************************************************************************
  !> \brief Main subroutine for the integration
  !
  !> This is the main subroutine where the solution is advanced integrating with 
  !> a predictor-corrector scheme. 
  !>
  !> \date 23/12/2013
  !> @author 
  !> Mattia de' Michieli Vitturi
  !******************************************************************************

  SUBROUTINE plumerise

    ! external variables
    USE meteo_module, ONLY : rho_atm
    USE mixture_module, ONLY : gas_mass_fraction , rho_mix, mass_flow_rate
    USE particles_module, ONLY : n_part , mom0 , mom
    USE particles_module, ONLY : distribution_variable
    USE particles_module, ONLY : solid_partial_mass_fraction
    USE plume_module, ONLY: s , w , z , vent_height , r , mag_u
    USE solver_module, ONLY: ds, ds0, f, ftemp, rhs, rhstemp
    USE solver_module, ONLY: f_stepold
    USE variables, ONLY : verbose_level
    USE variables, ONLY : dakota_flag , hysplit_flag
    USE variables, ONLY : pi_g
    USE variables, ONLY : height_weight , height_obj , mu_weight ,    &
       mu_obj , sigma_weight , sigma_obj , skew_weight , skew_obj

    ! external procedures
    USE inpout, ONLY: write_column , write_dakota , write_hysplit
    USE meteo_module, ONLY: zmet, initialize_meteo
    USE mixture_module, ONLY: initialize_mixture
    USE particles_module, ONLY: initialize_particles
    USE plume_module, ONLY: initialize_plume
    USE solver_module, ONLY: rate, lump, marching, unlump


    IMPLICIT NONE

    CHARACTER(LEN=20) :: description

    CHARACTER(len=8) :: x1 ! format descriptor

    INTEGER :: i_part

    REAL*8 :: mu(4)

    REAL*8 :: k1 , k2

    REAL*8 :: mu_phi , sigma_phi , skew_phi

    REAL*8 :: mass_fract

    REAL*8 :: solid_mass_flux , solid_mass_flux0

    REAL*8 :: solid_mass_flux_change

    REAL*8 :: obj_function

    REAL*8 :: w_old , w_oldold
    REAL*8 :: w_minrel , w_maxrel
    REAL*8 :: w_maxabs

    REAL*8 :: check_sb
    REAL*8 :: eps_sb

    REAL*8 :: z_array(10000000)
    REAL*8 :: w_array(10000000)

    REAL*8, ALLOCATABLE :: z_norm(:) , w_norm(:)
    REAL*8, ALLOCATABLE :: first_der_right(:) , first_der_left(:)
    REAL*8, ALLOCATABLE :: sec_der(:) , k(:)
    REAL*8, ALLOCATABLE :: first_der_central(:)

    REAL*8 :: k_max
    
    REAL*8 :: column_regime

    INTEGER :: idx , max_idx

    REAL*8 :: delta_rho

    REAL*8 :: plume_height
    REAL*8 :: height_nbl       
    REAL*8 :: deltarho_min
    REAL*8 :: rho_nbl

    REAL*8 :: deltarho , deltarho_old

    !
    ! ... Set initial conditions at the release height
    !
    CALL initialize_plume

    CALL initialize_meteo

    CALL initialize_particles

    !
    ! ... Get meteo variables at release height
    !
    CALL zmet

    CALL initialize_mixture
    
    description = 'Initial MFR'
    
    CALL WRITE_DAKOTA(description,mass_flow_rate)

    w_old = w
    w_oldold = w

    w_maxabs = w

    w_minrel = w
    w_maxrel = w

    idx = 1

    z_array(idx) = z
    w_array(idx) = w

    delta_rho = rho_mix - rho_atm

    DO i_part=1,n_part

       WRITE(x1,'(I2.2)') i_part ! converting integer to string using a 'internal file'

       description = 'Mean Diameter '//trim(x1)

       CALL write_dakota(description,mom0(i_part,1)/mom0(i_part,0))

       description = 'Sau. Mean Diam. '//trim(x1)

       CALL write_dakota(description,mom0(i_part,3)/mom0(i_part,2))

       !------ Convert from raw moments to central moments ----------------------
       !------ http://mathworld.wolfram.com/CentralMoment.html ------------------

       mu(1) = mom0(i_part,1)/mom0(i_part,0)

       mu(2) = -(mom0(i_part,1)/mom0(i_part,0))**2 + (mom0(i_part,2) /          &
            mom0(i_part,0))

       mu(3) = 2.D0 * (mom0(i_part,1)/mom0(i_part,0))**3                        &
            - 3.D0 * (mom0(i_part,1)/mom0(i_part,0)) * (mom0(i_part,2) /        &
            mom0(i_part,0))    &
            + (mom0(i_part,3)/mom0(i_part,0))

       mu(4) = -3.D0 * (mom0(i_part,1)/mom0(i_part,0))**4                       &
            + 6.D0 * (mom0(i_part,1)/mom0(i_part,0))**2 * (mom0(i_part,2) /     &
            mom0(i_part,0)) - 4.D0 * (mom0(i_part,1)/mom0(i_part,0)) *          &
            (mom0(i_part,3)/mom0(i_part,0)) + (mom0(i_part,4)/mom0(i_part,0)) 

       ! WRITE(*,*) 'mu',mu(1:4)
       ! READ(*,*)

       !       description = 'Variance '//trim(x1)
       !
       !       CALL write_dakota(description,mu(2))
       ! 
       !       description = 'Skewness '//trim(x1)
       !
       !       CALL write_dakota( description , mu(3) / mu(2)**(3.D0/2.D0) )
       ! 
       !       description = 'Kurtosis '//trim(x1)
       !
       !       CALL write_dakota( description , mu(4)/ ( mu(2)**2 ) - 3.D0 )

    END DO

    DO i_part=1,n_part

       WRITE(x1,'(I2.2)') i_part ! converting integer to string using a 'internal file'

       IF ( distribution_variable .EQ. 'particles_number' ) THEN

          k1 = log( 1.D3 * mom0(i_part,1) / mom0(i_part,0) )
          k2 = log( 1.D3 * mom0(i_part,3) / mom0(i_part,2) )

          mu_phi = - 0.25D0 * ( 5*k2 - k1 ) / log(2.D0)

          sigma_phi = sqrt( 0.5d0 * (k2-k1) ) / log(2.D0)

       ELSEIF ( distribution_variable .EQ. 'mass_fraction' ) THEN

          mu_phi = mom0(i_part,1)/mom0(i_part,0)

          sigma_phi = SQRT( -(mom0(i_part,1)/mom0(i_part,0))**2 +               &
               (mom0(i_part,2)/mom0(i_part,0)) )

          skew_phi = ( 2.D0 * (mom0(i_part,1)/mom0(i_part,0))**3                &
               - 3.D0 * mu_phi * ( mom0(i_part,2) / mom0(i_part,0) )            &
               + mom0(i_part,3) / mom0(i_part,0) ) / sigma_phi**3

       END IF

       description = 'Init Avg Diam '//trim(x1)

       CALL write_dakota(description,mu_phi)

       description = 'Init Var Diam '//trim(x1)

       CALL write_dakota(description,sigma_phi)

       description = 'Init Skw Diam '//trim(x1)

       CALL write_dakota(description,skew_phi)

       mass_fract = solid_partial_mass_fraction(i_part) * ( 1.D0 -              &
            gas_mass_fraction)

       description = 'Init Mass Fract '//trim(x1)

       CALL write_dakota(description,mass_fract)

       solid_mass_flux0 = solid_partial_mass_fraction(i_part) * ( 1.D0 -        &
            gas_mass_fraction) * rho_mix * pi_g * r**2 * mag_u

       description = 'Init Solid Flux '//trim(x1)

       CALL write_dakota(description,solid_mass_flux0)

    END DO

    !
    ! ... Lump physical variables 
    !
    CALL lump(f)
    !
    ! ----------------------------------------------------
    ! ... assign initial stepping length and
    ! ... start plume rise marching loop
    ! ----------------------------------------------------
    !
    ds = ds0

    IF ( .NOT.dakota_flag ) CALL write_column

    IF ( hysplit_flag ) CALL write_hysplit(.FALSE.)

    deltarho_min = 1000.D0

    deltarho_old = 0.D0

    main_loop: DO

       f_stepold = f

       !
       ! ---------------------------------------------------------
       ! ... compute meteo conditions and the rhs of the equations 
       ! ... for this location
       ! ---------------------------------------------------------
       !
       IF ( verbose_level .GE. 2 ) THEN

          WRITE(*,*)
          WRITE(*,*) '**************** BEFORE PREDICTOR STEP *****************'
          WRITE(*,*)

       END IF

       CALL unlump(f)

       w_oldold = w_old
       w_old = w

       CALL rate(rhs)
       !
       ! ... predictor step (compute temporary quantities)
       CALL marching(f,ftemp,rhs) 

       IF ( verbose_level .GE. 2 ) THEN

          WRITE(*,*)
          WRITE(*,*) '**************** BEFORE CORRECTOR STEP *****************'
          WRITE(*,*)

       END IF

       CALL unlump(ftemp)
       CALL rate(rhstemp)
       !
       ! ... corrector step
       rhs = 0.5D0 * ( rhstemp - rhs )
       CALL marching(ftemp,f,rhs)

       IF ( verbose_level .GE. 2 ) THEN

          WRITE(*,*)
          WRITE(*,*) '**************** AFTER CORRECTOR STEP *****************'
          WRITE(*,*)

       END IF


       CALL unlump(f)

       !
       ! Reduce step-size condition
       !
       IF ( w .LE. 0.D0) THEN

          ds = 0.5D0 * ds
          f = f_stepold

       ELSE

          idx = idx + 1

          z_array(idx) = z
          w_array(idx) = w

          IF ( ( w_old .LT. w ) .AND. ( w_old .LT. w_oldold ) )  THEN

             w_minrel = w_old

             !WRITE(*,*) 'minrel',w_oldold,w_old,w

          END IF
          
          IF ( w .GT. w_maxabs ) w_maxabs = w

          IF ( ( w_old .GT. w ) .AND. ( w_old .GT. w_oldold ) )  THEN

             w_maxrel = w_old

             !WRITE(*,*) 'maxrel',w_oldold,w_old,w

          END IF

          delta_rho = MIN( delta_rho , rho_mix - rho_atm )

          ! used to define the neutral buoyancy level 
          deltarho =  rho_mix - rho_atm


          IF ( deltarho * deltarho_old .LT. 0.D0 ) THEN

             rho_nbl = rho_mix
             height_nbl = z - vent_height

          END IF

          s = s + ds

          deltarho_old = deltarho

          IF ( verbose_level .GE. 2 ) THEN

             DO i_part=1,n_part

                WRITE(*,*) '**',mom(i_part,1)/mom(i_part,0)

             END DO

             READ(*,*)

          END IF

          IF ( .NOT.dakota_flag ) CALL write_column

          IF ( hysplit_flag ) CALL write_hysplit(.FALSE.)

          ! ... Exit condition

          IF ( w .LE. 1.D-5 ) THEN

             EXIT main_loop

          END IF

       END IF

    END DO main_loop

    IF ( hysplit_flag ) CALL write_hysplit(.TRUE.)


    max_idx = idx

    ALLOCATE( z_norm(max_idx) , w_norm(idx) )
    ALLOCATE( first_der_right(max_idx-1) , first_der_left(max_idx-1) )
    ALLOCATE( sec_der(max_idx-2) , k(max_idx-2) )
    ALLOCATE( first_der_central(max_idx-2) )

    z_norm = ( z_array(1:max_idx) - MINVAL( z_array(1:max_idx) ) )              &
         / ( MAXVAL( z_array(1:max_idx) ) - MINVAL( z_array(1:max_idx) ) )

    w_norm = w_array(1:max_idx) / MAXVAL( w_array(1:max_idx) )

    first_der_right(1:max_idx-2) = ( z_norm(3:max_idx) - z_norm(2:max_idx-1) )  &
         / ( w_norm(3:max_idx) - w_norm(2:max_idx-1) )

    first_der_left(1:max_idx-2) = ( z_norm(2:max_idx-1) - z_norm(1:max_idx-2) ) &
         / ( w_norm(2:max_idx-1) - w_norm(1:max_idx-2) )

    sec_der(1:max_idx-2) = ( first_der_right(1:max_idx-2) - first_der_left(1:max_idx-2) ) &
         / ( 0.5 * ( w_norm(3:max_idx) - w_norm(1:max_idx-2 ) ) )

    first_der_central = ( z_norm(3:max_idx) - z_norm(1:max_idx-2) )  &
         / ( w_norm(3:max_idx) - w_norm(1:max_idx-2) )

    k(1:max_idx-2) = sec_der(1:max_idx-2) / ( ( 1.0 + first_der_central(1:max_idx-2)**2 )**(1.5d0) )

    k_max = MAXVAL( k(1:max_idx-2) )

    check_sb = ( w_maxrel - w_minrel ) / w_maxabs

    eps_sb = 0.05

!!$    IF ( check_sb .GT. eps_sb ) THEN
!!$
!!$       WRITE(*,*) 'w_minrel,w_maxrel,w_maxabs',w_minrel,w_maxrel,w_maxabs
!!$
!!$       WRITE(*,*) 'superbuoyant'
!!$
!!$       column_regime = 2
!!$
!!$    ELSE
!!$
!!$       WRITE(*,*) 'k_max',k_max
!!$
!!$       IF ( k_max < 1.d0 ) THEN
!!$
!!$          WRITE(*,*) 'collapsing'
!!$
!!$          column_regime = 3
!!$
!!$       ELSE
!!$
!!$          WRITE(*,*) 'buoyant'
!!$
!!$          column_regime = 1
!!$
!!$       END IF


       ! modified criterium for collapsing

    IF ( delta_rho .GT. 0.d0 ) THEN
       
       WRITE(*,*) 'collapsing'
       
       column_regime = 3
       
    ELSE
       
       IF ( check_sb .GT. eps_sb ) THEN
          
          !WRITE(*,*) 'w_minrel,w_maxrel,w_maxabs',w_minrel,w_maxrel,w_maxabs
          
          WRITE(*,*) 'superbuoyant'
          
          column_regime = 2
          
       ELSE
          
          WRITE(*,*) 'buoyant'
          
          column_regime = 1
          
       END IF
          
    END IF


    description = 'Column regime'
    
    CALL WRITE_DAKOTA(description,column_regime)

    description = 'NBL height (atv)'

    CALL WRITE_DAKOTA(description,height_nbl)
    
    description = 'Plume height (atv)'
    plume_height = z - vent_height

    CALL WRITE_DAKOTA(description,plume_height)

    DO i_part=1,n_part

       WRITE(x1,'(I2.2)') i_part ! convert int to string using an 'internal file'

       IF ( distribution_variable .EQ. 'particles_number' ) THEN

          k1 = log( 1.D3 * mom(i_part,1) / mom(i_part,0) )
          k2 = log( 1.D3 * mom(i_part,3) / mom(i_part,2) )

          mu_phi = - 0.25D0 * ( 5*k2 - k1 ) / log(2.D0)

          sigma_phi = sqrt( 0.5d0 * (k2-k1) ) / log(2.D0)

       ELSEIF ( distribution_variable .EQ. 'mass_fraction' ) THEN

          mu_phi = mom(i_part,1)/mom(i_part,0)

          sigma_phi = SQRT( -(mom(i_part,1)/mom(i_part,0))**2 +                 &
               (mom(i_part,2)/mom(i_part,0)) )

          skew_phi = ( 2.D0 * (mom0(i_part,1)/mom0(i_part,0))**3                &
               - 3.D0 * mu_phi * ( mom0(i_part,2) / mom0(i_part,0) )            &
               + mom0(i_part,3) / mom0(i_part,0) ) / sigma_phi**3

       END IF

       description = 'Final Avg Diam '//trim(x1)

       CALL write_dakota(description,mu_phi)

       description = 'Final Var Diam '//trim(x1)

       CALL write_dakota(description,sigma_phi)

       description = 'Final Skw Diam '//trim(x1)

       CALL write_dakota(description,skew_phi)

       mass_fract = solid_partial_mass_fraction(i_part) * ( 1.D0 -              &
            gas_mass_fraction)

       description = 'Final Mass Fract '//trim(x1)

       CALL write_dakota(description,mass_fract)

       solid_mass_flux = solid_partial_mass_fraction(i_part) * ( 1.D0 -         &
            gas_mass_fraction) * rho_mix * pi_g * r**2 * mag_u

       description = 'Final Mass Flux '//trim(x1)

       CALL write_dakota(description,solid_mass_flux)

       solid_mass_flux_change = 1.D0 - solid_mass_flux / solid_mass_flux0

       description = 'Solid Flux Lost '//trim(x1)

       CALL write_dakota(description,solid_mass_flux_change)

            
    END DO

    description = 'Objective_function'
    
    obj_function = mu_weight * ( ( mu_phi-mu_obj ) /                            &
         MAX( ABS(mu_phi),ABS(mu_obj) ) )**2                                    &
         + sigma_weight * ( ( sigma_phi-sigma_obj ) /                           &
         MAX( ABS(sigma_phi),ABS(sigma_obj) ) )**2                              &
         + skew_weight * ( ( skew_phi-skew_obj ) /                              &
         MAX( ABS(skew_phi),ABS(skew_obj) ) )**2                                &
         + height_weight * ( ( plume_height-height_obj ) /                      &
         MAX( ABS(plume_height),ABS(height_obj) ) )**2      

    CALL WRITE_DAKOTA(description,obj_function)
    
    WRITE(*,*) 'plume_height,obj_function' , plume_height , obj_function
    WRITE(*,*) 'height_nbl',height_nbl


    RETURN

  END SUBROUTINE plumerise

END MODULE rise
!----------------------------------------------------------------------
