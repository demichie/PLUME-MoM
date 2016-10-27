 !********************************************************************************
!> \brief Solver module
!
!> This module contains the procedures to evaluate the terms of the differential 
!> equations of the model.
!> \date 22/10/2013
!> @author 
!> Mattia de' Michieli Vitturi
!********************************************************************************

MODULE solver_module
  USE moments_module, ONLY: n_mom
  USE particles_module, ONLY : n_part , distribution_variable
  !
  IMPLICIT NONE

  !> Right-Hand Side (rhs)
  REAL*8, ALLOCATABLE, DIMENSION(:) :: rhstemp

  !> Right-Hand Side (rhs)
  REAL*8, ALLOCATABLE, DIMENSION(:) :: rhs

  !> Integrated variables
  REAL*8, ALLOCATABLE, DIMENSION(:) :: ftemp

  !> Integrated variables
  REAL*8, ALLOCATABLE, DIMENSION(:) :: f

  !> Integrated variables
  REAL*8, ALLOCATABLE, DIMENSION(:) :: f_stepold

  !> Total number of equations
  INTEGER :: itotal

  !> Integration step
  REAL*8 :: ds                 

  !> Initial integration step
  REAL*8 :: ds0                
  !
  SAVE

CONTAINS

  !******************************************************************************
  !> \brief Solver variables allocation
  !
  !> This subroutine allocate the variables defining the terms on the two sides 
  !> of equations of the model (f and rhs).
  !> \date 22/10/2013
  !> @author 
  !> Mattia de' Michieli Vitturi
  !******************************************************************************
  
  SUBROUTINE allocate_matrix

    IMPLICIT NONE
    !
    itotal = n_part * n_mom + 9
    !
    ALLOCATE(f(itotal))
    ALLOCATE(f_stepold(itotal))
    ALLOCATE(rhs(itotal))
    ALLOCATE(rhstemp(itotal))
    ALLOCATE(ftemp(itotal))

    !
    f = 0.D0
    f_stepold = 0.D0
    ftemp = 0.D0
    rhs = 0.D0
    rhstemp = 0.D0
    !
    RETURN
  END SUBROUTINE allocate_matrix

  !******************************************************************************
  !> \brief Compute the right-hand side of the equations
  !
  !> This subroutine compute the right-hand side of the equations. 
  !> \param[out]    rhs_     right-hand side
  !> \date 22/10/2013
  !> @author 
  !> Mattia de' Michieli Vitturi
  !******************************************************************************
  
  SUBROUTINE rate(rhs_)

    USE meteo_module, ONLY: u_atm , rho_atm , ta , duatm_dz , cpair , rair ,    &
         cos_theta , sin_theta

    USE mixture_module, ONLY: rho_mix , tp , cpmix , rgasmix , rho_gas ,        &
         gas_mass_fraction

    USE particles_module, ONLY: mom , set_rhop_mom , set_cp_rhop_mom , set_mom ,&
         set_cp_mom , solid_volume_fraction , solid_mass_fraction

    USE plume_module, ONLY: s , r , u , w , mag_u , phi , alpha_inp , beta_inp ,&
         rp , prob_factor , particles_loss , r0

    USE variables, ONLY: gi

    !
    REAL*8, DIMENSION(:), INTENT(OUT) :: rhs_
    !
    INTEGER :: i
    INTEGER :: i_part

    REAL*8 :: ueps

    REAL*8 :: cos_phi
    REAL*8 :: factor0

    REAL*8 :: solid_term , cp_solid_term

    REAL*8 :: Ri

    REAL*8 :: alpha_p , beta_p

    REAL*8 :: A , C

    REAL*8 :: a_poly , b_poly , c_poly , d_poly
    REAL*8 :: s_star

    REAL*8 :: a_10 , a_10_deriv

    !WRITE(*,*) 'mag_u',mag_u

    cos_phi = u / mag_u

    IF ( alpha_inp .LE. 0.d0 ) THEN 

       s_star = s / r0
       
       IF ( s_star .LE. 10 ) THEN
          
          ! value and slope at s_star=10 from equation 3.4 Carazzo et al. 2006
          a_10 = 2.45 - 1.05* exp( -4.65e-3*10.d0)
          
          a_10_deriv = - 1.05* exp( -4.65e-3*10.d0) * ( -4.65e-3)
          
          ! coefficients for the 3rd order polynomial defining A as a function of z/D
          
          a_poly = 1.1
          
          b_poly = 0
                    
          d_poly = - ( a_10 - 1.1 - a_10_deriv ) / 500
          
          c_poly = ( a_10 - 1.1 - 1000*d_poly ) / 100

          ! Equation 12 Carazzo et al. 2008
          A = a_poly + b_poly*s_star + c_poly*s_star**2 + d_poly*s_star**3

       ELSE
        
          ! Equation 3.4 Carazzo et al. 2006
          A = 2.45d0 - 1.05d0* exp( -4.65e-3*s_star )
        
       END IF
    
       C = 0.135d0

       Ri = gi * ( rho_atm - rho_mix ) * r / ( rho_atm * mag_u**2 )

       alpha_p = MAX( 0.d0 , 0.5D0 * C + ( 1.d0 - 1.d0 / A ) * Ri )

       !WRITE(*,*) 's_star,Ri,alpha_p',s_star,Ri,alpha_p

    ELSE

       alpha_p = alpha_inp

    END IF

    IF ( particles_loss ) THEN

       !---- Probability of particle loss (Eq. 15 PlumeMoM - GMD) 
       factor0 = ( 1.D0 + 6.D0 / 5.D0 * alpha_p )** 2
       prob_factor = ( factor0 - 1.D0 ) / ( factor0 + 1.D0 ) 

    ELSE

       prob_factor = 0.D0

    END IF


     IF ( beta_inp .LE. 0.d0 ) THEN 

       beta_p = 0.0D0

    ELSE

       beta_p = beta_inp

    END IF
   
    !---- Entrainment velocity (Eq. 20 PlumeMoM - GMD) 
    ueps = alpha_p * ABS( mag_u - u_atm * cos_phi ) + beta_p * ABS( u_atm *     &
         SIN(phi))

    solid_term = 0.D0

    IF ( distribution_variable .EQ. "particles_number" ) THEN

       solid_term = SUM( solid_volume_fraction(1:n_part) *                      &
            set_rhop_mom(1:n_part,3) )

    ELSEIF ( distribution_variable .EQ. "mass_fraction" ) THEN

       solid_term = rho_mix * SUM( solid_mass_fraction(1:n_part) *              &
            set_mom(1:n_part,0) )

    END IF

    !---- Mass conservation of the plume  (Eq. 20 PlumeMoM - GMD)

    rhs_(1) = 2.D0 * r * rho_atm * ueps - prob_factor * 2.D0 * r *              &
         solid_term 

    !---- Horizontal momentum conservation   (Eq. 21 PlumeMoM - GMD)
    rhs_(2) = - r**2 * rho_mix * w * duatm_dz - u * prob_factor * 2.D0 * r *    &
         solid_term 

    !---- Vertical momentum conservation   (Eq. 22 PlumeMoM - GMD)  
    rhs_(3) = gi * r**2 * ( rho_atm - rho_mix ) - w * prob_factor * 2.D0 * r *  &
         solid_term

    !---- Mixture specific heat integration 

    cp_solid_term = 0.D0

    IF ( distribution_variable .EQ. "particles_number" ) THEN

       cp_solid_term = SUM( solid_volume_fraction(1:n_part) *                   &
            set_cp_rhop_mom(1:n_part,3) )

    ELSEIF ( distribution_variable .EQ. "mass_fraction" ) THEN

       cp_solid_term = rho_mix * SUM( solid_mass_fraction(1:n_part) *           &
            set_cp_mom(1:n_part,0) )

    END IF 

    !----  Mixture specific heat equation RHS term (Eq. 27 PlumeMoM - GMD)
    rhs_(4) = 1.D0 / ( rho_mix * mag_u * r**2 ) * ( - cpmix * rhs_(1)           &
         + cpair * 2.D0 * r * rho_atm * ueps - prob_factor * 2.D0 * r *         &
         cp_solid_term )
 
    !---- Energy conservation    (Eq.7 Bursik 2001) (Eq. 23 PlumeMoM - GMD)
    rhs_(5) = 2.D0 * r * ueps * rho_atm * cpair * ta - ( r**2 ) * w * rho_atm * &
         gi - tp * prob_factor * 2.D0 * r * cp_solid_term - rp * r *            &
         ( tp**4 - ta**4 )

    !---- Gas constant for the mixture integration  (Eq. 29 PlumeMoM - GMD)
    rhs_(6) = ( rair - rgasmix ) / ( rho_mix * gas_mass_fraction * mag_u*r**2 ) &
         * ( 2.D0 * r * rho_atm * ueps )
      
    !WRITE(*,*) '+++++++++++++++++++++++++++++++++++'
    !WRITE(*,*) 'rgasmix',rgasmix,ueps,r,mag_u
    !WRITE(*,*) 'rhs_(6)',rhs_(6)
    !WRITE(*,*) '+++++++++++++++++++++++++++++++++++'
    !READ(*,*) 

    !---- Z integration   (Eq. 30 PlumeMoM - GMD)
    rhs_(7) = w / mag_u

    !---- X integration   (Eq. 30 PlumeMoM - GMD)
    rhs_(8) = cos_phi * cos_theta

    !---- Y integration   (Eq. 30 PlumeMoM - GMD)
    rhs_(9) = cos_phi * sin_theta

    !---- Moments equations
    DO i_part=1,n_part

       DO i=0,n_mom-1

          IF ( distribution_variable .EQ. "particles_number" ) THEN

             !---- Momentum equation RHS term (Eq. 16 PlumeMoM - GMD)
             rhs_(10+i+(i_part-1)*n_mom) = - prob_factor * 2.D0 * r *           &
                  set_mom(i_part,i) * mom(i_part,i) 

          ELSEIF ( distribution_variable .EQ. "mass_fraction" ) THEN

             !---- Momentum equation RHS term (Eq. 32 PlumeMoM - GMD)
             rhs_(10+i+(i_part-1)*n_mom) = - prob_factor * 2.D0 * r *           &
                  rho_mix * set_mom(i_part,i) * mom(i_part,i)

          END IF

       END DO
    
    END DO

    RETURN

  END SUBROUTINE rate


  !******************************************************************************
  !> \brief Calculate the lumped variables
  !
  !> This subroutine calculates the lumped variables f_ (left-hand side of the 
  !> equations of the model).
  !> \param[out]    f_     lumped variables
  !> \date 22/10/2013
  !> @author 
  !> Mattia de' Michieli Vitturi
  !******************************************************************************
 
  SUBROUTINE lump(f_)

    USE mixture_module, ONLY: rgasmix , cpmix , rho_mix , tp
    USE meteo_module, ONLY: u_atm
    USE particles_module, ONLY: mom
    USE plume_module, ONLY: x , z , y , r , u , w , mag_u
    !
    IMPLICIT NONE
    REAL*8, DIMENSION(:), INTENT(OUT) :: f_
    INTEGER :: i_mom
    INTEGER :: i_part
    INTEGER :: idx

    f_(1) = rho_mix * mag_u * r**2
    f_(2) = f_(1) * ( u - u_atm )
    f_(3) = f_(1) * w
    f_(4) = cpmix
    f_(5) = f_(1) * cpmix * tp
    f_(6) = rgasmix
    f_(7) = z
    f_(8) = x
    f_(9) = y

    DO i_part=1,n_part

       DO i_mom=0,n_mom-1

          idx = 10+i_mom+(i_part-1)*n_mom

          IF ( distribution_variable .EQ. "particles_number" ) THEN

             f_(idx) = mag_u * r**2 * mom(i_part,i_mom)

          ELSEIF ( distribution_variable .EQ. "mass_fraction" ) THEN

             f_(idx) = mag_u * r**2 * mom(i_part,i_mom) * rho_mix

          END IF
             
       ENDDO

    END DO

    RETURN

  END SUBROUTINE lump

  !******************************************************************************
  !> \brief Marching s one step
  !
  !> This subroutine update the solution of the model from s to s+ds as 
  !> fnew=fold+ds*rate.
  !> \param[in]    fold    old lumped variables
  !> \param[in]    rate     rate of change of the lumped variables
  !> \param[out]   fnew    new lumped variables
  !> \date 22/10/2013
  !> @author 
  !> Mattia de' Michieli Vitturi
  !******************************************************************************

  SUBROUTINE marching(fold,fnew,rate)

    IMPLICIT NONE
    REAL*8, DIMENSION(:), INTENT(IN) :: fold, rate
    REAL*8, DIMENSION(:), INTENT(OUT) :: fnew
    INTEGER :: i
    INTEGER :: i_part, i_mom

    DO i=1,itotal

       fnew(i) = fold(i) + rate(i) * ds

    END DO

    DO i_part=1,n_part

       DO i_mom=0,n_mom-1

          IF ( fnew(10+i_mom+(i_part-1)*n_mom) .LE. 0.D0 ) THEN

             IF ( distribution_variable .EQ. 'particle_moments' ) THEN

                WRITE(*,*) 'WARNING: negative moment, part',i_part,'mom',i_mom
                READ(*,*)
                
             END IF

          END IF

       ENDDO

    END DO

    RETURN

  END SUBROUTINE marching

  !******************************************************************************
  !> \brief Calculate physical variables from lumped variables
  !
  !> This subroutine calculates a set of physical variables from lumped variables
  !> f_.
  !> \param[in]    f_     lumped variables
  !> \date 22/10/2013
  !> @author 
  !> Mattia de' Michieli Vitturi
  !******************************************************************************
  
  SUBROUTINE unlump(f_)

    USE meteo_module, ONLY: u_atm , rair , pa

    USE mixture_module, ONLY: rho_gas , rgasmix , cpmix , rho_mix , tp ,        &
         gas_volume_fraction , solid_tot_volume_fraction , gas_mass_fraction ,  &
         atm_mass_fraction , wvapour_mass_fraction , rwvapour

    USE particles_module, ONLY : mom , solid_partial_mass_fraction ,            &
         solid_partial_volume_fraction , solid_volume_fraction , distribution , &
         distribution_variable , solid_mass_fraction

    USE plume_module, ONLY: x , z , y , r , u , w , mag_u , phi

    USE moments_module, ONLY: n_nodes 

    USE variables, ONLY : pi_g , verbose_level

    USE meteo_module, ONLY : zmet

    USE moments_module, ONLY : moments_correction_wright
    USE moments_module, ONLY : moments_correction, wheeler_algorithm 
 
    USE particles_module, ONLY : eval_particles_moments 
    USE particles_module, ONLY : particles_density
   
    IMPLICIT NONE

    REAL*8, DIMENSION(:), INTENT(INOUT) :: f_
    !
    REAL*8, DIMENSION(n_part,n_nodes) :: xi , wi , wi_temp
    REAL*8, DIMENSION(n_part,n_nodes) :: part_dens_array

    REAL*8 :: rho_wvapour

    REAL*8 :: rhoB_solid_U_r2(n_part)
    REAL*8 :: rhoB_solid_tot_U_r2

    REAL*8 :: alfa_s_u_r2(1:n_part)
    REAL*8 :: alfa_g_u_r2

    REAL*8 :: alfa_g_wvapour
    REAL*8 :: alfa_g_atm

    REAL*8 :: atm_volume_fraction
    REAL*8 :: wvapour_volume_fraction

    REAL*8 :: u_r2

    REAL*8 :: rho_solid_avg(n_part)
    REAL*8 :: rho_solid_tot_avg

    INTEGER :: j
    INTEGER :: i_part
    INTEGER :: iter

    INTEGER :: idx1 , idx2
    

    REAL*8 :: rho_atm_tp

    phi = ATAN(w/u)
    z = f_(7)
    x = f_(8)
    y = f_(9)

    ! ---- evaluate the new atmospheric density ad u and temperature at z -------

    CALL zmet

    u = u_atm + f_(2)/f_(1)

    w = f_(3)/f_(1)

    mag_u = SQRT( u*u + w*w ) 

    cpmix = f_(4)

    tp = f_(5) / ( f_(1) * cpmix )

    rgasmix = f_(6)

    rho_gas = pa / ( rgasmix * tp )

    rho_wvapour = pa / ( rwvapour * tp )

    rho_atm_tp = pa / ( rair * tp )

    alfa_g_wvapour = ( rho_gas - rho_atm_tp ) / ( rho_wvapour - rho_atm_tp )

    alfa_g_atm = 1.D0 - alfa_g_wvapour

    IF ( verbose_level .GT. 2 ) THEN

       WRITE(*,*) '************** UNLUMP ***************'
       WRITE(*,*) 'rgasmix',rgasmix
    
       WRITE(*,*) 'rho_gas,rho_wvapour,rho_atm_tp',rho_gas,rho_wvapour,         &
            rho_atm_tp

       WRITE(*,*) 'alfa_g_wvapour,alfa_g_atm',alfa_g_wvapour,alfa_g_atm

    END IF


    DO i_part=1,n_part

       idx1 = 10 + 0 + n_mom * ( i_part - 1 )
       idx2 = 10 + n_mom - 1 + n_mom * ( i_part - 1 )

       IF ( distribution_variable .EQ. 'particles_number' ) THEN

          IF ( distribution .EQ. 'constant' ) THEN

             CALL moments_correction( f_(idx1:idx2) , iter )
!             CALL moments_correction_wright( f_(idx1:idx2) )
          ELSE

             CALL moments_correction( f_(idx1:idx2) , iter )

          END IF

       END IF

       IF ( distribution .EQ. 'constant' ) THEN

          CALL wheeler_algorithm( f_(idx1:idx1+1) , distribution , xi(i_part,:),&
               wi_temp(i_part,:) )

       ELSE

          CALL wheeler_algorithm( f_(idx1:idx2) , distribution , xi(i_part,:) , &
               wi_temp(i_part,:) )
 
       END IF

       DO j=1,n_nodes

          part_dens_array(i_part,j) = particles_density( i_part , xi(i_part,j) )

       END DO
       
       
       IF ( distribution_variable .EQ. 'particles_number' ) THEN

          rhoB_solid_U_r2(i_part) = SUM( part_dens_array(i_part,:) *            &
            wi_temp(i_part,:) * xi(i_part,:)**3 ) / 6.D0 * pi_g


          rho_solid_avg(i_part) = SUM( part_dens_array(i_part,:) *              &
               wi_temp(i_part,:) * xi(i_part,:)**3 ) /                          &
               SUM( wi_temp(i_part,:) * xi(i_part,:)**3 )

       ELSEIF ( distribution_variable .EQ. 'mass_fraction' ) THEN

          rhoB_solid_U_r2(i_part) = f_(idx1)

          rho_solid_avg(i_part) = 1.D0 / ( SUM( wi_temp(i_part,:) /             &
               part_dens_array(i_part,:) ) / SUM( wi_temp(i_part,:) ) )

       END IF

       IF ( verbose_level .GT. 0 ) THEN

          WRITE(*,*) 'rhoB_solid_U_r2',idx1,rhoB_solid_U_r2(i_part)
          WRITE(*,*) 'part_dens_array(i_part,:)',part_dens_array(i_part,:)
          WRITE(*,*) 'xi(i_part,:)',xi(i_part,:)
          WRITE(*,*) 'i_part,rho_solid_avg',i_part, rho_solid_avg(i_part)

       END IF

    END DO

    alfa_s_u_r2(1:n_part) = rhoB_solid_U_r2(1:n_part) / rho_solid_avg(1:n_part)

    rhoB_solid_tot_u_r2 = SUM( rhoB_solid_U_r2(1:n_part) )

    alfa_g_u_r2 = ( f_(1) - rhoB_solid_tot_U_r2 ) / rho_gas 

    u_r2 = SUM( alfa_s_u_r2(1:n_part) ) + alfa_g_u_r2

    r = DSQRT( u_r2 / mag_u )

    rho_mix = f_(1) / u_r2 

    IF ( verbose_level .GE. 1 ) THEN

       WRITE(*,*) 'rho_gas',rho_gas
       WRITE(*,*) '*********** SUM(alfa_s_u_r2(1:n_part))',SUM(alfa_s_u_r2(1:n_part))
       WRITE(*,*) 'u_r2,r,mag_u',u_r2,r,mag_u
       WRITE(*,*) 'f_(1),rho_gas',f_(1),rho_gas
       WRITE(*,*) 'rhoB_solid_tot_U_r2',rhoB_solid_tot_U_r2
       WRITE(*,*) 'rho_mix',rho_mix
       WRITE(*,*) ' alfa_g', alfa_g_u_r2/ u_r2
       READ(*,*)

    END IF

    ! -------- solid fractions --------------------------------------------------

    solid_volume_fraction(1:n_part) = alfa_s_u_r2(1:n_part) / u_r2

    solid_tot_volume_fraction = SUM( solid_volume_fraction(1:n_part) )

    solid_partial_volume_fraction(1:n_part) = solid_volume_fraction(1:n_part) / &
         solid_tot_volume_fraction

    solid_partial_mass_fraction = solid_partial_volume_fraction * rho_solid_avg &
         /  SUM( solid_partial_volume_fraction * rho_solid_avg )

    solid_mass_fraction(1:n_part) = solid_volume_fraction(1:n_part) *           &
         rho_solid_avg(1:n_part) / rho_mix

    rho_solid_tot_avg = SUM( solid_partial_volume_fraction * rho_solid_avg )

    ! --------- gas fractions ---------------------------------------------------

    gas_mass_fraction = ( f_(1) - rhoB_solid_tot_u_r2 ) / f_(1)

    gas_volume_fraction = 1.D0 - solid_tot_volume_fraction

    atm_volume_fraction = gas_volume_fraction * alfa_g_atm

    wvapour_volume_fraction = gas_volume_fraction * alfa_g_wvapour

    atm_mass_fraction = atm_volume_fraction * rho_atm_tp / rho_mix

    ! WRITE(*,*) 'SOLVER_RISE: atm_mass_fraction',atm_mass_fraction
    ! READ(*,*)


    wvapour_mass_fraction = wvapour_volume_fraction * rho_wvapour / rho_mix

    ! -------- update the moments -----------------------------------------------

    DO i_part=1,n_part

       idx1 = 10 + 0 + n_mom * ( i_part - 1 )
       idx2 = 10 + n_mom - 1 + n_mom * ( i_part - 1 )

       IF ( distribution_variable .EQ. 'particles_number' ) THEN

          mom(i_part,0:n_mom-1) = f_(idx1:idx2) / ( u_r2 ) 

       ELSEIF ( distribution_variable .EQ. 'mass_fraction' ) THEN

          mom(i_part,0:n_mom-1) = f_(idx1:idx2) / ( rho_mix * u_r2 ) 

       END IF

       IF ( verbose_level .GE. 2 ) THEN

          WRITE(*,*) 'i_part',i_part
          WRITE(*,*) 'f',f_(idx1:idx2)
          WRITE(*,*) 'mom',mom(i_part,0:n_mom-1)

       END IF

       IF ( distribution .EQ. 'constant' ) THEN

          CALL wheeler_algorithm( mom(i_part,0:1) , distribution , xi(i_part,:),&
               wi(i_part,:) )

       ELSE

          CALL wheeler_algorithm(  mom(i_part,0:n_mom-1) , distribution ,       &
               xi(i_part,:) , wi(i_part,:) )

       END IF

    END DO

    ! ------- evaluate the new averaged terms -----------------------------------

    CALL eval_particles_moments( xi , wi )

    IF ( verbose_level .GE. 2 ) THEN
       
       WRITE(*,*) ''
       WRITE(*,*) '************** UNLUMPED VARIABLES **************'
       WRITE(*,*) 'cpmix',cpmix
       WRITE(*,*) 'pres',pa
       WRITE(*,*) 'rgasmix',rgasmix
       WRITE(*,*) 'tp',tp
       WRITE(*,*) 'u,w',u,w
       WRITE(*,*) 'r',r

       WRITE(*,*) ''
       WRITE(*,*) '************** DENSITIES **************'
       WRITE(*,*) 'rho_gas',rho_gas       
       WRITE(*,*) 'rho solids',rho_solid_avg
       WRITE(*,*) 'rho solid tot avg',rho_solid_tot_avg
       WRITE(*,*) 'rho_mix',rho_mix

       WRITE(*,*) ''
       WRITE(*,*) '************** VOLUME FRACTIONS **************'
       WRITE(*,*) 'gas volume fraction',gas_volume_fraction
       WRITE(*,*) 'solid_tot_volume_fraction',solid_tot_volume_fraction       
       WRITE(*,*) 'partial solid volume fraction',solid_partial_volume_fraction 

       WRITE(*,*) ''
       WRITE(*,*) '************** MASS FRACTIONS **************'
       WRITE(*,*) 'solid partial mass fraction',solid_partial_mass_fraction 
       WRITE(*,*) 'solid mass fraction',solid_mass_fraction 
       WRITE(*,*) 'solid_tot_mass_fraction',1.D0 - gas_mass_fraction
       WRITE(*,*) 'gas mass fraction',gas_mass_fraction
       WRITE(*,*) 'atm_mass_fraction',atm_mass_fraction
       WRITE(*,*) 'wvapour_mass_fraction',wvapour_mass_fraction
       WRITE(*,*) ''
       READ(*,*)

    END IF

  END SUBROUTINE unlump

END MODULE solver_module
