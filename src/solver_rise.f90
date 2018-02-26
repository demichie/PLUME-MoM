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
  USE mixture_module, ONLY : n_gas
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

  !> Rate of change of volcanic gases
  REAL*8, ALLOCATABLE, DIMENSION(:) :: volcgas_rate

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
    itotal = n_part * n_mom + 9 + n_gas
    !
    ALLOCATE(f(itotal))
    ALLOCATE(f_stepold(itotal))
    ALLOCATE(rhs(itotal))
    ALLOCATE(rhstemp(itotal))
    ALLOCATE(ftemp(itotal))
    ALLOCATE(volcgas_rate(n_gas))

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
         cos_theta , sin_theta , sphu_atm , c_wv , h_wv0

    USE mixture_module, ONLY: rho_mix , tp , rgasmix , rho_gas ,                &
         gas_mass_fraction , cpvolcgas , rvolcgas

    USE particles_module, ONLY: mom , set_rhop_mom , set_cp_rhop_mom , set_mom ,&
         set_cp_mom , solid_volume_fraction , solid_mass_fraction

    USE plume_module, ONLY: s , r , u , w , mag_u , phi , alpha_inp , beta_inp ,&
         rp , prob_factor , particles_loss , r0 , z

    USE variables, ONLY: gi

    !
    REAL*8, DIMENSION(:), INTENT(OUT) :: rhs_
    !
    INTEGER :: i
    INTEGER :: i_part
    INTEGER :: i_gas

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

       ! See Eq. 34 PlumeMoM- GMD
       solid_term = rho_mix * SUM( set_mom(1:n_part,0) * mom(1:i_part,0) )

    END IF

    DO i_gas=1,n_gas

       volcgas_rate(i_gas) = 0.D0

    END DO

    !---- Mass conservation of the plume  (Eq. 20 PlumeMoM - GMD)

    rhs_(1) = 2.D0 * r * rho_atm * ueps - prob_factor * 2.D0 * r *              &
         solid_term + SUM( volcgas_rate(1:n_gas) )

    !WRITE(*,*) 'SUM( volcgas_rate(1:n_gas) )',  SUM( volcgas_rate(1:n_gas) )
    !READ(*,*)


    !---- Horizontal momentum conservation   (Eq. 21 PlumeMoM - GMD)
    rhs_(2) = - r**2 * rho_mix * w * duatm_dz - u * prob_factor * 2.D0 * r *    &
         solid_term + u * SUM( volcgas_rate(1:n_gas) )

    !---- Vertical momentum conservation   (Eq. 22 PlumeMoM - GMD)  
    rhs_(3) = gi * r**2 * ( rho_atm - rho_mix ) - w * prob_factor * 2.D0 * r *  &
         solid_term + w * SUM( volcgas_rate(1:n_gas) ) 

    !---- Mixture specific heat integration 

    cp_solid_term = 0.D0

    IF ( distribution_variable .EQ. "particles_number" ) THEN

       cp_solid_term = SUM( solid_volume_fraction(1:n_part) *                   &
            set_cp_rhop_mom(1:n_part,3) )

    ELSEIF ( distribution_variable .EQ. "mass_fraction" ) THEN

       cp_solid_term = rho_mix * SUM( solid_mass_fraction(1:n_part) *           &
            set_cp_mom(1:n_part,0) )

    END IF 

    !---- Enthalpy conservation (Eq.7 Bursik 2001)(Eq. 23 PlumeMoM-GMD)
    rhs_(4) = 2.D0 * r * ueps * rho_atm * ( cpair * ta + gi * z                 &
         + 0.5D0 * ueps**2 ) - ( r**2 ) * w * rho_atm * gi                      &
         - tp * prob_factor * 2.D0 * r * cp_solid_term                          &
         - rp * r * ( tp**4 - ta**4 )                                           &
         + tp * SUM( cpvolcgas(1:n_gas) * volcgas_rate(1:n_gas) )

    !WRITE(*,*) 'SUM( cpvolcgas(1:n_gas) * volcgas_rate(1:n_gas) '               &
     !	,SUM( cpvolcgas(1:n_gas) * volcgas_rate(1:n_gas))
    !READ(*,*)   

    !---- Energy conservation    (Eq.2d Folch 2016) 
    rhs_(4) = 2.D0 * r * ueps * rho_atm * ( cpair * ta * ( 1.D0 - sphu_atm )    &
         + sphu_atm * ( h_wv0 - c_wv * ta ) + gi * z                            &
         + 0.5D0 * ueps**2 ) - tp * prob_factor * 2.D0 * r * cp_solid_term      &
         + tp * SUM( cpvolcgas(1:n_gas) * volcgas_rate(1:n_gas) )   


    !---- Z integration   (Eq. 30 PlumeMoM - GMD)
    rhs_(5) = w / mag_u

    !---- X integration   (Eq. 30 PlumeMoM - GMD)
    rhs_(6) = cos_phi * cos_theta

    !---- Y integration   (Eq. 30 PlumeMoM - GMD)
    rhs_(7) = cos_phi * sin_theta

    !---- Moments equations
    DO i_part=1,n_part

       DO i=0,n_mom-1

          IF ( distribution_variable .EQ. "particles_number" ) THEN

             !---- Momentum equation RHS term (Eq. 16 PlumeMoM - GMD)
             rhs_(8+i+(i_part-1)*n_mom) = - prob_factor * 2.D0 * r *            &
                  set_mom(i_part,i) * mom(i_part,i) 

          ELSEIF ( distribution_variable .EQ. "mass_fraction" ) THEN

             !---- Momentum equation RHS term (Eq. 32 PlumeMoM - GMD)
             rhs_(8+i+(i_part-1)*n_mom) = - prob_factor * 2.D0 * r *            &
                  rho_mix * set_mom(i_part,i) * mom(i_part,i)

          END IF

       END DO
    
    END DO

    ! ---- Equations for entrained dry air
    rhs_(n_part*n_mom+8) =  ( 2.D0 * r * rho_atm * ueps ) * ( 1.D0 - sphu_atm )

    ! ---- Equations for H20 (volcanic+entrained)
    rhs_(n_part*n_mom+9) =  ( 2.D0 * r * rho_atm * ueps ) * ( sphu_atm )

    ! ---- Equations for additional volcanic gases 
    DO i_gas=1,n_gas

       rhs_(9+n_part*n_mom+i_gas) = volcgas_rate(i_gas)

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

    USE mixture_module, ONLY: rgasmix , rho_mix , tp
    USE mixture_module, ONLY: volcgas_mass_fraction , volcgas_mix_mass_fraction
    USE mixture_module, ONLY: atm_mass_fraction , mixture_enthalpy
    USE mixture_module, ONLY: dry_air_mass_fraction , water_mass_fraction
    USE mixture_module, ONLY: cpvolcgas_mix  , solid_tot_mass_fraction
    USE mixture_module, ONLY: liquid_water_mass_fraction , water_vapor_mass_fraction

    USE meteo_module, ONLY: u_atm , c_lw , c_wv , cpair , h_lw0 , h_wv0 , T_ref

    USE particles_module, ONLY: mom , cpsolid

    USE plume_module, ONLY: x , z , y , r , u , w , mag_u

    USE variables, ONLY: gi


    !
    IMPLICIT NONE
    REAL*8, DIMENSION(:), INTENT(OUT) :: f_
    INTEGER :: i_mom
    INTEGER :: i_part
    INTEGER :: idx
    INTEGER :: i_gas
    
    f_(1) = rho_mix * mag_u * r**2
    f_(2) = f_(1) * ( u - u_atm )
    f_(3) = f_(1) * w

    mixture_enthalpy = dry_air_mass_fraction * cpair * tp                       &
         + solid_tot_mass_fraction * cpsolid * tp                               & 
         + water_vapor_mass_fraction * ( h_wv0 + c_wv * ( tp - T_ref ) )        &
         + liquid_water_mass_fraction * ( h_lw0 + c_lw * ( tp - T_ref ) )       &
         + volcgas_mix_mass_fraction * cpvolcgas_mix * tp 


    ! ---- Total energy flow rate
    f_(4) = f_(1) * ( mixture_enthalpy + gi * z + 0.5D0 * mag_u**2 ) 
   
    f_(5) = z
    f_(6) = x
    f_(7) = y

    DO i_part=1,n_part

       DO i_mom=0,n_mom-1

          idx = 8+i_mom+(i_part-1)*n_mom

          IF ( distribution_variable .EQ. "particles_number" ) THEN

             f_(idx) = mag_u * r**2 * mom(i_part,i_mom)

          ELSEIF ( distribution_variable .EQ. "mass_fraction" ) THEN

             f_(idx) = mag_u * r**2 * mom(i_part,i_mom) * rho_mix

          END IF
             
       ENDDO

    END DO

    f_(n_part*n_mom+8) = ( rho_mix * dry_air_mass_fraction ) * mag_u * r**2 

    f_(n_part*n_mom+9) = ( rho_mix * water_mass_fraction ) * mag_u * r**2

    DO i_gas=1,n_gas

       f_(9+n_part*n_mom+i_gas) = ( rho_mix * volcgas_mass_fraction(i_gas) )    &
            * mag_u * r**2 

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

          IF ( fnew(8+i_mom+(i_part-1)*n_mom) .LE. 0.D0 ) THEN

             IF ( distribution_variable .EQ. 'particle_moments' ) THEN

                WRITE(*,*) 'WARNING: negative moment, part',i_part,'mom',i_mom
                
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

    USE meteo_module, ONLY: u_atm , rair , pa , cpair , rwv , rho_atm

    USE mixture_module, ONLY: rho_gas , rgasmix , rho_mix , tp ,                &
         gas_volume_fraction , solid_tot_volume_fraction , gas_mass_fraction ,  &
         atm_mass_fraction , rhovolcgas_mix , rvolcgas_mix ,                    &
         volcgas_mass_fraction , volcgas_mix_mass_fraction , cpvolcgas_mix ,    &
         rvolcgas , cpvolcgas , dry_air_mass_fraction , water_mass_fraction ,   &
         solid_tot_mass_fraction , liquid_water_mass_fraction ,                 &
         water_vapor_mass_fraction

    USE particles_module, ONLY : mom , solid_partial_mass_fraction ,            &
         solid_partial_volume_fraction , solid_volume_fraction , distribution , &
         distribution_variable , solid_mass_fraction , cp_rhop_mom , cp_mom ,   &
         rhop_mom , cpsolid

    USE plume_module, ONLY: x , z , y , r , u , w , mag_u , phi

    USE moments_module, ONLY: n_nodes 

    USE variables, ONLY : gi , pi_g , verbose_level

    USE meteo_module, ONLY : zmet

    USE mixture_module, ONLY : eval_temp

    USE moments_module, ONLY : moments_correction_wright
    USE moments_module, ONLY : moments_correction, wheeler_algorithm 
 
    USE particles_module, ONLY : eval_particles_moments 
    USE particles_module, ONLY : particles_density
    
   
    IMPLICIT NONE

    REAL*8, DIMENSION(:), INTENT(INOUT) :: f_
    !
    REAL*8, DIMENSION(n_part,n_nodes) :: xi , wi , wi_temp
    REAL*8, DIMENSION(n_part,n_nodes) :: part_dens_array

    REAL*8 :: rhoB_volcgas_U_r2(n_gas)

    REAL*8 :: rhoB_solid_U_r2(n_part)
    REAL*8 :: rhoB_solid_tot_U_r2

    REAL*8 :: alfa_s_u_r2(1:n_part)
    REAL*8 :: alfa_g_u_r2
    REAL*8 :: alfa_lw_u_r2

    REAL*8 :: atm_volume_fraction
    REAL*8 :: volcgas_mix_volume_fraction

    REAL*8 :: u_r2

    REAL*8 :: rho_solid_avg(n_part)
    REAL*8 :: rho_solid_tot_avg

    INTEGER :: j
    INTEGER :: i_part
    INTEGER :: iter
    INTEGER :: i_gas

    INTEGER :: idx1 , idx2
    
    REAL*8 :: enth

    REAL*8 :: gas_mix_volume_fraction

    ! Mass fraction of water vapor in the mixture
    REAL*8 :: wv_mf

    ! Density of liquid water in the mixture

    REAL*8 :: rho_lw

    ! Volume fraction of liquid water in the mixture

    REAL*8 :: liquid_water_volume_fraction

   
    rho_lw = 1000


    phi = ATAN(w/u)

    
    z = f_(5)
    x = f_(6)
    y = f_(7)

    ! Mass fractions of volcanic gases (H2O excluded ) in mixture of volc. gases
    volcgas_mass_fraction(1:n_gas) = f_(9+n_part*n_mom+1:9+n_part*n_mom+n_gas)  &
         / f_(1) 

    ! Sum of addional gas (H2O excluded) mass fractions
    volcgas_mix_mass_fraction = SUM( volcgas_mass_fraction(1:n_gas) )

    rvolcgas_mix = 0.D0
    cpvolcgas_mix = 0.D0

    ! Properties of the mixture of volcanic gases (H2O excluded)
    IF ( n_gas .GT. 0 ) THEN

       DO i_gas = 1,n_gas
       
          rvolcgas_mix = rvolcgas_mix + volcgas_mass_fraction(i_gas)            &
               * rvolcgas(i_gas)
       
          cpvolcgas_mix = cpvolcgas_mix + volcgas_mass_fraction(i_gas)          &
               * cpvolcgas(i_gas)
        END DO
        rvolcgas_mix = rvolcgas_mix / SUM(volcgas_mass_fraction(1:n_gas)) 
        cpvolcgas_mix = cpvolcgas_mix / SUM(volcgas_mass_fraction(1:n_gas)) 
        
        IF ( verbose_level .GE. 1 ) THEN
           
           WRITE(*,*) 'rvolcgas_mix :', rvolcgas_mix
           WRITE(*,*) 'cpvolcgas_mix :', cpvolcgas_mix
           
        END IF

    ELSE
        
       rvolcgas_mix=0 
       cpvolcgas_mix=0

    END IF
    


    ! mass fraction of dry air in the mixture
    dry_air_mass_fraction = f_(8+n_part*n_mom) / f_(1) 

    ! mass fraction of water in the mixture
    water_mass_fraction = f_(9+n_part*n_mom) / f_(1)

    ! solid mass fraction in the mixture
    solid_tot_mass_fraction = 1.D0- dry_air_mass_fraction - water_mass_fraction &
         - volcgas_mix_mass_fraction
    
    DO i_part=1,n_part

       idx1 = 8 + 0 + n_mom * ( i_part - 1 )
       idx2 = 8 + n_mom - 1 + n_mom * ( i_part - 1 )

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

       IF ( verbose_level .GT. 2 ) THEN

          WRITE(*,*) 'rhoB_solid_U_r2',idx1,rhoB_solid_U_r2(i_part)
          WRITE(*,*) 'part_dens_array(i_part,:)',part_dens_array(i_part,:)
          WRITE(*,*) 'xi(i_part,:)',xi(i_part,:)
          WRITE(*,*) 'i_part,rho_solid_avg',i_part, rho_solid_avg(i_part)

       END IF

    END DO


    ! ---- evaluate the new atmospheric density ad u and temperature at z -------

    CALL zmet

    u = u_atm + f_(2)/f_(1)

    w = f_(3)/f_(1)

    mag_u = SQRT( u*u + w*w ) 

    IF ( distribution_variable .EQ. "particles_number" ) THEN

       cpsolid = ( SUM( rhoB_solid_U_r2(1:n_part) * cp_rhop_mom(1:n_part,3)     &
            / rhop_mom(1:n_part,3) ) ) / ( SUM( rhoB_solid_U_r2(1:n_part) ) ) 

    ELSEIF ( distribution_variable .EQ. "mass_fraction" ) THEN

       cpsolid = ( SUM( rhoB_solid_U_r2(1:n_part) * cp_mom(1:n_part,0) ) )      &
            / ( SUM( rhoB_solid_U_r2(1:n_part) ) ) 

    END IF 


    enth =  f_(4) / f_(1) - gi * z - 0.5D0 * mag_u**2 

    ! --- Compute  water vapor mass fraction from other variables --------------
    CALL eval_temp(enth,pa,cpsolid,tp,wv_mf)

    ! mass fraction of water vapor in the mixture
    water_vapor_mass_fraction = wv_mf
    
    ! mass fraction of liquid water in the mixture    
    liquid_water_mass_fraction = water_mass_fraction - wv_mf
    
    ! constant for mixture of dry air + water vapor + other volcanic gases 
    rgasmix = ( f_(8+n_part*n_mom) * rair + wv_mf * f_(1) * rwv                 &
         + volcgas_mix_mass_fraction * f_(1) * rvolcgas_mix )                   &
         / ( f_(8+n_part*n_mom) + f_(1) * ( wv_mf + volcgas_mix_mass_fraction ) )

    ! density of mixture of dry air + water vapor + other volcanic gases 
    rho_gas = pa / ( rgasmix * tp )

    ! density of mixture of other volcanic gases (no H2O)
    rhovolcgas_mix = pa / ( rvolcgas_mix * tp )

    IF ( verbose_level .GT. 2 ) THEN

       WRITE(*,*) '************** UNLUMP ***************'
       WRITE(*,*) 'rgasmix',rgasmix
    
       WRITE(*,*) 'rho_gas,rhovolcgas_mix',rho_gas,rhovolcgas_mix

    END IF

    alfa_s_u_r2(1:n_part) = rhoB_solid_U_r2(1:n_part) / rho_solid_avg(1:n_part)

    rhoB_solid_tot_u_r2 = SUM( rhoB_solid_U_r2(1:n_part) )

    alfa_g_u_r2 = ( f_(1) * ( 1.D0 - liquid_water_mass_fraction ) -             &
         rhoB_solid_tot_U_r2 ) / rho_gas 

    alfa_lw_u_r2 = f_(1) * liquid_water_mass_fraction / rho_lw

    u_r2 = SUM( alfa_s_u_r2(1:n_part) ) + alfa_g_u_r2 + alfa_lw_u_r2 

    r = DSQRT( u_r2 / mag_u )

    rho_mix = f_(1) / u_r2 

    IF ( verbose_level .GE. 2 ) THEN

       IF ( liquid_water_mass_fraction .GT. 0.d0 ) THEN
          
          WRITE(*,*) '*********** SUM(alfa_s(1:n_part))' ,                      &
               SUM(alfa_s_u_r2(1:n_part))/u_r2
          WRITE(*,*) ' alfa_g', alfa_g_u_r2/ u_r2
          WRITE(*,*) ' alfa_lw', alfa_lw_u_r2/ u_r2
          WRITE(*,*) ( alfa_lw_u_r2 + alfa_g_u_r2 + SUM(alfa_s_u_r2(1:n_part)) )&
               / u_r2
          
       END IF


       WRITE(*,*) 'rho_gas',rho_gas
       WRITE(*,*) '*********** SUM(alfa_s_u_r2(1:n_part))' ,                    &
            SUM(alfa_s_u_r2(1:n_part))
       WRITE(*,*) 'u_r2,r,mag_u',u_r2,r,mag_u
       WRITE(*,*) 'f_(1),rho_gas',f_(1),rho_gas
       WRITE(*,*) 'rhoB_solid_tot_U_r2',rhoB_solid_tot_U_r2
       WRITE(*,*) 'rho_mix',rho_mix
       WRITE(*,*) 'alfa_g', alfa_g_u_r2 / u_r2
       READ(*,*)

    END IF

 
    ! --------- liquid fractions ------------------------------------------------

    liquid_water_volume_fraction = liquid_water_mass_fraction * ( rho_mix       &
         / rho_lw)

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
    ! --------- mixture of dry air + water vapor + other volcanic gases ---------
    
    gas_mass_fraction = ( f_(1)  * ( 1.D0 - liquid_water_mass_fraction ) -      &
         rhoB_solid_tot_u_r2 ) / f_(1)

    gas_volume_fraction = 1.D0 - solid_tot_volume_fraction -                    &
         liquid_water_volume_fraction

    volcgas_mix_volume_fraction = volcgas_mix_mass_fraction * ( rho_mix /       &
         rhovolcgas_mix )

    ! -------- update the moments -----------------------------------------------

    DO i_part=1,n_part

       idx1 = 8 + 0 + n_mom * ( i_part - 1 )
       idx2 = 8 + n_mom - 1 + n_mom * ( i_part - 1 )

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

    IF ( verbose_level .GE. 1 ) THEN
       
       WRITE(*,*) ''
       WRITE(*,*) '************** UNLUMPED VARIABLES **************'
       WRITE(*,*) 'pres',pa
       WRITE(*,*) 'rgasmix',rgasmix
       WRITE(*,*) 'tp',tp
       WRITE(*,*) 'u,w',u,w
       WRITE(*,*) 'r',r
       WRITE(*,*) 'z',z
       WRITE(*,*) 'mass_flow_rate', pi_g * rho_mix * mag_u * (r**2)
       
       WRITE(*,*) ''
       WRITE(*,*) '************** DENSITIES **************'
       WRITE(*,*) 'rho_gas',rho_gas       
       WRITE(*,*) 'rho solids',rho_solid_avg
       WRITE(*,*) 'rho solid tot avg',rho_solid_tot_avg
       WRITE(*,*) 'rho_atm',rho_atm
       WRITE(*,*) 'rho_mix',rho_mix

       WRITE(*,*) ''
       WRITE(*,*) '************** VOLUME FRACTIONS **************'
       WRITE(*,*) 'solid partial volume fractions',solid_partial_volume_fraction 
       WRITE(*,*) 'solid tot volume fraction',solid_tot_volume_fraction       
       WRITE(*,*) 'liquid water volume fraction',liquid_water_volume_fraction
       WRITE(*,*) 'gas volume fraction',gas_volume_fraction
       WRITE(*,*) 'sum of previous three volume fractions',                     &
            solid_tot_volume_fraction + liquid_water_volume_fraction +          &
            gas_volume_fraction

       WRITE(*,*) ''
       WRITE(*,*) '************** MASS FRACTIONS **************'
       WRITE(*,*) 'solid partial mass fractions',solid_partial_mass_fraction 
       WRITE(*,*) 'solid tot mass fraction',solid_tot_mass_fraction 
       WRITE(*,*) 'liquid water mass fraction',liquid_water_mass_fraction
       WRITE(*,*) 'gas mass fraction',gas_mass_fraction
       WRITE(*,*) 'sum of previous three mass fractions',                       &
            solid_tot_mass_fraction + liquid_water_mass_fraction +              &
            gas_mass_fraction

       WRITE(*,*) 

       WRITE(*,*) 'volcgas_mass_fractions',volcgas_mass_fraction
       WRITE(*,*) 'volcgas_mix_mass_fraction',volcgas_mix_mass_fraction
       WRITE(*,*) 'water vapor mass fraction',water_vapor_mass_fraction
       WRITE(*,*) 'dry air mass fraction',dry_air_mass_fraction
       WRITE(*,*) 'sum of previous three mass fractions',                       &
            volcgas_mix_mass_fraction + water_vapor_mass_fraction +             &
            dry_air_mass_fraction
       WRITE(*,*) 'liquid water mass fraction',liquid_water_mass_fraction
       
       WRITE(*,*) ''
       READ(*,*)

    END IF

  END SUBROUTINE unlump

END MODULE solver_module
