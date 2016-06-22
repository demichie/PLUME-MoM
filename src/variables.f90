!*****************************************************************************
!>\brief Global variables
!
!> This module contains global variables used in the other modules. 
!> \date 23/11/2008
!> @author 
!> Mattia de' Michieli Vitturi
!*****************************************************************************   

MODULE variables

  IMPLICIT NONE

  !> Gravity acceleration 
  REAL*8 :: gi          

  !> Greek pi  
  REAL*8 :: pi_g        

  !> Level of verbose output (0 = minimal output on screen)
  INTEGER :: verbose_level

  !> Flag for dakota run (less files on output)
  LOGICAL :: dakota_flag

  !> Flag for hysplit run 
  LOGICAL :: hysplit_flag

  !> Flag for hysplit output\n
  !> - '.TRUE.'          => last point of emission at neutral bouyancy level
  !> - '.FALSE.'         => last point of emission at maximum plume height
  !> .
  LOGICAL :: nbl_stop

  !> Maximum number of particle phases
  INTEGER, PARAMETER :: max_n_part = 50

  LOGICAL :: inversion_flag

  REAL*8 :: height_weight 
  REAL*8 :: height_obj 
  REAL*8 :: mu_weight
  REAL*8 :: mu_obj
  REAL*8 :: sigma_weight 
  REAL*8 :: sigma_obj 
  REAL*8 :: skew_weight 
  REAL*8 :: skew_obj

END MODULE variables
