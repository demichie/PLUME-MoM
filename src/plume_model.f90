!********************************************************************************
!> \mainpage   PLUME-MoM  
!> PLUME-MoM is a numerical code for the steady-state plume model, describing  
!> the rise in the atmosphere of a mixture of gas and volcanic ash during an  
!> eruption. The system of equation is formally the same of that presented in 
!> Barsotti et al. 2008, i.e. the equations we formulate describe the same 
!> conservation principles (conservation of mass, momentum and energy).
!> In this model, instead of assuming a finite number of particle sizes, a 
!> continuous distribution of particles size is considered and the mothod of
!> moments (Marchisio and Fox, 2013) is adopted to track the moments of the 
!> particle size distribution along the eruptive column.
!> Version 1.0:\n
!> 
!> Github project page: http://demichie.github.io/PLUME-MoM/
!> \n
!> \author Mattia de' Michieli Vitturi (*) \n
!> (*) Istituto Nazionale di Geofisica e vulcanologia, sezione di Pisa\n
!>     Via della Faggiola, 36\n
!>     I-56126 Pisa, Italy \n
!>     E-mail: mattia.demichielivitturi@ingv.it \n
!*********************************************************************

!> \brief Main Program 


PROGRAM plume_model
  
  USE inpout, ONLY: initialize , read_inp , check_hysplit
  
  USE inpout, ONLY: open_file_units , close_file_units

  USE rise, ONLY: plumerise

  USE solver_module, ONLY: allocate_matrix
  
  IMPLICIT NONE
 
  REAL*8 :: t1 , t2
 
  CALL cpu_time(t1)

  !
  !***  Initialize the input variables
  !
  CALL initialize
  !
  !***  Read from file the input parameters
  !
  CALL read_inp
  !
  !***  Open the units for output
  !
  CALL open_file_units
  !
  !***  Allocate varaibles for the colum model
  !
  CALL allocate_matrix
  
  !***  Solve the plume model
  CALL plumerise
  
  CALL close_file_units

  CALL check_hysplit
  
  CALL cpu_time(t2)

  WRITE(*,*) 'Time taken by the code was',t2-t1,'seconds'

END PROGRAM plume_model
