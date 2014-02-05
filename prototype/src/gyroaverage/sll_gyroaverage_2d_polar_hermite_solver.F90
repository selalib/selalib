!**************************************************************
!  Copyright INRIA
!  Authors : 
!     CALVI project team
!  
!  This code SeLaLib (for Semi-Lagrangian-Library) 
!  is a parallel library for simulating the plasma turbulence 
!  in a tokamak.
!  
!  This software is governed by the CeCILL-B license 
!  under French law and abiding by the rules of distribution 
!  of free software.  You can  use, modify and redistribute 
!  the software under the terms of the CeCILL-B license as 
!  circulated by CEA, CNRS and INRIA at the following URL
!  "http://www.cecill.info". 
!**************************************************************

module sll_module_gyroaverage_2d_polar_hermite_solver
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
use sll_module_gyroaverage_2d_base
use sll_gyroaverage_2d_polar
implicit none


  type,extends(sll_gyroaverage_2d_base) :: gyroaverage_2d_polar_hermite_solver     
  
    type(sll_plan_gyroaverage_polar), pointer                   :: gyro
    sll_int32 :: hermite_case
	! hermite_case
	! 1 : hermite standard
	! 2 : hermite C^1
    ! 3 : hermite C^1 with precomputation
    ! 4 : hermite C^1 with invariance

    contains
      procedure, pass(gyroaverage) :: initialize => &
        initialize_gyroaverage_2d_polar_hermite_solver
      procedure, pass(gyroaverage) :: compute_gyroaverage => &
        compute_gyroaverage_2d_polar_hermite
           
  end type gyroaverage_2d_polar_hermite_solver

contains
  function new_gyroaverage_2d_polar_hermite_solver( &
    eta_min, &
    eta_max, &
    Nc, &
    N_points, &
    interp_degree, &
    hermite_case) &     
    result(gyroaverage)
      
    type(gyroaverage_2d_polar_hermite_solver),pointer :: gyroaverage
    sll_real64, intent(in) :: eta_min(2)
    sll_real64, intent(in) :: eta_max(2)
    sll_int32, intent(in)  :: Nc(2)
    sll_int32, intent(in)  :: N_points  
    sll_int32, intent(in)  :: interp_degree(2)
    sll_int32, optional    :: hermite_case
    sll_int32 :: ierr
      
    SLL_ALLOCATE(gyroaverage,ierr)
    call initialize_gyroaverage_2d_polar_hermite_solver( &
      gyroaverage, &
      eta_min, &
      eta_max, &
      Nc, &
      N_points, &
      interp_degree, &
      hermite_case)
    
  end function new_gyroaverage_2d_polar_hermite_solver
  
  
  subroutine initialize_gyroaverage_2d_polar_hermite_solver( &
    gyroaverage, &
    eta_min, &
    eta_max, &
    Nc, &
    N_points, &
    interp_degree, &
    hermite_case)
    class(gyroaverage_2d_polar_hermite_solver) :: gyroaverage
    sll_real64, intent(in) :: eta_min(2)
    sll_real64, intent(in) :: eta_max(2)
    sll_int32, intent(in)  :: Nc(2)
    sll_int32, intent(in)  :: N_points  
    sll_int32, intent(in)  :: interp_degree(2)
    sll_int32              :: deriv_size
    sll_int32, optional    :: hermite_case
    sll_int32 :: ierr

    
    if(.not.(present(hermite_case)))then
      gyroaverage%hermite_case = 2
    else   
      gyroaverage%hermite_case = hermite_case  
    endif

	deriv_size=9
	select case(gyroaverage%hermite_case)
       case (2,3,4)
			deriv_size=4
	end select		
    
    select case(gyroaverage%hermite_case)
       case (1,2,3,4)
          gyroaverage%gyro => new_plan_gyroaverage_polar_hermite( &
          eta_min, &
          eta_max, &
          Nc, &
          N_points, &
          interp_degree, &
          deriv_size)
       case default
          print *,'#bad value of hermite_case=', gyroaverage%hermite_case
          print *,'#not implemented'
          print *,'#in initialize_gyroaverage_2d_polar_hermite_solver'
          stop
    end select   
           
  end subroutine initialize_gyroaverage_2d_polar_hermite_solver
  

  subroutine compute_gyroaverage_2d_polar_hermite( gyroaverage, larmor_rad, f )
    class(gyroaverage_2d_polar_hermite_solver), target :: gyroaverage
    sll_real64, intent(in) :: larmor_rad
    sll_real64,dimension(:,:),intent(inout) :: f

    select case(gyroaverage%hermite_case)
      case (1)
        call compute_gyroaverage_points_polar_hermite(gyroaverage%gyro,f,larmor_rad)           
      case (2)
        call compute_gyroaverage_points_polar_hermite_c1(gyroaverage%gyro,f,larmor_rad)
      case (3)
        call pre_compute_gyroaverage_polar_hermite_c1(gyroaverage%gyro,larmor_rad)
        call compute_gyroaverage_pre_compute_polar_hermite_c1(gyroaverage%gyro,f)
      case (4)
        call compute_gyroaverage_points_polar_with_invar_hermite_c1(gyroaverage%gyro,f,larmor_rad)
      case default
        print *,'#bad value of hermite_case=', gyroaverage%hermite_case
        print *,'#not implemented'
        print *,'compute_gyroaverage_2d_polar_hermite'
        stop
     end select
    
  end subroutine compute_gyroaverage_2d_polar_hermite
  
end module sll_module_gyroaverage_2d_polar_hermite_solver
