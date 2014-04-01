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



module sll_module_gyroaverage_2d_polar_pade_solver
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
use sll_module_gyroaverage_2d_base
use sll_gyroaverage_2d_polar
implicit none

  type,extends(sll_gyroaverage_2d_base) :: gyroaverage_2d_polar_pade_solver     
  
    type(sll_plan_gyroaverage_polar), pointer        :: gyro

    contains
      procedure, pass(gyroaverage) :: initialize => &
        initialize_gyroaverage_2d_polar_pade_solver
      procedure, pass(gyroaverage) :: compute_gyroaverage => &
        compute_gyroaverage_2d_polar_pade
           
  end type gyroaverage_2d_polar_pade_solver

contains
  function new_gyroaverage_2d_polar_pade_solver( &
    eta_min, &
    eta_max, &
    Nc) &     
    result(gyroaverage)
      
    type(gyroaverage_2d_polar_pade_solver),pointer :: gyroaverage
    sll_real64, intent(in) :: eta_min(2)
    sll_real64, intent(in) :: eta_max(2)
    sll_int32, intent(in)  :: Nc(2)
    sll_int32 :: ierr
      
    SLL_ALLOCATE(gyroaverage,ierr)
    call initialize_gyroaverage_2d_polar_pade_solver( &
      gyroaverage, &
      eta_min, &
      eta_max, &
      Nc)
    
  end function new_gyroaverage_2d_polar_pade_solver
  
  
  subroutine initialize_gyroaverage_2d_polar_pade_solver( &
    gyroaverage, &
    eta_min, &
    eta_max, &
    Nc)
    class(gyroaverage_2d_polar_pade_solver) :: gyroaverage
    sll_real64, intent(in) :: eta_min(2)
    sll_real64, intent(in) :: eta_max(2)
    sll_int32, intent(in)  :: Nc(2)
    sll_int32 :: ierr
    
    gyroaverage%gyro => new_plan_gyroaverage_polar_pade( &
    eta_min, &
    eta_max, &
    Nc) 
           
  end subroutine initialize_gyroaverage_2d_polar_pade_solver
  

  subroutine compute_gyroaverage_2d_polar_pade( gyroaverage, larmor_rad, f)
    class(gyroaverage_2d_polar_pade_solver), target :: gyroaverage
    sll_real64, intent(in) :: larmor_rad
    sll_real64,dimension(:,:),intent(inout) :: f

    call compute_gyroaverage_pade_polar(gyroaverage%gyro,f,larmor_rad)
    
  end subroutine compute_gyroaverage_2d_polar_pade
  
end module sll_module_gyroaverage_2d_polar_pade_solver
