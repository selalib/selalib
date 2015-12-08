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



module sll_m_gyroaverage_2d_polar_pade_solver
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_gyroaverage_2d_base, only: &
    sll_gyroaverage_2d_base

  use sll_m_gyroaverage_2d_polar, only: &
    compute_gyroaverage_pade_high_order_polar, &
    compute_gyroaverage_pade_polar, &
    new_plan_gyroaverage_polar_pade, &
    sll_plan_gyroaverage_polar

  implicit none

  public :: &
    new_gyroaverage_2d_polar_pade_solver

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type,extends(sll_gyroaverage_2d_base) :: gyroaverage_2d_polar_pade_solver     
  
    type(sll_plan_gyroaverage_polar), pointer                   :: gyro
     sll_int32 :: pade_case
	! pade_case
	! (/0,2/)
	! (/0,4/)
    ! (/2,4/)

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
    Nc, &
    pade_case) &     
    result(gyroaverage)
      
    type(gyroaverage_2d_polar_pade_solver),pointer :: gyroaverage
    sll_real64, intent(in) :: eta_min(2)
    sll_real64, intent(in) :: eta_max(2)
    sll_int32, intent(in)  :: Nc(2)
    sll_int32, optional    :: pade_case(2)
    sll_int32 :: ierr
      
    SLL_ALLOCATE(gyroaverage,ierr)
    call initialize_gyroaverage_2d_polar_pade_solver( &
      gyroaverage, &
      eta_min, &
      eta_max, &
      Nc, &
      pade_case)
    
  end function new_gyroaverage_2d_polar_pade_solver
  
  
  subroutine initialize_gyroaverage_2d_polar_pade_solver( &
    gyroaverage, &
    eta_min, &
    eta_max, &
    Nc, &
    pade_case)
    class(gyroaverage_2d_polar_pade_solver) :: gyroaverage
    sll_real64, intent(in) :: eta_min(2)
    sll_real64, intent(in) :: eta_max(2)
    sll_int32, intent(in)  :: Nc(2)
    sll_int32, optional    :: pade_case(2)
    
    if(.not.(present(pade_case)))then
      gyroaverage%pade_case = 0
    elseif ((pade_case(1)==0).and.(pade_case(2)==2)) then
      gyroaverage%pade_case = 1  
    elseif ((pade_case(1)==0).and.(pade_case(2)==4)) then
      gyroaverage%pade_case = 2  
    elseif ((pade_case(1)==2).and.(pade_case(2)==4)) then
      gyroaverage%pade_case = 3  
    else
      print *,'#bad value of pade_case=', gyroaverage%pade_case
      print *,'#not implemented'
      print *,'#in initialize_gyroaverage_2d_polar_pade_solver'
      stop  
    endif   

        gyroaverage%gyro => new_plan_gyroaverage_polar_pade( &
        eta_min, &
        eta_max, &
        Nc) 
      
           
  end subroutine initialize_gyroaverage_2d_polar_pade_solver
  

  subroutine compute_gyroaverage_2d_polar_pade( gyroaverage, larmor_rad, f)
    class(gyroaverage_2d_polar_pade_solver), target :: gyroaverage
    sll_real64, intent(in) :: larmor_rad
    sll_real64,dimension(:,:),intent(inout) :: f


    select case(gyroaverage%pade_case)
      case (0)
        call compute_gyroaverage_pade_polar(gyroaverage%gyro,f,larmor_rad)
      case (1)
        call compute_gyroaverage_pade_high_order_polar(gyroaverage%gyro,f,larmor_rad,(/0,2/))
      case (2)
        call compute_gyroaverage_pade_high_order_polar(gyroaverage%gyro,f,larmor_rad,(/0,4/))
      case (3)
        call compute_gyroaverage_pade_high_order_polar(gyroaverage%gyro,f,larmor_rad,(/2,4/))
      case default
        print *,'#bad value of pade_case=', gyroaverage%pade_case
        print *,'#not implemented'
        print *,'compute_gyroaverage_2d_polar_pade'
        stop
     end select
    
  end subroutine compute_gyroaverage_2d_polar_pade
  
end module sll_m_gyroaverage_2d_polar_pade_solver
