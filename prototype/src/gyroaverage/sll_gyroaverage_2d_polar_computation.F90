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



module sll_module_gyroaverage_2d_polar_computation
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
use sll_module_gyroaverage_2d_base
use sll_gyroaverage_2d_polar
implicit none

  sll_int32, parameter :: SLL_GYROAVERAGE_PADE = 0
  
  sll_int32, parameter :: SLL_GYROAVERAGE_HERMITE = 10
  sll_int32, parameter :: SLL_GYROAVERAGE_HERMITE_C1 = 11
  sll_int32, parameter :: SLL_GYROAVERAGE_HERMITE_C1_PRECOMPUTE = 12
  sll_int32, parameter :: SLL_GYROAVERAGE_HERMITE_C1_WITH_INVARIANCE = 13
  
  sll_int32, parameter :: SLL_GYROAVERAGE_SPLINES = 20
  sll_int32, parameter :: SLL_GYROAVERAGE_SPLINES_PRECOMPUTE = 21
  sll_int32, parameter :: SLL_GYROAVERAGE_SPLINES_WITH_INVARIANCE = 22
  sll_int32, parameter :: SLL_GYROAVERAGE_SPLINES_PRECOMPUTE_WITH_FFT = 23

  type,extends(sll_gyroaverage_2d_base) :: gyroaverage_2d_polar_computation     
  
    type(sll_plan_gyroaverage_polar), pointer                   :: gyro
    sll_int32 :: gyroaverage_case
    sll_real64  :: eta_min(2)
    sll_real64 :: eta_max(2)
    sll_int32 :: Nc(2)
    sll_int32 :: N_points  
    sll_int32 :: interp_degree(2)

    contains
      procedure, pass(gyroaverage) :: initialize => &
        initialize_gyroaverage_2d_polar_computation
      procedure, pass(gyroaverage) :: compute_gyroaverage => &
        compute_gyroaverage_2d_polar
           
  end type gyroaverage_2d_polar_computation

contains
  function new_gyroaverage_2d_polar_computation( &
    eta_min, &
    eta_max, &
    Nc, &
    N_points, &
    interp_degree, &
    gyroaverage_case) &     
    result(gyroaverage)
      
    type(gyroaverage_2d_polar_computation),pointer :: gyroaverage
    sll_real64, intent(in) :: eta_min(2)
    sll_real64, intent(in) :: eta_max(2)
    sll_int32, intent(in)  :: Nc(2)
    sll_int32, optional    :: N_points  
    sll_int32, optional    :: interp_degree(2)
    sll_int32, optional    :: gyroaverage_case
    sll_int32 :: ierr
      
    SLL_ALLOCATE(gyroaverage,ierr)
    call initialize_gyroaverage_2d_polar_computation( &
      gyroaverage, &
      eta_min, &
      eta_max, &
      Nc, &
      N_points, &
      interp_degree, &
      gyroaverage_case)
    
  end function new_gyroaverage_2d_polar_computation
  
  
  subroutine initialize_gyroaverage_2d_polar_computation( &
    gyroaverage, &
    eta_min, &
    eta_max, &
    Nc, &
    N_points, &
    interp_degree, &
    gyroaverage_case)
    class(gyroaverage_2d_polar_computation) :: gyroaverage
    sll_real64, intent(in) :: eta_min(2)
    sll_real64, intent(in) :: eta_max(2)
    sll_int32, intent(in)  :: Nc(2)
    sll_int32, optional    :: N_points  
    sll_int32, optional    :: interp_degree(2)
    sll_int32, optional    :: gyroaverage_case
    sll_int32 :: ierr
    
    
    gyroaverage%eta_min=eta_min
    gyroaverage%eta_max=eta_max
    gyroaverage%Nc=Nc
    
    if(.not.(present(N_points)))then
      gyroaverage%N_points = 4
    else   
      gyroaverage%N_points = N_points  
    endif
    
    
    if(.not.(present(interp_degree)))then
      gyroaverage%interp_degree(1) = 3
      gyroaverage%interp_degree(2) = 3
    else   
      gyroaverage%interp_degree = interp_degree 
    endif
    
    
    if(.not.(present(gyroaverage_case)))then
      gyroaverage%gyroaverage_case = SLL_GYROAVERAGE_HERMITE
    else      
      select case(gyroaverage_case)
        case (SLL_GYROAVERAGE_PADE)
        case (SLL_GYROAVERAGE_HERMITE)
        case (SLL_GYROAVERAGE_HERMITE_C1)
        case (SLL_GYROAVERAGE_HERMITE_C1_PRECOMPUTE)
        case (SLL_GYROAVERAGE_HERMITE_C1_WITH_INVARIANCE)
        case (SLL_GYROAVERAGE_SPLINES)
        case (SLL_GYROAVERAGE_SPLINES_PRECOMPUTE)
        case (SLL_GYROAVERAGE_SPLINES_WITH_INVARIANCE)
        case (SLL_GYROAVERAGE_SPLINES_PRECOMPUTE_WITH_FFT)   
        case default
          print *,'#bad value of gyroaverage_case=', gyroaverage%gyroaverage_case
          print *,'#not implemented'
          print *,'#in initialize_gyroaverage_2d_polar_computation'
          stop
      end select   
      gyroaverage%gyroaverage_case = gyroaverage_case
    endif
    
    
    
       
  end subroutine initialize_gyroaverage_2d_polar_computation
  

  subroutine compute_gyroaverage_2d_polar( gyroaverage, larmor_rad, f, Jf  )
    class(gyroaverage_2d_polar_computation), target :: gyroaverage
    sll_real64, intent(in) :: larmor_rad
    sll_real64,dimension(:,:),intent(in) :: f
    sll_real64,dimension(:,:),intent(out) :: Jf

    select case(gyroaverage%gyroaverage_case)
!      case (SLL_GYROAVERAGE_PADE)
!        call compute_gyroaverage_pade_polar(gyroaverage%gyro,f,larmor_rad)
!      case (SLL_GYROAVERAGE_HERMITE)
!        call compute_gyroaverage_points_polar_hermite(gyroaverage%gyro,f,larmor_rad)           
!      case (SLL_GYROAVERAGE_HERMITE_C1)
!        call compute_gyroaverage_points_polar_hermite_c1(gyroaverage%gyro,f,larmor_rad)
!      case (SLL_GYROAVERAGE_HERMITE_C1_PRECOMPUTE)
!        call pre_compute_gyroaverage_polar_hermite_c1(gyroaverage%gyro,larmor_rad)
!        call compute_gyroaverage_pre_compute_polar_hermite_c1(gyroaverage%gyro,f)
!      case (SLL_GYROAVERAGE_HERMITE_C1_WITH_INVARIANCE)
!        call compute_gyroaverage_points_polar_with_invar_hermite_c1(gyroaverage%gyro,f,larmor_rad)
!      case (SLL_GYROAVERAGE_SPLINES)
!        call compute_gyroaverage_points_polar_spl(gyroaverage%gyro,f,larmor_rad)
!      case (SLL_GYROAVERAGE_SPLINES_PRECOMPUTE)
!        call pre_compute_gyroaverage_polar_spl(gyroaverage%gyro,larmor_rad)
!        call compute_gyroaverage_pre_compute_polar_spl(gyroaverage%gyro,f)
!      case (SLL_GYROAVERAGE_SPLINES_WITH_INVARIANCE)
!        call compute_gyroaverage_points_polar_with_invar_spl(gyroaverage%gyro,f,larmor_rad)
!      case (SLL_GYROAVERAGE_SPLINES_PRECOMPUTE_WITH_FFT)
!        call pre_compute_gyroaverage_polar_spl_FFT(gyroaverage%gyro,larmor_rad)
!        call compute_gyroaverage_pre_compute_polar_spl_FFT(gyroaverage%gyro,f)
      case default
        print *,'#bad value of gyroaverage_case=', gyroaverage%gyroaverage_case
        print *,'#not implemented'
        print *,'in compute_gyroaverage_2d_polar'
        stop
     end select 
     
     Jf=f  
    
  end subroutine compute_gyroaverage_2d_polar
  
end module sll_module_gyroaverage_2d_polar_computation
