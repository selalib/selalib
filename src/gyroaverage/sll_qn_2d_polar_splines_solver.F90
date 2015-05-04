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



module sll_module_qn_2d_polar_splines_solver
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
use sll_timer
use sll_module_qn_2d_base
use sll_qn_2d_polar
implicit none


  type,extends(sll_qn_2d_base) :: qn_2d_polar_splines_solver     
  
    type(sll_plan_qn_polar), pointer                   :: quasineutral

    contains
      procedure, pass(qn) :: initialize => &
        initialize_qn_2d_polar_splines_solver
      procedure, pass(qn) :: precompute_qn => &
        precompute_qn_2d_polar_splines
      procedure, pass(qn) :: solve_qn => &
        solve_qn_2d_polar_splines
           
  end type qn_2d_polar_splines_solver

contains
  function new_qn_2d_polar_splines_solver( &
    eta_min, &
    eta_max, &
    Nc, &
    N_points, &
    lambda) &     
    result(qn)
      
    type(qn_2d_polar_splines_solver),pointer :: qn
    sll_real64, intent(in) :: eta_min(2)
    sll_real64, intent(in) :: eta_max(2)
    sll_int32, intent(in)  :: Nc(2)
    sll_int32, optional    :: N_points  
    sll_real64, dimension(:), intent(in)  :: lambda
    sll_int32 :: ierr
      
    SLL_ALLOCATE(qn,ierr)
    call initialize_qn_2d_polar_splines_solver( &
      qn, &
      eta_min, &
      eta_max, &
      Nc, &
      N_points, &
      lambda)
    
  end function new_qn_2d_polar_splines_solver
  
  
  subroutine initialize_qn_2d_polar_splines_solver( &
    qn, &
    eta_min, &
    eta_max, &
    Nc, &
    N_points, &
    lambda)
    class(qn_2d_polar_splines_solver) :: qn
    sll_real64, intent(in) :: eta_min(2)
    sll_real64, intent(in) :: eta_max(2)
    sll_int32, intent(in)  :: Nc(2)
    sll_int32, intent(in)  :: N_points  
    sll_real64, dimension(:), intent(in)  :: lambda
    sll_int32 :: ierr


          qn%quasineutral => new_plan_qn_polar_splines( &
          eta_min, &
          eta_max, &
          Nc, &
          N_points, &
          lambda)
 

  end subroutine initialize_qn_2d_polar_splines_solver
  
  
  
  
  subroutine precompute_qn_2d_polar_splines( qn, mu_points, mu_weights, N_mu)
    class(qn_2d_polar_splines_solver), target :: qn
    sll_int32,intent(in) :: N_mu
    sll_real64,dimension(1:N_mu),intent(in) :: mu_points
    sll_real64,dimension(1:N_mu),intent(in) :: mu_weights
    sll_real64,dimension(1:N_mu) :: rho_points
    sll_int32 :: i
    
    
    !call precompute_qn_polar_splines(qn%quasineutral,mu_points,mu_weights,N_mu)
    
do i=1,N_mu
  rho_points(i)=sqrt(2._f64*mu_points(i))
enddo
     
!call precompute_double_gyroaverage_matrix_polar_splines(qn%quasineutral,rho_points,N_mu)
!call precompute_inverse_qn_matrix_polar_splines(qn%quasineutral,mu_points,mu_weights,N_mu)  

!    print *, '#Start time mark t0'
!    call sll_set_time_mark(t0)
!print *,"# STEP 1 : pre_compute_double_gyroaverage_matrix_polar_splines"
!        call precompute_double_gyroaverage_matrix_polar_splines(qn%quasineutral,rho_points,N_mu)!
!        time1 = sll_time_elapsed_since(t0)
!        print *, '#time elapsed for STEP 1 : ',time1
!print *,"# STEP 2 : compute_double_gyroaverage_matrix_polar_splines"
!        call compute_double_gyroaverage_matrix_polar_splines(qn%quasineutral,N_mu)
!        print *, '#reached end of STEP 2'
!        time2 = sll_time_elapsed_since(t0)
!        print *, '#time elapsed for STEP 2 : ',time2-time1
!print *,"# STEP 3 : test => double gyroaverage"
!        call test_double_gyroaverage(qn%quasineutral,rho_points,N_mu)  
!print *,"# STEP 3 bis : test : print => profiles"
!        call print_profiles(qn%quasineutral,mu_points,N_mu)     
!print *,"# STEP 4 : test => qn solver"   
!        call test_solve_qn_polar(qn%quasineutral,mu_points,mu_weights,N_mu)   
!print *,"# STEP 5 : solve_qn_polar"
!        call solve_qn_polar(qn%quasineutral,phi,mu_points,mu_weights,N_mu)
!        time5 = sll_time_elapsed_since(t0)
!        print *, '#time elapsed for STEP 5 : ',time5-time1


  end subroutine precompute_qn_2d_polar_splines
  
  
  
  

  subroutine solve_qn_2d_polar_splines( qn, phi)
    class(qn_2d_polar_splines_solver), target :: qn
    sll_real64,dimension(:,:),intent(inout) :: phi
    
    call solve_qn_polar_new(qn%quasineutral,phi)
 
    
!    print *, '#Start time mark t0'
!    call sll_set_time_mark(t0)
!print *,"# STEP 1 : pre_compute_double_gyroaverage_matrix_polar_splines"
!        call pre_compute_double_gyroaverage_matrix_polar_splines(qn%quasineutral,rho_points,N_mu)
!        time1 = sll_time_elapsed_since(t0)
!        print *, '#time elapsed for STEP 1 : ',time1
!print *,"# STEP 2 : compute_double_gyroaverage_matrix_polar_splines"
!        call compute_double_gyroaverage_matrix_polar_splines(qn%quasineutral,N_mu)
!        print *, '#reached end of STEP 2'
!        time2 = sll_time_elapsed_since(t0)
!        print *, '#time elapsed for STEP 2 : ',time2-time1
!print *,"# STEP 3 : test => double gyroaverage"
!        call test_double_gyroaverage(qn%quasineutral,rho_points,N_mu)  
!print *,"# STEP 3 bis : test : print => profiles"
!        call print_profiles(qn%quasineutral,mu_points,N_mu)     
!print *,"# STEP 4 : test => qn solver"   
!        call test_solve_qn_polar(qn%quasineutral,mu_points,mu_weights,N_mu)   
!print *,"# STEP 5 : solve_qn_polar"
!        call solve_qn_polar(qn%quasineutral,phi,mu_points,mu_weights,N_mu)
!        time5 = sll_time_elapsed_since(t0)
!        print *, '#time elapsed for STEP 5 : ',time5-time1


  end subroutine solve_qn_2d_polar_splines
  
end module sll_module_qn_2d_polar_splines_solver