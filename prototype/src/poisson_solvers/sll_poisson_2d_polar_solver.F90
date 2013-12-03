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


!solves \sum_{i,j=1}^2 A_{i,j}\partial_{i,j} phi
!       +\sum_{i=1}^2B_i\partial_i phi
!       +C \phi = rho
!in polar coordinates
!this leads when A_{1,2}=A_{2,1}=0 and B_2 = 0
! A_11\partial_{1,e1}\hat{phi}+B_1\partial_{1}\hat{phi}+(C+A_{2,2}k^2)\hat{phi} = \hat{rho}


module sll_module_poisson_2d_polar_solver
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
!use sll_boundary_condition_descriptors
use sll_module_poisson_2d_base
use sll_poisson_2d_polar
implicit none

  type,extends(sll_poisson_2d_base) :: poisson_2d_polar_solver     
  
  !type(sll_plan_poisson_polar), pointer                   :: poiss
  
  contains
    procedure, pass(poisson) :: initialize => &
      initialize_poisson_2d_polar_solver
    procedure, pass(poisson) :: compute_phi_from_rho => &
      compute_phi_from_rho_2d_polar
    procedure, pass(poisson) :: compute_E_from_rho => &
      compute_E_from_rho_2d_polar
!    procedure, pass(poisson) :: compute_E_from_phi => &
!      compute_E_from_phi_2d_polar
      
  end type poisson_2d_polar_solver

contains
  function new_poisson_2d_polar_solver( &
    eta1_min, &
    eta1_max, &
    nc_eta1, &
    nc_eta2, &
    bc, &
    dlog_density, &
    inv_Te) &     
    result(poisson)
      
    type(poisson_2d_polar_solver),pointer :: poisson
    sll_real64, intent(in) :: eta1_min
    sll_real64, intent(in) :: eta1_max
    sll_int32, intent(in) :: nc_eta1
    sll_int32, intent(in) :: nc_eta2
    sll_int32, intent(in) :: bc(2)
    sll_real64, dimension(:), intent(in), optional :: dlog_density
    sll_real64, dimension(:), intent(in), optional :: inv_Te
    sll_int32 :: ierr
      
    SLL_ALLOCATE(poisson,ierr)
    call initialize_poisson_2d_polar_solver( &
      poisson, &
      eta1_min, &
      eta1_max, &
      nc_eta1, &
      nc_eta2, &
      bc, &
      dlog_density, &
      inv_Te)
    
  end function new_poisson_2d_polar_solver
  
  
  subroutine initialize_poisson_2d_polar_solver( &
    poisson, &
    eta1_min, &
    eta1_max, &
    nc_eta1, &
    nc_eta2, &
    bc, &
    dlog_density, &
    inv_Te)
    class(poisson_2d_polar_solver) :: poisson
    sll_real64, intent(in) :: eta1_min
    sll_real64, intent(in) :: eta1_max
    sll_int32, intent(in) :: nc_eta1
    sll_int32, intent(in) :: nc_eta2
    sll_int32, intent(in) :: bc(2)
    sll_real64, dimension(:), intent(in), optional :: dlog_density
    sll_real64, dimension(:), intent(in), optional :: inv_Te
    sll_int32 :: ierr
    sll_real64 :: delta_eta
    
!    delta_eta = (eta1_max-eta1_min)/real(nc_eta1,f64)
!    
!    poisson%poiss => new_plan_poisson_polar( &
!      delta_eta,& 
!      eta1_min, &
!      nc_eta1, &
!      nc_eta2, &
!      bc, &
!      dlog_density, &
!      inv_Te)

        
  end subroutine initialize_poisson_2d_polar_solver
  
  ! solves -\Delta phi = rho in 2d
  subroutine compute_phi_from_rho_2d_polar( poisson, phi, rho )
    class(poisson_2d_polar_solver) :: poisson
    sll_real64,dimension(:,:),intent(in) :: rho
    sll_real64,dimension(:,:),intent(out) :: phi
    
    !call solve( poisson%poiss, rho, phi)
    
  end subroutine compute_phi_from_rho_2d_polar

    ! solves E = -\nabla Phi in 2d
!    subroutine compute_E_from_phi_2d_fft( poisson, phi, E1, E2 )
!      class(poisson_2d_fft_solver) :: poisson
!      sll_real64,dimension(:,:),intent(in) :: phi
!      sll_real64,dimension(:,:),intent(out) :: E1
!      sll_real64,dimension(:,:),intent(out) :: E2
!    end subroutine compute_E_from_phi_2d_fft

    ! solves E = -\nabla Phi with -\Delta phi = rho in 2d 
    subroutine compute_E_from_rho_2d_polar( poisson, rho, E1, E2 )
      class(poisson_2d_polar_solver) :: poisson
      sll_real64,dimension(:,:),intent(in) :: rho
      sll_real64,dimension(:,:),intent(out) :: E1
      sll_real64,dimension(:,:),intent(out) :: E2
      
      print *,'#compute_E_from_rho_2d_polar'      
      print *,'#not implemented for the moment'
      stop
      
      !call solve( poisson%poiss, E1, E2, rho)
      
    end subroutine compute_E_from_rho_2d_polar
  
  
  
  
end module sll_module_poisson_2d_polar_solver
