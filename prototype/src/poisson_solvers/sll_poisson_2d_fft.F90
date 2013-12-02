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

module sll_module_poisson_2d_fft
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
!use sll_boundary_condition_descriptors
use sll_module_poisson_2d_base
use sll_poisson_2D_periodic
implicit none

  type,extends(sll_poisson_2d_base) :: poisson_2d_fft_solver     
  
  type(poisson_2d_periodic), pointer                   :: poiss
  
  contains
    procedure, pass(poisson) :: initialize => &
      initialize_poisson_2d_fft_solver
    procedure, pass(poisson) :: compute_phi_from_rho => &
      compute_phi_from_rho_2d_fft
    procedure, pass(poisson) :: compute_E_from_rho => &
      compute_E_from_rho_2d_fft
!    procedure, pass(poisson) :: compute_E_from_phi => &
!      compute_E_from_phi_2d_fft
      
  end type poisson_2d_fft_solver

contains
  function new_poisson_2d_fft_solver( &
    eta1_min, &
    eta1_max, &
    nc_eta1, &
    eta2_min, &
    eta2_max, &
    nc_eta2) &     
    result(poisson)
      
    type(poisson_2d_fft_solver),pointer :: poisson
    sll_real64 :: eta1_min
    sll_real64 :: eta1_max
    sll_int32 :: nc_eta1
    sll_real64 :: eta2_min
    sll_real64 :: eta2_max
    sll_int32 :: nc_eta2
    sll_int32 :: ierr
      
    SLL_ALLOCATE(poisson,ierr)
    call initialize_poisson_2d_fft_solver( &
    poisson, &
    eta1_min, &
    eta1_max, &
    nc_eta1, &
    eta2_min, &
    eta2_max, &
    nc_eta2)     
    
  end function new_poisson_2d_fft_solver
  
  
  subroutine initialize_poisson_2d_fft_solver( &
    poisson, &
    eta1_min, &
    eta1_max, &
    nc_eta1, &
    eta2_min, &
    eta2_max, &
    nc_eta2)     
    class(poisson_2d_fft_solver) :: poisson
    sll_real64 :: eta1_min
    sll_real64 :: eta1_max
    sll_int32 :: nc_eta1
    sll_real64 :: eta2_min
    sll_real64 :: eta2_max
    sll_int32 :: nc_eta2
    sll_int32 :: ierr
    
    SLL_ALLOCATE(poisson%poiss,ierr)
    
    call initialize( &
      poisson%poiss, &
      eta1_min, &
      eta1_max, &
      nc_eta1, &
      eta2_min, &
      eta2_max, &
      nc_eta2, &
      ierr) 

        
        
  end subroutine initialize_poisson_2d_fft_solver
  
  ! solves -\Delta phi = rho in 2d
  subroutine compute_phi_from_rho_2d_fft( poisson, phi, rho )
    class(poisson_2d_fft_solver) :: poisson
    sll_real64,dimension(:,:),intent(in) :: rho
    sll_real64,dimension(:,:),intent(out) :: phi
    
    call solve( poisson%poiss, phi, rho)
    
  end subroutine compute_phi_from_rho_2d_fft

    ! solves E = -\nabla Phi in 2d
!    subroutine compute_E_from_phi_2d_fft( poisson, phi, E1, E2 )
!      class(poisson_2d_fft_solver) :: poisson
!      sll_real64,dimension(:,:),intent(in) :: phi
!      sll_real64,dimension(:,:),intent(out) :: E1
!      sll_real64,dimension(:,:),intent(out) :: E2
!    end subroutine compute_E_from_phi_2d_fft

    ! solves E = -\nabla Phi with -\Delta phi = rho in 2d 
    subroutine compute_E_from_rho_2d_fft( poisson, rho, E1, E2 )
      class(poisson_2d_fft_solver) :: poisson
      sll_real64,dimension(:,:),intent(in) :: rho
      sll_real64,dimension(:,:),intent(out) :: E1
      sll_real64,dimension(:,:),intent(out) :: E2
      
      call solve( poisson%poiss, E1, E2, rho)
      
    end subroutine compute_E_from_rho_2d_fft
  
  
  
  
end module sll_module_poisson_2d_fft
