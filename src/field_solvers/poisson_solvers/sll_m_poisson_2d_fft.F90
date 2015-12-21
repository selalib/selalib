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

!> @ingroup poisson_solvers
!> @brief
!> Regular cartesian two dimensional mesh with periodic bounday conditions.
!> @details
!> Numerical method uses Fast Fourier Transform and periodic
!> boundary conditions.
!> @snippet poisson_solvers/test_poisson_2d_fft.F90 example
module sll_m_poisson_2d_fft
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_poisson_2d_base, only: &
    sll_c_poisson_2d_base, &
    sll_f_function_of_position

#ifdef FFTW
  use sll_m_poisson_2d_periodic_fftw, only: &
    sll_o_initialize, &
    sll_t_poisson_2d_periodic_fftw, &
    sll_o_solve

#define poisson_2d_periodic sll_t_poisson_2d_periodic_fftw
#else
use sll_m_poisson_2d_periodic_fftpack, only: &
    sll_o_initialize, &
    sll_t_poisson_2d_periodic_fftpack, &
    sll_o_solve

#define poisson_2d_periodic sll_t_poisson_2d_periodic_fftpack
#endif
  implicit none

  public :: &
    sll_f_new_poisson_2d_fft_solver, &
    sll_t_poisson_2d_fft_solver

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type, extends(sll_c_poisson_2d_base) :: sll_t_poisson_2d_fft_solver

    type(poisson_2d_periodic), private, pointer :: solver

  contains

    !> Create the Poisson solver
    procedure, public, pass(poisson) :: initialize => initialize_poisson_2d_fft_solver
    !> Compute potential solving the Poisson equation
    procedure, public, pass(poisson) :: compute_phi_from_rho => compute_phi_from_rho_2d_fft
    !> Compute electric fields solving the Poisson equation
    procedure, public, pass(poisson) :: compute_E_from_rho => compute_E_from_rho_2d_fft

    !> Compute the squared L_2 for given coefficients
    procedure :: &
         l2norm_squared => l2norm_squared_2d_fft
    !> Compute the right hand side from a given function
    procedure :: &
         compute_rhs_from_function => compute_rhs_from_function_2d_fft
    
      
  end type sll_t_poisson_2d_fft_solver


contains

  !> @returns a pointer to the derived type sll_t_poisson_2d_fft_solver.
  function sll_f_new_poisson_2d_fft_solver( &
    eta1_min, &
    eta1_max, &
    nc_eta1, &
    eta2_min, &
    eta2_max, &
    nc_eta2) &     
    result(poisson)
      
    type(sll_t_poisson_2d_fft_solver),pointer :: poisson
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
    
  end function sll_f_new_poisson_2d_fft_solver
  
  
  subroutine initialize_poisson_2d_fft_solver( &
    poisson, &
    eta1_min, &
    eta1_max, &
    nc_eta1, &
    eta2_min, &
    eta2_max, &
    nc_eta2)     
    class(sll_t_poisson_2d_fft_solver) :: poisson
    sll_real64 :: eta1_min
    sll_real64 :: eta1_max
    sll_int32 :: nc_eta1
    sll_real64 :: eta2_min
    sll_real64 :: eta2_max
    sll_int32 :: nc_eta2
    sll_int32 :: ierr
    
    SLL_ALLOCATE(poisson%solver,ierr)
    
    call sll_o_initialize( &
      poisson%solver, &
      eta1_min, &
      eta1_max, &
      nc_eta1, &
      eta2_min, &
      eta2_max, &
      nc_eta2, &
      ierr) 

  end subroutine initialize_poisson_2d_fft_solver
  
  !> solves \f$ -\Delta phi(x,y) = rho (x,y) \f$
  subroutine compute_phi_from_rho_2d_fft( poisson, phi, rho )
    class(sll_t_poisson_2d_fft_solver), target :: poisson
    sll_real64,dimension(:,:),intent(in) :: rho
    sll_real64,dimension(:,:),intent(out) :: phi
    
    call sll_o_solve( poisson%solver, phi, rho)
    
  end subroutine compute_phi_from_rho_2d_fft

  !> @brief
  !> sll_o_solve Poisson equation to compute electric fields
  !> @details
  !> solves 
  !> \f[ 
  !> E(x,y) = -\nabla \phi(x,y) \\
  !> -\Delta \phi(x,y) = \rho(x,y)
  !> \f]
  subroutine compute_E_from_rho_2d_fft( poisson, E1, E2, rho )
    class(sll_t_poisson_2d_fft_solver) :: poisson
    sll_real64,dimension(:,:),intent(in) :: rho
    sll_real64,dimension(:,:),intent(out) :: E1
    sll_real64,dimension(:,:),intent(out) :: E2
      
    call sll_o_solve( poisson%solver, E1, E2, rho)
      
  end subroutine compute_E_from_rho_2d_fft


  function l2norm_squared_2d_fft(poisson, coefs_dofs) result(r)
    class( sll_t_poisson_2d_fft_solver) , intent(in)        :: poisson !< Poisson solver object.
    sll_real64, intent(in)                                     :: coefs_dofs(:,:) !< Values of the coefficient vectors for each DoF
    sll_real64                                     :: r
    
    r = sum(coefs_dofs**2)*poisson%solver%dx* poisson%solver%dy
    
  end function l2norm_squared_2d_fft
  
  subroutine compute_rhs_from_function_2d_fft(poisson, func, coefs_dofs)
    class( sll_t_poisson_2d_fft_solver)                    :: poisson !< Maxwell solver object.
    procedure(sll_f_function_of_position)          :: func !< Function to be projected.
    sll_real64, intent(out)                        :: coefs_dofs(:) !< Coefficients of the projection.
    
    print*, 'compute_rhs_from_function not implemented for poisson_2d_fft.'
    
  end subroutine compute_rhs_from_function_2d_fft



  
end module sll_m_poisson_2d_fft
