#ifndef DOXYGEN_SHOULD_SKIP_THIS

!> @ingroup poisson_solvers
module sll_m_poisson_2d_periodic_solver
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_poisson_2d_base, only: &
    sll_c_poisson_2d_base

#ifdef FFTW
  use sll_m_poisson_2d_periodic_fftw, only: &
    sll_o_new, &
    sll_t_poisson_2d_periodic_fftw, &
    sll_o_solve

#define poisson_2d_periodic sll_t_poisson_2d_periodic_fftw
#else
use sll_m_poisson_2d_periodic_fftpack, only: &
    sll_o_new, &
    sll_t_poisson_2d_periodic_fftpack, &
    sll_o_solve

#define poisson_2d_periodic sll_t_poisson_2d_periodic_fftpack
#endif
  implicit none

  public :: &
    sll_f_new_poisson_2d_periodic_solver

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type,extends(sll_c_poisson_2d_base) :: poisson_2d_periodic_solver

    type(poisson_2d_periodic), pointer :: poiss

  contains

    procedure, pass(poisson) :: sll_o_initialize => &
      initialize_poisson_2d_periodic_solver
    procedure, pass(poisson) :: compute_phi_from_rho => &
      compute_phi_from_rho_2d_periodic
    procedure, pass(poisson) :: compute_E_from_rho => &
      compute_E_from_rho_2d_periodic
!    procedure, pass(poisson) :: compute_E_from_phi => &
!      compute_E_from_phi_2d_polar
      
  end type poisson_2d_periodic_solver

contains

  function sll_f_new_poisson_2d_periodic_solver( &
    eta1_min, &
    eta1_max, &
    nc_eta1, &
    eta2_min, &
    eta2_max, &
    nc_eta2) &     
    result(poisson)
      
    type(poisson_2d_periodic_solver),pointer :: poisson
    sll_real64, intent(in) :: eta1_min
    sll_real64, intent(in) :: eta1_max
    sll_int32,  intent(in) :: nc_eta1
    sll_real64, intent(in) :: eta2_min
    sll_real64, intent(in) :: eta2_max
    sll_int32,  intent(in) :: nc_eta2
    sll_int32 :: ierr
      
    SLL_ALLOCATE(poisson,ierr)
    call initialize_poisson_2d_periodic_solver( &
      poisson, &
      eta1_min, &
      eta1_max, &
      nc_eta1, &
      eta2_min, &
      eta2_max, &
      nc_eta2)     
    
  end function sll_f_new_poisson_2d_periodic_solver
  
  
  subroutine initialize_poisson_2d_periodic_solver( &
    poisson, &
    eta1_min, &
    eta1_max, &
    nc_eta1, &
    eta2_min, &
    eta2_max, &
    nc_eta2)     
    class(poisson_2d_periodic_solver) :: poisson
    sll_real64, intent(in) :: eta1_min
    sll_real64, intent(in) :: eta1_max
    sll_int32,  intent(in) :: nc_eta1
    sll_real64, intent(in) :: eta2_min
    sll_real64, intent(in) :: eta2_max
    sll_int32,  intent(in) :: nc_eta2
    sll_int32 :: ierr

    
    poisson%poiss => sll_o_new( &
      eta1_min, &
      eta1_max, &
      nc_eta1, &
      eta2_min, &
      eta2_max, &
      nc_eta2, &
      ierr)

  end subroutine initialize_poisson_2d_periodic_solver
  
  ! solves -\Delta phi = rho in 2d
  subroutine compute_phi_from_rho_2d_periodic( poisson, phi, rho )
    class(poisson_2d_periodic_solver), target :: poisson
    sll_real64,dimension(:,:),intent(in) :: rho
    sll_real64,dimension(:,:),intent(out) :: phi
    
    call sll_o_solve(poisson%poiss, phi, rho)
    
    
  end subroutine compute_phi_from_rho_2d_periodic

  ! solves E = -\nabla Phi with -\Delta phi = rho in 2d 
  subroutine compute_E_from_rho_2d_periodic( poisson, E1, E2, rho )
    class(poisson_2d_periodic_solver) :: poisson
    sll_real64,dimension(:,:),intent(in) :: rho
    sll_real64,dimension(:,:),intent(out) :: E1
    sll_real64,dimension(:,:),intent(out) :: E2
      
    call sll_o_solve(poisson%poiss, E1, E2, rho)
           
  end subroutine compute_E_from_rho_2d_periodic
  
end module sll_m_poisson_2d_periodic_solver
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
