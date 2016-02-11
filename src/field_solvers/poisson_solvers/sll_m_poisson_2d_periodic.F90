#ifndef DOXYGEN_SHOULD_SKIP_THIS

!> @ingroup poisson_solvers
module sll_m_poisson_2d_periodic
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_poisson_2d_base, only: &
    sll_c_poisson_2d_base, &
    sll_f_function_of_position

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
    sll_f_new_poisson_2d_periodic

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

    
    !> Compute the squarred L_2 for given coefficients
    procedure :: &
         l2norm_squared => l2norm_squarred_2d_periodic
    !> Compute the right hand side from a given function
    procedure :: &
         compute_rhs_from_function => compute_rhs_from_function_2d_periodic
    !> Destructor
    procedure :: free => delete_2d_periodic
      
  end type poisson_2d_periodic_solver

contains

  function sll_f_new_poisson_2d_periodic( &
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
    
  end function sll_f_new_poisson_2d_periodic
  
  
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
  
  function l2norm_squarred_2d_periodic(poisson, coefs_dofs) result(r)
    class( poisson_2d_periodic_solver) , intent(in)        :: poisson !< Poisson solver object.
    sll_real64   , intent(in)                                  :: coefs_dofs(:,:) !< Values of the coefficient vectors for each DoF
    sll_real64                                     :: r
    
    r = 0.0_f64
    print*, 'l2norm_squared not implemented for poisson_2d_periodic_solver.'
    
  end function l2norm_squarred_2d_periodic
  
  subroutine compute_rhs_from_function_2d_periodic(poisson, func, coefs_dofs)
    class( poisson_2d_periodic_solver)                    :: poisson !< Maxwell solver object.
    procedure(sll_f_function_of_position)          :: func !< Function to be projected.
    sll_real64, intent(out)                        :: coefs_dofs(:) !< Coefficients of the projection.
    
    print*, 'compute_rhs_from_function not implemented for poisson_2d_periodic_solver.'
    
  end subroutine compute_rhs_from_function_2d_periodic
  
  subroutine delete_2d_periodic(poisson)
    class( poisson_2d_periodic_solver)                    :: poisson !< Maxwell solver object.
  end subroutine delete_2d_periodic


end module sll_m_poisson_2d_periodic
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
