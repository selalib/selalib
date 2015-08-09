!> @ingroup maxwell_solvers
!> @brief
!> Module interface to solve Maxwell's equations in 1D
!> @details
!> Contains the abstract class to create a Maxwell solver in 1D.

module sll_m_maxwell_1d_fem
#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_utilities.h"

  use sll_m_maxwell_1d_base

  implicit none
  private
  
  public :: sll_new_maxwell_1d_fem

  type, public, extends(sll_maxwell_1d_base) :: sll_maxwell_1d_fem

     sll_real64 :: domain(2)

   contains
     procedure :: &
          compute_E_from_B => compute_E_from_B_1d_fem!< Solve E and B part of Amperes law with B constant in time
     procedure :: &
          compute_B_from_E => compute_B_from_E_1d_fem!< Solve Faraday equation with E constant in time
     procedure :: &
          compute_E_from_rho => compute_E_from_rho_1d_fem!< Solve E from rho using Poisson
  end type sll_maxwell_1d_fem

contains

  subroutine compute_E_from_B_1d_fem(this, delta_t, field_in, field_out)
    class(sll_maxwell_1d_fem) :: this
    sll_real64, intent(in)     :: delta_t
    sll_real64, intent(in)     :: field_in(:)
    sll_real64, intent(inout)  :: field_out(:)
  end subroutine compute_E_from_B_1d_fem

   subroutine compute_B_from_E_1d_fem(this, delta_t, field_in, field_out)
    class(sll_maxwell_1d_fem) :: this
    sll_real64, intent(in)     :: delta_t
    sll_real64, intent(in)     :: field_in(:)
    sll_real64, intent(inout)  :: field_out(:)
   end subroutine compute_B_from_E_1d_fem

  
   subroutine compute_E_from_rho_1d_fem(this, E, rho )       
     class(sll_maxwell_1d_fem) :: this
     sll_real64,dimension(:),intent(in) :: rho
     sll_real64,dimension(:),intent(out) :: E

     E = 0.0_f64

   end subroutine compute_E_from_rho_1d_fem

   function sll_new_maxwell_1d_fem(domain) result(this)
     sll_real64 :: domain(2)
     class(sll_maxwell_1d_fem), pointer :: this

     ! local variables
     sll_int32 :: ierr

     SLL_ALLOCATE(this, ierr)

     this%domain = domain

   end function sll_new_maxwell_1d_fem


end module sll_m_maxwell_1d_fem
