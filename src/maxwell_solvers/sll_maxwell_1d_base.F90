!> @ingroup maxwell_solvers
!> @brief
!> Module interface to solve Maxwell's equations in 1D
!> @details
!> Contains the abstract class to create a Maxwell solver in 1D.

module sll_m_maxwell_1d_base
#include "sll_working_precision.h"

  implicit none
  private
  
  type, public, abstract :: sll_maxwell_1d_base

   contains
     procedure(compute_field1_from_field2), deferred :: &
          compute_E_from_B !< Solve E and B part of Amperes law with B constant in time
     procedure(compute_field1_from_field2), deferred :: &
          compute_B_from_E !< Solve Faraday equation with E constant in time
     procedure(signature_compute_E_from_rho_1d), deferred :: &
          compute_E_from_rho !< Solve E from rho using Poisson
     !procedure(signature_solve), deferred :: &
     !     solve !< Solve Amperes law and Faraday equation
  end type sll_maxwell_1d_base

  abstract interface 
     subroutine compute_field1_from_field2(this, delta_t, field_in, field_out)
     use sll_working_precision
     import sll_maxwell_1d_base
     
     class(sll_maxwell_1d_base) :: this
     sll_real64, intent(in)     :: delta_t
     sll_real64, intent(in)     :: field_in(:)
     sll_real64, intent(inout)  :: field_out(:)
   end subroutine compute_field1_from_field2
  end interface

  abstract interface    
    subroutine signature_compute_E_from_rho_1d(this, E, rho )
      use sll_working_precision
      import sll_maxwell_1d_base       
      class(sll_maxwell_1d_base) :: this
      sll_real64,dimension(:),intent(in) :: rho
      sll_real64,dimension(:),intent(out) :: E
    end subroutine signature_compute_E_from_rho_1d
  end interface


end module sll_m_maxwell_1d_base
