!> @ingroup poisson_solvers
!> @brief
!> Parent class for Poisson solvers in 1 dimension
!> @details
!> This decribes two methods
!>   - compute_phi_from_rho
!>   - compute_e_from_rho
module sll_module_poisson_1d_base
#include "sll_working_precision.h"
  implicit none
  
  type, abstract :: sll_poisson_1d_base 
  contains
    !> Solve Poisson equation and compute electric potential from charge density, 
    procedure(signature_compute_phi_from_rho_1d), deferred, pass(poisson) :: &
      compute_phi_from_rho
    !> Solve Poisson equation and compute electric field from charge density
    procedure(signature_compute_E_from_rho_1d), deferred, pass(poisson) :: &
      compute_E_from_rho
  end type sll_poisson_1d_base


#ifndef DOXYGEN_SHOULD_SKIP_THIS

  abstract interface
    !> solves \f[ -\Delta phi = rho \f] in 1d 
    subroutine signature_compute_phi_from_rho_1d( poisson, phi, rho )
      use sll_working_precision
      import sll_poisson_1d_base      
      class(sll_poisson_1d_base), target :: poisson
      sll_real64,dimension(:),intent(in) :: rho
      sll_real64,dimension(:),intent(out) :: phi
    end subroutine signature_compute_phi_from_rho_1d
  end interface

  abstract interface    
    !> solves E = -\nabla Phi with -\Delta phi = rho in 1d 
    subroutine signature_compute_E_from_rho_1d( poisson, E, rho )
      use sll_working_precision
      import sll_poisson_1d_base       
      class(sll_poisson_1d_base) :: poisson
      sll_real64,dimension(:),intent(in) :: rho
      sll_real64,dimension(:),intent(out) :: E
    end subroutine signature_compute_E_from_rho_1d
  end interface

#endif /* DOXYGEN_SHOULD_SKIP_THIS */

end module sll_module_poisson_1d_base

