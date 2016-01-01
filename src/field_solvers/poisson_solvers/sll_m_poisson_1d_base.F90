!> @ingroup poisson_solvers
!> @brief
!> Module interface to solve Poisson equation in 1D.
!> @details
!> This decribes two methods
!>   - compute_phi_from_rho
!>   - compute_e_from_rho
module sll_m_poisson_1d_base
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"

  implicit none

  public :: &
    sll_c_poisson_1d_base

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  type, abstract :: sll_c_poisson_1d_base 
  contains
    !> Solve Poisson equation and compute electric potential from charge density, 
    !> \f[ -\Delta \phi_i = \rho_i \f]
    procedure(signature_compute_phi_from_rho_1d), deferred, pass(poisson) :: &
      compute_phi_from_rho
    !> Solve Poisson equation and compute electric field from charge density
    !> \f[ E_i = -\nabla \phi_i \f]  with \f[ -\Delta \phi_i = \rho_i \f].
    procedure(signature_compute_E_from_rho_1d), deferred, pass(poisson) :: &
      compute_E_from_rho
  end type sll_c_poisson_1d_base


#ifndef DOXYGEN_SHOULD_SKIP_THIS

  abstract interface
    subroutine signature_compute_phi_from_rho_1d( poisson, phi, rho )
      use sll_m_working_precision
      import sll_c_poisson_1d_base      
      class(sll_c_poisson_1d_base), target :: poisson
      sll_real64,dimension(:),intent(in) :: rho
      sll_real64,dimension(:),intent(out) :: phi
    end subroutine signature_compute_phi_from_rho_1d
  end interface

  abstract interface    
    subroutine signature_compute_E_from_rho_1d( poisson, E, rho )
      use sll_m_working_precision
      import sll_c_poisson_1d_base       
      class(sll_c_poisson_1d_base) :: poisson
      sll_real64,dimension(:),intent(in) :: rho
      sll_real64,dimension(:),intent(out) :: E
    end subroutine signature_compute_E_from_rho_1d
  end interface

#endif /* DOXYGEN_SHOULD_SKIP_THIS */

end module sll_m_poisson_1d_base

