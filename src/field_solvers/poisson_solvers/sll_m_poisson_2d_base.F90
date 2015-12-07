!> @ingroup poisson_solvers
!> @brief
!> Module interface to solve Poisson equation in 2D
!> @details
!> Contains the abstract class to create a Poisson solver in 2D.
module sll_m_poisson_2d_base
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"

  implicit none

  public :: &
    sll_poisson_2d_base

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  !> PLEASE ADD DOCUMENTATION
  type, abstract :: sll_poisson_2d_base 

  contains

    !> solves \f$ -\Delta \phi_{ij} = \rho_{ij} \f$
    procedure(signature_compute_phi_from_rho_2d), deferred, pass(poisson) :: &
      compute_phi_from_rho

    !> solves \f$ -\Delta \phi_{ij} = \rho_{ij} \f$ and \f$ E_{ij} = \nabla  \phi_{ij} \f$
    procedure(signature_compute_E_from_rho_2d), deferred, pass(poisson) :: &
      compute_E_from_rho

  end type sll_poisson_2d_base

#ifndef DOXYGEN_SHOULD_SKIP_THIS

  abstract interface

    ! solves -\Delta phi = rho in 2d

    subroutine signature_compute_phi_from_rho_2d( poisson, phi, rho )

      use sll_m_working_precision
      import sll_poisson_2d_base      

      class(sll_poisson_2d_base), target     :: poisson
      sll_real64,dimension(:,:), intent(in)  :: rho
      sll_real64,dimension(:,:), intent(out) :: phi

    end subroutine signature_compute_phi_from_rho_2d

  end interface

  abstract interface    
    ! solves E = -\nabla Phi with -\Delta phi = rho in 2d 
    subroutine signature_compute_E_from_rho_2d( poisson, E1, E2, rho )

      use sll_m_working_precision
      import sll_poisson_2d_base       

      class(sll_poisson_2d_base)              :: poisson
      sll_real64, dimension(:,:), intent(in)  :: rho
      sll_real64, dimension(:,:), intent(out) :: E1
      sll_real64, dimension(:,:), intent(out) :: E2

    end subroutine signature_compute_E_from_rho_2d
  end interface

#endif /* DOXYGEN_SHOULD_SKIP_THIS */

end module sll_m_poisson_2d_base
