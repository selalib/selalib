!> @ingroup poisson_solvers
!> @brief
!> Module interface to solve Poisson equation in 3D
!> @details
!> Contains the abstract class to create a Poisson solver in 3D.
module sll_m_poisson_3d_base
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"

  implicit none

  public :: &
    sll_c_poisson_3d_base

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  !> Abstract class for Poisson solver in 3 dimensions
  type, abstract :: sll_c_poisson_3d_base 
  contains
    !> PLEASE ADD DOCUMENTATION
    procedure(signature_compute_phi_from_rho_3d), deferred, pass(poisson) :: &
      compute_phi_from_rho
!    procedure(signature_compute_E_from_phi_2d), deferred, pass(poisson) :: &
!      compute_E_from_phi  
    !> PLEASE ADD DOCUMENTATION
    procedure(signature_compute_E_from_rho_3d), deferred, pass(poisson) :: &
      compute_E_from_rho
  end type sll_c_poisson_3d_base

#ifndef DOXYGEN_SHOULD_SKIP_THIS

  abstract interface
    ! solves -\Delta phi = rho in 2d or similar thing
    subroutine signature_compute_phi_from_rho_3d( poisson, phi, rho )
      use sll_m_working_precision
      import sll_c_poisson_3d_base      
      class(sll_c_poisson_3d_base), target :: poisson
      sll_real64,dimension(:,:,:),intent(in) :: rho
      sll_real64,dimension(:,:,:),intent(out) :: phi
    end subroutine signature_compute_phi_from_rho_3d
  end interface

!  abstract interface  
!    ! solves E = -\nabla Phi in 2d
!    subroutine signature_compute_E_from_phi_2d( poisson, phi, E1, E2 )
!      use sll_m_working_precision
!      import sll_poisson_2d_base       
!      class(sll_poisson_2d_base) :: poisson
!      sll_real64,dimension(:,:),intent(in) :: phi
!      sll_real64,dimension(:,:),intent(out) :: E1
!      sll_real64,dimension(:,:),intent(out) :: E2
!    end subroutine signature_compute_E_from_phi_2d
!  end interface
  
  abstract interface    
    ! solves E = -\nabla Phi with -\Delta phi = rho in 2d 
    subroutine signature_compute_E_from_rho_3d( poisson, E1, E2, E3, rho )
      use sll_m_working_precision
      import sll_c_poisson_3d_base       
      class(sll_c_poisson_3d_base) :: poisson
      sll_real64,dimension(:,:,:),intent(in) :: rho
      sll_real64,dimension(:,:,:),intent(out) :: E1
      sll_real64,dimension(:,:,:),intent(out) :: E2
      sll_real64,dimension(:,:,:),intent(out) :: E3
    end subroutine signature_compute_E_from_rho_3d
  end interface

#endif /* DOXYGEN_SHOULD_SKIP_THIS */

end module sll_m_poisson_3d_base
