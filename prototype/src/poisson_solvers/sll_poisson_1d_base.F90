module sll_module_poisson_1d_base
#include "sll_working_precision.h"
  implicit none
  
  ! For computing E and phi form rho, using Poisson equation
  type, abstract :: sll_poisson_1d_base 
  contains
    procedure(signature_compute_phi_from_rho_1d), deferred, pass(poisson) :: &
      compute_phi_from_rho
!    procedure(signature_compute_E_from_phi_2d), deferred, pass(poisson) :: &
!      compute_E_from_phi  
    procedure(signature_compute_E_from_rho_1d), deferred, pass(poisson) :: &
      compute_E_from_rho
  end type sll_poisson_1d_base

  abstract interface
    ! solves -\Delta phi = rho in 1d
    subroutine signature_compute_phi_from_rho_1d( poisson, phi, rho )
      use sll_working_precision
      import sll_poisson_1d_base      
      class(sll_poisson_1d_base), target :: poisson
      sll_real64,dimension(:),intent(in) :: rho
      sll_real64,dimension(:),intent(out) :: phi
    end subroutine signature_compute_phi_from_rho_1d
  end interface

!  abstract interface  
!    ! solves E = -\nabla Phi in 2d
!    subroutine signature_compute_E_from_phi_2d( poisson, phi, E1, E2 )
!      use sll_working_precision
!      import sll_poisson_2d_base       
!      class(sll_poisson_2d_base) :: poisson
!      sll_real64,dimension(:,:),intent(in) :: phi
!      sll_real64,dimension(:,:),intent(out) :: E1
!      sll_real64,dimension(:,:),intent(out) :: E2
!    end subroutine signature_compute_E_from_phi_2d
!  end interface
  
  abstract interface    
    ! solves E = -\nabla Phi with -\Delta phi = rho in 1d 
    subroutine signature_compute_E_from_rho_1d( poisson, E, rho )
      use sll_working_precision
      import sll_poisson_1d_base       
      class(sll_poisson_1d_base) :: poisson
      sll_real64,dimension(:),intent(in) :: rho
      sll_real64,dimension(:),intent(out) :: E
    end subroutine signature_compute_E_from_rho_1d
  end interface

end module sll_module_poisson_1d_base

