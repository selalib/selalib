module sll_m_qn_2d_base
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"

  implicit none

  public :: &
    sll_qn_2d_base

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type, abstract :: sll_qn_2d_base 
  contains
    procedure(signature_precompute_qn_2d), deferred, pass(qn) :: &
      precompute_qn
    procedure(signature_solve_qn_2d), deferred, pass(qn) :: &
      solve_qn  
  end type sll_qn_2d_base

  abstract interface

    subroutine signature_precompute_qn_2d( qn, mu_points, mu_weights , N_mu)
      use sll_m_working_precision
      import sll_qn_2d_base      
      class(sll_qn_2d_base), target :: qn
      sll_int32,intent(in) :: N_mu
      sll_real64,dimension(1:N_mu),intent(in) :: mu_points
      sll_real64,dimension(1:N_mu),intent(in) :: mu_weights
    end subroutine signature_precompute_qn_2d
  
    subroutine signature_solve_qn_2d( qn, phi)
      use sll_m_working_precision
      import sll_qn_2d_base      
      class(sll_qn_2d_base), target :: qn
      sll_real64,dimension(:,:),intent(inout) :: phi
    end subroutine signature_solve_qn_2d
  end interface

end module sll_m_qn_2d_base