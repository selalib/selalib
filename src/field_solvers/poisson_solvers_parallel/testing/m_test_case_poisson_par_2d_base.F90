module m_test_case_poisson_par_2d_base
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"

implicit none

public :: &
  c_test_case_poisson_2d_polar

private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type, abstract :: c_test_case_poisson_2d_polar

    sll_real64 :: rmin
    sll_real64 :: rmax
    sll_int32  :: bc_rmin
    sll_int32  :: bc_rmax

  contains
    ! 2D manufactured solution
    procedure( i_func_2d_real ), deferred :: phi_ex
    procedure( i_func_2d_real ), deferred :: phi_ex_diff1_r
    procedure( i_func_2d_real ), deferred :: phi_ex_diff2_r
    procedure( i_func_2d_real ), deferred :: phi_ex_diff2_th
    ! 2D right-hand side to solver
    procedure, non_overridable :: rhs => f_test_case__rhs

  end type c_test_case_poisson_2d_polar

  !-----------------------------------------------------------------------------
  abstract interface

    ! 2D polar profile, scalar real function
    pure function i_func_2d_real( self, r, th ) result( val )
      use sll_m_working_precision
      import c_test_case_poisson_2d_polar
      class( c_test_case_poisson_2d_polar ), intent(in) :: self
      sll_real64, intent(in) :: r
      sll_real64, intent(in) :: th
      sll_real64 :: val
    end function i_func_2d_real

  end interface

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
contains
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  pure function f_test_case__rhs( self, r, th ) result( rho )
    class(c_test_case_poisson_2d_polar), intent(in) :: self
    sll_real64, intent(in) :: r
    sll_real64, intent(in) :: th
    sll_real64 :: rho

    rho = - self%phi_ex_diff2_r( r, th ) - self%phi_ex_diff1_r( r, th )/r &
          - self%phi_ex_diff2_th( r, th )/r**2

  end function f_test_case__rhs

  
end module m_test_case_poisson_par_2d_base
