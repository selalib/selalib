!> @authors Yaman Güçlü, IPP Garching
!> @authors Edoardo Zoni, IPP Garching
!>
!> @brief
!> Method of manufactured solutions for Poisson's equation in 2D polar coordinates.
!>
!> @details
!> This module defines an abstract interface which requires subclasses to implement
!> an analytical function phi(r,theta) and its derivatives.
!> The analytical rho(r,theta) is calculated here according to the Poisson
!> equation, and it is not overridable.
!> In the numerical tests rho(r,theta) will be given to the Poisson solver, and
!> the resulting numerical phi will be compared to the exact solution.
!>
module m_test_poisson_2d_polar_base
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"

implicit none

public :: &
  c_test_poisson_2d_polar_base

private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type, abstract :: c_test_poisson_2d_polar_base

  contains
    ! Get domain limits and boundary conditions
    procedure( i_func_get_rlim), deferred :: get_rlim
    procedure( i_func_get_bcs ), deferred :: get_bcs
    ! 2D manufactured solution
    procedure( i_func_2d_real ), deferred :: phi_ex
    procedure( i_func_2d_real ), deferred :: phi_ex_diff1_r
    procedure( i_func_2d_real ), deferred :: phi_ex_diff2_r
    procedure( i_func_2d_real ), deferred :: phi_ex_diff2_th
    ! 2D right-hand side to solver
    procedure, non_overridable :: rhs => f_test_case__rhs

  end type c_test_poisson_2d_polar_base

  !-----------------------------------------------------------------------------
  abstract interface

    ! Get domain limits
    pure function i_func_get_rlim( self ) result( rlim )
      use sll_m_working_precision
      import c_test_poisson_2d_polar_base
      class( c_test_poisson_2d_polar_base ), intent(in) :: self
      sll_real64 :: rlim(2)
    end function i_func_get_rlim

    ! Get boundary conditions
    pure function i_func_get_bcs( self ) result( bcs )
      use sll_m_working_precision
      import c_test_poisson_2d_polar_base
      class( c_test_poisson_2d_polar_base ), intent(in) :: self
      sll_int32 :: bcs(2)
    end function i_func_get_bcs

    ! 2D polar profile, scalar real function
    pure function i_func_2d_real( self, r, th ) result( val )
      use sll_m_working_precision
      import c_test_poisson_2d_polar_base
      class( c_test_poisson_2d_polar_base ), intent(in) :: self
      sll_real64, intent(in) :: r
      sll_real64, intent(in) :: th
      sll_real64 :: val
    end function i_func_2d_real

  end interface

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
contains
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  pure function f_test_case__rhs( self, r, th ) result( rho )
    class(c_test_poisson_2d_polar_base), intent(in) :: self
    sll_real64                         , intent(in) :: r
    sll_real64                         , intent(in) :: th
    sll_real64 :: rho

    rho = - self%phi_ex_diff2_r( r, th ) - self%phi_ex_diff1_r( r, th )/r &
          - self%phi_ex_diff2_th( r, th )/r**2

  end function f_test_case__rhs

  
end module m_test_poisson_2d_polar_base
