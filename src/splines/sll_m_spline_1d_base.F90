!> @ingroup splines 
!> @brief   Abstract base class for 1D splines (uniform and non-uniform)
!> @author  Yaman Güçlü, IPP Garching

module sll_m_spline_1d_base
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"

  use sll_m_working_precision, only: &
    f64

  implicit none

  public :: &
    sll_c_spline_1d

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> Working precision
  integer, parameter :: wp = f64

  !> Abstract type for 1D splines
  type, abstract :: sll_c_spline_1d
  contains
    procedure(i_sub_free                ), deferred :: free
    procedure(i_sub_compute_interpolant ), deferred :: compute_interpolant
    procedure(i_fun_evaluate_at_x       ), deferred :: eval
    procedure(i_fun_evaluate_at_x       ), deferred :: eval_deriv
    procedure(i_sub_evaluate_at_x_array ), deferred :: eval_array       ! could be implemented here
    procedure(i_sub_evaluate_at_x_array ), deferred :: eval_array_deriv ! could be implemented here
    procedure(i_fun_get_pointer_to_coeff), deferred :: get_coeff
  end type sll_c_spline_1d


  abstract interface

  !-----------------------------------------------------------------------------
  subroutine i_sub_free( self )
   import sll_c_spline_1d
    class(sll_c_spline_1d), intent(inout) :: self
  end subroutine i_sub_free

  !-----------------------------------------------------------------------------
  subroutine i_sub_compute_interpolant( self, gtau, derivs_xmin, derivs_xmax )
   import sll_c_spline_1d, wp
    class(sll_c_spline_1d), intent(inout) :: self
    real(wp)              , intent(in   ) :: gtau(:)
    real(wp), optional    , intent(in   ) :: derivs_xmin(:)
    real(wp), optional    , intent(in   ) :: derivs_xmax(:)
  end subroutine i_sub_compute_interpolant

  !-----------------------------------------------------------------------------
  SLL_PURE function i_fun_evaluate_at_x( self, x ) result( y )
   import sll_c_spline_1d, wp
    class(sll_c_spline_1d), intent(in) :: self
    real(wp)              , intent(in) :: x
    real(wp) :: y
  end function i_fun_evaluate_at_x

  !-----------------------------------------------------------------------------
  SLL_PURE subroutine i_sub_evaluate_at_x_array( self, x, y )
   import sll_c_spline_1d, wp
    class(sll_c_spline_1d), intent(in   ) :: self
    real(wp)              , intent(in   ) :: x(:)
    real(wp)              , intent(  out) :: y(size(x))
  end subroutine i_sub_evaluate_at_x_array

  !-----------------------------------------------------------------------------
  function i_fun_get_pointer_to_coeff( self ) result( ptr )
   import sll_c_spline_1d, wp
    class(sll_c_spline_1d), target, intent(in) :: self
    real(wp), pointer :: ptr(:)
  end function i_fun_get_pointer_to_coeff

  !-----------------------------------------------------------------------------
  end interface

end module sll_m_spline_1d_base
