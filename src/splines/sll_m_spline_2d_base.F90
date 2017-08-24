!> @ingroup splines 
!> @brief   Abstract base class for 1D splines (uniform and non-uniform)
!> @author  Yaman Güçlü, IPP Garching

module sll_m_spline_2d_base
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"

  use sll_m_working_precision, only: &
    f64

  implicit none

  public :: &
    sll_c_spline_2d, &
    sll_t_spline_2d_boundary_data

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> Working precision
  integer, parameter :: wp = f64

  !> Abstract type for 1D splines
  type, abstract :: sll_c_spline_2d
  contains
    procedure(i_sub_free                ), deferred :: free
    procedure(i_sub_compute_interpolant ), deferred :: compute_interpolant
    procedure(i_fun_evaluate_at_x       ), deferred :: eval
    procedure(i_fun_evaluate_at_x       ), deferred :: eval_deriv_x1
    procedure(i_fun_evaluate_at_x       ), deferred :: eval_deriv_x2
    procedure(i_sub_evaluate_at_x_array ), deferred :: eval_array          ! could be implemented here
    procedure(i_sub_evaluate_at_x_array ), deferred :: eval_array_deriv_x1 ! could be implemented here
    procedure(i_sub_evaluate_at_x_array ), deferred :: eval_array_deriv_x2 ! could be implemented here
!    procedure(i_sub_get_interp_points   ), deferred :: get_interp_points
  end type sll_c_spline_2d

  !> Container for boundary condition data
  !>
  !>  x2_max  ____________
  !>         |            |
  !>         | c        d |
  !>         |            |
  !>         |            |
  !>         |            |
  !>         | a        b |
  !>  x2_min |____________|
  !>       x1_min       x1_max
  !>
  type :: sll_t_spline_2d_boundary_data
    real(wp), allocatable :: derivs_x1_min (:,:)
    real(wp), allocatable :: derivs_x1_max (:,:)
    real(wp), allocatable :: derivs_x2_min (:,:)
    real(wp), allocatable :: derivs_x2_max (:,:)
    real(wp), allocatable :: mixed_derivs_a(:,:)
    real(wp), allocatable :: mixed_derivs_b(:,:)
    real(wp), allocatable :: mixed_derivs_c(:,:)
    real(wp), allocatable :: mixed_derivs_d(:,:)
  end type sll_t_spline_2d_boundary_data


  abstract interface

  !-----------------------------------------------------------------------------
  subroutine i_sub_free( self )
   import sll_c_spline_2d
    class(sll_c_spline_2d), intent(inout) :: self
  end subroutine i_sub_free

  !-----------------------------------------------------------------------------
  subroutine i_sub_compute_interpolant( self, gtau, boundary_data )
   import sll_c_spline_2d, wp, sll_t_spline_2d_boundary_data
    class(sll_c_spline_2d)             , intent(inout)          :: self
    real(wp)                           , intent(in   )          :: gtau(:,:)
    type(sll_t_spline_2d_boundary_data), intent(in   ), optional:: boundary_data
  end subroutine i_sub_compute_interpolant

  !-----------------------------------------------------------------------------
  SLL_PURE function i_fun_evaluate_at_x( self, x1, x2 ) result( y )
   import sll_c_spline_2d, wp
    class(sll_c_spline_2d), intent(in) :: self
    real(wp)              , intent(in) :: x1
    real(wp)              , intent(in) :: x2
    real(wp) :: y
  end function i_fun_evaluate_at_x

  !-----------------------------------------------------------------------------
  SLL_PURE subroutine i_sub_evaluate_at_x_array( self, x1, x2, y )
   import sll_c_spline_2d, wp
    class(sll_c_spline_2d), intent(in   ) :: self
    real(wp)              , intent(in   ) :: x1(:,:)
    real(wp)              , intent(in   ) :: x2(:,:)
    real(wp)              , intent(  out) :: y (:,:)
!    real(wp)              , intent(in   ) :: x2(size(x1,1),size(x1,2))
!    real(wp)              , intent(  out) :: y (size(x1,1),size(x1,2))
  end subroutine i_sub_evaluate_at_x_array

  !-----------------------------------------------------------------------------
  subroutine i_sub_get_interp_points( self, tau1, tau2 )
   import sll_c_spline_2d, wp
    class(sll_c_spline_2d), intent(in   ) :: self
    real(wp),  allocatable, intent(  out) :: tau1(:)
    real(wp),  allocatable, intent(  out) :: tau2(:)
  end subroutine i_sub_get_interp_points

  !-----------------------------------------------------------------------------
  end interface

end module sll_m_spline_2d_base
