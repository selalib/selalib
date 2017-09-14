!> @brief  Define abstract base class for B-splines (uniform and non-uniform)
!> @author Yaman Güçlü  - IPP Garching
!> @author Edoardo Zoni - IPP Garching

module sll_m_bsplines_base
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"

  use sll_m_working_precision, only: f64

  implicit none

  public :: &
    sll_c_bsplines

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> Working precision
  integer, parameter :: wp = f64

  !> Abstract type, B-splines
  type, abstract :: sll_c_bsplines

    integer :: degree
    logical :: periodic
    logical :: uniform
    integer :: ncells
    integer :: nbasis
    integer :: offset

    real(wp) :: xmin
    real(wp) :: xmax

    real(wp), allocatable :: knots(:) ! Only used by non-uniform B-splines

  contains
    procedure(i_fun_find_cell              ), deferred :: find_cell
    procedure(i_sub_eval_basis             ), deferred :: eval_basis
    procedure(i_sub_eval_deriv             ), deferred :: eval_deriv
    procedure(i_sub_eval_basis_and_n_derivs), deferred :: eval_basis_and_n_derivs
    procedure(i_sub_free                   ), deferred :: free

  end type sll_c_bsplines

  abstract interface

    !---------------------------------------------------------------------------
    !> @brief     Find which grid cell contains the given point
    !> @param[in] self  B-splines object
    !> @param[in] x     point of interest
    !> results          cell index
    !---------------------------------------------------------------------------
    SLL_PURE function i_fun_find_cell( self, x ) result( icell )
     import sll_c_bsplines, wp
      class(sll_c_bsplines), intent(in) :: self
      real(wp)             , intent(in) :: x
      integer :: icell
    end function i_fun_find_cell

    !---------------------------------------------------------------------------
    !> Evaluate value at x of all basis functions with support in local cell
    !> values[j] = B_j(x) for jmin <= j <= jmin+degree
    !>
    !> @param[in]  self    B-splines object
    !> @param[in]  x       evaluation point
    !> @param[out] values  array of B-splines' values
    !> @param[out] jmin    index of first non-zero B-spline
    !---------------------------------------------------------------------------
    SLL_PURE subroutine i_sub_eval_basis( self, x, values, jmin )
     import sll_c_bsplines, wp
      class(sll_c_bsplines), intent(in   ) :: self
      real(wp)             , intent(in   ) :: x
      real(wp)             , intent(  out) :: values(:)
      integer              , intent(  out) :: jmin
    end subroutine i_sub_eval_basis

    !---------------------------------------------------------------------------
    !> Evaluate derivative at x of all basis functions with support in local cell
    !> derivs[j] = B_j'(x) for jmin <= j <= jmin+degree
    !>
    !> @param[in]  self    B-splines object
    !> @param[in]  x       evaluation point
    !> @param[out] derivs  array of B-splines' derivatives
    !> @param[out] jmin    index of first non-zero B-spline
    !---------------------------------------------------------------------------
    SLL_PURE subroutine i_sub_eval_deriv( self, x, derivs, jmin )
     import sll_c_bsplines, wp
      class(sll_c_bsplines), intent(in   ) :: self
      real(wp)             , intent(in   ) :: x
      real(wp)             , intent(  out) :: derivs(:)
      integer              , intent(  out) :: jmin
    end subroutine i_sub_eval_deriv

    !---------------------------------------------------------------------------
    !> Evaluate value and n derivatives at x of all basis functions with support in local cell
    !> derivs[i,j] = (d/dx)^i B_j(x) for 0 <= i <= n and jmin <= j <= jmin+degree
    !>
    !> @param[in]  self    B-splines object
    !> @param[in]  x       evaluation point
    !> @param[in]  n       number of required derivatives
    !> @param[out] derivs  array of B-splines' (multiple) derivatives
    !> @param[out] jmin    index of first non-zero B-spline
    !---------------------------------------------------------------------------
    SLL_PURE subroutine i_sub_eval_basis_and_n_derivs( self, x, n, derivs, jmin )
     import sll_c_bsplines, wp
      class(sll_c_bsplines), intent(in   ) :: self
      real(wp)             , intent(in   ) :: x
      integer              , intent(in   ) :: n
      real(wp)             , intent(  out) :: derivs(:,:)
      integer              , intent(  out) :: jmin
    end subroutine i_sub_eval_basis_and_n_derivs

    !---------------------------------------------------------------------------
    !> @brief        Free storage
    !> @param[inout] self  B-splines object
    !---------------------------------------------------------------------------
    subroutine i_sub_free( self )
     import sll_c_bsplines
      class(sll_c_bsplines), intent(inout) :: self
    end subroutine i_sub_free

  end interface

end module sll_m_bsplines_base
