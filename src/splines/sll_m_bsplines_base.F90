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

  contains
    procedure(i_sub_eval_basis             ), deferred :: eval_basis
    procedure(i_sub_eval_deriv             ), deferred :: eval_deriv
    procedure(i_sub_eval_basis_and_n_derivs), deferred :: eval_basis_and_n_derivs
    procedure(i_sub_free                   ), deferred :: free

  end type sll_c_bsplines

  abstract interface

    !> Evaluate value at x of all basis functions with support in local cell
    !> (jmin identifies index of basis functions)
    SLL_PURE subroutine i_sub_eval_basis( self, x, values, jmin )
     import sll_c_bsplines, wp
      class(sll_c_bsplines), intent(in   ) :: self
      real(wp)             , intent(in   ) :: x
      real(wp)             , intent(  out) :: values(:)
      integer              , intent(  out) :: jmin
    end subroutine i_sub_eval_basis

    !> Evaluate derivative at x of all basis functions with support in local cell
    SLL_PURE subroutine i_sub_eval_deriv( self, x, derivs, jmin )
     import sll_c_bsplines, wp
      class(sll_c_bsplines), intent(in   ) :: self
      real(wp)             , intent(in   ) :: x
      real(wp)             , intent(  out) :: derivs(:)
      integer              , intent(  out) :: jmin
    end subroutine i_sub_eval_deriv

    !> Evaluate value and n derivatives at x of all basis functions with support in local cell
    SLL_PURE subroutine i_sub_eval_basis_and_n_derivs( self, x, n, derivs, jmin )
     import sll_c_bsplines, wp
      class(sll_c_bsplines), intent(in   ) :: self
      real(wp)             , intent(in   ) :: x
      integer              , intent(in   ) :: n
      real(wp)             , intent(  out) :: derivs(:,:)
      integer              , intent(  out) :: jmin
    end subroutine i_sub_eval_basis_and_n_derivs

    !> Free storage
    subroutine i_sub_free( self )
     import sll_c_bsplines
      class(sll_c_bsplines), intent(inout) :: self
    end subroutine i_sub_free

  end interface

end module sll_m_bsplines_base
