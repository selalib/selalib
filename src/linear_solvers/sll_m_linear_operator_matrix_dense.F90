module sll_m_linear_operator_matrix_dense
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_errors.h"

  use sll_m_working_precision, only: f64

  use sll_m_linear_operator_base, only: sll_c_linear_operator

  use sll_m_vector_space_base, only: sll_c_vector_space

  use sll_m_vector_space_real_arrays, only: &
    sll_t_vector_space_real_1d, &
    sll_t_vector_space_real_2d

  implicit none

  public :: sll_t_linear_operator_matrix_dense

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! Working precision
  integer, parameter :: wp = f64

  type, extends(sll_c_linear_operator) :: sll_t_linear_operator_matrix_dense

    ! Pointer to matrix A
    type(sll_t_vector_space_real_2d), pointer :: A => null()

  contains

    procedure :: init      => s_linear_operator_matrix_dense__init
    procedure :: get_shape => f_linear_operator_matrix_dense__get_shape
    procedure :: dot       => s_linear_operator_matrix_dense__dot
    procedure :: free      => s_linear_operator_matrix_dense__free

  end type sll_t_linear_operator_matrix_dense

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! Initialize linear operator
  subroutine s_linear_operator_matrix_dense__init( self, A )
    class(sll_t_linear_operator_matrix_dense), intent(inout) :: self
    type(sll_t_vector_space_real_2d), target , intent(in   ) :: A

    self % A => A

  end subroutine s_linear_operator_matrix_dense__init

  ! Get shape of linear operator
  function f_linear_operator_matrix_dense__get_shape( self ) result( s )
    class(sll_t_linear_operator_matrix_dense), intent(in   ) :: self
    integer :: s(2)

    s = shape( self % A % array )

  end function f_linear_operator_matrix_dense__get_shape

  ! Implement Ax=y, with A dense matrix, x and y real 1D arrays
  subroutine s_linear_operator_matrix_dense__dot( self, x, y )
    class(sll_t_linear_operator_matrix_dense), intent(in   ) :: self
    class(sll_c_vector_space)                , intent(in   ) :: x
    class(sll_c_vector_space)                , intent(inout) :: y ! already constructed

    integer :: i, j, nx(1), ny(1), n(2)

    character(len=*), parameter :: this_sub_name = "sll_t_linear_operator_matrix_dense % dot"
    character(len=64) :: err_msg

    n = self % get_shape()

    ! Make sure to work with 1D real arrays
    select type ( x )

    type is ( sll_t_vector_space_real_1d )

      ! Check if A and x are compatible for multiplication
      nx = shape( x % array )
      SLL_ASSERT( n(2) == nx(1) )

      select type ( y )

      type is ( sll_t_vector_space_real_1d )

        ! Check if y and Ax are compatible for multiplication
        ny = shape( y % array )
        SLL_ASSERT( n(1) == ny(1) )

        y % array(:) = 0.0_wp
        do i = 1, n(1)
          do j = 1, n(2)
            y % array(i) = y % array(i) + self % A % array(i,j) * x % array(j)
          end do
        end do

      class default
        err_msg = "y must be of type sll_t_vector_space_real_1d"
        SLL_ERROR( this_sub_name, err_msg )

      end select

    class default
      err_msg = "x must be of type sll_t_vector_space_real_1d"
      SLL_ERROR( this_sub_name, err_msg )

    end select

  end subroutine s_linear_operator_matrix_dense__dot

  ! Free objects
  subroutine s_linear_operator_matrix_dense__free( self )
    class(sll_t_linear_operator_matrix_dense), intent(inout) :: self

    self % A => null()

  end subroutine s_linear_operator_matrix_dense__free

end module sll_m_linear_operator_matrix_dense
