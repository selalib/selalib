module sll_m_linear_operator_matrix_stencil_to_dense_new
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_errors.h"

  use sll_m_working_precision, only: f64

  use sll_m_linear_operator_base, only: sll_c_linear_operator

  use sll_m_vector_space_base, only: sll_c_vector_space

  use sll_m_vector_space_real_array_1d, only: sll_t_vector_space_real_array_1d

  use sll_m_vector_space_real_array_2d, only: sll_t_vector_space_real_array_2d

  implicit none

  public :: sll_t_linear_operator_matrix_stencil_to_dense_new

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! Working precision
  integer, parameter :: wp = f64

  type, extends(sll_c_linear_operator) :: sll_t_linear_operator_matrix_stencil_to_dense_new

    real(wp), allocatable :: A (:,:)

    integer :: n1
    integer :: n2

  contains

    procedure :: init      => s_linear_operator_matrix_stencil_to_dense_new__init
    procedure :: get_shape => f_linear_operator_matrix_stencil_to_dense_new__get_shape
    procedure :: dot       => s_linear_operator_matrix_stencil_to_dense_new__dot
    procedure :: dot_incr  => s_linear_operator_matrix_stencil_to_dense_new__dot_incr
    procedure :: to_array  => s_linear_operator_matrix_stencil_to_dense_new__to_array
    procedure :: free      => s_linear_operator_matrix_stencil_to_dense_new__free

  end type sll_t_linear_operator_matrix_stencil_to_dense_new

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! Initialize linear operator
  subroutine s_linear_operator_matrix_stencil_to_dense_new__init( self, n1, n2 )
    class(sll_t_linear_operator_matrix_stencil_to_dense_new), intent(inout) :: self
    integer                                                 , intent(in   ) :: n1
    integer                                                 , intent(in   ) :: n2

    allocate( self % A( n1, n2 ) )

    self % n1 = n1
    self % n2 = n2

  end subroutine s_linear_operator_matrix_stencil_to_dense_new__init

  ! Get shape of linear operator
  function f_linear_operator_matrix_stencil_to_dense_new__get_shape( self ) result( s )
    class(sll_t_linear_operator_matrix_stencil_to_dense_new), intent(in) :: self
    integer :: s(2)

    s = shape( self % A  )

  end function f_linear_operator_matrix_stencil_to_dense_new__get_shape

  ! Implement y=Ax, with A dense matrix (2D), x stencil vector (2D) and y dense vector (1D)
  subroutine s_linear_operator_matrix_stencil_to_dense_new__dot( self, x, y )
    class(sll_t_linear_operator_matrix_stencil_to_dense_new), intent(in   ) :: self
    class(sll_c_vector_space)                               , intent(in   ) :: x ! stencil
    class(sll_c_vector_space)                               , intent(inout) :: y ! dense

    integer :: j1, j2, i, j

    character(len=*), parameter :: this_sub_name = "sll_t_linear_operator_matrix_stencil_to_dense_new % dot"
    character(len=64) :: err_msg

    select type ( x )

    type is ( sll_t_vector_space_real_array_2d )

      associate( p1 => -lbound( x % array, 1 ) + 1, &
                 p2 => -lbound( x % array, 2 ) + 1 )

      associate( nx1 => ubound( x % array, 1 ) - p1, &
                 nx2 => ubound( x % array, 2 ) - p2 )

        ! Check dimensions
        SLL_ASSERT( self % n2 == nx1*nx2 )

        select type ( y )

        type is ( sll_t_vector_space_real_array_1d )

          ! Check dimensions
          SLL_ASSERT( self  % n1 == size( y % array ) )

          y % array = 0.0_wp
          do i = 1, size( y % array )
            do j2 = 1, nx2
              do j1 = 1, nx1
                j = (j1-1) * nx2 + j2
                y % array(i) = y % array(i) + self % A(i,j) * x % array(j1,j2)
              end do
            end do
          end do

        class default
          err_msg = "y must be of type sll_t_vector_space_real_array_1d"
          SLL_ERROR( this_sub_name, err_msg )

        end select

      end associate

      end associate

    class default
      err_msg = "x must be of type sll_t_vector_space_real_array_2d"
      SLL_ERROR( this_sub_name, err_msg )

    end select

  end subroutine s_linear_operator_matrix_stencil_to_dense_new__dot

  ! Implement y=y+Ax, with A dense matrix (2D), x stencil vector (2D) and y dense vector (1D)
  subroutine s_linear_operator_matrix_stencil_to_dense_new__dot_incr( self, x, y )
    class(sll_t_linear_operator_matrix_stencil_to_dense_new), intent(in   ) :: self
    class(sll_c_vector_space)                               , intent(in   ) :: x ! stencil
    class(sll_c_vector_space)                               , intent(inout) :: y ! dense

    integer :: j1, j2, i, j

    character(len=*), parameter :: this_sub_name = "sll_t_linear_operator_matrix_stencil_to_dense_new % dot"
    character(len=64) :: err_msg

    select type ( x )

    type is ( sll_t_vector_space_real_array_2d )

      associate( p1 => -lbound( x % array, 1 ) + 1, &
                 p2 => -lbound( x % array, 2 ) + 1 )

      associate( nx1 => ubound( x % array, 1 ) - p1, &
                 nx2 => ubound( x % array, 2 ) - p2 )

        ! Check dimensions
        SLL_ASSERT( self % n2 == nx1*nx2 )

        select type ( y )

        type is ( sll_t_vector_space_real_array_1d )

          ! Check dimensions
          SLL_ASSERT( self  % n1 == size( y % array ) )

          do i = 1, size( y % array )
            do j2 = 1, nx2
              do j1 = 1, nx1
                j = (j1-1) * nx2 + j2
                y % array(i) = y % array(i) + self % A(i,j) * x % array(j1,j2)
              end do
            end do
          end do

        class default
          err_msg = "y must be of type sll_t_vector_space_real_array_1d"
          SLL_ERROR( this_sub_name, err_msg )

        end select

      end associate

      end associate

    class default
      err_msg = "x must be of type sll_t_vector_space_real_array_2d"
      SLL_ERROR( this_sub_name, err_msg )

    end select

  end subroutine s_linear_operator_matrix_stencil_to_dense_new__dot_incr

  ! Convert dense matrix to array (trivial)
  subroutine s_linear_operator_matrix_stencil_to_dense_new__to_array( self, A )
    class(sll_t_linear_operator_matrix_stencil_to_dense_new), intent(in   ) :: self
    real(wp)                                                , intent(inout) :: A(:,:)

    SLL_ASSERT( size( A, 1 ) == self % n1 )
    SLL_ASSERT( size( A, 2 ) == self % n2 )

    A = self % A

  end subroutine s_linear_operator_matrix_stencil_to_dense_new__to_array

  ! Free objects
  subroutine s_linear_operator_matrix_stencil_to_dense_new__free( self )
    class(sll_t_linear_operator_matrix_stencil_to_dense_new), intent(inout) :: self

    deallocate( self % A )

  end subroutine s_linear_operator_matrix_stencil_to_dense_new__free

end module sll_m_linear_operator_matrix_stencil_to_dense_new
