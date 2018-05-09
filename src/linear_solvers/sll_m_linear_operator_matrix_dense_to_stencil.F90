module sll_m_linear_operator_matrix_dense_to_stencil
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_errors.h"

  use sll_m_working_precision, only: f64

  use sll_m_linear_operator_base, only: sll_c_linear_operator

  use sll_m_vector_space_base, only: sll_c_vector_space

  use sll_m_vector_space_real_array_1d, only: sll_t_vector_space_real_array_1d

  use sll_m_vector_space_real_array_2d, only: sll_t_vector_space_real_array_2d

  implicit none

  public :: sll_t_linear_operator_matrix_dense_to_stencil

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! Working precision
  integer, parameter :: wp = f64

  type, extends(sll_c_linear_operator) :: sll_t_linear_operator_matrix_dense_to_stencil

    real(wp), allocatable :: A (:,:)

    integer :: s1
    integer :: s2

  contains

    procedure :: init      => s_linear_operator_matrix_dense_to_stencil__init
    procedure :: get_shape => f_linear_operator_matrix_dense_to_stencil__get_shape
    procedure :: dot       => s_linear_operator_matrix_dense_to_stencil__dot
    procedure :: dot_incr  => s_linear_operator_matrix_dense_to_stencil__dot_incr
    procedure :: to_array  => s_linear_operator_matrix_dense_to_stencil__to_array
    procedure :: free      => s_linear_operator_matrix_dense_to_stencil__free

  end type sll_t_linear_operator_matrix_dense_to_stencil

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! Initialize linear operator
  subroutine s_linear_operator_matrix_dense_to_stencil__init( self, s1, s2 )
    class(sll_t_linear_operator_matrix_dense_to_stencil), intent(inout) :: self
    integer                                             , intent(in   ) :: s1
    integer                                             , intent(in   ) :: s2

    allocate( self % A( s1, s2 ) )

    self % s1 = s1
    self % s2 = s2

  end subroutine s_linear_operator_matrix_dense_to_stencil__init

  ! Get shape of linear operator
  function f_linear_operator_matrix_dense_to_stencil__get_shape( self ) result( s )
    class(sll_t_linear_operator_matrix_dense_to_stencil), intent(in) :: self
    integer :: s(2)

    s = shape( self % A  )

  end function f_linear_operator_matrix_dense_to_stencil__get_shape

  ! Implement y=Ax, with A dense matrix (2D), x dense vector (1D) and y stencil vector (2D)
  subroutine s_linear_operator_matrix_dense_to_stencil__dot( self, x, y )
    class(sll_t_linear_operator_matrix_dense_to_stencil), intent(in   ) :: self
    class(sll_c_vector_space)                           , intent(in   ) :: x ! dense
    class(sll_c_vector_space)                           , intent(inout) :: y ! stencil

    integer :: j1, j2, j, k

    character(len=*), parameter :: this_sub_name = "sll_t_linear_operator_matrix_dense_to_stencil % dot"
    character(len=64) :: err_msg

    select type ( x )

    type is ( sll_t_vector_space_real_array_1d )

      ! Check dimensions
      SLL_ASSERT( self % s2 == size( x % array ) )

      select type ( y )

      type is ( sll_t_vector_space_real_array_2d )

        associate( p1 => -lbound( y % array, 1 ) + 1, &
                   p2 => -lbound( y % array, 2 ) + 1 )

        associate( ny1 => ubound( y % array, 1 ) - p1, &
                   ny2 => ubound( y % array, 2 ) - p2 )

          ! Check dimensions
          SLL_ASSERT( self % s1 == ny1*ny2 )

          y % array = 0.0_wp
          do j2 = 1, ny2
            do j1 = 1, ny1
              j = (j1-1) * ny2 + j2
              do k = 1, size( x % array )
                y % array(j1,j2) = y % array(j1,j2) + self % A(j,k) * x % array(k)
              end do
            end do
          end do

          ! Update buffer regions
          y % array(1-p1:0      ,:) = y % array(ny1-p1+1:ny1,:)
          y % array(ny1+1:ny1+p1,:) = y % array(1:p1        ,:)
          y % array(:,1-p2:0      ) = y % array(:,ny2-p2+1:ny2)
          y % array(:,ny2+1:ny2+p2) = y % array(:,1:p2        )

        end associate

        end associate

      class default
        err_msg = "y must be of type sll_t_vector_space_real_array_2d"
        SLL_ERROR( this_sub_name, err_msg )

      end select

    class default
      err_msg = "x must be of type sll_t_vector_space_real_array_1d"
      SLL_ERROR( this_sub_name, err_msg )

    end select

  end subroutine s_linear_operator_matrix_dense_to_stencil__dot

  ! Implement y=y+Ax, with A dense matrix (2D), x dense vector (1D) and y stencil vector (2D)
  subroutine s_linear_operator_matrix_dense_to_stencil__dot_incr( self, x, y )
    class(sll_t_linear_operator_matrix_dense_to_stencil), intent(in   ) :: self
    class(sll_c_vector_space)                           , intent(in   ) :: x ! dense
    class(sll_c_vector_space)                           , intent(inout) :: y ! stencil

    integer :: j1, j2, j, k

    character(len=*), parameter :: this_sub_name = "sll_t_linear_operator_matrix_dense_to_stencil % dot_incr"
    character(len=64) :: err_msg

    select type ( x )

    type is ( sll_t_vector_space_real_array_1d )

      ! Check dimensions
      SLL_ASSERT( self % s2 == size( x % array ) )

      select type ( y )

      type is ( sll_t_vector_space_real_array_2d )

        associate( p1 => -lbound( y % array, 1 ) + 1, &
                   p2 => -lbound( y % array, 2 ) + 1 )

        associate( ny1 => ubound( y % array, 1 ) - p1, &
                   ny2 => ubound( y % array, 2 ) - p2 )

          ! Check dimensions
          SLL_ASSERT( self % s1 == ny1*ny2 )

          do j2 = 1, ny2
            do j1 = 1, ny1
              j = (j1-1) * ny2 + j2
              do k = 1, size( x % array )
                y % array(j1,j2) = y % array(j1,j2) + self % A(j,k) * x % array(k)
              end do
            end do
          end do

          ! Update buffer regions
          y % array(1-p1:0      ,:) = y % array(ny1-p1+1:ny1,:)
          y % array(ny1+1:ny1+p1,:) = y % array(1:p1        ,:)
          y % array(:,1-p2:0      ) = y % array(:,ny2-p2+1:ny2)
          y % array(:,ny2+1:ny2+p2) = y % array(:,1:p2        )

        end associate

        end associate

      class default
        err_msg = "y must be of type sll_t_vector_space_real_array_2d"
        SLL_ERROR( this_sub_name, err_msg )

      end select

    class default
      err_msg = "x must be of type sll_t_vector_space_real_array_1d"
      SLL_ERROR( this_sub_name, err_msg )

    end select

  end subroutine s_linear_operator_matrix_dense_to_stencil__dot_incr

  ! Convert dense matrix to array (trivial)
  subroutine s_linear_operator_matrix_dense_to_stencil__to_array( self, A )
    class(sll_t_linear_operator_matrix_dense_to_stencil), intent(in   ) :: self
    real(wp)                                            , intent(inout) :: A(:,:)

    SLL_ASSERT( size( A, 1 ) == self % s1 )
    SLL_ASSERT( size( A, 2 ) == self % s2 )

    A = self % A

  end subroutine s_linear_operator_matrix_dense_to_stencil__to_array

  ! Free objects
  subroutine s_linear_operator_matrix_dense_to_stencil__free( self )
    class(sll_t_linear_operator_matrix_dense_to_stencil), intent(inout) :: self

    deallocate( self % A )

  end subroutine s_linear_operator_matrix_dense_to_stencil__free

end module sll_m_linear_operator_matrix_dense_to_stencil
