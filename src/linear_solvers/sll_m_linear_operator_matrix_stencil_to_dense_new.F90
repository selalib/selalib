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

    real(wp), allocatable :: A (:,:,:)

    integer :: s1
    integer :: s2
    integer :: s3

  contains

    procedure :: init      => s_linear_operator_matrix_stencil_to_dense_new__init
    procedure :: get_shape => f_linear_operator_matrix_stencil_to_dense_new__get_shape
    procedure :: dot       => s_linear_operator_matrix_stencil_to_dense_new__dot
    procedure :: dot_incr  => s_linear_operator_matrix_stencil_to_dense_new__dot_incr
    procedure :: free      => s_linear_operator_matrix_stencil_to_dense_new__free

  end type sll_t_linear_operator_matrix_stencil_to_dense_new

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! Initialize linear operator
  subroutine s_linear_operator_matrix_stencil_to_dense_new__init( self, s1, s2, s3 )
    class(sll_t_linear_operator_matrix_stencil_to_dense_new), intent(inout) :: self
    integer                                                 , intent(in   ) :: s1
    integer                                                 , intent(in   ) :: s2
    integer                                                 , intent(in   ) :: s3

    allocate( self % A( s1, s2, s3 ) )

    self % s1 = s1
    self % s2 = s2
    self % s3 = s3

  end subroutine s_linear_operator_matrix_stencil_to_dense_new__init

  ! Get shape of linear operator
  function f_linear_operator_matrix_stencil_to_dense_new__get_shape( self ) result( s )
    class(sll_t_linear_operator_matrix_stencil_to_dense_new), intent(in) :: self
    integer :: s(2)

    s(1) = size( self % A, 1 )
    s(2) = size( self % A, 2 ) * size( self % A, 3 )

  end function f_linear_operator_matrix_stencil_to_dense_new__get_shape

  ! Implement y=Ax, with A dense matrix (2D), x stencil vector (2D) and y dense vector (1D)
  subroutine s_linear_operator_matrix_stencil_to_dense_new__dot( self, x, y )
    class(sll_t_linear_operator_matrix_stencil_to_dense_new), intent(in   ) :: self
    class(sll_c_vector_space)                               , intent(in   ) :: x ! stencil
    class(sll_c_vector_space)                               , intent(inout) :: y ! dense

    integer  :: j1, j2, i
    real(wp) :: temp

    character(len=*), parameter :: this_sub_name = "sll_t_linear_operator_matrix_stencil_to_dense_new % dot"
    character(len=64) :: err_msg

    select type ( x )

    type is ( sll_t_vector_space_real_array_2d )

      associate( p1 => -lbound( x % array, 1 ) + 1, &
                 p2 => -lbound( x % array, 2 ) + 1 )

      associate( nx1 => ubound( x % array, 1 ) - p1, &
                 nx2 => ubound( x % array, 2 ) - p2 )

        select type ( y )

        type is ( sll_t_vector_space_real_array_1d )

          ! Check dimensions
          SLL_ASSERT( self%s1 == size(y%array) )

          y % array = 0.0_wp

          !$OMP PARALLEL DO PRIVATE(temp)
          do i = 1, size( y % array )

            temp = 0.0_wp

            do j2 = 1, nx2
              do j1 = 1, p1
                temp = temp + self % A(i,j1,j2) * x % array(j1,j2)
              end do
            end do

            y % array(i) = y % array(i) + temp

          end do
          !$OMP END PARALLEL DO

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

    integer  :: j1, j2, i
    real(wp) :: temp

    character(len=*), parameter :: this_sub_name = "sll_t_linear_operator_matrix_stencil_to_dense_new % dot"
    character(len=64) :: err_msg

    select type ( x )

    type is ( sll_t_vector_space_real_array_2d )

      associate( p1 => -lbound( x % array, 1 ) + 1, &
                 p2 => -lbound( x % array, 2 ) + 1 )

      associate( nx1 => ubound( x % array, 1 ) - p1, &
                 nx2 => ubound( x % array, 2 ) - p2 )

        select type ( y )

        type is ( sll_t_vector_space_real_array_1d )

          ! Check dimensions
          SLL_ASSERT( self%s1 == size(y%array) )

          !$OMP PARALLEL DO PRIVATE(temp)
          do i = 1, size( y % array )

            temp = 0.0_wp

            do j2 = 1, nx2
              do j1 = 1, p1
                temp = temp + self % A(i,j1,j2) * x % array(j1,j2)
              end do
            end do

            y % array(i) = y % array(i) + temp

          end do
          !$OMP END PARALLEL DO

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

  ! Free objects
  subroutine s_linear_operator_matrix_stencil_to_dense_new__free( self )
    class(sll_t_linear_operator_matrix_stencil_to_dense_new), intent(inout) :: self

    deallocate( self % A )

  end subroutine s_linear_operator_matrix_stencil_to_dense_new__free

end module sll_m_linear_operator_matrix_stencil_to_dense_new
