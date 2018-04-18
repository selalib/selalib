module sll_m_linear_operator_matrix_stencil
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_errors.h"

  use sll_m_working_precision, only: f64

  use sll_m_linear_operator_base, only: sll_c_linear_operator

  use sll_m_vector_space_base, only: sll_c_vector_space

  use sll_m_vector_space_real_array_2d, only: sll_t_vector_space_real_array_2d

  implicit none

  public :: sll_t_linear_operator_matrix_stencil

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! Working precision
  integer, parameter :: wp = f64

  type, extends(sll_c_linear_operator) :: sll_t_linear_operator_matrix_stencil

    private

    real(wp), pointer :: A(:,:,:,:) => null()

    integer :: n1
    integer :: n2
    integer :: p1
    integer :: p2

  contains

    procedure :: init      => s_linear_operator_matrix_stencil__init
    procedure :: get_shape => f_linear_operator_matrix_stencil__get_shape
    procedure :: dot       => s_linear_operator_matrix_stencil__dot
    procedure :: to_array  => s_linear_operator_matrix_stencil__to_array
    procedure :: free      => s_linear_operator_matrix_stencil__free

  end type sll_t_linear_operator_matrix_stencil

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! Initialize linear operator
  subroutine s_linear_operator_matrix_stencil__init( self, l1, l2, A )
    class(sll_t_linear_operator_matrix_stencil), intent(inout) :: self
    integer                                    , intent(in   ) :: l1
    integer                                    , intent(in   ) :: l2
    real(wp), target                           , intent(in   ) :: A(l1:,l2:,:,:)

    self % A => A

    self % n1 = size( A, 3 )
    self % n2 = size( A, 4 )
    self % p1 = ubound( A, 1 )
    self % p2 = ubound( A, 2 )

    SLL_ASSERT( self % p1 == -l1 )
    SLL_ASSERT( self % p2 == -l2 )

  end subroutine s_linear_operator_matrix_stencil__init

  ! Get shape of linear operator
  function f_linear_operator_matrix_stencil__get_shape( self ) result( s )
    class(sll_t_linear_operator_matrix_stencil), intent(in) :: self
    integer :: s(2)

    s(1) = size( self % A, 3 ) * size( self % A, 4 )
    s(2) = s(1)

  end function f_linear_operator_matrix_stencil__get_shape

  ! Implement Ax=y, with A stencil matrix, x and y real 2D arrays
  subroutine s_linear_operator_matrix_stencil__dot( self, x, y )
    class(sll_t_linear_operator_matrix_stencil), intent(in   ) :: self
    class(sll_c_vector_space)                  , intent(in   ) :: x
    class(sll_c_vector_space)                  , intent(inout) :: y ! already constructed

    integer :: nx(2), ny(2)
    integer :: i1, i2, j1, j2, k1, k2

    character(len=*), parameter :: this_sub_name = "sll_t_linear_operator_matrix_stencil % dot"
    character(len=64) :: err_msg

    associate( n1 => self % n1, &
               n2 => self % n2, &
               p1 => self % p1, &
               p2 => self % p2 )

      ! Make sure to work with 1D real arrays
      select type ( x )

      type is ( sll_t_vector_space_real_array_2d )

        ! TODO: check dimensions
        nx = shape( x % array )

        select type ( y )

        type is ( sll_t_vector_space_real_array_2d )

          ! TODO: check dimensions
          ny = shape( y % array )

          y % array = 0.0_wp
          do i2 = 1, n2
            do i1 = 1, n1
              do k2 = -p2, p2
                do k1 = -p1, p1
                  j1 = modulo( i1 + k1, n1 ) + 1
                  j2 = modulo( i2 + k2, n2 ) + 1
                  y % array(i1,i2) = y % array(i1,i2) + self % A(k1,k2,i1,i2) * x % array(j1,j2)
                end do
              end do
            end do
          end do

          ! Update buffer regions
          y % array(1-p1:0    ,:) = y % array(n1-p1+1:n1,:)
          y % array(n1+1:n1+p1,:) = y % array(1:p1      ,:)
          y % array(:,1-p2:0    ) = y % array(:,n2-p2+1:n2)
          y % array(:,n2+1:n2+p2) = y % array(:,1:p2      )

        class default
          err_msg = "y must be of type sll_t_vector_space_real_array_2d"
          SLL_ERROR( this_sub_name, err_msg )

        end select

      class default
        err_msg = "x must be of type sll_t_vector_space_real_array_2d"
        SLL_ERROR( this_sub_name, err_msg )

      end select

    end associate

  end subroutine s_linear_operator_matrix_stencil__dot

  ! Convert stencil matrix to dense matrix
  subroutine s_linear_operator_matrix_stencil__to_array( self, A )
    class(sll_t_linear_operator_matrix_stencil), intent(in   ) :: self
    real(wp)                                   , intent(inout) :: A(:,:)

    integer :: i, j, i1, i2, j1, j2, k1, k2

    associate( n1 => self % n1, &
               n2 => self % n2, &
               p1 => self % p1, &
               p2 => self % p2 )

      SLL_ASSERT( size( A, 1 ) == n1*n2 )
      SLL_ASSERT( size( A, 2 ) == n1*n2 )

      do i2 = 1, n2
        do i1 = 1, n1
          do k2 = -p2, p2
            do k1 = -p1, p1
              j1 = modulo( i1 + k1, n1 ) + 1
              j2 = modulo( i2 + k2, n2 ) + 1
              i  = (i1-1) * n2 + i2
              j  = (j1-1) * n2 + j2
              A(i,j) = self % A(k1,k2,i1,i2)
            end do
          end do
        end do
      end do

    end associate

  end subroutine s_linear_operator_matrix_stencil__to_array

  ! Free objects
  subroutine s_linear_operator_matrix_stencil__free( self )
    class(sll_t_linear_operator_matrix_stencil), intent(inout) :: self

    self % A  => null()

  end subroutine s_linear_operator_matrix_stencil__free

end module sll_m_linear_operator_matrix_stencil
