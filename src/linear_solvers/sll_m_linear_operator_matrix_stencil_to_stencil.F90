module sll_m_linear_operator_matrix_stencil_to_stencil
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_errors.h"

  use sll_m_working_precision, only: f64

  use sll_m_linear_operator_base, only: sll_c_linear_operator

  use sll_m_vector_space_base, only: sll_c_vector_space

  use sll_m_vector_space_real_array_2d, only: sll_t_vector_space_real_array_2d

  implicit none

  public :: sll_t_linear_operator_matrix_stencil_to_stencil

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! Working precision
  integer, parameter :: wp = f64

  type, extends(sll_c_linear_operator) :: sll_t_linear_operator_matrix_stencil_to_stencil

    real(wp), allocatable :: A(:,:,:,:)

    integer :: s1
    integer :: s2
    integer :: p1
    integer :: p2

  contains

    procedure :: init      => s_linear_operator_matrix_stencil_to_stencil__init
    procedure :: get_shape => f_linear_operator_matrix_stencil_to_stencil__get_shape
    procedure :: dot       => s_linear_operator_matrix_stencil_to_stencil__dot
    procedure :: dot_incr  => s_linear_operator_matrix_stencil_to_stencil__dot_incr
    procedure :: to_array  => s_linear_operator_matrix_stencil_to_stencil__to_array
    procedure :: free      => s_linear_operator_matrix_stencil_to_stencil__free

  end type sll_t_linear_operator_matrix_stencil_to_stencil

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! Initialize linear operator
  subroutine s_linear_operator_matrix_stencil_to_stencil__init( self, s1, s2, p1, p2 )
    class(sll_t_linear_operator_matrix_stencil_to_stencil), intent(inout) :: self
    integer                                               , intent(in   ) :: s1
    integer                                               , intent(in   ) :: s2
    integer                                               , intent(in   ) :: p1
    integer                                               , intent(in   ) :: p2

    allocate( self % A( -p1:p1, -p2:p2, 1:s1, 1:s2 ) )

    self % s1 = s1
    self % s2 = s2
    self % p1 = p1
    self % p2 = p2

  end subroutine s_linear_operator_matrix_stencil_to_stencil__init

  ! Get shape of linear operator
  function f_linear_operator_matrix_stencil_to_stencil__get_shape( self ) result( s )
    class(sll_t_linear_operator_matrix_stencil_to_stencil), intent(in) :: self
    integer :: s(2)

    s(1) = size( self % A, 3 ) * size( self % A, 4 )
    s(2) = s(1)

  end function f_linear_operator_matrix_stencil_to_stencil__get_shape

  ! Implement y=Ax, with A stencil matrix (4D), x and y stencil vectors (2D)
  subroutine s_linear_operator_matrix_stencil_to_stencil__dot( self, x, y )
    class(sll_t_linear_operator_matrix_stencil_to_stencil), intent(in   ) :: self
    class(sll_c_vector_space)                             , intent(in   ) :: x ! stencil
    class(sll_c_vector_space)                             , intent(inout) :: y ! stencil

    integer :: i1, i2, j1, j2, k1, k2
    real(wp) :: temp

    character(len=*), parameter :: this_sub_name = "sll_t_linear_operator_matrix_stencil_to_stencil % dot"
    character(len=64) :: err_msg

    associate( s1 => self % s1, &
               s2 => self % s2, &
               p1 => self % p1, &
               p2 => self % p2 )

      select type ( x )

      type is ( sll_t_vector_space_real_array_2d )

        ! Check dimensions
        SLL_ASSERT( self % p1 == -lbound( x % array, 1 ) + 1 )
        SLL_ASSERT( self % p2 == -lbound( x % array, 2 ) + 1 )
        SLL_ASSERT( self % s1 ==  ubound( x % array, 1 ) + lbound( x % array, 1 ) - 1 )
        SLL_ASSERT( self % s2 ==  ubound( x % array, 2 ) + lbound( x % array, 2 ) - 1 )

        select type ( y )

        type is ( sll_t_vector_space_real_array_2d )

          ! Check dimensions
          SLL_ASSERT( self % p1 == -lbound( y % array, 1 ) + 1 )
          SLL_ASSERT( self % p2 == -lbound( y % array, 2 ) + 1 )
          SLL_ASSERT( self % s1 ==  ubound( y % array, 1 ) + lbound( y % array, 1 ) - 1 )
          SLL_ASSERT( self % s2 ==  ubound( y % array, 2 ) + lbound( y % array, 2 ) - 1 )

          y % array = 0.0_wp

          !$OMP PARALLEL DO PRIVATE(j1,j2,temp)
          do i2 = 1, s2
            do i1 = 1, s1

              temp = 0.0_wp

              do k2 = -p2, p2
                do k1 = -p1, p1

                  j1 = i1+k1
                  j2 = i2+k2

                  ! Hand-made modulo operation
                  if (j2 < 1) then
                    j2 = j2+s2
                  else if (j2 > s2) then
                    j2 = j2-s2
                  end if

                  temp = temp + self % A(k1,k2,i1,i2) * x % array(j1,j2)

                end do
              end do

              y % array(i1,i2) = y % array(i1,i2) + temp

            end do
          end do
          !$OMP END PARALLEL DO

          ! Update buffer regions
          y % array(1-p1:0    ,:) = y % array(s1-p1+1:s1,:)
          y % array(s1+1:s1+p1,:) = y % array(1:p1      ,:)
          y % array(:,1-p2:0    ) = y % array(:,s2-p2+1:s2)
          y % array(:,s2+1:s2+p2) = y % array(:,1:p2      )

        class default
          err_msg = "y must be of type sll_t_vector_space_real_array_2d"
          SLL_ERROR( this_sub_name, err_msg )

        end select

      class default
        err_msg = "x must be of type sll_t_vector_space_real_array_2d"
        SLL_ERROR( this_sub_name, err_msg )

      end select

    end associate

  end subroutine s_linear_operator_matrix_stencil_to_stencil__dot

  ! Implement y=y+Ax, with A stencil matrix (4D), x and y stencil vectors (2D)
  subroutine s_linear_operator_matrix_stencil_to_stencil__dot_incr( self, x, y )
    class(sll_t_linear_operator_matrix_stencil_to_stencil), intent(in   ) :: self
    class(sll_c_vector_space)                             , intent(in   ) :: x ! stencil
    class(sll_c_vector_space)                             , intent(inout) :: y ! stencil

    integer :: i1, i2, j1, j2, k1, k2
    real(wp) :: temp

    character(len=*), parameter :: this_sub_name = "sll_t_linear_operator_matrix_stencil_to_stencil % dot"
    character(len=64) :: err_msg

    associate( s1 => self % s1, &
               s2 => self % s2, &
               p1 => self % p1, &
               p2 => self % p2 )

      select type ( x )

      type is ( sll_t_vector_space_real_array_2d )

        ! Check dimensions
        SLL_ASSERT( self % p1 == -lbound( x % array, 1 ) + 1 )
        SLL_ASSERT( self % p2 == -lbound( x % array, 2 ) + 1 )
        SLL_ASSERT( self % s1 ==  ubound( x % array, 1 ) + lbound( x % array, 1 ) - 1 )
        SLL_ASSERT( self % s2 ==  ubound( x % array, 2 ) + lbound( x % array, 2 ) - 1 )

        select type ( y )

        type is ( sll_t_vector_space_real_array_2d )

          ! Check dimensions
          SLL_ASSERT( self % p1 == -lbound( y % array, 1 ) + 1 )
          SLL_ASSERT( self % p2 == -lbound( y % array, 2 ) + 1 )
          SLL_ASSERT( self % s1 ==  ubound( y % array, 1 ) + lbound( y % array, 1 ) - 1 )
          SLL_ASSERT( self % s2 ==  ubound( y % array, 2 ) + lbound( y % array, 2 ) - 1 )

          !$OMP PARALLEL DO PRIVATE(j1,j2,temp)
          do i2 = 1, s2
            do i1 = 1, s1

              temp = 0.0_wp

              do k2 = -p2, p2
                do k1 = -p1, p1

                  j1 = i1+k1
                  j2 = i2+k2

                  ! Hand-made modulo operation
                  if (j2 < 1) then
                    j2 = j2+s2
                  else if (j2 > s2) then
                    j2 = j2-s2
                  end if

                  temp = temp + self % A(k1,k2,i1,i2) * x % array(j1,j2)

                end do
              end do

              y % array(i1,i2) = y % array(i1,i2) + temp

            end do
          end do
          !$OMP END PARALLEL DO

          ! Update buffer regions
          y % array(1-p1:0    ,:) = y % array(s1-p1+1:s1,:)
          y % array(s1+1:s1+p1,:) = y % array(1:p1      ,:)
          y % array(:,1-p2:0    ) = y % array(:,s2-p2+1:s2)
          y % array(:,s2+1:s2+p2) = y % array(:,1:p2      )

        class default
          err_msg = "y must be of type sll_t_vector_space_real_array_2d"
          SLL_ERROR( this_sub_name, err_msg )

        end select

      class default
        err_msg = "x must be of type sll_t_vector_space_real_array_2d"
        SLL_ERROR( this_sub_name, err_msg )

      end select

    end associate

  end subroutine s_linear_operator_matrix_stencil_to_stencil__dot_incr

  ! Convert stencil matrix to dense matrix
  subroutine s_linear_operator_matrix_stencil_to_stencil__to_array( self, A )
    class(sll_t_linear_operator_matrix_stencil_to_stencil), intent(in   ) :: self
    real(wp)                                              , intent(inout) :: A(:,:)

    integer :: i, j, i1, i2, j1, j2, k1, k2

    associate( s1 => self % s1, &
               s2 => self % s2, &
               p1 => self % p1, &
               p2 => self % p2 )

      SLL_ASSERT( size( A, 1 ) == s1*s2 )
      SLL_ASSERT( size( A, 2 ) == s1*s2 )

      do i2 = 1, s2
        do i1 = 1, s1
          do k2 = -p2, p2
            do k1 = -p1, p1
              j1 = modulo( i1 - 1 + k1, s1 ) + 1
              j2 = modulo( i2 - 1 + k2, s2 ) + 1
              i  = (i1-1) * s2 + i2
              j  = (j1-1) * s2 + j2
              A(i,j) = self % A(k1,k2,i1,i2)
            end do
          end do
        end do
      end do

    end associate

  end subroutine s_linear_operator_matrix_stencil_to_stencil__to_array

  ! Free objects
  subroutine s_linear_operator_matrix_stencil_to_stencil__free( self )
    class(sll_t_linear_operator_matrix_stencil_to_stencil), intent(inout) :: self

    deallocate( self % A )

  end subroutine s_linear_operator_matrix_stencil_to_stencil__free

end module sll_m_linear_operator_matrix_stencil_to_stencil
