module sll_m_linear_operator_matrix_c1_block
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_errors.h"

   use sll_m_working_precision, only: f64

   use sll_m_linear_operator_base, only: sll_c_linear_operator

   use sll_m_linear_operator_matrix_dense_to_dense, only: sll_t_linear_operator_matrix_dense_to_dense

   use sll_m_linear_operator_matrix_dense_to_stencil, only: sll_t_linear_operator_matrix_dense_to_stencil

   use sll_m_linear_operator_matrix_stencil_to_dense, only: sll_t_linear_operator_matrix_stencil_to_dense

   use sll_m_linear_operator_matrix_stencil_to_stencil, only: sll_t_linear_operator_matrix_stencil_to_stencil

   use sll_m_vector_space_base, only: sll_c_vector_space

   use sll_m_vector_space_c1_block, only: sll_t_vector_space_c1_block

   implicit none

   public :: sll_t_linear_operator_matrix_c1_block

   private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   ! Working precision
   integer, parameter :: wp = f64

   type, extends(sll_c_linear_operator) :: sll_t_linear_operator_matrix_c1_block

      integer :: n1(4)
      integer :: n2(4)
      integer :: p1
      integer :: p2

      type(sll_t_linear_operator_matrix_dense_to_dense) :: block1
      type(sll_t_linear_operator_matrix_stencil_to_dense) :: block2
      type(sll_t_linear_operator_matrix_dense_to_stencil) :: block3
      type(sll_t_linear_operator_matrix_stencil_to_stencil) :: block4

   contains

      procedure :: init => s_linear_operator_matrix_c1_block__init
      procedure :: get_shape => f_linear_operator_matrix_c1_block__get_shape
      procedure :: dot => s_linear_operator_matrix_c1_block__dot
      procedure :: dot_incr => s_linear_operator_matrix_c1_block__dot_incr
      procedure :: to_array => s_linear_operator_matrix_c1_block__to_array
      procedure :: free => s_linear_operator_matrix_c1_block__free

   end type sll_t_linear_operator_matrix_c1_block

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   ! Initialize linear operator
   subroutine s_linear_operator_matrix_c1_block__init(self, n1, n2, p1, p2)
      class(sll_t_linear_operator_matrix_c1_block), intent(inout) :: self
      integer, intent(in) :: n1(4)
      integer, intent(in) :: n2(4)
      integer, intent(in) :: p1
      integer, intent(in) :: p2

      self%n1 = n1
      self%n2 = n2
      self%p1 = p1
      self%p2 = p2

      call self%block1%init(n1(1), n2(1))
      call self%block2%init(n1(2), n2(2))
      call self%block3%init(n1(3), n2(3))
      call self%block4%init(n1(4), n2(4), p1, p2)

   end subroutine s_linear_operator_matrix_c1_block__init

   ! Get shape of linear operator
   function f_linear_operator_matrix_c1_block__get_shape(self) result(s)
      class(sll_t_linear_operator_matrix_c1_block), intent(in) :: self
      integer :: s(2)

      s = self%block1%get_shape() + self%block2%get_shape() + &
          self%block3%get_shape() + self%block4%get_shape()

   end function f_linear_operator_matrix_c1_block__get_shape

   ! Implement y=Ax, with A C1 block matrix, x and y C1 block vectors
   subroutine s_linear_operator_matrix_c1_block__dot(self, x, y)
      class(sll_t_linear_operator_matrix_c1_block), intent(in) :: self
      class(sll_c_vector_space), intent(in) :: x
      class(sll_c_vector_space), intent(inout) :: y

      character(len=*), parameter :: this_sub_name = "sll_t_linear_operator_matrix_c1_block % dot"
      character(len=64) :: err_msg

      select type (x)

      type is (sll_t_vector_space_c1_block)

         select type (y)

         type is (sll_t_vector_space_c1_block)

            call self%block1%dot(x%vd, y%vd)
            call self%block2%dot_incr(x%vs, y%vd)
            call self%block3%dot(x%vd, y%vs)
            call self%block4%dot_incr(x%vs, y%vs)

         class default
            err_msg = "y must be of type sll_t_vector_space_c1_block"
            SLL_ERROR(this_sub_name, err_msg)

         end select

      class default
         err_msg = "x must be of type sll_t_vector_space_c1_block"
         SLL_ERROR(this_sub_name, err_msg)

      end select

   end subroutine s_linear_operator_matrix_c1_block__dot

   ! Implement y=y+Ax, with A C1 block matrix, x and y C1 block vectors
   subroutine s_linear_operator_matrix_c1_block__dot_incr(self, x, y)
      class(sll_t_linear_operator_matrix_c1_block), intent(in) :: self
      class(sll_c_vector_space), intent(in) :: x
      class(sll_c_vector_space), intent(inout) :: y

      character(len=*), parameter :: this_sub_name = "sll_t_linear_operator_matrix_c1_block % dot_incr"
      character(len=64) :: err_msg

      select type (x)

      type is (sll_t_vector_space_c1_block)

         select type (y)

         type is (sll_t_vector_space_c1_block)

            call self%block1%dot_incr(x%vd, y%vd)
            call self%block2%dot_incr(x%vs, y%vd)
            call self%block3%dot_incr(x%vd, y%vs)
            call self%block4%dot_incr(x%vs, y%vs)

         class default
            err_msg = "y must be of type sll_t_vector_space_c1_block"
            SLL_ERROR(this_sub_name, err_msg)

         end select

      class default
         err_msg = "x must be of type sll_t_vector_space_c1_block"
         SLL_ERROR(this_sub_name, err_msg)

      end select

   end subroutine s_linear_operator_matrix_c1_block__dot_incr

   ! Convert C1 block matrix to dense matrix
   subroutine s_linear_operator_matrix_c1_block__to_array(self, A)
      class(sll_t_linear_operator_matrix_c1_block), intent(in) :: self
      real(wp), intent(inout) :: A(:, :)

      integer :: i, j, i1, i2, j1, j2, k1, k2

      ! Indices 1,2,3 correspond to dense blocks
      SLL_ASSERT(size(A, 1) == self%n1(1) + self%n1(3))
      SLL_ASSERT(size(A, 2) == self%n2(1) + self%n2(2))

      associate (n1 => self%n1, &
                 n2 => self%n2, &
                 p1 => self%p1, &
                 p2 => self%p2)

         A(1:n1(1), 1:n2(1)) = self%block1%A(:, :)
         A(1:n1(1), n2(1) + 1:n2(1) + n2(2)) = self%block2%A(:, :)
         A(n1(1) + 1:n1(1) + n1(3), 1:n2(1)) = self%block3%A(:, :)

         do i2 = 1, n2(4)
            do i1 = 1, n1(4)
               do k2 = -p2, p2
                  do k1 = -p1, p1
                     j1 = modulo(i1 - 1 + k1, n1(4)) + 1
                     j2 = modulo(i2 - 1 + k2, n2(4)) + 1
                     i = (i1 - 1)*n2(4) + i2
                     j = (j1 - 1)*n2(4) + j2
                     A(i + 3, j + 3) = self%block4%A(k1, k2, i1, i2)
                  end do
               end do
            end do
         end do

      end associate

   end subroutine s_linear_operator_matrix_c1_block__to_array

   ! Free objects
   subroutine s_linear_operator_matrix_c1_block__free(self)
      class(sll_t_linear_operator_matrix_c1_block), intent(inout) :: self

      call self%block1%free()
      call self%block2%free()
      call self%block3%free()
      call self%block4%free()

   end subroutine s_linear_operator_matrix_c1_block__free

end module sll_m_linear_operator_matrix_c1_block
