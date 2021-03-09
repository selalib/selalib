module sll_m_linear_operator_matrix_dense_to_dense
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_errors.h"

   use sll_m_working_precision, only: f64

   use sll_m_linear_operator_base, only: sll_c_linear_operator

   use sll_m_vector_space_base, only: sll_c_vector_space

   use sll_m_vector_space_real_array_1d, only: sll_t_vector_space_real_array_1d

   implicit none

   public :: sll_t_linear_operator_matrix_dense_to_dense

   private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   ! Working precision
   integer, parameter :: wp = f64

   type, extends(sll_c_linear_operator) :: sll_t_linear_operator_matrix_dense_to_dense

      real(wp), allocatable :: A(:, :)
      real(wp), allocatable :: At(:, :)

      ! Logical flag telling whether A was given transposed or not
      logical :: transposed = .false.

      integer :: s1
      integer :: s2

   contains

      procedure :: init => s_linear_operator_matrix_dense_to_dense__init
      procedure :: get_shape => f_linear_operator_matrix_dense_to_dense__get_shape
      procedure :: dot => s_linear_operator_matrix_dense_to_dense__dot
      procedure :: dot_incr => s_linear_operator_matrix_dense_to_dense__dot_incr
      procedure :: free => s_linear_operator_matrix_dense_to_dense__free

   end type sll_t_linear_operator_matrix_dense_to_dense

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   ! Initialize linear operator
   subroutine s_linear_operator_matrix_dense_to_dense__init(self, s1, s2, transposed)
      class(sll_t_linear_operator_matrix_dense_to_dense), intent(inout) :: self
      integer, intent(in) :: s1
      integer, intent(in) :: s2
      logical, optional, intent(in) :: transposed

      if (present(transposed)) self%transposed = transposed

      if (self%transposed) then
         allocate (self%At(s1, s2))
      else
         allocate (self%A(s1, s2))
      end if

      self%s1 = s1
      self%s2 = s2

   end subroutine s_linear_operator_matrix_dense_to_dense__init

   ! Get shape of linear operator
   function f_linear_operator_matrix_dense_to_dense__get_shape(self) result(s)
      class(sll_t_linear_operator_matrix_dense_to_dense), intent(in) :: self
      integer :: s(2)

      if (self%transposed) then
         s = shape(self%At)
      else
         s = shape(self%A)
      end if

   end function f_linear_operator_matrix_dense_to_dense__get_shape

   ! Implement y=Ax, with A dense matrix (2D), x and y dense vectors (1D)
   subroutine s_linear_operator_matrix_dense_to_dense__dot(self, x, y)
      class(sll_t_linear_operator_matrix_dense_to_dense), intent(in) :: self
      class(sll_c_vector_space), intent(in) :: x ! dense
      class(sll_c_vector_space), intent(inout) :: y ! dense

      integer :: nx(1), ny(1), n(2)

      character(len=*), parameter :: this_sub_name = "sll_t_linear_operator_matrix_dense_to_dense % dot"
      character(len=64) :: err_msg

      n = self%get_shape()

      select type (x)

      type is (sll_t_vector_space_real_array_1d)

         ! Check dimensions
         nx = shape(x%array)
         SLL_ASSERT(n(2) == nx(1))

         select type (y)

         type is (sll_t_vector_space_real_array_1d)

            ! Check dimensions
            ny = shape(y%array)
            SLL_ASSERT(n(1) == ny(1))

            if (self%transposed) then
               y%array = matmul(x%array, self%At)
            else
               y%array = matmul(self%A, x%array)
            end if

         class default
            err_msg = "y must be of type sll_t_vector_space_real_array_1d"
            SLL_ERROR(this_sub_name, err_msg)

         end select

      class default
         err_msg = "x must be of type sll_t_vector_space_real_array_1d"
         SLL_ERROR(this_sub_name, err_msg)

      end select

   end subroutine s_linear_operator_matrix_dense_to_dense__dot

   ! Implement y=y+Ax, with A dense matrix (2D), x and y dense vectors (1D)
   subroutine s_linear_operator_matrix_dense_to_dense__dot_incr(self, x, y)
      class(sll_t_linear_operator_matrix_dense_to_dense), intent(in) :: self
      class(sll_c_vector_space), intent(in) :: x ! dense
      class(sll_c_vector_space), intent(inout) :: y ! dense

      integer :: nx(1), ny(1), n(2)

      character(len=*), parameter :: this_sub_name = "sll_t_linear_operator_matrix_dense_to_dense % dot"
      character(len=64) :: err_msg

      n = self%get_shape()

      select type (x)

      type is (sll_t_vector_space_real_array_1d)

         ! Check dimensions
         nx = shape(x%array)
         SLL_ASSERT(n(2) == nx(1))

         select type (y)

         type is (sll_t_vector_space_real_array_1d)

            ! Check dimensions
            ny = shape(y%array)
            SLL_ASSERT(n(1) == ny(1))

            if (self%transposed) then
               y%array = y%array + matmul(x%array, self%At)
            else
               y%array = y%array + matmul(self%A, x%array)
            end if

         class default
            err_msg = "y must be of type sll_t_vector_space_real_array_1d"
            SLL_ERROR(this_sub_name, err_msg)

         end select

      class default
         err_msg = "x must be of type sll_t_vector_space_real_array_1d"
         SLL_ERROR(this_sub_name, err_msg)

      end select

   end subroutine s_linear_operator_matrix_dense_to_dense__dot_incr

   ! Free objects
   subroutine s_linear_operator_matrix_dense_to_dense__free(self)
      class(sll_t_linear_operator_matrix_dense_to_dense), intent(inout) :: self

      if (self%transposed) then
         deallocate (self%At)
      else
         deallocate (self%A)
      end if

   end subroutine s_linear_operator_matrix_dense_to_dense__free

end module sll_m_linear_operator_matrix_dense_to_dense
