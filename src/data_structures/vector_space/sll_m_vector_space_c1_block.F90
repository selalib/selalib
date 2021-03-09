!> @ingroup vector_space
!> @authors Yaman Güçlü    - <yaman.guclu@gmail.com>
!> @authors Edoardo Zoni   - <edoardo.zoni@ipp.mpg.de>
!> @brief   Vector space for wrapping 2D Fortran real arrays
!
module sll_m_vector_space_c1_block
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_errors.h"

   use sll_m_working_precision, only: f64

   use sll_m_vector_space_base, only: sll_c_vector_space

   use sll_m_vector_space_real_array_1d, only: sll_t_vector_space_real_array_1d

   use sll_m_vector_space_real_array_2d, only: sll_t_vector_space_real_array_2d

   implicit none

   public :: sll_t_vector_space_c1_block

   private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   !> Working precision
   integer, parameter :: wp = f64

   type, extends(sll_c_vector_space) :: sll_t_vector_space_c1_block

      integer :: n1
      integer :: n2
      integer :: p1
      integer :: p2

      type(sll_t_vector_space_real_array_1d) :: vd ! dense component
      type(sll_t_vector_space_real_array_2d) :: vs ! stencil component

   contains

      procedure :: init => s_vector_space_c1_block__init

      !> @name Basic operations (overloading the abstract methods)
      procedure :: copy => s_vector_space_c1_block__copy
      procedure :: incr => s_vector_space_c1_block__incr
      procedure :: scal => s_vector_space_c1_block__scal

      !> @name Additional operations
      procedure :: add => s_vector_space_c1_block__add
      procedure :: mult => s_vector_space_c1_block__mult
      procedure :: mult_add => s_vector_space_c1_block__mult_add
      procedure :: incr_mult => s_vector_space_c1_block__incr_mult
      procedure :: lcmb => s_vector_space_c1_block__lcmb
      procedure :: incr_lcmb => s_vector_space_c1_block__incr_lcmb

      !> @name Optional subroutines and functions
      procedure :: norm => f_vector_space_c1_block__norm
      procedure :: inner => f_vector_space_c1_block__inner
      !> @}

   end type sll_t_vector_space_c1_block

   ! Error messages
   character(len=*), parameter :: wrong_type_x = "x not of type 'sll_t_vector_space_c1_block'"
   character(len=*), parameter :: wrong_type_y = "y not of type 'sll_t_vector_space_c1_block'"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   subroutine s_vector_space_c1_block__init(self, n1, n2, p1, p2)
      class(sll_t_vector_space_c1_block), intent(inout) :: self
      integer, intent(in) :: n1
      integer, intent(in) :: n2
      integer, intent(in) :: p1
      integer, intent(in) :: p2

      self%n1 = n1
      self%n2 = n2
      self%p1 = p1
      self%p2 = p2

   end subroutine s_vector_space_c1_block__init

   !-----------------------------------------------------------------------------
   subroutine s_vector_space_c1_block__copy(self, x)
      class(sll_t_vector_space_c1_block), intent(inout) :: self
      class(sll_c_vector_space), intent(in) :: x

      character(len=*), parameter :: this_sub_name = "sll_t_vector_space_c1_block % copy"

      select type (x)

      class is (sll_t_vector_space_c1_block)

         call self%vd%copy(x%vd)
         call self%vs%copy(x%vs)

      class default
         SLL_ERROR(this_sub_name, wrong_type_x)

      end select

   end subroutine s_vector_space_c1_block__copy

   !-----------------------------------------------------------------------------
   subroutine s_vector_space_c1_block__incr(self, x)
      class(sll_t_vector_space_c1_block), intent(inout) :: self
      class(sll_c_vector_space), intent(in) :: x

      character(len=*), parameter :: this_sub_name = "sll_t_vector_space_c1_block % incr"

      select type (x)

      class is (sll_t_vector_space_c1_block)

         call self%vd%incr(x%vd)
         call self%vs%incr(x%vs)

      class default
         SLL_ERROR(this_sub_name, wrong_type_x)

      end select

   end subroutine s_vector_space_c1_block__incr

   !-----------------------------------------------------------------------------
   subroutine s_vector_space_c1_block__scal(self, a)
      class(sll_t_vector_space_c1_block), intent(inout) :: self
      real(wp), intent(in) :: a

      call self%vd%scal(a)
      call self%vs%scal(a)

   end subroutine s_vector_space_c1_block__scal

   !-----------------------------------------------------------------------------
   subroutine s_vector_space_c1_block__add(self, x, y)
      class(sll_t_vector_space_c1_block), intent(inout) :: self
      class(sll_c_vector_space), intent(in) :: x
      class(sll_c_vector_space), intent(in) :: y

      character(len=*), parameter :: this_sub_name = "sll_t_vector_space_c1_block % add"

      select type (x)

      class is (sll_t_vector_space_c1_block)

         select type (y)

         class is (sll_t_vector_space_c1_block)

            call self%vd%add(x%vd, y%vd)
            call self%vs%add(x%vs, y%vs)

         class default

            SLL_ERROR(this_sub_name, wrong_type_y)

         end select

      class default

         SLL_ERROR(this_sub_name, wrong_type_x)

      end select

   end subroutine s_vector_space_c1_block__add

   !-----------------------------------------------------------------------------
   subroutine s_vector_space_c1_block__mult(self, a, x)
      class(sll_t_vector_space_c1_block), intent(inout) :: self
      real(wp), intent(in) :: a
      class(sll_c_vector_space), intent(in) :: x

      character(len=*), parameter :: this_sub_name = "sll_t_vector_space_c1_block % mult"

      select type (x)

      class is (sll_t_vector_space_c1_block)

         call self%vd%mult(a, x%vd)
         call self%vs%mult(a, x%vs)

      class default

         SLL_ERROR(this_sub_name, wrong_type_x)

      end select

   end subroutine s_vector_space_c1_block__mult

   !-----------------------------------------------------------------------------
   subroutine s_vector_space_c1_block__mult_add(self, a, x, y)
      class(sll_t_vector_space_c1_block), intent(inout) :: self
      real(wp), intent(in) :: a
      class(sll_c_vector_space), intent(in) :: x
      class(sll_c_vector_space), intent(in) :: y

      character(len=*), parameter :: this_sub_name = "sll_t_vector_space_c1_block % mult_add"

      select type (x)

      class is (sll_t_vector_space_c1_block)

         select type (y)

         class is (sll_t_vector_space_c1_block)

            call self%vd%mult_add(a, x%vd, y%vd)
            call self%vs%mult_add(a, x%vs, y%vs)

         class default

            SLL_ERROR(this_sub_name, wrong_type_y)

         end select

      class default

         SLL_ERROR(this_sub_name, wrong_type_x)

      end select

   end subroutine s_vector_space_c1_block__mult_add

   !-----------------------------------------------------------------------------
   subroutine s_vector_space_c1_block__incr_mult(self, a, x)
      class(sll_t_vector_space_c1_block), intent(inout) :: self
      real(wp), intent(in) :: a
      class(sll_c_vector_space), intent(in) :: x

      character(len=*), parameter :: this_sub_name = "sll_t_vector_space_c1_block % incr_mult"

      select type (x)

      class is (sll_t_vector_space_c1_block)

         call self%vd%incr_mult(a, x%vd)
         call self%vs%incr_mult(a, x%vs)

      class default

         SLL_ERROR(this_sub_name, wrong_type_x)

      end select

   end subroutine s_vector_space_c1_block__incr_mult

   !----------------------------------------------------------------------------
   subroutine s_vector_space_c1_block__lcmb(self, a, x)
      class(sll_t_vector_space_c1_block), intent(inout) :: self
      real(wp), intent(in) :: a(:)
      class(sll_c_vector_space), intent(in) :: x(:)

      character(len=*), parameter :: this_sub_name = "sll_t_vector_space_c1_block % lcmb"

      select type (x)

      class is (sll_t_vector_space_c1_block)

         call self%vd%lcmb(a, x%vd)
         call self%vs%lcmb(a, x%vs)

      class default

         SLL_ERROR(this_sub_name, wrong_type_x)

      end select

   end subroutine s_vector_space_c1_block__lcmb

   !-----------------------------------------------------------------------------
   subroutine s_vector_space_c1_block__incr_lcmb(self, a, x)
      class(sll_t_vector_space_c1_block), intent(inout) :: self
      real(wp), intent(in) :: a(:)
      class(sll_c_vector_space), intent(in) :: x(:)

      character(len=*), parameter :: this_sub_name = "sll_t_vector_space_c1_block % incr_lcmb"

      select type (x)

      class is (sll_t_vector_space_c1_block)

         call self%vd%incr_lcmb(a, x%vd)
         call self%vs%incr_lcmb(a, x%vs)

      class default

         SLL_ERROR(this_sub_name, wrong_type_x)

      end select

   end subroutine s_vector_space_c1_block__incr_lcmb

   !-----------------------------------------------------------------------------
   function f_vector_space_c1_block__norm(self) result(res)
      class(sll_t_vector_space_c1_block), intent(in) :: self
      real(wp) :: res

      res = sqrt(self%inner(self))

   end function f_vector_space_c1_block__norm

   !-----------------------------------------------------------------------------
   function f_vector_space_c1_block__inner(self, x) result(res)
      class(sll_t_vector_space_c1_block), intent(in) :: self
      class(sll_c_vector_space), intent(in) :: x
      real(wp) :: res

      character(len=*), parameter :: this_fun_name = "sll_t_vector_space_c1_block % inner"

      select type (x)

      class is (sll_t_vector_space_c1_block)

         associate (n1 => self%n1, n2 => self%n2)
            res = self%vd%inner(x%vd) + &
                  sum(self%vs%array(1:n1, 1:n2)*x%vs%array(1:n1, 1:n2))
         end associate

      class default

         SLL_ERROR(this_fun_name, wrong_type_x)

      end select

   end function f_vector_space_c1_block__inner

end module sll_m_vector_space_c1_block
