!> @ingroup vector_space
!> @authors Yaman Güçlü    - <yaman.guclu@gmail.com>
!> @authors Marco Restelli - <marco.restelli@gmail.com>
!> @authors Edoardo Zoni   - <edoardo.zoni@ipp.mpg.de>
!> @brief   Vector space for wrapping 3D Fortran real arrays
!
module sll_m_vector_space_real_array_3d
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_errors.h"

   use sll_m_working_precision, only: f64

   use sll_m_vector_space_base, only: sll_c_vector_space

   implicit none

   public :: sll_t_vector_space_real_array_3d

   private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   !> Working precision
   integer, parameter :: wp = f64

   !> @brief   Vector space for wrapping 3D Fortran real arrays
   !> @details Concrete derived type providing the efficient implementation
   !>          of all basic operations on a vector space that wraps a single
   !>          3D Fortran real array
   type, extends(sll_c_vector_space) :: sll_t_vector_space_real_array_3d

      real(wp), allocatable :: array(:, :, :)

   contains

      !> @name Basic operations (overloading the abstract methods)
      procedure :: copy => copy__real
      procedure :: incr => incr__real
      procedure :: scal => scal__real

      !> @name Additional operations
      procedure :: add => add__real
      procedure :: mult => mult__real
      procedure :: mult_add => mult_add__real
      procedure :: incr_mult => incr_mult__real
      procedure :: lcmb => lcmb__real
      procedure :: incr_lcmb => incr_lcmb__real

      !> @name Optional subroutines and functions
      procedure :: norm => norm__real
      procedure :: inner => inner__real
      !> @}

   end type sll_t_vector_space_real_array_3d

   ! Error messages
   character(len=*), parameter :: wrong_type_x = "x not of type 'sll_t_vector_space_real_array_3d'"
   character(len=*), parameter :: wrong_type_y = "y not of type 'sll_t_vector_space_real_array_3d'"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   subroutine copy__real(self, x)
      class(sll_t_vector_space_real_array_3d), intent(inout) :: self
      class(sll_c_vector_space), intent(in) :: x

      character(len=*), parameter :: this_sub_name = "sll_t_vector_space_real_array_3d % copy"

      select type (x)

      class is (sll_t_vector_space_real_array_3d)

         ! Automatic allocation
         self%array(:, :, :) = x%array

      class default

         SLL_ERROR(this_sub_name, wrong_type_x)

      end select

   end subroutine copy__real

   !-----------------------------------------------------------------------------
   subroutine incr__real(self, x)
      class(sll_t_vector_space_real_array_3d), intent(inout) :: self
      class(sll_c_vector_space), intent(in) :: x

      character(len=*), parameter :: this_sub_name = "sll_t_vector_space_real_array_3d % incr"

      select type (x)

      class is (sll_t_vector_space_real_array_3d)

         self%array(:, :, :) = self%array(:, :, :) + x%array(:, :, :)

      class default

         SLL_ERROR(this_sub_name, wrong_type_x)

      end select

   end subroutine incr__real

   !-----------------------------------------------------------------------------
   subroutine scal__real(self, a)
      class(sll_t_vector_space_real_array_3d), intent(inout) :: self
      real(wp), intent(in) :: a

      self%array(:, :, :) = self%array(:, :, :)*a

   end subroutine scal__real

   !-----------------------------------------------------------------------------
   subroutine add__real(self, x, y)
      class(sll_t_vector_space_real_array_3d), intent(inout) :: self
      class(sll_c_vector_space), intent(in) :: x
      class(sll_c_vector_space), intent(in) :: y

      character(len=*), parameter :: this_sub_name = "sll_t_vector_space_real_array_3d % add"

      select type (x)

      class is (sll_t_vector_space_real_array_3d)

         select type (y)

         class is (sll_t_vector_space_real_array_3d)

            self%array(:, :, :) = x%array(:, :, :) + y%array(:, :, :)

         class default

            SLL_ERROR(this_sub_name, wrong_type_y)

         end select

      class default

         SLL_ERROR(this_sub_name, wrong_type_x)

      end select

   end subroutine add__real

   !-----------------------------------------------------------------------------
   subroutine mult__real(self, a, x)
      class(sll_t_vector_space_real_array_3d), intent(inout) :: self
      real(wp), intent(in) :: a
      class(sll_c_vector_space), intent(in) :: x

      character(len=*), parameter :: this_sub_name = "sll_t_vector_space_real_array_3d % mult"

      select type (x)

      class is (sll_t_vector_space_real_array_3d)

         self%array(:, :, :) = a*x%array(:, :, :)

      class default

         SLL_ERROR(this_sub_name, wrong_type_x)

      end select

   end subroutine mult__real

   !-----------------------------------------------------------------------------
   subroutine mult_add__real(self, a, x, y)
      class(sll_t_vector_space_real_array_3d), intent(inout) :: self
      real(wp), intent(in) :: a
      class(sll_c_vector_space), intent(in) :: x
      class(sll_c_vector_space), intent(in) :: y

      character(len=*), parameter :: this_sub_name = "sll_t_vector_space_real_array_3d % mult_add"

      select type (x)

      class is (sll_t_vector_space_real_array_3d)

         select type (y)

         class is (sll_t_vector_space_real_array_3d)

            self%array(:, :, :) = a*x%array(:, :, :) + y%array(:, :, :)

         class default

            SLL_ERROR(this_sub_name, wrong_type_y)

         end select

      class default

         SLL_ERROR(this_sub_name, wrong_type_x)

      end select

   end subroutine mult_add__real

   !-----------------------------------------------------------------------------
   subroutine incr_mult__real(self, a, x)
      class(sll_t_vector_space_real_array_3d), intent(inout) :: self
      real(wp), intent(in) :: a
      class(sll_c_vector_space), intent(in) :: x

      character(len=*), parameter :: this_sub_name = "sll_t_vector_space_real_array_3d % incr_mult"

      select type (x)

      class is (sll_t_vector_space_real_array_3d)

         self%array(:, :, :) = self%array(:, :, :) + a*x%array(:, :, :)

      class default

         SLL_ERROR(this_sub_name, wrong_type_x)

      end select

   end subroutine incr_mult__real

   !----------------------------------------------------------------------------
   subroutine lcmb__real(self, a, x)
      class(sll_t_vector_space_real_array_3d), intent(inout) :: self
      real(wp), intent(in) :: a(:)
      class(sll_c_vector_space), intent(in) :: x(:)

      character(len=*), parameter :: this_sub_name = "sll_t_vector_space_real_array_3d % lcmb"

      integer :: i

      select type (x)

      class is (sll_t_vector_space_real_array_3d)

         ! Construct linear combination incrementally
         self%array(:, :, :) = a(1)*x(1)%array(:, :, :)
         do i = 2, size(a)
            self%array(:, :, :) = self%array(:, :, :) + a(i)*x(i)%array(:, :, :)
         end do

      class default

         SLL_ERROR(this_sub_name, wrong_type_x)

      end select

   end subroutine lcmb__real

   !-----------------------------------------------------------------------------
   subroutine incr_lcmb__real(self, a, x)
      class(sll_t_vector_space_real_array_3d), intent(inout) :: self
      real(wp), intent(in) :: a(:)
      class(sll_c_vector_space), intent(in) :: x(:)

      character(len=*), parameter :: this_sub_name = "sll_t_vector_space_real_array_3d % incr_lcmb"

      integer :: i

      select type (x)

      class is (sll_t_vector_space_real_array_3d)

         ! Construct linear combination incrementally
         do i = 1, size(a)
            self%array(:, :, :) = self%array(:, :, :) + a(i)*x(i)%array(:, :, :)
         end do

      class default

         SLL_ERROR(this_sub_name, wrong_type_x)

      end select

   end subroutine incr_lcmb__real

   !-----------------------------------------------------------------------------
   function norm__real(self) result(res)
      class(sll_t_vector_space_real_array_3d), intent(in) :: self
      real(wp) :: res

#ifdef __PGI
      res = sqrt(self%inner(self))
#else
      res = norm2(self%array)
#endif

   end function norm__real

   !-----------------------------------------------------------------------------
   function inner__real(self, x) result(res)
      class(sll_t_vector_space_real_array_3d), intent(in) :: self
      class(sll_c_vector_space), intent(in) :: x
      real(wp) :: res

      character(len=*), parameter :: this_fun_name = "sll_t_vector_space_real_array_3d % inner"

      select type (x)

      class is (sll_t_vector_space_real_array_3d)

         res = sum(self%array*x%array)

      class default

         SLL_ERROR(this_fun_name, wrong_type_x)

      end select

   end function inner__real

end module sll_m_vector_space_real_array_3d
