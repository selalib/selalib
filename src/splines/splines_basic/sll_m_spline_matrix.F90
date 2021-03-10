!> @ingroup splines
!> @brief   Access point to matrix class providing factory function
!> @author  Yaman Güçlü  - IPP Garching
!> @author  Edoardo Zoni - IPP Garching

module sll_m_spline_matrix

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_errors.h"

   use sll_m_working_precision, only: f64

   use sll_m_spline_matrix_base, only: sll_c_spline_matrix
   use sll_m_spline_matrix_dense, only: sll_t_spline_matrix_dense
   use sll_m_spline_matrix_banded, only: sll_t_spline_matrix_banded
   use sll_m_spline_matrix_periodic_banded, only: sll_t_spline_matrix_periodic_banded

   implicit none

   public :: &
      sll_c_spline_matrix, &
      sll_s_spline_matrix_new

   private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   !> Working precision
   integer, parameter :: wp = f64

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   subroutine sll_s_spline_matrix_new(matrix, matrix_type, n, kl, ku)
      class(sll_c_spline_matrix), allocatable, intent(out) :: matrix
      character(len=*), intent(in) :: matrix_type
      integer, intent(in) :: n
      integer, intent(in) :: kl
      integer, intent(in) :: ku

      character(len=*), parameter :: this_sub_name = "sll_s_spline_matrix_new"
      character(len=256) :: msg

      ! Allocate correct linear solver type
      select case (matrix_type)

      case ("dense")
         allocate (sll_t_spline_matrix_dense  :: matrix)

      case ("banded")
         allocate (sll_t_spline_matrix_banded :: matrix)

      case ("periodic_banded")
         allocate (sll_t_spline_matrix_periodic_banded :: matrix)

      case default
         msg = "Unrecognized matrix type: "//trim(matrix_type)
         SLL_ERROR(this_sub_name, msg)

      end select

      ! Initialize matrix (and linear solver)
      select type (matrix)

      type is (sll_t_spline_matrix_dense)
         call matrix%init(n)

      type is (sll_t_spline_matrix_banded)
         call matrix%init(n, kl, ku)

      type is (sll_t_spline_matrix_periodic_banded)
         if (kl /= ku) then
            msg = "Periodic banded matrix: Schur complement solver requires kl=ku, using their maximum instead."
            SLL_WARNING(this_sub_name, msg)
            call matrix%init(n, max(kl, ku), max(kl, ku))
         else
            call matrix%init(n, kl, ku)
         end if

      end select

   end subroutine sll_s_spline_matrix_new

end module sll_m_spline_matrix
