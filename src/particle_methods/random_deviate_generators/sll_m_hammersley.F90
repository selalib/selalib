!> hamm pseudo-random generator
!>
!>\author
!>\date created:
module sll_m_hammersley
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

   implicit none

   public :: &
      sll_f_suite_hamm

   private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

contains

!! !> @brief Returns
!! !> @param[in]
   function sll_f_suite_hamm(n, b)
      sll_int32, intent(in) :: n
      sll_int32, intent(in) :: b
      sll_real64 :: sll_f_suite_hamm

      sll_int32  :: m, u
      sll_real64 :: h, k

      k = 0.0_f64
      h = 0.0_f64
      m = n
      do while (m > 0)
         k = k + 1
         u = m/b
         h = h + b**(-k)*(m - u*b)
         m = u
      end do
      sll_f_suite_hamm = h
   end function sll_f_suite_hamm

   function vandercorput(n, p1, p2)!  p1>p2  !
      sll_int32, intent(in) :: n
      sll_int32, intent(in) :: p1
      sll_int32, intent(in) :: p2
      sll_real64 :: vandercorput

      sll_real64 :: s
      sll_int32  :: m

      m = n
      vandercorput = 0._f64
      s = 1.0_f64
      do while (m > 0)
         s = s/p1
         vandercorput = vandercorput + s*real(mod(p2*mod(m, p1), p1), f64)
         m = int(m/p1)
      end do
   end function vandercorput

end module sll_m_hammersley
