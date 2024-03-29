module sll_m_gyroaverage_2d_base
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"

   implicit none

   public :: &
      sll_c_gyroaverage_2d_base

   private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   type, abstract :: sll_c_gyroaverage_2d_base
   contains
      procedure(signature_compute_gyroaverage_2d), deferred, pass(gyroaverage) :: &
         compute_gyroaverage
   end type sll_c_gyroaverage_2d_base

   !> Compute Jf = gyroaverage of f with the Larmor radius larmor_rad
   abstract interface
      subroutine signature_compute_gyroaverage_2d(gyroaverage, larmor_rad, f)
         use sll_m_working_precision
         import sll_c_gyroaverage_2d_base
         class(sll_c_gyroaverage_2d_base), target :: gyroaverage
         sll_real64, intent(in) :: larmor_rad
         sll_real64, dimension(:, :), intent(inout) :: f
      end subroutine signature_compute_gyroaverage_2d
   end interface

end module sll_m_gyroaverage_2d_base
