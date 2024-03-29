program test_gaussian_deviate
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"

   use sll_m_gaussian

   implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   sll_int32, parameter :: num_particles = 1000000

   sll_int32  :: j
   sll_real64, allocatable :: r(:, :)
   sll_real64, allocatable :: z(:)
   sll_real64 :: z_mean, z_variance
   sll_real64 :: r1_mean, r1_variance
   sll_real64 :: r2_mean, r2_variance

   allocate (z(num_particles))
   allocate (r(2, num_particles))
   do j = 1, num_particles

      z(j) = sll_f_gaussian_deviate()
      call sll_s_gaussian_deviate_2d(r(:, j))

   end do

   z_mean = sum(z)/real(num_particles, f64)
   z_variance = sum((z - z_mean)*(z - z_mean))/real(num_particles, f64)
   r1_mean = sum(r(1, :))/real(num_particles, f64)
   r2_mean = sum(r(2, :))/real(num_particles, f64)
   r1_variance = sum((r(1, :) - r1_mean)*(r(1, :) - r1_mean))/real(num_particles, f64)
   r2_variance = sum((r(2, :) - r2_mean)*(r(2, :) - r2_mean))/real(num_particles, f64)

   print *, ' z mean       = ', abs(z_mean)
   print *, ' z variance   = ', abs(z_variance - 1.0_f64)
   print *, ' r1 mean      = ', abs(r1_mean)
   print *, ' r1 variance  = ', abs(r1_variance - 1.0_f64)
   print *, ' r2 mean      = ', abs(r2_mean)
   print *, ' r2 variance  = ', abs(r2_variance - 1.0_f64)
   write (*, *) "-------"
   write (*, *) "End of Program 'unit_test' for 'deviators' "
   if (abs(z_mean) < 1d-2 .and. abs(z_variance - 1.0_f64) < 1d-2 .and. &
       abs(r1_mean) < 1d-2 .and. abs(r1_variance - 1.0_f64) < 1d-2 .and. &
       abs(r2_mean) < 1d-2 .and. abs(r2_variance - 1.0_f64) < 1d-2) then
      print *, 'PASSED'
   else
      print *, 'FAILED'
   end if
   deallocate (z)

end program test_gaussian_deviate
