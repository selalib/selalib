program test_hammersley
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"

   use sll_m_hammersley, only: sll_f_suite_hamm

   implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#define num_particles 50000_i32

   sll_int32  :: j
   sll_real64, allocatable :: z(:)
   sll_real64 :: z_mean, z_variance

   allocate (z(num_particles))

   do j = 1, num_particles

      z(j) = sll_f_suite_hamm(j, 2)

   end do

   z_mean = sum(z)/real(num_particles, f64)
   z_variance = 6*sum((z - z_mean)*(z - z_mean))/real(num_particles, f64)

   if (abs(z_mean - 0.5_f64) < 1d-3 .and. abs(z_variance - 0.5_f64) < 1d-3) then
      print *, 'PASSED'
   else
      print *, 'FAILED'
   end if
end program test_hammersley
