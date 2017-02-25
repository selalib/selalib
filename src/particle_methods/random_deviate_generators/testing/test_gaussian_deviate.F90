program test_gaussian_deviate
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"

use sll_m_gaussian

implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#define Xmin  0._f64
#define Xmax  5._f64
#define num_particles 8000_i32
#define the_mean  -1._f64
#define the_spatial_variance  4._f64


sll_int32  :: j, incr_index
sll_real64 :: z, res(2)
sll_real64 :: shift
logical    :: is_passed = .TRUE.

!open(152,file='test_Hamm_Gauss.dat')

shift = 0.0734_f64
incr_index = 0
do j=1, num_particles

   z = sll_f_gaussian_deviate() * the_spatial_variance + (the_mean)
   ! This gives a Gaussian centered in 'the_mean'
   call sll_s_gaussian_deviate_2d(res)

enddo

print*, z 
write(*,*) "-------"
write(*,*) "End of Program 'unit_test' for 'deviators' "
if (is_passed) then
   print *, 'PASSED'
else
   print *, 'FAILED'
end if

end program test_gaussian_deviate
