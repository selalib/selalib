program sequential_advection

#include "sll_working_precision.h"
#include "sll_memory.h"
use sll_m_gnuplot

use sll_m_interpolators_2d_base
use sll_m_arbitrary_degree_spline_interpolator_2d

implicit none

class(sll_interpolator_2d_base), pointer                  :: interp
type(sll_arbitrary_degree_spline_interpolator_2d), target :: spl

sll_real64, dimension(:,:),  pointer       :: f

sll_int32  :: i_step
sll_real64 :: tcpu1
sll_real64 :: tcpu2
sll_int32  :: error
sll_real64 :: eta1
sll_real64 :: eta2

!#########################################################
!Simulation parameters and geometry sizes                !
                                                         !
sll_int32,  parameter :: nc_eta1 = 64                    !
sll_int32,  parameter :: nc_eta2 = 64                    !
sll_real64, parameter :: eta1_min = - 5.0_f64            !
sll_real64, parameter :: eta1_max =   5.0_f64            !
sll_real64, parameter :: eta2_min = - 5.0_f64            !
sll_real64, parameter :: eta2_max = + 5.0_f64            !
sll_real64 :: delta_eta1 = (eta1_max-eta1_min)/nc_eta1   !
sll_real64 :: delta_eta2 = (eta2_max-eta2_min)/nc_eta2   !
sll_real64, parameter :: delta_t = 0.05_f64              !
sll_int32,  parameter :: n_step  = 200                   !
sll_real64 :: alpha1
sll_real64 :: alpha2
                                                         !
!#########################################################

sll_int32  :: i, j

call cpu_time(tcpu1)

call spl%initialize(nc_eta1+1,    &
                    nc_eta2+1,    &
                    eta1_min,     &
                    eta1_max,     &
                    eta2_min,     &
                    eta2_max,     &
                    SLL_PERIODIC, &
                    SLL_PERIODIC, &
                    SLL_PERIODIC, &
                    SLL_PERIODIC, &
                    5, 5)

interp => spl

SLL_CLEAR_ALLOCATE(f(1:nc_eta1+1,1:nc_eta2+1),error)

do j=1,nc_eta2+1
do i=1,nc_eta1+1

   eta1  = eta1_min+(i-1)*delta_eta1
   eta2  = eta2_min+(j-1)*delta_eta2

   f(i,j)=exp(-.5*(eta1*eta1+eta2*eta2))

end do
end do

do i_step=1, 3*n_step

  call interp%compute_interpolants(f)
  alpha1 = delta_t
  alpha2 = delta_t
  do j=1,nc_eta2+1
     eta2 = eta2_min + (j-1)*delta_eta2 - alpha1
     do i = 1, nc_eta1+1
        eta1 = eta1_min + (i-1)*delta_eta1 - alpha2
        f(i,j) = interp%interpolate_value(eta1,eta2)
     end do
  end do

  call sll_gnuplot_2d( eta1_min, eta1_max, nc_eta1+1, &
                       eta2_min, eta2_max, nc_eta2+1, &
                       f, 'f_sequential', i_step, error)

  write(*,"(10x,' time = ', G15.3, ' sec' )") i_step*delta_t

end do

call cpu_time(tcpu2)
write(*,"(//10x,' Wall time = ', G15.3, ' sec' )") (tcpu2-tcpu1)

SLL_DEALLOCATE_ARRAY(f, error)


end program sequential_advection
