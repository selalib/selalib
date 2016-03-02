program bsl_advection_using_nufft

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

use sll_m_fft
use sll_m_boundary_condition_descriptors, only: sll_p_periodic
!use sll_m_gnuplot, only: sll_o_gnuplot_2d
use sll_m_constants
use sll_m_interpolators_2d_base, only: sll_c_interpolator_2d
use sll_m_nufft_interpolator_2d

implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

sll_real64, dimension(:,:),  allocatable :: f, g, h
sll_real64, dimension(:,:),  allocatable :: x, y

sll_int32  :: i_step
sll_real64 :: tcpu1
sll_real64 :: tcpu2
sll_int32  :: error
sll_real64 :: eta1
sll_real64 :: eta2

!#########################################################
!Simulation parameters and geometry sizes                !
                                                         !
sll_int32,  parameter :: nc_eta1 = 128                   !
sll_int32,  parameter :: nc_eta2 = 128                   !
sll_real64, parameter :: eta1_min =  0.0_f64             !
sll_real64, parameter :: eta1_max =   2*sll_p_pi         !
sll_real64, parameter :: eta2_min = - 0.0_f64            !
sll_real64, parameter :: eta2_max = + 2*sll_p_pi         !
sll_real64 :: delta_eta1 = (eta1_max-eta1_min)/nc_eta1   !
sll_real64 :: delta_eta2 = (eta2_max-eta2_min)/nc_eta2   !
sll_real64, parameter :: delta_t = 0.1_f64               !
                                                         !
!#########################################################

sll_int32  :: i, j
sll_real64 :: time
sll_int32  :: n_step  

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class(sll_c_interpolator_2d),     pointer :: interp2d
type(sll_t_nufft_interpolator_2d), target :: nufft2d
  
call cpu_time(tcpu1)

n_step = int(2.0_f64*sll_p_pi/delta_t)
  
call nufft2d%initialize( &
     nc_eta1+1,          &
     nc_eta2+1,          &
     eta1_min,           &
     eta1_max,           &
     eta2_min,           &
     eta2_max            )

interp2d => nufft2d

SLL_ALLOCATE(f(1:nc_eta1+1,1:nc_eta2+1),error)
SLL_ALLOCATE(g(1:nc_eta1+1,1:nc_eta2+1),error)
SLL_ALLOCATE(x(1:nc_eta1+1,1:nc_eta2+1),error)
SLL_ALLOCATE(y(1:nc_eta1+1,1:nc_eta2+1),error)

do j=1,nc_eta2+1
do i=1,nc_eta1+1

   eta1  = eta1_min+(i-1)*delta_eta1-0.5*(eta1_min+eta2_max)
   eta2  = eta2_min+(j-1)*delta_eta2-0.5*(eta1_min+eta1_max)

   f(i,j)=exp(-.5*(eta1*eta1+eta2*eta2)*10.)

end do
end do

h = f

time = 0.0_f64
do i_step=1, n_step

  time = time+delta_t

  do j = 1, nc_eta2+1
    do i = 1, nc_eta1+1
      eta1   = (i-1)*delta_eta1
      eta2   = (j-1)*delta_eta2 
      x(i,j) = eta1_min+modulo(eta1-delta_t,eta1_max-eta1_min)
      y(i,j) = eta2_min+modulo(eta2-delta_t,eta2_max-eta2_min)
    enddo
  enddo

  call interp2d%interpolate_array( nc_eta1+1, &
                                   nc_eta2+1, &
                                   f,         &
                                   x,         &
                                   y,         &
                                   g          )

!  call sll_o_gnuplot_2d( eta1_min, eta1_max, nc_eta1, &
!                         eta2_min, eta2_max, nc_eta2, &
!                         f, 'f_bsl', i_step, error)
end do

if(sum((h-g)**2)/real(nc_eta1*nc_eta2, f64) < 1d-2) then
  print*, '#PASSED'
end if

call cpu_time(tcpu2)
write(*,"(//10x,' Wall time = ', G15.3, ' sec' )") (tcpu2-tcpu1)

deallocate(f)

end program bsl_advection_using_nufft
