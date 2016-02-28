program advection_nufft_2d

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

use sll_m_fft
use sll_m_nufft_interpolation
use sll_m_boundary_condition_descriptors, only: sll_p_periodic
use sll_m_constants

implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

sll_real64, dimension(:,:),  allocatable :: f
sll_real64, dimension(:,:),  allocatable :: fe

sll_int32  :: i_step
sll_real64 :: tcpu1
sll_real64 :: tcpu2
sll_int32  :: error
sll_real64 :: eta1
sll_real64 :: eta2

!#########################################################
!Simulation parameters and geometry sizes                !
                                                         !
sll_int32,  parameter :: nc_eta1 = 64                   !
sll_int32,  parameter :: nc_eta2 = 128                   !
sll_real64, parameter :: eta1_min = - 5.0_f64            !
sll_real64, parameter :: eta1_max = + 5.0_f64            !
sll_real64, parameter :: eta2_min = - 5.0_f64            !
sll_real64, parameter :: eta2_max = + 5.0_f64            !
sll_real64 :: delta_eta1 = (eta1_max-eta1_min)/nc_eta1   !
sll_real64 :: delta_eta2 = (eta2_max-eta2_min)/nc_eta2   !
sll_real64, parameter :: delta_t = 0.1_f64               !
sll_int32,  parameter :: n_step  = 50                    !
                                                         !
!#########################################################

sll_int32  :: i, j
sll_real64 :: time, x1c, x2c
  
type(sll_t_nufft_2d) :: interp_2d

call cpu_time(tcpu1)

SLL_ALLOCATE(f(1:nc_eta1,1:nc_eta2),error)
SLL_ALLOCATE(fe(1:nc_eta1,1:nc_eta2),error)

do j=1,nc_eta2
do i=1,nc_eta1

   eta1  = eta1_min+(i-1)*delta_eta1
   eta2  = eta2_min+(j-1)*delta_eta2
   x1c   = (eta1_min+eta1_max)*0.5 + (eta1_max-eta1_min)*0.25
   x2c   = (eta1_min+eta1_max)*0.5 

   f(i,j) = exp(-.5*((eta1-x1c)**2+(eta2-x2c)**2)*10.)

end do
end do

call sll_s_nufft_2d_init( interp_2d,                   &
                          nc_eta1, eta1_min, eta1_max, &
                          nc_eta2, eta2_min, eta2_max)

time = 0.0_f64
do i_step=1, n_step

  time = time+delta_t

  call sll_s_nufft_2d_rotation( interp_2d, f, delta_t )

  do j=1,nc_eta2
    do i = 1, nc_eta1
      eta1    = eta1_min + (i-1)*delta_eta1 
      eta2    = eta2_min + (j-1)*delta_eta2 
      x1c     = cos(-time)*0.25*(eta1_max-eta1_min)+(eta1_min+eta1_max)*0.5
      x2c     = sin(-time)*0.25*(eta2_max-eta2_min)+(eta2_min+eta2_max)*0.5 
      fe(i,j) = exp(-.5*((eta1-x1c)**2+(eta2-x2c)**2)*10.)
    enddo
  enddo


end do

call cpu_time(tcpu2)
write(*,"(//10x,' Wall time = ', G15.3, ' sec' )") (tcpu2-tcpu1)


call sll_s_nufft_2d_free( interp_2d )

if (sum((f-fe)**2)/real(nc_eta1*nc_eta2,f64) < 1d-12) then
  print*,"#PASSED"
end if
deallocate(fe)
deallocate(f)

end program advection_nufft_2d
