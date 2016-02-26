program sequential_advection

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

use sll_m_fft
use sll_m_boundary_condition_descriptors, only: sll_p_periodic
use sll_m_gnuplot, only: sll_o_gnuplot_2d
use sll_m_constants

implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

sll_real64, dimension(:,:),  allocatable :: f, fe
sll_comp64, dimension(:,:),  allocatable :: fcmplx
sll_comp64, dimension(:,:),  allocatable :: ftilde

sll_int32  :: i_step
sll_real64 :: tcpu1
sll_real64 :: tcpu2
sll_int32  :: error
sll_real64 :: eta1
sll_real64 :: eta2

!#########################################################
!Simulation parameters and geometry sizes                !
                                                         !
sll_int32,  parameter :: nc_eta1 = 8                   !
sll_int32,  parameter :: nc_eta2 = 8                   !
sll_real64, parameter :: eta1_min =  0.0_f64
sll_real64, parameter :: eta1_max =  2*sll_p_pi           !
sll_real64, parameter :: eta2_min = - 0.0_f64
sll_real64, parameter :: eta2_max = + 2*sll_p_pi           !
sll_real64 :: delta_eta1 = (eta1_max-eta1_min)/nc_eta1   !
sll_real64 :: delta_eta2 = (eta2_max-eta2_min)/nc_eta2   !
sll_real64, parameter :: delta_t = 0.01_f64              !
sll_int32,  parameter :: n_step  = 500                   !
sll_real64 :: alpha1, x1
sll_real64 :: alpha2, x2
                                                         !
!#########################################################

sll_int32  :: i, j
sll_int32  :: m, n, p

sll_comp64, allocatable  :: ftemp_1d(:), f1d(:)
type(sll_t_fft) :: plan

sll_int32, allocatable  :: ntrace(:)
sll_int32, allocatable  :: mtrace(:)
sll_real64, allocatable :: x_array(:)
sll_real64, allocatable :: y_array(:)

sll_real64, parameter :: epsnufft=1.0d-9
sll_real64            :: time, x1c, x2c

call cpu_time(tcpu1)

SLL_ALLOCATE(x_array(1:nc_eta1*nc_eta2),error)
SLL_ALLOCATE(y_array(1:nc_eta1*nc_eta2),error)
SLL_ALLOCATE(mtrace(1:nc_eta1*nc_eta2),error)
SLL_ALLOCATE(ntrace(1:nc_eta2*nc_eta2),error)
SLL_ALLOCATE(f(1:nc_eta1,1:nc_eta2),error)
SLL_ALLOCATE(fe(1:nc_eta1,1:nc_eta2),error)
SLL_ALLOCATE(ftilde(1:nc_eta1,1:nc_eta2),error)
SLL_ALLOCATE(fcmplx(1:nc_eta1,1:nc_eta2),error)
SLL_ALLOCATE(f1d(1:nc_eta1),error)

do j=1,nc_eta2
do i=1,nc_eta1

   eta1  = eta1_min+(i-1)*delta_eta1-3*sll_p_pi*0.5
   eta2  = eta2_min+(j-1)*delta_eta2-sll_p_pi

   f(i,j)=exp(-.5*(eta1*eta1+eta2*eta2)*10.)

end do
end do

call sll_s_fft_init_c2c_2d(plan,nc_eta1,nc_eta2,fcmplx,ftilde,sll_p_fft_forward)

time = 0.0
do i_step=1, n_step

  time = time+delta_t

  fcmplx = cmplx(f,0.,f64)
  call sll_s_fft_exec_c2c_2d(plan, fcmplx, ftilde)
  ftilde = ftilde / cmplx(nc_eta1*nc_eta2,0.,f64)

  m=nc_eta1/2
  do n=1,nc_eta1
    f1d=ftilde(n,:)
    ftilde(n,:)=f1d((/ (/(p,p=m+1,nc_eta1)/), (/(p,p=1,m)/) /) )
  enddo
  do n=1,nc_eta2
    f1d=ftilde(:,n)
    ftilde(:,n)=f1d((/ (/(p,p=m+1,nc_eta1)/), (/(p,p=1,m)/) /) )
  enddo

  alpha1 = delta_t
  alpha2 = delta_t
  p=0
  do j=1,nc_eta2
    do i = 1, nc_eta1
      eta1 = eta1_min + (i-1)*delta_eta1 -sll_p_pi
      eta2 = eta2_min + (j-1)*delta_eta2 -sll_p_pi
      !x1   = eta1 - (eta2-sll_p_pi) * delta_t !cos(delta_t)*eta1-sin(delta_t)*eta2
      !x2   = eta2 + (eta1-sll_p_pi) * delta_t!sin(delta_t)*eta1+cos(delta_t)*eta2
      x1   = cos(delta_t)*eta1-sin(delta_t)*eta2 + sll_p_pi
      x2   = sin(delta_t)*eta1+cos(delta_t)*eta2 + sll_p_pi
     ! x1 = eta1_min + modulo(x1-eta1_min,eta1_max-eta1_min)
     ! x2 = eta2_min + modulo(x2-eta2_min,eta2_max-eta2_min)
      if ( eta1_min < x1 .and. x1 < eta1_max .and. &
           eta2_min < x2 .and. x2 < eta2_max ) then
        p=p+1
        ntrace(p)  = i
        mtrace(p)  = j
        x_array(p) = x1
        y_array(p) = x2
      else
        f(i,j) = 0.0_f64
      end if
      x1c   = cos(-time)*0.5*sll_p_pi
      x2c   = sin(-time)*0.5*sll_p_pi 
      fe(i,j)=exp(-.5*((eta1-x1c)**2+(eta2-x2c)**2)*10.)
    enddo
  enddo

  SLL_ALLOCATE(ftemp_1d(1:p), error)

  call nufft2d2f90(p,             &
                   x_array(1:p),  &
                   y_array(1:p),  &
                   ftemp_1d,      &
                         1,       &
                       epsnufft,  &
                      nc_eta1,    &
                      nc_eta2,    &
                       ftilde,    &
                       error)
  do n=1,p
    f(ntrace(n),mtrace(n))=real(ftemp_1d(n))
  enddo

  DEALLOCATE(ftemp_1d)

  call sll_o_gnuplot_2d( eta1_min, eta1_max, nc_eta1, &
                         eta2_min, eta2_max, nc_eta2, &
                         f, 'f_sequential', i_step, error)

  write(*,"(10x,' error = ', G15.3, i6 )") sum((f-fe)**2)/(nc_eta1*nc_eta2), nc_eta1
  

end do

call cpu_time(tcpu2)
write(*,"(//10x,' Wall time = ', G15.3, ' sec' )") (tcpu2-tcpu1)

call sll_s_fft_free(plan)
deallocate(f)
deallocate(ftilde)



end program sequential_advection
