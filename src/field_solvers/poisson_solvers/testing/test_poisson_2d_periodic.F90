program test_poisson_2d_periodic
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

use sll_m_constants
use sll_m_poisson_2d_base
use sll_m_poisson_2d_periodic

implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

sll_int32   :: nc_eta1
sll_int32   :: nc_eta2
sll_int32   :: info
sll_real64  :: eta1_max
sll_real64  :: eta1_min
sll_real64  :: eta2_max
sll_real64  :: eta2_min
sll_real64  :: error

sll_real64, dimension(:,:), allocatable :: ex
sll_real64, dimension(:,:), allocatable :: ey
sll_real64, dimension(:,:), allocatable :: ex_exact
sll_real64, dimension(:,:), allocatable :: ey_exact
sll_real64, dimension(:,:), allocatable :: rho
sll_real64, dimension(:,:), allocatable :: rhs
sll_real64, dimension(:,:), allocatable :: phi
sll_real64, dimension(:,:), allocatable :: phi_exact

type(sll_t_poisson_2d_periodic_fft) :: poisson

sll_real64                :: x1
sll_real64                :: x2
sll_int32                 :: mode
sll_int32                 :: i
sll_int32                 :: j

class(sll_c_poisson_2d_base), pointer :: poisson_class 

sll_real64 :: x1_min
sll_real64 :: x1_max
sll_real64 :: x2_min
sll_real64 :: x2_max
sll_int32  :: Nc_x1
sll_int32  :: Nc_x2

sll_real64, dimension(:,:), allocatable :: E1
sll_real64, dimension(:,:), allocatable :: E2

sll_int32 :: ierr
  
eta1_min = .0_f64; eta1_max = 2.0_f64*sll_p_pi
eta2_min = .0_f64; eta2_max = 2.0_f64*sll_p_pi

nc_eta1 = 128; nc_eta2 = 128

SLL_CLEAR_ALLOCATE(ex(nc_eta1+1,nc_eta2+1),info)
SLL_CLEAR_ALLOCATE(ey(nc_eta1+1,nc_eta2+1),info)
SLL_CLEAR_ALLOCATE(ex_exact(nc_eta1+1,nc_eta2+1),info)
SLL_CLEAR_ALLOCATE(ey_exact(nc_eta1+1,nc_eta2+1),info)
SLL_CLEAR_ALLOCATE(rhs(nc_eta1+1,nc_eta2+1),info)
SLL_CLEAR_ALLOCATE(rho(nc_eta1+1,nc_eta2+1),info)
SLL_CLEAR_ALLOCATE(phi(nc_eta1+1,nc_eta2+1),info)
SLL_CLEAR_ALLOCATE(phi_exact(nc_eta1+1,nc_eta2+1),info)

write(*,*) " eta1_min, eta1_max, nc_eta1 ", eta1_min, eta1_max, nc_eta1
write(*,*) " eta2_min, eta2_max, nc_eta2 ", eta2_min, eta2_max, nc_eta2

call sll_o_initialize( poisson, eta1_min, eta1_max, nc_eta1, &
                 eta2_min, eta2_max, nc_eta2, info) 

open(14, file="test_poisson_2d_rho.dat")
mode = 2
do i = 1, nc_eta1+1
   do j = 1, nc_eta2+1
      x1 = (i-1)*(eta1_max-eta1_min)/nc_eta1
      x2 = (j-1)*(eta2_max-eta2_min)/nc_eta2
      phi_exact(i,j) = real(mode,f64) * sin(mode*x1) * cos(mode*x2)
      ex_exact(i,j)  =  1._f64*real(mode,f64)**2*cos(mode*x1)*cos(mode*x2)
      ey_exact(i,j)  = -1._f64*real(mode,f64)**2*sin(mode*x1)*sin(mode*x2)
      rho(i,j) = -2._f64 * real(mode,f64)**3 * sin(mode*x1)*cos(mode*x2)
      write(14,*) x1, x2, rho(i,j)
   end do
end do

rhs = rho
call sll_o_solve( poisson, phi, rhs)
error =  maxval(abs(phi_exact+phi))
write(*,*) " Po Error = " , error
if (error > 1e-13) stop 'FAILED'
rhs = rho
call sll_o_solve( poisson, phi, rhs)
error = maxval(abs(phi_exact+phi))
write(*,*) " Po Error = " ,  error
if (error > 1e-13) stop 'FAILED'
rhs = rho
call sll_o_solve( poisson, ex, ey, rhs)
error = maxval(abs(ex_exact-ex))
write(*,*) " Ex Error = " , error
if (error > 1e-13) stop 'FAILED'
error = maxval(abs(ey_exact-ey))
write(*,*) " Ey Error = " , error
if (error > 1e-13) stop 'FAILED'
rhs = rho
call sll_o_solve( poisson, ex, ey, rhs)
error = maxval(abs(ex_exact-ex))
write(*,*) " Ex Error = " , error
if (error > 1e-13) stop 'FAILED'
error = maxval(abs(ey_exact-ey))
write(*,*) " Ey Error = " , error
if (error > 1e-13) stop 'FAILED'

deallocate(rho, phi)

x1_min = 0._f64
x1_max = 1._f64

x2_min = 0._f64
x2_max = 1._f64

Nc_x1 = 32
Nc_x2 = 64

SLL_ALLOCATE(phi(Nc_x1+1,Nc_x2+1),ierr)
SLL_ALLOCATE(E1(Nc_x1+1,Nc_x2+1),ierr)
SLL_ALLOCATE(E2(Nc_x1+1,Nc_x2+1),ierr)
SLL_ALLOCATE(rho(Nc_x1+1,Nc_x2+1),ierr)

rho = 1._f64

poisson_class => sll_f_new_poisson_2d_periodic(x1_min,x1_max,Nc_x1,x2_min,x2_max,Nc_x2)

call poisson_class%compute_phi_from_rho( phi, rho )

call poisson_class%compute_E_from_rho( E1, E2, rho )

if( maxval(phi) == minval(phi)) then
 print *, '#PASSED'
else
 stop '#FAILED'
end if

end program test_poisson_2d_periodic
