program periodic_interpolation
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"

  use sll_m_boundary_condition_descriptors, only: &
    sll_p_periodic

  use sll_m_constants, only: &
    sll_p_pi

  use sll_m_interpolators_1d_base, only: &
    sll_c_interpolator_1d

  use sll_m_lagrange_interpolator_1d, only: &
    sll_t_lagrange_interpolator_1d

  use sll_m_periodic_interp, only: &
    sll_p_lagrange, &
    sll_p_spline

  use sll_m_periodic_interpolator_1d, only: &
    sll_t_periodic_interpolator_1d

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
type(sll_t_periodic_interpolator_1d), target :: interp_per
type(sll_t_lagrange_interpolator_1d),     target :: interp_lagrange
sll_int32, parameter                       :: N0 = 16
sll_real64                                 :: u(16*N0+1)
sll_real64                                 :: u_exact(16*N0+1)
sll_real64                                 :: u_out(16*N0+1)
class(sll_c_interpolator_1d), pointer   :: interp
sll_real64, parameter :: xmin = 0.0_f64, xmax=3.0_f64   
sll_real64 :: alpha, error, old_error, L, xi,mode
sll_int32 :: i, p, N

alpha = 0.05_f64
mode = 3.0_f64
L = xmax - xmin

error = 1.0_f64
print*, 'Testing order of periodic interpolation (sll_p_spline)'
N = N0
do p=1,4
   N= 2*N 
   ! Interpolate non trivial smooth periodic function
   do i=0, N
     xi = xmin+real(i,f64)*L/real(N,f64)
     u(i+1) = 1.0_f64 / (2.0_f64 + sin(mode*2._f64*sll_p_pi/L*xi))
     u_exact(i+1) =  1.0_f64 / (2.0_f64 + sin(mode*2._f64*sll_p_pi/L*(xi-alpha)))
!    u(i+1) = cos(mode*twopi*i/N)
!    u_exact(i+1) = cos(mode*twopi*(i-alpha)/N)
   end do
   print*, 'p=', p
   call interp_per%initialize( N+1, xmin, xmax, sll_p_spline, 12)
   interp => interp_per
   call interp%interpolate_array_disp( N+1, u(1:N+1), -alpha, u_out(1:N+1))
   old_error = error
   error = maxval(abs(u_out(1:N+1)-u_exact(1:N+1)))
   print *, "error =", error
   print*, 'N=',N, 'error=', error, 'numerical order=', log(old_error/error)/log(2.0_f64)
end do

error = 1.0_f64
print*, 'Testing order of periodic interpolation (sll_p_lagrange)'
N = N0
do p=1,4
   N= 2*N 
   ! Interpolate non trivial smooth periodic function
   do i=0, N
     xi = xmin+real(i,f64)*L/real(N,f64)
     u(i+1) = 1.0_f64 / (2.0_f64 + sin(mode*2._f64*sll_p_pi/L*xi))
     u_exact(i+1) =  1.0_f64 / (2.0_f64 + sin(mode*2._f64*sll_p_pi/L*(xi-alpha)))
!    u(i+1) = cos(mode*twopi*i/N)
!    u_exact(i+1) = cos(mode*twopi*(i-alpha)/N)
   end do
   print*, 'p=', p
   call interp_per%initialize( N+1, xmin, xmax, sll_p_lagrange, 12)
   interp => interp_per
   call interp%interpolate_array_disp(N+1, u(1:N+1), -alpha, u_out(1:N+1))
   old_error = error
   error = maxval(abs(u_out(1:N+1)-u_exact(1:N+1)))
   print *, "error =", error
   print*, 'N=',N, 'error=', error, 'numerical order=', log(old_error/error)/log(2.0_f64)
end do


  error = 1.0_8
  print*,""
  print*, 'Testing order of sll_p_lagrange interpolation (periodic)'
  ! loop on N 
  N = N0
  do p=1,4
     N= 2*N 
     
     ! Interpolate non trivial smooth periodic function
     do  i=0, N
        xi = xmin+real(i,f64)*L/real(N,f64)
        u(i+1) = 1.0_f64 / (2.0_f64 + sin(mode*2._f64*sll_p_pi/L*xi))
        u_exact(i+1) =  1.0_f64 / (2.0_f64 + sin(mode*2._f64*sll_p_pi/L*(xi-alpha)))
!        u(i+1) = cos(mode*twopi*i/N)
!        u_exact(i+1) = cos(mode*twopi*(i-alpha)/N)
     end do

     call interp_lagrange%initialize( N+1,xmin,xmax,sll_p_periodic,6)
     interp => interp_lagrange
     call interp%interpolate_array_disp(N+1, u(1:N+1), -alpha, u_out(1:N+1))

     old_error = error
     error = maxval(abs(u_out(1:N+1)-u_exact(1:N+1)))
     print*, 'N=',N, 'error=', error, 'numerical order=', log(old_error/error)/log(2.0_f64) 

  end do

print*, 'PASSED'
end program periodic_interpolation
