program test_periodic_interp
#include "sll_working_precision.h"
#ifndef STDF95
  use sll_module_interpolators_1d_base
#endif
  !use util_constants
  use numeric_constants
  use sll_periodic_interpolator_1d
  use sll_lagrange_interpolator_1d
  implicit none 
  
  type(per_1d_interpolator), target       :: interp_per
  type(lagrange_1d_interpolator), target       :: interp_lagrange
  sll_int32, parameter    :: N0 = 16
  sll_real64               :: u(16*N0+1), u_exact(16*N0+1), u_out(16*N0+1)
#ifdef STDF95
  type(per_1d_interpolator), pointer     :: interp
#else
  class(sll_interpolator_1d_base), pointer     :: interp
#endif
  sll_real64, parameter :: xmin = 0., xmax=3.   
 ! sll_real64, parameter :: xmin = 0., xmax=2.   
  sll_real64 :: alpha, error, old_error, L, xi
  sll_int32 :: i, p, N, i0, mode ,j

  alpha = 0.05_8
 ! alpha = 0.1_8
 ! alpha = 0.0001_8

  L = xmax - xmin
  error = 0.0_8
  print*, 'Testing order of periodic interpolation'
  ! loop on N 
  N = N0
  do p=1,4
!!!rajout raf
u_exact=0.0_f64
u_out=0.0_f64
     N= 2*N 
     
     ! Interpolate non trivial smooth periodic function
     mode = 3
     do  i=0, N
        xi = xmin+i*L/N
        u(i+1) = 1.0_8 / (2 + sin(mode*2._f64*sll_pi/L*xi))
        u_exact(i+1) =  1.0_8 / (2 + sin(mode*2._f64*sll_pi/L*(xi-alpha)))
!        u(i+1) = cos(mode*twopi*i/N)
!        u_exact(i+1) = cos(mode*twopi*(i-alpha)/N)
!        u(i+1)=2.0_f64*xi*xi+1.0_f64
!        u_exact(i+1)=2.0_f64*(xi-alpha)*(xi-alpha)+1.0_f64
!        u(i+1)= 2.0_f64*(sin(xi) + 2.5_f64 + cos(xi))
!        u_exact(i+1)=2.0_f64*(sin(xi-alpha) + 2.5_f64 + cos(xi-alpha))
     end do
#ifdef STDF95
     call periodic_initialize(interp, N+1, xmin, xmax, SPLINE, 12)
#else
     call interp_per%initialize( N+1, xmin, xmax, SPLINE, 12)
#endif
     interp => interp_per
#ifdef STDF95
     u_out(1:N)=per_interpolate_array_at_displacement(interp, N+1, u(1:N+1), alpha)
#else
     u_out(1:N+1)=interp%interpolate_array_disp(N+1, u(1:N+1), alpha)
#endif
!end periodic_interp

     
     old_error = error
     !error = maxval(abs(u_out(1:N+1)-u_exact(1:N+1)))
     error = maxval(abs((u_out(1:N+1)-u_exact(1:N+1))/u_exact(1:N+1)))
     print*, 'N=',N, 'error=', error, 'numerical order=', log(old_error/error)/log(2.0_8) 
  end do



  error = 0.0_8
  print*, 'Testing order of lagrange interpolation'
  ! loop on N 
  N = N0
  do p=1,4
     N= 2*N 
     
     ! Interpolate non trivial smooth periodic function
     mode = 3
     do  i=0, N
        xi = xmin+i*L/N
        u(i+1) = 1.0_8 / (2 + sin(mode*2._f64*sll_pi/L*xi))
        u_exact(i+1) =  1.0_8 / (2 + sin(mode*2._f64*sll_pi/L*(xi-alpha)))
!        u(i+1) = cos(mode*twopi*i/N)
!        u_exact(i+1) = cos(mode*twopi*(i-alpha)/N)
!        u(i+1)=2.0_f64*xi*xi+1.0_f64
!        u_exact(i+1)=2.0_f64*(xi-alpha)*(xi-alpha)+1.0_f64
!        u(i+1)= 2.0_f64*(sin(xi) + 2.5_f64 + cos(xi))
!        u_exact(i+1)=2.0_f64*(sin(xi-alpha) + 2.5_f64 + cos(xi-alpha))
     end do

#ifdef STDF95
     call interp_lagrange_intialize(interp, N+1, xmin, xmax, PERIODIC_LAGRANGE,6)
#else
     call interp_lagrange%initialize( N+1,xmin,xmax,PERIODIC_LAGRANGE,6)
#endif
     interp => interp_lagrange
     u_out(1:N+1)=interp%interpolate_array_disp(N+1, u(1:N+1), -alpha)

     old_error = error
     !error = maxval(abs(u_out(1:N+1)-u_exact(1:N+1)))
     error = maxval(abs((u_out(1:N+1)-u_exact(1:N+1))/u_exact(1:N+1)))
     print*, 'N=',N, 'error=', error, 'numerical order=', log(old_error/error)/log(2.0_8) 
  end do


end program test_periodic_interp

