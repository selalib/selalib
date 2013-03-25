program test_periodic_interp
#include "sll_working_precision.h"
#ifndef STDF95
  use sll_module_interpolators_1d_base
#endif
  !use util_constants
  use numeric_constants
  use sll_periodic_interpolator_1d
  use sll_lagrange_interpolator_1d
!!!!ajout
  use sll_lagrange_interpolation2
  implicit none 
  
  type(per_1d_interpolator), target       :: interp_per
  type(lagrange_1d_interpolator), target       :: interp_lagrange
!!!ajout
  type( sll_lagrange_interpolation_1D2), pointer       :: l_i2
  !sll_int32, parameter    :: N0 = 16
  sll_int32, parameter    :: N0 = 4
  sll_real64               :: u(16*N0+1), u_exact(16*N0+1), u_out(16*N0+1)
#ifdef STDF95
  type(per_1d_interpolator), pointer     :: interp
#else
  class(sll_interpolator_1d_base), pointer     :: interp
#endif
  sll_real64, parameter :: xmin = 0.0_f64, xmax=3.0_f64   
  sll_real64 :: alpha, error, old_error, L, xi,mode
  sll_int32 :: i, p, N, i0, j

  alpha = 0.05_f64
  mode = 3.0_f64
  L = xmax - xmin

  error = 0.0_f64
  print*, 'Testing order of periodic interpolation'
  ! loop on N 
  N = N0
  do p=1,4
     N= 2*N 

     ! Interpolate non trivial smooth periodic function
     do  i=0, N
        xi = xmin+real(i,f64)*L/real(N,f64)
        u(i+1) = 1.0_f64 / (2.0_f64 + sin(mode*2._f64*sll_pi/L*xi))
        u_exact(i+1) =  1.0_f64 / (2.0_f64 + sin(mode*2._f64*sll_pi/L*(xi-alpha)))
!        u(i+1) = cos(mode*twopi*i/N)
!        u_exact(i+1) = cos(mode*twopi*(i-alpha)/N)
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
     
     old_error = error
     error = maxval(abs(u_out(1:N+1)-u_exact(1:N+1)))
     print*, 'N=',N, 'error=', error, 'numerical order=', log(old_error/error)/log(2.0_f64)
  end do



  error = 0.0_8
  print*,""
  print*, 'Testing order of lagrange interpolation'
  ! loop on N 
  N = N0
  do p=1,1!4
     N= 2*N 
     
     ! Interpolate non trivial smooth periodic function
     do  i=0, N
        xi = xmin+real(i,f64)*L/real(N,f64)
        u(i+1) = 1.0_f64 / (2.0_f64 + sin(mode*2._f64*sll_pi/L*xi))
        u_exact(i+1) =  1.0_f64 / (2.0_f64 + sin(mode*2._f64*sll_pi/L*(xi-alpha)))
!        u(i+1) = cos(mode*twopi*i/N)
!        u_exact(i+1) = cos(mode*twopi*(i-alpha)/N)
     end do

!!!modification
     call interp_lagrange%initialize( N+1,xmin,xmax,0,6)!PERIODIC_LAGRANGE,6)
     interp => interp_lagrange
     u_out(1:N+1)=interp%interpolate_array_disp(N+1, u(1:N+1), -alpha)

     old_error = error
     error = maxval(abs(u_out(1:N+1)-u_exact(1:N+1)))
     print*, 'N=',N, 'error=', error, 'numerical order=', log(old_error/error)/log(2.0_f64) 

if(p==1)then
print*,"u-ue"
print*,u_out(1:N+1)-u_exact(1:N+1)
end if
  end do


!!!!!!!!!!!!partie test avec lagrange2
  error = 0.0_8
  print*,""
  print*, 'Testing order of lagrange interpolation TESTTTT'
  ! loop on N 
  N = N0
  do p=1,1!4
     N= 2*N 
     
     ! Interpolate non trivial smooth periodic function
     do  i=0, N
        xi = xmin+real(i,f64)*L/real(N,f64)
        u(i+1) = 1.0_f64 / (2.0_f64 + sin(mode*2._f64*sll_pi/L*xi))
        u_exact(i+1) =  1.0_f64 / (2.0_f64 + sin(mode*2._f64*sll_pi/L*(xi-alpha)))
!        u(i+1) = cos(mode*twopi*i/N)
!        u_exact(i+1) = cos(mode*twopi*(i-alpha)/N)
     end do

	l_i2 =>new_lagrange_interpolation_1D2(N+1,xmin,xmax,0,6)
	call compute_lagrange_interpolation_1D2(u(1:N+1),-alpha,L_i2)
	call interpolate_array_values2(l_i2)
	u_out=l_i2%data_out

     old_error = error
     error = maxval(abs(u_out(1:N+1)-u_exact(1:N+1)))
     print*, 'N=',N, 'error=', error, 'numerical order=', log(old_error/error)/log(2.0_f64) 

if(p==1)then
print*,"u-ue"
print*,u_out(1:N+1)-u_exact(1:N+1)
end if
  end do


end program test_periodic_interp
