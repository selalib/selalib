program test_periodic_interp
#include "sll_working_precision.h"
  use sll_periodic_interpolator_1d
#ifndef STDF95
  use sll_module_interpolators_1d_base
#endif
  implicit none 
  
  sll_int32, parameter    :: N0 = 16
  sll_real64               :: u(16*N0), u_exact(16*N0), u_out(16*N0)
  type(periodic_1d_interpolator), target       :: interp_per
#ifdef STDF95
  type(periodic_1d_interpolator), pointer     :: interp
#else
  class(sll_interpolator_1d_base), pointer     :: interp
#endif
  sll_real64, parameter :: xmin = 1., xmax=3.   
  sll_real64 :: alpha, error, old_error, L, xi
  sll_int32 :: i, p, N, i0, mode 

  error = 0.0_8
  L = xmax - xmin
  print*, 'Testing order of periodic interpolation'
  ! loop on N 
  N = N0
  do p=1,4
     N= 2*N 
     alpha = 0.05_8
     
     ! Interpolate non trivial smooth periodic function
     mode = 3
     do  i=0, N-1
        xi = xmin+i*L/N
        u(i+1) = 1.0_8 / (2 + sin(mode*twopi/L*xi))
        u_exact(i+1) =  1.0_8 / (2 + sin(mode*twopi/L*(xi-alpha)))
        !u(i+1) = cos(mode*twopi*i/N)
        !u_exact(i+1) = cos(mode*twopi*(i-alpha)/N)
     end do
#ifdef STDF95
     call periodic_initialize(interp, N, xmin, xmax, SPLINE, 12)
#else
     call interp_per%initialize( N, xmin, xmax, SPLINE, 12)
#endif
     interp => interp_per
#ifdef STDF95
     u_out(1:N)=per_interpolate_array_at_displacement(interp, N, u, alpha)
#else
     u_out(1:N)=interp%interpolate_array_disp(N, u, alpha)
#endif
     
     old_error = error
     error = maxval(abs(u_out(1:N)-u_exact(1:N)))
     print*, 'N=',N, 'error=', error, 'numerical order=', log(old_error/error)/log(2.0_8) 
  end do
 

end program test_periodic_interp

