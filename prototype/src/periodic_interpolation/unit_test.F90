program test_periodic_interp
#include "sll_working_precision.h"
  use periodic_interp_module
  use sll_fft
  implicit none 
  
  sll_int32, parameter    :: N0 = 16
  sll_real64               :: u(16*N0), u_exact(16*N0), u_out(16*N0)
  type(periodic_interp_work), pointer :: interp
  sll_real64 :: alpha, error, old_error
  sll_int32 :: i, p, N, i0, mode 

  error = 0.0_8
  print*, 'Testing order of periodic interpolation'
  ! loop on N 
  N = N0
  do p=1,4
     N= 2*N 
     alpha = 0.05_8
     
     ! Interpolate non trivial smooth periodic function
     mode = 3
     do  i=0, N-1
        u(i+1) = 1.0_8 / (2 + sin(mode*twopi*i/N))
        u_exact(i+1) =  1.0_8 / (2 + sin(mode*twopi*(i-alpha)/N))
        !u(i+1) = cos(mode*twopi*i/N)
        !u_exact(i+1) = cos(mode*twopi*(i-alpha)/N)
     end do

     
     call initialize_periodic_interp(interp, N, SPLINE, 8)
     !call initialize_periodic_interp(interp, N, TRIGO_FFT_SELALIB, 8)
     !call initialize_periodic_interp(interp, N, TRIGO, 8)
     !call initialize_periodic_interp(interp, N, LAGRANGE, 16)
     call periodic_interp(interp, u_out,  u, alpha)
     
     old_error = error
     error = maxval(abs(u_out(1:N)-u_exact(1:N)))
     
     
     if (p>1) then
        print*, 'N=',N, 'error=', error, 'numerical order=', &
             log(old_error/error)/log(2.0_8) 
     endif
  end do
 

end program test_periodic_interp

