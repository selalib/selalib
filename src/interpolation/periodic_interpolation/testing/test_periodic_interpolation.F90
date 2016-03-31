program test_periodic_interpolation
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"

  use sll_m_constants, only: &
    sll_p_twopi

  use sll_m_periodic_interp, only: &
    sll_s_initialize_periodic_interp, &
    sll_s_periodic_interp, &
    sll_t_periodic_interp_work, &
    sll_p_spline

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  sll_int32, parameter    :: N0 = 16
  sll_real64               :: u(16*N0), u_exact(16*N0), u_out(16*N0)
  type(sll_t_periodic_interp_work), pointer :: interp
  sll_real64 :: alpha, error, old_error
  sll_int32 :: i, p, N
  !sll_int32 :: i0
  sll_int32 :: mode 

  error = 0.0_f64
  print*, 'Testing order of periodic interpolation'
  ! loop on N 
  N = N0
  do p=1,4
     N= 2*N 
     alpha = 0.05_f64
     
     ! Interpolate non trivial smooth periodic function
     mode = 3
     do  i=0, N-1
        u(i+1)       = 1.0_f64 / (2.0_f64 + sin(mode*sll_p_twopi*i/N))
        u_exact(i+1) = 1.0_f64 / (2.0_f64 + sin(mode*sll_p_twopi*(i-alpha)/N))
        !u(i+1) = cos(mode*twopi*i/N)
        !u_exact(i+1) = cos(mode*twopi*(i-alpha)/N)
     end do

     
     call sll_s_initialize_periodic_interp(interp, N, sll_p_spline, 8)
     !call sll_s_initialize_periodic_interp(interp, N, sll_p_trigo_fft_selalib, 8)
     !call sll_s_initialize_periodic_interp(interp, N, sll_p_trigo, 8)
     !call sll_s_initialize_periodic_interp(interp, N, sll_p_lagrange, 16)
     call sll_s_periodic_interp(interp, u_out,  u, alpha)
     
     old_error = error
     error = maxval(abs(u_out(1:N)-u_exact(1:N)))
     
     
     if (p>1) then
        print*, 'N=',N, 'error=', error, 'numerical order=', &
             log(old_error/error)/log(2.0_f64)
     endif
  end do
 
  print*, 'PASSED'
end program test_periodic_interpolation

