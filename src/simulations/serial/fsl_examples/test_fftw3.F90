program test_fftw3
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_fftw.h"

  use iso_c_binding, only: &
    c_double_complex, &
    c_f_pointer, &
    c_ptr, &
    c_size_t

#ifdef FFTW_F2003
  use sll_m_fftw3, only: &
    fftw_alloc_complex, &
    fftw_destroy_plan, &
    fftw_execute_dft_c2r, &
    fftw_execute_dft_r2c, &
    fftw_free, &
    fftw_measure, &
    fftw_plan_dft_c2r_2d, &
    fftw_plan_dft_r2c_2d
#endif

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#ifdef FFTW_F2003

type(C_PTR) :: fw, bw
complex(C_DOUBLE_COMPLEX), dimension(:,:), pointer :: rhot
integer(C_SIZE_T) :: sz_rhot
type(C_PTR) :: p_rhot
integer :: nc_x, nc_y
real(8), dimension(:,:), allocatable :: rho
integer :: i
real(8) :: t0, t1


call cpu_time(t0)
nc_x = 128
nc_y = 128
allocate(rho(nc_x,nc_y))

sz_rhot = int((nc_x/2+1)*nc_y,C_SIZE_T)
p_rhot = fftw_alloc_complex(sz_rhot)
call c_f_pointer(p_rhot, rhot, [nc_x/2+1,nc_y])

fw = fftw_plan_dft_r2c_2d(nc_y,nc_x,rho,rhot,FFTW_MEASURE)
bw = fftw_plan_dft_c2r_2d(nc_y,nc_x,rhot,rho,FFTW_MEASURE)

do i = 1, 10
  call fftw_execute_dft_r2c(fw, rho, rhot)
  call fftw_execute_dft_c2r(bw, rhot, rho)
end do

call fftw_free(p_rhot)
call fftw_destroy_plan(fw)
call fftw_destroy_plan(bw)

call cpu_time(t1)

write(*,"('CPU time =', g15.3)") t1-t0

#endif

end program test_fftw3
