!
!  Contact : Pierre Navaro http://wwww-irma.u-strasbg.fr/~navaro
!
program test_maxwell_3d_pstd
!-------------------------------------------------------------------
!  test 2D Maxwell solver based on FFT
!-------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_maxwell_solvers_macros.h"

use sll_m_constants, only: sll_p_pi
use sll_m_maxwell_3d_pstd
use sll_m_utilities

implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

sll_real64 :: x_max, x_min
sll_real64 :: y_max, y_min
sll_real64 :: z_max, z_min
sll_real64 :: delta_x, delta_y, delta_z

sll_int32  :: nc_x, nc_y, nc_z
sll_int32  :: error

type(sll_t_maxwell_3d_pstd)               :: maxwell
sll_int32                                 :: i, j, k
sll_real64                                :: w1, w2
sll_real64                                :: time
sll_int32                                 :: istep, nstep
sll_real64                                :: dt
sll_real64                                :: cfl = 0.5_f64
sll_real64, dimension(:,:,:), allocatable :: x
sll_real64, dimension(:,:,:), allocatable :: y
sll_real64, dimension(:,:,:), allocatable :: z
sll_real64, dimension(:,:,:), allocatable :: hx
sll_real64, dimension(:,:,:), allocatable :: hy
sll_real64, dimension(:,:,:), allocatable :: hz
sll_real64, dimension(:,:,:), allocatable :: ex
sll_real64, dimension(:,:,:), allocatable :: ey
sll_real64, dimension(:,:,:), allocatable :: ez
sll_real64 :: tstart, tend
sll_real64 :: pi = sll_p_pi

call cpu_time(tstart)


x_min = .0_f64; x_max = 1.0_f64
y_min = .0_f64; y_max = 1.0_f64
z_min = .0_f64; z_max = 1.0_f64

nc_x = 64; nc_y = 64; nc_z = 64

delta_x = (x_max-x_min)/nc_x
delta_y = (y_max-y_min)/nc_y
delta_z = (z_max-z_min)/nc_z

SLL_CLEAR_ALLOCATE(x(1:nc_x+1,1:nc_y+1,1:nc_z+1),error)
SLL_CLEAR_ALLOCATE(y(1:nc_x+1,1:nc_y+1,1:nc_z+1),error)
SLL_CLEAR_ALLOCATE(z(1:nc_x+1,1:nc_y+1,1:nc_z+1),error)

do k = 1, nc_z+1
  do j = 1, nc_y+1
    do i = 1, nc_x+1
      x(i,j,k) = x_min + (i-1)*delta_x
      y(i,j,k) = y_min + (j-1)*delta_y
      z(i,j,k) = z_min + (k-1)*delta_z
    end do
  end do
end do

dt    = cfl/sqrt(1./delta_x**2+1./delta_y**2+1./delta_z**2)
nstep = 10

SLL_CLEAR_ALLOCATE(hx(1:nc_x+1,1:nc_y+1,1:nc_z+1), error)
SLL_CLEAR_ALLOCATE(hy(1:nc_x+1,1:nc_y+1,1:nc_z+1), error)
SLL_CLEAR_ALLOCATE(hz(1:nc_x+1,1:nc_y+1,1:nc_z+1), error)
SLL_CLEAR_ALLOCATE(ex(1:nc_x+1,1:nc_y+1,1:nc_z+1), error)
SLL_CLEAR_ALLOCATE(ey(1:nc_x+1,1:nc_y+1,1:nc_z+1), error)
SLL_CLEAR_ALLOCATE(ez(1:nc_x+1,1:nc_y+1,1:nc_z+1), error)

call sll_s_maxwell_3d_pstd_init(maxwell,            &
                                x_min, x_max, nc_x, &
                                y_min, y_max, nc_y, & 
                                z_min, z_max, nc_z  )

w1 = + 2.0_f64*sqrt(3.0_f64)*pi
w2 = - sqrt(3.0_f64)

time = -0.5_f64*dt

hx = + w2 * cos(2*pi*(x+y+z)-w1*time)
hy = 0.0_f64
hz = - w2 * cos(2*pi*(x+y+z)-w1*time)

time = 0.0_f64

ex = cos(2*pi*(x+y+z)-w1*time)
ey = -2* ex
ez = cos(2*pi*(x+y+z)-w1*time)

do istep = 1, nstep !*** Loop over time

  write(*,"(10x,' istep = ',I6)",advance="no") istep
  write(*,"(', time  = ',g12.3,' sec')") time

  call sll_s_maxwell_3d_pstd_faraday(maxwell, ex, ey, ez, hx, hy, hz, dt)   
  time = time + 0.5*dt
  !call plot_field("hx",hx, w2*cos(2*pi*(x+y+z)-w1*time))
  !call plot_field("hz",hz,-w2*cos(2*pi*(x+y+z)-w1*time))
  call sll_s_maxwell_3d_pstd_ampere(maxwell, hx, hy, hz, ex, ey, ez, dt) 
  time = time + 0.5*dt
  !call plot_field("ex",ex,cos(2*pi*(x+y+z)-w1*time))
  !call plot_field("ez",ez,cos(2*pi*(x+y+z)-w1*time))

end do ! next time step

print*, "ex error :", maxval(abs(ex-cos(2*pi*(x+y+z)-w1*time)))
print*, "ey error :", maxval(abs(ey+2.*cos(2*pi*(x+y+z)-w1*time)))
print*, "ez error :", maxval(abs(ex-cos(2*pi*(x+y+z)-w1*time)))
time = time - 0.5_f64*dt
print*, "hx error :", maxval(abs(hx-w2*cos(2*pi*(x+y+z)-w1*time)))
print*, "hy error :", maxval(abs(hy))
print*, "hz error :", maxval(abs(hz+w2*cos(2*pi*(x+y+z)-w1*time)))



call cpu_time(tend)
print"('CPU time : ',g15.3)", tend-tstart
print*,'PASSED'

deallocate(hx,hy,hz,ex,ey,ez)

call sll_s_maxwell_3d_pstd_free(maxwell)

contains

subroutine plot_field( prefix, field1, field2 )

  character(len=*) :: prefix
  sll_real64       :: field1(:,:,:)
  sll_real64       :: field2(:,:,:)
  sll_int32        :: n1, n2, n3
  sll_int32        :: i, j
  character(len=4) :: cstep

  n1 = size(field1,1)
  n2 = size(field1,2)
  n3 = size(field1,3)

  call sll_s_int2string(istep, cstep)
  open(11, file=prefix//".gnu", position="append")
  if (istep == 1) rewind(11)
  write(11,"(a)") "splot '"//prefix//"_"//cstep//".dat' w l, &
                  & '"//prefix//"_"//cstep//".dat' u 1:2:4 w l"
  close(11)
  open(10, file=prefix//"_"//cstep//".dat")
  do i =1,n1
    do j=1,n2
      write(10,*) (i-1)*delta_x, (j-1)*delta_y, &
                  sngl(field1(i,j,n3/2)), sngl(field2(i,j,n3/2))
    end do
    write(10,*) 
  end do
  close(10)

end subroutine plot_field


end program test_maxwell_3d_pstd
