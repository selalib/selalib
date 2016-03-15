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
sll_real64                                :: omega
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
sll_real64  :: pi = sll_p_pi

call cpu_time(tstart)

!w  =   sqrt(6)*pi
!ex =   cos(pi*x)*sin(pi*y)*sin(-2*pi*z)*cos(omega*t)
!ey =   sin(pi*x)*cos(pi*y)*sin(-2*pi*z)*cos(omega*t)
!ez =   sin(pi*x)*sin(pi*y)*cos(-2*pi*z)*cos(omega*t)
!hx = - sin(pi*x)*cos(pi*y)*cos(-2*pi*z)*sin(omega*t)*3*pi/omega
!hy =   cos(pi*x)*sin(pi*y)*cos(-2*pi*z)*sin(omega*t)*3*pi/omega
!hz = + 0

x_min = .0_f64; x_max = 1.0_f64
y_min = .0_f64; y_max = 1.0_f64
z_min = .0_f64; z_max = 1.0_f64

nc_x = 64; nc_y = 64; nc_z = 64

delta_x = (x_max-x_min)/nc_x
delta_y = (y_max-y_min)/nc_y
delta_z = (z_max-z_min)/nc_z

SLL_ALLOCATE(x(nc_x+1,nc_y+1,nc_z+1),error)
SLL_ALLOCATE(y(nc_x+1,nc_y+1,nc_z+1),error)
SLL_ALLOCATE(z(nc_x+1,nc_y+1,nc_z+1),error)

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
omega = sqrt(6.0_f64)*pi
nstep = int(2*pi / dt,i32)

print*, nstep*dt- 2*pi
nstep = nstep +1
dt    = 2.0_f64*pi / (nstep)
print*, nstep*dt- 2*pi

SLL_CLEAR_ALLOCATE(hx(1:nc_x+1,1:nc_y+1,1:nc_z+1), error)
SLL_CLEAR_ALLOCATE(hy(1:nc_x+1,1:nc_y+1,1:nc_z+1), error)
SLL_CLEAR_ALLOCATE(hz(1:nc_x+1,1:nc_y+1,1:nc_z+1), error)
SLL_CLEAR_ALLOCATE(ex(1:nc_x+1,1:nc_y+1,1:nc_z+1), error)
SLL_CLEAR_ALLOCATE(ey(1:nc_x+1,1:nc_y+1,1:nc_z+1), error)
SLL_CLEAR_ALLOCATE(ez(1:nc_x+1,1:nc_y+1,1:nc_z+1), error)

call sll_o_create(maxwell,            &
                  x_min, x_max, nc_x, &
                  y_min, y_max, nc_y, & 
                  z_min, z_max, nc_z  )

time = -0.5_f64*dt
hx   = - sin(pi*x)*cos(pi*y)*cos(-2.0_f64*pi*z)*sin(omega*time)*3.0_f64*pi/omega
hy   =   cos(pi*x)*sin(pi*y)*cos(-2.0_f64*pi*z)*sin(omega*time)*3.0_f64*pi/omega
hz   =   0.0_f64

time =  0.0_f64
ex   =   cos(pi*x)*sin(pi*y)*sin(-2.0_f64*pi*z)*cos(omega*time)
ey   =   sin(pi*x)*cos(pi*y)*sin(-2.0_f64*pi*z)*cos(omega*time)
ez   =   sin(pi*x)*sin(pi*y)*cos(-2.0_f64*pi*z)*cos(omega*time)

do istep = 1, nstep !*** Loop over time

  write(*,"(10x,' istep = ',I6)",advance="no") istep
  write(*,"(', time  = ',g12.3,' sec')",advance="no") time
  write(*,"(', error = ',g15.5)") &
  maxval(abs(hx+sin(pi*x)*cos(pi*y)*cos(-2*pi*z)*sin(omega*time)*3.0_f64*pi/omega)) 

  call plot_field( 'ex.dat', ex )
  call plot_field( 'ey.dat', ey )
  call plot_field( 'ez.dat', ez )
  call ampere(maxwell, ex, ey, ez, hx, hy, hz, dt) 
  time = time + 0.5_f64*dt
  call plot_field( 'hx.dat', hx )
  call plot_field( 'hy.dat', hy )
  call plot_field( 'hz.dat', hz )
  hx = - sin(pi*x)*cos(pi*y)*cos(-2.0_f64*pi*z)*sin(omega*time)*3.0_f64*pi/omega
  hy =   cos(pi*x)*sin(pi*y)*cos(-2.0_f64*pi*z)*sin(omega*time)*3.0_f64*pi/omega
  hz =   0.0_f64
  call plot_field( 'hx.ref', hx )
  call plot_field( 'hy.ref', hy )
  call plot_field( 'hz.ref', hz )
  stop

  call faraday(maxwell, hx, hy, hz, ex, ey, ez, dt)   
  time = time + 0.5_f64*dt

  ex   =   cos(pi*x)*sin(pi*y)*sin(-2.0_f64*pi*z)*cos(omega*time)
  ey   =   sin(pi*x)*cos(pi*y)*sin(-2.0_f64*pi*z)*cos(omega*time)
  ez   =   sin(pi*x)*sin(pi*y)*cos(-2.0_f64*pi*z)*cos(omega*time)
  call plot_field( 'ex.ref', ex )
  call plot_field( 'ey.ref', ey )
  call plot_field( 'ez.ref', ez )


end do ! next time step

call cpu_time(tend)
print"('CPU time : ',g15.3)", tend-tstart
print*,'PASSED'

deallocate(hx,hy,hz,ex,ey,ez)

contains

subroutine plot_field( filename, field )

  character(len=*) :: filename
  sll_real64       :: field(:,:,:)
  sll_int32        :: n1, n2, n3
  sll_int32        :: i, j

  n1 = size(field,1)
  n2 = size(field,2)
  n3 = size(field,3)

  open(10, file=filename)
  do i =1,n1
    do j=1,n2
      write(10,*) i, j, sngl(field(i,j,n3/2))
    end do
    write(10,*) 
  end do
  close(10)

end subroutine plot_field


end program test_maxwell_3d_pstd
