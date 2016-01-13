program test_pic2d
#include "sll_working_precision.h"
#include "sll_memory.h"

use zone
use particules
use diagno

use sll_m_poisson_2d_periodic_fftpack

implicit none

type(tm_mesh_fields) :: f
type(particle)       :: p

sll_real64 :: time
sll_real64 :: xmin
sll_real64 :: xmax
sll_real64 :: ymin
sll_real64 :: ymax

sll_int32  :: istep
sll_int32  :: iplot
sll_int32  :: iargc
sll_int32  :: n
sll_int32  :: i
sll_int32  :: j
sll_int32  :: error

sll_real64 :: aux1, aux2

character(len=272) :: argv
type(sll_t_poisson_2d_periodic_fftpack) :: poisson

n = iargc()
if (n == 0) stop 'Usage: ./bin/test_pic2d fichier-de-donnees.nml'
do i = 1, n
  call getarg( i, argv)
  write(*,'(i2, 1x, a)') i, argv
end do

call readin( trim(argv) )

SLL_CLEAR_ALLOCATE(f%ex(0:nx,0:ny),  error)
SLL_CLEAR_ALLOCATE(f%ey(0:nx,0:ny),  error)
SLL_CLEAR_ALLOCATE(f%bz(0:nx,0:ny), error)
SLL_CLEAR_ALLOCATE(f%jx(0:nx,0:ny),  error)
SLL_CLEAR_ALLOCATE(f%jy(0:nx,0:ny),  error)
SLL_CLEAR_ALLOCATE(f%r0(0:nx,0:ny),    error) !rho au temps n
SLL_CLEAR_ALLOCATE(f%r1(0:nx,0:ny),    error) !rho au temps n+1

time  = 0.d0
iplot = 0

if( nstep > nstepmax ) nstep = nstepmax

!********************************************************************

istep = 1

do i=0,nx-1
  aux1 = alpha/kx * sin(kx*x(i))
  aux2 = alpha * cos(kx*x(i))
  do j=0,ny
    f%ex(i,j) = aux1
    f%r1(i,j) = aux2
  enddo
enddo
      
xmin = 0.0_f64; xmax = dimx
ymin = 0.0_f64; ymax = dimy

call plasma( p, time ) 

call calcul_rho( p, f )

call sll_o_gnuplot_2d(xmin, xmax, nx+1, &
                      ymin, ymax, ny+1, &
                      f%r0, 'rho', 1, error)

call sll_o_initialize( poisson, xmin, xmax, nx, &
                       ymin, ymax, ny, error) 

call sll_o_solve( poisson, f%ex, f%ey, f%r0)

call sll_o_gnuplot_2d(xmin, xmax, nx+1, &
                      ymin, ymax, ny+1, &
                      f%ex, 'ex', 1, error)

do istep = 1, nstep

  call interpol_eb( f, p )
  call avancee_vitesse( p )
  call avancee_part( p, 1.d0 )
  call sortie_part( p )
  call calcul_rho( p, f )
  call sll_o_solve( poisson, f%ex, f%ey, f%r0)

  time = time + dt

  call modeE( f, istep, time )
  write(*,"('istep = ', i6, ' time = ')", advance='no') istep

end do

print*,'PASSED'

end program test_pic2d
