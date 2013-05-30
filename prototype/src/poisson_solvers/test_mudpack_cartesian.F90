!
!     a sample program/test driver for mudpack is listed below.  it
!     can be executed as an initial test.  the output is listed for
!     the test case described.
!
!     test the driver below by solving the separable elliptic pde
!
!     pxx + pyy = r(x,y)
!
!     on the square [0:2pi]x[0:2pi] with periodic boundary conditions
!

program test_mudpack_cartesian
#include "sll_working_precision.h"
#include "sll_utilities.h"
#include "sll_file_io.h"
#include "sll_constants.h"

use sll_mudpack_cartesian

implicit none

integer, parameter :: nc_eta1 = 64, nc_eta2 = 64
real(8) :: phi(nc_eta1+1,nc_eta2+1)
real(8) :: rhs(nc_eta1+1,nc_eta2+1)
real(8) :: eta1(nc_eta1+1,nc_eta2+1)
real(8) :: eta2(nc_eta1+1,nc_eta2+1)

real(8) :: eta1_min, eta1_max, eta2_min, eta2_max
real(8) :: delta_eta1, delta_eta2
real(8) :: pe, errmax

integer :: i, j, error
integer :: nx, ny, mode = 2

!set end points of solution rectangle in (x,y) space
eta1_min = 0.0
eta1_max = 2.*sll_pi
eta2_min = 0.0
eta2_max = 2.*sll_pi

nx = nc_eta1+1
ny = nc_eta2+1

!set mesh increments
delta_eta1 = (eta1_max-eta1_min)/float(nc_eta1)
delta_eta2 = (eta2_max-eta2_min)/float(nc_eta2)

!set right hand side in rhs and initialize phi to zero
do i=1,nx
   do j=1,ny
      eta1(i,j) = eta1_min+float(i-1)*delta_eta1
      eta2(i,j) = eta2_min+float(j-1)*delta_eta2
      rhs(i,j) = -2_f64 * mode**3 * sin(mode*eta1(i,j))*cos(mode*eta2(i,j))
      phi(i,j) = 0.0
   end do
end do

call sll_gnuplot_corect_2d(eta1_min, eta1_max, nx, &
                           eta2_min, eta2_max, ny, &
                           rhs, "rhs", 1, error)

call initialize_mudpack_cartesian(phi, rhs, &
                                  eta1_min, eta1_max, nc_eta1, &
                                  eta2_min, eta2_max, nc_eta2, &
                                  PERIODIC, PERIODIC, PERIODIC, PERIODIC)

call solve_mudpack_cartesian(phi, rhs)

call sll_gnuplot_corect_2d(eta1_min, eta1_max, nx, &
                           eta2_min, eta2_max, ny, &
                           phi, "phi", 1, error)


write (*,108) error
if (error > 0) call exit(0)

if (error <= 0) then

   !compute and print maximum norm of error
   errmax = 0.0
   do j=1,ny
      do i=1,nx
         call exact(eta1(i,j),eta2(i,j),mode,pe)
         errmax = dmax1(errmax,abs((phi(i,j)-pe)))
      end do
   end do
   write(*,201) errmax

end if

! attempt to improve approximation to fourth order

print*,"PASSED"

108 format(/'mudpack test ', ' Error = ',i2)
201 format(' maximum error  =  ',e10.3)
     
end program test_mudpack_cartesian

!> set an exact solution for testing purpose
subroutine exact(x,y,mode,pe)
implicit none
real(8)  :: x,y,pe
integer  :: mode
pe  = mode * sin(mode*x) * cos(mode*y)
return
end
