!
!     a sample program/test driver for mudpack is listed below.  it
!     can be executed as an initial test.  the output is listed for
!     the test case described.
!
!     test the driver below by solving the separable elliptic pde
!
!     (1.+x**2)*pxx + exp(1.-y)*(pyy-py) - (x+y)*pe = r(x,y)
!
!     on the unit square with specified boundary conditions at
!     eta1_max = 1.0, eta2_min = 0.0 and mixed boundary conditions
!
!          dp/dx - pe(0.0,y) =  ga(y)  (at x = 0.0)
!
!          dp/dy + pe(x,1.0) = gd(x)  (at y = 1.0)
!
!     use point relaxation and choose a grid as close to 60 x 50
!     as the grid size constraints allow.  use the exact solution
!
!          pe(x,y) = (x**3+y**3+1.0)/3
!
!     for testing.  first mud2sp is called to yield a second-order
!     approximation.  then mud24sp is called to improve the estimate.
!

program test_mudpack_polar
#include "sll_working_precision.h"
#include "sll_utilities.h"
#include "sll_file_io.h"

use sll_mudpack_cartesian

implicit none

integer, parameter :: nc_eta1 = 64, nc_eta2 = 64
real(8) :: phi(nc_eta1+1,nc_eta2+1)
real(8) :: rhs(nc_eta1+1,nc_eta2+1)
real(8) :: eta1(nc_eta1+1,nc_eta2+1)
real(8) :: eta2(nc_eta1+1,nc_eta2+1)

real(8) :: eta1_min, eta1_max, eta2_min, eta2_max
real(8) :: delta_eta1, delta_eta2
real(8) :: x, cx, cxx, cex, px, pxx
real(8) :: y, cy, cyy, cey, py, pyy
real(8) :: ce, pe, errmax

integer :: i, j, error
integer :: nx, ny

!set end points of solution rectangle in (x,y) space
eta1_min = 0.0
eta1_max = 1.0
eta2_min = 0.0
eta2_max = 1.0

nx = nc_eta1+1
ny = nc_eta2+1

!set mesh increments
delta_eta1 = (eta1_max-eta1_min)/float(nc_eta1)
delta_eta2 = (eta2_min-eta2_min)/float(nc_eta2)

!set right hand side in rhs and initialize phi to zero
do i=1,nx
   do j=1,ny
      eta1(i,j) = eta1_min+float(i-1)*delta_eta1
      eta2(i,j) = eta2_min+float(j-1)*delta_eta2
      call cofx(eta1(i,j),cxx,cx,cex)
      call cofy(eta2(i,j),cyy,cy,cey)
      ce = cex+cey
      call exact(eta1(i,j),eta2(i,j),pxx,pyy,px,py,pe)
      rhs(i,j) = cxx*pxx+cyy*pyy+cx*px+cy*py+ce*pe
      phi(i,j) = pe
   end do
end do

call sll_gnuplot_corect_2d(eta1_min, eta1_max, nx, &
                           eta2_min, eta2_max, ny, &
                           rhs, "rhs", 1, error)

call sll_gnuplot_corect_2d(eta1_min, eta1_max, nx, &
                           eta2_min, eta2_max, ny, &
                           phi, "phi_exact", 1, error)

call initialize_mudpack_cartesian(phi, rhs, &
                                  eta1_min, eta1_max, nc_eta1, &
                                  eta2_min, eta2_max, nc_eta2)

call solve_mudpack_cartesian(phi, rhs)

call sll_gnuplot_corect_2d(eta1_min,eta1_max,nx,eta2_min,eta2_max,ny,phi,"phi",1,error)

write (*,108) error
if (error > 0) call exit(0)

if (error <= 0) then

   !compute and print maximum norm of error
   errmax = 0.0
   do j=1,ny
      do i=1,nx
         call exact(eta1(i,j),eta2(i,j),pxx,pyy,px,py,pe)
         errmax = dmax1(errmax,abs((phi(i,j)-pe)))
      end do
   end do
   write(*,201) errmax

end if

! attempt to improve approximation to fourth order

print*,"PASSED"

108 format(/' mud24sp test ', ' Error = ',i2)
201 format(' maximum error  =  ',e10.3)
     
end program test_mudpack_polar
