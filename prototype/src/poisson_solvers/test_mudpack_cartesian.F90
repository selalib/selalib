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
#include "sll_memory.h"

use sll_mudpack_cartesian

implicit none

integer :: nc_eta1
integer :: nc_eta2
type(mudpack_2d) :: poisson
real(8), allocatable :: phi(:,:)
real(8), allocatable :: rhs(:,:)
real(8), allocatable :: eta1(:,:)
real(8), allocatable :: eta2(:,:)

real(8) :: eta1_min, eta1_max, eta2_min, eta2_max
real(8) :: delta_eta1, delta_eta2

integer :: i, j, error, istep
real, parameter :: mode_1 = 1.0
real, parameter :: mode_2 = 1.0
real(8) :: avg, dt = 0.1

nc_eta1 = 128
nc_eta2 = 128

SLL_CLEAR_ALLOCATE(eta1(1:nc_eta1+1,1:nc_eta2+1),error)
SLL_CLEAR_ALLOCATE(eta2(1:nc_eta1+1,1:nc_eta2+1),error)

!set end points of solution rectangle in (x,y) space
eta1_min = 0.0
eta1_max = 4.*sll_pi
eta2_min = 0.0
eta2_max = 4.*sll_pi

!set mesh increments
delta_eta1 = (eta1_max-eta1_min)/float(nc_eta1)
delta_eta2 = (eta2_max-eta2_min)/float(nc_eta2)

!set right hand side in rhs and initialize phi to zero
do i=1,nc_eta1+1
   do j=1,nc_eta2+1
      eta1(i,j) = eta1_min+float(i-1)*delta_eta1
      eta2(i,j) = eta2_min+float(j-1)*delta_eta2
   end do
end do

call initialize_mudpack_cartesian(poisson,                     &
                                  eta1_min, eta1_max, nc_eta1, &
                                  eta2_min, eta2_max, nc_eta2, &
                                  SLL_PERIODIC, SLL_PERIODIC,  &
                                  SLL_PERIODIC, sll_PERIODIC)

SLL_CLEAR_ALLOCATE(phi( 1:nc_eta1+1,1:nc_eta2+1),error)
SLL_CLEAR_ALLOCATE(rhs( 1:nc_eta1+1,1:nc_eta2+1),error)


do istep = 0, 9

   rhs = cos(mode_1*eta1)*cos(mode_2*eta2)
   rhs  = -2*mode_1**3 * sin(mode_1*eta1)*cos(mode_2*eta2)*cos(sll_pi*istep*dt)

   call sll_gnuplot_corect_2d(eta1_min, eta1_max, nc_eta1+1, &
                              eta2_min, eta2_max, nc_eta2+1, &
                              rhs, "rhs", istep+1, error)

   
   call solve_mudpack_cartesian(poisson, phi, rhs)

   avg  = sum(phi)/((nc_eta1+1)*(nc_eta2+1))
   rhs = rhs - phi
   
   call sll_gnuplot_corect_2d(eta1_min, eta1_max, nc_eta1+1, &
                              eta2_min, eta2_max, nc_eta2+1, &
                              phi, "phi", istep+1, error)

end do

!compute and print maximum norm of error
!write(*,201) maxval(abs((phi-mode*sin(mode*eta1)*cos(mode*eta2))))

print*,"PASSED"

!201 format(' maximum error  =  ',e10.3)
     
end program test_mudpack_cartesian
