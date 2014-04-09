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
type(mudpack_2d) :: periodic
type(mudpack_2d) :: dirichlet
real(8), allocatable :: sol(:,:)
real(8), allocatable :: phi(:,:)
real(8), allocatable :: rhs(:,:)
real(8), allocatable :: eta1(:,:)
real(8), allocatable :: eta2(:,:)

real(8) :: eta1_min, eta1_max, eta2_min, eta2_max
real(8) :: delta_eta1, delta_eta2

integer :: i, j, error

nc_eta1 = 128
nc_eta2 = 128

SLL_CLEAR_ALLOCATE(eta1(1:nc_eta1+1,1:nc_eta2+1),error)
SLL_CLEAR_ALLOCATE(eta2(1:nc_eta1+1,1:nc_eta2+1),error)

!set end points of solution rectangle in (x,y) space
eta1_min = 0.0
eta1_max = 4.0
eta2_min = 0.0
eta2_max = 4.0

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

SLL_CLEAR_ALLOCATE(sol( 1:nc_eta1+1,1:nc_eta2+1),error)
SLL_CLEAR_ALLOCATE(phi( 1:nc_eta1+1,1:nc_eta2+1),error)
SLL_CLEAR_ALLOCATE(rhs( 1:nc_eta1+1,1:nc_eta2+1),error)


!Poisson periodic

call initialize_mudpack_cartesian(periodic,                    &
                                  eta1_min, eta1_max, nc_eta1, &
                                  eta2_min, eta2_max, nc_eta2, &
                                  SLL_PERIODIC, SLL_PERIODIC,  &
                                  SLL_PERIODIC, SLL_PERIODIC)


sol  = sin(2*sll_pi*eta1)*sin(2*sll_pi*eta2)
rhs  = -8*sll_pi**2 * sol + 1.

call solve_mudpack_cartesian(periodic, phi, rhs)

call sll_gnuplot_corect_2d(eta1_min, eta1_max, nc_eta1+1, &
                           eta2_min, eta2_max, nc_eta2+1, &
                           phi, "sinsin", 1, error)

!compute and print maximum norm of error
write(*,201) maxval(abs(phi-sol))

print*,"PASSED"
     
call initialize_mudpack_cartesian(dirichlet,                    &
                                  eta1_min, eta1_max, nc_eta1,  &
                                  eta2_min, eta2_max, nc_eta2,  &
                                  SLL_DIRICHLET, SLL_DIRICHLET, &
                                  SLL_DIRICHLET, SLL_DIRICHLET)


sol = sin(2*sll_pi*eta1)*cos(2*sll_pi*eta2)
rhs = -8*sll_pi**2 * sol
phi = sol

call solve_mudpack_cartesian(dirichlet, phi, rhs)

call sll_gnuplot_corect_2d(eta1_min, eta1_max, nc_eta1+1, &
                           eta2_min, eta2_max, nc_eta2+1, &
                           phi, "sincos", 2, error)

!compute and print maximum norm of error
write(*,201) maxval(abs(phi-sol))

print*,"PASSED"

201 format(' maximum error  =  ',e10.3)
     
end program test_mudpack_cartesian
