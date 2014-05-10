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

nc_eta1 = 64
nc_eta2 = 64

SLL_CLEAR_ALLOCATE(eta1(1:nc_eta1+1,1:nc_eta2+1),error)
SLL_CLEAR_ALLOCATE(eta2(1:nc_eta1+1,1:nc_eta2+1),error)
SLL_CLEAR_ALLOCATE(sol( 1:nc_eta1+1,1:nc_eta2+1),error)
SLL_CLEAR_ALLOCATE(phi( 1:nc_eta1+1,1:nc_eta2+1),error)
SLL_CLEAR_ALLOCATE(rhs( 1:nc_eta1+1,1:nc_eta2+1),error)

!set end points of solution rectangle in (x,y) space
eta1_min = 0.0
eta1_max = 4.0
eta2_min = 0.0
eta2_max = 4.0

delta_eta1 = (eta1_max-eta1_min)/float(nc_eta1)
delta_eta2 = (eta2_max-eta2_min)/float(nc_eta2)
do i=1,nc_eta1+1
   do j=1,nc_eta2+1
      eta1(i,j) = eta1_min+float(i-1)*delta_eta1
      eta2(i,j) = eta2_min+float(j-1)*delta_eta2
   end do
end do

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

eta1_min = -5.0
eta1_max =  5.0
eta2_min = -5.0
eta2_max =  5.0

delta_eta1 = (eta1_max-eta1_min)/float(nc_eta1)
delta_eta2 = (eta2_max-eta2_min)/float(nc_eta2)
do i=1,nc_eta1+1
   do j=1,nc_eta2+1
      eta1(i,j) = eta1_min+float(i-1)*delta_eta1
      eta2(i,j) = eta2_min+float(j-1)*delta_eta2
   end do
end do
     
call initialize_mudpack_cartesian(dirichlet,                    &
                                  eta1_min, eta1_max, nc_eta1,  &
                                  eta2_min, eta2_max, nc_eta2,  &
                                  SLL_DIRICHLET, SLL_DIRICHLET, &
                                  SLL_DIRICHLET, SLL_DIRICHLET)


sol = exp(-(eta1*eta1+eta2*eta2))

call sll_gnuplot_corect_2d(eta1_min, eta1_max, nc_eta1+1, &
                           eta2_min, eta2_max, nc_eta2+1, &
                           sol, "sol_dirichlet", 1, error)


do j=2,nc_eta2
   do i=2,nc_eta1
      rhs(i,j) = (sol(i-1,j)-2.*sol(i,j)+sol(i+1,j))/(delta_eta1*delta_eta1) &
               + (sol(i,j-1)-2.*sol(i,j)+sol(i,j+1))/(delta_eta2*delta_eta2)
   end do
end do

!rhs = 4 * sol * (eta1*eta1 + eta2*eta2 - 1)

call sll_gnuplot_corect_2d(eta1_min, eta1_max, nc_eta1+1, &
                           eta2_min, eta2_max, nc_eta2+1, &
                           rhs, "rhs_dirichlet", 1, error)

!rhs = 4.0_f64
phi(:,1) = sol(:,1)
phi(:,nc_eta2+1) = sol(:,nc_eta2+1)
phi(1,:) = sol(1,:)
phi(nc_eta2+1,:) = sol(nc_eta2+1,:)

call solve_mudpack_cartesian(dirichlet, phi, rhs)

call sll_gnuplot_corect_2d(eta1_min, eta1_max, nc_eta1+1, &
                           eta2_min, eta2_max, nc_eta2+1, &
                           phi, "phi_dirichlet", 1, error)

!compute and print maximum norm of error
write(*,201) maxval(abs(phi-sol))

print*,"PASSED"

201 format(' maximum error  =  ',e10.3)
     
end program test_mudpack_cartesian
