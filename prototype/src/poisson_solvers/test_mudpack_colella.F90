program test_mudpack_colella
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_constants.h"

use sll_mudpack_colella

implicit none

sll_real64, dimension(:,:), allocatable :: rhs
sll_real64, dimension(:,:), allocatable :: phi
sll_real64, dimension(:,:), allocatable :: phi_sol
sll_real64, dimension(:),   allocatable :: eta1
sll_real64, dimension(:),   allocatable :: eta2

sll_int32  :: i, j
sll_int32  :: nc_eta1
sll_int32  :: nc_eta2
sll_real64 :: eta1_min, eta1_max, delta_eta1
sll_real64 :: eta2_min, eta2_max, delta_eta2,alpha
sll_int32  :: error

sll_int32, parameter  :: n = 4

type(mudpack_2d) :: poisson

integer :: mode = 2

eta1_min   = 0.5_f64
eta1_max   = 1.0_f64 !*sll_pi

eta2_min  = 0.0_f64
eta2_max  = 1.0_f64  !* sll_pi
alpha = 0.1_f64

nc_eta1 = 64
nc_eta2 = 64
delta_eta1 = (eta1_max-eta1_min)/real(nc_eta1,f64)
delta_eta2 = 2.0_f64*sll_pi/real(nc_eta2,f64)

SLL_CLEAR_ALLOCATE(rhs(1:nc_eta1+1,1:nc_eta2+1),error)
SLL_CLEAR_ALLOCATE(phi(1:nc_eta1+1,1:nc_eta2+1),error)
SLL_CLEAR_ALLOCATE(phi_sol(1:nc_eta1+1,1:nc_eta2+1),error)
SLL_CLEAR_ALLOCATE(eta1(1:nc_eta1+1),error)
SLL_CLEAR_ALLOCATE(eta2(1:nc_eta2+1),error)

do i = 1, nc_eta1+1
   eta1(i)=eta1_min+(i-1)*delta_eta1
end do
do j = 1, nc_eta2+1
   eta2(j)=eta2_min+(j-1)*delta_eta2
end do

open(10,file="phi_colella_sol.dat")
do j=1,nc_eta2+1
   do i=1,nc_eta1+1
      phi_sol(i,j) = cos(eta1(i)+alpha*sin(2*sll_pi*eta1(i))*sin(2*sll_pi*eta2(j)))* &
                     cos(eta2(j)+alpha*sin(2*sll_pi*eta1(i))*sin(2*sll_pi*eta2(j)))
      write(10,*) sngl(eta1(i)+alpha*sin(2*sll_pi*eta1(i))*sin(2*sll_pi*eta2(j))), &
                  sngl(eta2(j)+alpha*sin(2*sll_pi*eta1(i))*sin(2*sll_pi*eta2(j))), &
                  sngl(phi_sol(i,j))
   end do
   write(10,*)
end do
close(10)

!tol   = 1.0e-14_f64

do i =1,nc_eta1+1
   do j=1,nc_eta2+1
      rhs(i,j) = 0.0_f64 !f_cos(eta1(i), eta2(j))
   end do
end do

call initialize_poisson_colella_mudpack(poisson, phi, rhs, &
                                  eta1_min, eta1_max, &
                                  eta2_min, eta2_max, &
                                  nc_eta1,nc_eta2)



print*,'PASS 1'
call solve_poisson_colella_mudpack(poisson,phi, rhs)
print*,'PASS 2'
!changer de polaire -> colella
open(11,file="phi_colella.dat")
do j = 1, nc_eta2+1
   do i = 1, nc_eta1+1
      write(11,*) sngl(eta1(i)+alpha*sin(2*sll_pi*eta1(i))*sin(2*sll_pi*eta2(j))), &
                  sngl(eta2(j)+alpha*sin(2*sll_pi*eta1(i))*sin(2*sll_pi*eta2(j))), &
                  sngl(phi(i,j)),      &
                  sngl(phi_sol(i,j))
   end do
   write(11,*)
end do
close(11)



contains

sll_real64 function f_cos( r, eta2 )

   !sage: assume(r>=1)
   !sage: assume(r<=2)
   !sage: phi = (r-eta1_min)*(r-eta1_max)*r*cos(n*eta2)
   !sage: diff(r*diff(phi,r),r)/r + diff(phi,eta2,eta2)/(r*r)

   sll_real64 :: r
   sll_real64 :: eta2

   f_cos = 0.0_8


end function f_cos

end program test_mudpack_colella
