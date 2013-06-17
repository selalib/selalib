program test_mudpack_colella
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_constants.h"

use sll_mudpack_colella

implicit none

sll_real64, dimension(:,:), allocatable :: rhs
sll_real64, dimension(:,:), allocatable :: phi
sll_real64, dimension(:,:), allocatable :: phi_cos
sll_real64, dimension(:,:), allocatable :: phi_sin
sll_real64, dimension(:),   allocatable :: eta1
sll_real64, dimension(:),   allocatable :: eta2

sll_int32  :: i, j
sll_int32  :: nc_eta1
sll_int32  :: nc_eta2
sll_real64 :: eta1_min, eta1_max, delta_eta1
sll_real64 :: eta2_min, eta2_max, delta_eta2
sll_int32  :: error
sll_real64 :: tol,l1,l2,linf

sll_int32, parameter  :: n = 4

eta1_min   = 0.0_f64
eta1_max   = 2.0_f64*sll_pi

eta2_min  = 0.0_f64
eta2_max  = 2.0_f64 * sll_pi

nc_eta1 = 64
nc_eta2 = 64
delta_eta1 = (eta1_max-eta1_min)/real(nc_eta1,f64)
delta_eta2 = 2.0_f64*sll_pi/real(nc_eta2,f64)

SLL_CLEAR_ALLOCATE(rhs(1:nc_eta1+1,1:nc_eta2+1),error)
SLL_CLEAR_ALLOCATE(phi(1:nc_eta1+1,1:nc_eta2+1),error)
SLL_CLEAR_ALLOCATE(phi_cos(1:nc_eta1+1,1:nc_eta2+1),error)
SLL_CLEAR_ALLOCATE(phi_sin(1:nc_eta1+1,1:nc_eta2+1),error)
SLL_CLEAR_ALLOCATE(eta1(1:nc_eta1+1),error)
SLL_CLEAR_ALLOCATE(eta2(1:nc_eta2+1),error)

do i = 1, nc_eta1+1
   eta1(i)=eta1_min+(i-1)*delta_eta1
end do
do j = 1, nc_eta2+1
   eta2(j)=eta2_min+(j-1)*delta_eta2
end do

open(10,file="phi_colella.dat")
do j=1,nc_eta2+1
   do i=1,nc_eta1+1
      phi_cos(i,j) = 0.0 !(eta1(i)-eta1_min)*(eta1(i)-eta1_max)*cos(n*theta(j))*eta1(i)
      phi_sin(i,j) = (eta1(i)-eta1_min)*(eta1(i)-eta1_max)*cos(n*eta2(j))*eta1(i)
      write(10,*) sngl(eta1(i)*cos(eta2(j))), &
                  sngl(eta1(i)*sin(eta2(j))), &
                  sngl(phi_cos(i,j)),      &
                  sngl(phi_sin(i,j))
   end do
   write(10,*)
end do
close(10)

tol   = 1.0e-14_f64

do i =1,nc_eta1+1
   do j=1,nc_eta2+1
      rhs(i,j) = 0 !f_cos(eta1(i), eta2(j))
   end do
end do

call initialize_poisson_colella_mudpack(phi_cos, rhs, &
                                      eta1_min, eta1_max, &
                                      eta2_min, eta2_max, &
                                      nc_eta1, nc_eta2)
print*,'PASS 1'
call solve_poisson_colella_mudpack(phi_cos, rhs)
print*,'PASS 2'
!changer de polaire -> colella
do j = 1, nc_eta2+1
   do i = 1, nc_eta1+1
      write(11,*) sngl(eta1(i)*cos(eta2(j))), &
                  sngl(eta1(i)*sin(eta2(j))), &
                  sngl(phi_cos(i,j)),      &
                  sngl(phi_sin(i,j))
   end do
   write(11,*)
end do

l1   = 0.0_f64
l2   = 0.0_f64
linf = 0.0_f64

if (l1>tol .or. l2>tol .or. linf>tol) then
   print*,'FAILED',tol,l1,l2,linf
   stop
end if

do i =1,nc_eta1+1
   do j=1,nc_eta2+1
      rhs(i,j)= f_sin(eta1(i),eta2(j))
   end do
end do

if (l1>tol .or. l2>tol .or. linf>tol) then
   print*,'FAILED',tol,l1,l2,linf
   stop
end if

print*,'PASSED'

contains

sll_real64 function f_cos( r, eta2 )

   !sage: assume(r>=1)
   !sage: assume(r<=2)
   !sage: phi = (r-eta1_min)*(r-eta1_max)*r*cos(n*eta2)
   !sage: diff(r*diff(phi,r),r)/r + diff(phi,eta2,eta2)/(r*r)

   sll_real64 :: r
   sll_real64 :: eta2

   f_cos = -(r-eta1_max)*(r-eta1_min)*n*n*cos(n*eta2)/r &
           + ((r-eta1_max)*(r-eta1_min)*cos(n*eta2)  &
           + (r-eta1_max)*r*cos(n*eta2) + (r-eta1_min)*r*cos(n*eta2) &
           + 2*((r-eta1_max)*cos(n*eta2) + (r-eta1_min)*cos(n*eta2) &
           + r*cos(n*eta2))*r)/r


end function f_cos

sll_real64 function f_sin( r, eta2)

   !sage: assume(r>=1)
   !sage: assume(r<=2)
   !sage: phi = (r-eta1_min)*(r-eta1_max)*r*sin(n*eta2)
   !sage: diff(r*diff(phi,r),r)/r + diff(phi,eta2,eta2)/(r*r)

   sll_real64 :: r
   sll_real64 :: eta2
   
   f_sin = -(r-eta1_max)*(r-eta1_min)*n*n*sin(n*eta2)/r &
         + ((r-eta1_max)*(r-eta1_min)*sin(n*eta2) &
         + (r-eta1_max)*r*sin(n*eta2) + (r-eta1_min)*r*sin(n*eta2) &
         + 2*((r-eta1_max)*sin(n*eta2) + (r-eta1_min)*sin(n*eta2)  &
         + r*sin(n*eta2))*r)/r

end function f_sin

end program test_mudpack_colella
