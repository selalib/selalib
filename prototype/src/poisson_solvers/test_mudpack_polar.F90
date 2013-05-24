program test_poisson_polar
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_constants.h"

use sll_mudpack_polar

implicit none

sll_real64, dimension(:,:), allocatable :: rhs
sll_real64, dimension(:,:), allocatable :: phi
sll_real64, dimension(:,:), allocatable :: phi_cos
sll_real64, dimension(:,:), allocatable :: phi_sin
sll_real64, dimension(:),   allocatable :: r
sll_real64, dimension(:),   allocatable :: theta

sll_int32  :: i, j
sll_int32  :: nr
sll_int32  :: ntheta
sll_real64 :: r_min, r_max, delta_r
sll_real64 :: theta_min, theta_max, delta_theta
sll_int32  :: error
sll_real64 :: tol,l1,l2,linf
sll_int32, parameter  :: n = 4

print*,'Testing the Poisson solver in 2D, polar coordinate'

r_min   = 1.0_f64
r_max   = 2.0_f64

theta_min  = 0.0_f64
theta_max  = 2.0_f64 * sll_pi

nr     = 33
ntheta = 129
delta_r     = (r_max-r_min)/real(nr-1,f64)
delta_theta = 2.0_f64*sll_pi/real(ntheta-1,f64)

SLL_CLEAR_ALLOCATE(rhs(1:nr,1:ntheta),error)
SLL_CLEAR_ALLOCATE(phi(1:nr,1:ntheta),error)
SLL_CLEAR_ALLOCATE(phi_cos(1:nr,1:ntheta),error)
SLL_CLEAR_ALLOCATE(phi_sin(1:nr,1:ntheta),error)
SLL_CLEAR_ALLOCATE(r(1:nr),error)
SLL_CLEAR_ALLOCATE(theta(1:ntheta),error)

print*, r_min
do i = 1, nr
   r(i)=r_min+(i-1)*delta_r
end do
print*, r
do j = 1, ntheta
   theta(j)=(j-1)*delta_theta
end do

open(10,file="phi_polar.dat")
do j=1,ntheta
   do i=1,nr
      phi_cos(i,j) = 0.0 !(r(i)-r_min)*(r(i)-r_max)*cos(n*theta(j))*r(i)
      phi_sin(i,j) = (r(i)-r_min)*(r(i)-r_max)*cos(n*theta(j))*r(i)
      write(10,*) sngl(r(i)*cos(theta(j))), &
                  sngl(r(i)*sin(theta(j))), &
                  sngl(phi_cos(i,j)),      &
                  sngl(phi_sin(i,j))
   end do
   write(10,*)
end do
close(10)

tol   = 1.0e-14_f64

do i =1,nr
   do j=1,ntheta
      rhs(i,j) = f_cos(r(i), theta(j))
   end do
end do

call solve_poisson_polar_mudpack(phi_cos, rhs, &
                                 r_min, r_max, &
                                 theta_min, theta_max, &
                                 nr, ntheta)

do j = 1, ntheta
   do i = 1, nr
      write(11,*) sngl(r(i)*cos(theta(j))), &
                  sngl(r(i)*sin(theta(j))), &
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

do i =1,nr
   do j=1,ntheta
      rhs(i,j)= f_sin(r(i),theta(j))
   end do
end do

if (l1>tol .or. l2>tol .or. linf>tol) then
   print*,'FAILED',tol,l1,l2,linf
   stop
end if

print*,'PASSED'

contains

sll_real64 function f_cos( r, theta )

   !sage: assume(r>=1)
   !sage: assume(r<=2)
   !sage: phi = (r-r_min)*(r-r_max)*r*cos(n*theta)
   !sage: diff(r*diff(phi,r),r)/r + diff(phi,theta,theta)/(r*r)

   sll_real64 :: r
   sll_real64 :: theta

   f_cos = -(r-r_max)*(r-r_min)*n*n*cos(n*theta)/r &
           + ((r-r_max)*(r-r_min)*cos(n*theta)  &
           + (r-r_max)*r*cos(n*theta) + (r-r_min)*r*cos(n*theta) &
           + 2*((r-r_max)*cos(n*theta) + (r-r_min)*cos(n*theta) &
           + r*cos(n*theta))*r)/r


end function f_cos

sll_real64 function f_sin( r, theta)

   !sage: assume(r>=1)
   !sage: assume(r<=2)
   !sage: phi = (r-r_min)*(r-r_max)*r*sin(n*theta)
   !sage: diff(r*diff(phi,r),r)/r + diff(phi,theta,theta)/(r*r)

   sll_real64 :: r
   sll_real64 :: theta
   
   f_sin = -(r-r_max)*(r-r_min)*n*n*sin(n*theta)/r &
         + ((r-r_max)*(r-r_min)*sin(n*theta) &
         + (r-r_max)*r*sin(n*theta) + (r-r_min)*r*sin(n*theta) &
         + 2*((r-r_max)*sin(n*theta) + (r-r_min)*sin(n*theta)  &
         + r*sin(n*theta))*r)/r

end function f_sin

end program test_poisson_polar
