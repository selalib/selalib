program test_mudpack_polar
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_constants.h"

use sll_mudpack_polar

implicit none

type(mudpack_2d) :: poisson
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
sll_real64 :: tolmax,errmax

sll_int32, parameter  :: n = 4

print*,'Testing the Poisson solver in 2D, polar coordinate'

r_min   = 1.0_f64
r_max   = 2.0_f64

theta_min = 0.0_f64
theta_max = 2.0_f64 * sll_pi

nr     = 33
ntheta = 129
delta_r     = (r_max-r_min)/(nr-1)
delta_theta = 2.0_f64*sll_pi/(ntheta-1)

SLL_CLEAR_ALLOCATE(rhs(1:nr,1:ntheta),error)
SLL_CLEAR_ALLOCATE(phi(1:nr,1:ntheta),error)
SLL_CLEAR_ALLOCATE(phi_cos(1:nr,1:ntheta),error)
SLL_CLEAR_ALLOCATE(phi_sin(1:nr,1:ntheta),error)
SLL_CLEAR_ALLOCATE(r(1:nr),error)
SLL_CLEAR_ALLOCATE(theta(1:ntheta),error)

do i = 1, nr
   r(i)=r_min+(i-1)*delta_r
end do
do j = 1, ntheta
   theta(j)=(j-1)*delta_theta
end do

do j=1,ntheta
   do i=1,nr
      phi_cos(i,j) = (r(i)-r_min)*(r(i)-r_max)*cos(n*theta(j))*r(i)
      phi_sin(i,j) = (r(i)-r_min)*(r(i)-r_max)*sin(n*theta(j))*r(i)
   end do
end do

tolmax   = 1.0e-4_f64

call initialize_poisson_polar_mudpack(poisson,  &
                                      r_min, r_max, nr, &
                                      theta_min, theta_max, ntheta, &
                                      DIRICHLET, DIRICHLET, PERIODIC, PERIODIC )
do i =1,nr
   do j=1,ntheta
      rhs(i,j) = f_cos(r(i), theta(j))
   end do
end do

poisson%iguess = 0 ! no initial guess

call solve_poisson_polar_mudpack(poisson, phi, rhs)

call plot_field( "mudpack_polar_cos.dat", phi )

errmax = maxval(abs(phi_cos-phi))
write(*,201) errmax
if ( errmax > tolmax ) then
   print*,'FAILED'
   stop
else
   print*,'PASSED'
end if

do i =1,nr
   do j=1,ntheta
      rhs(i,j) = f_sin(r(i), theta(j))
   end do
end do

poisson%iguess = 0 ! no initial guess

call solve_poisson_polar_mudpack(poisson, phi, rhs)

call plot_field( "mudpack_polar_sin.dat", phi )

errmax = maxval(abs(phi_sin-phi))
write(*,201) errmax
if ( errmax > tolmax ) then
   print*,'FAILED'
   stop
else
   print*,'PASSED'
end if

201 format(' maximum error  =  ',e10.3)

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

subroutine plot_field( filename, field )

   character(len=*) :: filename
   sll_real64       :: field(:,:)
   open(10, file=filename)
   do i =1,nr
      do j=1,ntheta
         write(10,*) sngl(r(i)*cos(theta(j))), &
                     sngl(r(i)*sin(theta(j))), &
                     sngl(field(i,j))
      end do
      write(10,*) 
   end do
   close(10)

end subroutine plot_field

end program test_mudpack_polar
