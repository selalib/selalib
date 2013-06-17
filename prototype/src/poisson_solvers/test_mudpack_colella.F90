program test_mudpack_colella
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_constants.h"

use sll_mudpack_colella

implicit none

type(mudpack_2d) :: poisson
sll_real64, dimension(:,:), allocatable :: rhs
sll_real64, dimension(:,:), allocatable :: phi
sll_real64, dimension(:,:), allocatable :: err
sll_real64, dimension(:),   allocatable :: eta1
sll_real64, dimension(:),   allocatable :: eta2

sll_int32  :: i, j
sll_int32  :: nc_eta1
sll_int32  :: nc_eta2
sll_real64 :: eta1_min, eta1_max, delta_eta1
sll_real64 :: eta2_min, eta2_max, delta_eta2
sll_int32  :: error
sll_real64 :: tol,l1,l2,linf, x, y

#define alpha 0.00
#define mode 2

<<<<<<< HEAD
eta1_min   = 0.0_f64
eta1_max   = 2.0_f64*sll_pi
=======
sll_int32, parameter  :: n = 4
>>>>>>> origin/prototype-devel

eta1_min  = 0.0_f64
eta1_max  = 1.0_f64
eta2_min  = 0.0_f64
eta2_max  = 1.0_f64

<<<<<<< HEAD
nc_eta1 = 64
nc_eta2 = 64
=======
nc_eta1 = 32
nc_eta2 = 32
>>>>>>> origin/prototype-devel
delta_eta1 = (eta1_max-eta1_min)/real(nc_eta1,f64)
delta_eta2 = (eta2_max-eta2_min)/real(nc_eta2,f64)

SLL_CLEAR_ALLOCATE(rhs(1:nc_eta1+1,1:nc_eta2+1),error)
SLL_CLEAR_ALLOCATE(phi(1:nc_eta1+1,1:nc_eta2+1),error)
SLL_CLEAR_ALLOCATE(eta1(1:nc_eta1+1),error)
SLL_CLEAR_ALLOCATE(eta2(1:nc_eta2+1),error)

do i = 1, nc_eta1+1
   eta1(i)=eta1_min+(i-1)*delta_eta1
end do
do j = 1, nc_eta2+1
   eta2(j)=eta2_min+(j-1)*delta_eta2
end do

do j=1,nc_eta2+1
   do i=1,nc_eta1+1
<<<<<<< HEAD
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
=======
      x = eta1(i) + alpha*sin(2*sll_pi*eta1(i))*sin(2*sll_pi*eta2(j))
      y = eta2(j) + alpha*sin(2*sll_pi*eta1(i))*sin(2*sll_pi*eta2(j))
      phi(i,j) = sin(2*sll_pi*mode*x)*sin(2*sll_pi*mode*y)
      rhs(i,j) = -8*sll_pi*sll_pi*mode*mode*phi(i,j)
      write(11,*) sngl(x), sngl(y), sngl(phi(i,j))
>>>>>>> origin/prototype-devel
   end do
   write(11,*)
end do

call initialize_poisson_colella_mudpack(poisson, phi, rhs, &
                                        eta1_min, eta1_max, &
                                        eta2_min, eta2_max, &
                                        nc_eta1, nc_eta2)

call solve_poisson_colella_mudpack(poisson, phi, rhs)

SLL_CLEAR_ALLOCATE(err(1:nc_eta1+1,1:nc_eta2+1),error)
do j=1,nc_eta2+1
   do i=1,nc_eta1+1
      x = eta1(i) + alpha*sin(eta1(i))*sin(eta2(j))
      y = eta2(j) + alpha*sin(eta1(i))*sin(eta2(j))
      err(i,j) = phi(i,j) - exact(x,y)
      write(12,*) sngl(x),sngl(y),sngl(phi(i,j))
   end do
   write(12,*)
end do


write(*,201) maxval(abs(err))
print*,'PASSED'
201 format(' maximum error  =  ',e10.3)

contains

real(8) function exact(x, y)
   real(8) :: x, y
   exact = sin(2*sll_pi*mode*x)*sin(2*sll_pi*mode*y)
   return

end function exact

end program test_mudpack_colella
!> input pde coefficients at any grid point (x,y) in the solution region
!> (xa.le.x.le.xb,yc.le.y.le.yd) to mud2cr
subroutine coef(x,y,cxx,cxy,cyy,cx,cy,ce)
implicit none
real(8) :: c1, c2, s1, s2, beta
real(8) :: c12, c13, c14, c15, c16
real(8) :: c22, c23, c24, c25, c26
real(8) :: s12, s13, s14, s15, s16
real(8) :: s22, s23, s24, s25, s26
real(8) :: pi2, pi3, pi4
real(8) :: alpha2, alpha3
real(8) :: x,y,cxx,cxy,cyy,cx,cy,ce,pi
pi = 4*atan(1.)
pi2 = pi*pi
pi3 = pi2*pi
pi4 = pi3*pi

c1 = cos(pi*x); c12=c1*c1; c13=c12*c1; c14=c13*c1; c15=c14*c1; c16=c15*c1
c2 = cos(pi*y); c22=c2*c2; c23=c22*c2; c14=c23*c2; c25=c24*c2; c26=c25*c2
s1 = sin(pi*x); s12=s1*s1; s13=s12*s1; s14=s13*s1; s15=s14*s1; s16=s15*s1
s2 = sin(pi*y); s22=s2*s2; s23=s22*s2; s14=s23*s2; s25=s24*s2; s26=s25*s2

alpha2 = alpha*alpha; alpha3 = alpha2*alpha

beta=512*alpha3*c16*c23*pi3*s23 + 1536*alpha3*c15*c24*pi3*s1*s22 &
+ 1536*alpha3*c14*c25*pi3*s12*s2 + 512*alpha3*c13*c26*pi3*s13 - &
768*alpha3*c15*c22*pi3*s1*s22 - 1536*alpha3*c14*c23*pi3*s12*s2 - &
768*alpha3*c14*c23*pi3*s23 - 768*alpha3*c13*c24*pi3*s13 - &
1536*alpha3*c13*c24*pi3*s1*s22 - 768*alpha3*c12*c25*pi3*s12*s2 + &
192*pi2*alpha2*c14*c22*s22 + 384*pi2*alpha2*c13*c23*s1*s2 + &
192*pi2*alpha2*c12*c24*s12 + 384*alpha3*c14*c2*pi3*s12*s2 + &
384*alpha3*c13*c22*pi3*s13 + 768*alpha3*c13*c22*pi3*s1*s22 + &
768*alpha3*c12*c23*pi3*s12*s2 + 384*alpha3*c12*c23*pi3*s23 + &
384*alpha3*c1*c24*pi3*s1*s22 - 192*pi2*alpha2*c13*c2*s1*s2 - &
192*pi2*alpha2*c12*c22*s12 - 192*pi2*alpha2*c12*c22*s22 - &
192*pi2*alpha2*c1*c23*s1*s2 - 64*alpha3*c13*pi3*s13 - &
192*alpha3*c12*c2*pi3*s12*s2 - 192*alpha3*c1*c22*pi3*s1*s22 - &
64*alpha3*c23*pi3*s23 + 48*pi2*alpha2*c12*s12 + &
96*pi2*alpha2*c1*c2*s1*s2 + 48*pi2*alpha2*c22*s22 + &
24*pi*alpha*c12*c2*s2 + 24*pi*alpha*c1*c22*s1 - 12*pi*alpha*c1*s1 - &
12*pi*alpha*c2*s2 + 1

cxx=-1024*alpha3*c16*c25*pi3*s2 - 1024*alpha3*c15*c26*pi3*s1 + &
1024*alpha3*c16*c23*pi3*s2 + 1536*alpha3*c15*c24*pi3*s1 + &
1536*alpha3*c14*c25*pi3*s2 + 1024*alpha3*c13*c26*pi3*s1 - &
256*pi2*alpha2*c14*c24 + 128*pi2*alpha2*c13*c23*s1*s2 - &
256*alpha3*c16*c2*pi3*s2 - 768*alpha3*c15*c22*pi3*s1 - &
1536*alpha3*c14*c23*pi3*s2 - 1536*alpha3*c13*c24*pi3*s1 - &
512*alpha3*c12*c25*pi3*s2 + 256*pi2*alpha2*c14*c22 - &
64*pi2*alpha2*c13*c2*s1*s2 + 256*pi2*alpha2*c12*c24 - &
64*pi2*alpha2*c1*c23*s1*s2 + 128*alpha3*c15*pi3*s1 + &
384*alpha3*c14*c2*pi3*s2 + 768*alpha3*c13*c22*pi3*s1 + &
512*alpha3*c12*c23*pi3*s2 - 64*pi2*alpha2*c14 - &
256*pi2*alpha2*c12*c22 + 32*pi2*alpha2*c1*c2*s1*s2 - &
128*alpha3*c13*pi3*s1 - 128*alpha3*c12*c2*pi3*s2 + &
64*pi2*alpha2*c12 + 8*pi*alpha*c12*c2*s2 + 24*pi*alpha*c1*c22*s1 - &
12*pi*alpha*c1*s1 - 4*pi*alpha*c2*s2 + 1

cyy=-1024*alpha3*c16*c25*pi3*s2 - 1024*alpha3*c15*c26*pi3*s1 + &
1024*alpha3*c16*c23*pi3*s2 + 1536*alpha3*c15*c24*pi3*s1 + &
1536*alpha3*c14*c25*pi3*s2 + 1024*alpha3*c13*c26*pi3*s1 - &
256*pi2*alpha2*c14*c24 + 128*pi2*alpha2*c13*c23*s1*s2 - &
512*alpha3*c15*c22*pi3*s1 - 1536*alpha3*c14*c23*pi3*s2 - &
1536*alpha3*c13*c24*pi3*s1 - 768*alpha3*c12*c25*pi3*s2 - &
256*alpha3*c1*c26*pi3*s1 + 256*pi2*alpha2*c14*c22 - &
64*pi2*alpha2*c13*c2*s1*s2 + 256*pi2*alpha2*c12*c24 - &
64*pi2*alpha2*c1*c23*s1*s2 + 512*alpha3*c13*c22*pi3*s1 + &
768*alpha3*c12*c23*pi3*s2 + 384*alpha3*c1*c24*pi3*s1 + &
128*alpha3*c25*pi3*s2 - 256*pi2*alpha2*c12*c22 + &
32*pi2*alpha2*c1*c2*s1*s2 - 64*pi2*alpha2*c24 - &
128*alpha3*c1*c22*pi3*s1 - 128*alpha3*c23*pi3*s2 + &
64*pi2*alpha2*c22 + 24*pi*alpha*c12*c2*s2 + 8*pi*alpha*c1*c22*s1 - &
4*pi*alpha*c1*s1 - 12*pi*alpha*c2*s2 + 1 

cxy=2048*alpha3*c16*c25*pi3*s2 + 2048*alpha3*c15*c26*pi3*s1 - &
2048*alpha3*c16*c23*pi3*s2 - 3072*alpha3*c15*c24*pi3*s1 - &
3072*alpha3*c14*c25*pi3*s2 - 2048*alpha3*c13*c26*pi3*s1 + &
256*pi2*alpha2*c14*c24 - 512*pi2*alpha2*c13*c23*s1*s2 + &
512*alpha3*c16*c2*pi3*s2 + 1024*alpha3*c15*c22*pi3*s1 + &
3072*alpha3*c14*c23*pi3*s2 + 3072*alpha3*c13*c24*pi3*s1 + &
1024*alpha3*c12*c25*pi3*s2 + 512*alpha3*c1*c26*pi3*s1 - &
256*pi2*alpha2*c14*c22 + 256*pi2*alpha2*c13*c2*s1*s2 - &
256*pi2*alpha2*c12*c24 + 256*pi2*alpha2*c1*c23*s1*s2 - &
768*alpha3*c14*c2*pi3*s2 - 1024*alpha3*c13*c22*pi3*s1 - &
1024*alpha3*c12*c23*pi3*s2 - 768*alpha3*c1*c24*pi3*s1 + &
32*pi2*alpha2*c14 + 256*pi2*alpha2*c12*c22 - &
128*pi2*alpha2*c1*c2*s1*s2 + 32*pi2*alpha2*c24 + &
256*alpha3*c12*c2*pi3*s2 + 256*alpha3*c1*c22*pi3*s1 - &
32*pi2*alpha2*c12 - 32*pi2*alpha2*c22 - 16*pi*alpha*c12*c2*s2 - &
16*pi*alpha*c1*c22*s1 + 8*pi*alpha*c1*s1 + 8*pi*alpha*c2*s2 

cx=512*pi4*alpha3*c15*c2*s1*s2 + 512*pi4*alpha3*c1*c25*s1*s2 - &
512*pi4*alpha3*c13*c2*s1*s2 - 512*pi4*alpha3*c1*c23*s1*s2 + &
256*pi4*alpha3*c1*c2*s1*s2 + 32*pi2*alpha*c1*c2*s1*s2 + &
64*alpha2*c13*pi3*s1 + 64*alpha2*c23*pi3*s2 - 32*alpha2*c1*pi3*s1 - &
32*alpha2*c2*pi3*s2 

cy=512*pi4*alpha3*c15*c2*s1*s2 + 512*pi4*alpha3*c1*c25*s1*s2 - &
512*pi4*alpha3*c13*c2*s1*s2 - 512*pi4*alpha3*c1*c23*s1*s2 + &
256*pi4*alpha3*c1*c2*s1*s2 + 32*pi2*alpha*c1*c2*s1*s2 + &
64*alpha2*c13*pi3*s1 + 64*alpha2*c23*pi3*s2 - 32*alpha2*c1*pi3*s1 - &
32*alpha2*c2*pi3*s2

cxx = cxx / beta
cyy = cyy / beta
cx = cx / beta
cy = cy / beta

ce  = 0.0 
return
end subroutine

!> at upper y boundary
subroutine bnd(kbdy,xory,alfa,beta,gama,gbdy)
implicit none
integer  :: kbdy
real(8)  :: xory,alfa,beta,gama,gbdy

!! Set bounday condition value

return
end subroutine
