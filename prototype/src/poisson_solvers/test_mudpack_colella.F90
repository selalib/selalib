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

#define alpha 0.2
#define mode 2

sll_int32, parameter  :: n = 4

eta1_min  = 0.0_f64
eta1_max  = 2*sll_pi
eta2_min  = 0.0_f64
eta2_max  = 2*sll_pi

nc_eta1 = 32
nc_eta2 = 32
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
      x = eta1(i) !+ alpha*sin(eta1(i))*sin(eta2(j))
      y = eta2(j) !+ alpha*sin(eta1(i))*sin(eta2(j))
      phi(i,j) = sin(mode*x)*sin(mode*y)
      rhs(i,j) = -2*mode*mode*phi(i,j)
      write(11,*) sngl(x), sngl(y), sngl(phi(i,j))
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
      err(i,j) = phi(i,j) - sin(mode*x)*sin(mode*y)
      write(12,*) sngl(x),sngl(y),sngl(phi(i,j))
   end do
   write(12,*)
end do


write(*,201) maxval(abs(err))
print*,'PASSED'
201 format(' maximum error  =  ',e10.3)

end program test_mudpack_colella

!> input pde coefficients at any grid point (x,y) in the solution region
!> (xa.le.x.le.xb,yc.le.y.le.yd) to mud2cr
subroutine coef(eta_1,eta_2,cxx,cxy,cyy,cx,cy,ce)
implicit none
real(8) :: eta_1,eta_2,cxx,cxy,cyy,cx,cy,ce,pi, beta
pi = 4*atan(1.)

beta = (512*pi**3*alpha**3*sin(pi*eta_1)**3*cos(pi*eta_1)**3*cos(pi*eta_2)**6       &
       - 64*pi**3*alpha**3*sin(pi*eta_1)**3*cos(pi*eta_1)**3 - &
       48*pi**2*alpha**2*sin(pi*eta_1)**2*cos(pi*eta_1)**2 + &
       64*(8*pi**3*alpha**3*sin(pi*eta_1)**6 - 12*pi**3*alpha**3*sin(pi*eta_1)**4 + &
       6*pi**3*alpha**3*sin(pi*eta_1)**2 - &
       pi**3*alpha**3)*sin(pi*eta_2)**3*cos(pi*eta_2)**3 - &
       192*(4*pi**3*alpha**3*sin(pi*eta_1)**3*cos(pi*eta_1)**3 + &
       pi**2*alpha**2*sin(pi*eta_1)**2*cos(pi*eta_1)**2)*cos(pi*eta_2)**4 - &
       12*pi*alpha*sin(pi*eta_1)*cos(pi*eta_1) + &
       48*(8*(4*pi**3*alpha**3*sin(pi*eta_1)**5*cos(pi*eta_1) - &
       4*pi**3*alpha**3*sin(pi*eta_1)**3*cos(pi*eta_1) + &
       pi**3*alpha**3*sin(pi*eta_1)*cos(pi*eta_1))*cos(pi*eta_2)**4 - &
       (16*pi**3*alpha**3*sin(pi*eta_1)**5*cos(pi*eta_1) - &
       16*pi**3*alpha**3*sin(pi*eta_1)**3*cos(pi*eta_1) + &
       4*pi**3*alpha**3*sin(pi*eta_1)*cos(pi*eta_1) + &
       4*pi**2*alpha**2*sin(pi*eta_1)**4 - 4*pi**2*alpha**2*sin(pi*eta_1)**2 + &
       pi**2*alpha**2)*cos(pi*eta_2)**2)*sin(pi*eta_2)**2 + &
       24*(16*pi**3*alpha**3*sin(pi*eta_1)**3*cos(pi*eta_1)**3 + &
       8*pi**2*alpha**2*sin(pi*eta_1)**2*cos(pi*eta_1)**2 + &
       pi*alpha*sin(pi*eta_1)*cos(pi*eta_1))*cos(pi*eta_2)**2 + &
       12*(64*(2*pi**3*alpha**3*sin(pi*eta_1)**4*cos(pi*eta_1)**2 - &
       pi**3*alpha**3*sin(pi*eta_1)**2*cos(pi*eta_1)**2)*cos(pi*eta_2)**5 - &
       16*(8*pi**3*alpha**3*sin(pi*eta_1)**4*cos(pi*eta_1)**2 - &
       4*pi**3*alpha**3*sin(pi*eta_1)**2*cos(pi*eta_1)**2 + &
       2*pi**2*alpha**2*sin(pi*eta_1)**3*cos(pi*eta_1) - &
       pi**2*alpha**2*sin(pi*eta_1)*cos(pi*eta_1))*cos(pi*eta_2)**3 + &
       (32*pi**3*alpha**3*sin(pi*eta_1)**4*cos(pi*eta_1)**2 + &
       16*pi**2*alpha**2*sin(pi*eta_1)**3*cos(pi*eta_1) - &
       8*pi**2*alpha**2*sin(pi*eta_1)*cos(pi*eta_1) - &
       2*(8*pi**3*alpha**3*cos(pi*eta_1)**2 - pi*alpha)*sin(pi*eta_1)**2 - &
       pi*alpha)*cos(pi*eta_2))*sin(pi*eta_2) - 1)


cxx = 0.
cxy = 0.
     
cyy = 0.
cx  = 0.0 
cy  = 0.0 
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
