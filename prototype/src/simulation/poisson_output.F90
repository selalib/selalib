program poisson_output
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"

  use poisson_polar
  use sll_constants
  implicit none

  type(sll_plan_poisson_polar), pointer :: plan
  sll_real64 :: rmin,rmax,dr,dtheta,r1,r2,r,theta,x,y
  sll_int32 :: nr,ntheta,i,j
  sll_real64, dimension(:,:), allocatable :: f,phi0

  rmin=1.0_f64
  rmax=10.0_f64

  nr=512
  ntheta=64

  dr=(rmax-rmin)/real(nr,f64)
  dtheta=2.0_f64*sll_pi/real(ntheta,f64)

  SLL_ALLOCATE(f(nr+1,ntheta+1),i)
  SLL_ALLOCATE(phi0(nr+1,ntheta+1),i)

  r1=4.0_f64
  r2=5.0_f64
  do i=1,nr+1
     r=rmin+real(i-1,f64)*dr
     if (r>=r1 .and. r<=r2) then
        f(i,:)=1._f64
     end if
  end do

  plan => new_plan_poisson_polar(dr,rmin,nr,ntheta)
  call poisson_solve_polar(plan,f,phi0)

  open(21,file='phi0.dat')
  do i=1,nr+1
     r=rmin+real(i-1,f64)*dr
     do j=1,ntheta+1
        theta=real(j-1,f64)*dtheta
        x=r*cos(theta)
        y=r*sin(theta)
        write(21,*)r,theta,x,y,phi0(i,j)
     end do
  end do

end program poisson_output
