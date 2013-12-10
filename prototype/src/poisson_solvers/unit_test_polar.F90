program test_poisson
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"

  use polar_operators
  use sll_poisson_2d_polar
  use sll_constants
  implicit none

  !this code should become part of the unit_test for the Poisson solver in polar coordinates
  !the following is for compilation
  !copy it in the CMakeListe.txt :
  !  add_executable(test_poisson
  !  test_poisson.F90
  !  sll_poisson_2d_polar.F90
  !  )
  !  target_link_libraries(test_poisson
  !  sll_fft
  !  sll_tridiagonal
  !  sll_constants
  !  sll_memory
  !  sll_assertion
  !  sll_utilities
  !  #sll_poisson_2d
  !  )
  
  sll_real64, dimension(:,:), allocatable :: f, phi
  type(sll_plan_poisson_polar), pointer :: plan
  sll_int32 :: nr,ntheta,i,j,err,bc(2)
  sll_real64 :: dr,dtheta,rmin,rmax
  logical :: test
  sll_real64 :: tol,r,theta,a,l1,l2,linf
  sll_int32 :: mod

  print*,'#Testing the Poisson solver in 2D, polar coordinate'

  rmin=1.0_f64
  rmax=2.0_f64
  nr=512
  ntheta=256
  dr=(rmin-rmax)/real(nr,f64)
  dtheta=2.0_f64*sll_pi/real(ntheta,f64)

  SLL_ALLOCATE(f(nr+1,ntheta+1),err)
  SLL_ALLOCATE(phi(nr+1,ntheta+1),err)

  plan => new_plan_poisson_polar(dr,rmin,nr,ntheta,bc)

  tol=1.0e-14_f64
  test= .true.
  mod=0
  bc(1)=SLL_DIRICHLET
  bc(2)=SLL_DIRICHLET

  do while (test .and. mod<ntheta/2)
     do i =1,nr+1
        r=rmin+real(i-1,f64)*dr
        do j=1,ntheta
           theta=real(j-1,f64)*dtheta
           if (bc(1)+bc(2)==5) then
              a=rmax
              f(i,j)=-(2.0_f64+(2.0_f64*r-rmin-a)/r-(real(mod,f64)/r)**2*(r-rmin)*(r-a))*cos(real(mod,f64)*theta)
           else if (bc(1)+bc(2)==6) then
              a=2*rmax-rmin
              f(i,j)=-(2.0_f64+(2.0_f64*r-rmin-a)/r-(real(mod,f64)/r)**2*(r-rmin)*(r-a))*cos(real(mod,f64)*theta)
           end if
        end do
     end do
     call poisson_solve_polar(plan,f,phi)
     
     
     l1=0.0_f64
     l2=0.0_f64
     linf=0.0_f64

     do j=1,ntheta
        theta=real(j-1,f64)*dtheta
        do i=2,nr
           r=rmin+real(i-1,f64)*dr
           linf=max(linf,abs(phi(i,j)-(r-rmin)*(r-a)*cos(real(mod,f64)*theta)*r))
           l1=l1+abs(phi(i,j)-(r-rmin)*(r-a)*cos(real(mod,f64)*theta))*r
           l2=l2+(phi(i,j)-(r-rmin)*(r-a)*cos(real(mod,f64)*theta))**2*r
        end do
        linf=max(linf,abs(phi(1,j)))
        r=rmax
        linf=max(linf,abs(phi(i,j)-(r-rmin)*(r-a)*cos(real(mod,f64)*theta)*r))
        l1=l1+abs(phi(nr+1,j)-(r-rmin)*(r-a)*cos(real(mod,f64)*theta))/2.0_f64*r
        l2=l2+(phi(nr+1,j)-(r-rmin)*(r-a)*cos(real(mod,f64)*theta))**2/2.0_f64*r
     end do
     linf=linf*dr*dtheta
     l1=l1*dr*dtheta
     l2=sqrt(l2*dr*dtheta)

     if (l1>tol .or. l2>tol .or. linf>tol) then
        test=.false.
     end if

     mod=mod+1
     if (mod==ntheta/2 .and. (bc(1)+bc(2)==5)) then
        bc(1)=SLL_DIRICHLET
        bc(2)=SLL_NEUMANN
        mod=0
     end if
  end do

  do while (test .and. mod<ntheta/2)
     do i =1,nr+1
        r=rmin+real(i-1,f64)*dr
        do j=1,ntheta
           theta=real(j-1,f64)*dtheta
           if (bc(1)+bc(2)==5) then
              a=rmax
              f(i,j)=-(2.0_f64+(2.0_f64*r-rmin-a)/r-(real(mod,f64)/r)**2*(r-rmin)*(r-a))*sin(real(mod,f64)*theta)
           else if (bc(1)+bc(2)==6) then
              a=2*rmax-rmin
              f(i,j)=-(2.0_f64+(2.0_f64*r-rmin-a)/r-(real(mod,f64)/r)**2*(r-rmin)*(r-a))*sin(real(mod,f64)*theta)
           end if
        end do
     end do
     call poisson_solve_polar(plan,f,phi)

     l1=0.0_f64
     l2=0.0_f64
     linf=0.0_f64

     do j=1,ntheta
        theta=real(j-1,f64)*dtheta
        do i=2,nr
           r=rmin+real(i-1,f64)*dr
           linf=max(linf,abs(phi(i,j)-(r-rmin)*(r-a)*sin(real(mod,f64)*theta)*r))
           l1=l1+abs(phi(i,j)-(r-rmin)*(r-a)*sin(real(mod,f64)*theta))*r
           l2=l2+(phi(i,j)-(r-rmin)*(r-a)*sin(real(mod,f64)*theta))**2*r
        end do
        linf=max(linf,abs(phi(1,j)))
        r=rmax
        linf=max(linf,abs(phi(i,j)-(r-rmin)*(r-a)*sin(real(mod,f64)*theta)*r))
        l1=l1+abs(phi(nr+1,j)-(r-rmin)*(r-a)*sin(real(mod,f64)*theta))/2.0_f64*r
        l2=l2+(phi(nr+1,j)-(r-rmin)*(r-a)*sin(real(mod,f64)*theta))**2/2.0_f64*r
     end do
     linf=linf*dr*dtheta
     l1=l1*dr*dtheta
     l2=sqrt(l2*dr*dtheta)

     if (l1>tol .or. l2>tol .or. linf>tol) then
        test=.false.
     end if

     mod=mod+1
     if (mod==ntheta/2 .and. bc(1)+bc(2)==5) then
        bc(2)=SLL_NEUMANN
        bc(1)=SLL_DIRICHLET
        mod=0
     end if
  end do

  if (test) then
     print*,'PASS'
  else
     print*,'FAIL',tol,l1,l2,linf
  end if

end program test_poisson

