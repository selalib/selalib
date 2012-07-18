program vlaspois
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"

  use sll_fft
  use polar_advection
  use poisson_polar
  use numeric_constants
  implicit none

  sll_int32 :: i, j, step, err
  sll_int32 :: nr, ntheta, nb_step
  sll_int32 :: fcase, scheme
  sll_real64 :: dr, dtheta, rmin, rmax, r, theta, dt, tf, x, y, r1, r2
  sll_real64 :: w0, w, l10, l1, l20, l2, e, e0
  sll_real64, dimension(:,:), pointer :: f, phi ,fdemi
  sll_real64, dimension(:,:,:), pointer :: grad_phi
  type(sll_fft_plan), pointer ::pfwd, pinv

  sll_int32 :: mod
  sll_real64 :: mode

  !for fcase=3
  !modes is used to test the fft with f(r)*cos(mode*theta)
  !namelist /modes/ mod
  mod=1
  !read(*,NML=modes)
  mode=real(mod,f64)

  rmin=1.0_f64
  rmax=10.0_f64

  ! number of step in r and theta directions
  ! /= of number of points
  nr=128
  ntheta=128

  dr=real(rmax-rmin,f64)/real(nr,f64)
  dtheta=2.0_f64*sll_pi/real(ntheta,f64)

  !choose the way to define dt, tf and nb_step
  !the tree ways are equivalent
  !we should have dt<=0.1*dr
  !default
  tf=1.0_f64
  dt=0.1_f64*dr
  nb_step=floor(tf/dt)

!!$  !definition of dt=tf/nb_step
!!$  tf=1.0_f64
!!$  nb_step=-1
!!$  dt=tf/real(nb_step,f64)

  !definition of nb_step=tf/dt
  dt=0.05_f64*dr
  tf=5.0_f64
  nb_step=floor(tf/dt)

!!$  !definition of tf=dt*nb_step
!!$  nb_step=1000
!!$  dt=0.1_f64
!!$  tf=dt*real(nb_step,f64)

  print*,'# nb_step =',nb_step,' dt =',dt,'tf =',real(nb_step,f64)*dt

  SLL_ALLOCATE(f(nr+1,ntheta+1),err)
  SLL_ALLOCATE(phi(nr+1,ntheta+1),err)
  SLL_ALLOCATE(grad_phi(2,nr+1,ntheta+1),err)

  !initialization of FFT
  pfwd => fft_new_plan(ntheta,f(1,1:ntheta),f(1,1:ntheta),FFT_FORWARD,FFT_NORMALIZE)
  pinv => fft_new_plan(ntheta,phi(1,1:ntheta),phi(1,1:ntheta),FFT_INVERSE)

  phi=0.0_f64

  !distribution function
  ! 1 : gaussienne in r, constant in theta
  ! 2 : f(r,theta)=1[r1,r2](r)*cos(theta)
  ! 3 : test distribution for poisson solver
  ! 4 : (gaussienne in r)*cos(theta)
  fcase=4

  !chose the way to calcul
  ! 1 : Semi-Lagrangien scheme
  ! 2 : Semi-Lagrangien scheme with control
  scheme=1
  if (scheme==2) then
     SLL_ALLOCATE(fdemi(nr+1,ntheta+1),err)
  end if

  if (fcase==1) then
     do i=1,nr+1
        r=rmin+real(i-1,f64)*dr
        f(i,:)=1.0_f64/(1.0_f64*sqrt(2.0_f64*sll_pi))*exp(-(r-(real(rmax-rmin)/2.0_f64))**2/(2.0_f64*1.0_f64**2))
     end do

  else if (fcase==2) then
     r1=4.0_f64
     r2=5.0_f64
     do i=1,nr+1
        r=rmin+real(i-1,f64)*dr
        if (r>=r1 .and. r<=r2) then
           do j=1,ntheta+1
              theta=real(j-1,f64)*dtheta
              f(i,j)=cos(3.0_f64*theta)
           end do
        end if
     end do

  else if (fcase==3) then
     do i=1,nr+1
        r=rmin+real(i-1,f64)*dr
        do j=1,ntheta+1
           theta=real(j-1,f64)*dtheta
           f(i,j)=-(r-rmin)*(r-rmax)/r**2*((36.0_f64-mode**2)*r**4+(2.0_f64*mode**2-39.0_f64)*r**3*(rmin+rmax) &
                & +(9.0_f64-mode**2)*r**2*(rmin**2+rmax**2)+(30.0_f64-4.0_f64*mode**2)*r**2*rmin*rmax &
                & +(2.0_f64*mode**2-3.0_f64)*r*rmin*rmax*(rmin+rmax)-mode**2*rmin**2*rmax**2) &
                & *sin(mode*theta)
        end do
     end do

  else if (fcase==4) then
     do i=1,nr+1
        r=rmin+real(i-1,f64)*dr
        do j=1,ntheta+1
           theta=real(j-1,f64)*dtheta
           f(i,j)=1.0_f64/(5.0_f64*sqrt(2.0_f64*sll_pi))*exp(-(r-(real(rmax-rmin)/2.0_f64))**2/50.0_f64)*cos(theta)
        end do
     end do

  else
     print*,"f is not defined"
     print*,'see line 76 of file selalib/prototype/src/simulation/VlasovPoisson.F90'
     print*,'can not go any further'
     print*,'exiting...'
     stop
  end if

  !write f in a file before calculations
  open (unit=20,file='CGinit.dat')
  do i=1,nr+1
     r=rmin+real(i-1,f64)*dr
     do j=1,ntheta+1
        theta=real(j-1,f64)*dtheta
        x=r*cos(theta)
        y=r*sin(theta)
        write(20,*)r,theta,x,y,f(i,j)
     end do
     write(20,*)' '
  end do
  close(20)

  open(unit=23,file='thdiag.dat')
  write(23,*)'#tf = ',tf,'  nb_step = ',nb_step,'  dt = ',dt
  write(23,*)'#   t   //   w   //   l1   //   l2   //   e' 
  w0=0.0_f64
  l10=0.0_f64
  l20=0.0_f64
  e0=0.0_f64
  call poisson_solve_polar(f,rmin,dr,nr,ntheta,pfwd,pinv,phi)
  phi(:,ntheta+1)=phi(:,1)
  call compute_grad_field(nr,ntheta,dr,dtheta,rmin,rmax,phi,grad_phi)
  do i=1,nr
     r=rmin+real(i-1,f64)*dr
     do j=1,ntheta
        w0=w0+r*f(i,j)
        l10=l10+r*abs(f(i,j))
        l20=l20+r*f(i,j)**2
        e0=e0+r*(grad_phi(1,i,j)**2+grad_phi(2,i,j)**2)
     end do
  end do

  w0=w0*dr*dtheta
  l10=l10*dr*dtheta
  l20=sqrt(l20*dr*dtheta)
  e0=e0*dr*dtheta
  write(23,*)'#',w0,l10,l20,e0
  write(23,*)0.0_f64,1.0_f64,1.0_f64,1.0_f64,0.0_f64
  !print*,'e=',e0, log(e0)

  do step=1,nb_step
     !initialisation of weight (w), l1, l2 and energy (e)
     w=0.0_f64
     l1=0.0_f64
     l2=0.0_f64
     e=0.0_f64

     if (scheme==1) then
        !classical semi-Lagrangian scheme
        call SL_classic(dt,dr,dtheta,nr,ntheta,rmin,rmax,pfwd,pinv,f,phi,grad_phi)

     else if (scheme==2) then
        !semi-Lagrangian scheme with control
        call SL_controlled(dt,dr,dtheta,nr,ntheta,rmin,rmax,pfwd,pinv,f,fdemi,phi,grad_phi)

     else
        print*,'no scheme define'
        print*,"the program won't do anything"
        print*,'see line 85 of file selalib/prototype/src/simulation to solve the probleme'
        print*,'exiting the loop'
        exit
     end if

     f(:,ntheta+1)=f(:,1)
     grad_phi(:,:,ntheta+1)=grad_phi(:,:,1)
     phi(:,ntheta+1)=phi(:,1)

     do i=1,nr
        r=rmin+real(i-1,f64)*dr
        do j=1,ntheta
           w=w+r*f(i,j)
           l1=l1+r*abs(f(i,j))
           l2=l2+r*f(i,j)**2
           e=e+r*(grad_phi(1,i,j)**2+grad_phi(2,i,j)**2)
        end do
     end do
     w=w*dr*dtheta
     l1=l1*dr*dtheta
     l2=sqrt(l2*dr*dtheta)
     e=e*dr*dtheta
     write(23,*)dt*real(step,f64),w/w0,l1/l10,l2/l20,e-e0

     if ((step/10)*10==step) then
        print*,'#step',step
     end if
  end do
  close(23)

  !write the final f in a file
  !w0=0.0_f64
  !w=0.0_f64
  open (unit=21,file='CGfinal.dat')
  do i=1,nr+1
     r=rmin+real(i-1,f64)*dr
     do j=1,ntheta+1
        theta=real(j-1,f64)*dtheta
        x=r*cos(theta)
        y=r*sin(theta)
        !w0=max(w0,abs(phi(i,j)))
        !w=max(w,abs(phi(i,j)-(r-rmin)**3*(r-rmax)**3*sin(mode*theta)))
        write(21,*)r,theta,x,y,f(i,j)!,phi(i,j),(r-rmin)**3*(r-rmax)**3*cos(mode*theta)
     end do
     write(21,*)' '
  end do
  close(21)
  !print*,dr,w0,w,w/w0,'#dr, w0, w,w/w0'

  call fft_delete_plan(pinv)
  call fft_delete_plan(pfwd)

  SLL_DEALLOCATE(f,err)
  SLL_DEALLOCATE(phi,err)
  SLL_DEALLOCATE(grad_phi,err)
  if (scheme==2) then
     SLL_DEALLOCATE(fdemi,err)
  end if

end program vlaspois
