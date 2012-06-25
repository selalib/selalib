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
  sll_real64, dimension(:,:), pointer :: f, phi ,fdemi
  sll_real64, dimension(:,:,:), pointer :: grad_phi
  type(sll_fft_plan), pointer ::pfwd, pinv

  rmin=0.2_f64
  rmax=0.8_f64

  ! number of step in r and theta directions
  ! /= of number of points
  nr=64
  ntheta=64

  dr=real(rmax-rmin,f64)/real(nr,f64)
  dtheta=2.0_f64*sll_pi/real(ntheta,f64)

  tf=1.0_f64
  nb_step=20
  dt=tf/real(nb_step,f64)

  SLL_ALLOCATE(f(nr+1,ntheta+1),err)
  SLL_ALLOCATE(phi(nr+1,ntheta+1),err)
  SLL_ALLOCATE(grad_phi(2,nr+1,ntheta+1),err)

  !initialization of FFT
  pfwd => fft_new_plan(ntheta,f(1,1:ntheta),f(1,1:ntheta),FFT_FORWARD,FFT_NORMALIZE)
  pinv => fft_new_plan(ntheta,phi(1,1:ntheta),phi(1,1:ntheta),FFT_INVERSE,FFT_NORMALIZE)

  phi=0.0_f64

  !distribution function
  ! 1 : gaussienne in r, constant in theta
  ! 2 : f(r,theta)=1[r1,r2](r)*cos(theta)
  fcase=2

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
        f(i,:)=1.0_f64/(5.0_f64*sqrt(2.0_f64*sll_pi))*exp(-(r-(real(rmax-rmin)/2.0_f64))**2/50.0_f64)
     end do
  else if (fcase==2) then
     r1=rmin+(rmax-rmin)/3.0_f64
     r2=rmin+2.0_f64*(rmax-rmin)/3.0_f64
     do i=1,nr+1
        r=rmin+real(i-1,f64)*dr
        if (r>=r1 .and. r<=r2) then
           do j=1,ntheta+1
              theta=real(j-1,f64)*dtheta
              f(i,j)=cos(theta)
           end do
        end if
     end do
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
  end do
  close(20)

  do step=1,nb_step

     if (scheme==1) then
        !classical semi-Lagrangian scheme
        call SL_classic(dt,dr,dtheta,nr,ntheta,rmin,rmax,pfwd,pinv,f,phi,grad_phi)

     else if (scheme==2) then
        !semi-Lagrangian scheme with control
        call SL_controlled(dt,dr,dtheta,nr,ntheta,rmin,rmax,pfwd,pinv,f,fdemi,phi,grad_phi)
     end if

     f(:,ntheta+1)=f(:,1)
     grad_phi(:,:,ntheta+1)=grad_phi(:,:,1)
     phi(:,ntheta+1)=phi(:,1)
  
  end do

  !write the final f in a file
  open (unit=21,file='CGfinal.dat')
  do i=1,nr+1
     r=rmin+real(i-1,f64)*dr
     do j=1,ntheta+1
        theta=real(j-1,f64)*dtheta
        x=r*cos(theta)
        y=r*sin(theta)
        write(21,*)r,theta,x,y,f(i,j)
     end do
  end do
  close(21)

  call fft_delete_plan(pinv)
  call fft_delete_plan(pfwd)

end program vlaspois
