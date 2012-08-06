program cg_polar
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"

  use sll_timer
  use polar_kind
  use sll_fft
  use polar_advection
  use poisson_polar
  use numeric_constants
  implicit none

  type(polar_data), pointer :: data
  type(polar_VP_data), pointer :: adv
  type(polar_VP_rk4), pointer :: rk
  type(time_mark), pointer :: t1,t2,t3
  sll_real64, dimension (:,:), allocatable :: div
  sll_int32 :: i, j, step,fin
  sll_int32 :: nr, ntheta, nb_step
  sll_int32 :: fcase, scheme
  sll_real64 :: dr, dtheta, rmin, rmax, r, theta, dt, tf, x, y, r1, r2
  sll_real64 :: w0, w, l10, l1, l20, l2, e, e0
  character (len=40) :: cgf, thd
  sll_int32 :: mod
  sll_real64 :: mode,temps
  integer :: hh,min,ss
  integer, dimension(3) :: time

  !python script for fcase=3
  !modes is used to test the fft with f(r)*cos(mode*theta)
  !namelist /modes/ mod
  mod=3
  !read(*,NML=modes)
  mode=real(mod,f64)

  t1 => new_time_mark()
  t2 => new_time_mark()
  t3 => new_time_mark()

  rmin=1.0_f64
  rmax=10.0_f64

  ! number of step in r and theta directions
  ! /= of number of points
  nr=512
  ntheta=128

  dr=real(rmax-rmin,f64)/real(nr,f64)
  dtheta=2.0_f64*sll_pi/real(ntheta,f64)
  print*,'#dr=',dr,'dtheta=',dtheta

  !choose the way to define dt, tf and nb_step
  !the tree ways are equivalent
  !we should have dt<=0.1*dr
  !default
  tf=1.0_f64
  dt=0.1_f64*dr
  nb_step=ceiling(tf/dt)

!!$  !definition of dt=tf/nb_step
!!$  tf=5.0_f64
!!$  nb_step=5690
!!$  dt=tf/real(nb_step,f64)

!!$  !definition of nb_step=tf/dt
!!$  dt=0.05_f64*dr
!!$  tf=5.0_f64
!!$  nb_step=ceiling(tf/dt)
!!$
  !definition of tf=dt*nb_step
  nb_step=1
  dt=0.05_f64*dr
  tf=dt*real(nb_step,f64)

  tf=real(nb_step,f64)*dt
  fin=floor(tf+0.5_f64)
  print*,'# nb_step =',nb_step,' dt =',dt,'tf =',tf

  data => new_polar_data(nb_step,dt,rmin,rmax,nr,ntheta)
  adv => new_vp_data(data)
  rk => new_polar_vp_rk4(nr,ntheta)
  SLL_ALLOCATE(div(nr+1,ntheta+1),i)

!!$  adv%rr=1.0_f64
!!$  do i=1,nr+1
!!$     adv%rr(i)=rmax
!!$  end do

  adv%phi=0.0_f64
  adv%grad_phi=0.0_f64
  adv%f=0.0_f64

  !distribution function
  ! 1 : gaussienne in r, constant in theta
  ! 2 : f(r,theta)=1[r1,r2](r)*cos(theta)
  ! 3 : test distribution for poisson solver
  ! 4 : (gaussienne in r)*cos(theta)
  fcase=3

  !chose the way to calcul
  ! 1 : Semi-Lagrangien scheme
  ! 2 : Semi-Lagrangien scheme order 2
  ! 3 : ?jump-sheep? scheme
  scheme=2

  call filename('CGfinal',len('CGfinal'),fcase,scheme,mod,fin,cgf)
  call filename('thdiag',len('thdiag'),fcase,scheme,mod,fin,thd)

  if (fcase==1) then
     do i=1,nr+1
        r=rmin+real(i-1,f64)*dr
        adv%f(i,:)=1.0_f64/(1.0_f64*sqrt(2.0_f64*sll_pi))*exp(-(r-(real(rmax-rmin)/2.0_f64))**2/(2.0_f64*1.0_f64**2))
     end do

  else if (fcase==2) then
     r1=4.0_f64
     r2=5.0_f64
     do i=1,nr+1
        r=rmin+real(i-1,f64)*dr
        if (r>=r1 .and. r<=r2) then
           do j=1,ntheta+1
              theta=real(j-1,f64)*dtheta
               adv%f(i,j)=cos(3.0_f64*theta)
           end do
        end if
     end do

  else if (fcase==3) then
     do i=1,nr+1
        r=rmin+real(i-1,f64)*dr
        do j=1,ntheta+1
           theta=real(j-1,f64)*dtheta
            adv%f(i,j)=(r-rmin)*(r-rmax)/r**2*((36.0_f64-mode**2)*r**4+(2.0_f64*mode**2-39.0_f64)*r**3*(rmin+rmax) &
                & +(9.0_f64-mode**2)*r**2*(rmin**2+rmax**2)+(30.0_f64-4.0_f64*mode**2)*r**2*rmin*rmax &
                & +(2.0_f64*mode**2-3.0_f64)*r*rmin*rmax*(rmin+rmax)-mode**2*rmin**2*rmax**2) &
                & *cos(mode*theta)
        end do
     end do

  else if (fcase==4) then
     do i=1,nr+1
        r=rmin+real(i-1,f64)*dr
        do j=1,ntheta+1
           theta=real(j-1,f64)*dtheta
            adv%f(i,j)=1.0_f64/(0.5_f64*sqrt(2.0_f64*sll_pi))*exp(-(r-(real(rmax-rmin)/2.0_f64))**2/(2*0.5_f64**2))*sin(mode*theta)
        end do
     end do

  else
     print*,"f is not defined"
     print*,'see variable "fcase" in file selalib/prototype/src/simulation/CG_polar.F90'
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
        write(20,*)r,theta,x,y, adv%f(i,j)
     end do
     write(20,*)' '
  end do
  close(20)

  open(unit=23,file=thd)
  write(23,*)'#tf = ',tf,'  nb_step = ',nb_step,'  dt = ',dt
  write(23,*)'#   t   //   w   //   l1 rel  //   l2  rel //   e' 
  call poisson_solve_polar(adv)
  call compute_grad_field(adv)
  w0=0.0_f64
  l10=0.0_f64
  l20=0.0_f64
  e0=0.0_f64
  do j=1,ntheta
     w0=w0+(adv%f(1,j)*rmin+adv%f(nr+1,j)*rmax)/2.0_f64
     l10=l10+abs(adv%f(1,j)*rmin+adv%f(nr+1,j)*rmax)/2.0_f64
     l20=l20+(adv%f(1,j)/2.0_f64)**2*rmin+(adv%f(nr+1,j)/2.0_f64)**2*rmax
     e0=e0+rmin*(adv%grad_phi(1,1,j)/2.0_f64)**2+rmax*(adv%grad_phi(1,nr+1,j)/2.0_f64)**2+ &
          & rmin*(adv%grad_phi(2,1,j)/2.0_f64)**2+rmax*(adv%grad_phi(2,nr+1,j)/2.0_f64)**2
     do i=2,nr
        r=rmin+real(i-1,f64)*dr
        w0=w0+r*adv%f(i,j)
        l10=l10+r*abs(adv%f(i,j))
        l20=l20+r*adv%f(i,j)**2
        e0=e0+r*(adv%grad_phi(1,i,j)**2+adv%grad_phi(2,i,j)**2)
     end do
  end do
  w0=w0*dr*dtheta
  l10=l10*dr*dtheta
  l20=sqrt(l20*dr*dtheta)
  e0=e0*dr*dtheta
  write(23,*)'#t=0',w0,l10,l20,e0
  write(23,*)0.0_f64,w0,1.0_f64,1.0_f64,0.0_f64

  t1 => start_time_mark(t1)
  do step=1,nb_step

     if (step==101) then
        t2 => start_time_mark(t2)
        temps=time_elapsed_between(t1,t2)
        temps=temps/100*real(nb_step,f32)
        hh=floor(temps/3600.0d0)
        min=floor((temps-3600.0d0*real(hh))/60.0d0)
        ss=floor(temps-3600.0d0*real(hh)-60.0d0*real(min))
        print*,'# temps de calcul estimmé : ',hh,'h',min,'min',ss,'s'
        call itime(time)
        time(3)=time(3)+ss
        time(2)=time(2)+floor(real(time(3))/60.0)
        time(3)=time(3)-60*floor(real(time(3))/60.0)
        time(2)=time(2)+min
        time(1)=time(1)+floor(real(time(2))/60.0)
        time(2)=time(2)-60*floor(real(time(2))/60.0)
        time(1)=time(1)+hh
        print*,'#fin estimmée à',time(1),'h',time(2),"'",time(3),'"'
     end if

     if (scheme==1) then
        !classical semi-Lagrangian scheme
        call SL_classic(adv,rk)

     else if (scheme==2) then
        !semi-Lagrangian predictiv-correctiv scheme
        call SL_ordre_2(adv,rk)

     else if (scheme==3) then
        !?jump-sheep scheme?
        if (step==1) then
           call SL_ordre_2(adv,rk)
        else 
           call poisson_solve_polar(adv)
           call compute_grad_field(adv)
           adv%f_fft=adv%f
           adv%f=adv%fdemi
           adv%fdemi=adv%f_fft
           call advect_CG_polar(adv,rk)
        end if

     else
        print*,'no scheme define'
        print*,"the program won't do anything"
        print*,'see variable scheme in file selalib/prototype/src/simulation to solve the probleme'
        print*,'exiting the loop'
        exit
     end if

     adv%f(:,ntheta+1)=adv%f(:,1)
     adv%grad_phi(:,:,ntheta+1)=adv%grad_phi(:,:,1)
     adv%phi(:,ntheta+1)=adv%phi(:,1)

     !computation of mass (w), l1, l2 and energy (e)
     w=0.0_f64
     l1=0.0_f64
     l2=0.0_f64
     e=0.0_f64
     do j=1,ntheta
        w=w+(adv%f(1,j)*rmin+adv%f(nr+1,j)*rmax)/2.0_f64
        l1=l1+abs(adv%f(1,j)*rmin+adv%f(nr+1,j)*rmax)/2.0_f64
        l2=l2+(adv%f(1,j)/2.0_f64)**2*rmin+(adv%f(nr+1,j)/2.0_f64)**2*rmax
        e=e+rmin*(adv%grad_phi(1,1,j)/2.0_f64)**2+rmax*(adv%grad_phi(1,nr+1,j)/2.0_f64)**2+ &
             & rmin*(adv%grad_phi(2,1,j)/2.0_f64)**2+rmax*(adv%grad_phi(2,nr+1,j)/2.0_f64)**2
        do i=2,nr
           r=rmin+real(i-1,f64)*dr
           w=w+r*adv%f(i,j)
           l1=l1+r*abs(adv%f(i,j))
           l2=l2+r*adv%f(i,j)**2
           e=e+r*(adv%grad_phi(1,i,j)**2+adv%grad_phi(2,i,j)**2)
        end do
     end do
     w=w*dr*dtheta
     l1=l1*dr*dtheta
     l2=sqrt(l2*dr*dtheta)
     e=e*dr*dtheta
     write(23,*)dt*real(step,f64),w,l1/l10,l2/l20,e-e0

     if ((step/100)*100==step) then
        print*,'#step',step
     end if
  end do
  close(23)

  t3 => start_time_mark(t3)
  temps=time_elapsed_between(t1,t3)
  hh=floor(temps/3600.0d0)
  min=floor((temps-3600.0d0*real(hh))/60.0d0)
  ss=floor(temps-3600.0d0*real(hh)-60.0d0*real(min))
  print*,'# temps pour faire la boucle en temps : ',hh,'h',min,'min',ss,'s'

  !checking divergence of field
  call poisson_solve_polar(adv)
  call compute_grad_field(adv)
  call divergence_ortho_field(adv,div)

  !write the final f in a file
  open (unit=21,file=cgf)
  do i=1,nr+1
     r=adv%rr(i)
     do j=1,ntheta+1
     !j=1
        theta=adv%ttheta(j)
        x=r*cos(theta)
        y=r*sin(theta)
        !<for fase=3, checking the poisson solveur>
        !w0=max(w0,abs(phi(i,j)))
        !w=max(w,abs(phi(i,j)-(r-rmin)**3*(r-rmax)**3*sin(mode*theta)))
!!$        write(21,*)r,theta,x,y,adv%phi(i,j),(r-rmin)**3*(r-rmax)**3*cos(mode*theta)
        !</for fase=3, checking the poisson solveur>
        !write(21,*)r,theta,x,y,adv%f(i,j),div(i,j)
        write(21,*)r,theta,x,y,adv%grad_phi(1,i,j),adv%grad_phi(2,i,j),3.0_f64*(r-rmin)**2*(r-rmax)**2*(2.0_f64*r-rmin-rmax)*cos(mode*theta), &
             & -mode*(r-rmin)**3*(r-rmax)**3*sin(mode*theta)/r
     end do
     write(21,*)' '
  end do
!!$  write(21,*)' '
!!$  do j=1,ntheta+1
!!$     i=55
!!$     r=adv%rr(i)
!!$     theta=adv%ttheta(j)
!!$     x=r*cos(theta)
!!$     y=r*sin(theta)
!!$     write(21,*)r,theta,x,y,adv%grad_phi(1,i,j),adv%grad_phi(2,i,j)
!!$  end do
  close(21)

  SLL_DEALLOCATE_ARRAY(div,i)
  t1 => delete_time_mark(t1)
  t2 => delete_time_mark(t2)
  call vp_data_delete(adv)
  call vp_rk4_delete(rk)

contains

  subroutine filename(str,lenght,fcase,scheme,mod,tf,out)

    implicit none

    character (len=40), intent(out) :: out
    sll_int32, intent(in) :: scheme, mod,tf,fcase,lenght
    character (len=lenght) :: str

    integer :: i1,i2,i3
    character (len=2) :: mode,f,sch
    character (len=3) :: fin

    i1=mod/10
    i2=mod-i1
    mode=char(i1+48)//char(i2+48)
    i1=fcase/10
    i2=fcase-i1
    f=char(i1+48)//char(i2+48)
    i1=tf/100
    i2=(tf-100*i1)/10
    i3=tf-100*i1-10*i2
    fin=char(i1+48)//char(i2+48)//char(i3+48)
    i1=scheme/10
    i2=scheme-i1
    sch=char(i1+48)//char(i2+48)
    out=str//char(095)//'f'//f//char(095)//'s'//sch//char(095)//''//mode//char(095)//'i01'//char(095)//fin//'s.dat'

  end subroutine filename

end program cg_polar
