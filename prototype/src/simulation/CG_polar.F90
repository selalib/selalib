program cg_polar
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"

  use sll_timer
  use polar_operators
  use polar_advection
  use poisson_polar
  use numeric_constants
  implicit none

  type(sll_SL_polar), pointer :: plan_sl
  type(time_mark), pointer :: t1,t2,t3
  sll_real64, dimension (:,:), allocatable :: div,f,fp1,g
  sll_int32 :: i, j, step,fin,visustep
  sll_int32 :: nr, ntheta, nb_step
  sll_int32 :: fcase, scheme,carac,grad,visu
  sll_real64 :: dr, dtheta, rmin, rmax, r, theta, dt, tf, x, y, r1, r2
  sll_real64 :: w0, w, l10, l1, l20, l2, e, e0
  sll_int32 :: mod
  sll_real64 :: mode,temps,alpha
  sll_real64, dimension(2,2) :: dom
  integer :: hh,min,ss
  integer, dimension(3) :: time
  character (len=20) :: f_file

  !python script for fcase=3
  !modes is used to test the fft with f(r)*cos(mode*theta)
  !namelist /modes/ mod
  mod=3
  alpha = 1.e-3_f64
  !alpha = 1.e-3_f64
  !read(*,NML=modes)
  mode=real(mod,f64)

  t1 => new_time_mark()
  t2 => new_time_mark()
  t3 => new_time_mark()

  !>files 'CG_data.dat' and 'CG_case.dat' are included in directory selalib/prototype/src/simulation
  !>copy them in the same directory as the executable
  open(27,file='CG_data.dat',action="read")
  read(27,*)rmin
  read(27,*)rmax
  read(27,*)nr
  read(27,*)ntheta
  read(27,*)fin
  read(27,*)dt

  read(27,*)
  
  read(27,*)carac
  read(27,*)grad
  read(27,*)fcase
  read(27,*)scheme
  read(27,*)visu
  close(27)

  tf=real(fin,f64)

!!$  rmin=1.0_f64
!!$  rmax=10.0_f64
!!$
!!$  ! number of step in r and theta directions
!!$  ! /= of number of points
!!$  nr=256
!!$  ntheta=128

  dr=real(rmax-rmin,f64)/real(nr,f64)
  dtheta=2.0_f64*sll_pi/real(ntheta,f64)
  print*,'#dr=',dr,'dtheta=',dtheta
  dom(1,1)=rmin
  dom(1,2)=0.0_f64
  dom(2,1)=rmax
  dom(2,2)=2.0_f64*sll_pi

  !choose the way to define dt, tf and nb_step
  !the tree ways are equivalent
  !we should have dt<=0.1*dr

!!$  !definition of dt=tf/nb_step
!!$  tf=5.0_f64
!!$  nb_step=5690
!!$  dt=tf/real(nb_step,f64)

  !definition of nb_step=tf/dt
  dt=dt*dr
  !tf=100.0_f64
  nb_step=ceiling(tf/dt)

!!$  !definition of tf=dt*nb_step
!!$  nb_step=0
!!$  dt=0.05_f64*dr
!!$  tf=dt*real(nb_step,f64)

  tf=real(nb_step,f64)*dt
  print*,'# nb_step =',nb_step,' dt =',dt,'tf =',tf
  visustep=2000

!!$  !scheme to compute caracteristics
!!$  ! 1 : using explicit Euler method
!!$  ! 2 : rotation, rotation speed = -1
!!$  ! 3 : using symplectic Euler with linear interpolation
!!$  ! 4 : using symplectic Verlet with linear interpolation
!!$  ! 5 : using fixed point method
!!$  ! 6 : using modified symplectic Euler
!!$  ! 7 : using modified symplectic Verlet
!!$  ! 8 : using modified fixed point
!!$  carac=4
!!$
!!$  !scheme to compute gradian
!!$  ! 1 : final differencies in r and theta
!!$  ! 2 : fft in theta, final differencies in r
!!$  ! 3 : splines in r and theta
!!$  grad=3
!!$
!!$  !distribution function
!!$  ! 1 : gaussian in r, constant in theta
!!$  ! 2 : f(r,theta)=1_[r1,r2](r)*(1+cos(theta))
!!$  ! 3 : test distribution for poisson solver
!!$  ! 4 : (gaussian in r)*(1+cos(theta))
!!$  ! 5 : read f in a file with syntax : r theta x y f(i,j)
!!$  fcase=2
!!$  !f_file='CGfinal04.dat' !not working
!!$
!!$  !choose the way to compute
!!$  ! 1 : Semi-Lagrangian scheme order 1
!!$  ! 2 : Semi-Lagrangian scheme order 2
!!$  ! 3 : leap-frog scheme
!!$  scheme=2
!!$
!!$  !choose the visualization
!!$  ! 0 : gnuplot
!!$  ! 1 : vtk
!!$  visu=0

  plan_sl => new_SL(rmin,rmax,dr,dtheta,dt,nr,ntheta,grad,carac)
  SLL_ALLOCATE(div(nr+1,ntheta+1),i)
  SLL_ALLOCATE(f(nr+1,ntheta+1),i)
  SLL_ALLOCATE(g(ntheta+1,nr+1),i)
  SLL_ALLOCATE(fp1(nr+1,ntheta+1),i)

  step=0
  f=0.0_f64

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
              f(i,j)=1._f64+alpha*cos(mode*theta)
           end do
           g(j,i)=f(i,j)
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
                & *cos(mode*theta)
        end do
     end do

  else if (fcase==4) then
     do i=1,nr+1
        r=rmin+real(i-1,f64)*dr
        do j=1,ntheta+1
           theta=real(j-1,f64)*dtheta
           f(i,j)=1.0_f64/(0.5_f64*sqrt(2.0_f64*sll_pi))*exp(-(r-(real(rmax-rmin)/2.0_f64))**2/(2*0.5_f64**2))*(1.0_f64+alpha*sin(mode*theta))
        end do
     end do

  else if (fcase==5) then
     open(25,file='CGfinal.dat',action="read")
     read(25,*)
     read(25,*)
     read(25,*)
     do i=i,nr+1
        do j=1,ntheta+1
           print*,i,j
           read(25,'(2X,5(1F18.16,8X))')r,theta,x,y,f(i,j)
        end do
        !read(25,*)
     end do
     close(25)

  else
     print*,"f is not defined"
     print*,'see variable "fcase" in file selalib/prototype/src/simulation/CG_polar.F90'
     print*,'can not go any further'
     print*,'exiting...'
     stop
     print*,'so far so good'
  end if

  fp1=0.0_f64

  call poisson_solve_polar(plan_sl%poisson,f,plan_sl%phi)
  call compute_grad_field(plan_sl%grad,plan_sl%phi,plan_sl%adv%field)
  !do i=1,nr+1
  !   r=rmin+dr*real(i-1,f64)
  !   do j=1,ntheta+1
  !      temps=plan_sl%adv%field(1,i,j)/r
  !      plan_sl%adv%field(1,i,j)=-plan_sl%adv%field(2,i,j)
  !      plan_sl%adv%field(2,i,j)=temps
  !   end do
  !end do
  call divergence_scalar_field(plan_sl%grad,plan_sl%adv%field,div)

  !write f in a file before calculations
  call print2dper(dom,f(1:nr+1,1:ntheta),Nr+1,Ntheta,visu,step,"CG")
!!$  open (unit=20,file='CGinit.dat')
!!$    do i=1,nr+1
!!$     r=rmin+real(i-1,f64)*dr
!!$     do j=1,ntheta+1
!!$        theta=real(j-1,f64)*dtheta
!!$        x=r*cos(theta)
!!$        y=r*sin(theta)
!!$        write(20,*)r,theta,x,y,f(i,j),div(i,j)
!!$     end do
!!$     write(20,*)' '
!!$  end do
!!$  close(20)





!!$  call poisson_solve_polar(data,adv%f,adv%phi,adv%f_fft,adv%fk,adv%phik,adv%a,adv%cts,adv%ipiv,adv%pfwd,adv%pinv)
!!$  call compute_grad_field(data,adv%phi,adv%grad_phi,adv%spl_phi)
!!$  open (unit=21,file='test.dat')
!!$  do i=1,nr+1
!!$     r=adv%rr(i)
!!$     do j=1,ntheta+1
!!$        theta=adv%ttheta(j)
!!$        x=r*cos(theta)
!!$        y=r*sin(theta)
!!$        write(21,*)r,theta,x,y,adv%grad_phi(1,i,j),adv%grad_phi(2,i,j),adv%phi(i,j), &
!!$             & 3.0_f64*(r-rmin)**2*(r-rmax)**2*(2.0_f64*r-rmin-rmax)*cos(mode*theta), &
!!$             & -mode*(r-rmin)**3*(r-rmax)**3*sin(mode*theta)/r, (r-rmin)**3*(r-rmax)**3*cos(mode*theta)
!!$     end do
!!$     write(21,*)' '
!!$  end do
!!$
!!$  stop





  open(unit=23,file='thdiag.dat')
  write(23,*)'#fcase',fcase,'scheme',scheme,'mode',mode,'grad',grad,'carac',carac
  write(23,*)'#nr',nr,'ntheta',ntheta
  write(23,*)'#tf = ',tf,'  nb_step = ',nb_step,'  dt = ',dt
  write(23,*)'#   t   //   w   //   l1 rel  //   l2  rel //   e' 

  do i=1,nr+1
     r=rmin+real(i-1,f64)*dr
     plan_sl%adv%field(2,i,:)=plan_sl%adv%field(2,i,:)/r
  end do

  w0=0.0_f64
  l10=0.0_f64
  l20=0.0_f64
  e0=0.0_f64
  do j=1,ntheta
     w0=w0+(f(1,j)*rmin+f(nr+1,j)*rmax)/2.0_f64
     l10=l10+abs(f(1,j)*rmin+f(nr+1,j)*rmax)/2.0_f64
     l20=l20+(f(1,j)/2.0_f64)**2*rmin+(f(nr+1,j)/2.0_f64)**2*rmax
     e0=e0+rmin*(plan_sl%adv%field(1,1,j)/2.0_f64)**2+rmax*(plan_sl%adv%field(1,nr+1,j)/2.0_f64)**2+ &
          & rmin*(plan_sl%adv%field(2,1,j)/2.0_f64)**2+rmax*(plan_sl%adv%field(2,nr+1,j)/2.0_f64)**2
     do i=2,nr
        r=rmin+real(i-1,f64)*dr
        w0=w0+r*f(i,j)
        l10=l10+r*abs(f(i,j))
        l20=l20+r*f(i,j)**2
        e0=e0+r*(plan_sl%adv%field(1,i,j)**2+plan_sl%adv%field(2,i,j)**2)
     end do
  end do
  w0=w0*dr*dtheta
  l10=l10*dr*dtheta
  l20=sqrt(l20*dr*dtheta)
  e0=e0*dr*dtheta/2.0_f64
  write(23,*)'#t=0',w0,l10,l20,e0
  write(23,*)0.0_f64,w0,1.0_f64,1.0_f64,0.0_f64,e0

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
        !classical semi-Lagrangian scheme (order 1)
        call SL_classic(plan_sl,f,fp1)

     else if (scheme==2) then
        !semi-Lagrangian predictive-corrective scheme
        call SL_ordre_2(plan_sl,f,fp1)

!!$     else if (scheme==3) then
!!$        !leap-frog scheme
!!$        if (step==1) then
!!$           call SL_ordre_2()
!!$        else 
!!$           call poisson_solve_polar()
!!$           call compute_grad_field()
!!$           call advect_CG_polar()
!!$        end if

     else
        print*,'no scheme define'
        print*,"the program won't do anything"
        print*,'see variable scheme in file selalib/prototype/src/simulation to solve the probleme'
        print*,'exiting the loop'
        exit
     end if

     fp1(:,ntheta+1)=fp1(:,1)
     f=fp1

     !computation of mass (w), l1, l2 and energy (e)
     do i=1,nr+1
        r=rmin+real(i-1,f64)*dr
        plan_sl%adv%field(2,i,:)=plan_sl%adv%field(2,i,:)/r
     end do
     w=0.0_f64
     l1=0.0_f64
     l2=0.0_f64
     e=0.0_f64
     do j=1,ntheta
        w=w+(f(1,j)*rmin+f(nr+1,j)*rmax)/2.0_f64
        l1=l1+abs(f(1,j)*rmin+f(nr+1,j)*rmax)/2.0_f64
        l2=l2+(f(1,j)/2.0_f64)**2*rmin+(f(nr+1,j)/2.0_f64)**2*rmax
        e=e+rmin*(plan_sl%adv%field(1,1,j)/2.0_f64)**2+rmax*(plan_sl%adv%field(1,nr+1,j)/2.0_f64)**2+ &
             & rmin*(plan_sl%adv%field(2,1,j)/2.0_f64)**2+rmax*(plan_sl%adv%field(2,nr+1,j)/2.0_f64)**2
        do i=2,nr
           r=rmin+real(i-1,f64)*dr
           w=w+r*f(i,j)
           l1=l1+r*abs(f(i,j))
           l2=l2+r*f(i,j)**2
           e=e+r*(plan_sl%adv%field(1,i,j)**2+plan_sl%adv%field(2,i,j)**2)
        end do
     end do
     w=w*dr*dtheta
     l1=l1*dr*dtheta
     l2=sqrt(l2*dr*dtheta)
     e=e*dr*dtheta/2.0_f64
     write(23,*)dt*real(step,f64),w,l1/l10,l2/l20,e-e0,e

     if ((step/500)*500==step) then
        print*,'#step',step
     end if
     
     if (step/visustep*visustep==step) then
        call print2dper(dom,f(1:nr+1,1:ntheta),Nr+1,Ntheta,visu,step,"CG")
     end if

!!$     if (abs(real(step)*dt-125.)<=1e-3) then
!!$        open(24,file='125s.dat')
!!$        do i=1,nr+1
!!$           r=rmin+real(i-1,f64)*dr
!!$           do j=1,ntheta+1
!!$              theta=real(j-1,f64)*dtheta
!!$              x=r*cos(theta)
!!$              y=r*sin(theta)
!!$              write(24,*)r,theta,x,y,f(i,j)
!!$           end do
!!$        end do
!!$        close(24)
!!$     end if

  end do
  close(23)

  t3 => start_time_mark(t3)
  temps=time_elapsed_between(t1,t3)
  hh=floor(temps/3600.0d0)
  min=floor((temps-3600.0d0*real(hh))/60.0d0)
  ss=floor(temps-3600.0d0*real(hh)-60.0d0*real(min))
  print*,'# temps pour faire la boucle en temps : ',hh,'h',min,'min',ss,'s'

  !checking divergence of field
  call poisson_solve_polar(plan_sl%poisson,f,plan_sl%phi)
  call compute_grad_field(plan_sl%grad,plan_sl%phi,plan_sl%adv%field)
  do i=1,nr+1
     r=rmin+dr*real(i-1,f64)
     do j=1,ntheta+1
        temps=plan_sl%adv%field(1,i,j)/r
        plan_sl%adv%field(1,i,j)=-plan_sl%adv%field(2,i,j)
        plan_sl%adv%field(2,i,j)=temps
     end do
  end do
  call divergence_scalar_field(plan_sl%grad,plan_sl%adv%field,div)

  !write the final f in a file
  call print2dper(dom,f(1:nr+1,1:ntheta),Nr+1,Ntheta,visu,step,"CG")
  open (unit=21,file='CGfinal.dat')
  write(21,*)'#fcase',fcase,'scheme',scheme,'mode',mode,'grad',grad,'carac',carac
  write(21,*)'#nr',nr,'ntheta',ntheta
  write(21,*)'#tf = ',tf,'  nb_step = ',nb_step,'  dt = ',dt
  do i=1,nr+1
     r=rmin+real(i-1,f64)*dr
     do j=1,ntheta+1
        theta=real(j-1,f64)*dtheta
        x=r*cos(theta)
        y=r*sin(theta)
        write(21,*)r,theta,x,y,f(i,j)!,div(i,j),plan_sl%phi(i,j),&
        !& plan_sl%adv%field(1,i,j),plan_sl%adv%field(2,i,j)
     end do
     !write(21,*)' '
  end do
  close(21)

  SLL_DEALLOCATE_ARRAY(div,i)
  SLL_DEALLOCATE_ARRAY(f,i)
  t1 => delete_time_mark(t1)
  t2 => delete_time_mark(t2)
  t3 => delete_time_mark(t3)
  call delete_SL_polar(plan_sl)

end program cg_polar
