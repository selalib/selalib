program cg_polar
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"

  use sll_timer 
  use polar_operators
  use polar_advection
  use sll_fft
  !use poisson_polar
  use sll_constants
  implicit none

  type(sll_SL_polar), pointer :: plan_sl
  type(sll_time_mark)  :: t1,t2,t3
  sll_real64, dimension (:,:), allocatable :: div,f,fp1,g,phi_ref
  sll_real64, dimension (:)  , allocatable :: int_r
  sll_int32  :: i, j, step, visustep, hh, min, ss,ii
  sll_int32  :: nr, ntheta, nb_step
  sll_int32  :: fcase, scheme, carac, grad, visu, interp_case, PPM_order
  sll_int32  :: ierr_poiss
  sll_real64 :: dr, dtheta, rmin, rmax, r, theta, dt, tf, r1, r2
  sll_real64 :: w0, w, l10, l1, l20, l2, e, e0
  sll_int32  :: mod, bc(2)
  sll_real64 :: mode, temps, temps_mode, alpha, tmp, err_loc
  sll_real64, dimension(2,2) :: dom
  character (len=16) :: f_file!,bctop,bcbot
  !used for testing poisson with fcase=2
  sll_real64 :: c1, c2, c3, k1, k2, k3, x, y
  sll_real64 :: c1_mode, c2_mode, c3_mode, k1_mode, k2_mode, k3_mode
  sll_real64 :: mode_slope(1:8),time_mode(1:8)

  !>files 'CG_data.dat'is included in directory selalib/prototype/src/simulation
  !>copy it in the same directory as the executable
  open(27,file='CG_data.txt',action="read")
  read(27,*)rmin
  read(27,*)rmax
  read(27,*)nr
  read(27,*)ntheta
  read(27,*)r1
  read(27,*)r2
  read(27,*)alpha
  read(27,*)mod
  read(27,*)nb_step
  read(27,*)dt
  read(27,*)visustep
  read(27,*)
  read(27,*)carac
  read(27,*)grad
  read(27,*)fcase
  read(27,*)scheme
  read(27,*)interp_case
  read(27,*)PPM_order
  read(27,*)visu
  read(27,*)f_file
  read(27,*)
  read(27,*)bc(1)
  read(27,*)bc(2)
  close(27)

  mode   = real(mod,f64)
  dr     = real(rmax-rmin,f64)/real(nr,f64)
  dtheta = 2.0_f64*sll_pi/real(ntheta,f64)

  call print_defaultfftlib()
  dom(1,1) = rmin
  dom(1,2) = 0.0_f64
  dom(2,1) = rmax
  dom(2,2) = 2.0_f64*sll_pi

  !--> computation of the final time
  tf = dt*real(nb_step,f64)
  print*,'# nb_step =',nb_step,' dt =',dt,'tf =',tf

  plan_sl => new_SL(rmin,rmax,dr,dtheta,dt,nr,ntheta,grad,carac,bc)
  SLL_ALLOCATE(div(nr+1,ntheta+1),i)
  SLL_ALLOCATE(f(nr+1,ntheta+1),i)
  SLL_ALLOCATE(g(ntheta+1,nr+1),i)
  SLL_ALLOCATE(fp1(nr+1,ntheta+1),i)
  SLL_ALLOCATE(int_r(ntheta),i)

  SLL_ALLOCATE(phi_ref(nr+1,ntheta+1),i)

  step = 0
  f    = 0.0_f64

  if (fcase==1) then
    do i = 1,nr+1
      r      = rmin+real(i-1,f64)*dr
      f(i,:) = 1.0_f64/(1.0_f64*sqrt(2.0_f64*sll_pi)) * &
        exp(-(r-(real(rmax-rmin)/2.0_f64))**2/(2.0_f64*1.0_f64**2))
    end do

  else if (fcase==2) then
    do i = 1,nr+1
      r = rmin+real(i-1,f64)*dr
      if (r>=r1 .and. r<=r2) then
        do j = 1,ntheta+1
          theta  = real(j-1,f64)*dtheta
          f(i,j) = 1._f64+alpha*cos(mode*theta)
        end do
      end if
      do j = 1,ntheta+1
        theta = real(j-1,f64)*dtheta
      enddo
    end do
!f=10._f64

  else if (fcase==3) then
    do i = 1,nr+1
      r = rmin+real(i-1,f64)*dr
      do j = 1,ntheta+1
        theta  = real(j-1,f64)*dtheta
        f(i,j) = -(r-rmin)*(r-rmax)/r**2 * &
          ((36.0_f64-mode**2)*r**4 + &
          (2.0_f64*mode**2-39.0_f64)*r**3*(rmin+rmax) + &
          (9.0_f64-mode**2)*r**2*(rmin**2+rmax**2) + &
          (30.0_f64-4.0_f64*mode**2)*r**2*rmin*rmax + &
          (2.0_f64*mode**2-3.0_f64)*r*rmin*rmax*(rmin+rmax) - &
          mode**2*rmin**2*rmax**2)*cos(mode*theta)
      end do
    end do

  else if (fcase==4) then
    do i = 1,nr+1
      r = rmin+real(i-1,f64)*dr
      do j = 1,ntheta+1
        theta  = real(j-1,f64)*dtheta
        f(i,j) = 1.0_f64/(0.5_f64*sqrt(2.0_f64*sll_pi)) * &
          exp(-(r-(real(rmax-rmin)/2.0_f64))**2/(2*0.5_f64**2)) * &
          (1.0_f64+alpha*sin(mode*theta))
      end do
    end do

  else if (fcase==5) then
    open(25,file=f_file,action="read")
    read(25,*)f
    close(25)

  else if (fcase==6) then
    !essai
    r1 = 4.0_f64
    r2 = 5.0_f64
    do j = 1,ntheta+1
      theta = real(j-1,f64)*dtheta
      do i = 1,nr+1
        r = rmin+real(i-1,f64)*dr
        if (r>=r1 .and. r<=r2) then
          f(i,j)=1._f64+alpha*mode**2*cos(mode*theta)/r
        else
          f(i,j)=alpha*mode**2*cos(mode*theta)/r
        end if
      end do
    end do

  else
    print*,"f is not defined"
    print*,'see variable "fcase" in file selalib/prototype/src/simulation/CG_polar.F90'
    print*,'can not go any further'
    print*,'exiting...'
    stop
  end if

  fp1 = 0.0_f64

  !call solve_poisson_polar(plan_sl%poisson,f,plan_sl%phi)!,ierr_poiss)

  call poisson_solve_polar(plan_sl%poisson,f,plan_sl%phi,ierr_poiss)
  if (ierr_poiss.ne.0) then
    print*, 'WARNING: poisson error is larger than 1.e-12 for ', &
      ierr_poiss, 'values of poloidal modes for step 0'
  end if
  call compute_grad_field(plan_sl%grad,plan_sl%phi, &
    plan_sl%adv%field)

  !write f in a file before calculations
  call print2d(dom,f(1:(nr+1),1:(ntheta+1)),Nr,Ntheta, &
    visu,step,"CG")

  print *,'#bc(1)=',bc(1)
  print *,'#bc(2)=',bc(2)
  print *,'#SLL_DIRICHLET=',SLL_DIRICHLET
  print *,'#SLL_NEUMANN=',SLL_NEUMANN
  print *,'#SLL_NEUMANN_MODE_0=',SLL_NEUMANN_MODE_0
  

  !mode=1.5
  if(bc(1)==SLL_DIRICHLET)then
    k1 = (r1**2-r2**2+2.0_f64*r1**2*log(rmax/r1) + &
      2.0_f64*r2**2*log(r2/rmax))/(4.0_f64*log(rmin/rmax))
    k2 = (r1**2-r2**2+2.0_f64*r1**2*log(rmin/r1) + &
      2.0_f64*r2**2*log(r2/rmax))/(4.0_f64*log(rmin/rmax))
    k3 = (r1**2-r2**2+2.0_f64*r1**2*log(rmin/r1) + &
      2.0_f64*r2**2*log(r2/rmin))/(4.0_f64*log(rmin/rmax))
    c1 = (2.0_f64*r1**2*log(rmax/r1)+2.0_f64*r2**2*log(r2/rmax) + &
      r1**2-r2**2)*log(rmin)/(-4.0_f64*log(rmin/rmax))
    c2 = (2.0_f64*r2**2*log(rmin)*log(r2/rmax) + &
      2.0_f64*r1**2*log(rmax)*log(rmin/r1)+r1**2*log(rmax) - &
      r2**2*log(rmin))/(-4.0_f64*log(rmin/rmax))
    c3 = (r1**2-r2**2+2.0_f64*r2**2*log(r2/rmin) + &
      2.0_f64*r1**2*log(rmin/r1))*log(rmax) / &
      (-4.0_f64*log(rmin/rmax))
  endif
  if((bc(1)==SLL_NEUMANN).or.(bc(1)==SLL_NEUMANN_MODE_0))then
    k1 = 0._f64
    k2 = r1**2/2._f64
    k3 = r1**2/2._f64-r2**2/2._f64
    c1 = r1**2*log(r1)/2._f64+r2**2/4._f64-r2**2*log(r2)/2._f64
    c1 = c1-r1**2*log(rmax)/2._f64+r2**2*log(rmax)/2._f64 - &
      r1**2/4._f64
    c2 = r2**2/4._f64-r2**2*log(r2)/2._f64-r1**2 * &
      log(rmax)/2._f64+r2**2*log(rmax)/2._f64
    c3 = -log(rmax)*(r1**2-r2**2)/2._f64
  endif

  if (mode==0) then 
    if(bc(1)==SLL_DIRICHLET)then
      k1_mode = (r1**2-r2**2+2.0_f64*r1**2*log(rmax/r1) + &
        2.0_f64*r2**2*log(r2/rmax))/(4.0_f64*log(rmin/rmax))
      k2_mode = (r1**2-r2**2+2.0_f64*r1**2*log(rmin/r1) + &
        2.0_f64*r2**2*log(r2/rmax))/(4.0_f64*log(rmin/rmax))
      k3_mode = (r1**2-r2**2+2.0_f64*r1**2*log(rmin/r1) + &
        2.0_f64*r2**2*log(r2/rmin))/(4.0_f64*log(rmin/rmax))
      c1_mode = (2.0_f64*r1**2*log(rmax/r1) + &
        2.0_f64*r2**2*log(r2/rmax)+r1**2-r2**2)*log(rmin) / &
        (-4.0_f64*log(rmin/rmax))
      c2_mode = (2.0_f64*r2**2*log(rmin)*log(r2/rmax) + &
        2.0_f64*r1**2*log(rmax)*log(rmin/r1) + &
        r1**2*log(rmax)-r2**2*log(rmin))/(-4.0_f64*log(rmin/rmax))
      c3_mode = (r1**2-r2**2+2.0_f64*r2**2*log(r2/rmin) + &
        2.0_f64*r1**2*log(rmin/r1))*log(rmax) / &
        (-4.0_f64*log(rmin/rmax))
    endif
    if((bc(1)==SLL_NEUMANN).or.(bc(1)==SLL_NEUMANN_MODE_0))then
      k1_mode = 0._f64
      k2_mode = r1**2/2._f64
      k3_mode = r1**2/2._f64-r2**2/2._f64
      c1_mode = r1**2*log(r1)/2._f64+r2**2/4._f64-r2**2 * &
        log(r2)/2._f64
      c1_mode = c1_mode-r1**2*log(rmax)/2._f64+r2**2 * &
        log(rmax)/2._f64-r1**2/4._f64
      c2_mode = r2**2/4._f64-r2**2*log(r2)/2._f64-r1**2 * &
        log(rmax)/2._f64+r2**2*log(rmax)/2._f64
      c3_mode = -log(rmax)*(r1**2-r2**2)/2._f64
    endif
  endif

  if(mode==1)then
    if((bc(1)==SLL_NEUMANN))then
      k1_mode = (-3._f64*r1*rmax**2-r2**3 + &
        3*rmax**2*r2+r1**3)/(6._f64*(rmin**2+rmax**2))
      k2_mode = (3._f64*r2*rmax**2-r2**3 + &
        3*rmin**2*r1+r1**3)/(6._f64*(rmin**2+rmax**2))
      k3_mode = (-3._f64*r2*rmin**2-r2**3 + &
        3*rmin**2*r1+r1**3)/(6._f64*(rmin**2+rmax**2))
      c1_mode = (-3._f64*r1*rmax**2-r2**3 + &
        3*rmax**2*r2+r1**3)*rmin**2/(6._f64*(rmin**2+rmax**2))
      c2_mode = -(3._f64*r1*rmin**2*rmax**2 + &
        r2**3*rmin**2-3*rmin**2*rmax**2*r2+rmax**2*r1**3) / &
        (6._f64 *(rmin**2+rmax**2))
      c3_mode = -(-3._f64*r2*rmin**2-r2**3 + &
        3*rmin**2*r1+r1**3)*rmax**2/(6._f64*(rmin**2+rmax**2))
    endif

    if((bc(1)==SLL_DIRICHLET).or.(bc(1)==SLL_NEUMANN_MODE_0))then
      k1_mode = (-3._f64*r1*rmax**2-r2**3 + &
        3*rmax**2*r2+r1**3)/(6._f64*(-rmin**2+rmax**2))
      k2_mode = (3._f64*r2*rmax**2-r2**3 - &
        3*rmin**2*r1+r1**3)/(6._f64*(-rmin**2+rmax**2))
      k3_mode = (3._f64*r2*rmin**2-r2**3 - &
        3*rmin**2*r1+r1**3)/(6._f64*(-rmin**2+rmax**2))
      c1_mode = (-3._f64*r1*rmax**2-r2**3 + &
        3*rmax**2*r2+r1**3)*rmin**2/(6._f64*(rmin**2-rmax**2))
      c2_mode = (-3._f64*r1*rmin**2*rmax**2 - &
        r2**3*rmin**2+3*rmin**2*rmax**2*r2+rmax**2*r1**3) / &
        (6._f64 *(rmin**2-rmax**2))
      c3_mode = (-3._f64*r1*rmin**2-r2**3 + &
        3*rmin**2*r2+r1**3)*rmax**2/(6._f64*(rmin**2-rmax**2))
    endif
  endif

  if(mode==3)then
    if((bc(1)==SLL_DIRICHLET).or.(bc(1)==SLL_NEUMANN_MODE_0))then
      k1_mode = (r1*r2*(r2**5-r1**5) - &
        5._f64*rmax**6*(r2-r1))/(30._f64*r2*r1*(rmin**6-rmax**6))
      k2_mode = (r1*r2*(r2**5-r1**5) - &
        5._f64*(rmin**6*r2-rmax**6*r1))/(30._f64*r2*r1 * &
        (rmin**6-rmax**6))    
      k3_mode = (-r1*r2*(r1**5-r2**5) - &
        5._f64*rmin**6*(r2-r1))/(30._f64*r2*r1*(rmin**6-rmax**6))
      c1_mode = (-r1*r2*(r2**5-r1**5) + &
        5._f64*rmax**6*(r2-r1))/(30._f64*r2*r1*(rmin**6-rmax**6))
      c2_mode = (-r1*r2*(rmin**6*r2**5-rmax**6*r1**5) + &
        5._f64*(rmin*rmax)**6*(r2-r1)) / &
        (30._f64*r2*r1*(rmin**6-rmax**6))
      c3_mode = rmax**6*(-r1*r2*(r2**5-r1**5) + &
        5._f64*rmin**6*(r2-r1))/(30._f64*r2*r1*(rmin**6-rmax**6))
    endif
  endif

  if(mode==7)then
    if((bc(1)==SLL_DIRICHLET).or.(bc(1)==SLL_NEUMANN_MODE_0))then
      k1_mode = -5._f64*r1**5*r2**14 + &
        5._f64*r1**14*r2**5-9._f64*r1**5*rmax**14 + &
        9._f64*r2**5*rmax**14
      k1_mode = k1/(630._f64*r1**5*r2**5*(rmin**14+rmax**14))
      k2_mode = -5._f64*r1**5*r2**14 + &
        5._f64*r1**14*r2**5-9._f64*r2**5*rmin**14 - &
        9._f64*r1**5*rmax**14
      k2_mode = k2/(630._f64*r1**5*r2**5*(rmin**14+rmax**14))
      k3_mode = -5._f64*r1**5*r2**14 + &
        5._f64*r1**14*r2**5-9._f64*r2**5*rmin**14 + &
        9._f64*r1**5*rmin**14
      k3_mode = k3/(630._f64*r1**5*r2**5*(rmin**14+rmax**14))
      c1_mode = -5._f64*r1**5*r2**14 + &
        5._f64*r1**14*r2**5-9._f64*r1**5*rmax**14 + &
        9._f64*r2**5*rmax**14
      c1_mode = rmin**14*c1 / &
        (630._f64*r1**5*r2**5*(rmin**14+rmax**14))
      c2_mode = 5._f64*rmin**14*r1**5*r2**14 + &
        9._f64*rmax**14*rmin**14*r1**5
      c2_mode = c2_mode-9._f64*rmax**14*r2**5*rmin**14 + &
        5._f64*rmax**14*r1**14*r2**5
      c2_mode = -c2_mode/(630._f64*r1**5*r2**5*(rmin**14+rmax**14))
      c3_mode = -5._f64*r1**5*r2**14 + &
        5._f64*r1**14*r2**5+9._f64*r1**5*rmin**14 - &
        9._f64*r2**5*rmin**14
      c3_mode = -rmax**14*c3_mode / &
        (630._f64*r1**5*r2**5*(rmin**14+rmax**14))
    endif
  endif

  tmp = 0._f64
  l1  = 0.0_f64
  l2  = 0.0_f64

  open (unit=20,file='CGinit.dat')
  do i = 1,nr+1
    r = rmin+real(i-1,f64)*dr
    if(mode==0)then
      if (r<r1) then
        temps_mode = k1_mode*log(r)+c1_mode
      else if (r>r2) then
        temps_mode = k3_mode*log(r)+c3_mode
      else
        temps_mode = k2_mode*log(r)+c2_mode-r**2/4.0_f64
      end if
    end if

    if((mode>=3).or.(mode==1))then
      if (r<r1) then
        temps_mode = k1_mode*r**mode+c1_mode/r**(mode)
      else if (r>r2) then
        temps_mode = k3_mode*r**mode+c3_mode/r**mode
      else
        temps_mode = k2_mode*r**mode+c2_mode/r**mode + &
          r**2/(mode**2-4._f64)
      end if
    end if


    if (r<r1) then
      temps = k1*log(r)+c1
    else if (r>r2) then
      temps = k3*log(r)+c3
    else
      temps = k2*log(r)+c2-r**2/4.0_f64
    end if

    temps_mode = temps_mode*alpha
    do j = 1,ntheta+1
      theta   = real(j-1,f64)*dtheta
      x       = r*cos(theta)
      y       = r*sin(theta)
      err_loc = abs(plan_sl%phi(i,j) - &
        temps_mode*cos(mode*theta)-temps)
      phi_ref(i,j) = temps_mode*cos(mode*theta)+temps

      write(20,*) r, theta, x, y, plan_sl%phi(i,j), &
        temps+temps_mode*cos(mode*theta)
      tmp = max(tmp,err_loc)
      if (i==1 .or. i==nr+1) then
        l1 = l1+err_loc*r/2.0_f64
        l2 = l2+err_loc**2*r/2.0_f64
      else
        l1 = l1+err_loc*r
        l2 = l2+err_loc**2*r
      end if
    end do
    write(20,*)' '
  end do
  close(20)
  l1 = l1*dr*dtheta
  l2 = sqrt(l2*dr*dtheta)
  print*,"#error for phi in initialization", &
    nr,dr,tmp/(1._f64+abs(alpha)), &
    l1/(1._f64+abs(alpha)),l2/(1._f64+abs(alpha))

  open(unit=23,file='thdiag.dat')
  write(23,*)'#fcase',fcase,'scheme',scheme, &
    'mode',mode,'grad',grad,'carac',carac
  write(23,*)'#nr',nr,'ntheta',ntheta,'alpha',alpha
  write(23,*)'#tf = ',tf,'  nb_step = ',nb_step,'  dt = ',dt
  write(23,*)'#   t   //   w   //   l1 rel  //   l2  rel //   e   //   re   //   im'

  do i = 1,nr+1
    r = rmin+real(i-1,f64)*dr
    plan_sl%adv%field(2,i,:) = plan_sl%adv%field(2,i,:)/r
  end do

  w0    = 0.0_f64
  l10   = 0.0_f64
  l20   = 0.0_f64
  e0    = 0.0_f64
  int_r = 0.0_f64
  do j = 1,ntheta
    w0  = w0+(f(1,j)*rmin+f(nr+1,j)*rmax)/2.0_f64
    l10 = l10+abs(f(1,j)*rmin+f(nr+1,j)*rmax)/2.0_f64
    l20 = l20+(f(1,j)/2.0_f64)**2*rmin+(f(nr+1,j)/2.0_f64)**2*rmax
    e0  = e0+rmin*(plan_sl%adv%field(1,1,j))**2/2.0_f64 + &
      rmax*(plan_sl%adv%field(1,nr+1,j))**2/2.0_f64 + &
      rmin*(plan_sl%adv%field(2,1,j))**2/2.0_f64 + &
      rmax*(plan_sl%adv%field(2,nr+1,j))**2/2.0_f64
    int_r(j) = (f(1,j)*rmin+f(nr+1,j)*rmax)/2.0_f64
    do i = 2,nr
      r   = rmin+real(i-1,f64)*dr
      w0  = w0+r*f(i,j)
      !w0  = w0+f(i,j)
      l10 = l10+r*abs(f(i,j))
      l20 = l20+r*f(i,j)**2
      e0  = e0+r*(plan_sl%adv%field(1,i,j)**2 + &
        plan_sl%adv%field(2,i,j)**2)
      int_r(j) = int_r(j)+f(i,j)*r
    end do
  end do

  w0    = w0*dr*dtheta
  l10   = l10*dr*dtheta
  l20   = sqrt(l20*dr*dtheta)
  e0    = e0*dr*dtheta/2.0_f64
  int_r = int_r*dr
  call fft_apply_plan(plan_sl%poisson%pfwd,int_r,int_r)
  do ii=1,8
    time_mode(ii) = real(fft_get_mode(plan_sl%poisson%pfwd,int_r,ii-1))**2 &
     &+aimag(fft_get_mode(plan_sl%poisson%pfwd,int_r,ii-1))**2
  enddo
  !mode_slope=time_mode
  write(23,*)'#t=0',w0,l10,l20,e0
  write(23,*)0.0_f64,w0,1.0_f64,1.0_f64,0.0_f64,e0, &
     time_mode(1:8)!,time_mode(1:8)
!    real(fft_get_mode(plan_sl%poisson%pfwd,int_r,0)), &
!    aimag(fft_get_mode(plan_sl%poisson%pfwd,int_r,0)), &
!    real(fft_get_mode(plan_sl%poisson%pfwd,int_r,1)), &
!    aimag(fft_get_mode(plan_sl%poisson%pfwd,int_r,1)), &
!    real(fft_get_mode(plan_sl%poisson%pfwd,int_r,2)), &
!    aimag(fft_get_mode(plan_sl%poisson%pfwd,int_r,2)), &
!    real(fft_get_mode(plan_sl%poisson%pfwd,int_r,3)), &
!    aimag(fft_get_mode(plan_sl%poisson%pfwd,int_r,3)), &
!    real(fft_get_mode(plan_sl%poisson%pfwd,int_r,4)), &
!    aimag(fft_get_mode(plan_sl%poisson%pfwd,int_r,4)), &
!    real(fft_get_mode(plan_sl%poisson%pfwd,int_r,7)), &
!    aimag(fft_get_mode(plan_sl%poisson%pfwd,int_r,7)), &
!    real(fft_get_mode(plan_sl%poisson%pfwd,int_r,ntheta-1)), &
!    aimag(fft_get_mode(plan_sl%poisson%pfwd,int_r,ntheta-1)), &
!    real(fft_get_mode(plan_sl%poisson%pfwd,int_r,ntheta-2)), &
!    aimag(fft_get_mode(plan_sl%poisson%pfwd,int_r,ntheta-2)), &
!    real(fft_get_mode(plan_sl%poisson%pfwd,int_r,ntheta-3)), &
!    aimag(fft_get_mode(plan_sl%poisson%pfwd,int_r,ntheta-3)), &
!    real(fft_get_mode(plan_sl%poisson%pfwd,int_r,ntheta-4)), &
!    aimag(fft_get_mode(plan_sl%poisson%pfwd,int_r,ntheta-4)), &
!    real(fft_get_mode(plan_sl%poisson%pfwd,int_r,ntheta-7)), &
!    aimag(fft_get_mode(plan_sl%poisson%pfwd,int_r,ntheta-7))

  




  call sll_set_time_mark(t1)

  do step = 1,nb_step
    if (step==101) then
      call sll_set_time_mark(t2)
      temps = sll_time_elapsed_between(t1,t2)
      temps = temps/100*real(nb_step,f32)
      hh    = floor(temps/3600.0d0)
      min   = floor((temps-3600.0d0*real(hh))/60.0d0)
      ss    = floor(temps-3600.0d0*real(hh)-60.0d0*real(min))
      print*,'# CPU time estimation : ',hh,'h',min,'min',ss,'s'
    end if

    if(scheme==0)then
      call poisson_solve_polar(plan_sl%poisson,f,plan_sl%phi)
      call advect_CG_polar2(plan_sl%adv,f,fp1,plan_sl%phi)
    endif

    if (scheme==1) then
      !classical semi-Lagrangian scheme (order 1)
      call SL_classic(plan_sl,f,fp1)

    else if (scheme==10) then
      !semi-Lagrangian scheme remap (Euler, order 1)
      call SL_remap(plan_sl,f,fp1,interp_case,PPM_order)

    else if (scheme==2) then
      !semi-Lagrangian predictive-corrective scheme
      call SL_ordre_2(plan_sl,f,fp1)

    else if (scheme==20) then
      !semi-Lagrangian scheme remap (predictive-corrective, order 2)
      call SL_remap_ordre_2(plan_sl,f,fp1,interp_case,PPM_order)

    else if (scheme==3) then
      !leap-frog scheme
      if (step==1) then
        call SL_ordre_2(plan_sl,f,fp1)
        plan_sl%adv%dt=2.0_f64*dt
      else 
        call poisson_solve_polar(plan_sl%poisson,f,plan_sl%phi)
        call compute_grad_field(plan_sl%grad,plan_sl%phi, &
          plan_sl%adv%field)
        call advect_CG_polar(plan_sl%adv,g,fp1)
      end if
      g=f

    else
      print*,'no scheme define'
      print*,"the program won't do anything"
      print*,'see variable scheme in file selalib/prototype/src/simulation to solve the probleme'
      print*,'exiting the loop'
      exit
    end if

    !---> periodicity
    fp1(:,ntheta+1) = fp1(:,1)
    f = fp1

    !---> poisson solving
    call poisson_solve_polar(plan_sl%poisson,f,plan_sl%phi,ierr_poiss)
    if (ierr_poiss.ne.0) then
      print*, 'WARNING: poisson error is larger than 1.e-12 for ', &
        ierr_poiss, 'values of poloidal modes for step ', step
    end if

    !---> Computation of the field
    call compute_grad_field(plan_sl%grad,plan_sl%phi, &
      plan_sl%adv%field)

    do i = 1,nr+1
      r = rmin+real(i-1,f64)*dr
      plan_sl%adv%field(2,i,:) = plan_sl%adv%field(2,i,:)/r
    end do
    w     = 0.0_f64
    l1    = 0.0_f64
    l2    = 0.0_f64
    e     = 0.0_f64
    int_r = 0.0_f64
    do j = 1,ntheta
      w  = w+(f(1,j)*rmin+f(nr+1,j)*rmax)/2.0_f64
      l1 = l1+abs(f(1,j)*rmin+f(nr+1,j)*rmax)/2.0_f64
      l2 = l2+(f(1,j)/2.0_f64)**2*rmin+(f(nr+1,j)/2.0_f64)**2*rmax
      e  = e+rmin*(plan_sl%adv%field(1,1,j))**2/2.0_f64 + &
        rmax*(plan_sl%adv%field(1,nr+1,j))**2/2.0_f64 + &
        rmin*(plan_sl%adv%field(2,1,j))**2/2.0_f64 + &
        rmax*(plan_sl%adv%field(2,nr+1,j))**2/2.0_f64
      int_r(j) = (f(1,j)*rmin+f(nr+1,j)*rmax)/2.0_f64
      int_r(j) = (plan_sl%phi(1,j)*rmin + &
        plan_sl%phi(nr+1,j)*rmax)/2.0_f64
      do i = 2,nr
        r  = rmin+real(i-1,f64)*dr
        w  = w+r*f(i,j)
        !w  = w+f(i,j)
        l1 = l1+r*abs(f(i,j))
        l2 = l2+r*f(i,j)**2
        e  = e+r*(plan_sl%adv%field(1,i,j)**2 + &
          plan_sl%adv%field(2,i,j)**2)
        int_r(j) = int_r(j)+plan_sl%phi(i,j)*r
      end do
    end do
    w     = w*dr*dtheta
    l1    = l1*dr*dtheta
    l2    = sqrt(l2*dr*dtheta)
    e     = e*dr*dtheta/2.0_f64
    int_r = int_r*dr
    call fft_apply_plan(plan_sl%poisson%pfwd,int_r,int_r)
    do ii=1,8
      mode_slope(ii) = time_mode(ii)
      time_mode(ii) = real(fft_get_mode(plan_sl%poisson%pfwd,int_r,ii-1))**2 &
       &+aimag(fft_get_mode(plan_sl%poisson%pfwd,int_r,ii-1))**2
       mode_slope(ii) = (log(time_mode(ii))-log(mode_slope(ii)))/dt
    enddo
    !mode_slope=(time_mode-mode_slope)/dt
    
    
    write(23,*) dt*real(step,f64),w,l1/l10,l2/l20,e-e0,e, &
      time_mode,mode_slope
!      real(fft_get_mode(plan_sl%poisson%pfwd,int_r,0)), &   !$7
!      aimag(fft_get_mode(plan_sl%poisson%pfwd,int_r,0)), &  !$8
!      real(fft_get_mode(plan_sl%poisson%pfwd,int_r,1)), &   !$9
!      aimag(fft_get_mode(plan_sl%poisson%pfwd,int_r,1)), &  !$10
!      real(fft_get_mode(plan_sl%poisson%pfwd,int_r,2)), &   !$11
!      aimag(fft_get_mode(plan_sl%poisson%pfwd,int_r,2)), &  !$12
!      real(fft_get_mode(plan_sl%poisson%pfwd,int_r,3)), &   !$13
!      aimag(fft_get_mode(plan_sl%poisson%pfwd,int_r,3)), &  !$14
!      real(fft_get_mode(plan_sl%poisson%pfwd,int_r,4)), &   !$15
!      aimag(fft_get_mode(plan_sl%poisson%pfwd,int_r,4)), &  !$16
!      real(fft_get_mode(plan_sl%poisson%pfwd,int_r,7)), &   !$17
!      aimag(fft_get_mode(plan_sl%poisson%pfwd,int_r,7)), &  !$18
!      real(fft_get_mode(plan_sl%poisson%pfwd,int_r,ntheta-1)), & 
!      aimag(fft_get_mode(plan_sl%poisson%pfwd,int_r,ntheta-1)), &
!      real(fft_get_mode(plan_sl%poisson%pfwd,int_r,ntheta-2)), &
!      aimag(fft_get_mode(plan_sl%poisson%pfwd,int_r,ntheta-2)), &
!      real(fft_get_mode(plan_sl%poisson%pfwd,int_r,ntheta-3)), &
!      aimag(fft_get_mode(plan_sl%poisson%pfwd,int_r,ntheta-3)), &
!      real(fft_get_mode(plan_sl%poisson%pfwd,int_r,ntheta-4)), &
!      aimag(fft_get_mode(plan_sl%poisson%pfwd,int_r,ntheta-4)), &
!      real(fft_get_mode(plan_sl%poisson%pfwd,int_r,ntheta-7)), &
!      aimag(fft_get_mode(plan_sl%poisson%pfwd,int_r,ntheta-7))

    if ((step/500)*500==step) then
      print*,'#step',step
    end if

#ifndef NOHDF5
    if (step==1 .or. step/visustep*visustep==step) then
      call plot_f(step/visustep)
    end if
#endif
  end do
  write(23,*)' '
  write(23,*)' '
  close(23)

  call sll_set_time_mark(t3)
  temps = sll_time_elapsed_between(t1,t3)
  hh    = floor(temps/3600.0d0)
  min   = floor((temps-3600.0d0*real(hh))/60.0d0)
  ss    = floor(temps-3600.0d0*real(hh)-60.0d0*real(min))
  print*,'# CPU time for the loop in time : ', &
    hh,'h',min,'min',ss,'s'


  !---> write the final f in a file
  call print2d(dom,f(1:(nr+1),1:(ntheta+1)) , &
    Nr,Ntheta,visu,step,"CG")
  open (unit=21,file='CGrestart.dat')
  write(21,*)f
  close(21)

  SLL_DEALLOCATE_ARRAY(div,i)
  SLL_DEALLOCATE_ARRAY(f,i)
  SLL_DEALLOCATE_ARRAY(fp1,i)
  SLL_DEALLOCATE_ARRAY(g,i)
  call delete_SL_polar(plan_sl)


#ifndef NOHDF5
!*********************
contains
!*********************

  !---------------------------------------------------
  ! Save the mesh structure
  !---------------------------------------------------
  subroutine plot_f(iplot)
    use sll_xdmf
    use sll_hdf5_io
    sll_int32 :: file_id
    sll_int32 :: error
    sll_real64, dimension(:,:), allocatable :: x1
    sll_real64, dimension(:,:), allocatable :: x2
    sll_int32 :: i, j
    sll_int32, intent(in) :: iplot
    character(len=4)      :: cplot
    sll_int32             :: nnodes_x1, nnodes_x2

    nnodes_x1 = nr+1
    nnodes_x2 = ntheta+1

    if (iplot == 1) then

      SLL_ALLOCATE(x1(nnodes_x1,nnodes_x2), error)
      SLL_ALLOCATE(x2(nnodes_x1,nnodes_x2), error)
      do j = 1,nnodes_x2
        do i = 1,nnodes_x1
          r       = rmin+real(i-1,f32)*dr
          theta   = real(j-1,f32)*dtheta
          x1(i,j) = r*cos(theta)
          x2(i,j) = r*sin(theta)
        end do
      end do
      call sll_hdf5_file_create("polar_mesh-x1.h5",file_id,error)
      call sll_hdf5_write_array(file_id,x1,"/x1",error)
      call sll_hdf5_file_close(file_id, error)
      call sll_hdf5_file_create("polar_mesh-x2.h5",file_id,error)
      call sll_hdf5_write_array(file_id,x2,"/x2",error)
      call sll_hdf5_file_close(file_id, error)
      deallocate(x1)
      deallocate(x2)

    end if

    call int2string(iplot,cplot)
    call sll_xdmf_open("f"//cplot//".xmf","polar_mesh",nnodes_x1,nnodes_x2,file_id,error)
    call sll_xdmf_write_array("f"//cplot,f,"values",error,file_id,"Node")
    call sll_xdmf_close(file_id,error)
  end subroutine plot_f

#endif
 
end program cg_polar
