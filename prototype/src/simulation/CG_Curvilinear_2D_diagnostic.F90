module module_diagnostic

#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "selalib.h"

   use module_cg_curvi_structure
   
contains


subroutine diagnostic_1(f,plan_sl,phi_ref,int_r,bc,rmin,rmax,dr,dtheta,nr &
       & ,ntheta,nb_step,fcase,scheme,carac,grad,mode,&
       & l10,l20,e0,dt,alpha,r1,r2)

implicit none

type(sll_SL_curvilinear)     , pointer          :: plan_sl
sll_real64, dimension (:,:)  , intent(inout)    :: phi_ref
sll_real64, dimension (:,:)  , intent(in)       :: f
sll_real64, dimension (:)    , intent(inout)    :: int_r
sll_int32,  intent(in)                          :: bc(2)
sll_real64, intent(in)                          :: r1,r2,rmax,rmin,dr,dtheta,alpha,dt
sll_int32,  intent(in)                          :: fcase, scheme, carac, grad,nr, ntheta
sll_real64, intent(out)                         :: l10,l20,e0

sll_real64 :: r,theta,x,y
sll_real64 :: k1,k2,k3,c1,c2,c3 
sll_real64 :: k1_mode,k2_mode,k3_mode,c1_mode,c2_mode,c3_mode
sll_int32  :: NEUMANN,NEUMANN_MODE0,DIRICHLET
sll_real64 :: w0,w, l1, l2, e, re, im,tmp,temps
sll_int32  :: mode,nb_step,i,j
sll_real64 :: temps_mode,err_loc


 !mode=1.5
DIRICHLET=1
NEUMANN=2
NEUMANN_MODE0=3

  if(bc(1)==DIRICHLET)then
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
  if((bc(1)==NEUMANN).or.(bc(1)==NEUMANN_MODE0))then
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
    if(bc(1)==DIRICHLET)then
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
    if((bc(1)==NEUMANN).or.(bc(1)==NEUMANN_MODE0))then
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
    if((bc(1)==NEUMANN))then
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

    if((bc(1)==DIRICHLET).or.(bc(1)==NEUMANN_MODE0))then
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
    if((bc(1)==DIRICHLET).or.(bc(1)==NEUMANN_MODE0))then
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
    if((bc(1)==DIRICHLET).or.(bc(1)==NEUMANN_MODE0))then
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
    l1/(1._f64+abs(alpha)),l2/(1._f64+abs(alpha)),alpha,tmp,l1

  open(unit=23,file='diagMF.dat')
  write(23,*)'#fcase',fcase,'scheme',scheme, &
    'mode',mode,'grad',grad,'carac',carac
  write(23,*)'#nr',nr,'ntheta',ntheta,'alpha',alpha
  write(23,*)'#nb_step = ',nb_step,'  dt = ',dt
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
  write(23,*)'#t=0',w0,l10,l20,e0
  write(23,*)0.0_f64,w0,1.0_f64,1.0_f64,0.0_f64,e0, &
    real(fft_get_mode(plan_sl%poisson%pfwd,int_r,0)), &
    aimag(fft_get_mode(plan_sl%poisson%pfwd,int_r,0)), &
    real(fft_get_mode(plan_sl%poisson%pfwd,int_r,1)), &
    aimag(fft_get_mode(plan_sl%poisson%pfwd,int_r,1)), &
    real(fft_get_mode(plan_sl%poisson%pfwd,int_r,2)), &
    aimag(fft_get_mode(plan_sl%poisson%pfwd,int_r,2)), &
    real(fft_get_mode(plan_sl%poisson%pfwd,int_r,3)), &
    aimag(fft_get_mode(plan_sl%poisson%pfwd,int_r,3)), &
    real(fft_get_mode(plan_sl%poisson%pfwd,int_r,4)), &
    aimag(fft_get_mode(plan_sl%poisson%pfwd,int_r,4)), &
    real(fft_get_mode(plan_sl%poisson%pfwd,int_r,7)), &
    aimag(fft_get_mode(plan_sl%poisson%pfwd,int_r,7)), &
    real(fft_get_mode(plan_sl%poisson%pfwd,int_r,ntheta-1)), &
    aimag(fft_get_mode(plan_sl%poisson%pfwd,int_r,ntheta-1)), &
    real(fft_get_mode(plan_sl%poisson%pfwd,int_r,ntheta-2)), &
    aimag(fft_get_mode(plan_sl%poisson%pfwd,int_r,ntheta-2)), &
    real(fft_get_mode(plan_sl%poisson%pfwd,int_r,ntheta-3)), &
    aimag(fft_get_mode(plan_sl%poisson%pfwd,int_r,ntheta-3)), &
    real(fft_get_mode(plan_sl%poisson%pfwd,int_r,ntheta-4)), &
    aimag(fft_get_mode(plan_sl%poisson%pfwd,int_r,ntheta-4)), &
    real(fft_get_mode(plan_sl%poisson%pfwd,int_r,ntheta-7)), &
    aimag(fft_get_mode(plan_sl%poisson%pfwd,int_r,ntheta-7))

 do i = 1,nr+1
    r = rmin+real(i-1,f64)*dr
    plan_sl%adv%field(2,i,:) = plan_sl%adv%field(2,i,:)*r
 end do
end subroutine diagnostic_1

!!******************
subroutine diagnostic_2(f,plan_sl,phi_ref,int_r,rmin,rmax,dr,dtheta,nr &
       & ,ntheta,step,l10,l20,e0,mode,dt,alpha)
implicit none

type(sll_SL_curvilinear)     ,  pointer       :: plan_sl
sll_real64, dimension (:,:)  ,  intent(inout) :: phi_ref
sll_real64, dimension (:,:)  ,  intent(in)    :: f
sll_real64, dimension (:)    ,  intent(inout) :: int_r
sll_real64, intent(in)                        :: dr,dtheta,alpha,dt,rmax,rmin
sll_int32 , intent(in)                        :: nr,ntheta,step,mode
sll_real64   :: r,theta
sll_real64   :: w, l10, l1, l20, l2, e, e0, re
sll_int32    :: i,j


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
    write(23,*) dt*real(step,f64),w,l1/l10,l2/l20,e-e0,e, &
      real(fft_get_mode(plan_sl%poisson%pfwd,int_r,0)), &   !$7
      aimag(fft_get_mode(plan_sl%poisson%pfwd,int_r,0)), &  !$8
      real(fft_get_mode(plan_sl%poisson%pfwd,int_r,1)), &   !$9
      aimag(fft_get_mode(plan_sl%poisson%pfwd,int_r,1)), &  !$10
      real(fft_get_mode(plan_sl%poisson%pfwd,int_r,2)), &   !$11
      aimag(fft_get_mode(plan_sl%poisson%pfwd,int_r,2)), &  !$12
      real(fft_get_mode(plan_sl%poisson%pfwd,int_r,3)), &   !$13
      aimag(fft_get_mode(plan_sl%poisson%pfwd,int_r,3)), &  !$14
      real(fft_get_mode(plan_sl%poisson%pfwd,int_r,4)), &   !$15
      aimag(fft_get_mode(plan_sl%poisson%pfwd,int_r,4)), &  !$16
      real(fft_get_mode(plan_sl%poisson%pfwd,int_r,7)), &   !$17
      aimag(fft_get_mode(plan_sl%poisson%pfwd,int_r,7)), &  !$18
      real(fft_get_mode(plan_sl%poisson%pfwd,int_r,ntheta-1)), & 
      aimag(fft_get_mode(plan_sl%poisson%pfwd,int_r,ntheta-1)), &
      real(fft_get_mode(plan_sl%poisson%pfwd,int_r,ntheta-2)), &
      aimag(fft_get_mode(plan_sl%poisson%pfwd,int_r,ntheta-2)), &
      real(fft_get_mode(plan_sl%poisson%pfwd,int_r,ntheta-3)), &
      aimag(fft_get_mode(plan_sl%poisson%pfwd,int_r,ntheta-3)), &
      real(fft_get_mode(plan_sl%poisson%pfwd,int_r,ntheta-4)), &
      aimag(fft_get_mode(plan_sl%poisson%pfwd,int_r,ntheta-4)), &
      real(fft_get_mode(plan_sl%poisson%pfwd,int_r,ntheta-7)), &
      aimag(fft_get_mode(plan_sl%poisson%pfwd,int_r,ntheta-7))

do i = 1,nr+1
    r = rmin+real(i-1,f64)*dr
    plan_sl%adv%field(2,i,:) = plan_sl%adv%field(2,i,:)*r
  end do
end subroutine diagnostic_2
!!******************************************************************************

end module module_diagnostic
