!**************************************************************
!  Copyright INRIA
!  Authors : 
!     CALVI project team
!  
!  This code SeLaLib (for Semi-Lagrangian-Library) 
!  is a parallel library for simulating the plasma turbulence 
!  in a tokamak.
!  
!  This software is governed by the CeCILL-B license 
!  under French law and abiding by the rules of distribution 
!  of free software.  You can  use, modify and redistribute 
!  the software under the terms of the CeCILL-B license as 
!  circulated by CEA, CNRS and INRIA at the following URL
!  "http://www.cecill.info". 
!**************************************************************
!> @author
!> Michel Mehrenberger (mehrenbe@math.unistra.fr)
!> Nicolas Crouseilles
!> Pierre Glanc
!> Eric Madaule
!> @brief 
!> basis for the advection in CG_polar and vp4d_dk
!> and CG_polar
!> should be obsolete and replaced/translated in advection, characteristics modules

module polar_advection
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"

  use polar_operators
  use sll_poisson_2d_polar
  use sll_constants
  use sll_cubic_splines
  implicit none

  !>type sll_plan_adv_polar
  !>type for advection with center-guide equations
  !>the field and other needed data/object are within
  type sll_plan_adv_polar
     sll_real64 :: rmin,rmax,dr,dtheta,dt
     sll_int32 :: nr,ntheta
     type(sll_cubic_spline_2D), pointer :: spl_f
     sll_int32 :: time_scheme
     sll_real64, dimension(:,:,:), pointer :: field,carac
  end type sll_plan_adv_polar

  !>type sll_SL_polar
  !>type for semi Lagrangian
  !>contains other types for the routines called in SL routines
  type sll_SL_polar
     type(sll_plan_adv_polar), pointer :: adv
     type(plan_polar_op), pointer :: grad
     type(sll_plan_poisson_polar), pointer :: poisson
     sll_real64, dimension(:,:), pointer :: phi
  end type sll_SL_polar

contains

!===================================
!  creation of sll_plan_adv_polar
!===================================

  !>function new_plan_adv_polar(rmin,rmax,dr,dtheta,dt,nr,ntheta,time_scheme)
  !>rmin : interior adius
  !>dr, dtheta and dt : size of step in direction r and theta and in time
  !>nr and ntheta : number of space in direction r and theta
  !>time_scheme : integer to choose the scheme for advection
  !>                   1 : using explicit Euler method
  !>                   2 : rotation, rotation speed = -1
  !>                   3 : using symplectic Euler with linear interpolation
  !>                   4 : using symplectic Verlet with linear interpolation
  !>                   5 : using fixed point method
  !>                   6 : using modified symplectic Euler
  !>                   7 : using modified symplectic Verlet
  !>                   8 : using modified fixed point
  function new_plan_adv_polar(rmin,rmax,dr,dtheta,dt,nr,ntheta,time_scheme) result(this)

    type(sll_plan_adv_polar), pointer :: this
    sll_real64, intent(in) :: rmin,rmax,dr,dtheta,dt
    sll_int32, intent(in) :: nr,ntheta
    sll_int32, intent(in) :: time_scheme

    sll_int32 :: err

    SLL_ALLOCATE(this,err)
    SLL_ALLOCATE(this%field(2,nr+1,ntheta+1),err)
    SLL_ALLOCATE(this%carac(2,nr+1,ntheta+1),err)

    this%field=0.0_f64
    this%rmin=rmin
    this%rmax=rmax
    this%dr=dr
    this%dtheta=dtheta
    this%dt=dt
    this%nr=nr
    this%ntheta=ntheta
    this%time_scheme=time_scheme

    this%spl_f => new_spline_2D(nr+1,ntheta+1,rmin,rmax,0._f64, 2._f64*sll_pi, &
         & SLL_HERMITE, SLL_PERIODIC,const_slope_x1_min = 0._f64,const_slope_x1_max = 0._f64)

   ! this%spl_f => new_spline_2D(nr+1,ntheta+1,rmin,rmax,0._f64, 2._f64*sll_pi, &
   !      & SLL_PERIODIC, SLL_PERIODIC)

  end function new_plan_adv_polar

!===================================
!  deletion of sll_plan_adv_polar
!===================================

  !>delete_plan_adv_polar(this)
  !>deletion of sll_plan_adv_polar object
  subroutine delete_plan_adv_polar(this)

    implicit none

    type(sll_plan_adv_polar), pointer :: this

    sll_int32 :: err

    if (associated(this)) then
       SLL_DEALLOCATE_ARRAY(this%field,err)
       call delete_spline_2d(this%spl_f)
       SLL_DEALLOCATE(this,err)
       this=>null()
    end if

  end subroutine delete_plan_adv_polar

!==================================
!  creation of sll_SL_polar type
!==================================

  !>function new_SL(rmin,dr,dtheta,dt,nr,ntheta,grad_case,time_scheme)
  !>creation of sll_SL_polar object for semi Lagrangian scheme in polar coordinates
  !>rmin : interior adius
  !>dr, dtheta and dt : size of step in direction r and theta and in time
  !>nr and ntheta : number of space in direction r and theta
  !>grad_case : sll_int32, see function new_polar_op
  !>time_scheme : sll_int32, see function new_plan_adv_polar
  function new_SL(rmin,rmax,dr,dtheta,dt,nr,ntheta,grad_case,time_scheme,bc) result(this)

    type(sll_SL_polar), pointer :: this
    sll_real64 :: rmin,rmax,dr,dtheta,dt
    sll_int32 :: nr,ntheta,bc(2)
    sll_int32, intent(in) :: grad_case,time_scheme

    sll_int32 :: err

    SLL_ALLOCATE(this,err)
    SLL_ALLOCATE(this%phi(nr+1,ntheta+1),err)

    this%poisson => new_plan_poisson_polar(dr,rmin,nr,ntheta,bc)
    this%grad => new_polar_op(rmin,rmax,dr,dtheta,nr,ntheta,grad_case)
    this%adv => new_plan_adv_polar(rmin,rmax,dr,dtheta,dt,nr,ntheta,time_scheme)

  end function new_SL

!=============================
!  deletion of sll_SL_polar
!=============================

  !>subroutine delete_SL_polar(this)
  !>deletion of sll_SL_polar object
  subroutine delete_SL_polar(this)

    implicit none

    type(sll_SL_polar), pointer :: this

    sll_int32 :: err

    if (associated(this)) then
       call delete_plan_adv_polar(this%adv)
       call delete_plan_poisson_polar(this%poisson)
       call delete_plan_polar_op(this%grad)
       SLL_DEALLOCATE(this,err)
    end if

  end subroutine delete_SL_polar

!==============
!  advection
!==============


  subroutine advect_CG_polar2(plan,fn,fnp1,phi)

    implicit none

    type(sll_plan_adv_polar), pointer       :: plan
    sll_real64, dimension(:,:), intent(in)  :: fn,phi
    sll_real64, dimension(:,:), intent(out) :: fnp1

    sll_int32 :: nr, ntheta
    sll_real64 :: dt, dr, dtheta, rmin, rmax
    sll_int32 :: i,j,kr,k
    sll_real64 :: r,theta

    nr=plan%nr
    ntheta=plan%ntheta
    dt=plan%dt
    dr=plan%dr
    dtheta=plan%dtheta
    rmin=plan%rmin
    rmax=plan%rmax

    !construction of spline coefficients for f
    call compute_spline_2D(fn,plan%spl_f)

    !if (plan%time_scheme==1) then
       !explicit Euler
       fnp1=fn
       do i=1,nr
          
          do j=1,ntheta
             theta=real(j-1,f64)*dtheta
             r=rmin+real(i-1,f64)*dr
             !theta=theta-dt*plan%field(1,i,j)/r
             !r=r+dt*plan%field(2,i,j)/r

             theta=theta-dt*(phi(i+1,j)-phi(i,j))/(r*dr)
             r=r+dt*(phi(i,j+1)-phi(i,j))/(r*dtheta)


             !call correction_r(r,rmin,rmax)
             !call correction_theta2(theta)
             
             
             
             r=(r-rmin)/(rmax-rmin)
             if(r>=1._f64)then
               r=1._f64
             endif
             if(r<0._f64)then
               r=0._f64
             endif             
             
             r=r*real(nr,f64)
             kr=floor(r)+1
             r=r-real(kr-1,f64)

             theta=theta/(2._f64*sll_pi)
             
             do while(theta>=1._f64)
               theta=theta-1._f64
             enddo
             do while(theta<0._f64)
               theta=theta+1._f64
             enddo
             
             theta=theta*real(ntheta,f64)
             k=floor(theta)+1
             theta=theta-real(k-1,f64)
             
             if(k==ntheta+1)then
               k=ntheta
               theta=1._f64               
             endif
             if(kr==nr+1)then
               kr=nr
               r=1._f64               
             endif
             
             if(k>=ntheta+1)then
               print *,'k=',k,theta
             endif
             if(kr>=nr+1)then
               print *,'kr=',kr,r
             endif


             fnp1(i,j)=(1._f64-theta)*((1.-r)*fn(kr,k)+r*fn(kr+1,k))&
             +theta*((1.-r)*fn(kr,k+1)+r*fn(kr+1,k+1))
             
             r=rmin+(real(kr,f64)-1._f64+r)*dr
             theta=(real(k,f64)-1._f64+theta)*dtheta
             fnp1(i,j)=interpolate_value_2D(r,theta,plan%spl_f)

          end do
       end do
  fnp1(:,ntheta+1)=fnp1(:,1)
  
end subroutine advect_CG_polar2


  !>subroutine advect_CG_polar(plan,in,out)
  !>compute step for Center-Guide equation
  !>plan : sll_plan_adv_polar object
  !>in : distribution function at time t_n, size (nr+1)*(ntheta+1)
  !>out : distribution function at time t_(n+1), size (nr+1)*(ntheta+1)
  subroutine advect_CG_polar(plan,fn,fnp1)

    implicit none

    type(sll_plan_adv_polar),pointer        :: plan
    sll_real64, dimension(:,:), intent(in)  :: fn
    sll_real64, dimension(:,:), intent(out) :: fnp1

    sll_int32 :: nr, ntheta
    sll_real64 :: dt, dr, dtheta, rmin, rmax
    sll_int32 :: i,j,maxiter,iter,kr,k
    sll_real64 :: r,theta,rr,rrn,ttheta,tthetan,tolr,tolth,ar,atheta

    nr=plan%nr
    ntheta=plan%ntheta
    dt=plan%dt
    dr=plan%dr
    dtheta=plan%dtheta
    rmin=plan%rmin
    rmax=plan%rmax

    !construction of spline coefficients for f
    call compute_spline_2D(fn,plan%spl_f)

    if (plan%time_scheme==1) then
       !explicit Euler
       fnp1=fn
       do i=1,nr
          
          do j=1,ntheta
             theta=real(j-1,f64)*dtheta
             r=rmin+real(i-1,f64)*dr
             theta=theta-dt*plan%field(1,i,j)/r
             r=r+dt*plan%field(2,i,j)/r

             !theta=theta-dt*(plan%phi(i+1,j)-plan%phi(i,j))/(r*dr)
             !r=r+dt*(plan%phi(i,j+1)-plan%phi(i,j))/(r*dtheta)


             !call correction_r(r,rmin,rmax)
             !call correction_theta2(theta)
             
             
             
             r=(r-rmin)/(rmax-rmin)
             if(r>=1._f64)then
               r=1._f64
             endif
             if(r<0._f64)then
               r=0._f64
             endif             
             
             r=r*real(nr,f64)
             kr=floor(r)+1
             r=r-real(kr-1,f64)

             theta=theta/(2._f64*sll_pi)
             
             do while(theta>=1._f64)
               theta=theta-1._f64
             enddo
             do while(theta<0._f64)
               theta=theta+1._f64
             enddo
             
             theta=theta*real(ntheta,f64)
             k=floor(theta)+1
             theta=theta-real(k-1,f64)
             
             if(k==ntheta+1)then
               k=ntheta
               theta=1._f64               
             endif
             if(kr==nr+1)then
               kr=nr
               r=1._f64               
             endif
             
             if(k>=ntheta+1)then
               print *,'k=',k,theta
             endif
             if(kr>=nr+1)then
               print *,'kr=',kr,r
             endif


             fnp1(i,j)=(1._f64-theta)*((1.-r)*fn(kr,k)+r*fn(kr+1,k))&
             +theta*((1.-r)*fn(kr,k+1)+r*fn(kr+1,k+1))
             
             r=rmin+(real(kr,f64)-1._f64+r)*dr
             theta=(real(k,f64)-1._f64+theta)*dtheta
             fnp1(i,j)=interpolate_value_2D(r,theta,plan%spl_f)

          end do
       end do

    else if (plan%time_scheme==2) then
       !rotation
       !independant of field

       do i=1,nr+1
          do j=1,ntheta+1
             theta=real(j-1,f64)*dtheta-dt
             call correction_theta(theta)
             fnp1(i,j)=interpolate_value_2D(r,theta,plan%spl_f)
          end do
       end do

    else if (plan%time_scheme==3) then
       !using symplectic Euler with linear interpolation
       !we fix the tolerance and the maximum of iteration
       tolr=dr/5.0_f64
       tolr=1e-14
       maxiter=1000

       do j=1,ntheta
          do i=1,nr+1
             !initialization for r interpolation
             rr=rmin+real(i-1,f64)*dr+dt*plan%field(2,i,j)/(rmin+real(i-1,f64)*dr)
             r=0.0_f64
             iter=0

             call correction_r(rr,rmin,rmax)
             do while (iter<maxiter .and. abs(rrn-rr)>tolr)
                r=(rr-rmin)/(rmax-rmin)
                r=r*real(nr,f64)
                k=floor(r)+1
                r=r-real(k-1,f64)
                rrn=rr
                if (k==nr+1) then
                   !r=0
                   rr=rmin+real(i-1,f64)*dr+dt*plan%field(2,k,j)/rr
                else if (k<nr+1 .and. k>=1) then
                   rr=rmin+real(i-1,f64)*dr+dt*((1.0_f64-r)*plan%field(2,k,j)/rr+r*plan%field(2,k+1,j)/rr)
                else
                   print*,'k is outside of boundaries : error'
                   print*,'exiting the program...'
                   stop
                end if
                call correction_r(rr,rmin,rmax)
                iter=iter+1
             end do
             if (iter==maxiter .and. abs(rrn-rr)>tolr) then
                print*,'not enought iterations for r in symplectic Euler',i,j,rr,rrn
                stop
             end if
             r=(rr-rmin)/(rmax-rmin)
             r=r*real(nr,f64)
             k=floor(r)+1
             r=r-real(k-1,f64)

             if (k/=nr+1) then
                theta=real(j-1,f64)*dtheta-dt*((1.0_f64-r)*plan%field(1,k,j)/rr+r*plan%field(1,k+1,j)/rr)
             else
                theta=real(j-1,f64)*dtheta-dt*plan%field(1,k,j)/rr
             end if
             call correction_theta(theta)

             fnp1(i,j)=interpolate_value_2d(rr,theta,plan%spl_f)
          end do
       end do

    else if (plan%time_scheme==4) then
       !using symplectic Verlet with linear interpolation

       !we fix the tolerance and the maximum of iteration
       tolr=1e-12
       tolth=1e-12
       !tolr=1e-4
       !tolth=1e-4
       maxiter=1000

       do j=1,ntheta
          do i=1,nr+1
             !initialization for r interpolation
             rr=rmin+real(i-1,f64)*dr+dt/2.0_f64*plan%field(2,i,j)/(rmin+real(i-1,f64)*dr)
             rrn=0.0_f64
             r=0.0_f64
             kr=1
             iter=0

             call correction_r(rr,rmin,rmax)
             do while (iter<maxiter .and. abs(rrn-rr)>tolr)
                r=(rr-rmin)/(rmax-rmin)
                r=r*real(nr,f64)
                kr=floor(r)+1
                r=r-real(kr-1,f64)
                rrn=rr
                if (kr==nr+1) then
                   rr=rmin+real(i-1,f64)*dr+0.5_f64*dt*plan%field(2,kr,j)/rr
                else if (kr>0 .and. kr<nr+1) then
                   rr=rmin+real(i-1,f64)*dr+0.5_f64*dt*((1.0_f64-r)*plan%field(2,kr,j)/rr+r*plan%field(2,kr+1,j)/rr)
                else
                   print*,kr
                   print*,'error : kr is not in range'
                   print*,'exiting'
                   stop
                end if
                call correction_r(rr,rmin,rmax)

                iter=iter+1
             end do
             if (iter==maxiter .and. abs(rrn-rr)>tolr) then
                print*,'not enought iterations for r in symplectic Verlet',i,j,kr,rr,rrn
                stop
             end if
             r=(rr-rmin)/(rmax-rmin)
             r=r*real(nr,f64)
             kr=floor(r)+1
             r=r-real(kr-1,f64)

             !initialization for theta interpolation
             ttheta=real(j-1,f64)*dtheta-dt*plan%field(1,i,j)/(rmin+real(i-1,f64)*dr)
             tthetan=3.0_f64*sll_pi
             theta=0.0_f64
             k=1
             iter=0

             call correction_theta(theta)
             do while (iter<maxiter .and. abs(tthetan-ttheta)>tolth .and. &
                  & abs(tthetan+2.0_f64*sll_pi-ttheta)>tolth .and.  abs(tthetan-ttheta-2.0_f64*sll_pi)>tolth)
                theta=ttheta/(2.0_f64*sll_pi)
                theta=theta-real(floor(theta),f64)
                theta=theta*real(ntheta,f64)
                k=floor(theta)+1
                theta=theta-real(k-1,f64)
                if (k==ntheta+1) then
                   k=1
                   theta=0.0_f64
                end if
                tthetan=ttheta
                if (kr==nr+1) then
                   ttheta=real(j-1,f64)*dtheta-0.5_f64*dt*((1.0_f64-theta)*plan%field(1,kr,k)/rr+theta*plan%field(1,kr,k+1)/rr)
                   ttheta=ttheta-0.5_f64*dt*plan%field(1,kr,j)/rr
                else
                   ttheta=real(j-1,f64)*dtheta-0.5_f64*dt*((1.0_f64-theta)*((1.0_f64-r)*plan%field(1,kr,k)/rr+r*plan%field(1,kr+1,k)/rr) &
                        & +theta*((1.0_f64-r)*plan%field(1,kr,k+1)/rr+r*plan%field(1,kr+1,k+1)/rr))
                   ttheta=ttheta-0.5_f64*dt*((1.0_f64-r)*plan%field(1,kr,j)/rr+r*plan%field(1,kr+1,j)/rr)
                end if
                call correction_theta(ttheta)

                iter=iter+1
             end do
             if (iter==maxiter .and. abs(tthetan-ttheta)>tolth .and. abs(tthetan+2.0_f64*sll_pi-ttheta)>tolth &
                  & .and.abs(tthetan-ttheta-2.0_f64*sll_pi)>tolth) then
                print*,'not enought iterations for theta in symplectic Verlet',i,j,k,ttheta,tthetan
                stop
             end if
             theta=ttheta/(2.0_f64*sll_pi)
             theta=theta-real(floor(theta),f64)
             theta=theta*real(ntheta,f64)
             k=floor(theta)+1
             theta=theta-real(k-1,f64)
             if (k==ntheta+1) then
               k=1
               theta=0.0_f64
             end if
             if (kr==nr+1) then
                rr=rr+0.5_f64*dt*((1.0_f64-theta)*plan%field(2,kr,k)/rr+theta*plan%field(2,kr,k+1)/rr)
             else
                rr=rr+0.5_f64*dt*((1.0_f64-theta)*((1.0_f64-r)*plan%field(2,kr,k)/rr+r*plan%field(2,kr+1,k)/rr) &
                     & +theta*((1.0_f64-r)*plan%field(2,kr,k+1)/rr+r*plan%field(2,kr+1,k+1)/rr))
             end if
             call correction_r(rr,rmin,rmax)

             fnp1(i,j)=interpolate_value_2d(rr,ttheta,plan%spl_f)

          end do
       end do

    else if (plan%time_scheme==5) then
       !using fixed point method

       !initialization
       maxiter=10 !10
       tolr=1e-12

       do j=1,ntheta
          do i=1,nr+1
             rr=rmin+real(i-1,f64)*dr
             ttheta=real(j-1,f64)*dtheta
             r=0.0_f64
             kr=i
             ar=0.0_f64
             atheta=0.0_f64
             iter=0

             do while (iter<maxiter .and. abs((rrn-rr)+(tthetan-ttheta))>tolr .and. abs((rrn-rr)+(tthetan+2.0_f64*sll_pi-ttheta))>tolr &
                  & .and. abs((rrn-rr)+(tthetan-ttheta-2.0_f64*sll_pi))>tolr)
                r=(rr-rmin)/(rmax-rmin)
                kr=floor(r)+1
                r=r-real(kr-1,f64)
                theta=ttheta/(2.0_f64*sll_pi)
                theta=theta-real(floor(theta),f64)
                theta=theta*real(ntheta,f64)
                k=floor(theta)+1
                theta=theta-real(k-1,f64)
                if (k==ntheta+1) then
                   k=1
                   theta=0.0_f64
                end if
                if (kr==nr+1) then
                   ar=0.5_f64*dt*((1.0_f64-theta)*((1.0_f64-r)*plan%field(2,kr,k)/rr+theta*((1.0_f64-r)*plan%field(2,kr,k+1)/rr)))
                   atheta=-0.5_f64*dt*((1.0_f64-theta)*((1.0_f64-r)*plan%field(1,kr,k)/rr+theta*((1.0_f64-r)*plan%field(1,kr,k+1)/rr)))
                else
                   ar=-0.5_f64*dt*((1.0_f64-theta)*((1.0_f64-r)*plan%field(2,kr,k)/rr+r*plan%field(2,kr+1,k)/rr) &
                        & +theta*((1.0_f64-r)*plan%field(2,kr,k+1)/rr+r*plan%field(2,kr+1,k+1)/rr))
                   atheta=0.5_f64*dt*((1.0_f64-theta)*((1.0_f64-r)*plan%field(1,kr,k)/rr+r*plan%field(1,kr+1,k)/rr) &
                        & +theta*((1.0_f64-r)*plan%field(1,kr,k+1)/rr+r*plan%field(1,kr+1,k+1)/rr))
                end if
                rrn=rr
                tthetan=ttheta
                rr=rmin+real(i-1,f64)*dr-ar
                ttheta=real(j-1,f64)*dtheta-atheta

                iter=iter+1
             end do
             if (iter==maxiter .and. abs(rrn-rr)+abs(tthetan-ttheta)>tolr) then
                !print*,'no convergence in fixed point method',i,j,&
                !&abs(rrn-rr)+abs(tthetan-ttheta),tolr
             end if

             rr=rmin+real(i-1,f64)*dr-2.0_f64*ar
             ttheta=real(j-1,f64)*dtheta-2.0_f64*atheta
             call correction_r(rr,rmin,rmax)
             call correction_theta(ttheta)
             fnp1(i,j)=interpolate_value_2d(rr,ttheta,plan%spl_f)
          end do
       end do

    else if (plan%time_scheme==6) then
       !using modified symplectic Euler with linear interpolation
       !we fix the tolerance and the maximum of iteration
       tolr=dr/5.0_f64
       tolr=1e-14
       maxiter=1000

       do j=1,ntheta
          do i=1,nr+1
             !initialization for r interpolation
             rr=rmin+real(i-1,f64)*dr-dt*plan%field(1,i,j)
             r=0.0_f64
             iter=0

             call correction_r(rr,rmin,rmax)
             do while (iter<maxiter .and. abs(rrn-rr)>tolr)
                r=(rr-rmin)/(rmax-rmin)
                r=r*real(nr,f64)
                k=floor(r)+1
                r=r-real(k-1,f64)
                rrn=rr
                if (k==nr+1) then
                   !r=0
                   rr=rmin+real(i-1,f64)*dr-dt*plan%field(2,k,j)/rr
                else if (k<nr+1 .and. k>=1) then
                   rr=rmin+real(i-1,f64)*dr-dt*((1.0_f64-r)*plan%field(2,k,j)/rr+r*plan%field(1,k+1,j)/rr)
                else
                   print*,'k is outside of boundaries : error'
                   print*,'exiting the program...'
                   stop
                end if
                call correction_r(rr,rmin,rmax)
                iter=iter+1
             end do
             if (iter==maxiter .and. abs(rrn-rr)>tolr) then
                print*,'not enought iterations for r in modified symplectic Euler',i,j,rr,rrn
                stop
             end if
             r=(rr-rmin)/(rmax-rmin)
             r=r*real(nr,f64)
             k=floor(r)+1
             r=r-real(k-1,f64)

             if (k/=nr+1) then
                theta=real(j-1,f64)*dtheta*rr/(rmin+real(i-1,f64)*dr)+dt &
                     & *((1.0_f64-r)*(plan%field(1,k,j)/(rmin+real(k-1,f64)*dr)+real(j-1,f64)*dtheta/rr**2*plan%field(2,k,j)) &
                     & +r*(plan%field(1,k+1,j)/(rmin+real(k+1-1,f64)*dr)+real(j-1,f64)*dtheta/rr**2*plan%field(2,k+1,j)))
             else
                theta=real(j-1,f64)*dtheta*rr/(rmin+real(i-1,f64)*dr)+dt &
                     & *(plan%field(1,k,j)/(rmin+real(k-1,f64)*dr)+real(j-1,f64)*dtheta/(rr**2)*plan%field(2,k,j))
             end if
             call correction_theta(theta)

             fnp1(i,j)=interpolate_value_2d(rr,theta,plan%spl_f)
          end do
       end do

    else if (plan%time_scheme==7) then
       !using modified symplectic Verlet with linear interpolation

       !we fix the tolerance and the maximum of iteration
       tolr=1e-12
       tolth=1e-12
       tolr=1e-4
       tolth=1e-4
       maxiter=1000

       do j=1,ntheta
          do i=1,nr+1
             !initialization for r interpolation
             rr=rmin+real(i-1,f64)*dr+dt/2.0_f64*plan%field(2,i,j)/(rmin+real(i-1,f64)*dr)
             rrn=0.0_f64
             r=0.0_f64
             kr=1
             iter=0

             call correction_r(rr,rmin,rmax)
             do while (iter<maxiter .and. abs(rrn-rr)>tolr)
                r=(rr-rmin)/(rmax-rmin)
                r=r*real(nr,f64)
                kr=floor(r)+1
                r=r-real(kr-1,f64)
                rrn=rr
                if (kr==nr+1) then
                   rr=rmin+real(i-1,f64)*dr-0.5_f64*dt*plan%field(2,kr,j)/rr
                else if (kr>0 .and. kr<nr+1) then
                   rr=rmin+real(i-1,f64)*dr-0.5_f64*dt*((1.0_f64-r)*plan%field(2,kr,j)/rr+r*plan%field(2,kr+1,j)/rr)
                else
                   print*,kr
                   print*,'error : kr is not in range'
                   print*,'exiting'
                   stop
                end if
                call correction_r(rr,rmin,rmax)

                iter=iter+1
             end do
             if (iter==maxiter .and. abs(rrn-rr)>tolr) then
                print*,'not enought iterations for r in symplectic Verlet',i,j,kr,rr,rrn
                stop
             end if
             r=(rr-rmin)/(rmax-rmin)
             r=r*real(nr,f64)
             kr=floor(r)+1
             r=r-real(kr-1,f64)

             !initialization for theta interpolation
             ttheta=real(j-1,f64)*dtheta-dt*plan%field(2,i,j)
             tthetan=3.0_f64*sll_pi
             theta=0.0_f64
             k=1
             iter=0

             call correction_theta(theta)
             do while (iter<maxiter .and. abs(tthetan-ttheta)>tolth .and. &
                  & abs(tthetan+2.0_f64*sll_pi-ttheta)>tolth .and.  abs(tthetan-ttheta-2.0_f64*sll_pi)>tolth)
                theta=ttheta/(2.0_f64*sll_pi)
                theta=theta-real(floor(theta),f64)
                theta=theta*real(ntheta,f64)
                k=floor(theta)+1
                theta=theta-real(k-1,f64)
                if (k==ntheta+1) then
                   k=1
                   theta=0.0_f64
                end if
                tthetan=ttheta
                if (kr==nr+1) then

                   ttheta=real(j-1,f64)*dtheta*rr/(rmin+real(i-1,f64)*dr)+0.5_f64*dt &
                        & *((1.0_f64-theta)*(plan%field(1,kr,k)/(rmin+real(kr-1,f64)*dr)+real(k-1,f64)*dtheta/(rr**2)*plan%field(2,kr,k)) &
                        & +theta*(plan%field(1,kr,k+1)/(rmin+real(kr-1,f64)*dr)+real(k+1-1,f64)*dtheta/(rr**2)*plan%field(2,kr,k+1)))
                    ttheta=ttheta+0.5_f64*dt*(plan%field(1,kr,j)/(rmin+real(kr-1,f64)*dr)+real(j-1,f64)*dtheta/(rr**2)*plan%field(2,kr,j))

                else

                   ttheta=real(j-1,f64)*dtheta*rr/(rmin+real(i-1,f64)*dr)+0.5_f64*dt &
                        & *((1.0_f64-theta)*((1.0_f64-r)*(plan%field(1,kr,k)/(rmin+real(kr-1,f64)*dr)+real(k-1,f64)*dtheta/(rr**2)*plan%field(2,kr,k)) &
                        & +r*(plan%field(1,kr+1,k)/(rmin+real(kr+1-1,f64)*dr)+real(k-1,f64)*dtheta/rr**2*plan%field(2,kr+1,k))) &
                        & +theta*((1.0_f64-r)*(plan%field(1,kr,k+1)/(rmin+real(i-1,f64)*dr)+real(j+1-1,f64)*dtheta/(rr**2)*plan%field(2,kr,k+1)) &
                        & +r*(plan%field(1,kr+1,k+1)/(rmin+real(kr+1-1,f64)*dr)+real(k+1-1,f64)*dtheta/rr**2*plan%field(2,kr+1,k+1))))
                   ttheta=ttheta+0.5_f64*dt*((1.0_f64-r)*(plan%field(1,kr,j)/(rmin+real(k-1,f64)*dr)+real(j-1,f64)*dtheta/(rr**2)*plan%field(2,kr,j)) &
                        & +r*(plan%field(1,kr+1,j)/(rmin+real(kr+1-1,f64)*dr)+real(j-1,f64)*dtheta/rr**2*plan%field(2,kr+1,j)))

                end if
                call correction_theta(ttheta)

                iter=iter+1
             end do
             if (iter==maxiter .and. abs(tthetan-ttheta)>tolth .and. abs(tthetan+2.0_f64*sll_pi-ttheta)>tolth &
                  & .and.abs(tthetan-ttheta-2.0_f64*sll_pi)>tolth) then
                print*,'not enought iterations for theta in symplectic Verlet',i,j,k,ttheta,tthetan
                stop
             end if
             theta=ttheta/(2.0_f64*sll_pi)
             theta=theta-real(floor(theta),f64)
             theta=theta*real(ntheta,f64)
             k=floor(theta)+1
             theta=theta-real(k-1,f64)
             if (k==ntheta+1) then
               k=1
               theta=0.0_f64
             end if
             if (kr==nr+1) then
                rr=rr-0.5_f64*dt*((1.0_f64-theta)*plan%field(2,kr,k)/rr+theta*plan%field(2,kr,k+1)/rr)
             else
                rr=rr+0.5_f64*dt*((1.0_f64-theta)*((1.0_f64-r)*plan%field(2,kr,k)/rr+r*plan%field(2,kr+1,k)/rr) &
                     & +theta*((1.0_f64-r)*plan%field(2,kr,k+1)/rr+r*plan%field(2,kr+1,k+1)/rr))
             end if
             call correction_r(rr,rmin,rmax)

             fnp1(i,j)=interpolate_value_2d(rr,ttheta,plan%spl_f)

          end do
       end do

    else if (plan%time_scheme==8) then
       !using modified fixed point method

       !initialization
       maxiter=1000
       tolr=1e-10

       do j=1,ntheta
          do i=1,nr+1
             rr=rmin+real(i-1,f64)*dr
             ttheta=real(j-1,f64)*dtheta
             r=0.0_f64
             kr=i
             ar=0.0_f64
             atheta=0.0_f64
             iter=0

             do while (iter<maxiter .and. abs((rrn-rr)+(tthetan-ttheta))>tolr .and. abs((rrn-rr)+(tthetan+2.0_f64*sll_pi-ttheta))>tolr &
                  & .and. abs((rrn-rr)+(tthetan-ttheta-2.0_f64*sll_pi))>tolr)
                r=(rr-rmin)/(rmax-rmin)
                kr=floor(r)+1
                r=r-real(kr-1,f64)
                theta=ttheta/(2.0_f64*sll_pi)
                theta=theta-real(floor(theta),f64)
                theta=theta*real(ntheta,f64)
                k=floor(theta)+1
                theta=theta-real(k-1,f64)
                if (k==ntheta+1) then
                   k=1
                   theta=0.0_f64
                end if
                if (kr==nr+1) then

                   ar=-0.5_f64*dt*((1.0_f64-theta)*plan%field(2,kr,k)/rr+theta*plan%field(2,kr,k+1)/rr)
                   atheta=0.5_f64*dt*((1.0_f64-theta)*(plan%field(1,kr,k)/(rmin+real(kr-1,f64)*dr)+real(k-1,f64)*dtheta/rr*2*plan%field(2,kr,k)) &
                        & +theta*(plan%field(1,kr,k+1)/(rmin+real(kr-1,f64)*dr)+real(k+1-1,f64)*dtheta/rr*2*plan%field(2,kr,k+1)))

                else

                   ar=-0.5_f64*dt*((1.0_f64-theta)*((1.0_f64-r)*plan%field(2,kr,k)/rr+r*plan%field(1,kr+1,k)/rr) &
                        & +theta*((1.0_f64-r)*plan%field(1,kr,k+1)/rr+r*plan%field(1,kr+1,k+1)/rr))
                   atheta=0.5_f64*dt*((1.0_f64-theta)*((1.0_f64-r)*(plan%field(1,kr,k)/(rmin+real(kr-1,f64)*dr)+real(k-1,f64)*dtheta/rr**2*plan%field(2,kr,k)) &
                        & +r*(plan%field(1,kr+1,k)/(rmin+real(kr+1-1,f64)*dr)+real(k-1,f64)*dtheta/rr**2*plan%field(2,kr+1,k))) &
                        & +theta*((1.0_f64-r)*(plan%field(1,kr,k+1)/(rmin+real(kr-1,f64)*dr)+real(k+1-1,f64)*dtheta/rr**2*plan%field(2,kr,k+1)) &
                        & +r*(plan%field(1,kr+1,k+1)/(rmin+real(kr+1-1,f64)*dr)+real(k+1-1,f64)*dtheta/rr**2*plan%field(2,kr+1,k+1))))

                end if
                rrn=rr
                tthetan=ttheta
                rr=rmin+real(i-1,f64)*dr-ar
                ttheta=real(j-1,f64)*dtheta-atheta

                iter=iter+1
             end do
             if (iter==maxiter .and. (rrn-rr)+(tthetan-ttheta)>tolr) then
                print*,'no convergence in fixe point methode',i,j
             end if

             rr=rmin+real(i-1,f64)*dr-2.0_f64*ar
             ttheta=real(j-1,f64)*dtheta-2.0_f64*atheta
             call correction_r(rr,rmin,rmax)
             call correction_theta(ttheta)
             fnp1(i,j)=interpolate_value_2d(rr,ttheta,plan%spl_f)
          end do
       end do

    end if

    fnp1(:,ntheta+1)=fnp1(:,1)

  end subroutine advect_CG_polar

  subroutine compute_remap(plan,fn,fnp1,interp_case,PPM_order)

    implicit none

    type(sll_plan_adv_polar),pointer        :: plan
    sll_real64, dimension(:,:), intent(in)  :: fn
    sll_real64, dimension(:,:), intent(out) :: fnp1
    sll_int32,intent(in) :: interp_case,PPM_order
    !sll_real64, intent(in)  :: dt

    sll_int32 :: nr,ntheta
    sll_real64 :: dr, dtheta, rmin, rmax
    sll_int32 :: i,j,maxiter,iter,kr,k
    sll_real64 :: r,theta,rr,rrn,ttheta,tthetan,tolr,tolth,ar,atheta,geom_x(2,2),dt
    sll_real64,dimension(:,:,:),allocatable::carac
    sll_real64,dimension(:,:),allocatable::buf2d

    

    nr=plan%nr
    ntheta=plan%ntheta
    dt=plan%dt
    dr=plan%dr
    dtheta=plan%dtheta
    rmin=plan%rmin
    rmax=plan%rmax
    geom_x(1,1)=rmin
    geom_x(2,1)=rmax-rmin
    geom_x(1,2)=0._f64
    geom_x(2,2)=2._f64*sll_pi
    !interp_case=3
    !ppm_order=1

    allocate(carac(2,-1:nr+1,-1:ntheta))
    allocate(buf2d(0:nr,0:ntheta-1))
    call compute_carac(plan,nr,ntheta,dt,dr,dtheta,rmin,rmax,carac)
!do i=-1,nr
!do j=-1,ntheta
!print*,i,j,carac(1,i,j), carac(2,i,j)
!enddo
!enddo
    do i=1,nr+1
      fnp1(i,1:ntheta) = (rmin+real(i-1,f64)*(rmax-rmin)/real(nr,f64))*fn(i,1:ntheta)
    enddo
    call advect2d_CSL_CG(geom_x,fnp1(1:nr+1,1:ntheta),nr,ntheta,buf2d,interp_case,carac,PPM_order)
    !fnp1=fn

    do i=1,nr+1
      fnp1(i,1:ntheta) = fnp1(i,1:ntheta)/(rmin+real(i-1,f64)*(rmax-rmin)/real(nr,f64))
    enddo


!print*,fn
!stop

    deallocate(carac,buf2d)

  end subroutine compute_remap

  subroutine compute_carac(plan,nr,ntheta,dt,dr,dtheta,rmin,rmax,carac)

    implicit none

    type(sll_plan_adv_polar),pointer :: plan
    sll_int32,intent(in) :: nr, ntheta
    sll_real64,intent(in) :: dt, dr, dtheta, rmin, rmax
    sll_real64 :: r,theta,rr,rrn,ttheta,tthetan,tolr,tolth,ar,atheta
    sll_int32 :: i,j,maxiter,iter,kr,k
    sll_real64,dimension(2,-1:nr+1,-1:ntheta),intent(inout)::carac  


    !nr=plan%nr
    !ntheta=plan%ntheta
    !dt=plan%dt
    !dr=plan%dr
    !dtheta=plan%dtheta
    !rmin=plan%rmin

    if (plan%time_scheme==1) then
       !explicit Euler
       do i=1,nr+1
          do j=1,ntheta+1
             theta=real(j-1,f64)*dtheta
             r=rmin+real(i-1,f64)*dr
             theta=theta-dt*plan%field(1,i,j)/r
             r=r+dt*plan%field(2,i,j)/r
             carac(1,i-1,j-1)=r
             carac(2,i-1,j-1)=theta
!print*,plan%field(1,i,j)/r
!print*,plan%field(2,i,j)/r
          enddo
       enddo

    else if (plan%time_scheme==4) then
       !symplectic Verlet with linear interpolation

       !we fix the tolerance and the maximum of iteration
       tolr=1e-12
       tolth=1e-12
       !tolr=1e-4
       !tolth=1e-4
       maxiter=10!00

       do j=1,ntheta+1
          do i=1,nr+1
             !initialization for r interpolation
             rr=rmin+real(i-1,f64)*dr+dt/2.0_f64*plan%field(2,i,j)/(rmin+real(i-1,f64)*dr)
             rrn=0.0_f64
             r=0.0_f64
             kr=1
             iter=0

             !call correction_r(rr,rmin,rmax)
             do while (iter<maxiter .and. abs(rrn-rr)>tolr)
                r=(rr-rmin)/(rmax-rmin)
                r=r*real(nr,f64)
                kr=floor(r)+1
                r=r-real(kr-1,f64)
                rrn=rr
                if (kr>nr) then
                   rr=rmin+real(i-1,f64)*dr+0.5_f64*dt*plan%field(2,nr+1,j)/rr
                else if (kr>=1 .and. kr<=nr) then
                   rr=rmin+real(i-1,f64)*dr+0.5_f64*dt*((1.0_f64-r)*plan%field(2,kr,j)/rr+r*plan%field(2,kr+1,j)/rr)
                else
                   rr=rmin+real(i-1,f64)*dr+0.5_f64*dt*plan%field(2,1,j)/rr
                   !print*,kr
                   !print*,'error : kr is not in range'
                   !print*,'exiting'
                   !stop
                end if
                !call correction_r(rr,rmin,rmax)

                iter=iter+1
             end do
             !if (iter==maxiter .and. abs(rrn-rr)>tolr) then
             !   print*,'not enought iterations for r in symplectic Verlet',i,j,kr,rr,rrn
             !   stop
             !end if
             r=(rr-rmin)/(rmax-rmin)
             r=r*real(nr,f64)
             kr=floor(r)+1
             r=r-real(kr-1,f64)

             !initialization for theta interpolation
             ttheta=real(j-1,f64)*dtheta-dt*plan%field(1,i,j)/(rmin+real(i-1,f64)*dr)
             tthetan=3.0_f64*sll_pi
             theta=0.0_f64
             k=1
             iter=0

             !call correction_theta(theta)
             do while (iter<maxiter .and. abs(tthetan-ttheta)>tolth .and. &
                  & abs(tthetan+2.0_f64*sll_pi-ttheta)>tolth .and.  abs(tthetan-ttheta-2.0_f64*sll_pi)>tolth)
                theta=ttheta/(2.0_f64*sll_pi)
                theta=theta-real(floor(theta),f64)
                theta=theta*real(ntheta,f64)
                k=floor(theta)+1
                theta=theta-real(k-1,f64)
                if (k>ntheta) then
                   k=modulo(k-1,ntheta)+1
                !   theta=0.0_f64
                end if
                tthetan=ttheta
                if (kr>nr) then
                   ttheta=real(j-1,f64)*dtheta-0.5_f64*dt*((1.0_f64-theta)*plan%field(1,nr+1,k)/rr+theta*plan%field(1,nr+1,k+1)/rr)
                   ttheta=ttheta-0.5_f64*dt*plan%field(1,nr+1,j)/rr
                else if (kr>=1 .and. kr<=nr) then
                   ttheta=real(j-1,f64)*dtheta-0.5_f64*dt*((1.0_f64-theta)*((1.0_f64-r)*plan%field(1,kr,k)/rr+r*plan%field(1,kr+1,k)/rr) &
                        & +theta*((1.0_f64-r)*plan%field(1,kr,k+1)/rr+r*plan%field(1,kr+1,k+1)/rr))
                   ttheta=ttheta-0.5_f64*dt*((1.0_f64-r)*plan%field(1,kr,j)/rr+r*plan%field(1,kr+1,j)/rr)
                else
                   ttheta=real(j-1,f64)*dtheta-0.5_f64*dt*((1.0_f64-theta)*plan%field(1,1,k)/rr+theta*plan%field(1,1,k+1)/rr)
                   ttheta=ttheta-0.5_f64*dt*plan%field(1,1,j)/rr
                end if
                !call correction_theta(ttheta)

                iter=iter+1
             end do
             !if (iter==maxiter .and. abs(tthetan-ttheta)>tolth .and. abs(tthetan+2.0_f64*sll_pi-ttheta)>tolth &
             !     & .and.abs(tthetan-ttheta-2.0_f64*sll_pi)>tolth) then
             !   print*,'not enought iterations for theta in symplectic Verlet',i,j,k,ttheta,tthetan
             !   stop
             !end if
             theta=ttheta/(2.0_f64*sll_pi)
             theta=theta-real(floor(theta),f64)
             theta=theta*real(ntheta,f64)
             k=floor(theta)+1
             theta=theta-real(k-1,f64)
             if (k>ntheta) then
                k=modulo(k-1,ntheta)+1
             !   theta=0.0_f64
             end if
             if (kr>nr) then
                rr=rr+0.5_f64*dt*((1.0_f64-theta)*plan%field(2,nr+1,k)/rr+theta*plan%field(2,nr+1,k+1)/rr)
             else if (kr>=1 .and. kr<=nr) then
                rr=rr+0.5_f64*dt*((1.0_f64-theta)*((1.0_f64-r)*plan%field(2,kr,k)/rr+r*plan%field(2,kr+1,k)/rr) &
                     & +theta*((1.0_f64-r)*plan%field(2,kr,k+1)/rr+r*plan%field(2,kr+1,k+1)/rr))
             else
                rr=rr+0.5_f64*dt*((1.0_f64-theta)*plan%field(2,1,k)/rr+theta*plan%field(2,1,k+1)/rr)
             end if
             !call correction_r(rr,rmin,rmax)
!if ((i==1) .and. (j==1)) then
!print*,rr,ttheta
!endif
             carac(1,i-1,j-1)=rr
             carac(2,i-1,j-1)=ttheta

          end do
       end do

    endif


       do i=0,nr
          !carac(1,i,ntheta) = carac(1,i,0)
          !carac(2,i,ntheta) = carac(2,i,0)+ntheta*dtheta
          carac(1,i,-1) = carac(1,i,ntheta-1)
          carac(2,i,-1) = carac(2,i,ntheta-1)-ntheta*dtheta
       enddo
       do j=-1,ntheta
          !more suitable for periodic conditions, but here we are Dirichlet for r
          !carac(1,nr,j) = carac(1,0,j)+nr*dr
          !carac(2,nr,j) = carac(2,0,j)          
          !carac(1,-1,j) = carac(1,nr-1,j)-nr*dr
          !carac(2,-1,j) = carac(2,nr-1,j)      
          
          carac(1,nr+1,j) = carac(1,nr,j)+dr
          carac(2,nr+1,j) = carac(2,nr,j)          
          carac(1,-1,j) = carac(1,0,j)-dr
          carac(2,-1,j) = carac(2,0,j)      
       enddo
!print*,carac(2,0,:)
!stop

  end subroutine compute_carac

  subroutine compute_plan_carac(plan,nr,ntheta,dt,dr,dtheta,rmin,rmax)

    implicit none

    type(sll_plan_adv_polar),pointer :: plan
    sll_int32,intent(in) :: nr, ntheta
    sll_real64,intent(in) :: dt, dr, dtheta, rmin, rmax
    sll_real64 :: r,theta,rr,rrn,ttheta,tthetan,tolr,tolth,ar,atheta
    sll_int32 :: i,j,maxiter,iter,kr,k
    !sll_real64,dimension(2,-1:nr,-1:ntheta),intent(inout)::carac  


    !nr=plan%nr
    !ntheta=plan%ntheta
    !dt=plan%dt
    !dr=plan%dr
    !dtheta=plan%dtheta
    !rmin=plan%rmin

    if (plan%time_scheme==1) then
       !explicit Euler
       do i=1,nr+1
          do j=1,ntheta+1
             theta=real(j-1,f64)*dtheta
             r=rmin+real(i-1,f64)*dr
             theta=theta-dt*plan%field(1,i,j)/r
             r=r+dt*plan%field(2,i,j)/r
             plan%carac(1,i,j)=r
             plan%carac(2,i,j)=theta
!print*,plan%field(1,i,j)/r
!print*,plan%field(2,i,j)/r
          enddo
       enddo

    else if (plan%time_scheme==4) then
       !symplectic Verlet with linear interpolation

       !we fix the tolerance and the maximum of iteration
       tolr=1e-12
       tolth=1e-12
       !tolr=1e-4
       !tolth=1e-4
       maxiter=10!00

       do j=1,ntheta+1
          do i=1,nr+1
             !initialization for r interpolation
             rr=rmin+real(i-1,f64)*dr+dt/2.0_f64*plan%field(2,i,j)/(rmin+real(i-1,f64)*dr)
             rrn=0.0_f64
             r=0.0_f64
             kr=1
             iter=0

             !call correction_r(rr,rmin,rmax)
             do while (iter<maxiter .and. abs(rrn-rr)>tolr)
                r=(rr-rmin)/(rmax-rmin)
                r=r*real(nr,f64)
                kr=floor(r)+1
                r=r-real(kr-1,f64)
                rrn=rr
                if (kr>nr) then
                   rr=rmin+real(i-1,f64)*dr+0.5_f64*dt*plan%field(2,nr+1,j)/rr
                else if (kr>=1 .and. kr<=nr) then
                   rr=rmin+real(i-1,f64)*dr+0.5_f64*dt*((1.0_f64-r)*plan%field(2,kr,j)/rr+r*plan%field(2,kr+1,j)/rr)
                else
                   rr=rmin+real(i-1,f64)*dr+0.5_f64*dt*plan%field(2,1,j)/rr
                   !print*,kr
                   !print*,'error : kr is not in range'
                   !print*,'exiting'
                   !stop
                end if
                !call correction_r(rr,rmin,rmax)

                iter=iter+1
             end do
             !if (iter==maxiter .and. abs(rrn-rr)>tolr) then
             !   print*,'not enought iterations for r in symplectic Verlet',i,j,kr,rr,rrn
             !   stop
             !end if
             r=(rr-rmin)/(rmax-rmin)
             r=r*real(nr,f64)
             kr=floor(r)+1
             r=r-real(kr-1,f64)

             !initialization for theta interpolation
             ttheta=real(j-1,f64)*dtheta-dt*plan%field(1,i,j)/(rmin+real(i-1,f64)*dr)
             tthetan=3.0_f64*sll_pi
             theta=0.0_f64
             k=1
             iter=0

             !call correction_theta(theta)
             do while (iter<maxiter .and. abs(tthetan-ttheta)>tolth .and. &
                  & abs(tthetan+2.0_f64*sll_pi-ttheta)>tolth .and.  abs(tthetan-ttheta-2.0_f64*sll_pi)>tolth)
                theta=ttheta/(2.0_f64*sll_pi)
                theta=theta-real(floor(theta),f64)
                theta=theta*real(ntheta,f64)
                k=floor(theta)+1
                theta=theta-real(k-1,f64)
                if (k>ntheta) then
                   k=modulo(k-1,ntheta)+1
                !   theta=0.0_f64
                end if
                tthetan=ttheta
                if (kr>nr) then
                   ttheta=real(j-1,f64)*dtheta-0.5_f64*dt*((1.0_f64-theta)*plan%field(1,nr+1,k)/rr+theta*plan%field(1,nr+1,k+1)/rr)
                   ttheta=ttheta-0.5_f64*dt*plan%field(1,nr+1,j)/rr
                else if (kr>=1 .and. kr<=nr) then
                   ttheta=real(j-1,f64)*dtheta-0.5_f64*dt*((1.0_f64-theta)*((1.0_f64-r)*plan%field(1,kr,k)/rr+r*plan%field(1,kr+1,k)/rr) &
                        & +theta*((1.0_f64-r)*plan%field(1,kr,k+1)/rr+r*plan%field(1,kr+1,k+1)/rr))
                   ttheta=ttheta-0.5_f64*dt*((1.0_f64-r)*plan%field(1,kr,j)/rr+r*plan%field(1,kr+1,j)/rr)
                else
                   ttheta=real(j-1,f64)*dtheta-0.5_f64*dt*((1.0_f64-theta)*plan%field(1,1,k)/rr+theta*plan%field(1,1,k+1)/rr)
                   ttheta=ttheta-0.5_f64*dt*plan%field(1,1,j)/rr
                end if
                !call correction_theta(ttheta)

                iter=iter+1
             end do
             !if (iter==maxiter .and. abs(tthetan-ttheta)>tolth .and. abs(tthetan+2.0_f64*sll_pi-ttheta)>tolth &
             !     & .and.abs(tthetan-ttheta-2.0_f64*sll_pi)>tolth) then
             !   print*,'not enought iterations for theta in symplectic Verlet',i,j,k,ttheta,tthetan
             !   stop
             !end if
             theta=ttheta/(2.0_f64*sll_pi)
             theta=theta-real(floor(theta),f64)
             theta=theta*real(ntheta,f64)
             k=floor(theta)+1
             theta=theta-real(k-1,f64)
             if (k>ntheta) then
                k=modulo(k-1,ntheta)+1
             !   theta=0.0_f64
             end if
             if (kr>nr) then
                rr=rr+0.5_f64*dt*((1.0_f64-theta)*plan%field(2,nr+1,k)/rr+theta*plan%field(2,nr+1,k+1)/rr)
             else if (kr>=1 .and. kr<=nr) then
                rr=rr+0.5_f64*dt*((1.0_f64-theta)*((1.0_f64-r)*plan%field(2,kr,k)/rr+r*plan%field(2,kr+1,k)/rr) &
                     & +theta*((1.0_f64-r)*plan%field(2,kr,k+1)/rr+r*plan%field(2,kr+1,k+1)/rr))
             else
                rr=rr+0.5_f64*dt*((1.0_f64-theta)*plan%field(2,1,k)/rr+theta*plan%field(2,1,k+1)/rr)
             end if
             !call correction_r(rr,rmin,rmax)
!if ((i==1) .and. (j==1)) then
!print*,rr,ttheta
!endif
             plan%carac(1,i,j)=rr
             plan%carac(2,i,j)=ttheta

          end do
       end do

    endif



  end subroutine compute_plan_carac




  !> subroutine correction_r(r,rmin,rmax)
  !> correction of r to stay in the domain
  !> r : value tu correct
  !> rmin and rmax : boundaries of the domain
  subroutine correction_r(r,rmin,rmax)

    implicit none

    sll_real64, intent(in) :: rmin,rmax
    sll_real64, intent(inout) :: r

    if (r>rmax) then
       r=rmax
    else if (r<rmin) then
       r=rmin
    end if

  end subroutine correction_r


  !> subroutine correction_theta(theta)
  !> correction of theta to stay in [0;2pi]
  !> theta : angle to correct
  subroutine correction_theta2(theta)

    implicit none

    sll_real64, intent(inout) :: theta
    sll_real64 :: th

    th=theta
    if (theta>2.0_f64*sll_pi) then
       !theta=theta-real(floor(theta/(2.0_f64*sll_pi)),f64)*2.0_f64*sll_pi
       theta=modulo(theta,2.0_f64*sll_pi)
    end if
    if (theta>2.0_f64*sll_pi) then
       print*,'je ne sais pas calculer!'
       print*,th
       print*,theta
    end if
    if (theta<0.0_f64) then
       theta=-theta
       theta=modulo(theta,2.0_f64*sll_pi)
       theta=-theta+2.0_f64*sll_pi
    end if
    if (theta<0.0_f64) then
       print*,'je ne sais pas calculer!'
       print*,th
       print*,theta
    end if

  end subroutine correction_theta2

  subroutine correction_theta(theta)

    implicit none

    sll_real64, intent(inout) :: theta
    sll_real64 :: th

    th=theta
    theta=theta-real(floor(theta/(2.0_f64*sll_pi)),f64)*2.0_f64*sll_pi
    if (theta>2.0_f64*sll_pi) then
       theta=theta-2.0_f64*sll_pi
    end if
    if (theta<0.0_f64) then
       theta=theta+2.0_f64*sll_pi
    end if

  end subroutine correction_theta

  !>subroutine SL_remap(plan,in,out,interp_case,PPM_order)
  !>plan : sll_SL_polar object, contains plan for Poisson, gradient and advection
  !>in : distribution function at time n, size (nr+1)*(ntheta+1)
  !>out : distribution function at time n+1, size (nr+1)*(ntheta+1)
  subroutine SL_remap(plan,in,out,interp_case,PPM_order)
 
    implicit none

    type(sll_SL_polar), pointer               :: plan
    sll_real64, dimension(:,:), intent(inout) :: in
    sll_real64, dimension(:,:), intent(out)   :: out
    sll_int32,intent(in) :: interp_case,PPM_order
    sll_real64 :: dt

    dt=plan%adv%dt
    call poisson_solve_polar(plan%poisson,in,plan%phi)
    call compute_grad_field(plan%grad,plan%phi,plan%adv%field)    
    call compute_remap(plan%adv,in,out,interp_case,PPM_order)

    !print *,sum(abs(plan%adv%field(1,1,:)))

  end subroutine SL_remap


  !>subroutine SL_remap(plan,in,out,interp_case,PPM_order)
  !>plan : sll_SL_polar object, contains plan for Poisson, gradient and advection
  !>in : distribution function at time n, size (nr+1)*(ntheta+1)
  !>out : distribution function at time n+1, size (nr+1)*(ntheta+1)
  subroutine SL_remap_ordre_2(plan,in,out,interp_case,PPM_order)
 
    implicit none

    type(sll_SL_polar), pointer               :: plan
    sll_real64, dimension(:,:), intent(inout) :: in
    sll_real64, dimension(:,:), intent(out)   :: out
    sll_real64, dimension(:,:), allocatable   :: aux
    sll_int32,intent(in) :: interp_case,PPM_order
    sll_real64 :: dt
    sll_int32 :: nr,ntheta

    nr=plan%adv%nr
    ntheta=plan%adv%ntheta
    allocate(aux(1:nr+1,1:ntheta+1))
    dt=plan%adv%dt/2.0_f64
    plan%adv%dt=dt
    aux=in

    call poisson_solve_polar(plan%poisson,in,plan%phi)
    call compute_grad_field(plan%grad,plan%phi,plan%adv%field)    
    call compute_remap(plan%adv,aux,out,interp_case,PPM_order)

    call poisson_solve_polar(plan%poisson,out,plan%phi)
    call compute_grad_field(plan%grad,plan%phi,plan%adv%field)

    dt=2._f64*dt
    plan%adv%dt=dt
    call compute_remap(plan%adv,in,out,interp_case,PPM_order)

    deallocate(aux)
    !print *,sum(abs(plan%adv%field(1,1,:)))

  end subroutine SL_remap_ordre_2


  !>subroutine SL_classic(plan,in,out)
  !>computes the classic semi-Lagrangian scheme for Vlasov-Poisson equation
  !>plan : sll_SL_polar object, contains plan for Poisso, gradient and advection
  !>in : distribution function at time n, size (nr+1)*(ntheta+1)
  !>out : distribution function at time n+1, size (nr+1)*(ntheta+1)
  subroutine SL_classic(plan,in,out)

    implicit none

    type(sll_SL_polar), pointer               :: plan
    sll_real64, dimension(:,:), intent(inout) :: in
    sll_real64, dimension(:,:), intent(out)   :: out

    call poisson_solve_polar(plan%poisson,in,plan%phi)
    call compute_grad_field(plan%grad,plan%phi,plan%adv%field)
    
        
    call advect_CG_polar(plan%adv,in,out)

    !print *,sum(abs(plan%adv%field(1,1,:)))

  end subroutine SL_classic


  !>subroutine SL_ordre_2(plan,in,out)
  !>computes the semi-Lagrangian scheme order 2
  !>plan : sll_SL_polar object, contains plan for Poisso, gradient and advection
  !>in : distribution function at time n, size (nr+1)*(ntheta+1)
  !>out : distribution function at time n+1, size (nr+1)*(ntheta+1)
  subroutine SL_ordre_2(plan,in,out)

    implicit none

    type(sll_SL_polar), pointer               :: plan
    sll_real64, dimension(:,:), intent(inout) :: in
    sll_real64, dimension(:,:), intent(out)   :: out

    sll_real64 :: dt

    dt          = plan%adv%dt
    plan%adv%dt = dt/2.0_f64

    call poisson_solve_polar(plan%poisson,in,plan%phi)
    call compute_grad_field(plan%grad,plan%phi,plan%adv%field)
    call advect_CG_polar(plan%adv,in,out)
    !we just obtained f^(n+1/2)
    call poisson_solve_polar(plan%poisson,out,plan%phi)
    call compute_grad_field(plan%grad,plan%phi,plan%adv%field)
    !we just obtained E^(n+1/2)
    plan%adv%dt = dt
    call advect_CG_polar(plan%adv,in,out)

  end subroutine SL_ordre_2







  subroutine print2d(dom,ftab,Nx,Ny,visucase,step,filename)
    sll_int32,intent(in)::Nx,Ny,visucase,step
    sll_real64,dimension(0:1,0:1),intent(in)::dom
    sll_real64,dimension(0:Nx,0:Ny),intent(in)::ftab
    character(len=*),intent(in)::filename

    if(visucase==0)then
       !gnuplot
       call printgp2d(dom,ftab,Nx,Ny,step,filename)
    endif
    if(visucase==1)then
       !vtk
       call printvtk2d(dom,ftab,Nx,Ny,step,filename)
    endif
  end subroutine print2d



  subroutine print2dper(dom,ftab,Nx,Ny,visucase,step,filename)
    sll_int32,intent(in)::Nx,Ny,visucase,step
    sll_real64,dimension(0:1,0:1),intent(in)::dom
    sll_real64,dimension(0:Nx-1,0:Ny-1),intent(in)::ftab
    character(len=*),intent(in)::filename

    if(visucase==0)then
       !gnuplot
       call printgp2dper(dom,ftab,Nx,Ny,step,filename)
    endif
    if(visucase==1)then
       !vtk
       call printvtk2dper(dom,ftab,Nx,Ny,step,filename)
    endif
  end subroutine print2dper

  subroutine printgp2dper(dom,ftab,Nx,Ny,step,filename)
    sll_int32,intent(in)::Nx,Ny,step
    sll_real64,dimension(0:1,0:1),intent(in)::dom
    sll_real64,dimension(0:Nx-1,0:Ny-1),intent(in)::ftab
    sll_int32::i,j
    sll_real64::z(0:1),dz(0:1)
    character(len=*),intent(in)::filename
    character(len=80)::str,str2
    write(str2,*)step
    str=trim(adjustl((filename)))//trim(adjustl((str2)))//'.dat'

    dz(0)=(dom(1,0)-dom(0,0))/real(Nx,f64);dz(1)=(dom(1,1)-dom(0,1))/real(Ny,f64)
    open(unit=900,file=str)
    do j=0,Ny-1
       do i=0,Nx-1
          z(0)=dom(0,0)+real(i,f64)*dz(0)
          z(1)=dom(0,1)+real(j,f64)*dz(1)
          write(900,*) z(0),z(1),ftab(i,j)
       enddo
       i=Nx
       z(0)=dom(0,0)+real(i,f64)*dz(0)
       z(1)=dom(0,1)+real(j,f64)*dz(1)
       write(900,*) z(0),z(1),ftab(0,j)      
       write(900,*) ''      
    enddo
    j=Ny
    do i=0,Nx-1
       z(0)=dom(0,0)+real(i,f64)*dz(0)
       z(1)=dom(0,1)+real(j,f64)*dz(1)
       write(900,*) z(0),z(1),ftab(i,0)
    enddo
    i=Nx
    z(0)=dom(0,0)+real(i,f64)*dz(0)
    z(1)=dom(0,1)+real(j,f64)*dz(1)
    write(900,*)z(0),z(1),ftab(0,0)
    write(900,*)''
    close(900)  
  end subroutine printgp2dper


  subroutine printgp2d(dom,ftab,Nx,Ny,step,filename)
    sll_int32,intent(in)::Nx,Ny,step
    sll_real64,dimension(0:1,0:1),intent(in)::dom
    sll_real64,dimension(0:Nx,0:Ny),intent(in)::ftab
    sll_int32::i,j
    sll_real64::z(0:1),dz(0:1)
    character(len=*),intent(in)::filename
    character(len=80)::str,str2
    write(str2,*)step
    str=trim(adjustl((filename)))//trim(adjustl((str2)))//'.dat'

    dz(0)=(dom(1,0)-dom(0,0))/real(Nx,f64);dz(1)=(dom(1,1)-dom(0,1))/real(Ny,f64)
    open(unit=900,file=str)
    do j=0,Ny
       do i=0,Nx
          z(0)=dom(0,0)+real(i,f64)*dz(0)
          z(1)=dom(0,1)+real(j,f64)*dz(1)
          write(900,*) z(0),z(1),ftab(i,j)
       enddo
       write(900,*) ''      
    enddo
    close(900)  
  end subroutine printgp2d


  subroutine printvtk2dper(dom,ftab,Nx,Ny,step,filename)
    sll_int32,intent(in)::Nx,Ny
    sll_real64,dimension(0:1,0:1),intent(in)::dom
    sll_real64,dimension(0:Nx-1,0:Ny-1),intent(in)::ftab
    sll_int32::i,j
    sll_int32,intent(in):: step
    sll_real64::z(0:1),dz(0:1)
    character(len=*),intent(in)::filename
    character(len=80)::str,str2
    write(str2,*)step
    !write(str,*) 'mv f.dat f'//trim(adjustl((str2)))//'.dat';call system(str)
    write(str,*) 'f'//trim(adjustl((filename)))//trim(adjustl((str2)))//'.vtk';!call system(str)
    str=trim(adjustl((filename)))//trim(adjustl((str2)))//'.vtk';!call system(str)
    dz(0)=(dom(1,0)-dom(0,0))/real(Nx,f64);dz(1)=(dom(1,1)-dom(0,1))/real(Ny,f64)
    !open(unit=900,file='f.vtk')
    open(unit=900,file=str,form='formatted')
    write(900,'(A)')                  '# vtk DataFile Version 2.0'
    write(900,'(A)')                  'Exemple'
    write(900,'(A)')                  'ASCII'
    write(900,'(A)')                  'DATASET STRUCTURED_POINTS'
    write(900,'(A,I0,A,I0,A,I0)') 'DIMENSIONS ', Nx+1,' ', Ny+1,' ', 1
    write(900,'(A,I0,A,I0,A,I0)') 'ORIGIN ', floor(dom(0,0)+0.1),' ' , floor(dom(0,1)+0.1),' ' , 0
    !write(900,'(A,F10.4,A,F10.4,A,F10.4)') 'SPACING ', dz(0),' ', dz(1),' ', 1. 
    write(900,*) 'SPACING ', dz(0),' ', dz(1),' ', 1. 
    write(900,*)
    write(900,'(A,I0)')           'POINT_DATA ',(Nx+1)*(Ny+1)
    write(900,'(A,I0)')           'SCALARS f float ',1
    write(900,'(A)')                  'LOOKUP_TABLE default'

    do j=0,Ny-1
       do i=0,Nx-1
          z(0)=dom(0,0)+real(i,f64)*dz(0)
          z(1)=dom(0,1)+real(j,f64)*dz(1)
          !write(900,'(F0.8)') ftab(i,j)
          write(900,*) ftab(i,j)
       enddo
       i=Nx
       z(0)=dom(0,0)+real(i,f64)*dz(0)
       z(1)=dom(0,1)+real(j,f64)*dz(1)
       !write(900,'(F0.8)') ftab(0,j)            
       write(900,*) ftab(0,j)            
    enddo
    j=Ny
    do i=0,Nx-1
       z(0)=dom(0,0)+real(i,f64)*dz(0)
       z(1)=dom(0,1)+real(j,f64)*dz(1)
       !write(900,'(F0.8)') ftab(i,0)
       write(900,*) ftab(i,0)
    enddo
    i=Nx
    z(0)=dom(0,0)+real(i,f64)*dz(0)
    z(1)=dom(0,1)+real(j,f64)*dz(1)
    !write(900,'(F0.8)') ftab(0,0)	  	       
    write(900,*) ftab(0,0)
    close(900)  
  end subroutine printvtk2dper

  subroutine printvtk2d(dom,ftab,Nx,Ny,step,filename)
    sll_int32,intent(in)::Nx,Ny
    sll_real64,dimension(0:1,0:1),intent(in)::dom
    sll_real64,dimension(0:Nx,0:Ny),intent(in)::ftab
    sll_int32::i,j
    sll_int32,intent(in):: step
    sll_real64::z(0:1),dz(0:1)
    character(len=*),intent(in)::filename
    character(len=80)::str,str2
    write(str2,*)step
    !write(str,*) 'mv f.dat f'//trim(adjustl((str2)))//'.dat';call system(str)
    write(str,*) 'f'//trim(adjustl((filename)))//trim(adjustl((str2)))//'.vtk';!call system(str)
    str=trim(adjustl((filename)))//trim(adjustl((str2)))//'.vtk';!call system(str)
    dz(0)=(dom(1,0)-dom(0,0))/real(Nx,f64);dz(1)=(dom(1,1)-dom(0,1))/real(Ny,f64)
    !open(unit=900,file='f.vtk')
    open(unit=900,file=str,form='formatted')
    write(900,'(A)')                  '# vtk DataFile Version 2.0'
    write(900,'(A)')                  'Exemple'
    write(900,'(A)')                  'ASCII'
    write(900,'(A)')                  'DATASET STRUCTURED_POINTS'
    write(900,'(A,I0,A,I0,A,I0)') 'DIMENSIONS ', Nx+1,' ', Ny+1,' ', 1
    write(900,'(A,I0,A,I0,A,I0)') 'ORIGIN ', floor(dom(0,0)+0.1),' ' , floor(dom(0,1)+0.1),' ' , 0
    !write(900,'(A,F10.4,A,F10.4,A,F10.4)') 'SPACING ', dz(0),' ', dz(1),' ', 1. 
    write(900,*) 'SPACING ', dz(0),' ', dz(1),' ', 1. 
    write(900,*)
    write(900,'(A,I0)')           'POINT_DATA ',(Nx+1)*(Ny+1)
    write(900,'(A,I0)')           'SCALARS f float ',1
    write(900,'(A)')                  'LOOKUP_TABLE default'

    do j=0,Ny
       do i=0,Nx
          z(0)=dom(0,0)+real(i,f64)*dz(0)
          z(1)=dom(0,1)+real(j,f64)*dz(1)
          !write(900,'(F0.8)') ftab(i,j)
          write(900,*) ftab(i,j)
       enddo
    enddo
    close(900)  
  end subroutine printvtk2d

  subroutine advect2d_CSL_CG(dom,f,N0,N1,buf2d,interp_case,carac,ppm_order) !conservative 2d remapping algorithm
    sll_real64,dimension(0:1,0:1),intent(in)::dom
    sll_int32,intent(in)::N0,N1,interp_case,ppm_order
    sll_real64,dimension(0:N0,0:N1-1)::buf2d
    sll_real64,dimension(1:N0+1,1:N1)::f
    sll_real64,dimension(2,-1:N0+1,-1:N1)::carac
    sll_real64::xx(4),yy(4),xA,yA,xB,yB,res,xx0,yy0,x,xxn,y,yyn,dx,dy
    sll_real64::xA_loc,yA_loc,xB_loc,yB_loc
    sll_int32::i,j,ii(4),jj(4),im1,jm1,i0,j0,i1,j1,s,k,sx,sy,minfl,maxfl,ix,iii,iy,iiii,jjjj
    sll_real64,dimension(:),allocatable::intx,inty
    sll_int32,dimension(:,:),allocatable::tnbr
    sll_int32,dimension(:,:,:),allocatable::cell
    sll_real64,dimension(:,:),allocatable::tt,tcell,dir,aretesh,aretesv,sommets,aretesvg,aretesvd,areteshb,areteshh,&
    sommetsbg,sommetsbd,sommetshg,sommetshd
    sll_real64,dimension(:,:,:,:),allocatable::tpts
    sll_int32::nbx(4),nby(4),nbmax,dirx,diry,ell,ell1,imin,jmin,imax,jmax,i0_loc,j0_loc
    sll_real64::xx1,xx2,yy1,yy2,w00,w10,w01,w20,w02,w11,w21,&
		w12,w22,c00,c10,c01,c20,c02,c11,c12,c21,c22,&
		fij,fim1jm1,fim1j,fim1jp1,fijm1,fijp1,fip1jm1,fip1j,fip1jp1,xxx,yyy
    sll_int32::im2,ib,ip1,ip2,jm2,jb,jp1,jp2   

    nbmax=30
    !print*,f(1,33)
    !return
    allocate(tt(nbmax,2),cell(2,nbmax,4),tcell(nbmax,4),intx(0:nbmax),inty(0:nbmax),dir(nbmax,2))
    do j=0,N1-1
      do i=0,N0
        buf2d(i,j)=f(i+1,j+1)
      enddo
    enddo
    dx=dom(1,0)/real(N0,f64)
    dy=dom(1,1)/real(N1,f64)


    if ((interp_case==3) .or. (interp_case==6) .or. (interp_case==7)) then
      allocate(aretesh(0:N0,0:N1-1),aretesv(0:N0,0:N1-1),sommets(0:N0,0:N1-1))
      allocate(aretesvg(0:N0,0:N1-1),aretesvd(0:N0,0:N1-1),areteshb(0:N0,0:N1-1),areteshh(0:N0,0:N1-1),&
      sommetsbg(0:N0,0:N1-1),sommetsbd(0:N0,0:N1-1),sommetshg(0:N0,0:N1-1),sommetshd(0:N0,0:N1-1))
      call aux(N0,N1,buf2d,aretesh,aretesv,sommets,ppm_order,dom)
    endif
    
    if (interp_case==4) then
      allocate(aretesvg(0:N0,0:N1-1),aretesvd(0:N0,0:N1-1),areteshb(0:N0,0:N1-1),areteshh(0:N0,0:N1-1),&
      sommetsbg(0:N0,0:N1-1),sommetsbd(0:N0,0:N1-1),sommetshg(0:N0,0:N1-1),sommetshd(0:N0,0:N1-1))
      call aux2(N0,N1,buf2d,areteshb,areteshh,aretesvg,aretesvd,sommetsbg,sommetsbd,sommetshg,sommetshd,ppm_order,carac,dom)
    endif

    do j=0,N1-1
      do i=0,N0

	f(i+1,j+1)=0._f64

        !im1=modulo(i-1,N0)
	!i1=modulo(i+1,N0)
	!jm1=modulo(j-1,N1)
	!j1=modulo(j+1,N1)
        
	!computation of the feet of the characteristics
        xx(1)=0.25_f64*(carac(1,i-1,j-1)+carac(1,i,j-1)+carac(1,i,j)+carac(1,i-1,j))
        xx(2)=0.25_f64*(carac(1,i,j-1)+carac(1,i+1,j-1)+carac(1,i+1,j)+carac(1,i,j))
        xx(3)=0.25_f64*(carac(1,i,j)+carac(1,i+1,j)+carac(1,i+1,j+1)+carac(1,i,j+1))
        xx(4)=0.25_f64*(carac(1,i-1,j)+carac(1,i,j)+carac(1,i,j+1)+carac(1,i-1,j+1))

        yy(1)=0.25_f64*(carac(2,i-1,j-1)+carac(2,i,j-1)+carac(2,i,j)+carac(2,i-1,j))
        yy(2)=0.25_f64*(carac(2,i,j-1)+carac(2,i+1,j-1)+carac(2,i+1,j)+carac(2,i,j))
        yy(3)=0.25_f64*(carac(2,i,j)+carac(2,i+1,j)+carac(2,i+1,j+1)+carac(2,i,j+1))
        yy(4)=0.25_f64*(carac(2,i-1,j)+carac(2,i,j)+carac(2,i,j+1)+carac(2,i-1,j+1))      

!if ((i==0) .and. (j==32)) then
!print*,'ok',i,j,carac(2,i,j)
!stop
!endif
!print*,i,j,xx(1),yy(1),xx(2),yy(2),xx(3),yy(3),xx(4),yy(4)

        !normalization
	xx(1)=(xx(1)-dom(0,0))/dom(1,0)*real(N0,f64)
	xx(2)=(xx(2)-dom(0,0))/dom(1,0)*real(N0,f64)
	xx(3)=(xx(3)-dom(0,0))/dom(1,0)*real(N0,f64)
	xx(4)=(xx(4)-dom(0,0))/dom(1,0)*real(N0,f64)
	
	yy(1)=(yy(1)-dom(0,1))/dom(1,1)*real(N1,f64)
	yy(2)=(yy(2)-dom(0,1))/dom(1,1)*real(N1,f64)
	yy(3)=(yy(3)-dom(0,1))/dom(1,1)*real(N1,f64)
	yy(4)=(yy(4)-dom(0,1))/dom(1,1)*real(N1,f64)

        xx=xx+0.5_f64
        yy=yy+0.5_f64 

!if ((i==0) .and. (j==32)) then
!print*,i,j,xx(1),yy(1),xx(2),yy(2),xx(3),yy(3),xx(4),yy(4)
!endif

	
        
        ii(1)=floor(xx(1))
        ii(2)=floor(xx(2))
        ii(3)=floor(xx(3))
        ii(4)=floor(xx(4))
	
        jj(1)=floor(yy(1))
        jj(2)=floor(yy(2))
        jj(3)=floor(yy(3))
        jj(4)=floor(yy(4))
	
	imin=min(ii(1),ii(2),ii(3),ii(4))
	jmin=min(jj(1),jj(2),jj(3),jj(4))
	imax=max(ii(1),ii(2),ii(3),ii(4))
	jmax=max(jj(1),jj(2),jj(3),jj(4))

!        allocate(tnbr(imin:imax,jmin:jmax))
!        allocate(tpts(imin:imax,jmin:jmax,100,2))

!        tnbr=0    
	!Computation of external edges

!#ifdef COMPUTE_EDGE  

    do ell=1,4

	ell1=ell+1;if(ell1==5)ell1=1  
	xA=xx(ell);yA=yy(ell);xB=xx(ell1);yB=yy(ell1)
	i0=ii(ell);j0=jj(ell);i1=ii(ell1);j1=jj(ell1)

	s=1;
	if(i0<i1)then
	  do k=i0+1,i1
	    tt(s,1)=(real(k,f64)-xA)/(xB-xA)
	    s=s+1
	  enddo	  
	  dirx=1
	endif
	if(i0>i1)then
	  do k=i0,i1+1,-1
	    tt(s,1)=(real(k,f64)-xA)/(xB-xA)
	    s=s+1
	  enddo
	  dirx=-1
	endif
        nbx(ell)=s-1;
	s=1;
	if(j0<j1)then
	  do k=j0+1,j1
	    tt(s,2)=(real(k,f64)-yA)/(yB-yA)
	    s=s+1
	  enddo
	  diry=1
	endif
	if(j0>j1)then
	  do k=j0,j1+1,-1
	    tt(s,2)=(real(k,f64)-yA)/(yB-yA)
	    s=s+1
	  enddo
	  diry=-1
	endif
	nby(ell)=s-1
	
	cell(1,1,ell)=i0
	cell(2,1,ell)=j0
	tcell(1,ell)=0._f64
	sx=1;sy=1
	s=1

	do while((sx<=nbx(ell)).and.(sy<=nby(ell)))
	  if(tt(sx,1)<tt(sy,2))then
	    s=s+1
	    cell(1,s,ell)=cell(1,s-1,ell)+dirx
	    cell(2,s,ell)=cell(2,s-1,ell)
	    tcell(s,ell)=tt(sx,1)
	    intx(2*cell(1,s-1,ell)+(dirx-1)/2-2*imin)=yA+tt(sx,1)*(yB-yA)
	    dir(s,1)=1
	    dir(s,2)=dirx  
	    sx=sx+1	
	  else
	    s=s+1
	    cell(1,s,ell)=cell(1,s-1,ell)
	    cell(2,s,ell)=cell(2,s-1,ell)+diry
	    tcell(s,ell)=tt(sy,2)
	    inty(2*cell(2,s-1,ell)+(diry-1)/2-2*jmin)=xA+tt(sy,2)*(xB-xA)
	    dir(s,1)=2
	    dir(s,2)=diry	  
	    sy=sy+1	    
	  endif
	enddo
	do while(sx<=nbx(ell))
	  s=s+1
	  cell(1,s,ell)=cell(1,s-1,ell)+dirx
	  cell(2,s,ell)=cell(2,s-1,ell)
	  tcell(s,ell)=tt(sx,1)
	  intx(2*cell(1,s-1,ell)+(dirx-1)/2-2*imin)=yA+tt(sx,1)*(yB-yA)
	  dir(s,1)=1
	  dir(s,2)=dirx
	  sx=sx+1
	enddo  
	do while(sy<=nby(ell))
	  s=s+1
	  cell(1,s,ell)=cell(1,s-1,ell)
	  cell(2,s,ell)=cell(2,s-1,ell)+diry
	  tcell(s,ell)=tt(sy,2)
	  inty(2*cell(2,s-1,ell)+(diry-1)/2-2*jmin)=xA+tt(sy,2)*(xB-xA)
	  dir(s,1)=2
	  dir(s,2)=diry	  
	  sy=sy+1
	enddo        
!#endif

!Computation of extern edges

!#ifdef COMPUTE

if ((ell==1) .or. (ell==4) .or. ((ell==2) .and. (i==N0)) .or. ((ell==3) .and. (j==N1-1))) then
        xB_loc=xx(ell)
        yB_loc=yy(ell)
        do k=2,nbx(ell)+nby(ell)+2
          xA_loc=xB_loc
          yA_loc=yB_loc
          i0_loc=cell(1,k-1,ell)
          j0_loc=cell(2,k-1,ell)
          if(k==nbx(ell)+nby(ell)+2)then
            xB_loc=xx(ell1)
            yB_loc=yy(ell1)
          else
            if(dir(k,1)==1)then
              xB_loc=real(cell(1,k,ell)+(1-dir(k,2))/2,f64)
              yB_loc=yy(ell)+tcell(k,ell)*(yy(ell1)-yy(ell))        	  
            else
              xB_loc=xx(ell)+tcell(k,ell)*(xx(ell1)-xx(ell))					
              yB_loc=real(cell(2,k,ell)+(1-dir(k,2))/2,f64)       	  
            endif
          endif

          call calcule_coeff(N0,N1,buf2d,i0_loc,j0_loc,xA_loc,yA_loc,xB_loc,yB_loc,res,aretesh,aretesv,sommets,areteshb,areteshh,&
          aretesvg,aretesvd,sommetsbg,sommetsbd,sommetshg,sommetshd,interp_case)

	  f(i+1,j+1)=f(i+1,j+1)+res

          if ((i-1>=0) .and. (ell==4)) then
            f(i,j+1)=f(i,j+1)-res
          endif
          if ((j-1>=0) .and. (ell==1)) then
            f(i+1,j)=f(i+1,j)-res
          endif
        enddo
endif
        
      end do  	

!#endif

!Computation of vertical intern edges

!#ifdef COMPUTE_V
	
	do ell=0,imax-imin-1
	minfl=min(floor(intx(2*ell)),floor(intx(2*ell+1)))
	maxfl=max(floor(intx(2*ell)),floor(intx(2*ell+1)))

        i0_loc=imin+ell
        yB_loc=min(intx(2*ell),intx(2*ell+1))
        do k=0,maxfl-minfl
          yA_loc=yB_loc
          j0_loc=minfl+k  
          if(k==maxfl-minfl)then
            yB_loc=max(intx(2*ell),intx(2*ell+1))
          else
            yB_loc=real(minfl+k+1,f64)
          endif

          call calcule_coeffv(N0,N1,buf2d,i0_loc,j0_loc,yA_loc,yB_loc,res,aretesh,aretesv,sommets,areteshb,areteshh,&
	  aretesvg,aretesvd,sommetsbg,sommetsbd,sommetshg,sommetshd,interp_case)
          f(i+1,j+1)=f(i+1,j+1)+res

	  enddo
	enddo

!#endif

!Computation of horizontal intern edges

!#ifdef COMPUTE_X
	
	do ell=0,jmax-jmin-1
	minfl=min(floor(inty(2*ell)),floor(inty(2*ell+1)))
	maxfl=max(floor(inty(2*ell)),floor(inty(2*ell+1)))

        j0_loc=jmin+ell
        xA_loc=min(inty(2*ell),inty(2*ell+1))
        do k=0,maxfl-minfl
          xB_loc=xA_loc
          i0_loc=minfl+k  
          if(k==maxfl-minfl)then
            xA_loc=max(inty(2*ell),inty(2*ell+1))
          else
            xA_loc=real(minfl+k+1,f64)
          endif
          
          call calcule_coeffh(N0,N1,buf2d,i0_loc,j0_loc,xA_loc,xB_loc,res,aretesh,aretesv,sommets,areteshb,areteshh,&
	  aretesvg,aretesvd,sommetsbg,sommetsbd,sommetshg,sommetshd,interp_case)
	  f(i+1,j+1)=f(i+1,j+1)+res

	  enddo
	enddo

!#endif
!        deallocate(tnbr)
!        deallocate(tpts)

      enddo

    enddo
   
!#ifdef DIAG_TIME
!      f=buf2d
!#endif

    deallocate(tt,cell,tcell,intx,inty,dir)

    if ((interp_case==3) .or. (interp_case==6) .or. (interp_case==7)) then
      deallocate(aretesh,aretesv,sommets)
      deallocate(areteshb,areteshh,aretesvg,aretesvd,sommetsbg,sommetsbd,sommetshg,sommetshd)
    endif
    if (interp_case==4) then
      deallocate(areteshb,areteshh,aretesvg,aretesvd,sommetsbg,sommetsbd,sommetshg,sommetshd)
    endif

   

  end subroutine advect2d_CSL_CG

	subroutine aux(N0,N1,f,aretesh,aretesv,sommets,ordre,dom) !used in PPM case

		sll_int32,intent(in)::N0,N1,ordre
                sll_real64,dimension(0:1,0:1),intent(in)::dom
		real(f64),dimension(0:N0,0:N1-1),intent(in)::f
            !    sll_real64,dimension(2,-1:N0+1,-1:N1+1),intent(in)::carac
		sll_int32::i,j,im3,im2,im1,ib,ip1,ip2,jm3,jm2,jm1,jb,jp1,jp2
		real(f64),dimension(0:N0,0:N1-1),intent(inout)::aretesh,aretesv,sommets
!		print*,N0,N1,ordre
!stop

		do i=0,N0
			do j=0,N1-1
				im3=i-3
				im2=i-2
				im1=i-1
				ip1=i+1
				ip2=i+2
				if(i-3<=0)im3=0;
				if(i-2<=0)im2=0;
				if(i-1<=0)im1=0;
				ib=i;
				if(i+1>=N0)ip1=N0;
				if(i+2>=N0)ip2=N0;
				
				!im3=modulo(i-3,N0)
				!im2=modulo(i-2,N0)
				!im1=modulo(i-1,N0)
				!ib=modulo(i,N0)
				!ip1=modulo(i+1,N0)
				!ip2=modulo(i+2,N0)
				jm3=modulo(j-3,N1)
				jm2=modulo(j-2,N1)
				jm1=modulo(j-1,N1)
				jb=modulo(j,N1)
				jp1=modulo(j+1,N1)
				jp2=modulo(j+2,N1)

				if (ordre==1) then !PPM1
					aretesv(i,j)=7._f64/12._f64*(f(im1,jb)+f(ib,jb)) &
					-1._f64/12._f64*(f(im2,jb)+f(ip1,jb))
					aretesh(i,j)=7._f64/12._f64*(f(ib,jm1)+f(ib,jb)) &	
					-1._f64/12._f64*(f(ib,jm2)+f(ib,jp1))
				else if (ordre==2) then !PPM2
					aretesv(i,j)=1._f64/60._f64*(f(ip2,jb)+f(im3,jb)) &
					-8._f64/60._f64*(f(ip1,jb)+f(im2,jb)) &
					+37._f64/60._f64*(f(ib,jb)+f(im1,jb))
					aretesh(i,j)=1._f64/60._f64*(f(ib,jp2)+f(ib,jm3)) &
					-8._f64/60._f64*(f(ib,jp1)+f(ib,jm2)) &
					+37._f64/60._f64*(f(ib,jb)+f(ib,jm1))
				else if (ordre==0) then !PPM0
					aretesv(i,j)=1._f64/2._f64*(f(ib,jb)+f(im1,jb))
					aretesh(i,j)=1._f64/2._f64*(f(ib,jb)+f(ib,jm1))
				end if
			end do
		end do	

		do i=0,N0
			do j=0,N1-1
				!im3=modulo(i-3,N0)
				!im2=modulo(i-2,N0)
				!im1=modulo(i-1,N0)
				!ib=modulo(i,N0)
				!ip1=modulo(i+1,N0)
				!ip2=modulo(i+2,N0)
				im3=i-3
				im2=i-2
				im1=i-1
				ip1=i+1
				ip2=i+2
				if(i-3<=0)im3=0;
				if(i-2<=0)im2=0;
				if(i-1<=0)im1=0;
				ib=i;
				if(i+1>=N0)ip1=N0;
				if(i+2>=N0)ip2=N0;
				
				
				
				jm3=modulo(j-3,N1)
				jm2=modulo(j-2,N1)
				jm1=modulo(j-1,N1)
				jb=modulo(j,N1)
				jp1=modulo(j+1,N1)
				jp2=modulo(j+2,N1)

				if (ordre==1) then !PPM1
					sommets(i,j)=7._f64/12._f64*(aretesv(ib,jm1)+aretesv(ib,jb)) &
					-1._f64/12._f64*(aretesv(ib,jm2)+aretesv(ib,jp1))
				else if (ordre==2) then !PPM2
					sommets(i,j)=1._f64/60._f64*(aretesv(ib,jp2)+aretesv(ib,jm3)) &
					-8._f64/60._f64*(aretesv(ib,jp1)+aretesv(ib,jm2)) &
					+37._f64/60._f64*(aretesv(ib,jb)+aretesv(ib,jm1))
				else if (ordre==0) then !PPM0
					sommets(i,j)=1._f64/2._f64*(aretesv(ib,jb)+aretesv(ib,jm1))
				end if
			end do
		end do

	end subroutine aux


	subroutine aux2(N0,N1,f,areteshb,areteshh,aretesvg,aretesvd,sommetsbg,sommetsbd,sommetshg,sommetshd,ordre,carac,dom) !used in PPM case

		integer,intent(in)::N0,N1,ordre
                real(f64),dimension(0:1,0:1),intent(in)::dom
		real(f64),dimension(0:N0,0:N1-1),intent(in)::f
                real(f64),dimension(2,-1:N0+1,-1:N1+1),intent(in)::carac
		integer::i,j,im3,im2,im1,ib,ip1,ip2,jm3,jm2,jm1,jb,jp1,jp2,i3
		real(f64),dimension(0:N0,0:N1-1),intent(inout)::areteshb,areteshh,aretesvg,aretesvd,&
                sommetsbg,sommetsbd,sommetshg,sommetshd
    real(f64) ::w(-ordre:ordre+1),tmp,ww(-ordre:ordre)
    
    !f2py intent(in)::buf,f
    integer::r,s,ii,d   
!		print*,N0,N1,ordre
!stop

    d=ordre
    r=-d
    s=d+1
    
    !maple code for generation of w
    !for k from r to -1 do
    !  C[k]:=product((k-j),j=r..k-1)*product((k-j),j=k+1..s):
    !  C[k]:=1/C[k]*product((-j),j=r..k-1)*product((-j),j=k+1..-1)*product((-j),j=1..s):
    !od:
    !for k from 1 to s do
    !  C[k]:=product((k-j),j=r..k-1)*product((k-j),j=k+1..s):
    !  C[k]:=1/C[k]*product((-j),j=r..-1)*product((-j),j=1..k-1)*product((-j),j=k+1..s):
    !od:
    !C[0]:=-add(C[k],k=r..-1)-add(C[k],k=1..s):
    
    do i=r,-1
      tmp=1._f64
      do j=r,i-1
        tmp=tmp*real(i-j,f64)
      enddo
      do j=i+1,s
        tmp=tmp*real(i-j,f64)
      enddo
      tmp=1._f64/tmp
      do j=r,i-1
        tmp=tmp*real(-j,f64)
      enddo
      do j=i+1,-1
        tmp=tmp*real(-j,f64)
      enddo
      do j=1,s
        tmp=tmp*real(-j,f64)
      enddo
      w(i)=tmp      
    enddo

!    do i=r,-1
!      tmp=1._f64
!      !do j=r,i-1
!      !  tmp=tmp*real(i-j,f64)
!      !enddo
!      !do j=i+1,s
!      !  tmp=tmp*real(i-j,f64)
!      !enddo
!      !tmp=1._f64/tmp
!      do j=r,i-1 !-j/(i-j)=j/(j-i)=1/(1-i/j)
!        tmp=tmp*(1._f64-real(i,f64)/real(j,f64))
!      enddo
!      do j=i+1,-1
!        tmp=tmp*(1._f64-real(i,f64)/real(j,f64))
!      enddo
!      do j=1,s
!        tmp=tmp*(1._f64-real(i,f64)/real(j,f64))
!      enddo
!      tmp=tmp*real(i,f64)
!      w(i)=1._f64/tmp      
!    enddo



    do i=1,s
      tmp=1._f64
      do j=r,i-1
        tmp=tmp*real(i-j,f64)
      enddo
      do j=i+1,s
        tmp=tmp*real(i-j,f64)
      enddo
      tmp=1._f64/tmp
      do j=r,-1
        tmp=tmp*real(-j,f64)
      enddo
      do j=1,i-1
        tmp=tmp*real(-j,f64)
      enddo
      do j=i+1,s
        tmp=tmp*real(-j,f64)
      enddo
      w(i)=tmp      
    enddo

    tmp=0._f64
    do i=r,-1
      tmp=tmp+w(i)
    enddo
    do i=1,s
      tmp=tmp+w(i)
    enddo
    w(0)=-tmp
    
    
    
    !print *,'w',w
    !do ii=r,s
    !  print *,ii,w(r+s-ii)
    !enddo
    
    !compute now ww
    !maple code
    !#for conservative formulation
    !tmp:=0:
    !for k from r to -1 do
    !tmp:=tmp+C[k]:
    !CC[k]:=-tmp:
    !od:
    !tmp:=0:
    !for k from s to 1 by -1 do
    !  tmp:=tmp+C[k]:
    !  CC[k-1]:=tmp:
    !od:
    !seq(CC[k],k=r..s-1);
    !evalf(%);

    tmp=0._f64
    do i=r,-1
      tmp=tmp+w(i)
      ww(i)=-tmp
    enddo
    tmp=0._f64
    do i=s,1,-1
      tmp=tmp+w(i)
      ww(i-1)=tmp
    enddo

    !print *,'ww',ww
    !do ii=r,s-1
    !  print *,ii,ww(r+s-1-ii)
    !enddo
    !stop
    do j=0,N1-1
      do i=0,N0
        tmp=0._f64
        do ii=r,s-1
          i3=i+ii-1;if(i3<=0)i3=0;if(i3>=N0)i3=N0          
          tmp=tmp+ww(r+s-1-ii)*f(i3,j)
        enddo
        aretesvg(i,j)=tmp
        tmp=0._f64
        do ii=r,s-1
          i3=i+ii;if(i3<=0)i3=0;if(i3>=N0)i3=N0          
          tmp=tmp+ww(ii)*f(i3,j)
        enddo
        aretesvd(i,j)=tmp
        tmp=0._f64
        do ii=r,s-1
          tmp=tmp+ww(ii)*f(i,modulo(j+ii,N1))
        enddo
        areteshh(i,j)=tmp
        tmp=0._f64
        do ii=r,s-1
          tmp=tmp+ww(r+s-1-ii)*f(i,modulo(j+ii-1,N1))
        enddo
        areteshb(i,j)=tmp
      enddo
    enddo

    do j=0,N1-1
      do i=0,N0
        tmp=0._f64
        do ii=r,s-1
          i3=i+ii;if(i3<=0)i3=0;if(i3>=N0)i3=N0          
          tmp=tmp+ww(ii)*areteshh(i3,j)
        enddo
        sommetshd(i,j)=tmp
        tmp=0._f64
        do ii=r,s-1
          i3=i+ii-1;if(i3<=0)i3=0;if(i3>=N0)i3=N0          
          tmp=tmp+ww(r+s-1-ii)*areteshh(i3,j)
        enddo
        sommetshg(i,j)=tmp
        tmp=0._f64
        do ii=r,s-1
          i3=i+ii;if(i3<=0)i3=0;if(i3>=N0)i3=N0          
          tmp=tmp+ww(ii)*areteshb(i3,j)
        enddo
        sommetsbd(i,j)=tmp
        tmp=0._f64
        do ii=r,s-1
          i3=i+ii-1;if(i3<=0)i3=0;if(i3>=N0)i3=N0          
          tmp=tmp+ww(r+s-1-ii)*areteshb(i3,j)
        enddo
        sommetsbg(i,j)=tmp
      enddo
    enddo

	end subroutine aux2
	
	
	subroutine calcule_coeff(N0,N1,a_moyenne,i,j,x1,y1,x2,y2,res,aretesh,aretesv,sommets,areteshb,areteshh,aretesvg,aretesvd,&
                sommetsbg,sommetsbd,sommetshg,sommetshd,cas)

	        sll_int32,intent(in)::N0,N1,i,j,cas
		real(f64),intent(in)::x1,y1,x2,y2
		real(f64),dimension(0:N0,0:N1-1),intent(in)::a_moyenne,aretesh,aretesv,sommets,areteshb,areteshh,aretesvg,aretesvd,&
                sommetsbg,sommetsbd,sommetshg,sommetshd
		real(f64)::xx1,xx2,yy1,yy2,aux,dax,day,bx,by,c,w00,w10,w01,w20,w02,w11,w21,w12,w22,&
                c00,c10,c01,c20,c02,c11,c12,c21,c22,sij,sip1j,sijp1,sip1jp1,avij,avip1j,ahij,ahijp1,&
		fij,fim1jm1,fim1j,fim1jp1,fijm1,fijp1,fip1jm1,fip1j,fip1jp1
		sll_int32::im2,im1,ib,ip1,ip2,jm2,jm1,jb,jp1,jp2
		real(f64),intent(out)::res
		
		if(cas==0)then
		  return
		endif


		!im2=modulo(i-2,N0)
		!im1=modulo(i-1,N0)
		!ib=modulo(i,N0)
		!ip1=modulo(i+1,N0)
		!ip2=modulo(i+2,N0)
		im2=i-2
		im1=i-1
		ib=i
		ip1=i+1
		ip2=i+2
		
		if(i-2<=0)im2=0;
		if(i-1<=0)im1=0;
		if(i<=0)ib=0;
		if(i+1<=0)ip1=0;
		if(i+2<=0)ip2=0;

		if(i-2>=N0)im2=N0;
		if(i-1>=N0)im1=N0;		
		if(i>=N0)ib=N0;
		if(i+1>=N0)ip1=N0;
		if(i+2>=N0)ip2=N0;

		
		jm2=modulo(j-2,N1)
		jm1=modulo(j-1,N1)
		jb=modulo(j,N1)
		jp1=modulo(j+1,N1)
		jp2=modulo(j+2,N1)

!1: Lauritzen, 2: Lag3, 3: PPM

		xx1=x1-real(i,f64)
		xx2=x2-real(i,f64)
		yy1=y1-real(j,f64)
		yy2=y2-real(j,f64)

		w00=1._f64/2._f64*(xx2+xx1)*(yy2-yy1)
		w10=1._f64/6._f64*(xx2**2+xx2*xx1+xx1**2)*(yy2-yy1)
		w01=-1._f64/6._f64*(yy2**2+yy2*yy1+yy1**2)*(xx2-xx1)
		w20=1._f64/12._f64*(xx2+xx1)*(xx2**2+xx1**2)*(yy2-yy1)
		w02=-1._f64/12._f64*(yy2+yy1)*(yy2**2+yy1**2)*(xx2-xx1)
		w11=1._f64/24._f64*(yy2*(3._f64*xx2**2+2._f64*xx2*xx1+xx1**2)+yy1*(xx2**2+2._f64*xx2*xx1+3._f64*xx1**2))*(yy2-yy1)
		w21=-1._f64/60._f64*(6._f64*xx1**2*yy1**2+6._f64*xx2**2*yy2**2+xx1**2*yy2**2+xx2**2*yy1**2+3._f64*xx1**2*yy1*yy2+ &
		3._f64*xx2**2*yy1*yy2+3._f64*xx1*xx2*yy1**2+3._f64*xx1*xx2*yy2**2+4._f64*xx1*xx2*yy1*yy2)*(xx2-xx1)
		w12=1._f64/60._f64*(6._f64*xx1**2*yy1**2+6._f64*xx2**2*yy2**2+xx1**2*yy2**2+xx2**2*yy1**2+3._f64*xx1**2*yy1*yy2+ &
		3._f64*xx2**2*yy1*yy2+3._f64*xx1*xx2*yy1**2+3._f64*xx1*xx2*yy2**2+4._f64*xx1*xx2*yy1*yy2)*(yy2-yy1)
		w22=-1._f64/180._f64*(10._f64*xx1**2*yy1**3+10._f64*xx2**2*yy2**3+xx2**2*yy1**3+xx1**2*yy2**3+4._f64*xx1*xx2*yy1**3 &
		+4._f64*xx1*xx2*yy2**3+6._f64*xx1**2*yy1**2*yy2+3._f64*xx1**2*yy1*yy2**2+3._f64*xx2**2*yy1**2*yy2+6._f64*xx2**2*yy1*yy2**2 &
		+6._f64*xx1*xx2*yy1**2*yy2+6._f64*xx1*xx2*yy1*yy2**2)*(xx2-xx1)

	if (cas==1) then !Lauritzen

		dax=1._f64/12._f64*(-a_moyenne(ip2,jb)+8._f64*a_moyenne(ip1,jb)-8._f64*a_moyenne(im1,jb)+a_moyenne(im2,jb))
		day=1._f64/12._f64*(-a_moyenne(ib,jp2)+8._f64*a_moyenne(ib,jp1)-8._f64*a_moyenne(ib,jm1)+a_moyenne(ib,jm2))
		bx=1._f64/4._f64*(a_moyenne(ip2,jb)-6._f64*a_moyenne(ip1,jb)+10._f64*a_moyenne(ib,jb)-6._f64*a_moyenne(im1,jb)+a_moyenne(im2,jb))
		by=1._f64/4._f64*(a_moyenne(ib,jp2)-6._f64*a_moyenne(ib,jp1)+10._f64*a_moyenne(ib,jb)-6._f64*a_moyenne(ib,jm1)+a_moyenne(ib,jm2))
		c=1._f64/4._f64*(a_moyenne(ip1,jp1)-a_moyenne(im1,jp1)-a_moyenne(ip1,jm1)+a_moyenne(im1,jm1))

		c00=a_moyenne(ib,jb)-0.5_f64*dax-0.5_f64*day-1._f64/6._f64*bx-1._f64/6._f64*by+1._f64/4._f64*c
		c10=dax+bx-0.5_f64*c
		c01=day+by-0.5_f64*c
		c20=-bx
		c02=-by
		c11=c

		res=w00*c00+w10*c10+w01*c01+w20*c20+w02*c02+w11*c11

	else if (cas==2) then !Lagrange 3

		fim1jm1=a_moyenne(im1,jm1)
		fim1j=a_moyenne(im1,jb)
		fim1jp1=a_moyenne(im1,jp1)
		fijm1=a_moyenne(ib,jm1)
		fij=a_moyenne(ib,jb)
		fijp1=a_moyenne(ib,jp1)
		fip1jm1=a_moyenne(ip1,jm1)
		fip1j=a_moyenne(ip1,jb)
		fip1jp1=a_moyenne(ip1,jp1)

		c00=1._f64/9._f64*fim1jm1+5._f64/18._f64*fim1j-1._f64/18._f64*fim1jp1+5._f64/18._f64*fijm1+25._f64/36._f64*fij-5._f64/36._f64*fijp1 &
		-1._f64/18._f64*fip1jm1-5._f64/36._f64*fip1j+1._f64/36._f64*fip1jp1
		c10=-1._f64/3._f64*fim1jm1-5._f64/6._f64*fim1j+1._f64/6._f64*fim1jp1+1._f64/3._f64*fijm1+5._f64/6._f64*fij-1._f64/6._f64*fijp1
		c20=1._f64/6._f64*fim1jm1+5._f64/12._f64*fim1j-1._f64/12._f64*fim1jp1-1._f64/3._f64*fijm1-5._f64/6._f64*fij+1._f64/6._f64*fijp1 &
		+1._f64/6._f64*fip1jm1+5._f64/12._f64*fip1j-1._f64/12._f64*fip1jp1
		c01=-1._f64/3._f64*fim1jm1+1._f64/3._f64*fim1j-5._f64/6._f64*fijm1+5._f64/6._f64*fij+1._f64/6._f64*fip1jm1-1._f64/6._f64*fip1j
		c11=fim1jm1-fim1j-fijm1+fij
		c21=-1._f64/2._f64*fim1jm1+1._f64/2._f64*fim1j+fijm1-fij-1._f64/2._f64*fip1jm1+1._f64/2._f64*fip1j
		c02=1._f64/6._f64*fim1jm1-1._f64/3._f64*fim1j+1._f64/6._f64*fim1jp1+5._f64/12._f64*fijm1-5._f64/6._f64*fij+5._f64/12._f64*fijp1 &
		-1._f64/12._f64*fip1jm1+1._f64/6._f64*fip1j-1._f64/12._f64*fip1jp1
		c12=-1._f64/2._f64*fim1jm1+fim1j-1._f64/2._f64*fim1jp1+1._f64/2._f64*fijm1-fij+1._f64/2._f64*fijp1
		c22=1._f64/4._f64*fim1jm1-1._f64/2._f64*fim1j+1._f64/4._f64*fim1jp1-1._f64/2._f64*fijm1+fij-1._f64/2._f64*fijp1+1._f64/4._f64*fip1jm1 &
		-1._f64/2._f64*fip1j+1._f64/4._f64*fip1jp1

		res=w00*c00+w10*c10+w01*c01+w20*c20+w02*c02+w11*c11+w21*c21+w12*c12+w22*c22

	else if (cas==3) then !PPM

		sij=sommets(ib,jb)
		sip1j=sommets(ip1,jb)
		sijp1=sommets(ib,jp1)
		sip1jp1=sommets(ip1,jp1)
		avij=aretesv(ib,jb)
		avip1j=aretesv(ip1,jb)
		ahij=aretesh(ib,jb)
		ahijp1=aretesh(ib,jp1)
		fij=a_moyenne(ib,jb)

		!sij=49._f64/144._f64*(a_moyenne(im1,jm1)+a_moyenne(ib,jm1)+a_moyenne(im1,jb)+a_moyenne(ib,jb))&
		!-7._f64/144._f64*(a_moyenne(im2,jm1)+a_moyenne(ip1,jm1)+a_moyenne(im2,jb)+a_moyenne(ip1,jb))&
		!-7._f64/144._f64*(a_moyenne(im1,jm2)+a_moyenne(ib,jm2)+a_moyenne(im1,jp1)+a_moyenne(ib,jp1))&
		!+1._f64/144._f64*(a_moyenne(im2,jm2)+a_moyenne(ip1,jm2)+a_moyenne(im2,jp1)+a_moyenne(ip1,jp1))
		!sip1j=49._f64/144._f64*(a_moyenne(ib,jm1)+a_moyenne(ip1,jm1)+a_moyenne(ib,jb)+a_moyenne(ip1,jb))&
		!-7._f64/144._f64*(a_moyenne(im1,jm1)+a_moyenne(ip2,jm1)+a_moyenne(im1,jb)+a_moyenne(ip2,jb))&
		!-7._f64/144._f64*(a_moyenne(ib,jm2)+a_moyenne(ip1,jm2)+a_moyenne(ib,jp1)+a_moyenne(ip1,jp1))&
		!+1._f64/144._f64*(a_moyenne(im1,jm2)+a_moyenne(ip2,jm2)+a_moyenne(im1,jp1)+a_moyenne(ip2,jp1))
		!sijp1=49._f64/144._f64*(a_moyenne(im1,jb)+a_moyenne(ib,jb)+a_moyenne(im1,jp1)+a_moyenne(ib,jp1))&
		!-7._f64/144._f64*(a_moyenne(im2,jb)+a_moyenne(ip1,jb)+a_moyenne(im2,jp1)+a_moyenne(ip1,jp1))&
		!-7._f64/144._f64*(a_moyenne(im1,jm1)+a_moyenne(ib,jm1)+a_moyenne(im1,jp2)+a_moyenne(ib,jp2))&
		!+1._f64/144._f64*(a_moyenne(im2,jm1)+a_moyenne(ip1,jm1)+a_moyenne(im2,jp2)+a_moyenne(ip1,jp2))
		!sip1jp1=49._f64/144._f64*(a_moyenne(ib,jb)+a_moyenne(ip1,jb)+a_moyenne(ib,jp1)+a_moyenne(ip1,jp1))&
		!-7._f64/144._f64*(a_moyenne(im1,jb)+a_moyenne(ip2,jb)+a_moyenne(im1,jp1)+a_moyenne(ip2,jp1))&
		!-7._f64/144._f64*(a_moyenne(ib,jm1)+a_moyenne(ip1,jm1)+a_moyenne(ib,jp2)+a_moyenne(ip1,jp2))&
		!+1._f64/144._f64*(a_moyenne(im1,jm1)+a_moyenne(ip2,jm1)+a_moyenne(im1,jp2)+a_moyenne(ip2,jp2))
		!avij=7._f64/12._f64*(a_moyenne(im1,jb)+a_moyenne(ib,jb))-1._f64/12._f64*(a_moyenne(im2,jb)+a_moyenne(ip1,jb))
		!avip1j=7._f64/12._f64*(a_moyenne(ib,jb)+a_moyenne(ip1,jb))-1._f64/12._f64*(a_moyenne(im1,jb)+a_moyenne(ip2,jb))
		!ahij=7._f64/12._f64*(a_moyenne(ib,jm1)+a_moyenne(ib,jb))-1._f64/12._f64*(a_moyenne(ib,jm2)+a_moyenne(ib,jp1))
		!ahijp1=7._f64/12._f64*(a_moyenne(ib,jb)+a_moyenne(ib,jp1))-1._f64/12._f64*(a_moyenne(ib,jm1)+a_moyenne(ib,jp2))


		c00=sij
		c10=-4._f64*sij-2._f64*sip1j+6._f64*ahij
		c20=3._f64*sij+3._f64*sip1j-6._f64*ahij
		c01=-4._f64*sij-2._f64*sijp1+6._f64*avij
		c11=16._f64*sij+8._f64*sip1j-24._f64*ahij+8._f64*sijp1+4._f64*sip1jp1-12._f64*ahijp1-24._f64*avij-12._f64*avip1j+36._f64*fij
		c21=-12._f64*sij-12._f64*sip1j+24._f64*ahij-6._f64*sijp1-6*sip1jp1+12._f64*ahijp1+18._f64*avij+18._f64*avip1j-36._f64*fij
		c02=3._f64*sij+3._f64*sijp1-6._f64*avij
		c12=-12._f64*sij-6._f64*sip1j+18._f64*ahij-12._f64*sijp1-6._f64*sip1jp1+18._f64*ahijp1+24._f64*avij+12._f64*avip1j-36._f64*fij
		c22=9._f64*sij+9._f64*sip1j-18._f64*ahij+9._f64*sijp1+9._f64*sip1jp1-18._f64*ahijp1-18._f64*avij-18*avip1j+36._f64*fij

		res=w00*c00+w10*c10+w01*c01+w20*c20+w02*c02+w11*c11+w21*c21+w12*c12+w22*c22

	else if (cas==4) then

		sij=sommetshd(ib,jb)
		sip1j=sommetshg(ip1,jb)
		sijp1=sommetsbd(ib,jp1)
		sip1jp1=sommetsbg(ip1,jp1)
		avij=aretesvd(ib,jb)
		avip1j=aretesvg(ip1,jb)
		ahij=areteshh(ib,jb)
		ahijp1=areteshb(ib,jp1)
		fij=a_moyenne(ib,jb)

		c00=sij
		c10=-4._f64*sij-2._f64*sip1j+6._f64*ahij
		c20=3._f64*sij+3._f64*sip1j-6._f64*ahij
		c01=-4._f64*sij-2._f64*sijp1+6._f64*avij
		c11=16._f64*sij+8._f64*sip1j-24._f64*ahij+8._f64*sijp1+4._f64*sip1jp1-12._f64*ahijp1-24._f64*avij-12._f64*avip1j+36._f64*fij
		c21=-12._f64*sij-12._f64*sip1j+24._f64*ahij-6._f64*sijp1-6*sip1jp1+12._f64*ahijp1+18._f64*avij+18._f64*avip1j-36._f64*fij
		c02=3._f64*sij+3._f64*sijp1-6._f64*avij
		c12=-12._f64*sij-6._f64*sip1j+18._f64*ahij-12._f64*sijp1-6._f64*sip1jp1+18._f64*ahijp1+24._f64*avij+12._f64*avip1j-36._f64*fij
		c22=9._f64*sij+9._f64*sip1j-18._f64*ahij+9._f64*sijp1+9._f64*sip1jp1-18._f64*ahijp1-18._f64*avij-18*avip1j+36._f64*fij

		res=w00*c00+w10*c10+w01*c01+w20*c20+w02*c02+w11*c11+w21*c21+w12*c12+w22*c22

	  endif

	end subroutine calcule_coeff



	subroutine calcule_coeffh(N0,N1,a_moyenne,i,j,x1,x2,res,aretesh,aretesv,sommets,areteshb,areteshh,aretesvg,aretesvd,&
                sommetsbg,sommetsbd,sommetshg,sommetshd,cas)

	        sll_int32,intent(in)::N0,N1,i,j,cas
		real(f64),intent(in)::x1,x2
		real(f64),dimension(0:N0,0:N1-1),intent(in)::a_moyenne,aretesh,aretesv,sommets,areteshb,areteshh,aretesvg,aretesvd,&
                sommetsbg,sommetsbd,sommetshg,sommetshd
		real(f64)::xx1,xx2,yy,aux,day,by,w01,w02,w21,w22,w00,w10,w20,w11,w12,c01,c02,c21,c22,c00,c10,c20,c11,c12,&
		sij,sip1j,sijp1,sip1jp1,avij,avip1j,ahij,ahijp1,&
		fij,fim1jm1,fim1j,fim1jp1,fijm1,fijp1,fip1jm1,fip1j,fip1jp1,c,yy1,yy2
		sll_int32::im2,im1,ib,ip1,ip2,jm2,jm1,jb,jp1,jp2
		real(f64),intent(out)::res
		
		if(cas==0)then
		  return
		endif

		!im2=modulo(i-2,N0)
		!im1=modulo(i-1,N0)
		!ib=modulo(i,N0)
		!ip1=modulo(i+1,N0)
		!ip2=modulo(i+2,N0)
		im2=i-2
		im1=i-1
		ib=i
		ip1=i+1
		ip2=i+2
		
		if(i-2<=0)im2=0;
		if(i-1<=0)im1=0;
		if(i<=0)ib=0;
		if(i+1<=0)ip1=0;
		if(i+2<=0)ip2=0;

		if(i-2>=N0)im2=N0;
		if(i-1>=N0)im1=N0;		
		if(i>=N0)ib=N0;
		if(i+1>=N0)ip1=N0;
		if(i+2>=N0)ip2=N0;

		
		
		jm2=modulo(j-2,N1)
		jm1=modulo(j-1,N1)
		jb=modulo(j,N1)
		jp1=modulo(j+1,N1)
		jp2=modulo(j+2,N1)

!1: Lauritzen, 2: Lag3, 3: PPM

		xx1=x1-real(i,f64)
		xx2=x2-real(i,f64)

		w01=-1._f64/2._f64*(xx2-xx1)
		w02=-1._f64/3._f64*(xx2-xx1)
		w21=-1._f64/6._f64*(xx1**2+xx2**2+xx1*xx2)*(xx2-xx1)
		w22=-1._f64/9._f64*(xx1**2+xx2**2+xx1*xx2)*(xx2-xx1)

	if (cas==1) then !Lauritzen

		day=1._f64/12._f64*(-a_moyenne(ib,jp2)+8._f64*a_moyenne(ib,jp1)-8._f64*a_moyenne(ib,jm1)+a_moyenne(ib,jm2))
		by=1._f64/4._f64*(a_moyenne(ib,jp2)-6._f64*a_moyenne(ib,jp1)+10._f64*a_moyenne(ib,jb)-6._f64*a_moyenne(ib,jm1)+a_moyenne(ib,jm2))
		c=1._f64/4._f64*(a_moyenne(ip1,jp1)-a_moyenne(im1,jp1)-a_moyenne(ip1,jm1)+a_moyenne(im1,jm1))

		c01=day+by-0.5_f64*c
		c02=-by

		res=w01*c01+w02*c02

	else if (cas==2) then !Lagrange 3

		fim1jm1=a_moyenne(im1,jm1)
		fim1j=a_moyenne(im1,jb)
		fim1jp1=a_moyenne(im1,jp1)
		fijm1=a_moyenne(ib,jm1)
		fij=a_moyenne(ib,jb)
		fijp1=a_moyenne(ib,jp1)
		fip1jm1=a_moyenne(ip1,jm1)
		fip1j=a_moyenne(ip1,jb)
		fip1jp1=a_moyenne(ip1,jp1)

		c01=-1._f64/3._f64*fim1jm1+1._f64/3._f64*fim1j-5._f64/6._f64*fijm1+5._f64/6._f64*fij+1._f64/6._f64*fip1jm1-1._f64/6._f64*fip1j
		c21=-1._f64/2._f64*fim1jm1+1._f64/2._f64*fim1j+fijm1-fij-1._f64/2._f64*fip1jm1+1._f64/2._f64*fip1j
		c02=1._f64/6._f64*fim1jm1-1._f64/3._f64*fim1j+1._f64/6._f64*fim1jp1+5._f64/12._f64*fijm1-5._f64/6._f64*fij+5._f64/12._f64*fijp1 &
		-1._f64/12._f64*fip1jm1+1._f64/6._f64*fip1j-1._f64/12._f64*fip1jp1
		c22=1._f64/4._f64*fim1jm1-1._f64/2._f64*fim1j+1._f64/4._f64*fim1jp1-1._f64/2._f64*fijm1+fij-1._f64/2._f64*fijp1+1._f64/4._f64*fip1jm1 &
		-1._f64/2._f64*fip1j+1._f64/4._f64*fip1jp1

		res=w01*c01+w02*c02+w21*c21+w22*c22

	else if (cas==3) then !PPM

		sij=sommets(ib,jb)
		sip1j=sommets(ip1,jb)
		sijp1=sommets(ib,jp1)
		sip1jp1=sommets(ip1,jp1)
		avij=aretesv(ib,jb)
		avip1j=aretesv(ip1,jb)
		ahij=aretesh(ib,jb)
		ahijp1=aretesh(ib,jp1)
		fij=a_moyenne(ib,jb)

		!sij=49._f64/144._f64*(a_moyenne(im1,jm1)+a_moyenne(ib,jm1)+a_moyenne(im1,jb)+a_moyenne(ib,jb))&
		!-7._f64/144._f64*(a_moyenne(im2,jm1)+a_moyenne(ip1,jm1)+a_moyenne(im2,jb)+a_moyenne(ip1,jb))&
		!-7._f64/144._f64*(a_moyenne(im1,jm2)+a_moyenne(ib,jm2)+a_moyenne(im1,jp1)+a_moyenne(ib,jp1))&
		!+1._f64/144._f64*(a_moyenne(im2,jm2)+a_moyenne(ip1,jm2)+a_moyenne(im2,jp1)+a_moyenne(ip1,jp1))
		!sip1j=49._f64/144._f64*(a_moyenne(ib,jm1)+a_moyenne(ip1,jm1)+a_moyenne(ib,jb)+a_moyenne(ip1,jb))&
		!-7._f64/144._f64*(a_moyenne(im1,jm1)+a_moyenne(ip2,jm1)+a_moyenne(im1,jb)+a_moyenne(ip2,jb))&
		!-7._f64/144._f64*(a_moyenne(ib,jm2)+a_moyenne(ip1,jm2)+a_moyenne(ib,jp1)+a_moyenne(ip1,jp1))&
		!+1._f64/144._f64*(a_moyenne(im1,jm2)+a_moyenne(ip2,jm2)+a_moyenne(im1,jp1)+a_moyenne(ip2,jp1))
		!sijp1=49._f64/144._f64*(a_moyenne(im1,jb)+a_moyenne(ib,jb)+a_moyenne(im1,jp1)+a_moyenne(ib,jp1))&
		!-7._f64/144._f64*(a_moyenne(im2,jb)+a_moyenne(ip1,jb)+a_moyenne(im2,jp1)+a_moyenne(ip1,jp1))&
		!-7._f64/144._f64*(a_moyenne(im1,jm1)+a_moyenne(ib,jm1)+a_moyenne(im1,jp2)+a_moyenne(ib,jp2))&
		!+1._f64/144._f64*(a_moyenne(im2,jm1)+a_moyenne(ip1,jm1)+a_moyenne(im2,jp2)+a_moyenne(ip1,jp2))
		!sip1jp1=49._f64/144._f64*(a_moyenne(ib,jb)+a_moyenne(ip1,jb)+a_moyenne(ib,jp1)+a_moyenne(ip1,jp1))&
		!-7._f64/144._f64*(a_moyenne(im1,jb)+a_moyenne(ip2,jb)+a_moyenne(im1,jp1)+a_moyenne(ip2,jp1))&
		!-7._f64/144._f64*(a_moyenne(ib,jm1)+a_moyenne(ip1,jm1)+a_moyenne(ib,jp2)+a_moyenne(ip1,jp2))&
		!+1._f64/144._f64*(a_moyenne(im1,jm1)+a_moyenne(ip2,jm1)+a_moyenne(im1,jp2)+a_moyenne(ip2,jp2))
		!avij=7._f64/12._f64*(a_moyenne(im1,jb)+a_moyenne(ib,jb))-1._f64/12._f64*(a_moyenne(im2,jb)+a_moyenne(ip1,jb))
		!avip1j=7._f64/12._f64*(a_moyenne(ib,jb)+a_moyenne(ip1,jb))-1._f64/12._f64*(a_moyenne(im1,jb)+a_moyenne(ip2,jb))
		!ahij=7._f64/12._f64*(a_moyenne(ib,jm1)+a_moyenne(ib,jb))-1._f64/12._f64*(a_moyenne(ib,jm2)+a_moyenne(ib,jp1))
		!ahijp1=7._f64/12._f64*(a_moyenne(ib,jb)+a_moyenne(ib,jp1))-1._f64/12._f64*(a_moyenne(ib,jm1)+a_moyenne(ib,jp2))

		c01=-4._f64*sij-2._f64*sijp1+6._f64*avij
		c21=-12._f64*sij-12._f64*sip1j+24._f64*ahij-6._f64*sijp1-6*sip1jp1+12._f64*ahijp1+18._f64*avij+18._f64*avip1j-36._f64*fij
		c02=3._f64*sij+3._f64*sijp1-6._f64*avij
		c22=9._f64*sij+9._f64*sip1j-18._f64*ahij+9._f64*sijp1+9._f64*sip1jp1-18._f64*ahijp1-18._f64*avij-18*avip1j+36._f64*fij

		res=w01*c01+w02*c02+w21*c21+w22*c22

	else if (cas==4) then

		sij=sommetshd(ib,jb)
		sip1j=sommetshg(ip1,jb)
		sijp1=sommetsbd(ib,jp1)
		sip1jp1=sommetsbg(ip1,jp1)
		avij=aretesvd(ib,jb)
		avip1j=aretesvg(ip1,jb)
		ahij=areteshh(ib,jb)
		ahijp1=areteshb(ib,jp1)
		fij=a_moyenne(ib,jb)

		c01=-4._f64*sij-2._f64*sijp1+6._f64*avij
		c21=-12._f64*sij-12._f64*sip1j+24._f64*ahij-6._f64*sijp1-6*sip1jp1+12._f64*ahijp1+18._f64*avij+18._f64*avip1j-36._f64*fij
		c02=3._f64*sij+3._f64*sijp1-6._f64*avij
		c22=9._f64*sij+9._f64*sip1j-18._f64*ahij+9._f64*sijp1+9._f64*sip1jp1-18._f64*ahijp1-18._f64*avij-18*avip1j+36._f64*fij

		res=w01*c01+w02*c02+w21*c21+w22*c22

	  endif

	end subroutine calcule_coeffh



	subroutine calcule_coeffv(N0,N1,a_moyenne,i,j,y1,y2,res,aretesh,aretesv,sommets,areteshb,areteshh,aretesvg,aretesvd,&
                sommetsbg,sommetsbd,sommetshg,sommetshd,cas)

	        sll_int32,intent(in)::N0,N1,i,j,cas
		real(f64),intent(in)::y1,y2
		real(f64),dimension(0:N0,0:N1-1),intent(in)::a_moyenne,aretesh,aretesv,sommets,areteshb,areteshh,aretesvg,&
                aretesvd,sommetsbg,sommetsbd,sommetshg,sommetshd
		real(f64)::xx,yy1,yy2,aux,dax,bx,by,c,w00,w10,w20,w11,w12,w01,w02,w21,w22,&
		c00,c10,c20,c11,c12,c01,c02,c21,c22,sij,sip1j,sijp1,sip1jp1,avij,avip1j,ahij,ahijp1,&
		fij,fim1jm1,fim1j,fim1jp1,fijm1,fijp1,fip1jm1,fip1j,fip1jp1,day,xx1,xx2
		sll_int32::im2,im1,ib,ip1,ip2,jm2,jm1,jb,jp1,jp2
		real(f64),intent(out)::res
		
		if(cas==0)then
		  return
		endif

		!im2=modulo(i-2,N0)
		!im1=modulo(i-1,N0)
		!ib=modulo(i,N0)
		!ip1=modulo(i+1,N0)
		!ip2=modulo(i+2,N0)
		im2=i-2
		im1=i-1
		ib=i
		ip1=i+1
		ip2=i+2
		
		if(i-2<=0)im2=0;
		if(i-1<=0)im1=0;
		if(i<=0)ib=0;
		if(i+1<=0)ip1=0;
		if(i+2<=0)ip2=0;

		if(i-2>=N0)im2=N0;
		if(i-1>=N0)im1=N0;		
		if(i>=N0)ib=N0;
		if(i+1>=N0)ip1=N0;
		if(i+2>=N0)ip2=N0;

		
		
		
		jm2=modulo(j-2,N1)
		jm1=modulo(j-1,N1)
		jb=modulo(j,N1)
		jp1=modulo(j+1,N1)
		jp2=modulo(j+2,N1)

!1: Lauritzen, 2: Lag3, 3: PPM
		yy1=y1-real(j,f64)
		yy2=y2-real(j,f64)

		w00=yy2-yy1
		w10=1._f64/2._f64*(yy2-yy1)
		w20=1._f64/3._f64*(yy2-yy1)
		w11=1._f64/4._f64*(yy1+yy2)*(yy2-yy1)
		w12=1._f64/6._f64*(yy1**2+yy2**2+yy1*yy2)*(yy2-yy1)

	if (cas==1) then !Lauritzen

		dax=1._f64/12._f64*(-a_moyenne(ip2,jb)+8._f64*a_moyenne(ip1,jb)-8._f64*a_moyenne(im1,jb)+a_moyenne(im2,jb))
		day=1._f64/12._f64*(-a_moyenne(ib,jp2)+8._f64*a_moyenne(ib,jp1)-8._f64*a_moyenne(ib,jm1)+a_moyenne(ib,jm2))
		bx=1._f64/4._f64*(a_moyenne(ip2,jb)-6._f64*a_moyenne(ip1,jb)+10._f64*a_moyenne(ib,jb)-6._f64*a_moyenne(im1,jb)+a_moyenne(im2,jb))
		by=1._f64/4._f64*(a_moyenne(ib,jp2)-6._f64*a_moyenne(ib,jp1)+10._f64*a_moyenne(ib,jb)-6._f64*a_moyenne(ib,jm1)+a_moyenne(ib,jm2))
		c=1._f64/4._f64*(a_moyenne(ip1,jp1)-a_moyenne(im1,jp1)-a_moyenne(ip1,jm1)+a_moyenne(im1,jm1))

		c00=a_moyenne(ib,jb)-0.5_f64*dax-0.5_f64*day-1._f64/6._f64*bx-1._f64/6._f64*by+1._f64/4._f64*c
		c10=dax+bx-0.5_f64*c
		c20=-bx
		c11=c

		res=w00*c00+w10*c10+w20*c20+w11*c11

	else if (cas==2) then !Lagrange 3

		fim1jm1=a_moyenne(im1,jm1)
		fim1j=a_moyenne(im1,jb)
		fim1jp1=a_moyenne(im1,jp1)
		fijm1=a_moyenne(ib,jm1)
		fij=a_moyenne(ib,jb)
		fijp1=a_moyenne(ib,jp1)
		fip1jm1=a_moyenne(ip1,jm1)
		fip1j=a_moyenne(ip1,jb)
		fip1jp1=a_moyenne(ip1,jp1)

		c00=1._f64/9._f64*fim1jm1+5._f64/18._f64*fim1j-1._f64/18._f64*fim1jp1+5._f64/18._f64*fijm1+25._f64/36._f64*fij-5._f64/36._f64*fijp1 &
		-1._f64/18._f64*fip1jm1-5._f64/36._f64*fip1j+1._f64/36._f64*fip1jp1
		c10=-1._f64/3._f64*fim1jm1-5._f64/6._f64*fim1j+1._f64/6._f64*fim1jp1+1._f64/3._f64*fijm1+5._f64/6._f64*fij-1._f64/6._f64*fijp1
		c20=1._f64/6._f64*fim1jm1+5._f64/12._f64*fim1j-1._f64/12._f64*fim1jp1-1._f64/3._f64*fijm1-5._f64/6._f64*fij+1._f64/6._f64*fijp1 &
		+1._f64/6._f64*fip1jm1+5._f64/12._f64*fip1j-1._f64/12._f64*fip1jp1
		c11=fim1jm1-fim1j-fijm1+fij
		c12=-1._f64/2._f64*fim1jm1+fim1j-1._f64/2._f64*fim1jp1+1._f64/2._f64*fijm1-fij+1._f64/2._f64*fijp1

		res=w00*c00+w10*c10+w20*c20+w11*c11+w12*c12

	else if (cas==3) then !PPM

		sij=sommets(ib,jb)
		sip1j=sommets(ip1,jb)
		sijp1=sommets(ib,jp1)
		sip1jp1=sommets(ip1,jp1)
		avij=aretesv(ib,jb)
		avip1j=aretesv(ip1,jb)
		ahij=aretesh(ib,jb)
		ahijp1=aretesh(ib,jp1)
		fij=a_moyenne(ib,jb)

		!sij=49._f64/144._f64*(a_moyenne(im1,jm1)+a_moyenne(ib,jm1)+a_moyenne(im1,jb)+a_moyenne(ib,jb))&
		!-7._f64/144._f64*(a_moyenne(im2,jm1)+a_moyenne(ip1,jm1)+a_moyenne(im2,jb)+a_moyenne(ip1,jb))&
		!-7._f64/144._f64*(a_moyenne(im1,jm2)+a_moyenne(ib,jm2)+a_moyenne(im1,jp1)+a_moyenne(ib,jp1))&
		!+1._f64/144._f64*(a_moyenne(im2,jm2)+a_moyenne(ip1,jm2)+a_moyenne(im2,jp1)+a_moyenne(ip1,jp1))
		!sip1j=49._f64/144._f64*(a_moyenne(ib,jm1)+a_moyenne(ip1,jm1)+a_moyenne(ib,jb)+a_moyenne(ip1,jb))&
		!-7._f64/144._f64*(a_moyenne(im1,jm1)+a_moyenne(ip2,jm1)+a_moyenne(im1,jb)+a_moyenne(ip2,jb))&
		!-7._f64/144._f64*(a_moyenne(ib,jm2)+a_moyenne(ip1,jm2)+a_moyenne(ib,jp1)+a_moyenne(ip1,jp1))&
		!+1._f64/144._f64*(a_moyenne(im1,jm2)+a_moyenne(ip2,jm2)+a_moyenne(im1,jp1)+a_moyenne(ip2,jp1))
		!sijp1=49._f64/144._f64*(a_moyenne(im1,jb)+a_moyenne(ib,jb)+a_moyenne(im1,jp1)+a_moyenne(ib,jp1))&
		!-7._f64/144._f64*(a_moyenne(im2,jb)+a_moyenne(ip1,jb)+a_moyenne(im2,jp1)+a_moyenne(ip1,jp1))&
		!-7._f64/144._f64*(a_moyenne(im1,jm1)+a_moyenne(ib,jm1)+a_moyenne(im1,jp2)+a_moyenne(ib,jp2))&
		!+1._f64/144._f64*(a_moyenne(im2,jm1)+a_moyenne(ip1,jm1)+a_moyenne(im2,jp2)+a_moyenne(ip1,jp2))
		!sip1jp1=49._f64/144._f64*(a_moyenne(ib,jb)+a_moyenne(ip1,jb)+a_moyenne(ib,jp1)+a_moyenne(ip1,jp1))&
		!-7._f64/144._f64*(a_moyenne(im1,jb)+a_moyenne(ip2,jb)+a_moyenne(im1,jp1)+a_moyenne(ip2,jp1))&
		!-7._f64/144._f64*(a_moyenne(ib,jm1)+a_moyenne(ip1,jm1)+a_moyenne(ib,jp2)+a_moyenne(ip1,jp2))&
		!+1._f64/144._f64*(a_moyenne(im1,jm1)+a_moyenne(ip2,jm1)+a_moyenne(im1,jp2)+a_moyenne(ip2,jp2))
		!avij=7._f64/12._f64*(a_moyenne(im1,jb)+a_moyenne(ib,jb))-1._f64/12._f64*(a_moyenne(im2,jb)+a_moyenne(ip1,jb))
		!avip1j=7._f64/12._f64*(a_moyenne(ib,jb)+a_moyenne(ip1,jb))-1._f64/12._f64*(a_moyenne(im1,jb)+a_moyenne(ip2,jb))
		!ahij=7._f64/12._f64*(a_moyenne(ib,jm1)+a_moyenne(ib,jb))-1._f64/12._f64*(a_moyenne(ib,jm2)+a_moyenne(ib,jp1))
		!ahijp1=7._f64/12._f64*(a_moyenne(ib,jb)+a_moyenne(ib,jp1))-1._f64/12._f64*(a_moyenne(ib,jm1)+a_moyenne(ib,jp2))

		c00=sij
		c10=-4._f64*sij-2._f64*sip1j+6._f64*ahij
		c20=3._f64*sij+3._f64*sip1j-6._f64*ahij
		c11=16._f64*sij+8._f64*sip1j-24._f64*ahij+8._f64*sijp1+4._f64*sip1jp1-12._f64*ahijp1-24._f64*avij-12._f64*avip1j+36._f64*fij
		c12=-12._f64*sij-6._f64*sip1j+18._f64*ahij-12._f64*sijp1-6._f64*sip1jp1+18._f64*ahijp1+24._f64*avij+12._f64*avip1j-36._f64*fij

		res=w00*c00+w10*c10+w20*c20+w11*c11+w12*c12

	else if (cas==4) then

		sij=sommetshd(ib,jb)
		sip1j=sommetshg(ip1,jb)
		sijp1=sommetsbd(ib,jp1)
		sip1jp1=sommetsbg(ip1,jp1)
		avij=aretesvd(ib,jb)
		avip1j=aretesvg(ip1,jb)
		ahij=areteshh(ib,jb)
		ahijp1=areteshb(ib,jp1)
		fij=a_moyenne(ib,jb)
		c00=sij

		c10=-4._f64*sij-2._f64*sip1j+6._f64*ahij
		c20=3._f64*sij+3._f64*sip1j-6._f64*ahij
		c11=16._f64*sij+8._f64*sip1j-24._f64*ahij+8._f64*sijp1+4._f64*sip1jp1-12._f64*ahijp1-24._f64*avij-12._f64*avip1j+36._f64*fij
		c12=-12._f64*sij-6._f64*sip1j+18._f64*ahij-12._f64*sijp1-6._f64*sip1jp1+18._f64*ahijp1+24._f64*avij+12._f64*avip1j-36._f64*fij

		res=w00*c00+w10*c10+w20*c20+w11*c11+w12*c12

	  endif

	end subroutine calcule_coeffv



	subroutine init_random_seed()
		implicit none
		sll_int32 :: i, n, heure
		sll_int32, DIMENSION(:), ALLOCATABLE :: graine
        
		CALL RANDOM_SEED(size = n)
		ALLOCATE(graine(n))
		CALL SYSTEM_CLOCK(COUNT=heure)
		graine = heure + 37 * (/ (i - 1, i = 1, n) /)
!print*,heure
		CALL RANDOM_SEED(PUT = graine)
		DEALLOCATE(graine)

	end subroutine



end module polar_advection
