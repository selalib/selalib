module polar_advection
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"

  use polar_operators
  use poisson_polar
  use numeric_constants
  use sll_splines
  implicit none

  !>type sll_plan_adv_polar
  !>type for advection with center-guide equations
  !>the field and other needed data/object are within
  type sll_plan_adv_polar
     sll_real64 :: rmin,rmax,dr,dtheta,dt
     sll_int32 :: nr,ntheta
     type(sll_spline_2D), pointer :: spl_f
     sll_int32 :: time_scheme
     sll_real64, dimension(:,:,:), pointer :: field
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
         & HERMITE_SPLINE, PERIODIC_SPLINE,const_slope_x1_min = 0._f64,const_slope_x1_max = 0._f64)

   ! this%spl_f => new_spline_2D(nr+1,ntheta+1,rmin,rmax,0._f64, 2._f64*sll_pi, &
   !      & PERIODIC_SPLINE, PERIODIC_SPLINE)

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
  !>grad_case : integer, see function new_polar_op
  !>time_scheme : integer, see function new_plan_adv_polar
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
       maxiter=10
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


  !> subroutine correction_r(r,rmin,rmax)
  !> correction of r to stay in the domain
  !> r : value tu corret
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

    dt=plan%adv%dt
    plan%adv%dt=dt/2.0_f64

    call poisson_solve_polar(plan%poisson,in,plan%phi)
    call compute_grad_field(plan%grad,plan%phi,plan%adv%field)
    call advect_CG_polar(plan%adv,in,out)
    !we just obtained f^(n+1/2)
    call poisson_solve_polar(plan%poisson,out,plan%phi)
    call compute_grad_field(plan%grad,plan%phi,plan%adv%field)
    !we just obtained E^(n+1/2)
    plan%adv%dt=dt
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





end module polar_advection
