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
     type(sll_polar_data), pointer :: data
     type(sll_spline_2D), pointer :: spl_f
     sll_int32 :: interpolate_case
     sll_real64, dimension(:,:,:), allocatable :: field
  end type sll_plan_adv_polar

  !>type sll_SL_polar
  !>type for semi Lagrangian
  !>contains other types for the routines called in SL routines
  type sll_SL_polar
     type(sll_plan_adv_polar), pointer :: adv
     type(plan_polar_op), pointer :: grad
     type(sll_plan_poisson_polar), pointer :: poisson
     sll_real64, dimension(:,:), allocatable :: phi
  end type sll_SL_polar

contains

!===================================
!  creation of sll_plan_adv_polar
!===================================

  !>function new_plan_adv_polar(data,interpolate_case)
  !>data : sll_polar_data object, contains data about the domain
  !>interpolate_case : integer to choose the scheme for advection
  !>                   1 : using explicit Euler method
  !>                   2 : rotation, rotation speed = -1
  !>                   3 : using symplectic Euler with linear interpolation
  !>                   4 : using symplectic Verlet with linear interpolation
  !>                   5 : using fixed point method
  function new_plan_adv_polar(data,interpolate_case) result(this)

    type(sll_plan_adv_polar), pointer :: this
    type(sll_polar_data), intent(in), pointer :: data
    sll_int32, intent(in) :: interpolate_case

    sll_int32 :: err

    SLL_ALLOCATE(this,err)
    SLL_ALLOCATE(this%data,err)
    SLL_ALLOCATE(this%field(2,data%nr+1,data%ntheta+1),err)

    this%field=0.0_f64
    this%data=data
    this%interpolate_case=interpolate_case

    this%spl_f => new_spline_2D(data%nr+1,data%ntheta+1,data%rmin,data%rmax,0._f64, 2._f64*sll_pi, &
         & HERMITE_SPLINE, PERIODIC_SPLINE,const_slope_x1_min = 0._f64,const_slope_x1_max = 0._f64)

  end function new_plan_adv_polar

!===================================
!  deletion of sll_plan_adv_polar
!===================================

  !>delete_plan_adv_polar(this)
  !>deletion of sll_plan_adv_polar object
  subroutine delete_plan_adv_polar(this)

    implicit none

    type(sll_plan_adv_polar), intent(inout), pointer :: this

    sll_int32 :: err

    if (associated(this)) then
       SLL_DEALLOCATE(this%data,err)
       SLL_DEALLOCATE_ARRAY(this%field,err)
       call delete_spline_2d(this%spl_f)
       this%data=>null()
       SLL_DEALLOCATE(this,err)
       this=>null()
    end if

  end subroutine delete_plan_adv_polar

!==================================
!  creation of sll_SL_polar type
!==================================

  !>function new_SL(data,grad_case,interpolate_case)
  !>creation of sll_SL_polar object for semi Lagrangian scheme in polar coordinates
  !>data : sll_polar_data
  !>grad_case : integer, see function new_polar_op
  !>interpolate_case : integer, see function new_plan_adv_polar
  function new_SL(data,grad_case,interpolate_case) result(this)

    type(sll_SL_polar), pointer :: this
    type(sll_polar_data), intent(in), pointer :: data
    sll_int32, intent(in) :: grad_case,interpolate_case

    sll_int32 :: err

    SLL_ALLOCATE(this,err)
    SLL_ALLOCATE(this%phi(data%nr+1,data%ntheta+1),err)

    this%poisson => new_plan_poisson_polar(data)
    this%grad => new_polar_op(data,grad_case)
    this%adv => new_plan_adv_polar(data,interpolate_case)

  end function new_SL

!=============================
!  deletion of sll_SL_polar
!=============================

  !>subroutine delete_SL_polar(this)
  !>deletion of sll_SL_polar object
  subroutine delete_SL_polar(this)

    implicit none

    type(sll_SL_polar), intent(inout), pointer :: this

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

  !>subroutine advect_CG_polar(plan,in,out)
  !>compute step for Center-Guide equation
  !>plan : sll_plan_adv_polar object
  !>in : distribution function at time t_n, size (nr+1)*(ntheta+1)
  !>out : distribution function at time t_(n+1), size (nr+1)*(ntheta+1)
  subroutine advect_CG_polar(plan,fn,fnp1)

    implicit none

    type(sll_plan_adv_polar), intent(inout), pointer :: plan
    sll_real64, dimension(:,:), intent(in) :: fn
    sll_real64, dimension(:,:), intent(out) :: fnp1

    sll_int32 :: nr, ntheta
    sll_real64 :: dt, dr, dtheta, rmin, rmax
    sll_int32 :: i,j,maxiter,iter,kr,k
    sll_real64 :: r,theta,rr,rrn,ttheta,tthetan,tolr,tolth,ar,atheta

    nr=plan%data%nr
    ntheta=plan%data%ntheta
    dt=plan%data%dt
    dr=plan%data%dr
    dtheta=plan%data%dtheta
    rmin=plan%data%rmin
    rmax=plan%data%rmax

    !construction of spline coefficients for f
    call compute_spline_2D(fn,plan%spl_f)

    if (plan%interpolate_case==1) then
       !explicit Euler
       do i=1,nr+1
          r=rmin+real(i-1,f64)*dr
          do j=1,ntheta+1
             theta=real(j-1,f64)*dtheta

             theta=theta-dt*plan%field(2,i,j)
             r=r-dt*plan%field(1,i,j)

             call correction_r(r,rmin,rmax)
             call correction_theta(theta)
             fnp1(i,j)=interpolate_value_2D(r,theta,plan%spl_f)

          end do
       end do

    else if (plan%interpolate_case==2) then
       !rotation

       do i=1,nr+1
          do j=1,ntheta+1
             theta=real(j-1,f64)*dtheta-dt
             call correction_theta(theta)
             fnp1(i,j)=interpolate_value_2D(r,theta,plan%spl_f)
          end do
       end do

    else if (plan%interpolate_case==3) then
       !using symplectic Euler with linear interpolation
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
                   rr=rmin+real(i-1,f64)*dr-dt*plan%field(1,i,j)
                else if (k<nr+1 .and. k>=1) then
                   rr=rmin+real(i-1,f64)*dr-dt*((1.0_f64-r)*plan%field(1,k,j)+r*plan%field(1,k+1,j))
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
                theta=real(j-1,f64)*dtheta-dt*((1.0_f64-r)*plan%field(2,k,j)+r*plan%field(2,k+1,j))
             else
                theta=real(j-1,f64)*dtheta-dt*plan%field(2,k,j)
             end if
             call correction_theta(theta)

             fnp1(i,j)=interpolate_value_2d(rr,theta,plan%spl_f)
          end do
       end do

    else if (plan%interpolate_case==4) then
       !using symplectic Verlet with linear interpolation

       !we fix the tolerance and the maximum of iteration
       tolr=1e-12
       tolth=1e-12
       tolr=1e-4
       tolth=1e-4
       maxiter=1000

       do j=1,ntheta
          do i=1,nr+1
             !initialization for r interpolation
             rr=rmin+real(i-1,f64)*dr+dt/2.0_f64*plan%field(1,i,j)
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
                   rr=rmin+real(i-1,f64)*dr-0.5_f64*dt*plan%field(1,kr,j)
                else if (kr>0 .and. kr<nr+1) then
                   rr=rmin+real(i-1,f64)*dr-0.5_f64*dt*((1.0_f64-r)*plan%field(1,kr,j)+r*plan%field(1,kr+1,j))
                else
                   print*,kr
                   print*,'error : kr is not in range'
                   print*,'exiting'
                   stop
                end if
                call correction_r(rr,rmin,rmax)

                iter=iter+1
             end do
!!$             if (iter==maxiter .and. abs(rrn-rr)>tolr) then
!!$                print*,'not enought iterations for r in symplectic Verlet',i,j,kr,rr,rrn
!!$                stop
!!$             end if
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
                   ttheta=real(j-1,f64)*dtheta-0.5_f64*dt*((1.0_f64-theta)*plan%field(2,kr,k)+theta*plan%field(2,kr,k+1))
                   ttheta=ttheta-0.5_f64*dt*plan%field(2,kr,j)
                else
                   ttheta=real(j-1,f64)*dtheta-0.5_f64*dt*((1.0_f64-theta)*((1.0_f64-r)*plan%field(2,kr,k)+r*plan%field(2,kr+1,k)) &
                        & +theta*((1.0_f64-r)*plan%field(1,kr,k+1)+r*plan%field(2,kr+1,k+1)))
                   ttheta=ttheta-0.5_f64*dt*((1.0_f64-r)*plan%field(2,kr,j)+r*plan%field(2,kr+1,j))
                end if
                call correction_theta(ttheta)

                iter=iter+1
             end do
!!$             if (iter==maxiter .and. abs(tthetan-ttheta)>tolth .and. abs(tthetan+2.0_f64*sll_pi-ttheta)>tolth &
!!$                  & .and.abs(tthetan-ttheta-2.0_f64*sll_pi)>tolth) then
!!$                print*,'not enought iterations for theta in symplectic Verlet',i,j,k,ttheta,tthetan
!!$                stop
!!$             end if
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
                rr=rr+0.5_f64*dt*((1.0_f64-theta)*plan%field(1,kr,k)+theta*plan%field(1,kr,k+1))
             else
                rr=rr+0.5_f64*dt*((1.0_f64-theta)*((1.0_f64-r)*plan%field(1,kr,k)+r*plan%field(1,kr+1,k)) &
                     & +theta*((1.0_f64-r)*plan%field(1,kr,k+1)+r*plan%field(1,kr+1,k+1)))
             end if
             call correction_r(rr,rmin,rmax)

             fnp1(i,j)=interpolate_value_2d(rr,ttheta,plan%spl_f)

          end do
       end do

    else if (plan%interpolate_case==5) then
       !using fixed point method

       !initialization
       maxiter=10
       tolr=(dr+dtheta)/5.0_f64

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
                   ar=0.5_f64*dt*((1.0_f64-theta)*((1.0_f64-r)*plan%field(1,kr,k)+theta*((1.0_f64-r)*plan%field(1,kr,k+1))))
                   atheta=0.5_f64*dt*((1.0_f64-theta)*((1.0_f64-r)*plan%field(2,kr,k)+theta*((1.0_f64-r)*plan%field(2,kr,k+1))))
                else
                   ar=0.5_f64*dt*((1.0_f64-theta)*((1.0_f64-r)*plan%field(1,kr,k)+r*plan%field(1,kr+1,k)) &
                        & +theta*((1.0_f64-r)*plan%field(1,kr,k+1)+r*plan%field(1,kr+1,k+1)))
                   atheta=0.5_f64*dt*((1.0_f64-theta)*((1.0_f64-r)*plan%field(2,kr,k)+r*plan%field(2,kr+1,k)) &
                        & +theta*((1.0_f64-r)*plan%field(2,kr,k+1)+r*plan%field(1,kr+1,k+1)))
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
  subroutine correction_theta(theta)

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

  end subroutine correction_theta


  !>subroutine SL_classic(plan,in,out)
  !>computes the classic semi-Lagrangian scheme for Vlasov-Poisson equation
  !>plan : sll_SL_polar object, contains plan for Poisso, gradian and advection
  !>in : distribution function at time n, size (nr+1)*(ntheta+1)
  !>out : distribution function at time n+1, size (nr+1)*(ntheta+1)
  subroutine SL_classic(plan,in,out)

    implicit none

    type(sll_SL_polar), intent(inout), pointer :: plan
    sll_real64, dimension(:,:), intent(inout) :: in,out

    sll_int32 :: i,j
    sll_real64 :: temp,r

    call poisson_solve_polar(plan%poisson,in,plan%phi)
    call compute_grad_field(plan%grad,plan%phi,plan%adv%field)
    do i=1,plan%adv%data%nr+1
       r=plan%adv%data%rmin+plan%adv%data%dr*real(i-1,f64)
       do j=1,plan%adv%data%ntheta+1
          temp=plan%adv%field(1,i,j)/r
          plan%adv%field(1,i,j)=-plan%adv%field(2,i,j)
          plan%adv%field(2,i,j)=temp
       end do
    end do
    call advect_CG_polar(plan%adv,in,out)

  end subroutine SL_classic


  !>subroutine SL_ordre_2(plan,in,out)
  !>computes the semi-Lagrangian scheme order 2
  !>plan : sll_SL_polar object, contains plan for Poisso, gradian and advection
  !>in : distribution function at time n, size (nr+1)*(ntheta+1)
  !>out : distribution function at time n+1, size (nr+1)*(ntheta+1)
  subroutine SL_ordre_2(plan,in,out)

    implicit none

    type(sll_SL_polar), intent(inout), pointer :: plan
    sll_real64, dimension(:,:), intent(inout) :: in,out

    sll_int32 :: i,j
    sll_real64 :: dt,temp,r

    dt=plan%adv%data%dt
    plan%adv%data%dt=dt/2.0_f64

    call poisson_solve_polar(plan%poisson,in,plan%phi)
    call compute_grad_field(plan%grad,plan%phi,plan%adv%field)
    do i=1,plan%adv%data%nr+1
       r=plan%adv%data%rmin+plan%adv%data%dr*real(i-1,f64)
       do j=1,plan%adv%data%ntheta+1
          temp=plan%adv%field(1,i,j)/r
          plan%adv%field(1,i,j)=-plan%adv%field(2,i,j)
          plan%adv%field(2,i,j)=temp
       end do
    end do
    call advect_CG_polar(plan%adv,in,out)
    !we just obtained f^(n+1/2)
    call poisson_solve_polar(plan%poisson,out,plan%phi)
    call compute_grad_field(plan%grad,plan%phi,plan%adv%field)
    do i=1,plan%adv%data%nr+1
       r=plan%adv%data%rmin+plan%adv%data%dr*real(i-1,f64)
       do j=1,plan%adv%data%ntheta+1
          temp=plan%adv%field(1,i,j)/r
          plan%adv%field(1,i,j)=-plan%adv%field(2,i,j)
          plan%adv%field(2,i,j)=temp
       end do
    end do
    !we just obtained E^(n+1/2)
    plan%adv%data%dt=dt
    call advect_CG_polar(plan%adv,in,out)

  end subroutine SL_ordre_2

end module polar_advection
