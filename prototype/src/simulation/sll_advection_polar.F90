module polar_advection
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"

  use polar_kind
  use poisson_polar
  use numeric_constants
  use sll_splines
  implicit none

contains

  !>subroutine compute_grad_field(data,phi,grad_phi,spl_phi)
  !>compute grad(phi) for phi scalar field in polar coordinate
  !>data : polar_data, contains data about the mesh
  !>phi : scalar field, size (nr+1)*(ntheta+1)
  !>grad_phi : grad(phi), size 2*(nr+1)*(ntheta+1)
  !>           grad_phi(1,:,:)=d_r(phi), grad_phi(2,:,:)=d_theta(phi)/r
  !>spl_phi : sll_spline_2D object, spline of pji, computed if needed
  subroutine compute_grad_field(data,phi,grad_phi,spl_phi)

    implicit none

    type(polar_data), intent(inout), pointer :: data
    sll_real64, dimension(:,:), intent(in) :: phi
    sll_real64, dimension(:,:,:), intent(out) :: grad_phi
    type(sll_spline_2D), pointer :: spl_phi

    sll_int32 :: nr, ntheta
    sll_real64 :: dr, dtheta, rmin, rmax
    sll_int32 :: i,j,calculus
    sll_real64 :: r,theta

    nr=data%nr
    ntheta=data%ntheta
    dr=data%dr
    dtheta=data%dtheta
    rmin=data%rmin
    rmax=data%rmax

    ! way to compute a
    ! 1 : center for r and theta
    !     decenter on boundaries
    ! 2 : center for r, decenter on boundaries
    !     using fft for theta
    ! 3 : using splines
    calculus = 3

    if (calculus==1) then
       ! center formula for r end theta
       ! decenter for r on boundaries
       do i=2,nr
          r=rmin+real(i-1,f64)*dr
          do j=1,ntheta
             grad_phi(1,i,j)=(phi(i+1,j)-phi(i-1,j))/(2*dr)
             grad_phi(2,i,j)=(phi(i,modulo(j+1-1+ntheta,ntheta)+1)-phi(i,modulo(j-1-1+ntheta,ntheta)+1))/(2*r*dtheta)
          end do
       end do
       do j=1,ntheta
!!$          grad_phi(1,1,j)=phi(2,j)-0.0_f64)/(2.0_f64*dr)
!!$          grad_phi(1,nr+1,j)=(0.0_f64-phi(nr,j))/(dr*2.0_f64)
!!$          grad_phi(1,1,j)=0.0_f64
!!$          grad_phi(1,nr+1,j)=0.0_f64
          grad_phi(1,1,j)=(phi(2,j)-phi(1,j))/dr
          grad_phi(1,nr+1,j)=(phi(nr,j)-phi(nr+1,j))/dr
!!$          grad_phi(1,1,j)=-(1.5_f64*phi(1,j)-2.0_f64*phi(2,j)+0.5_f64*phi(3,j))/dr
!!$          grad_phi(1,nr+1,j)=-(1.5_f64*phi(nr+1,j)-2.0_f64*phi(nr,j)+0.5_f64*phi(nr-1,j))/dr
          grad_phi(2,1,j)=(phi(1,modulo(j+1-1+ntheta,ntheta)+1)-phi(1,modulo(j-1-1+ntheta,ntheta)+1))/(2*rmin*dtheta)
          grad_phi(2,nr+1,j)=(phi(nr+1,modulo(j+1-1+ntheta,ntheta)+1)-phi(nr+1,modulo(j-1-1+ntheta,ntheta)+1))/(2*rmax*dtheta)
       end do

    else if (calculus==2) then
       ! center formula for r, decenter on boundaries
       ! using fft for theta
       !not done

!!$       do i=2,nr
!!$          r=rmin+real(i-1,f64)*dr
!!$          grad_phi(1,i,:)=(phi(i+1,:)-phi(i-1,:))/(2*dr)
!!$       end do
!!$       do j=1,ntheta
!!$          grad_phi(1,1,j)=(phi(2,j)-phi(1,j))/dr
!!$          grad_phi(1,nr+1,j)=(phi(nr+1,j)-phi(nr,j))/dr
!!$       end do
!!$
!!$       call derivate_fft()

    else if (calculus==3) then
       ! using splines for r and theta

       call compute_spline_2D(phi,spl_phi)

       do j=1,ntheta
          theta=real(j-1,f64)*dtheta
          do i=1,nr+1
             r=rmin+real(i-1,f64)*dr
             grad_phi(1,i,j)=interpolate_x1_derivative_2D(r,theta,spl_phi)
             grad_phi(2,i,j)=interpolate_x2_derivative_2D(r,theta,spl_phi)/r
          end do
       end do

    else
       print*,'no choosen way to compute grad'
       print*,'see line 46 of file selalib/prototype/src/simulation/sll_advection_polar.F90'
       print*,'initializing grad to 0'
       grad_phi=0.0_f64
    end if

    grad_phi(:,:,ntheta+1)=grad_phi(:,:,1)

  end subroutine compute_grad_field


  !>subroutine advect_CG_polar(data,rk,f,f_fft,phi,grad_phi,spl_f,spl_phi,spl_a1,spl_a2)
  !>compute step for Center-Guide equation
  !>data : polar_data object, contains all datas about the domain
  !>rk : polar_vp_rk4 object, contains vectors for RK2 and RK4
  !>f : distribution function, size (nr+1)*(ntheta+1)
  !>f_fft and phi : used only with RK, size (nr+1)*(ntheta+1), optional
  !>f_fft : copy of f for Poisson in RK
  !>phi : solution of -Laplacian(phi)=f, it is already computed and useless if RK not used
  !>grad_phi : gradient of phi, size 2*(nr+1)*(ntheta+1)
  !>           grad_phi(1,:,:)=d_r(phi), grad_phi(2,:,:)=d_\theta(phi)/r
  !>spl_f : sll_splines_2D, spline of f
  !>spl_phi,spl_a1,spl_a2 : sll_splines_2D, optional, those are used only with RK
  !>                        respectively splines of phi, grad_phi(1,:,:), grad_phi(2,:,:)
  subroutine advect_CG_polar(data,rk,f,f_fft,phi,grad_phi,spl_f,spl_phi,spl_a1,spl_a2)

    implicit none

    type(polar_data), intent(inout), pointer :: data
    type(polar_vp_rk4), intent(inout), pointer :: rk
    sll_real64, dimension(:,:), intent(inout) :: f
    sll_real64, dimension(:,:), intent(inout), optional :: f_fft,phi
    sll_real64, dimension(:,:,:), intent(inout) :: grad_phi
    type(sll_spline_2D), pointer :: spl_f
    type(sll_spline_2D), optional, pointer :: spl_phi,spl_a1,spl_a2

    sll_int32 :: nr, ntheta
    sll_real64 :: dt, dr, dtheta, rmin, rmax
    sll_int32 :: interpolate_case
    sll_int32 :: i,j,maxiter,iter,kr,k
    sll_real64 :: r,theta,rr,rrn,ttheta,tthetan,tolr,tolth,ar,atheta

    nr=data%nr
    ntheta=data%ntheta
    dt=data%dt
    dr=data%dr
    dtheta=data%dtheta
    rmin=data%rmin
    rmax=data%rmax

    !interpolation
    ! 1 : using explicit Euler method
    ! 2 : rotation, this case ignore the field phi
    !     rotation speed = -1
    ! 3 : using RK4 // A REPRENDRE
    ! 4 : using RK2 // ib
    ! 5 : using symplectic Euler with linear interpolation
    ! 6 : using symplectic Verlet with linear interpolation
    ! 7 : using fixed point method
    interpolate_case=6
    !in grad_phi(2,:,:), the field is already divided by r, there is non need to do it here
    !hypothesis for 5, 6, 7 : field = 0 every where outside of the domain => grad_phi=0

    if (interpolate_case==1 .or. interpolate_case==2) then

       !construction of spline coefficients for f
       call compute_spline_2D(f,spl_f)

       do i=1,nr+1
          r=rmin+real(i-1,f64)*dr
          do j=1,ntheta+1
             theta=real(j-1,f64)*dtheta

             if (interpolate_case==1) then
                !Euler methode
                theta=theta-dt*grad_phi(1,i,j)/r
                r=r+dt*grad_phi(2,i,j)

             else if (interpolate_case==2) then
                !rotation
                theta=theta-dt
             end if

             call correction_r(r,rmin,rmax)
             call correction_theta(theta)
             f(i,j)=interpolate_value_2D(r,theta,spl_f)

          end do
       end do

    else if (interpolate_case==3) then
       !call rk4_polar_advect(adv,rk)!bad done

    else if (interpolate_case==4) then
       !call rk2_polar_advect(adv,rk)!bad done

    else if (interpolate_case==5) then
       !using symplectic Euler with linear interpolation
       !construction of spline coefficients
       call compute_spline_2D(f,spl_f)

       !we fix the tolerance and the maximum of iteration
       tolr=dr/5.0_f64
       tolr=1e-14
       maxiter=1000

       do j=1,ntheta
          do i=1,nr+1
             !initialization for r interpolation
             rr=rmin+real(i-1,f64)*dr+dt*grad_phi(2,i,j)
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
                   rr=rmin+real(i-1,f64)*dr+dt*grad_phi(2,k,j)
                else if (k<nr+1 .and. k>=1) then
                   rr=rmin+real(i-1,f64)*dr+dt*((1.0_f64-r)*grad_phi(2,k,j)+r*grad_phi(2,k+1,j))
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
                theta=real(j-1,f64)*dtheta-dt*((1.0_f64-r)*grad_phi(1,k,j)/(rmin+real(k-1,f64)*dr)+r*grad_phi(1,k+1,j)/(rmin+real(k+1-1,f64)*dr))
             else
                theta=real(j-1,f64)*dtheta-dt*grad_phi(1,k,j)/(rmin+real(k-1,f64)*dr)
             end if
             call correction_theta(theta)

             f(i,j)=interpolate_value_2d(rr,theta,spl_f)
          end do
       end do

    else if (interpolate_case==6) then
       !using symplectic Verlet with linear interpolation
       !construction of spline coefficients
       call compute_spline_2D(f,spl_f)

       !we fix the tolerance and the maximum of iteration
       tolr=1e-12
       tolth=1e-12
       maxiter=1000

       do j=1,ntheta
          do i=1,nr+1
             !initialization for r interpolation
             rr=rmin+real(i-1,f64)*dr+dt*0.5_f64*grad_phi(2,i,j)
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
                   rr=rmin+real(i-1,f64)*dr+0.5_f64*dt*grad_phi(2,kr,j)
                else if (kr>0 .and. kr<nr+1) then
                   rr=rmin+real(i-1,f64)*dr+0.5_f64*dt*((1.0_f64-r)*grad_phi(2,kr,j)+r*grad_phi(2,kr+1,j))
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
             ttheta=real(j-1,f64)*dtheta-dt*grad_phi(1,i,j)/(rmin+real(i-1,f64)*dr)
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
                   ttheta=real(j-1,f64)*dtheta-0.5_f64*dt*((1.0_f64-theta)*grad_phi(1,kr,k)/(rmin+real(kr-1,f64)*dr) &
                        & +theta*grad_phi(1,kr,k+1)/(rmin+real(kr-1,f64)*dr))
                   ttheta=ttheta-0.5_f64*dt*grad_phi(1,kr,j)/(rmin+real(kr-1,f64)*dr)
                else
                   ttheta=real(j-1,f64)*dtheta-0.5_f64*dt*((1.0_f64-theta) &
                        & *((1.0_f64-r)*grad_phi(1,kr,k)/(rmin+real(kr-1,f64)*dr) &
                        & +r*grad_phi(1,kr+1,k)/(rmin+real(kr+1-1,f64)*dr)) &
                        & +theta*((1.0_f64-r)*grad_phi(1,kr,k+1)/(rmin+real(kr-1,f64)*dr) &
                        & +r*grad_phi(1,kr+1,k+1)/(rmin+real(kr+1-1,f64)*dr)))
                   ttheta=ttheta-0.5_f64*dt*((1.0_f64-r)*grad_phi(1,kr,j)/(rmin+real(kr-1,f64)*dr)+r*grad_phi(1,kr+1,j)/(rmin+real(kr+1-1,f64)*dr))
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
                rr=rr+0.5_f64*dt*((1.0_f64-theta)*grad_phi(2,kr,k)+theta*grad_phi(2,kr,k+1))
             else
                rr=rr+0.5_f64*dt*((1.0_f64-theta)*((1.0_f64-r)*grad_phi(2,kr,k)+r*grad_phi(2,kr+1,k)) &
                     & +theta*((1.0_f64-r)*grad_phi(2,kr,k+1)+r*grad_phi(2,kr+1,k+1)))
             end if
             call correction_r(rr,rmin,rmax)

             f(i,j)=interpolate_value_2d(rr,ttheta,spl_f)

          end do
       end do

    else if (interpolate_case==7) then
       !using fixed point method
       !construction of spline coefficients
       call compute_spline_2D(f,spl_f)

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
                   ar=-0.5_f64*dt*((1.0_f64-theta)*((1.0_f64-r)*grad_phi(2,kr,k)+theta*((1.0_f64-r)*grad_phi(2,kr,k+1))))
                   atheta=0.5_f64*dt*((1.0_f64-theta)*((1.0_f64-r)*grad_phi(1,kr,k)/(rmin+real(kr-1,f64)*dr)+theta*((1.0_f64-r)*grad_phi(1,kr,k+1)/(rmin+real(kr-1,f64)*dr))))
                else
                   ar=-0.5_f64*dt*((1.0_f64-theta)*((1.0_f64-r)*grad_phi(2,kr,k)+r*grad_phi(2,kr+1,k)) &
                        & +theta*((1.0_f64-r)*grad_phi(2,kr,k+1)+r*grad_phi(2,kr+1,k+1)))
                   atheta=0.5_f64*dt*((1.0_f64-theta)*((1.0_f64-r)*grad_phi(1,kr,k)/(rmin+real(kr-1,f64)*dr)+r*grad_phi(1,kr+1,k)/(rmin+real(kr+1-1,f64)*dr)) &
                        & +theta*((1.0_f64-r)*grad_phi(1,kr,k+1)/(rmin+real(kr-1,f64)*dr)+r*grad_phi(1,kr+1,k+1)/(rmin+real(kr+1-1,f64)*dr)))
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
             f(i,j)=interpolate_value_2d(rr,ttheta,spl_f)
          end do
       end do

    else
       print*,'no way chosen to compute r and theta'
       print*,'nothing will be done'
       print*,'see variable interpolate_case, routine advect_VP_polar in file prototype/src/simulation/sll_advection_polar.F90'
    end if

    f(:,ntheta+1)=f(:,1)

  end subroutine advect_CG_polar


  !>subroutine rk4_polar_advect(adv,rk)
  !>RK4 for polar advection only
  !>adv : polar_vp_data object, all data are included in
  !>rk : polar_vp_rk4 object, array for rk4
  subroutine rk4_polar_advect(adv,rk)
!A REPRENDRE
    implicit none

    type(polar_vp_data), intent(inout), pointer :: adv
    type(polar_vp_rk4), intent(inout), pointer :: rk

    sll_real64 :: dt, dr, dtheta, rmin, rmax
    sll_int32 :: nr, ntheta
    sll_int32 :: i,j
    sll_real64 :: r,theta,rr

    nr=adv%data%nr
    ntheta=adv%data%ntheta
    dt=adv%data%dt
    dr=adv%data%dr
    dtheta=adv%data%dtheta
    rmin=adv%data%rmin
    rmax=adv%data%rmax

    !construction of spline coeficients for f
    call compute_spline_2D(adv%f,adv%spl_f)

    !first step of RK4
    do i=1,nr+1
       !r=rmin+real(i-1)*dr
       rk%r1(i,:)=-dt*adv%grad_phi(2,i,:)/adv%rr(i)
       rk%theta1(i,:)=dt*adv%grad_phi(1,i,:)/adv%rr(i)
    end do

    !2nd step of RK4
    do i=1,nr+1
       !r=rmin+real(i-1)*dr
       do j=1,ntheta+1
          theta=real(j-1,f64)*dtheta
          rk%r4(i,j)=adv%rr(i)-dt/2.0_f64*adv%grad_phi(2,i,j)/adv%rr(i)
          rk%theta4(i,j)=adv%ttheta(j)+dt/2.0_f64*adv%grad_phi(1,i,j)/adv%rr(i)

          call correction_r(rk%r4(i,j),rmin,rmax)
          call correction_theta(rk%theta4(i,j))

          adv%f(i,j)=interpolate_value_2D(rk%r4(i,j),rk%theta4(i,j),adv%spl_f)
       end do
    end do

    !call poisson_solve_polar(adv)
    !call compute_grad_field(adv)

    !construction of spline coeficients for a
    call compute_spline_2d(adv%grad_phi(1,:,:),adv%spl_a1)
    call compute_spline_2d(adv%grad_phi(2,:,:),adv%spl_a2)

    do i=1,nr+1
       !rr=1.0_f64!rmin+real(i-1)*dr
       do j=1,ntheta+1
          r=adv%rr(i)+rk%r1(i,j)/2.0_f64
          theta=adv%ttheta(j)+rk%theta1(i,j)/2.0_f64
          call correction_r(r,rmin,rmax)
          call correction_theta(theta)
          rk%r2(i,j)=-interpolate_value_2d(r,theta,adv%spl_a2)*dt/rk%r4(i,j)
          rk%theta2(i,j)=interpolate_value_2d(r,theta,adv%spl_a1)*dt/rk%r4(i,j)
       end do
    end do

    !3rd step of RK4
    do i=1,nr+1
       !rr=1.0_f64!rmin+real(i-1)*dr
       do j=1,ntheta+1
          r=adv%rr(i)+rk%r2(i,j)/2.0_f64
          theta=adv%ttheta(j)+rk%theta2(i,j)/2.0_f64
          call correction_r(r,rmin,rmax)
          call correction_theta(theta)
          rk%r4(i,j)=rk%r4(i,j)*adv%rr(i)

          rk%r3(i,j)=-interpolate_value_2d(r,theta,adv%spl_a2)*dt/rk%r4(i,j)
          rk%theta3(i,j)=interpolate_value_2d(r,theta,adv%spl_a1)*dt/rk%r4(i,j)
          !print*,'3',i,j
       end do
    end do

    !4th step of RK4
    do i=1,nr+1
       r=adv%rr(i)!rmin+real(i-1)*dr
       do j=1,ntheta+1
          theta=adv%ttheta(j)!real(j-i)*dtheta
          rk%r4(i,j)=r+rk%r1(i,j)
          rk%theta4(i,j)=theta+rk%theta1(i,j)
          call correction_r(rk%r4(i,j),rmin,rmax)
          call correction_theta(rk%theta4(i,j))

          adv%f(i,j)=interpolate_value_2D(rk%r4(i,j),rk%theta4(i,j),adv%spl_f)
       end do
    end do

    !call poisson_solve_polar(adv)
    !call compute_grad_field(adv)

    !construction of spline coeficients for a
    call compute_spline_2d(adv%grad_phi(1,:,:),adv%spl_a1)
    call compute_spline_2d(adv%grad_phi(2,:,:),adv%spl_a2)

    do i=1,nr+1
       rr=r*adv%rr(i)
       do j=1,ntheta+1
          r=adv%rr(i)+rk%r3(i,j)
          theta=adv%ttheta(j)+rk%theta3(i,j)
          call correction_r(r,rmin,rmax)
          call correction_theta(theta)

          rk%r4(i,j)=-interpolate_value_2d(r,theta,adv%spl_a2)*dt/rr
          rk%theta4(i,j)=interpolate_value_2d(r,theta,adv%spl_a1)*dt/rr
       end do
    end do

    !sommation
    do i=1,nr+1
       r=adv%rr(i)!rmin+real(i-1,f64)*dr
       do j=1,ntheta+1
          theta=adv%ttheta(j)!real(j-1,f64)*dtheta
          rk%r1(i,j)=r+rk%r1(i,j)/6.0_f64+rk%r2(i,j)/3.0_f64+rk%r3(i,j)/3.0_f64+rk%r4(i,j)/6.0_f64
          rk%theta1(i,j)=theta+rk%theta1(i,j)/6.0_f64+rk%theta2(i,j)/3.0_f64+rk%theta3(i,j)/3.0_f64+rk%theta4(i,j)/6.0_f64

          call correction_r(rk%r1(i,j),rmin,rmax)
          call correction_theta(rk%theta1(i,j))
       end do
    end do

    !updating the distribution function

  end subroutine rk4_polar_advect


  !>subroutine rk2_polar_advect(adv,rk)
  !>RK2 for polar advection only
  !>adv : polar_vp_data object
  !>rk : polar_vp_rk4 object, contains more array than necessary
  subroutine rk2_polar_advect(adv,rk)
!A REPRENDRE
    implicit none

    type(polar_vp_data), intent(inout), pointer :: adv
    type(polar_vp_rk4), intent(inout), pointer :: rk

    sll_real64 :: dt, dr, dtheta, rmin, rmax
    sll_int32 :: nr, ntheta
    sll_int32 :: i,j
    sll_real64 :: r,theta

    nr=adv%data%nr
    ntheta=adv%data%ntheta
    dt=adv%data%dt
    dr=adv%data%dr
    dtheta=adv%data%dtheta
    rmin=adv%data%rmin
    rmax=adv%data%rmax

    !construction of spline coeficients for f
    call compute_spline_2D(adv%f,adv%spl_f)

    do i=1,nr+1
       do j=1,ntheta+1
          rk%r1(i,j)=adv%rr(i)+dt/2.0_f64*adv%grad_phi(2,i,j)/adv%rr(i)
          rk%theta1(i,j)=adv%ttheta(j)-dt/2.0_f64*adv%grad_phi(1,i,j)/adv%rr(i)
          call correction_r(rk%r1(i,j),rmin,rmax)
          call correction_theta(rk%theta1(i,j))

          adv%f(i,j)=interpolate_value_2D(rk%r1(i,j),rk%theta1(i,j),adv%spl_f)
       end do
    end do

    !call poisson_solve_polar(adv)
    !call compute_grad_field(adv)

    !construction of spline coeficients for a
    do i=1,nr+1
       adv%grad_phi(:,i,:)=adv%grad_phi(:,i,:)*adv%rr(i)
    end do
    call compute_spline_2d(adv%grad_phi(1,:,:),adv%spl_a1)
    call compute_spline_2d(adv%grad_phi(2,:,:),adv%spl_a2)

    do i=1,nr+1
       do j=1,ntheta+1
          r=rk%r1(i,j)
          theta=rk%theta1(i,j)

          rk%r2(i,j)=interpolate_value_2d(r,theta,adv%spl_a2)/r
          rk%theta2(i,j)=-interpolate_value_2d(r,theta,adv%spl_a1)/r
       end do
    end do

    !sommation
    do i=1,nr+1
       do j=1,ntheta+1
          rk%r1(i,j)=adv%rr(i)+dt*rk%r2(i,j)
          rk%theta1(i,j)=adv%ttheta(j)+dt*rk%theta2(i,j)

          call correction_r(rk%r1(i,j),rmin,rmax)
          call correction_theta(rk%theta1(i,j))

          !updating distribution function f
          adv%f(i,j)=interpolate_value_2d(rk%r1(i,j),rk%theta1(i,j),adv%spl_f)
       end do
    end do

  end subroutine rk2_polar_advect


  !>subroutine SL_classic(data,rk,f,f_fft,phi,grad_phi,pfwd,pinv,fk,phik,spl_f,spl_phi,spl_a1,spl_a2,a,cts,ipiv)
  !>computes the classic semi-Lagrangian scheme for Vlasov-Poisson equation
  !>data: polar_data object, contains data about the domain
  !>rk : polar_vp_rk4 object, array for rk4
  !>f, f_fft, phi : distribution function, copy for Poisson and field solution of -Laplacian(phi)=f
  !>                size (nr+1)*(ntheta+1)
  !>grad_phi : gradient(phi), size 2*(nr+1)*(ntheta+1)
  !>pfwd and pinv : sll_fft_plan object for fft forward and inverse
  !>fk and phik : complex arrays for Poisson, size ntheta/2+1
  !>spl_f, spl_phi, spl_a1,spl_a2 : sll_spline_2D objects for spline
  !>spl_phi is needed only if the grad is computed with spline
  !>spl_a1 and spl_a2 are needed only if advect_CG_polar uses RK
  !>a, cts, ipiv : objects for the tridiagonal solver, see the tridiagonal solver section
  subroutine SL_classic(data,rk,f,f_fft,phi,grad_phi,pfwd,pinv,fk,phik,spl_f,spl_phi,spl_a1,spl_a2,a,cts,ipiv)

    implicit none

    type(polar_data), intent(inout), pointer :: data
    type(polar_vp_rk4), intent(inout), pointer :: rk
    sll_real64, dimension(:,:), intent(inout) :: f,f_fft,phi
    sll_real64, dimension(:,:,:), intent(inout) :: grad_phi
    type(sll_fft_plan), intent(inout), pointer :: pfwd,pinv
    sll_comp64, dimension(:), intent(inout) :: fk,phik
    type(sll_spline_2D), intent(inout), pointer :: spl_f
    type(sll_spline_2D), intent(inout), optional, pointer :: spl_phi, spl_a1, spl_a2
    sll_real64, dimension(:), intent(inout) :: a,cts
    sll_int32, dimension(:), intent(inout) :: ipiv


    if ( .not. present(spl_phi)) then
       spl_phi => NULL()
    end if
    if ( .not. present(spl_a1)) then
       spl_a1 => NULL()
    end if
    if ( .not. present(spl_a2)) then
       spl_a2 => NULL()
    end if

    call poisson_solve_polar(data,f,phi,f_fft,fk,phik,a,cts,ipiv,pfwd,pinv)
    call compute_grad_field(data,phi,grad_phi,spl_phi)
    call advect_CG_polar(data,rk,f,f_fft,phi,grad_phi,spl_f,spl_phi,spl_a1,spl_a2)

  end subroutine SL_classic


  !>subroutine SL_controlled(data,rk,f,fdemi,f_fft,phi,grad_phi,pfwd,pinv,spl_f,spl_phi,spl_a1,spl_a2,a,cts,ipiv)
  !>computes the semi-Lagrangian scheme order 2
  !>data: polar_data object, contains data about the domain
  !>rk : polar_vp_rk4 object, array for rk4
  !>f, fdemi, f_fft, phi : distribution function, copy for Poisson and field solution of -Laplacian(phi)=f
  !>                       size (nr+1)*(ntheta+1)
  !>grad_phi : gradient(phi), size 2*(nr+1)*(ntheta+1)
  !>pfwd and pinv : sll_fft_plan object for fft forward and inverse
  !>fk and phik : complex arrays for Poisson, size ntheta/2+1
  !>spl_f, spl_phi, spl_a1,spl_a2 : sll_spline_2D objects for spline
  !>spl_phi is needed only if the grad is computed with spline
  !>spl_a1 and spl_a2 are needed only if advect_CG_polar uses RK
  !>a, cts, ipiv : objects for the tridiagonal solver, see the tridiagonal solver section
  subroutine SL_ordre_2(data,rk,f,fdemi,f_fft,phi,grad_phi,pfwd,pinv,fk,phik,spl_f,spl_phi,spl_a1,spl_a2,a,cts,ipiv)

    implicit none

    type(polar_data), intent(inout), pointer :: data
    type(polar_vp_rk4), intent(inout), pointer :: rk
    sll_real64, dimension(:,:), intent(inout) :: f,fdemi,f_fft,phi
    sll_real64, dimension(:,:,:), intent(inout) :: grad_phi
    type(sll_fft_plan), intent(inout), pointer :: pfwd,pinv
    sll_comp64, dimension(:), intent(inout) :: fk,phik
    type(sll_spline_2D), intent(inout), pointer :: spl_f
    type(sll_spline_2D), intent(inout), optional, pointer :: spl_phi, spl_a1, spl_a2
    sll_real64, dimension(:), intent(inout) :: a,cts
    sll_int32, dimension(:), intent(inout) :: ipiv

    sll_real64 :: dt

    fdemi=f
    dt=data%dt
    data%dt=dt/2.0_f64

    if ( .not. present(spl_phi)) then
       spl_phi => NULL()
    end if
    if ( .not. present(spl_a1)) then
       spl_a1 => NULL()
    end if
    if ( .not. present(spl_a2)) then
       spl_a2 => NULL()
    end if

    call poisson_solve_polar(data,f,phi,f_fft,fk,phik,a,cts,ipiv,pfwd,pinv)
    call compute_grad_field(data,phi,grad_phi,spl_phi)
    call advect_CG_polar(data,rk,f,f_fft,phi,grad_phi,spl_f,spl_phi,spl_a1,spl_a2)
    !we just obtained f^(n+1/2)
    call poisson_solve_polar(data,f,phi,f_fft,fk,phik,a,cts,ipiv,pfwd,pinv)
    call compute_grad_field(data,phi,grad_phi,spl_phi)
    !we just obtained E^(n+1/2)
    data%dt=dt
    f=fdemi
    call advect_CG_polar(data,rk,f,f_fft,phi,grad_phi,spl_f,spl_phi,spl_a1,spl_a2)

  end subroutine SL_ordre_2


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


  !>subroutine divergence_ortho_field(data,grad_phi,div)
  !>compute divergence of the field used for CG equations
  !>this field is orthogonal to grad_phi
  !>data : polar_data object
  !>grad_phi : size 2*(nr+1)*(ntheta+1), grad_phi(1,:,:)=d_r(phi), grad_phi(2,:,:)=d_theta(phi)
  !>div : polar divergence of grad_phi at point (i,j), size (nr+1)*(ntheta+1)
  subroutine divergence_ortho_field(data,grad_phi,div)

    implicit none

    type(polar_data), intent(in), pointer :: data
    sll_real64, dimension(:,:,:), intent(inout) :: grad_phi
    sll_real64, dimension(:,:), intent(out) :: div

    sll_real64 :: dr,dtheta,rmin,rmax
    sll_int32 :: nr,ntheta
    sll_real64 :: r
    sll_int32 :: i,j

    nr=data%nr
    ntheta=data%ntheta
    dr=data%dr
    dtheta=data%dtheta
    rmin=data%rmin
    rmax=data%rmax

    div=-grad_phi(2,:,:)
    grad_phi(2,:,:)=grad_phi(1,:,:)
    grad_phi(1,:,:)=div
    div=0.0_f64

    do i=2,nr
       r=rmin+real(i-1,f64)*dr
       do j=1,ntheta
          div(i,j)=1/r*((grad_phi(1,i+1,j)*(r+dr)-grad_phi(1,i-1,j)*(r-dr))/(2*dr) &
               & +(grad_phi(2,i,modulo(j+1-1+ntheta,ntheta+1)+1)-grad_phi(2,i,modulo(j-1-1+ntheta,ntheta+1)+1))/(2*dtheta))
       end do
    end do
    do j=1,ntheta
       div(1,j)=1.0_f64/rmin*((grad_phi(1,2,j)*(rmin+dr)-grad_phi(1,1,j)*(rmin))/(dr) &
               & +(grad_phi(2,1,modulo(j+1-1+ntheta,ntheta+1)+1)-grad_phi(2,1,modulo(j-1-1+ntheta,ntheta+1)+1))/(2*dtheta))
       div(nr+1,j)=1.0_f64/r*((grad_phi(1,nr+1,j)*(rmax)-grad_phi(1,nr+1-1,j)*(rmax-dr))/(2*dr) &
               & +(grad_phi(2,nr+1,modulo(j+1-1+ntheta,ntheta+1)+1)-grad_phi(2,nr+1,modulo(j-1-1+ntheta,ntheta+1)+1))/(2*dtheta))
!!$       div(1,j)=1.0_f64/rmin*(-(1.5_f64*grad_phi(1,1,j)*rmin-2.0_f64*grad_phi(1,2,j)*(rmin+dr)+0.5_f64*grad_phi(1,3,j)*(rmin+2.0_f64*dr))/(dr) &
!!$               & +(grad_phi(2,1,modulo(j+1-1+ntheta,ntheta+1)+1)-grad_phi(2,1,modulo(j-1-1+ntheta,ntheta+1)+1))/(2*dtheta))
!!$       div(nr+1,j)=1.0_f64/rmax*(-(1.5_f64*grad_phi(1,nr+1,j)*rmin-2.0_f64*grad_phi(1,nr,j)*(rmin+dr)+0.5_f64*grad_phi(1,nr-1,j)*(rmin+2.0_f64*dr))/(2*dr) &
!!$               & +(grad_phi(2,nr+1,modulo(j+1-1+ntheta,ntheta+1)+1)-grad_phi(2,nr+1,modulo(j-1-1+ntheta,ntheta+1)+1))/(2*dtheta))
    end do


!!$          grad_phi(1,1,j)=-(1.5_f64*phi(1,j)-2.0_f64*phi(2,j)+0.5_f64*phi(3,j))/dr
!!$          grad_phi(1,nr+1,j)=-(1.5_f64*phi(nr+1,j)-2.0_f64*phi(nr,j)+0.5_f64*phi(nr-1,j))/dr

    div(:,ntheta+1)=div(:,1)

  end subroutine divergence_ortho_field

end module polar_advection
