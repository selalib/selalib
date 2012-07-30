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

  !>subroutine compute_grad_field(adv)
  !>compute a = grad(phi) for phi scalar field
  !>adv : polar_vp_data object, all datas are included in
  subroutine compute_grad_field(adv)

    implicit none

    type(polar_vp_data), intent(inout), pointer :: adv

    sll_int32 :: nr, ntheta
    sll_real64 :: dr, dtheta, rmin, rmax
    sll_int32 :: i,j,calculus
    sll_real64 :: r,theta

    nr=adv%data%nr
    ntheta=adv%data%ntheta
    dr=adv%data%dr
    dtheta=adv%data%dtheta
    rmin=adv%data%rmin
    rmax=adv%data%rmax

    ! way to compute a
    ! 1 : center for r and theta
    !     decenter on boundaries
    ! 2 : center for r, decenter on boundaries
    !     using fft for theta
    ! 3 : center for r, decenter on boundaries (not writen)
    !     using splines for theta
    calculus = 1

    if (calculus==1) then
       ! center formula for r end theta
       ! decenter for r on boundaries

       do i=2,nr
          r=rmin+real(i-1,f64)*dr
          do j=1,ntheta
             adv%grad_phi(1,i,j)=(adv%phi(i+1,j)-adv%phi(i-1,j))/(2*dr)
             adv%grad_phi(2,i,j)=(adv%phi(i,modulo(j+1-1+ntheta,ntheta)+1)-adv%phi(i,modulo(j-1-1+ntheta,ntheta)+1))/(2*r*dtheta)
          end do
       end do
       do j=1,ntheta
          adv%grad_phi(1,1,j)=(adv%phi(2,j)-adv%phi(1,j))/dr
          adv%grad_phi(1,nr+1,j)=(adv%phi(nr+1,j)-adv%phi(nr,j))/dr
          adv%grad_phi(2,1,j)=(adv%phi(1,modulo(j+1-1+ntheta,ntheta)+1)-adv%phi(1,modulo(j-1-1+ntheta,ntheta)+1))/(2*rmin*dtheta)
          adv%grad_phi(2,nr+1,j)=(adv%phi(nr+1,modulo(j+1-1+ntheta,ntheta)+1)-adv%phi(nr+1,modulo(j-1-1+ntheta,ntheta)+1))/(2*rmax*dtheta)
       end do

    else if (calculus==2) then
       ! center formula for r, decenter on boundaries
       ! using fft for theta

       do i=2,nr
          r=rmin+real(i-1,f64)*dr
          adv%grad_phi(1,i,:)=(adv%phi(i+1,:)-adv%phi(i-1,:))/(2*dr)
       end do
       do j=1,ntheta
          adv%grad_phi(1,1,j)=(adv%phi(2,j)-adv%phi(1,j))/dr
          adv%grad_phi(1,nr+1,j)=(adv%phi(nr+1,j)-adv%phi(nr,j))/dr
       end do

       call derivate_fft(adv)

    else if (calculus==3) then
       ! center formula for r, decenter on boundaries
       ! using splines for theta

       do i=1,nr+1
          r=rmin+real(i,f64)*dr
          adv%f_fft(i,:)=adv%phi(i,:)*r
       end do
       call compute_spline_2D(adv%f_fft,adv%spl_phi)

       do i=2,nr
          r=rmin+real(i-1,f64)*dr
          adv%grad_phi(1,i,:)=(adv%phi(i+1,:)-adv%phi(i-1,:))/(2*dr)
       end do
       do j=1,ntheta
          adv%grad_phi(1,1,j)=(adv%phi(2,j)-adv%phi(1,j))/dr
          adv%grad_phi(1,nr+1,j)=(adv%phi(nr+1,j)-adv%phi(nr,j))/dr
          theta=real(j-1,f64)*dtheta
          do i=1,nr+1
             r=rmin+real(i-1,f64)*dr
             adv%grad_phi(2,i,j)=interpolate_x2_derivative_2D(r,theta,adv%spl_phi)/r
          end do
       end do

    else
       print*,'no choosen way to compute grad'
       print*,'see line 46 of file selalib/prototype/src/simulation/sll_advection_polar.F90'
       print*,'initializing grad to 0'
       adv%grad_phi=0.0_f64
    end if

    adv%grad_phi(:,:,ntheta+1)=adv%grad_phi(:,:,1)

  end subroutine compute_grad_field


  !>subroutine advect_CG_polar(adv)
  !>compute step for Center-Guide equation
  !>adv : polar_vp_data object, all data are included in
  subroutine advect_CG_polar(adv,rk)

    implicit none

    type(polar_vp_data), intent(inout), pointer :: adv
    type(polar_vp_rk4), intent(inout), pointer :: rk

    sll_int32 :: nr, ntheta
    sll_real64 :: dt, dr, dtheta, rmin, rmax
    sll_int32 :: interpolate_case
    sll_int32 :: i,j
    sll_real64 :: r, theta

    nr=adv%data%nr
    ntheta=adv%data%ntheta
    dt=adv%data%dt
    dr=adv%data%dr
    dtheta=adv%data%dtheta
    rmin=adv%data%rmin
    rmax=adv%data%rmax

    !interpolation
    ! 1 : using explicit Eulerian scheme
    ! 2 : rotation, this case ignore the field phi
    !     rotation speed = -1
    ! 3 : using RK4
    ! 4 : using RK2
    interpolate_case=1

    if (interpolate_case<3 .and. interpolate_case>0) then

       !construction of spline coefficients for f
       do i=1,nr+1
          adv%f_fft(i,:)=adv%f(i,:)!*adv%rr(i)
       end do
       call compute_spline_2D(adv%f_fft,adv%spl_f)

       do i=1,nr+1
          r=adv%rr(i)
          do j=1,ntheta+1
             theta=adv%ttheta(j)

             if (interpolate_case==1) then
                !Eulerian scheme
                theta=theta+dt*adv%grad_phi(1,i,j)/r
                r=r-dt*adv%grad_phi(2,i,j)/r

             else if (interpolate_case==2) then
                !rotation
                theta=theta+dt
             end if

             call correction_r(r,rmin,rmax)
             call correction_theta(theta)
             adv%f(i,j)=interpolate_value_2D(r,theta,adv%spl_f)!/adv%rr(i)

          end do
       end do

    else if (interpolate_case==3) then
       call rk4_polar_advect(adv,rk)

    else if (interpolate_case==4) then
       call rk2_polar_advect(adv,rk)

    else
       print*,'no way chosen to compute r and theta'
       print*,'nothing will be done'
       print*,'see variable interpolate_case, routine advect_VP_polar in file prototype/src/simulation/sll_advection_polar.F90'
    end if

    adv%f(:,ntheta+1)=adv%f(:,1)

  end subroutine advect_CG_polar


  !>subroutine rk4_polar_advect(adv,rk)
  !>RK4 for polar advection only
  !>adv : polar_vp_data object, all data are included in
  !>rk : polar_vp_rk4 object, array for rk4
  subroutine rk4_polar_advect(adv,rk)

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
    do i=1,nr+1
       !r=1.0_f64!rmin+real(i,f64)*dr
       adv%f_fft(i,:)=adv%f(i,:)*adv%rr(i)
    end do
    call compute_spline_2D(adv%f_fft,adv%spl_f)

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

          adv%f(i,j)=interpolate_value_2D(rk%r4(i,j),rk%theta4(i,j),adv%spl_f)/adv%rr(i)
       end do
    end do

    call poisson_solve_polar(adv)
    call compute_grad_field(adv)

    !construction of spline coeficients for a
    do i=1,nr+1
       !r=1.0_f64!rmin+real(i,f64)*dr
       adv%grad_phi(:,i,:)=adv%grad_phi(:,i,:)*adv%rr(i)
    end do
    call compute_spline_2d(adv%grad_phi(1,:,:),adv%spl_a1)
    call compute_spline_2d(adv%grad_phi(2,:,:),adv%spl_a2)

    do i=1,nr+1
       !rr=1.0_f64!rmin+real(i-1)*dr
       do j=1,ntheta+1
          r=adv%rr(i)+rk%r1(i,j)/2.0_f64
          theta=adv%ttheta(j)+rk%theta1(i,j)/2.0_f64
          call correction_r(r,rmin,rmax)
          call correction_theta(theta)
          rk%r2(i,j)=-interpolate_value_2d(r,theta,adv%spl_a2)*dt/(rk%r4(i,j)*adv%rr(i))
          rk%theta2(i,j)=interpolate_value_2d(r,theta,adv%spl_a1)*dt/(rk%r4(i,j)*adv%rr(i))
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
!print*,'rk 3 ok'
    !4th step of RK4
    do i=1,nr+1
       r=adv%rr(i)!rmin+real(i-1)*dr
       do j=1,ntheta+1
          theta=adv%ttheta(j)!real(j-i)*dtheta
          rk%r4(i,j)=r+rk%r1(i,j)
          rk%theta4(i,j)=theta+rk%theta1(i,j)
          call correction_r(rk%r4(i,j),rmin,rmax)
          call correction_theta(rk%theta4(i,j))

          adv%f(i,j)=interpolate_value_2D(rk%r4(i,j),rk%theta4(i,j),adv%spl_f)/r
          !print*,'41',i,j
       end do
    end do

    call poisson_solve_polar(adv)
    call compute_grad_field(adv)

    !construction of spline coeficients for a
    do i=1,nr+1
       r=adv%rr(i)!1.0_f64!rmin+real(i,f64)*dr
       adv%grad_phi(:,i,:)=adv%grad_phi(:,i,:)*r
    end do
    call compute_spline_2d(adv%grad_phi(1,:,:),adv%spl_a1)
    call compute_spline_2d(adv%grad_phi(2,:,:),adv%spl_a2)
!print*,'42'
    do i=1,nr+1
       rr=r*adv%rr(i)
       do j=1,ntheta+1
          r=adv%rr(i)+rk%r3(i,j)
          theta=adv%ttheta(j)+rk%theta3(i,j)
          call correction_r(r,rmin,rmax)
          call correction_theta(theta)

          rk%r4(i,j)=-interpolate_value_2d(r,theta,adv%spl_a2)*dt/rr
          rk%theta4(i,j)=interpolate_value_2d(r,theta,adv%spl_a1)*dt/rr
          !!print*,'43',i,j
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
    call deposit_value_2D(rk%r1,rk%theta1,adv%spl_f,adv%f)
    do i=1,nr+1
       adv%f(i,:)=adv%f(i,:)/adv%rr(i)
    end do

  end subroutine rk4_polar_advect


  !>subroutine rk2_polar_advect(adv,rk)
  !>RK2 for polar advection only
  !>adv : polar_vp_data object
  !>rk : polar_vp_rk4 object, contains more array than necessary
  subroutine rk2_polar_advect(adv,rk)

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
    do i=1,nr+1
       !r=1.0_f64!rmin+real(i,f64)*dr
       adv%f_fft(i,:)=adv%f(i,:)*adv%rr(i)
    end do
    call compute_spline_2D(adv%f_fft,adv%spl_f)

    !first step of rk2
    do i=1,nr+1
       do j=1,ntheta+1
          rk%r1(i,j)=adv%rr(i)-dt/2.0_f64*adv%grad_phi(2,i,j)/adv%rr(i)
          rk%theta1(i,j)=adv%ttheta(j)+dt/2.0_f64*adv%grad_phi(1,i,j)/adv%rr(i)
       end do
    end do

    !2nd step of rk2
    do i=1,nr+1
       do j=1,ntheta+1
          theta=real(j-1,f64)*dtheta
          rk%r4(i,j)=adv%rr(i)-dt/2.0_f64*adv%grad_phi(2,i,j)/adv%rr(i)
          rk%theta4(i,j)=adv%ttheta(j)+dt/2.0_f64*adv%grad_phi(1,i,j)/adv%rr(i)

          call correction_r(rk%r4(i,j),rmin,rmax)
          call correction_theta(rk%theta4(i,j))

          adv%f(i,j)=interpolate_value_2D(rk%r4(i,j),rk%theta4(i,j),adv%spl_f)/adv%rr(i)
       end do
    end do

    call poisson_solve_polar(adv)
    call compute_grad_field(adv)

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
          call correction_r(r,rmin,rmax)
          call correction_theta(theta)

          rk%r2(i,j)=-interpolate_value_2d(r,theta,adv%spl_a2)*dt/(rk%r4(i,j)*adv%rr(i))
          rk%theta2(i,j)=interpolate_value_2d(r,theta,adv%spl_a1)*dt/(rk%r4(i,j)*adv%rr(i))
       end do
    end do

    !sommation
    do i=1,nr+1
       do j=1,ntheta+1
          rk%r1(i,j)=adv%rr(i)+dt*rk%r2(i,j)
          rk%theta1(i,j)=adv%ttheta(j)+dt*rk%theta2(i,j)

          call correction_r(rk%r1(i,j),rmin,rmax)
          call correction_theta(rk%theta1(i,j))
       end do
    end do

    !updating the distribution function
    call deposit_value_2D(rk%r1,rk%theta1,adv%spl_f,adv%f)
    do i=1,nr+1
       adv%f(i,:)=adv%f(i,:)/adv%rr(i)
    end do

  end subroutine rk2_polar_advect


  !>subroutine SL_classic(adv,rk)
  !>computes the classic semi-Lagrangian scheme for Vlasov-Poisson equation
  !>adv : polar_vp_data object, all data are included in
  !>rk : polar_vp_rk4 object, array for rk4
  subroutine SL_classic(adv,rk)

    implicit none

    type(polar_vp_data), intent(inout), pointer :: adv
    type(polar_vp_rk4), intent(inout), pointer :: rk

    call poisson_solve_polar(adv)
    call compute_grad_field(adv)
    call advect_CG_polar(adv,rk)

  end subroutine SL_classic


  !>subroutine SL_controlled(adv,rk)
  !>computes the semi-Lagrangian scheme with control for Vlasov-Poisson equation
  !>adv : polar_vp_data object, all data are included in
  !>rk : polar_vp_rk4 object, array for rk4
  subroutine SL_ordre_2(adv,rk)

    implicit none

    type(polar_vp_data), intent(inout), pointer :: adv
    type(polar_vp_rk4), intent(inout), pointer :: rk

    sll_real64 :: dt

    adv%fdemi=adv%f
    dt=adv%data%dt
    adv%data%dt=dt/2.0_f64

    call poisson_solve_polar(adv)
    call compute_grad_field(adv)
    call advect_CG_polar(adv,rk)
    !we just obtained f^(n+1/2)
    call poisson_solve_polar(adv)
    call compute_grad_field(adv)
    !we just obtained E^(n+1/2)
    adv%data%dt=dt
    adv%f=adv%fdemi
    call advect_CG_polar(adv,rk)

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
       print*,'POney'
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


  !>subroutine divergence_ortho_field(adv)
  !>compute divergence of the field used for CG equations
  !>this field is orthogonal to adv%grad_phi
  !>adv : polar_vp_data object
  subroutine divergence_ortho_field(adv,div)

    implicit none

    type(polar_vp_data), intent(in), pointer :: adv
    sll_real64, dimension(:,:), intent(out) :: div

    sll_real64 :: dr,dtheta,rmin,rmax
    sll_int32 :: nr,ntheta
    sll_real64 :: r
    sll_int32 :: i,j

    nr=adv%data%nr
    ntheta=adv%data%ntheta
    dr=adv%data%dr
    dtheta=adv%data%dtheta
    rmin=adv%data%rmin
    rmax=adv%data%rmax

    div=-adv%grad_phi(2,:,:)
    adv%grad_phi(2,:,:)=adv%grad_phi(1,:,:)
    adv%grad_phi(1,:,:)=div
    div=0.0_f64

    do i=2,nr
       r=rmin+real(i-1,f64)*dr
       do j=1,ntheta
          div(i,j)=1/r*((adv%grad_phi(1,i+1,j)*(r+dr)-adv%grad_phi(1,i-1,j)*(r-dr))/(2*dr) &
               & +(adv%grad_phi(2,i,modulo(j+1-1+ntheta,ntheta+1)+1)-adv%grad_phi(2,i,modulo(j-1-1+ntheta,ntheta+1)+1))/(2*dtheta))
       end do
    end do
    do j=1,ntheta
       div(1,j)=1/rmin*((adv%grad_phi(1,2,j)*(rmin+dr)-adv%grad_phi(1,i,j)*(rmin))/(dr) &
               & +(adv%grad_phi(2,1,modulo(j+1-1+ntheta,ntheta+1)+1)-adv%grad_phi(2,1,modulo(j-1-1+ntheta,ntheta+1)+1))/(2*dtheta))
       div(nr+1,j)=1/r*((adv%grad_phi(1,nr+1,j)*(rmax)-adv%grad_phi(1,nr+1-1,j)*(rmax-dr))/(2*dr) &
               & +(adv%grad_phi(2,nr+1,modulo(j+1-1+ntheta,ntheta+1)+1)-adv%grad_phi(2,nr+1,modulo(j-1-1+ntheta,ntheta+1)+1))/(2*dtheta))
    end do

    div(:,ntheta+1)=div(:,1)

  end subroutine divergence_ortho_field

end module polar_advection


