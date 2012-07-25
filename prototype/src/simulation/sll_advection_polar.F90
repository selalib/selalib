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

  !>subroutine compute_grad_field(nr,ntheta,dr,dtheta,rmin,rmax,phi,a)
  !>compute a = grad(phi) for phi scalar field
  !>nr and ntheta : numbers of steps in directions r and theta
  !>dr and dtheta : size of step in directions r and theta
  !>rmin : radius of the hole
  !>rmax : radius of the disc
  !>phi : known field, size (nr+1)X(ntheta+1)
  !>a = grad(phi), size 2X(nr+1)X(ntheta+1)
  !>a(1,i,j) = d_r(phi)(r_i,theta_j)
  !>a(2,i,j) = d_theta(phi)(r_i,theta_j)
  subroutine compute_grad_field(adv)

    implicit none

    type(polar_vp_data), intent(inout), pointer :: adv

    sll_int32 :: nr, ntheta
    sll_real64 :: dr, dtheta, rmin
    sll_int32 :: i,j,calculus
    sll_real64 :: r,theta

    nr=adv%data%nr
    ntheta=adv%data%ntheta
    dr=adv%data%dr
    dtheta=adv%data%dtheta
    rmin=adv%data%rmin

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
             adv%grad_phi(2,i,j)=(adv%phi(i,modulo(j+1-1+ntheta,ntheta)+1)-adv%phi(i,modulo(j-1-1+ntheta,ntheta)+1))/(2*dtheta)
          end do
       end do
       do j=1,ntheta
          adv%grad_phi(1,1,j)=(adv%phi(2,j)-adv%phi(1,j))/dr
          adv%grad_phi(1,nr+1,j)=(adv%phi(nr+1,j)-adv%phi(nr,j))/dr
          adv%grad_phi(2,1,j)=(adv%phi(1,modulo(j+1-1+ntheta,ntheta)+1)-adv%phi(1,modulo(j-1-1+ntheta,ntheta)+1))/(2*dtheta)
          adv%grad_phi(2,nr+1,j)=(adv%phi(nr+1,modulo(j+1-1+ntheta,ntheta)+1)-adv%phi(nr+1,modulo(j-1-1+ntheta,ntheta)+1))/(2*dtheta)
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

       call compute_spline_2D(adv%phi,adv%spl_phi)

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
             adv%grad_phi(2,i,j)=interpolate_x2_derivative_2D(r,theta,adv%spl_phi)
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


  !>subroutine advect_CG_polar(dt,dr,dtheta,nr,ntheta,rmin,rmax,phi,a,f,pfwd,pinv)
  !>compute step for Center-Guide equation
  !>dt : size of time step
  !>dr and dtheta : size of r- and theta-step
  !>nr and ntheta : number of steps in direction r and theta
  !>rmax : disc radius
  !>rmin : hole radius
  !>phi : field, size nrXntheta
  !>a : advection field, size 2XnrXntheta
  !>a = grad(phi), where phi is solution of laplacien(phi) = -f in Vlasov-Poisson equations
  !>f = distribution function, inout, size nrXntheta
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
    interpolate_case=3

    if (interpolate_case<3 .and. interpolate_case>0) then

       !construction of spline coefficients for f
       call compute_spline_2D(adv%f,adv%spl_f)

       do i=1,nr+1
          do j=1,ntheta+1
             r=rmin+real(i-1,f64)*dr
             theta=real(j-1,f64)*dtheta

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

             adv%f(i,j)=interpolate_value_2D(r,theta,adv%spl_f)

          end do
       end do

    else if (interpolate_case==3) then
       call rk4_polar_advect(adv,rk)

    else
       print*,'no way chosen to compute r'
       print*,'nothing will be done'
       print*,'see variable interpolate_case, routine advect_VP_polar in file prototype/src/simulation/sll_advection_polar.F90'
    end if

    adv%f(:,ntheta+1)=adv%f(:,1)

  end subroutine advect_CG_polar


  !>subroutine rk4_polar_advect(dt,dr,dtheta,nr,ntheta,rmin,rmax,f,phi,grad_phi,pfwd,pinv)
  !>RK4 for polar advection only
  !>dt, dr, dtheta : size of time step, step in direction r and theta
  !>nr, ntheta : number of steps in direction r and theta
  !>rmin and rmax : radius of the hole and of the disc
  !>f : distribution function, size (nr+1)Xntheta
  !>phi : field, size (nr+1)Xntheta
  !>grad_phi : gradient(phi), size (nr+1)Xntheta
  !>           grad_phi(1,i,j)=d_r(phi)(r_i,theta_j)
  !>           grad_phi(2,i,j)=d_theta(phi)(r_i,theta_j)
  subroutine rk4_polar_advect(adv,rk)

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

    !first step of RK4
    do i=1,nr+1
       r=rmin+real(i-1)*dr
       rk%r1(i,:)=-dt*adv%grad_phi(2,i,:)/r
       rk%theta1(i,:)=dt*adv%grad_phi(1,i,:)/r
    end do

    !2nd step of RK4
    do i=1,nr+1
       r=rmin+real(i-1)*dr
       do j=1,ntheta+1
          theta=real(j-1,f64)*dtheta
          rk%r4(i,j)=r-dt/2.0_f64*adv%grad_phi(2,i,j)/r
          rk%theta4(i,j)=theta+dt/2.0_f64*adv%grad_phi(1,i,j)/r

          call correction_r(rk%r4(i,j),rmin,rmax)
          call correction_theta(rk%theta4(i,j))

          adv%f(i,j)=interpolate_value_2D(rk%r4(i,j),rk%theta4(i,j),adv%spl_f)
       end do
    end do

    call poisson_solve_polar(adv)
    call compute_grad_field(adv)

    !construction of spline coeficients for a
    call compute_spline_2d(adv%grad_phi(1,:,:),adv%spl_a1)
    call compute_spline_2d(adv%grad_phi(2,:,:),adv%spl_a2)

    do i=1,nr+1
       rk%r4(i,:)=rmin+real(i-1,f64)*dr
       do j=1,ntheta+1
          r=rmin+real(i-1)*dr+rk%r1(i,j)/2.0_f64
          theta=real(j-1,f64)*dtheta+rk%theta1(i,j)/2.0_f64
          call correction_r(r,rmin,rmax)
          call correction_theta(theta)

          rk%r2(i,j)=-interpolate_value_2d(r,theta,adv%spl_a2)*dt/rk%r4(i,j)
          rk%theta2(i,j)=interpolate_value_2d(r,theta,adv%spl_a1)*dt/rk%r4(i,j)
       end do
    end do

    !3rd step of RK4
    do i=1,nr+1
       do j=1,ntheta+1
          r=rmin+real(i-1)*dr+rk%r2(i,j)/2.0_f64
          theta=real(j-1,f64)*dtheta+rk%theta2(i,j)/2.0_f64
          call correction_r(r,rmin,rmax)
          call correction_theta(theta) 

          rk%r3(i,j)=-interpolate_value_2d(r,theta,adv%spl_a2)*dt/rk%r4(i,j)
          rk%theta3(i,j)=interpolate_value_2d(r,theta,adv%spl_a1)*dt/rk%r4(i,j)
       end do
    end do

    !4th step of RK4
    do i=1,nr+1
       r=rmin+real(i-1)*dr
       do j=1,ntheta+1
          theta=real(j-i)*dtheta
          rk%r4(i,j)=r+rk%r1(i,j)
          rk%theta4(i,j)=theta+rk%theta1(i,j)
          call correction_r(rk%r4(i,j),rmin,rmax)
          call correction_theta(rk%theta4(i,j))

          adv%f(i,j)=interpolate_value_2D(rk%r4(i,j),rk%theta4(i,j),adv%spl_f)
       end do
    end do

    call poisson_solve_polar(adv)
    call compute_grad_field(adv)

    !construction of spline coeficients for a
    call compute_spline_2d(adv%grad_phi(1,:,:),adv%spl_a1)
    call compute_spline_2d(adv%grad_phi(2,:,:),adv%spl_a2)

    do i=1,nr+1
       do j=1,ntheta+1
          r=rmin+real(i-1)*dr+rk%r3(i,j)
          theta=real(j-1,f64)*dtheta+rk%theta3(i,j)
          call correction_r(r,rmin,rmax)
          call correction_theta(theta)

          rk%r4(i,j)=-interpolate_value_2d(r,theta,adv%spl_a2)*dt/r
          rk%theta4(i,j)=interpolate_value_2d(r,theta,adv%spl_a1)*dt/r
       end do
    end do

    !sommation
    do i=1,nr+1
       r=rmin+real(i-1,f64)*dr
       do j=1,ntheta+1
          theta=real(j-1,f64)*dtheta
          rk%r1(i,j)=r+rk%r1(i,j)/6.0_f64+rk%r2(i,j)/3.0_f64+rk%r3(i,j)/3.0_f64+rk%r4(i,j)/6.0_f64
          rk%theta1(i,j)=theta+rk%theta1(i,j)/6.0_f64+rk%theta2(i,j)/3.0_f64+rk%theta3(i,j)/3.0_f64+rk%theta4(i,j)/6.0_f64

          call correction_r(rk%r1(i,j),rmin,rmax)
          call correction_theta(rk%theta1(i,j))
       end do
    end do

    !updating the distribution function
    call deposit_value_2D(rk%r1,rk%theta1,adv%spl_f,adv%f)
!!$    do i=1,nr+1
!!$       do j=1,ntheta+1
!!$          adv%f(i,j)=interpolate_value_2d(rk%r1(i,j),rk%theta1(i,j),adv%spl_f)
!!$       end do
!!$    end do

  end subroutine rk4_polar_advect


  !>subroutine SL_classic(dt,dr,dtheta,nr,ntheta,rmin,rmax,pfwd,pinv,f,phi,grad_phi)
  !>computes the classic semi-Lagrangian scheme for Vlasov-Poisson equation
  !>dt, dr, dtheta : size of step in direction t (time), r and theta
  !>nr and ntheta : number ox step in direction r and theta
  !>                number of points are nr+1 and ntheta+1
  !>pfwd and pinv : initialized object for FFT forward and backward
  !>f : distribution function, size (nr+1)X(ntheta+1)
  !>phi : field, size (nr+1)X(ntheta+1)
  !>grad_phi : gradient(phi), size 2X(nr+1)X(ntheta+1)
  !>           grad_phi(1,i,j)=d_r(phi)(r_i,theta_j)
  !>           grad_phi(2,i,j)=d_theta(phi)(r_i,theta_j)
  subroutine SL_classic(adv,rk)

    implicit none

    type(polar_vp_data), intent(inout), pointer :: adv
    type(polar_vp_rk4), intent(inout), pointer :: rk

    call poisson_solve_polar(adv)
    call compute_grad_field(adv)
    call advect_CG_polar(adv,rk)

  end subroutine SL_classic


  !>subroutine SL_controlled(dt,dr,dtheta,nr,ntheta,rmin,rmax,pfwd,pinv,f,fdemi,phi,grad_phi)
  !>computes the semi-Lagrangian scheme with control for Vlasov-Poisson equation
  !>dt, dr, dtheta : size of step in direction t (time), r and theta
  !>nr and ntheta : number ox step in direction r and theta
  !>                number of points are nr+1 and ntheta+1
  !>pfwd and pinv : initialized object for FFT forward and backward
  !>f : distribution function, size (nr+1)X(ntheta+1)
  !>fdemi : distribution function at time t^(n+1/2)
  !>phi : field, size (nr+1)X(ntheta+1)
  !>grad_phi : gradient(phi), size 2X(nr+1)X(ntheta+1)
  !>           grad_phi(1,i,j)=d_r(phi)(r_i,theta_j)
  !>           grad_phi(2,i,j)=d_theta(phi)(r_i,theta_j)
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

    do while (theta>2.0_f64*sll_pi)
       theta=theta-2.0_f64*sll_pi
    end do
    do while (theta<0.0_f64)
       theta=theta+2.0_f64*sll_pi
    end do

  end subroutine correction_theta

end module polar_advection
