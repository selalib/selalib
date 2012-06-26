module polar_advection
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"

  use poisson_polar
  use numeric_constants
  use sll_splines
  implicit none

contains

  !compute a = grad(phi) for phi scalar field
  subroutine compute_advection(nr,ntheta,dr,dtheta,rmin,rmax,phi,a)
    ! nr and ntheta : numbers of steps in directions r and theta
    ! dr and dtheta : size of step in directions r and theta
    ! rmin : radius of the hole
    ! rmax : radius of the disc
    ! phi : known field, size (nr+1)X(ntheta+1)
    ! a = grad(phi), size 2X(nr+1)X(ntheta+1)
    ! a(1,i,j) = d_r(phi)(r_i,theta_j)
    ! a(2,i,j) = d_theta(phi)(r_i,theta_j)

    implicit none

    sll_int32, intent(in) :: nr, ntheta
    sll_real64, intent(in) :: dr, dtheta, rmin, rmax
    sll_real64, dimension(:,:), intent(in), pointer :: phi
    sll_real64, dimension(:,:,:), intent(out), pointer :: a

    sll_int32 :: i,j,calculus
    sll_real64 :: r,theta,x,y

    ! way to compute a
    ! 1 : center for r and theta
    !     decenter on boundaries
    calculus = 1

    if (calculus==1) then
       ! center formula for r end theta
       ! decenter for r on boundaries
       
       do i=2,nr
          r=rmin+real(i-1,f64)*dr
          do j=1,ntheta
             a(1,i,j)=(phi(i+1,j)-phi(i-1,j))/(2*dr)
             a(2,i,j)=(phi(i,modulo(j+1-1+ntheta,ntheta)+1)-phi(i,modulo(j-1-1+ntheta,ntheta)+1))/(2*dtheta)
          end do
       end do
       do j=1,ntheta
          a(1,1,j)=(phi(2,j)-phi(1,j))/(2*dr)
          a(1,nr+1,j)=(phi(nr,j)-phi(nr-1,j))/(2*dr)
          a(2,1,j)=(phi(1,modulo(j+1-1+ntheta,ntheta)+1)-phi(1,modulo(j-1-1+ntheta,ntheta)+1))/(2*dtheta)
          a(2,nr+1,j)=(phi(nr,modulo(j+1-1+ntheta,ntheta)+1)-phi(nr,modulo(j-1-1+ntheta,ntheta)+1))/(2*dtheta)
       end do
    end if

    a(:,:,ntheta+1)=a(:,:,1)

    open(unit=22,file='phi.dat')
    do i=1,nr+1
       r=rmin+real(i-1,f64)*dr
       do j=1,ntheta+1
          theta=real(j-1,f64)*2.0_f64*sll_pi
          x=r*cos(theta)
          y=r*sin(theta)
          write (22,*)x,y,a(1,i,j),a(2,i,j)
       end do
    end do
    close(22)

  end subroutine compute_advection


  !compute step for Vlasov-Poisson equation
  subroutine advect_VP_polar(dt,dr,dtheta,nr,ntheta,rmin,rmax,phi,a,f,pfwd,pinv)
    ! dt : size of time step
    ! dr and dtheta : size of r- and theta-step
    ! nr and ntheta : number of steps in direction r and theta
    ! rmax : disc radius
    ! rmin : hole radius
    ! phi : field, size nrXntheta
    ! a : advection field, size 2XnrXntheta
    ! a = grad(phi), where phi is solution of laplacien(phi) = -f in Vlasov-Poisson equations
    ! f = distribution function, inout, size nrXntheta

    implicit none

    sll_int32, intent(in) :: nr, ntheta
    sll_real64, intent(in) :: dt, dr, dtheta, rmin, rmax
    sll_real64, dimension(:,:,:), intent(in),pointer :: a
    sll_real64, dimension(:,:), intent(in), pointer :: phi
    sll_real64, dimension(:,:), intent(inout), pointer :: f
    type(sll_fft_plan), intent(in), pointer ::pfwd, pinv

    type(sll_spline_2D), pointer :: spl_f
    sll_int32 :: interpolate_case
    sll_int32 :: i,j
    sll_real64 :: r, theta

    ! creation spline 
    spl_f => new_spline_2D(Nr+1, Ntheta+1, &
         rmin, rmax, &
         0.0_f64, 2.0_f64*sll_pi, &
         HERMITE_SPLINE, PERIODIC_SPLINE,&
         const_slope_x1_min = 0.0_f64,const_slope_x1_max = 0.0_f64)

    !interpolation
    ! 1 : using explicite Eulerian scheme
    ! 2 : rotation, this case ignore the field phi
    !     roation speed = -1
    ! 3 : using RK4
    interpolate_case=1

    if (interpolate_case/=3) then

       !construction of spline coeficients for f
       call compute_spline_2D(f,spl_f)

       do i=1,nr+1
          do j=1,ntheta+1
             r=rmin+real(i-1,f64)*dr
             theta=real(j-1,f64)*dtheta

             if (interpolate_case==1) then
                !Eulerian scheme
                r=r-dt*a(2,i,j)/r
                theta=theta+dt*a(1,i,j)
             else if (interpolate_case==2) then
                !rotation
                theta=theta+dt
             end if

             !correction of r
             if (r>rmax) then
                r=rmax
             else if (r<rmin) then
                r=rmin
             end if
             !correction of theta
             if (theta<0.0_f64) then
                theta=theta+2*sll_pi
             else if (theta>2*sll_pi) then
                theta=theta-2*sll_pi
             end if

             f(i,j)=interpolate_value_2D(r,theta,spl_f)

          end do
       end do

    else if (interpolate_case==3) then
       call rk4_polar_advect(dt,dr,dtheta,nr,ntheta,rmin,rmax,f,phi,a,pfwd,pinv)
    end if

    call delete_spline_2d(spl_f)

  end subroutine advect_VP_polar


  !KR4 for polar advection only
  subroutine rk4_polar_advect(dt,dr,dtheta,nr,ntheta,rmin,rmax,f,phi,grad_phi,pfwd,pinv)
    ! dt, dr, dtheta : size of time step, step in direction r and theta
    ! nr, ntheta : number of steps in direction r and theta
    ! rmin and rmax : radius of the hole and of the disc
    ! f : distribution function, size (nr+1)Xntheta
    ! phi : field, size (nr+1)Xntheta
    ! grad_phi : gradient(phi), size (nr+1)Xntheta
    !            grad_phi(1,i,j)=d_r(phi)(r_i,theta_j)
    !            grad_phi(2,i,j)=d_theta(phi)(r_i,theta_j)

    implicit none

    sll_real64, intent(in) :: dt, dr, dtheta, rmin, rmax
    sll_int32, intent(in) :: nr, ntheta
    sll_real64, dimension(:,:), intent(in), pointer :: phi
    sll_real64, dimension(:,:,:), intent(in), pointer :: grad_phi
    sll_real64, dimension(:,:), intent(inout), pointer :: f
    type(sll_fft_plan), intent(in), pointer ::pfwd, pinv

    sll_int32 :: i,j,err
    sll_real64 :: r,theta
    sll_real64, dimension(:,:), pointer :: r1,r2,r3,r4, theta1,theta2,theta3,theta4
    sll_real64, dimension(:,:), pointer :: fcopie, phicopie
    sll_real64, dimension(:,:,:), pointer :: a
    type(sll_spline_2D), pointer :: spl_f, spl_a1, spl_a2

    SLL_ALLOCATE(fcopie(nr+1,ntheta+1),err)
    SLL_ALLOCATE(phicopie(nr+1,ntheta+1),err)
    SLL_ALLOCATE(r1(nr+1,ntheta+1),err)
    SLL_ALLOCATE(r2(nr+1,ntheta+1),err)
    SLL_ALLOCATE(r3(nr+1,ntheta+1),err)
    SLL_ALLOCATE(r4(nr+1,ntheta+1),err)
    SLL_ALLOCATE(theta1(nr+1,ntheta+1),err)
    SLL_ALLOCATE(theta2(nr+1,ntheta+1),err)
    SLL_ALLOCATE(theta3(nr+1,ntheta+1),err)
    SLL_ALLOCATE(theta4(nr+1,ntheta+1),err)
    SLL_ALLOCATE(a(2,nr+1,ntheta+1),err)

    ! creation spline 
    spl_f => new_spline_2D(nr+1,ntheta+1,rmin,rmax,0._f64, 2._f64*sll_pi, &
         HERMITE_SPLINE, PERIODIC_SPLINE,&
         const_slope_x1_min = 0._f64,const_slope_x1_max = 0._f64)
    spl_a1 => new_spline_2D(nr+1,ntheta+1,rmin,rmax,0._f64, 2._f64*sll_pi, &
         HERMITE_SPLINE, PERIODIC_SPLINE,&
         const_slope_x1_min = 0._f64,const_slope_x1_max = 0._f64)
    spl_a2 => new_spline_2D(nr+1,ntheta+1,rmin,rmax,0._f64, 2._f64*sll_pi, &
         HERMITE_SPLINE, PERIODIC_SPLINE,&
         const_slope_x1_min = 0._f64,const_slope_x1_max = 0._f64)

    fcopie=f
    phicopie=phi
    a=grad_phi

    !construction of spline coeficients for f
    call compute_spline_2D(f,spl_f)

    !first step of RK4
    do i=1,nr+1
       r=rmin+real(i-1)*dr
       r1(i,:)=dt*a(2,i,:)/r
    end do
    theta1(:,:)=-dt*a(1,:,:)

    !2nd step of RK4
    do i=1,nr+1
       r=rmin+real(i-1)*dr
       do j=1,ntheta+1
          theta=real(j-1,f64)*dtheta
          r3(i,j)=r-dt/2.0_f64*a(2,i,j)/r
          theta3(i,j)=theta+dt/2.0_f64*a(1,i,j)

          !correction of r
          if (r3(i,j)>rmax) then
             r3(i,j)=rmax
          else if (r3(i,j)<rmin) then
             r3(i,j)=rmin
          end if
          !correction of theta
          if (theta3(i,j)>2.0_f64*sll_pi) then
             theta3(i,j)=theta3(i,j)-2.0_f64*sll_pi
          else if (theta3(i,j)<0.0_f64) then
             theta3(i,j)=theta3(i,j)+2.0_f64*sll_pi
          end if

          fcopie(i,j)=interpolate_value_2D(r3(i,j),theta3(i,j),spl_f)
       end do
    end do

    call poisson_solve_polar(fcopie,rmin,dr,nr,ntheta-1,pfwd,pinv,phicopie)
    call compute_advection(nr,ntheta,dr,dtheta,rmin,rmax,phicopie,a)

    !construction of spline coeficients for a
    call compute_spline_2d(a(1,:,:),spl_a1)
    call compute_spline_2d(a(2,:,:),spl_a2)

    do i=1,nr+1
       do j=1,ntheta+1
          r2(i,j)=interpolate_value_2d(r3(i,j)+r1(i,j)/2.0_f64,theta3(i,j)+theta1(i,j)/2.0_f64,spl_a2)* &
               & dt/(r3(i,j)+r1(i,j)/2.0_f64)
          theta2(i,j)=-interpolate_value_2d(r3(i,j)+r1(i,j)/2.0_f64,theta3(i,j)+theta1(i,j)/2.0_f64,spl_a2)*dt
       end do
    end do

    !3rd step of RK4
    !fcopie, phicopie and a already at t+dt/2
    !r and theta already at t+dt/2
    do i=1,nr+1
       do j=1,ntheta+1
          r3(i,j)=interpolate_value_2d(r3(i,j)+r2(i,j)/2.0_f64,theta3(i,j)+theta2(i,j)/2.0_f64,spl_a2)* &
               & dt/(r3(i,j)+r2(i,j)/2.0_f64)
          theta3(i,j)=-interpolate_value_2d(r3(i,j)+r2(i,j)/2.0_f64,theta3(i,j)+theta2(i,j)/2.0_f64,spl_a2)*dt
       end do
    end do

    !4th step of RK4
    fcopie=f
    phicopie=phi
    a=grad_phi
    do i=1,nr+1
       r=rmin+real(i-1)*dr
       do j=1,ntheta+1
          r4(i,j)=r-dt*a(2,i,j)/r
          theta4(i,j)=theta+dt*a(1,i,j)

          !correction of r
          if (r4(i,j)>rmax) then
             r4(i,j)=rmax
          else if (r4(i,j)<rmin) then
             r4(i,j)=rmin
          end if
          !correction of theta
          if (theta4(i,j)>2.0_f64*sll_pi) then
             theta4(i,j)=theta4(i,j)-2.0_f64*sll_pi
          else if (theta4(i,j)<0.0_f64) then
             theta4(i,j)=theta4(i,j)+2.0_f64*sll_pi
          end if

          fcopie(i,j)=interpolate_value_2D(r4(i,j),theta4(i,j),spl_f)
       end do
    end do

    call poisson_solve_polar(fcopie,rmin,dr,nr,ntheta-1,pfwd,pinv,phicopie)
    call compute_advection(nr,ntheta,dr,dtheta,rmin,rmax,phicopie,a)

    !construction of spline coeficients for a
    call compute_spline_2d(a(1,:,:),spl_a1)
    call compute_spline_2d(a(2,:,:),spl_a2)

    do i=1,nr+1
       do j=1,ntheta+1
          r4(i,j)=interpolate_value_2d(r4(i,j)+r3(i,j),theta4(i,j)+theta3(i,j),spl_a2)* &
               & dt/(r4(i,j)+r3(i,j))
          theta3(i,j)=-interpolate_value_2d(r4(i,j)+r3(i,j),theta4(i,j)+theta3(i,j),spl_a2)*dt
       end do
    end do

    !sommation
    f=f+r1/6.0_f64+r2/3.0_f64+r3/3.0_f64+r4/6.0_f64

    call delete_spline_2d(spl_f)
    call delete_spline_2d(spl_a1)
    call delete_spline_2d(spl_a2)
    SLL_DEALLOCATE(fcopie,err)
    SLL_DEALLOCATE(phicopie,err)
    SLL_DEALLOCATE(a,err)
    SLL_DEALLOCATE(r1,err)
    SLL_DEALLOCATE(r2,err)
    SLL_DEALLOCATE(r3,err)
    SLL_DEALLOCATE(r4,err)
    SLL_DEALLOCATE(theta1,err)
    SLL_DEALLOCATE(theta2,err)
    SLL_DEALLOCATE(theta3,err)
    SLL_DEALLOCATE(theta4,err)

  end subroutine rk4_polar_advect


  !compute the classic semi-Lagrangian scheme
  subroutine SL_classic(dt,dr,dtheta,nr,ntheta,rmin,rmax,pfwd,pinv,f,phi,grad_phi)
    ! dt, dr, dtheta : size of step in direction t (time), r and theta
    ! nr and ntheta : number ox step in direction r and theta
    !                 number of points are nr+1 and ntheta+1
    ! pfwd and pinv : initialized object for FFT forward and backward
    ! f : distribution function, size (nr+1)X(ntheta+1)
    ! phi : field, size (nr+1)X(ntheta+1)
    ! grad_phi : gradient(phi), size 2X(nr+1)X(ntheta+1)
    !            grad_phi(1,i,j)=d_r(phi)(r_i,theta_j)
    !            grad_phi(2,i,j)=d_theta(phi)(r_i,theta_j)

    implicit none

    sll_real64, intent(in) :: dt,dr,dtheta,rmin,rmax
    sll_int32, intent(in) :: nr, ntheta
    type(sll_fft_plan), intent(inout), pointer ::pfwd, pinv
    sll_real64, dimension(:,:), intent(inout), pointer :: f, phi
    sll_real64, dimension(:,:,:), intent(inout), pointer :: grad_phi

    call poisson_solve_polar(f,rmin,dr,nr,ntheta,pfwd,pinv,phi)
    call compute_advection(nr,ntheta,dr,dtheta,rmin,rmax,phi,grad_phi)
    call advect_VP_polar(dt,dr,dtheta,nr,ntheta,rmin,rmax,phi,grad_phi,f,pfwd,pinv)
  end subroutine SL_classic


  !compute the semi-Lagrangian scheme with control
  subroutine SL_controlled(dt,dr,dtheta,nr,ntheta,rmin,rmax,pfwd,pinv,f,fdemi,phi,grad_phi)
    ! dt, dr, dtheta : size of step in direction t (time), r and theta
    ! nr and ntheta : number ox step in direction r and theta
    !                 number of points are nr+1 and ntheta+1
    ! pfwd and pinv : initialized object for FFT forward and backward
    ! f : distribution function, size (nr+1)X(ntheta+1)
    ! fdemi : distribution function at time t^(n+1/2)
    ! phi : field, size (nr+1)X(ntheta+1)
    ! grad_phi : gradient(phi), size 2X(nr+1)X(ntheta+1)
    !            grad_phi(1,i,j)=d_r(phi)(r_i,theta_j)
    !            grad_phi(2,i,j)=d_theta(phi)(r_i,theta_j)

    implicit none

    sll_real64, intent(in) :: dt,dr,dtheta,rmin,rmax
    sll_int32 :: nr, ntheta
    type(sll_fft_plan), intent(inout), pointer ::pfwd, pinv
    sll_real64, dimension(:,:), intent(inout), pointer :: f, fdemi, phi
    sll_real64, dimension(:,:,:), intent(inout), pointer :: grad_phi

    fdemi=f
    call poisson_solve_polar(fdemi,rmin,dr,nr,ntheta,pfwd,pinv,phi)
    call compute_advection(nr,ntheta,dr,dtheta,rmin,rmax,phi,grad_phi)
    call advect_VP_polar(dt/2.0_f64,dr,dtheta,nr,ntheta,rmin,rmax,phi,grad_phi,fdemi,pfwd,pinv)
    !we just obtained f^(n+1/2)
    call poisson_solve_polar(fdemi,rmin,dr,nr,ntheta,pfwd,pinv,phi)
    !we just obtained E^(n+1/2)
    call compute_advection(nr,ntheta,dr,dtheta,rmin,rmax,phi,grad_phi)
    call advect_VP_polar(dt,dr,dtheta,nr,ntheta,rmin,rmax,phi,grad_phi,f,pfwd,pinv)

  end subroutine SL_controlled

end module polar_advection
