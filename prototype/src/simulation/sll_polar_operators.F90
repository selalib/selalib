module polar_operators
! definition of gradient and divergence in polar coordinate
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"

  use sll_poisson_2d_polar
  use sll_fft
  use sll_cubic_splines
  implicit none

  !>type plan_polar_op
  !>data type for polar gradient and polar divergence
  type plan_polar_op
     sll_real64 :: rmin,rmax,dr,dtheta
     sll_int32 :: nr,ntheta
     sll_int32 :: grad_case
     type(sll_cubic_spline_2D), pointer :: spl_phi
     type(sll_fft_plan), pointer :: pfwd,pinv
     sll_comp64, dimension(:,:), pointer :: grad_fft
  end type plan_polar_op

contains

!===================================
!  construction of plan_polar_op
!===================================

  !>function new_polar_op(rmin,rmax,dr,dtheta,nr,ntheta,grad_case)
  !>creation of plan_polar_op
  !>rmin : interior radius
  !>rmax : exterior radius
  !>dr and dtheta : size of step in direction r and theta
  !>nr and ntheta : number of step in direction r and theta
  !>grad_case : integer, enable to choose the way to compute the gradient, optional
  !>            1 : finit differences center for r and theta, decenter on boundaries
  !>            2 : finit differences center for r, decenter on boundaries, using fft for theta
  !>                this case is not possible with fftpack
  !>            3 : using splines
  !>            default case is 3
  function new_polar_op(rmin,rmax,dr,dtheta,nr,ntheta,grad_case) result(this)

    type(plan_polar_op), pointer :: this
    sll_real64, intent(in) :: rmin,rmax,dr,dtheta
    sll_int32, intent(in) :: nr, ntheta
    sll_int32, intent(in), optional :: grad_case
    sll_int32 :: err
    sll_real64, dimension(:), allocatable :: bufr
    sll_comp64, dimension(:), allocatable :: bufc

    SLL_ALLOCATE(this,err)

    this%rmin=rmin
    this%rmax=rmax
    this%dr=dr
    this%dtheta=dtheta
    this%nr=nr
    this%ntheta=ntheta

    if (.not. present(grad_case)) then
       this%grad_case=3
    else
       this%grad_case=grad_case
    end if

    if (grad_case==2) then
       SLL_ALLOCATE(this%grad_fft(nr+1,ntheta/2+1),err)
    end if

    this%spl_phi => new_spline_2D(nr+1,ntheta+1,rmin,rmax,0._f64, 2._f64*sll_pi, &
         & SLL_HERMITE, SLL_PERIODIC,const_slope_x1_min = 0._f64,const_slope_x1_max = 0._f64)

    SLL_ALLOCATE(bufr(ntheta),err)
    SLL_ALLOCATE(bufc(ntheta/2+1),err)

    this%pfwd => fft_new_plan(ntheta,bufr,bufc,FFT_NORMALIZE)
    this%pinv => fft_new_plan(ntheta,bufc,bufr)
    SLL_DEALLOCATE_ARRAY(bufr,err)
    SLL_DEALLOCATE_ARRAY(bufc,err)

  end function new_polar_op

!===============================
!  deletion of plan_polar_op
!===============================

  !>subroutine delete_plan_polar_op(this)
  !>deletion of plan_polar_op object
  subroutine delete_plan_polar_op(this)

    implicit none

    type(plan_polar_op), pointer :: this
    sll_int32 :: err

    if (this%grad_case==2) then
       SLL_DEALLOCATE_ARRAY(this%grad_fft,err)
    end if
    call fft_delete_plan(this%pfwd)
    call fft_delete_plan(this%pinv)
    call sll_delete(this%spl_phi)
    SLL_DEALLOCATE(this,err)
    this=>null()

  end subroutine delete_plan_polar_op

!===========================================
!  computation of gradient and divergence
!===========================================

  !>subroutine compute_grad_field(plan,phi,grad_phi)
  !>compute grad(phi) for phi scalar field in cartesian coordinate.
  !>For computation in polar coordinate, juste divide the second coordinate by r
  !>plan : plan_polar_op object, data for grad and div in polar
  !>phi : scalar field, size (nr+1)*(ntheta+1), input
  !>grad_phi : grad(phi), size 2*(nr+1)*(ntheta+1), output
  !>           grad_phi(1,:,:)=d_r(phi), grad_phi(2,:,:)=d_theta(phi)/r
  subroutine compute_grad_field(plan,phi,grad_phi)

    implicit none

    type(plan_polar_op), pointer              :: plan
    sll_real64, dimension(:,:), intent(inout) :: phi
    !phi is inout because when I wrote the code it needed to be inout in fft_apply_plan
    sll_real64, dimension(:,:,:), intent(out) :: grad_phi

    sll_int32 :: nr, ntheta
    sll_real64 :: dr, dtheta, rmin, rmax
    sll_int32 :: i,j
    sll_real64 :: r,theta,tmp
    sll_comp64 :: temp

    nr=plan%nr
    ntheta=plan%ntheta
    dr=plan%dr
    dtheta=plan%dtheta
    rmin=plan%rmin
    rmax=plan%rmax
    
    !experimental
    tmp=0._f64
    do j=1,ntheta
      tmp=tmp+phi(1,j)
    enddo
    tmp=tmp/real(ntheta,f64)
    phi(1,:)=tmp

    tmp=0._f64
    do j=1,ntheta
      tmp=tmp+phi(nr+1,j)
    enddo
    tmp=tmp/real(ntheta,f64)
    phi(nr+1,:)=tmp

    
    
    
    if (plan%grad_case==1) then
       ! center formula for r end theta
       ! decenter for r on boundaries
       do i=2,nr
          do j=1,ntheta
             grad_phi(1,i,j)=(phi(i+1,j)-phi(i-1,j))/(2*dr)
             grad_phi(2,i,j)=(phi(i,modulo(j+1-1+ntheta,ntheta)+1)-phi(i,modulo(j-1-1+ntheta,ntheta)+1))/(2*dtheta)
          end do
       end do
       do j=1,ntheta
!!$          grad_phi(1,1,j)=phi(2,j)-0.0_f64)/(2.0_f64*dr)
!!$          grad_phi(1,nr+1,j)=(0.0_f64-phi(nr,j))/(dr*2.0_f64)
!!$          grad_phi(1,1,j)=0.0_f64
!!$          grad_phi(1,nr+1,j)=0.0_f64
          grad_phi(1,1,j)=(phi(2,j)-phi(1,j))/dr
          grad_phi(1,nr+1,j)=(phi(nr+1,j)-phi(nr,j))/dr
!!$          grad_phi(1,1,j)=-(1.5_f64*phi(1,j)-2.0_f64*phi(2,j)+0.5_f64*phi(3,j))/dr
!!$          grad_phi(1,nr+1,j)=-(1.5_f64*phi(nr+1,j)-2.0_f64*phi(nr,j)+0.5_f64*phi(nr-1,j))/dr
          grad_phi(2,1,j)=(phi(1,modulo(j+1-1+ntheta,ntheta)+1)-phi(1,modulo(j-1-1+ntheta,ntheta)+1))/(2*dtheta)
          grad_phi(2,nr+1,j)=(phi(nr+1,modulo(j+1-1+ntheta,ntheta)+1)-phi(nr+1,modulo(j-1-1+ntheta,ntheta)+1))/(2*dtheta)
       end do

    else if (plan%grad_case==2) then
       ! center formula for r, decenter on boundaries
       ! using fft for theta
       !not done

       do i=2,nr
          r=rmin+real(i-1,f64)*dr
          grad_phi(1,i,:)=(phi(i+1,:)-phi(i-1,:))/(2*dr)
       end do
       do j=1,ntheta
          grad_phi(1,1,j)=(phi(2,j)-phi(1,j))/dr
          grad_phi(1,nr+1,j)=(phi(nr+1,j)-phi(nr,j))/dr
       end do

       do i=1,nr+1
          !call fft_apply_plan_r2c_1d(plan%pfwd,phi(i,1:ntheta),plan%grad_fft(i,:))
          call fft_apply_plan(plan%pfwd,phi(i,1:ntheta),plan%grad_fft(i,:))
       end do

       do i=1,nr+1
          do j=0,ntheta/2
             temp=fft_get_mode(plan%pfwd,plan%grad_fft(i,:),j)
             temp=temp*cmplx(0.0_f64,real(j,f64),kind=f64)
             call fft_set_mode(plan%pinv,plan%grad_fft(i,:),temp,j)
          end do
       end do

       do i=1,nr+1
          call fft_apply_plan(plan%pinv,plan%grad_fft(i,:),grad_phi(2,i,1:ntheta))
          grad_phi(2,i,ntheta+1)=grad_phi(2,i,1)
       end do

    else if (plan%grad_case==3) then
       ! using splines for r and theta

       call compute_spline_2D(phi,plan%spl_phi)

       do j=1,ntheta
          theta=real(j-1,f64)*dtheta
          do i=1,nr+1
             r=rmin+real(i-1,f64)*dr
             grad_phi(1,i,j)=interpolate_x1_derivative_2D(r,theta,plan%spl_phi)
             grad_phi(2,i,j)=interpolate_x2_derivative_2D(r,theta,plan%spl_phi)
          end do
       end do

    else
       print*,'no choosen way to compute grad'
       print*,'initialization missing'
    end if

    grad_phi(:,:,ntheta+1)=grad_phi(:,:,1)

    !print *,sum(abs(phi(2,:)-phi(1,:)))
    !print *,'#',sum(abs(grad_phi(1,1,:))),sum(abs(phi(1,:)))
  end subroutine compute_grad_field


  !>subroutine divergence_scalar_field(plan,field,div)
  !>compute divergence of field in polar coordinate
  !>plan : plan_polar_op object
  !>field : size 2*(nr+1)*(ntheta+1), input
  !>div : polar divergence of field at point (i,j), size (nr+1)*(ntheta+1) output
  subroutine divergence_scalar_field(plan,field,div)

    implicit none

    type(plan_polar_op), pointer             :: plan
    sll_real64, dimension(:,:,:), intent(in) :: field
    sll_real64, dimension(:,:), intent(out)  :: div

    sll_real64 :: dr,dtheta,rmin,rmax
    sll_int32 :: nr,ntheta
    sll_real64 :: r
    sll_int32 :: i,j

    nr=plan%nr
    ntheta=plan%ntheta
    dr=plan%dr
    dtheta=plan%dtheta
    rmin=plan%rmin
    rmax=plan%rmax

    div=0.0_f64

    do i=2,nr
       r=rmin+real(i-1,f64)*dr
       do j=1,ntheta
          div(i,j)=1/r*((field(1,i+1,j)*(r+dr)-field(1,i-1,j)*(r-dr))/(2*dr) &
               & +(field(2,i,modulo(j+1-1+ntheta,ntheta+1)+1)-field(2,i,modulo(j-1-1+ntheta,ntheta+1)+1))/(2*dtheta))
       end do
    end do
    do j=1,ntheta
       div(1,j)=1/rmin*((field(1,2,j)*(rmin+dr)-field(1,i,j)*(rmin))/(dr) &
            & +(field(2,1,modulo(j+1-1+ntheta,ntheta+1)+1)-field(2,1,modulo(j-1-1+ntheta,ntheta+1)+1))/(2*dtheta))
       div(nr+1,j)=1/r*((field(1,nr+1,j)*(rmax)-field(1,nr+1-1,j)*(rmax-dr))/(2*dr) &
            & +(field(2,nr+1,modulo(j+1-1+ntheta,ntheta+1)+1)-field(2,nr+1,modulo(j-1-1+ntheta,ntheta+1)+1))/(2*dtheta))
    end do

    div(:,ntheta+1)=div(:,1)

  end subroutine divergence_scalar_field

end module polar_operators
