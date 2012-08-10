module polar_kind
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"

  use sll_fft
  use sll_splines
  use numeric_constants
  implicit none

  !>type polar_data
  !>generic type for problems in polar
  !>contains size of time and space steps, boundaries, number of steps in time and space,
  !>and the final time
  type polar_data
     sll_real64 :: dt, dr, dtheta
     sll_real64 :: tf,rmin, rmax
     sll_int32 :: nb_step,nr, ntheta
  end type polar_data

  !>type polar_VP_data
  !>contains most needed data for the resolution of Vlasov-Poisson center-line equations
  !>contains the distribution function, the field and its grad, //...
  type polar_VP_data
     type(polar_data), pointer :: data
     type(sll_fft_plan), pointer :: pfwd, pinv
     sll_real64, dimension(:), allocatable :: rr,ttheta
     sll_real64, dimension(:,:), allocatable :: f,phi,f_fft,fdemi
     sll_real64, dimension(:,:,:), allocatable :: grad_phi
     sll_comp64, dimension(:), allocatable :: fk,phik
     type(sll_spline_2D), pointer :: spl_f, spl_a1, spl_a2, spl_phi
     !for the tridiagonal solver
     sll_real64, dimension(:), allocatable :: cts,a
     sll_int32, dimension(:), allocatable :: ipiv
  end type polar_VP_data

  !>type polar_VP_rk4
  !>used for RK4 in Vlasov Poisson
  !>contains r1, r2, r3, r4, theta1, theta2, theta3, theta4 for all points
  type polar_VP_rk4
     sll_real64, dimension(:,:), allocatable :: r1,r2,r3,r4
     sll_real64, dimension(:,:), allocatable :: theta1, theta2, theta3, theta4
  end type polar_VP_rk4

!   //==============\\
!   ||  INTERFACES  ||
!   \\==============//

  !>new_polar_data
  !>build a polar_data object
  !>can be build knowing only two of the folowing data : final time(tf), number of step in time (nb_step),
  !>and size of time step (dt)
  !>Syntaxe to use the interface :
  !>1/ knowing tf and dt : new_polar_data(dt,tf,rmin,rmax,nr,ntheta)
  !>2/ knowing dt and nb_step : new_polar_data(nb_step,dt,rmin,rmax,nr,ntheta)
  !>3/ knowing tf and nb_step : new_polar_data(tf,rmin,rmax,nb_step,nr,ntheta)
  interface new_polar_data
     module procedure new_polar_data_dt_tf, new_polar_data_tf_nbstep, new_polar_data_dt_nbstep
  end interface new_polar_data

  !>new_VP_data
  !>build a polar_VP_data object
  !>can be build with a polar_data object or with the list of datas
  !>Syntaxe to use the interface :
  !>1/ with a polar_data object : new_VP_data(polar_data)
  !>2/ without a polar_data : see "new_polar_data" section, same syntaxe with new_VP_data
  interface new_VP_data
     module procedure new_VP_dat_from_polar_data, new_VP_dat_from_all_dt_nbstep, new_VP_dat_from_all_dt_tf, &
          & new_VP_dat_from_all_tf_nbstep
  end interface new_VP_data

contains

!=========================================
!  beginnig of creation of polar data
!=========================================

  function new_polar_data_dt_nbstep(nb_step,dt,rmin,rmax,nr,ntheta) result(this)

    implicit none

    type(polar_data), pointer :: this
    sll_real64, intent(in) :: dt
    sll_real64, intent(in) :: rmin, rmax
    sll_int32, intent(in) :: nb_step,nr, ntheta
    sll_int32 :: err

    SLL_ALLOCATE(this,err)
    this%dt=dt
    this%dr=(rmax-rmin)/real(nr,f64)
    this%dtheta=2.0_f64*sll_pi/real(ntheta,f64)
    this%rmin=rmin
    this%rmax=rmax
    this%nb_step=nb_step
    this%nr=nr
    this%ntheta=ntheta
    this%tf=dt*real(nb_step,f64)

  end function new_polar_data_dt_nbstep


  function new_polar_data_dt_tf(dt,tf,rmin,rmax,nr,ntheta) result(this)

    implicit none

    type(polar_data), pointer :: this
    sll_real64, intent(in) :: dt
    sll_real64, intent(in) :: tf,rmin, rmax
    sll_int32, intent(in) :: nr, ntheta
    sll_int32 :: err

    SLL_ALLOCATE(this,err)
    this%dt=dt
    this%dr=(rmax-rmin)/real(nr,f64)
    this%dtheta=2.0_f64*sll_pi/real(ntheta,f64)
    this%rmin=rmin
    this%rmax=rmax
    this%nb_step=ceiling(tf/dt)
    this%nr=nr
    this%ntheta=ntheta
    this%tf=tf

  end function new_polar_data_dt_tf


  function new_polar_data_tf_nbstep(tf,rmin,rmax,nb_step,nr,ntheta) result(this)

    implicit none

    type(polar_data), pointer :: this
    sll_real64, intent(in) :: tf,rmin, rmax
    sll_int32, intent(in) :: nb_step,nr, ntheta
    sll_int32 :: err

    SLL_ALLOCATE(this,err)
    this%dt=tf/real(nb_step,f64)
    this%dr=(rmax-rmin)/real(nr,f64)
    this%dtheta=2.0_f64*sll_pi/real(ntheta,f64)
    this%rmin=rmin
    this%rmax=rmax
    this%nb_step=nb_step
    this%nr=nr
    this%ntheta=ntheta
    this%tf=tf

  end function new_polar_data_tf_nbstep

!============================================
!  beginning of creation of polar_vp_data
!============================================

  function new_VP_dat_from_polar_data(data_pol) result(this)

    implicit none

    type(polar_data) :: data_pol
    type(polar_VP_data), pointer :: this

    sll_int32 :: err,i
    sll_real64,dimension(:),pointer :: buf

    SLL_ALLOCATE(this,err)
    SLL_ALLOCATE(buf(data_pol%ntheta),err)
    SLL_ALLOCATE(this%f(data_pol%nr+1,data_pol%ntheta+1),err)
    SLL_ALLOCATE(this%f_fft(data_pol%nr+1,data_pol%ntheta+1),err)
    SLL_ALLOCATE(this%fdemi(data_pol%nr+1,data_pol%ntheta+1),err)
    SLL_ALLOCATE(this%phi(data_pol%nr+1,data_pol%ntheta+1),err)
    SLL_ALLOCATE(this%grad_phi(2,data_pol%nr+1,data_pol%ntheta+1),err)
    SLL_ALLOCATE(this%a(3*(data_pol%nr+1)),err)
    SLL_ALLOCATE(this%cts(7*(data_pol%nr+1)),err)
    SLL_ALLOCATE(this%ipiv(data_pol%nr+1),err)
    SLL_ALLOCATE(this%fk(data_pol%nr+1),err)
    SLL_ALLOCATE(this%phik(data_pol%nr+1),err)
    SLL_ALLOCATE(this%rr(data_pol%nr+1),err)
    SLL_ALLOCATE(this%ttheta(data_pol%ntheta+1),err)

    this%pfwd => fft_new_plan(data_pol%ntheta,buf,buf,FFT_FORWARD,FFT_NORMALIZE)
    this%pinv => fft_new_plan(data_pol%ntheta,buf,buf,FFT_INVERSE)

    SLL_ALLOCATE(this%data,err)
    this%data=data_pol

    do i=1,this%data%nr+1
       this%rr(i)=this%data%rmin+real(i-1,f64)*this%data%dr
    end do
    do i=1,this%data%ntheta+1
       this%ttheta(i)=real(i-1,f64)*this%data%dtheta
    end do

    this%spl_f => new_spline_2D(data_pol%nr+1,data_pol%ntheta+1,data_pol%rmin,data_pol%rmax,0._f64, 2._f64*sll_pi, &
         HERMITE_SPLINE, PERIODIC_SPLINE,&
         const_slope_x1_min = 0._f64,const_slope_x1_max = 0._f64)
    this%spl_a1 => new_spline_2D(data_pol%nr+1,data_pol%ntheta+1,data_pol%rmin,data_pol%rmax,0._f64, 2._f64*sll_pi, &
         HERMITE_SPLINE, PERIODIC_SPLINE,&
         const_slope_x1_min = 0._f64,const_slope_x1_max = 0._f64)
    this%spl_a2 => new_spline_2D(data_pol%nr+1,data_pol%ntheta+1,data_pol%rmin,data_pol%rmax,0._f64, 2._f64*sll_pi, &
         HERMITE_SPLINE, PERIODIC_SPLINE,&
         const_slope_x1_min = 0._f64,const_slope_x1_max = 0._f64)
    this%spl_phi => new_spline_2D(data_pol%nr+1,data_pol%ntheta+1,data_pol%rmin,data_pol%rmax,0._f64, 2._f64*sll_pi, &
         HERMITE_SPLINE, PERIODIC_SPLINE,&
         const_slope_x1_min = 0._f64,const_slope_x1_max = 0._f64)

    SLL_DEALLOCATE(buf,err)

  end function new_VP_dat_from_polar_data


  function new_VP_dat_from_all_dt_nbstep(nb_step,dt,rmin,rmax,nr,ntheta) result(this)

    implicit none

    sll_real64, intent(in) :: dt
    sll_real64, intent(in) :: rmin, rmax
    sll_int32, intent(in) :: nb_step,nr, ntheta

    type(polar_VP_data), pointer :: this

    sll_int32 :: err,i
    sll_real64,dimension(:),pointer :: buf

    SLL_ALLOCATE(this,err)
    SLL_ALLOCATE(buf(ntheta),err)
    SLL_ALLOCATE(this%f(nr+1,ntheta+1),err)
    SLL_ALLOCATE(this%f_fft(nr+1,ntheta+1),err)
    SLL_ALLOCATE(this%fdemi(nr+1,ntheta+1),err)
    SLL_ALLOCATE(this%phi(nr+1,ntheta+1),err)
    SLL_ALLOCATE(this%grad_phi(2,nr+1,ntheta+1),err)
    SLL_ALLOCATE(this%a(3*(nr+1)),err)
    SLL_ALLOCATE(this%cts(7*(nr+1)),err)
    SLL_ALLOCATE(this%ipiv(nr+1),err)
    SLL_ALLOCATE(this%fk(nr+1),err)
    SLL_ALLOCATE(this%phik(nr+1),err)
    SLL_ALLOCATE(this%rr(nr+1),err)
    SLL_ALLOCATE(this%ttheta(ntheta+1),err)

    this%pfwd => fft_new_plan(ntheta,buf,buf,FFT_FORWARD,FFT_NORMALIZE)
    this%pinv => fft_new_plan(ntheta,buf,buf,FFT_INVERSE)

    SLL_ALLOCATE(this%data,err)
    this%data => new_polar_data(dt,rmin,rmax,nb_step,nr,ntheta)

    do i=1,nr+1
       this%rr(i)=rmin+real(i-1,f64)*this%data%dr
    end do
    do i=1,ntheta+1
       this%ttheta(i)=real(i-1,f64)*this%data%dtheta
    end do

    this%spl_f => new_spline_2D(nr+1,ntheta+1,rmin,rmax,0._f64, 2._f64*sll_pi, &
         HERMITE_SPLINE, PERIODIC_SPLINE,&
         const_slope_x1_min = 0._f64,const_slope_x1_max = 0._f64)
    this%spl_a1 => new_spline_2D(nr+1,ntheta+1,rmin,rmax,0._f64, 2._f64*sll_pi, &
         HERMITE_SPLINE, PERIODIC_SPLINE,&
         const_slope_x1_min = 0._f64,const_slope_x1_max = 0._f64)
    this%spl_a2 => new_spline_2D(nr+1,ntheta+1,rmin,rmax,0._f64, 2._f64*sll_pi, &
         HERMITE_SPLINE, PERIODIC_SPLINE,&
         const_slope_x1_min = 0._f64,const_slope_x1_max = 0._f64)
    this%spl_phi => new_spline_2D(nr+1,ntheta+1,rmin,rmax,0._f64, 2._f64*sll_pi, &
         HERMITE_SPLINE, PERIODIC_SPLINE,&
         const_slope_x1_min = 0._f64,const_slope_x1_max = 0._f64)

    SLL_DEALLOCATE(buf,err)

  end function new_VP_dat_from_all_dt_nbstep


  function new_VP_dat_from_all_tf_nbstep(tf,rmin,rmax,nb_step,nr,ntheta) result(this)

    implicit none

    sll_real64, intent(in) :: tf,rmin, rmax
    sll_int32, intent(in) :: nb_step,nr, ntheta

    type(polar_VP_data), pointer :: this

    sll_int32 :: err,i
    sll_real64,dimension(:),pointer :: buf

    SLL_ALLOCATE(this,err)
    SLL_ALLOCATE(buf(ntheta),err)
    SLL_ALLOCATE(this%f(nr+1,ntheta+1),err)
    SLL_ALLOCATE(this%f_fft(nr+1,ntheta+1),err)
    SLL_ALLOCATE(this%fdemi(nr+1,ntheta+1),err)
    SLL_ALLOCATE(this%phi(nr+1,ntheta+1),err)
    SLL_ALLOCATE(this%grad_phi(2,nr+1,ntheta+1),err)
    SLL_ALLOCATE(this%a(3*(nr+1)),err)
    SLL_ALLOCATE(this%cts(7*(nr+1)),err)
    SLL_ALLOCATE(this%ipiv(nr+1),err)
    SLL_ALLOCATE(this%fk(nr+1),err)
    SLL_ALLOCATE(this%phik(nr+1),err)
    SLL_ALLOCATE(this%rr(nr+1),err)
    SLL_ALLOCATE(this%ttheta(ntheta+1),err)

    this%pfwd => fft_new_plan(ntheta,buf,buf,FFT_FORWARD,FFT_NORMALIZE)
    this%pinv => fft_new_plan(ntheta,buf,buf,FFT_INVERSE)

    do i=1,nr+1
       this%rr(i)=rmin+real(i-1,f64)*this%data%dr
    end do
    do i=1,ntheta+1
       this%ttheta(i)=real(i-1,f64)*this%data%dtheta
    end do

    SLL_ALLOCATE(this%data,err)
    this%data => new_polar_data(tf,rmin,rmax,nb_step,nr,ntheta)

    this%spl_f => new_spline_2D(nr+1,ntheta+1,rmin,rmax,0._f64, 2._f64*sll_pi, &
         HERMITE_SPLINE, PERIODIC_SPLINE,&
         const_slope_x1_min = 0._f64,const_slope_x1_max = 0._f64)
    this%spl_a1 => new_spline_2D(nr+1,ntheta+1,rmin,rmax,0._f64, 2._f64*sll_pi, &
         HERMITE_SPLINE, PERIODIC_SPLINE,&
         const_slope_x1_min = 0._f64,const_slope_x1_max = 0._f64)
    this%spl_a2 => new_spline_2D(nr+1,ntheta+1,rmin,rmax,0._f64, 2._f64*sll_pi, &
         HERMITE_SPLINE, PERIODIC_SPLINE,&
         const_slope_x1_min = 0._f64,const_slope_x1_max = 0._f64)
    this%spl_phi => new_spline_2D(nr+1,ntheta+1,rmin,rmax,0._f64, 2._f64*sll_pi, &
         HERMITE_SPLINE, PERIODIC_SPLINE,&
         const_slope_x1_min = 0._f64,const_slope_x1_max = 0._f64)

    SLL_DEALLOCATE(buf,err)

  end function new_VP_dat_from_all_tf_nbstep


  function new_VP_dat_from_all_dt_tf(dt,tf,rmin,rmax,nr,ntheta) result(this)

    implicit none

    sll_real64, intent(in) :: dt
    sll_real64, intent(in) :: tf,rmin, rmax
    sll_int32, intent(in) :: nr, ntheta

    type(polar_VP_data), pointer :: this

    sll_int32 :: err,i
    sll_real64,dimension(:),pointer :: buf

    SLL_ALLOCATE(this,err)
    SLL_ALLOCATE(buf(ntheta),err)
    SLL_ALLOCATE(this%f(nr+1,ntheta+1),err)
    SLL_ALLOCATE(this%f_fft(nr+1,ntheta+1),err)
    SLL_ALLOCATE(this%fdemi(nr+1,ntheta+1),err)
    SLL_ALLOCATE(this%phi(nr+1,ntheta+1),err)
    SLL_ALLOCATE(this%grad_phi(2,nr+1,ntheta+1),err)
    SLL_ALLOCATE(this%a(3*(nr+1)),err)
    SLL_ALLOCATE(this%cts(7*(nr+1)),err)
    SLL_ALLOCATE(this%ipiv(nr+1),err)
    SLL_ALLOCATE(this%fk(nr+1),err)
    SLL_ALLOCATE(this%phik(nr+1),err)
    SLL_ALLOCATE(this%rr(nr+1),err)
    SLL_ALLOCATE(this%ttheta(ntheta+1),err)

    this%pfwd => fft_new_plan(ntheta,buf,buf,FFT_FORWARD,FFT_NORMALIZE)
    this%pinv => fft_new_plan(ntheta,buf,buf,FFT_INVERSE)

    SLL_ALLOCATE(this%data,err)
    this%data => new_polar_data(dt,tf,rmin,rmax,nr,ntheta)

    do i=1,nr+1
       this%rr(i)=rmin+real(i-1,f64)*this%data%dr
    end do
    do i=1,ntheta+1
       this%ttheta(i)=real(i-1,f64)*this%data%dtheta
    end do

    this%spl_f => new_spline_2D(nr+1,ntheta+1,rmin,rmax,0._f64, 2._f64*sll_pi, &
         HERMITE_SPLINE, PERIODIC_SPLINE,&
         const_slope_x1_min = 0._f64,const_slope_x1_max = 0._f64)
    this%spl_a1 => new_spline_2D(nr+1,ntheta+1,rmin,rmax,0._f64, 2._f64*sll_pi, &
         HERMITE_SPLINE, PERIODIC_SPLINE,&
         const_slope_x1_min = 0._f64,const_slope_x1_max = 0._f64)
    this%spl_a2 => new_spline_2D(nr+1,ntheta+1,rmin,rmax,0._f64, 2._f64*sll_pi, &
         HERMITE_SPLINE, PERIODIC_SPLINE,&
         const_slope_x1_min = 0._f64,const_slope_x1_max = 0._f64)

    SLL_DEALLOCATE(buf,err)

  end function new_VP_dat_from_all_dt_tf

!================================
!  deletion of polar_vp_data
!================================

  !>VP_data_delete(this)
  !>delete object of type polar_vp_data
  subroutine VP_data_delete(this)

    implicit none

    type(polar_VP_data), intent(inout), pointer :: this
    sll_int32 :: err
    if (associated(this)) then
       call fft_delete_plan(this%pfwd)
       call fft_delete_plan(this%pinv)
       SLL_DEALLOCATE_ARRAY(this%f,err)
       SLL_DEALLOCATE_ARRAY(this%f_fft,err)
       SLL_DEALLOCATE_ARRAY(this%fdemi,err)
       SLL_DEALLOCATE_ARRAY(this%phi,err)
       SLL_DEALLOCATE_ARRAY(this%grad_phi,err)
       SLL_DEALLOCATE_ARRAY(this%a,err)
       SLL_DEALLOCATE_ARRAY(this%cts,err)
       SLL_DEALLOCATE_ARRAY(this%ipiv,err)
       SLL_DEALLOCATE_ARRAY(this%rr,err)
       SLL_DEALLOCATE_ARRAY(this%ttheta,err)

       call delete_spline_2d(this%spl_f)
       call delete_spline_2d(this%spl_a1)
       call delete_spline_2d(this%spl_a2)
       this%data => null()
    end if

    SLL_DEALLOCATE(this,err)

  end subroutine VP_data_delete

!==================================
!  construction of polar_vp_rk4
!==================================

  !>function new_polar_vp_rk4(nr,ntheta)
  !>initialize a polar_vp_rk4 object
  function new_polar_vp_rk4(nr,ntheta) result(this)

    implicit none

    type(polar_vp_rk4), pointer :: this
    sll_int32, intent(in) :: nr,ntheta
    sll_int32 :: err

    SLL_ALLOCATE(this,err)
    SLL_ALLOCATE(this%r1(nr+1,ntheta+1),err)
    SLL_ALLOCATE(this%r2(nr+1,ntheta+1),err)
    SLL_ALLOCATE(this%r3(nr+1,ntheta+1),err)
    SLL_ALLOCATE(this%r4(nr+1,ntheta+1),err)
    SLL_ALLOCATE(this%theta1(nr+1,ntheta+1),err)
    SLL_ALLOCATE(this%theta2(nr+1,ntheta+1),err)
    SLL_ALLOCATE(this%theta3(nr+1,ntheta+1),err)
    SLL_ALLOCATE(this%theta4(nr+1,ntheta+1),err)
  end function new_polar_vp_rk4

!==============================
!  deletion of polar_vp_rk4
!==============================

  !>subroutine vp_rk4_delete(this)
  !>delete a polar_vp_rk4 object
  subroutine vp_rk4_delete(this)

    implicit none

    type(polar_vp_rk4),intent(inout), pointer :: this
    sll_int32 :: err

    if (associated(this)) then
       SLL_DEALLOCATE_ARRAY(this%r1,err)
       SLL_DEALLOCATE_ARRAY(this%r2,err)
       SLL_DEALLOCATE_ARRAY(this%r3,err)
       SLL_DEALLOCATE_ARRAY(this%r4,err)
       SLL_DEALLOCATE_ARRAY(this%theta1,err)
       SLL_DEALLOCATE_ARRAY(this%theta2,err)
       SLL_DEALLOCATE_ARRAY(this%theta3,err)
       SLL_DEALLOCATE_ARRAY(this%theta4,err)
    end if
    SLL_DEALLOCATE(this,err)

  end subroutine vp_rk4_delete

end module polar_kind
