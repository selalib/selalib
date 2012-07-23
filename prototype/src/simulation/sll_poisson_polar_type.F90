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
     sll_real64, dimension(:,:), pointer :: f,phi,f_fft,fdemi
     sll_real64, dimension(:,:,:), pointer :: grad_phi
     sll_comp64, dimension(:), pointer :: fk,phik
     type(sll_spline_2D), pointer :: spl_f, spl_a1, spl_a2, spl_phi
     !for the tridiagonal solver
     sll_real64, dimension(:), pointer :: cts,a
     sll_int32, dimension(:), pointer :: ipiv
  end type polar_VP_data

  !>type polar_VP_rk4
  !>used for RK4 in Vlasov Poisson
  !>contains r1, r2, r3, r4, theta1, theta2, theta3, theta4 for all points
  type polar_VP_rk4
     sll_real64, dimension(:,:), pointer :: r1,r2,r3,r4
     sll_real64, dimension(:,:), pointer :: theta1, theta2, theta3, theta4
  end type polar_VP_rk4

!>   //==============\\
!>   ||  INTERFACES  ||
!>   \\==============//

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
print*,'in new_polar_data_dt_nbstep',ntheta
    this%dt=dt
    this%dr=(rmax-rmin)/real(nr,f64)
    this%dtheta=2.0_f64*sll_pi/real(ntheta,f64)
    this%rmin=rmin
    this%rmax=rmax
    this%nb_step=nb_step
    this%nr=nr
    this%ntheta=ntheta
    this%tf=dt*real(nb_step,f64)
print*,'in new_polar_data_dt_nbstep',this%ntheta
  end function new_polar_data_dt_nbstep


  function new_polar_data_dt_tf(dt,tf,rmin,rmax,nr,ntheta) result(this)

    implicit none

    type(polar_data), pointer :: this
    sll_real64, intent(in) :: dt
    sll_real64, intent(in) :: tf,rmin, rmax
    sll_int32, intent(in) :: nr, ntheta

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

    sll_int32 :: err
    sll_real64,dimension(:),pointer :: buf
print*,'new_VP_dat_from_polar_data, entre'
print*,data_pol%ntheta
    SLL_ALLOCATE(buf(data_pol%ntheta),err)
print*,1
    SLL_ALLOCATE(this%f(data_pol%nr+1,data_pol%ntheta+1),err)
print*,2
    SLL_ALLOCATE(this%f_fft(data_pol%nr+1,data_pol%ntheta+1),err)
print*,3
    SLL_ALLOCATE(this%fdemi(data_pol%nr+1,data_pol%ntheta+1),err)
print*,4
    SLL_ALLOCATE(this%phi(data_pol%nr+1,data_pol%ntheta+1),err)
print*,5
    SLL_ALLOCATE(this%grad_phi(2,data_pol%nr+1,data_pol%ntheta+1),err)
print*,6
    SLL_ALLOCATE(this%a(3*data_pol%nr+1),err)
    SLL_ALLOCATE(this%cts(7*data_pol%nr+1),err)
    SLL_ALLOCATE(this%ipiv(data_pol%nr+1),err)
    SLL_ALLOCATE(this%fk(data_pol%ntheta),err)
    SLL_ALLOCATE(this%phik(data_pol%ntheta),err)
print*,'alloc ok'
    this%pfwd => fft_new_plan(data_pol%ntheta,buf,buf,FFT_FORWARD,FFT_NORMALIZE)
    this%pinv => fft_new_plan(data_pol%ntheta,buf,buf,FFT_INVERSE)

    this%data=data_pol
print*,'fft+data ok'
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
print*,'spl ok'
    SLL_DEALLOCATE(buf,err)

  end function new_VP_dat_from_polar_data


  function new_VP_dat_from_all_dt_nbstep(nb_step,dt,rmin,rmax,nr,ntheta) result(this)

    implicit none

    sll_real64, intent(in) :: dt
    sll_real64, intent(in) :: rmin, rmax
    sll_int32, intent(in) :: nb_step,nr, ntheta

    type(polar_VP_data), pointer :: this

    sll_int32 :: err
    sll_real64,dimension(:),pointer :: buf

    SLL_ALLOCATE(buf(ntheta),err)
    SLL_ALLOCATE(this%f(nr+1,ntheta+1),err)
    SLL_ALLOCATE(this%f_fft(nr+1,ntheta+1),err)
    SLL_ALLOCATE(this%fdemi(nr+1,ntheta+1),err)
    SLL_ALLOCATE(this%phi(nr+1,ntheta+1),err)
    SLL_ALLOCATE(this%grad_phi(2,nr+1,ntheta+1),err)
    SLL_ALLOCATE(this%a(3*(nr+1)),err)
    SLL_ALLOCATE(this%cts(7*(nr+1)),err)
    SLL_ALLOCATE(this%ipiv(nr+1),err)
    SLL_ALLOCATE(this%fk(ntheta),err)
    SLL_ALLOCATE(this%phik(ntheta),err)

    this%pfwd => fft_new_plan(ntheta,buf,buf,FFT_FORWARD,FFT_NORMALIZE)
    this%pinv => fft_new_plan(ntheta,buf,buf,FFT_INVERSE)

    this%data => new_polar_data(dt,rmin,rmax,nb_step,nr,ntheta)

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

    sll_int32 :: err
    sll_real64,dimension(:),pointer :: buf

    SLL_ALLOCATE(buf(ntheta),err)
    SLL_ALLOCATE(this%f(nr+1,ntheta+1),err)
    SLL_ALLOCATE(this%f_fft(nr+1,ntheta+1),err)
    SLL_ALLOCATE(this%fdemi(nr+1,ntheta+1),err)
    SLL_ALLOCATE(this%phi(nr+1,ntheta+1),err)
    SLL_ALLOCATE(this%grad_phi(2,nr+1,ntheta+1),err)
    SLL_ALLOCATE(this%a(3*(nr+1)),err)
    SLL_ALLOCATE(this%cts(7*(nr+1)),err)
    SLL_ALLOCATE(this%ipiv(nr+1),err)
    SLL_ALLOCATE(this%fk(ntheta),err)
    SLL_ALLOCATE(this%phik(ntheta),err)

    this%pfwd => fft_new_plan(ntheta,buf,buf,FFT_FORWARD,FFT_NORMALIZE)
    this%pinv => fft_new_plan(ntheta,buf,buf,FFT_INVERSE)

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

    sll_int32 :: err
    sll_real64,dimension(:),pointer :: buf

    SLL_ALLOCATE(buf(ntheta),err)
    SLL_ALLOCATE(this%f(nr+1,ntheta+1),err)
    SLL_ALLOCATE(this%f_fft(nr+1,ntheta+1),err)
    SLL_ALLOCATE(this%fdemi(nr+1,ntheta+1),err)
    SLL_ALLOCATE(this%phi(nr+1,ntheta+1),err)
    SLL_ALLOCATE(this%grad_phi(2,nr+1,ntheta+1),err)
    SLL_ALLOCATE(this%a(3*(nr+1)),err)
    SLL_ALLOCATE(this%cts(7*(nr+1)),err)
    SLL_ALLOCATE(this%ipiv(nr+1),err)
    SLL_ALLOCATE(this%fk(ntheta),err)
    SLL_ALLOCATE(this%phik(ntheta),err)

    this%pfwd => fft_new_plan(ntheta,buf,buf,FFT_FORWARD,FFT_NORMALIZE)
    this%pinv => fft_new_plan(ntheta,buf,buf,FFT_INVERSE)

    this%data => new_polar_data(dt,tf,rmin,rmax,nr,ntheta)

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

    call fft_delete(this%pfwd)
    call fft_delete(this%pinv)

    if (associated(this%f)) then
       SLL_DEALLOCATE(this%f,err)
    end if
    if (associated(this%f_fft)) then
       SLL_DEALLOCATE(this%f_fft,err)
    end if
    if (associated(this%phi)) then
       SLL_DEALLOCATE(this%phi,err)
    end if
    if (associated(this%grad_phi)) then
       SLL_DEALLOCATE(this%grad_phi,err)
    end if
    if (associated(this%a)) then
       SLL_DEALLOCATE(this%a,err)
    end if
    if (associated(this%cts)) then
       SLL_DEALLOCATE(this%cts,err)
    end if
    if (associated(this%ipiv)) then
       SLL_DEALLOCATE(this%ipiv,err)
    end if

    call delete_spline_2d(this%spl_f)
    call delete_spline_2d(this%spl_a1)
    call delete_spline_2d(this%spl_a2)

    this%data => null()
    this%f => null()
    this%phi => null()
    this%grad_phi => null()

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

    if (associated(this%r1)) then
       SLL_DEALLOCATE(this%r1,err)
    end if
    if (associated(this%r2)) then
       SLL_DEALLOCATE(this%r2,err)
    end if
    if (associated(this%r3)) then
       SLL_DEALLOCATE(this%r3,err)
    end if
    if (associated(this%r4)) then
       SLL_DEALLOCATE(this%r4,err)
    end if
    if (associated(this%theta1)) then
       SLL_DEALLOCATE(this%theta1,err)
    end if
    if (associated(this%theta2)) then
       SLL_DEALLOCATE(this%theta1,err)
    end if
    if (associated(this%theta3)) then
       SLL_DEALLOCATE(this%theta1,err)
    end if
    if (associated(this%theta4)) then
       SLL_DEALLOCATE(this%theta1,err)
    end if
    this => null()

  end subroutine vp_rk4_delete

end module polar_kind
