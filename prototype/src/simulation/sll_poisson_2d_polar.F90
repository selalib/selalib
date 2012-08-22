module poisson_polar
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"

  use sll_fft
  use sll_tridiagonal
  use numeric_constants
  implicit none

  !>type polar_data
  !>generic type for problems in polar
  !>contains size of time and space steps, boundaries, number of steps in time and space,
  !>and the final time
  type sll_polar_data
     sll_real64 :: dt, dr, dtheta
     sll_real64 :: tf,rmin, rmax
     sll_int32 :: nb_step,nr, ntheta
  end type sll_polar_data

  !>type sll_plan_poisson_polar
  !>type for the Poisson solver in polar coordinate
  type sll_plan_poisson_polar
     type(sll_polar_data), pointer :: data
     type(sll_fft_plan), pointer :: pfwd,pinv
     sll_real64, dimension(:,:), allocatable :: f_fft
     sll_comp64, dimension(:), allocatable :: fk,phik
     !for the tridiagonal solver
     sll_real64, dimension(:), allocatable :: a,cts
     sll_int32, dimension(:), allocatable :: ipiv
  end type sll_plan_poisson_polar

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

contains
  
!=========================================
!  beginnig of creation of polar data
!=========================================

  function new_polar_data_dt_nbstep(nb_step,dt,rmin,rmax,nr,ntheta) result(this)

    implicit none

    type(sll_polar_data), pointer :: this
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

    type(sll_polar_data), pointer :: this
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

    type(sll_polar_data), pointer :: this
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

!========================================
!  creation of sll_plan_poisson_polar
!========================================

  !>new_plan_poisson_polar(data)
  !>build a sll_plan_poisson_polar object for the Poisson solver in polar coordinate
  !>data : sll_polar_data object
  function new_plan_poisson_polar(data) result(this)

    implicit none

    type(sll_plan_poisson_polar), pointer :: this
    type(sll_polar_data), intent(in) :: data

    sll_int32 :: err
    sll_real64, dimension(:), allocatable :: buf

    SLL_ALLOCATE(this,err)
    SLL_ALLOCATE(buf(data%ntheta),err)
    SLL_ALLOCATE(this%f_fft(data%nr+1,data%ntheta+1),err)
    SLL_ALLOCATE(this%fk(data%nr+1),err)
    SLL_ALLOCATE(this%phik(data%nr+1),err)
    SLL_ALLOCATE(this%a(3*(data%nr+1)),err)
    SLL_ALLOCATE(this%cts(7*(data%nr+1)),err)
    SLL_ALLOCATE(this%ipiv(data%nr+1),err)

    this%pfwd => fft_new_plan(data%ntheta,buf,buf,FFT_FORWARD,FFT_NORMALIZE)
    this%pinv => fft_new_plan(data%ntheta,buf,buf,FFT_INVERSE)

    SLL_ALLOCATE(this%data,err)
    this%data=data

    SLL_DEALLOCATE_ARRAY(buf,err)

  end function new_plan_poisson_polar

!======================================
! deletion of sll_plan_poisson_polar
!======================================

  !>delete_plan_poisson_polar(plan)
  !>delete a sll_plan_poisson_polar object
  subroutine delete_plan_poisson_polar(this)

    implicit none

    type(sll_plan_poisson_polar), intent(inout), pointer :: this
    sll_int32 :: err
    if (associated(this)) then
       call fft_delete_plan(this%pfwd)
       call fft_delete_plan(this%pinv)
       SLL_DEALLOCATE_ARRAY(this%f_fft,err)
       SLL_DEALLOCATE_ARRAY(this%fk,err)
       SLL_DEALLOCATE_ARRAY(this%phik,err)
       SLL_DEALLOCATE_ARRAY(this%a,err)
       SLL_DEALLOCATE_ARRAY(this%cts,err)
       SLL_DEALLOCATE_ARRAY(this%ipiv,err)

       this%data => null()

       SLL_DEALLOCATE(this,err)
    end if

  end subroutine delete_plan_poisson_polar

!===================
!  Poisson solver
!===================

  !>subroutine poisson_solve_polar(plan,f,phi)
  !>poisson solver for polar system : -\Delta (phi)=f
  !>plan : sll_plan_poisson_polar, contains data for the solver
  !>f : distribution function, size (nr+1)*(ntheta+1), input
  !>phi : unknown field, size (nr+1)*(ntheta+1), output
  !>initialization must be done outside the solver
  subroutine poisson_solve_polar(plan,f,phi)

    implicit none

    type(sll_plan_poisson_polar), intent(inout), pointer :: plan
    sll_real64, dimension(plan%data%nr+1,plan%data%ntheta+1), intent(in) :: f
    sll_real64, dimension(plan%data%nr+1,plan%data%ntheta+1), intent(out) :: phi

    sll_real64 :: rmin,dr
    sll_int32 :: nr, ntheta

    sll_real64 :: r
    sll_int32::i,ind_k

    nr=plan%data%nr
    ntheta=plan%data%ntheta
    rmin=plan%data%rmin
    dr=plan%data%dr

    plan%f_fft=f

    do i=1,nr+1
       call fft_apply_plan(plan%pfwd,plan%f_fft(i,1:ntheta),plan%f_fft(i,1:ntheta))
    end do

   ! poisson solver
    do ind_k=0,ntheta/2
       do i=1,nr+1
          r=rmin+real(i-1,f64)*dr
          plan%a(3*i)=-1.0_f64/dr**2-1.0_f64/(2.0_f64*dr*r)
          plan%a(3*i-1)=2.0_f64/dr**2+(real(ind_k,f64)/r)**2
          plan%a(3*i-2)=-1.0_f64/dr**2+1.0_f64/(2.0_f64*dr*r)
          plan%fk(i)=fft_get_mode(plan%pfwd,plan%f_fft(i,1:ntheta),ind_k)
       enddo

       plan%a(1)=0.0_f64
       plan%a(3*nr+3)=0.0_f64
       !a(2)=1.0_f64
       !a(3*nr+2)=1.0_f64

       call setup_cyclic_tridiag(plan%a,nr+1,plan%cts,plan%ipiv)
       call solve_cyclic_tridiag(plan%cts,plan%ipiv,plan%fk,nr+1,plan%phik)

       do i=1,nr+1
          call fft_set_mode(plan%pinv,phi(i,1:ntheta),plan%phik(i),ind_k)
       end do
    end do

print*,'ok'
    ! FFT INVERSE
    do i=1,Nr+1
       call fft_apply_plan(plan%pinv,phi(i,1:ntheta),phi(i,1:ntheta))
    end do

    phi(1,:)=0.0_f64
    phi(nr+1,:)=0.0_f64
    phi(:,ntheta+1)=phi(:,1)

  end subroutine poisson_solve_polar

end module poisson_polar
