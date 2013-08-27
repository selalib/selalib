module sll_vp_cartesian_2d
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_field_2d.h"
#ifdef STDF95
  use sll_cubic_spline_interpolator_1d
#else
  use sll_module_interpolators_1d_base
  use sll_time_splitting
#endif
  use distribution_function
  use sll_poisson_1d_periodic
  implicit none

  type :: app_field_params
     sll_real64 :: Edrmax, tflat, tL, tR, twL, twR, t0
     sll_real64 :: kmode, omegadr
     logical    :: turn_drive_off, driven
  end type app_field_params

#ifdef STDF95
  type  :: vp_cartesian_2d
    sll_real64 :: current_time = 0.0_f64
    type(cubic_spline_1d_interpolator), pointer    :: interpx, interpv
#else
  type, extends(time_splitting) :: vp_cartesian_2d
     class(sll_interpolator_1d_base), pointer    :: interpx, interpv
#endif
     type(sll_distribution_function_2d), pointer   :: dist_func
     type(poisson_1d_periodic), pointer            :: poisson_1d
     sll_int32 :: Ncx, Ncv
     type(app_field_params)  :: params
#ifndef STDF95
   contains
     procedure, pass(this) :: operator1 => advection_x
     procedure, pass(this) :: operator2 => advection_v
#endif
  end type vp_cartesian_2d

contains

  subroutine vp_cartesian_2d_initialize(this, dist_func, poisson_1d, Ncx, Ncv, interpx, interpv, params)
    type(vp_cartesian_2d) :: this 
    type(sll_distribution_function_2d), target   :: dist_func
    type(poisson_1d_periodic), target            :: poisson_1d
    sll_int32 :: Ncx, Ncv
#ifdef STDF95
    type(cubic_spline_1d_interpolator), pointer    :: interpx, interpv
#else
    class(sll_interpolator_1d_base), pointer    :: interpx, interpv
#endif
    type(app_field_params)  :: params
    this%dist_func  => dist_func
    this%poisson_1d => poisson_1d
    this%Ncx = Ncx
    this%Ncv = Ncv
    this%interpx => interpx
    this%interpv => interpv
    this%params = params
  end subroutine vp_cartesian_2d_initialize

#ifdef STDF95
  subroutine vp_cartesian_2d_operator1(this, dt)
    type(vp_cartesian_2d) :: this 
#else
  subroutine advection_x(this, dt)
    class(vp_cartesian_2d) :: this 
#endif
    sll_real64, intent(in) :: dt
    ! local variables
    sll_real64, dimension(:), pointer :: f1d
    sll_real64 :: displacement
    sll_int32 :: j
    sll_real64 :: vmin, vmax, delta_v

#ifdef STDF95
    vmin = x2_node_discrete(this%dist_func%extend_type%mesh, 1,1)
    vmax = x2_node_discrete(this%dist_func%extend_type%mesh, 1,this%Ncv+1)
    delta_v = (vmax - vmin) /  this%dist_func%extend_type%mesh%mesh%num_cells2
    do j = 1, this%Ncv+1
       displacement = (vmin + (j-1) * delta_v) * dt
       f1d => FIELD_DATA(this%dist_func%extend_type) (:,j)
       f1d = cubic_spline_interpolate_array_at_displacement( this%interpx, this%Ncx+1, f1d, displacement)
    end do
#else    
    vmin = this%dist_func%mesh%x2_at_node(1,1)
    vmax = this%dist_func%mesh%x2_at_node(1,this%Ncv+1)
    delta_v = (vmax - vmin) /  this%dist_func%mesh%mesh%num_cells2
    do j = 1, this%Ncv+1
       displacement = (vmin + (j-1) * delta_v) * dt
       f1d => FIELD_DATA(this%dist_func) (:,j)
       f1d = this%interpx%interpolate_array_disp(this%Ncx+1, f1d, displacement)
    end do
#endif
  end subroutine 

#ifdef STDF95
  subroutine vp_cartesian_2d_operator2(this, dt)
    type(vp_cartesian_2d) :: this 
#else
  subroutine advection_v(this, dt)
    class(vp_cartesian_2d) :: this 
#endif
    sll_real64, intent(in) :: dt
    ! local variables
    sll_real64 :: time
    sll_real64, dimension(this%Ncx) :: rho, efield, e_app
    sll_real64, dimension(:), pointer :: f1d
    sll_real64 :: displacement
    sll_real64 :: adr
    sll_real64 :: arg
    sll_int32 :: i
    sll_real64 :: xmin, xmax, delta_x
    sll_real64 :: vmin, vmax, delta_v

    
    time = this%current_time
#ifdef STDF95
    xmin = x1_node_discrete(this%dist_func%extend_type%mesh, 1,1)
    xmax = x1_node_discrete(this%dist_func%extend_type%mesh, this%Ncx+1,1)
    delta_x = (xmax - xmin) /  this%dist_func%extend_type%mesh%mesh%num_cells1
    vmin = x2_node_discrete(this%dist_func%extend_type%mesh, 1,1)
    vmax = x2_node_discrete(this%dist_func%extend_type%mesh, 1,this%Ncv+1)
    delta_v = (vmax - vmin) /  this%dist_func%extend_type%mesh%mesh%num_cells2
#else
    xmin = this%dist_func%mesh%x1_at_node(1,1)
    xmax = this%dist_func%mesh%x1_at_node(this%Ncx+1,1)
    delta_x = (xmax - xmin) /  this%dist_func%mesh%mesh%num_cells1
    vmin = this%dist_func%mesh%x2_at_node(1,1)
    vmax = this%dist_func%mesh%x2_at_node(1,this%Ncv+1)
    delta_v = (vmax - vmin) /  this%dist_func%mesh%mesh%num_cells2
#endif
    
    ! compute electric field
    !-----------------------
#ifdef STDF95
    rho = 1.0_f64 - delta_v * sum(FIELD_DATA(this%dist_func%extend_type), DIM = 2)
#else
    rho = 1.0_f64 - delta_v * sum(FIELD_DATA(this%dist_func), DIM = 2)
#endif
    call solve(this%poisson_1d, efield, rho)
    if (this%params%driven) then
       call PFenvelope(adr, time, this%params)
       do i = 1, this%Ncx + 1
          arg = this%params%kmode * real(i-1,8) * delta_x - this%params%omegadr*time
          e_app(i) = this%params%Edrmax * adr * this%params%kmode * sin(arg)
       enddo
    endif
    ! do advection for given electric field
    do i = 1, this%Ncx+1
        displacement = -(efield(i)+e_app(i)) * 0.5_f64 * dt
#ifdef STDF95
        f1d => FIELD_DATA(this%dist_func%extend_type) (i,:) 
        f1d = cubic_spline_interpolate_array_at_displacement(this%interpv, this%Ncv+1, f1d, displacement)
#else
        f1d => FIELD_DATA(this%dist_func) (i,:) 
        f1d = this%interpv%interpolate_array_disp(this%Ncv+1, f1d, displacement)
#endif
     end do
  end subroutine

  elemental function f_equilibrium(v)
    sll_real64, intent(in) :: v
    sll_real64 :: f_equilibrium

    f_equilibrium = 1.0_f64/sqrt(2*sll_pi)*exp(-0.5_f64*v*v)
  end function f_equilibrium

  subroutine PFenvelope(S, t, params)

    ! DESCRIPTION
    ! -----------
    ! S: the wave form at a given point in time. This wave form is 
    !    not scaled (its maximum value is 1).
    ! t: the time at which the envelope is being evaluated
    ! tflat, tL, tR, twL, twR, tstart, t0: the parameters defining the
    !    envelope, defined in the main portion of this program.
    ! turn_drive_off: 1 if the drive should be turned off after a time
    !    tflat, and 0 otherwise

    sll_real64, intent(out) :: S
    sll_real64, intent(in) :: t
    type(app_field_params)  :: params
    ! local variables
    sll_real64 :: t0, twL, twR, tflat, tL, tR
    sll_real64 :: epsilon

    tflat = params%tflat
    tL = params%tR 
    twL = params%twL
    twR = params%twR 
    t0 = params%t0
    ! The envelope function is defined such that it is zero at t0,
    ! rises to 1 smoothly, stay constant for tflat, and returns
    ! smoothly to zero.
    if(params%turn_drive_off) then
       epsilon = 0.5*(tanh((t0-tL)/twL) - tanh((t0-tR)/twR))
       S = 0.5*(tanh((t-tL)/twL) - tanh((t-tR)/twR)) - epsilon
       S = S / (1-epsilon)
    else
       epsilon = 0.5*(tanh((t0-tL)/twL) + 1)
       S = 0.5*(tanh((t-tL)/twL) + 1) - epsilon
       S = S / (1-epsilon)
    endif
    if(S<0) then
       S = 0.
    endif
    return
  end subroutine PFenvelope

end module sll_vp_cartesian_2d
