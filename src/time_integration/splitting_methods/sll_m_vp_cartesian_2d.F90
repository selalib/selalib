#ifndef DOXYGEN_SHOULD_SKIP_THIS

module sll_m_vp_cartesian_2d
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_field_2d.h"

   use sll_m_cartesian_meshes, only: &
      sll_t_cartesian_mesh_2d

   use sll_m_constants, only: &
      sll_p_pi

   use sll_m_distribution_function, only: &
      sll_t_distribution_function_2d

   use sll_m_interpolators_1d_base, only: &
      sll_c_interpolator_1d

   use sll_m_poisson_1d_periodic, only: &
      sll_t_poisson_1d_periodic, &
      sll_o_solve

   use sll_m_time_splitting, only: &
      sll_c_time_splitting

   implicit none

   private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   type :: app_field_params
      sll_real64 :: Edrmax, tflat, tL, tR, twL, twR, t0
      sll_real64 :: kmode, omegadr
      logical    :: turn_drive_off, driven
   end type app_field_params

   type, extends(sll_c_time_splitting) :: vp_cartesian_2d
      class(sll_c_interpolator_1d), pointer    :: interpx, interpv
      type(sll_t_distribution_function_2d), pointer   :: dist_func
      type(sll_t_poisson_1d_periodic), pointer            :: poisson_1d
      sll_int32 :: Ncx, Ncv
      type(app_field_params)  :: params
   contains
      procedure, pass(this) :: operator1 => advection_x
      procedure, pass(this) :: operator2 => advection_v
   end type vp_cartesian_2d

contains

   subroutine vp_cartesian_2d_initialize(this, dist_func, poisson_1d, Ncx, Ncv, interpx, interpv, params)
      type(vp_cartesian_2d) :: this
      type(sll_t_distribution_function_2d), target   :: dist_func
      type(sll_t_poisson_1d_periodic), target            :: poisson_1d
      sll_int32 :: Ncx, Ncv
      class(sll_c_interpolator_1d), pointer    :: interpx, interpv
      type(app_field_params)  :: params
      this%dist_func => dist_func
      this%poisson_1d => poisson_1d
      this%Ncx = Ncx
      this%Ncv = Ncv
      this%interpx => interpx
      this%interpv => interpv
      this%params = params
   end subroutine vp_cartesian_2d_initialize

   subroutine advection_x(this, dt)
      class(vp_cartesian_2d) :: this
      sll_real64, intent(in) :: dt
      ! local variables
      sll_real64, dimension(:), pointer :: f1d
      sll_real64 :: displacement
      sll_int32 :: j
      sll_real64 :: vmin, vmax, delta_v

      associate (mesh => this%dist_func%transf%mesh)

         vmin = this%dist_func%transf%x2_at_node(1, 1)
         vmax = this%dist_func%transf%x2_at_node(1, this%Ncv + 1)
         delta_v = (vmax - vmin)/mesh%num_cells2

      end associate

      do j = 1, this%Ncv + 1
         displacement = -(vmin + (j - 1)*delta_v)*dt
         f1d => FIELD_DATA(this%dist_func) (:, j)
         call this%interpx%interpolate_array_disp_inplace(this%Ncx + 1, f1d, displacement)
      end do
   end subroutine

   subroutine advection_v(this, dt)
      class(vp_cartesian_2d) :: this
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

      associate (mesh => this%dist_func%transf%mesh)

         xmin = this%dist_func%transf%x1_at_node(1, 1)
         xmax = this%dist_func%transf%x1_at_node(this%Ncx + 1, 1)
         delta_x = (xmax - xmin)/mesh%num_cells1
         vmin = this%dist_func%transf%x2_at_node(1, 1)
         vmax = this%dist_func%transf%x2_at_node(1, this%Ncv + 1)
         delta_v = (vmax - vmin)/mesh%num_cells2

      end associate

      ! compute electric field
      !-----------------------

      rho = 1.0_f64 - delta_v*sum(FIELD_DATA(this%dist_func), DIM=2)
      call sll_o_solve(this%poisson_1d, efield, rho)
      if (this%params%driven) then
         call PF_envelope(adr, time, this%params)
         do i = 1, this%Ncx + 1
            arg = this%params%kmode*real(i - 1, 8)*delta_x - this%params%omegadr*time
            e_app(i) = this%params%Edrmax*adr*this%params%kmode*sin(arg)
         end do
      end if
      ! do advection for given electric field
      do i = 1, this%Ncx + 1
         displacement = (efield(i) + e_app(i))*0.5_f64*dt
         f1d => FIELD_DATA(this%dist_func) (i, :)
         call this%interpv%interpolate_array_disp_inplace(this%Ncv + 1, f1d, displacement)
      end do
   end subroutine

   elemental function f_equilibrium(v)
      sll_real64, intent(in) :: v
      sll_real64 :: f_equilibrium

      f_equilibrium = 1.0_f64/sqrt(2*sll_p_pi)*exp(-0.5_f64*v*v)
   end function f_equilibrium

   subroutine PF_envelope(S, t, params)

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
      type(app_field_params), intent(in)  :: params
      ! local variables
      sll_real64 :: t0, twL, twR, tflat, tL, tR
      sll_real64 :: epsilon

      tflat = params%tflat
      tR = params%tR
      tL = params%tL
      twL = params%twL
      twR = params%twR
      t0 = params%t0
      ! The envelope function is defined such that it is zero at t0,
      ! rises to 1 smoothly, stay constant for tflat, and returns
      ! smoothly to zero.
      if (params%turn_drive_off) then
         epsilon = 0.5*(tanh((t0 - tL)/twL) - tanh((t0 - tR)/twR))
         S = 0.5*(tanh((t - tL)/twL) - tanh((t - tR)/twR)) - epsilon
         S = S/(1 - epsilon)
      else
         epsilon = 0.5*(tanh((t0 - tL)/twL) + 1.0_f64)
         S = 0.5*(tanh((t - tL)/twL) + 1.0_f64) - epsilon
         S = S/(1.0_f64 - epsilon)
      end if
      if (S < 0) then
         S = 0.0_f64
      end if
      return
   end subroutine PF_envelope

end module sll_m_vp_cartesian_2d
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
