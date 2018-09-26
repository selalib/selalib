module sll_m_time_integrator
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_errors.h"

  use sll_m_working_precision, only: f64

  use sll_m_simulation_state, only: sll_t_simulation_state

  use sll_m_advector_2d_pseudo_cartesian, only: sll_t_advector_2d_pseudo_cartesian

  use sll_m_poisson_2d_fem_sps_stencil_new, only: sll_t_poisson_2d_fem_sps_stencil_new

  use sll_m_spline_interpolator_2d, only: sll_t_spline_interpolator_2d

	implicit none

	public :: sll_t_time_integrator

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! Working precision
  integer, parameter :: wp = f64

  type :: sll_t_time_integrator

    real(wp) :: dt

    type(sll_t_advector_2d_pseudo_cartesian)  , pointer :: advector       => null()
    type(sll_t_poisson_2d_fem_sps_stencil_new), pointer :: poisson_solver => null()
    type(sll_t_spline_interpolator_2d)        , pointer :: spline_interp_2d => null()

    type(sll_t_simulation_state), pointer :: sim_state
    type(sll_t_simulation_state)          :: sim_state_aux

  contains

    procedure :: init            => s_time_integrator__init
    procedure :: advance_in_time => s_time_integrator__advance_in_time
    procedure :: free            => s_time_integrator__free

  end type sll_t_time_integrator

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine s_time_integrator__init( self, dt, advector, poisson_solver, spline_interp_2d, sim_state )
    class(sll_t_time_integrator)                      , intent(inout) :: self
    real(wp)                                          , intent(in   ) :: dt
    type(sll_t_advector_2d_pseudo_cartesian)  , target, intent(in   ) :: advector
    type(sll_t_poisson_2d_fem_sps_stencil_new), target, intent(in   ) :: poisson_solver
    type(sll_t_spline_interpolator_2d)        , target, intent(in   ) :: spline_interp_2d
    type(sll_t_simulation_state)              , target, intent(in   ) :: sim_state

    self % dt = dt

    self % advector         => advector
    self % poisson_solver   => poisson_solver
    self % spline_interp_2d => spline_interp_2d

    self % sim_state => sim_state
    call self % sim_state % copy( self % sim_state_aux )

  end subroutine s_time_integrator__init

  !-----------------------------------------------------------------------------
  subroutine s_time_integrator__advance_in_time( self, sim_state, success )
    class(sll_t_time_integrator), intent(inout) :: self
    type(sll_t_simulation_state), intent(inout) :: sim_state
    logical                     , intent(  out) :: success

    integer :: ic

    associate( dt                => self % dt                           , &
               nc                => size( sim_state % point_charges )   , &
               ntau2             => size( sim_state % rho, 2 ) - 1      , &
               spline_interp_2d  => self % spline_interp_2d             , &
               spline_2d_rho     => sim_state % spline_2d_rho           , &
               spline_2d_phi     => sim_state % spline_2d_phi           , &
               point_charges     => sim_state % point_charges           , &
               point_charges_aux => self % sim_state_aux % point_charges, &
               rho_aux           => self % sim_state_aux % rho )

      ! Predictor: evolve density
      call self % advector % advect_distribution( -0.5_wp*dt, success, rho_aux )
      if ( .not. success ) return

      ! Predictor: evolve point charges
      if ( nc /= 0 ) then
        do ic = 1, nc
          point_charges_aux(ic) % location(:) = point_charges(ic) % location(:)
        end do
        call self % advector % advect_point_charges( 0.5_wp*dt, success, point_charges )
        if ( .not. success ) return
      end if

      call spline_interp_2d % compute_interpolant( spline_2d_rho, rho_aux(:,1:ntau2) )

      ! Solve Poisson equation
      call self % poisson_solver % reset_charge()
      call self % poisson_solver % accumulate_charge ( spline_2d_rho )
      if ( nc /= 0 ) then
        do ic = 1, nc
          associate( intensity => point_charges(ic)%intensity, &
                     location  => point_charges(ic)%location )
            call self % poisson_solver % accumulate_charge( intensity, location )
          end associate
        end do
      end if
      call self % poisson_solver % solve ( spline_2d_phi )

      call spline_interp_2d % compute_interpolant( spline_2d_rho, sim_state % rho(:,1:ntau2) )

      ! Corrector: evolve density
      call self % advector % advect_distribution( -dt, success, rho_aux )
      if ( .not. success ) return

      ! Corrector: evolve point charges
      if ( nc /= 0 ) then
        do ic = 1, nc
          point_charges(ic) % location(:) = point_charges_aux(ic) % location(:)
        end do
        call self % advector % advect_point_charges( dt, success, point_charges )
        if ( .not. success ) return
      end if

      call spline_interp_2d % compute_interpolant( spline_2d_rho, rho_aux(:,1:ntau2) )

      ! Solve Poisson equation
      call self % poisson_solver % reset_charge()
      call self % poisson_solver % accumulate_charge ( spline_2d_rho )
      if ( nc /= 0 ) then
        do ic = 1, nc
          associate( intensity => point_charges(ic)%intensity, &
                     location  => point_charges(ic)%location )
            call self % poisson_solver % accumulate_charge( intensity, location )
          end associate
        end do
      end if
      call self % poisson_solver % solve ( spline_2d_phi )

      sim_state % rho = rho_aux

    end associate

  end subroutine s_time_integrator__advance_in_time

  !-----------------------------------------------------------------------------
  subroutine s_time_integrator__free( self )
    class(sll_t_time_integrator), intent(inout) :: self

    nullify( self % advector )
    nullify( self % poisson_solver )
    nullify( self % spline_interp_2d )

  end subroutine s_time_integrator__free

end module sll_m_time_integrator
