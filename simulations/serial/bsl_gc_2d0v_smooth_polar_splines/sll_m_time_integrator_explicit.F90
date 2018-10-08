module sll_m_time_integrator_explicit
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"

  use sll_m_working_precision, only: f64

  use sll_m_simulation_state, only: sll_t_simulation_state

  use sll_m_spline_interpolator_2d, only: sll_t_spline_interpolator_2d

  use sll_m_polar_mapping_iga, only: sll_t_polar_mapping_iga

  use sll_m_poisson_2d_fem_sps_stencil_new, only: sll_t_poisson_2d_fem_sps_stencil_new

  use sll_m_time_integrator_base, only: sll_c_time_integrator

  implicit none

  public :: sll_t_time_integrator_explicit

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! Working precision
  integer, parameter :: wp = f64

  type, extends(sll_c_time_integrator) :: sll_t_time_integrator_explicit

    real(wp), allocatable :: eta_copy(:,:,:)

    type(sll_t_simulation_state) :: sim_state_copy

  contains

    procedure :: init            => s_time_integrator_explicit__init
    procedure :: predictor       => s_time_integrator_explicit__predictor
    procedure :: corrector       => s_time_integrator_explicit__corrector
    procedure :: advance_in_time => s_time_integrator_explicit__advance_in_time
    procedure :: free            => s_time_integrator_explicit__free

  end type sll_t_time_integrator_explicit

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine s_time_integrator_explicit__init( &
    self            , &
    dt              , &
    tau_eta1        , &
    tau_eta2        , &
    mapping_discrete, &
    spline_interp_2d, &
    sim_state       , &
    poisson_solver )
    class(sll_t_time_integrator_explicit)             , intent(inout) :: self
    real(wp)                                          , intent(in   ) :: dt
    real(wp)                                  , target, intent(in   ) :: tau_eta1(:)
    real(wp)                                  , target, intent(in   ) :: tau_eta2(:)
    type(sll_t_polar_mapping_iga)                     , intent(in   ) :: mapping_discrete
    type(sll_t_spline_interpolator_2d)        , target, intent(in   ) :: spline_interp_2d
    type(sll_t_simulation_state)              , target, intent(in   ) :: sim_state
    type(sll_t_poisson_2d_fem_sps_stencil_new), target, intent(in   ) :: poisson_solver

    self % dt = dt

    self % tau_eta1 => tau_eta1
    self % tau_eta2 => tau_eta2

    call self % jacobian_2d_pseudo_cartesian % init( mapping_discrete )

    self % poisson_solver   => poisson_solver
    self % spline_interp_2d => spline_interp_2d

    allocate( self % eta_copy( 2, size( tau_eta1 ), size( tau_eta2 ) ) )

    self % sim_state => sim_state
    call self % sim_state % copy( self % sim_state_copy )

  end subroutine s_time_integrator_explicit__init

  !-----------------------------------------------------------------------------
  subroutine s_time_integrator_explicit__predictor( self )
    class(sll_t_time_integrator_explicit), intent(inout) :: self

    integer  :: i1, i2, ic
    real(wp) :: x(2), u(2)

    associate( dt                 => self % dt                            , &
               nc                 => self % sim_state % nc                , &
               ntau1              => size( self % tau_eta1 )              , &
               ntau2              => size( self % tau_eta2 )              , &
               sim_state          => self % sim_state                     , &
               spline_2d_rho      => self % sim_state % spline_2d_rho     , &
               rho_copy           => self % sim_state_copy % rho          , &
               point_charges      => self % sim_state      % point_charges, &
               point_charges_copy => self % sim_state_copy % point_charges, &
               eta_copy           => self % eta_copy )

      ! Evolve density

      !$OMP PARALLEL DO PRIVATE(x,u)
      do i2 = 1, ntau2
        do i1 = 1, ntau1

          eta_copy(1,i1,i2) = self % tau_eta1(i1)
          eta_copy(2,i1,i2) = self % tau_eta2(i2)

          x = self % logical_to_pseudo_cartesian( eta_copy(:,i1,i2) )
          u = self % advection_field( sim_state , eta_copy(:,i1,i2) )
          x = x - dt * u
          eta_copy(:,i1,i2) = self % pseudo_cartesian_to_logical( x )

          rho_copy(i1,i2) = spline_2d_rho % eval( eta_copy(1,i1,i2), eta_copy(2,i1,i2) )

        end do
      end do
      !$OMP END PARALLEL DO

      ! Apply periodicity along theta
      rho_copy(:,ntau2+1) = rho_copy(:,1)

      ! Evolve point charges
      if ( self % sim_state % point_charges_present ) then

        do ic = 1, nc

          point_charges_copy(ic) % location = point_charges(ic) % location

          x = self % logical_to_pseudo_cartesian( point_charges_copy(ic) % location )
          u = self % advection_field( sim_state, point_charges_copy(ic) % location )
          x = x + dt * u
          point_charges_copy(ic) % location = self % pseudo_cartesian_to_logical( x )

        end do

      end if

    end associate

  end subroutine s_time_integrator_explicit__predictor

  !-----------------------------------------------------------------------------
  subroutine s_time_integrator_explicit__corrector( self )
    class(sll_t_time_integrator_explicit), intent(inout) :: self

    integer  :: i1, i2, ic
    real(wp) :: eta(2), x(2), u1(2), u2(2)

    associate( dt                 => self % dt                            , &
               nc                 => self % sim_state % nc                , &
               ntau1              => size( self % tau_eta1 )              , &
               ntau2              => size( self % tau_eta2 )              , &
               sim_state          => self % sim_state                     , &
               sim_state_copy     => self % sim_state_copy                , &
               spline_2d_rho      => self % sim_state % spline_2d_rho     , &
               rho                => self % sim_state % rho               , &
               point_charges      => self % sim_state      % point_charges, &
               point_charges_copy => self % sim_state_copy % point_charges, &
               eta_copy           => self % eta_copy )

      ! Evolve density

      !$OMP PARALLEL DO PRIVATE(eta,x,u1,u2)
      do i2 = 1, ntau2
        do i1 = 1, ntau1

          eta(1) = self % tau_eta1(i1)
          eta(2) = self % tau_eta2(i2)

          x  = self % logical_to_pseudo_cartesian( eta )
          u1 = self % advection_field( sim_state_copy, eta )
          u2 = self % advection_field( sim_state, eta_copy(:,i1,i2) )
          x  = x - 0.5_wp * dt * ( u1 + u2 )
          eta = self % pseudo_cartesian_to_logical( x )

          rho(i1,i2) = spline_2d_rho % eval( eta(1), eta(2) )

        end do
      end do
      !$OMP END PARALLEL DO

      ! Apply periodicity along theta
      rho(:,ntau2+1) = rho(:,1)

      ! Evolve point charges
      if ( self % sim_state % point_charges_present ) then

        do ic = 1, nc

          x  = self % logical_to_pseudo_cartesian( point_charges(ic) % location )
          u1 = self % advection_field( sim_state_copy, point_charges_copy(ic) % location )
          u2 = self % advection_field( sim_state, point_charges(ic) % location )
          x  = x + 0.5_wp * dt * ( u1 + u2 )
          point_charges(ic) % location = self % pseudo_cartesian_to_logical( x )

        end do

      end if

    end associate

  end subroutine s_time_integrator_explicit__corrector

  !-----------------------------------------------------------------------------
  subroutine s_time_integrator_explicit__advance_in_time( self )
    class(sll_t_time_integrator_explicit), intent(inout) :: self

    integer :: ic

    associate( nc                 => self % sim_state % nc                , &
               ntau2              => size( self % sim_state % rho, 2 ) - 1, &
               spline_interp_2d   => self % spline_interp_2d              , &
               spline_2d_rho      => self % sim_state      % spline_2d_rho, &
               spline_2d_rho_copy => self % sim_state_copy % spline_2d_rho, &
               spline_2d_phi      => self % sim_state      % spline_2d_phi, &
               spline_2d_phi_copy => self % sim_state_copy % spline_2d_phi, &
               rho                => self % sim_state      % rho          , &
               rho_copy           => self % sim_state_copy % rho          , &
               point_charges      => self % sim_state      % point_charges, &
               point_charges_copy => self % sim_state_copy % point_charges )

      call self % predictor()

      call spline_interp_2d % compute_interpolant( spline_2d_rho_copy, rho_copy(:,1:ntau2) )

      ! Solve Poisson equation
      call self % poisson_solver % reset_charge()
      call self % poisson_solver % accumulate_charge ( spline_2d_rho_copy )
      if ( self % sim_state % point_charges_present ) then
        do ic = 1, nc
          call self % poisson_solver % accumulate_charge( &
            point_charges_copy(ic) % intensity, &
            point_charges_copy(ic) % location )
        end do
      end if
      call self % poisson_solver % solve ( spline_2d_phi_copy )

      call self % corrector()

      call spline_interp_2d % compute_interpolant( spline_2d_rho, rho(:,1:ntau2) )

      ! Solve Poisson equation
      call self % poisson_solver % reset_charge()
      call self % poisson_solver % accumulate_charge ( spline_2d_rho )
      if ( self % sim_state % point_charges_present ) then
        do ic = 1, nc
          call self % poisson_solver % accumulate_charge( &
            point_charges(ic) % intensity, &
            point_charges(ic) % location )
        end do
      end if
      call self % poisson_solver % solve ( spline_2d_phi )

    end associate

  end subroutine s_time_integrator_explicit__advance_in_time

  !-----------------------------------------------------------------------------
  subroutine s_time_integrator_explicit__free( self )
    class(sll_t_time_integrator_explicit), intent(inout) :: self

    nullify( self % tau_eta1 )
    nullify( self % tau_eta2 )

    nullify( self % poisson_solver )
    nullify( self % spline_interp_2d )

    deallocate( self % eta_copy )

    call self % sim_state_copy % free()

  end subroutine s_time_integrator_explicit__free

end module sll_m_time_integrator_explicit
