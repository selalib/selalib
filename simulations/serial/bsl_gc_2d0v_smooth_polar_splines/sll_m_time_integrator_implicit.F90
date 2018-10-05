module sll_m_time_integrator_implicit
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"

  use sll_m_working_precision, only: f64

  use sll_m_simulation_state, only: sll_t_simulation_state

  use sll_m_spline_interpolator_2d, only: sll_t_spline_interpolator_2d

  use sll_m_polar_mapping_iga, only: sll_t_polar_mapping_iga

  use sll_m_point_charge, only: sll_t_point_charge

  use sll_m_poisson_2d_fem_sps_stencil_new, only: sll_t_poisson_2d_fem_sps_stencil_new

  use sll_m_time_integrator_base, only: sll_c_time_integrator

  implicit none

  public :: sll_t_time_integrator_implicit

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! Working precision
  integer, parameter :: wp = f64

  type, extends(sll_c_time_integrator) :: sll_t_time_integrator_implicit

    real(wp) :: abs_tol
    real(wp) :: rel_tol
    integer  :: maxiter

    type(sll_t_simulation_state) :: sim_state_copy

  contains

    procedure :: init                 => s_time_integrator_implicit__init
    procedure :: advect_single_coords => s_time_integrator_implicit__advect_single_coords
    procedure :: predictor            => s_time_integrator_implicit__predictor
    procedure :: corrector            => s_time_integrator_implicit__corrector
    procedure :: advance_in_time      => s_time_integrator_implicit__advance_in_time
    procedure :: free                 => s_time_integrator_implicit__free

  end type sll_t_time_integrator_implicit

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine s_time_integrator_implicit__init( &
    self            , &
    dt              , &
    tau_eta1        , &
    tau_eta2        , &
    mapping_discrete, &
    spline_interp_2d, &
    sim_state       , &
    poisson_solver  , &
    abs_tol         , &
    rel_tol         , &
    maxiter )
    class(sll_t_time_integrator_implicit)             , intent(inout) :: self
    real(wp)                                          , intent(in   ) :: dt
    real(wp)                                  , target, intent(in   ) :: tau_eta1(:)
    real(wp)                                  , target, intent(in   ) :: tau_eta2(:)
    type(sll_t_polar_mapping_iga)                     , intent(in   ) :: mapping_discrete
    type(sll_t_spline_interpolator_2d)        , target, intent(in   ) :: spline_interp_2d
    type(sll_t_simulation_state)              , target, intent(in   ) :: sim_state
    type(sll_t_poisson_2d_fem_sps_stencil_new), target, intent(in   ) :: poisson_solver
    real(wp)                                          , intent(in   ) :: abs_tol
    real(wp)                                          , intent(in   ) :: rel_tol
    integer                                           , intent(in   ) :: maxiter

    self % dt = dt

    self % tau_eta1 => tau_eta1
    self % tau_eta2 => tau_eta2

    call self % jacobian_2d_pseudo_cartesian % init( mapping_discrete )

    self % poisson_solver   => poisson_solver
    self % spline_interp_2d => spline_interp_2d

    self % sim_state => sim_state
    call self % sim_state % copy( self % sim_state_copy )

    self % abs_tol = abs_tol
    self % rel_tol = rel_tol
    self % maxiter = maxiter

  end subroutine s_time_integrator_implicit__init

  !-----------------------------------------------------------------------------
  recursive subroutine s_time_integrator_implicit__advect_single_coords( self, h, sim_state, success, eta )
    class(sll_t_time_integrator_implicit), intent(in   ) :: self
    real(wp)                             , intent(in   ) :: h
    type(sll_t_simulation_state)         , intent(in   ) :: sim_state
    logical                              , intent(  out) :: success
    real(wp)                             , intent(inout) :: eta(2)

    integer  :: j
    real(wp) :: tol_sqr, eta0(2), x0(2), x(2), dx(2), temp(2), k2(2), a0(2)

    character(len=*), parameter :: this_sub_name = "sll_t_time_integrator_implicit % advect_single_coords"

    success = .true.

    eta0 = eta

    x0 = self % logical_to_pseudo_cartesian( eta )

    tol_sqr = ( self % abs_tol + self % rel_tol * norm2( x0 ) )**2

    ! Pseudo-Cartesian components of advection field
    a0 = self % advection_field( sim_state, eta )

    ! First iteration
    dx  = h * a0
    x   = x0 + dx
    eta = self % pseudo_cartesian_to_logical( x )

    ! Successive iterations if first iteration did not converge
    if ( dot_product( dx, dx ) > tol_sqr ) then

      success = .false.

      associate( k1 => a0, h_half => 0.5_wp*h, dx_old => temp, error => temp )

        do j = 2, self % maxiter

          ! k2 = f(t,x_{i-1})
          k2 = self % advection_field( sim_state, eta )

          dx_old = dx
          dx     = h_half*(k1+k2)
          error  = dx_old - dx
          x      = x0 + dx
          eta    = self % pseudo_cartesian_to_logical( x )

          if ( dot_product( error, error ) <= tol_sqr ) then
            success = .true.
            return
          end if

        end do

      end associate

    end if

    ! Check if integrator converged
    if ( .not. success ) then

      ! Sub-stepping
      eta = eta0
      call self % advect_single_coords( h*0.5_wp, sim_state, success, eta )
      call self % advect_single_coords( h*0.5_wp, sim_state, success, eta )
      return

    end if

  end subroutine s_time_integrator_implicit__advect_single_coords

  !-----------------------------------------------------------------------------
  subroutine s_time_integrator_implicit__predictor( self, success )
    class(sll_t_time_integrator_implicit), intent(inout) :: self
    logical                              , intent(  out) :: success

    integer  :: i1, i2, ic
    real(wp) :: eta(2)

    associate( dt                 => self % dt                            , &
               nc                 => self % sim_state % nc                , &
               ntau1              => size( self % tau_eta1 )              , &
               ntau2              => size( self % tau_eta2 )              , &
               sim_state          => self % sim_state                     , &
               spline_2d_rho      => self % sim_state % spline_2d_rho     , &
               rho_copy           => self % sim_state_copy % rho          , &
               point_charges      => self % sim_state      % point_charges, &
               point_charges_copy => self % sim_state_copy % point_charges )

      ! Evolve density
      do i2 = 1, ntau2
        do i1 = 1, ntau1

          eta(1) = self % tau_eta1(i1)
          eta(2) = self % tau_eta2(i2)

          call self % advect_single_coords( -0.5_wp*dt, sim_state, success, eta )

          ! Check if integrator converged
          if ( success ) then
            rho_copy(i1,i2) = spline_2d_rho % eval( eta(1), eta(2) )
          else
            return
          end if

        end do
      end do
      ! Apply periodicity along theta
      rho_copy(:,ntau2+1) = rho_copy(:,1)

      ! Evolve point charges
      if ( self % sim_state % point_charges_present ) then

        do ic = 1, nc
          point_charges_copy(ic) % location = point_charges(ic) % location
          call self % advect_single_coords( 0.5_wp*dt, sim_state, success, point_charges_copy(ic) % location )
          if ( .not. success ) return
        end do

      end if

    end associate

  end subroutine s_time_integrator_implicit__predictor

  !-----------------------------------------------------------------------------
  subroutine s_time_integrator_implicit__corrector( self, success )
    class(sll_t_time_integrator_implicit), intent(inout) :: self
    logical                              , intent(  out) :: success

    integer  :: i1, i2, ic
    real(wp) :: eta(2)

    associate( dt             => self % dt                       , &
               nc             => self % sim_state % nc           , &
               ntau1          => size( self % tau_eta1 )         , &
               ntau2          => size( self % tau_eta2 )         , &
               sim_state_copy => self % sim_state_copy           , &
               spline_2d_rho  => self % sim_state % spline_2d_rho, &
               rho            => self % sim_state % rho          , &
               point_charges  => self % sim_state % point_charges )

      ! Evolve density
      do i2 = 1, ntau2
        do i1 = 1, ntau1

          eta(1) = self % tau_eta1(i1)
          eta(2) = self % tau_eta2(i2)

          call self % advect_single_coords( -dt, sim_state_copy, success, eta )

          ! Check if integrator converged
          if ( success ) then
            rho(i1,i2) = spline_2d_rho % eval( eta(1), eta(2) )
          else
            return
          end if

        end do
      end do
      ! Apply periodicity along theta
      rho(:,ntau2+1) = rho(:,1)

      ! Evolve point charges
      if ( self % sim_state % point_charges_present ) then

        do ic = 1, nc
          call self % advect_single_coords( dt, sim_state_copy, success, point_charges(ic) % location )
          if ( .not. success ) return
        end do

      end if

    end associate

  end subroutine s_time_integrator_implicit__corrector

  !-----------------------------------------------------------------------------
  subroutine s_time_integrator_implicit__advance_in_time( self, success )
    class(sll_t_time_integrator_implicit), intent(inout) :: self
    logical                              , intent(  out) :: success

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

      call self % predictor( success )
      if ( .not. success ) return

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

      call self % corrector( success )
      if ( .not. success ) return

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

  end subroutine s_time_integrator_implicit__advance_in_time

  !-----------------------------------------------------------------------------
  subroutine s_time_integrator_implicit__free( self )
    class(sll_t_time_integrator_implicit), intent(inout) :: self

    nullify( self % tau_eta1 )
    nullify( self % tau_eta2 )

    nullify( self % poisson_solver )
    nullify( self % spline_interp_2d )

    call self % sim_state_copy % free()

  end subroutine s_time_integrator_implicit__free

end module sll_m_time_integrator_implicit
