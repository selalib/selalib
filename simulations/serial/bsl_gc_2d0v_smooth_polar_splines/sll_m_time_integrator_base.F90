module sll_m_time_integrator_base
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"

  use sll_m_working_precision, only: f64

  use sll_m_constants, only: sll_p_twopi

  use sll_m_simulation_state, only: sll_t_simulation_state

  use sll_m_point_charge, only: sll_t_point_charge

  use sll_m_spline_interpolator_2d, only: sll_t_spline_interpolator_2d

  use sll_m_polar_mapping_iga, only: sll_t_polar_mapping_iga

  use sll_m_jacobian_2d_pseudo_cartesian, only: sll_t_jacobian_2d_pseudo_cartesian

  use sll_m_poisson_2d_fem_sps_stencil_new, only: sll_t_poisson_2d_fem_sps_stencil_new

	implicit none

	public :: sll_c_time_integrator

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! Working precision
  integer, parameter :: wp = f64

  type, abstract :: sll_c_time_integrator

    ! Time step
    real(wp) :: dt

    ! Interpolation points
    real(wp), pointer :: tau_eta1(:) => null()
    real(wp), pointer :: tau_eta2(:) => null()

    type(sll_t_jacobian_2d_pseudo_cartesian) :: jacobian_2d_pseudo_cartesian

    type(sll_t_poisson_2d_fem_sps_stencil_new), pointer :: poisson_solver   => null()
    type(sll_t_spline_interpolator_2d)        , pointer :: spline_interp_2d => null()

    type(sll_t_simulation_state), pointer :: sim_state

  contains

    ! Deferred procedures
    procedure(i_sub__advect) , deferred :: advect_single_coords
    procedure(i_sub__advance), deferred :: predictor
    procedure(i_sub__advance), deferred :: corrector
    procedure(i_sub__advance), deferred :: advance_in_time
    procedure(i_sub__free)   , deferred :: free

    ! Non-deferred procedures
    procedure :: advection_field             => f_time_integrator__advection_field
    procedure :: advect_distribution         => s_time_integrator__advect_distribution
    procedure :: advect_point_charges        => s_time_integrator__advect_point_charges
    procedure :: logical_to_pseudo_cartesian => f_time_integrator__logical_to_pseudo_cartesian
    procedure :: pseudo_cartesian_to_logical => f_time_integrator__pseudo_cartesian_to_logical

  end type sll_c_time_integrator

  abstract interface

    subroutine i_sub__advect( self, h, success, eta )
      import wp
      import sll_c_time_integrator
      class(sll_c_time_integrator), intent(in   ) :: self
      real(wp)                    , intent(in   ) :: h
      logical                     , intent(  out) :: success
      real(wp)                    , intent(inout) :: eta(2)
    end subroutine i_sub__advect

    subroutine i_sub__advance( self, success )
      import sll_c_time_integrator
      class(sll_c_time_integrator), intent(inout) :: self
      logical                     , intent(  out) :: success
    end subroutine i_sub__advance

    subroutine i_sub__free( self )
      import sll_c_time_integrator
      class(sll_c_time_integrator), intent(inout) :: self
    end subroutine i_sub__free

  end interface

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  SLL_PURE function f_time_integrator__advection_field( self, eta ) result( u )
    class(sll_c_time_integrator), intent(in) :: self
    real(wp)                    , intent(in) :: eta(2)
    real(wp) :: u(2)

    real(wp) :: jmat_comp(2,2), E(2)

    jmat_comp = self % jacobian_2d_pseudo_cartesian % eval( eta )

    ! Cartesian components of electric field
    E = self % sim_state % electric_field % eval( eta )

    ! Pseudo-Cartesian components of advection field
    u(1) = jmat_comp(1,1) * (-E(2)) + jmat_comp(1,2) * E(1)
    u(2) = jmat_comp(2,1) * (-E(2)) + jmat_comp(2,2) * E(1)

  end function f_time_integrator__advection_field

  !-----------------------------------------------------------------------------
  subroutine s_time_integrator__advect_distribution( self, h, success, rho_new )
    class(sll_c_time_integrator), intent(inout) :: self
    real(wp)                    , intent(in   ) :: h
    logical                     , intent(  out) :: success
    real(wp)                    , intent(inout) :: rho_new(:,:)

    integer  :: i1, i2, ntau1, ntau2
    real(wp) :: eta(2)

    ntau1 = size( self % tau_eta1 )
    ntau2 = size( self % tau_eta2 )

    success = .true.

    !$OMP PARALLEL DO PRIVATE(eta,success)
    do i2 = 1, ntau2
      do i1 = 1, ntau1

        eta(1) = self % tau_eta1(i1)
        eta(2) = self % tau_eta2(i2)

        call self % advect_single_coords( h, success, eta )

        ! Check if integrator converged
        if ( success ) then
          rho_new(i1,i2) = self % sim_state % spline_2d_rho % eval( eta(1), eta(2) )
        else
          exit
        end if

      end do
    end do
    !$OMP END PARALLEL DO

    ! Apply periodicity along theta
    if ( success ) rho_new(:,ntau2+1) = rho_new(:,1)

  end subroutine s_time_integrator__advect_distribution

  !-----------------------------------------------------------------------------
  subroutine s_time_integrator__advect_point_charges( self, h, success, point_charges )
    class(sll_c_time_integrator), intent(inout) :: self
    real(wp)                    , intent(in   ) :: h
    logical                     , intent(  out) :: success
    type(sll_t_point_charge)    , intent(inout) :: point_charges(:)

    integer :: ic

    SLL_ASSERT( size( point_charges ) > 0 )

    associate( nc => size( point_charges ) )

      do ic = 1, nc
        call self % advect_single_coords( h, success, point_charges(ic) % location )
        if ( .not. success ) exit
      end do

    end associate

  end subroutine s_time_integrator__advect_point_charges

  !-----------------------------------------------------------------------------
  SLL_PURE function f_time_integrator__logical_to_pseudo_cartesian( self, eta ) result( x )
    class(sll_c_time_integrator), intent(in) :: self
    real(wp)                    , intent(in) :: eta(2)
    real(wp) :: x(2)

    x(1) = eta(1) * cos( eta(2) )
    x(2) = eta(1) * sin( eta(2) )

  end function f_time_integrator__logical_to_pseudo_cartesian

  !-----------------------------------------------------------------------------
  SLL_PURE function f_time_integrator__pseudo_cartesian_to_logical( self, x ) result( eta )
    class(sll_c_time_integrator), intent(in) :: self
    real(wp)                    , intent(in) :: x(2)
    real(wp) :: eta(2)

    eta(1) = norm2( x )
    eta(2) = modulo( atan2( x(2), x(1) ), sll_p_twopi ) ! atan2 returns theta in [-pi,pi)

    ! For characteristics going out of the domain
    if ( eta(1) > 1.0_wp ) eta(1) = 1.0_wp

  end function f_time_integrator__pseudo_cartesian_to_logical

end module sll_m_time_integrator_base
