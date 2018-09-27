module sll_m_advector_2d_pseudo_cartesian
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_errors.h"

  use sll_m_working_precision, only: f64

  use sll_m_constants, only: &
    sll_p_pi, &
    sll_p_twopi

  use sll_m_polar_mapping_iga, only: sll_t_polar_mapping_iga

  use sll_m_jacobian_2d_pseudo_cartesian, only: sll_t_jacobian_2d_pseudo_cartesian

  use sll_m_point_charge, only: sll_t_point_charge

  use sll_m_simulation_state, only: sll_t_simulation_state

  implicit none

  public :: sll_t_advector_2d_pseudo_cartesian

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! Working precision
  integer, parameter :: wp = f64

  type :: sll_t_advector_2d_pseudo_cartesian

    real(wp) :: abs_tol
    real(wp) :: rel_tol
    integer  :: maxiter

    real(wp), pointer :: tau_eta1(:)
    real(wp), pointer :: tau_eta2(:)

    type(sll_t_polar_mapping_iga), pointer   :: mapping_discrete => null()
    type(sll_t_jacobian_2d_pseudo_cartesian) :: jacobian_2d_pseudo_cartesian

    type(sll_t_simulation_state), pointer :: sim_state

  contains

    procedure :: init                 => s_advector_2d_pseudo_cartesian__init
    procedure :: advect_field         => f_advector_2d_pseudo_cartesian__advect_field
    procedure :: advect_single_coords => s_advector_2d_pseudo_cartesian__advect_single_coords
    procedure :: advect_distribution  => s_advector_2d_pseudo_cartesian__advect_distribution
    procedure :: advect_point_charges => s_advector_2d_pseudo_cartesian__advect_point_charges
    procedure :: free                 => s_advector_2d_pseudo_cartesian__free

  end type sll_t_advector_2d_pseudo_cartesian

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine s_advector_2d_pseudo_cartesian__init( &
    self            , &
    tau_eta1        , &
    tau_eta2        , &
    mapping_discrete, &
    sim_state       , &
    abs_tol         , &
    rel_tol         , &
    maxiter )
    class(sll_t_advector_2d_pseudo_cartesian), intent(inout) :: self
    real(wp)                     , target    , intent(in   ) :: tau_eta1(:)
    real(wp)                     , target    , intent(in   ) :: tau_eta2(:)
    type(sll_t_polar_mapping_iga), target    , intent(in   ) :: mapping_discrete
    type(sll_t_simulation_state) , target    , intent(in   ) :: sim_state
    real(wp)                                 , intent(in   ) :: abs_tol
    real(wp)                                 , intent(in   ) :: rel_tol
    integer                                  , intent(in   ) :: maxiter

    self % tau_eta1 => tau_eta1
    self % tau_eta2 => tau_eta2

    self % mapping_discrete => mapping_discrete
    call self % jacobian_2d_pseudo_cartesian % init( mapping_discrete )

    self % sim_state => sim_state

    self % abs_tol = abs_tol
    self % rel_tol = rel_tol
    self % maxiter = maxiter

  end subroutine s_advector_2d_pseudo_cartesian__init

  !-----------------------------------------------------------------------------
  SLL_PURE function f_advector_2d_pseudo_cartesian__advect_field( self, eta ) result( u )
    class(sll_t_advector_2d_pseudo_cartesian), intent(in) :: self
    real(wp)                                 , intent(in) :: eta(2)
    real(wp) :: u(2)

    real(wp) :: jmat_comp(2,2), E(2)

    jmat_comp = self % jacobian_2d_pseudo_cartesian % eval( eta )

    ! Cartesian components of electric field
    E = self % sim_state % electric_field % eval( eta )

    ! Pseudo-Cartesian components of advection field
    u(1) = jmat_comp(1,1) * (-E(2)) + jmat_comp(1,2) * E(1)
    u(2) = jmat_comp(2,1) * (-E(2)) + jmat_comp(2,2) * E(1)

  end function f_advector_2d_pseudo_cartesian__advect_field

  !-----------------------------------------------------------------------------
  recursive subroutine s_advector_2d_pseudo_cartesian__advect_single_coords( self, h, success, eta )
    class(sll_t_advector_2d_pseudo_cartesian), intent(in   ) :: self
    real(wp)                                 , intent(in   ) :: h
    logical                                  , intent(  out) :: success
    real(wp)                                 , intent(inout) :: eta(2)

    integer  :: j
    real(wp) :: tol_sqr, eta0(2), x0(2), x(2), dx(2), temp(2), k2(2), a0(2)

    character(len=*), parameter :: this_sub_name = "sll_t_advector_2d_pseudo_cartesian % advect_single_coords"
    character(len=256) :: err_msg

    success = .true.

    eta0 = eta

    x0 = logical_to_pseudo_cartesian( eta )

    tol_sqr = ( self % abs_tol + self % rel_tol * norm2( x0 ) )**2

    ! Pseudo-Cartesian components of advection field
    a0 = self % advect_field( eta )

    ! First iteration
    dx  = h * a0
    x   = x0 + dx
    eta = pseudo_cartesian_to_logical( x )

    ! Successive iterations if first iteration did not converge
    if ( dot_product( dx, dx ) > tol_sqr ) then

      success = .false.

      associate( k1 => a0, h_half => 0.5_wp*h, dx_old => temp, error => temp )

        do j = 2, self % maxiter

          ! k2 = f(t,x_{i-1})
          k2 = self % advect_field( eta )

          dx_old = dx
          dx     = h_half*(k1+k2)
          error  = dx_old - dx
          x      = x0 + dx
          eta    = pseudo_cartesian_to_logical( x )

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
      call self % advect_single_coords( h*0.5_wp, success, eta )
      call self % advect_single_coords( h*0.5_wp, success, eta )
      return

    end if

  end subroutine s_advector_2d_pseudo_cartesian__advect_single_coords

  !-----------------------------------------------------------------------------
  subroutine s_advector_2d_pseudo_cartesian__advect_distribution( self, h, success, rho_new )
    class(sll_t_advector_2d_pseudo_cartesian), intent(inout) :: self
    real(wp)                                 , intent(in   ) :: h
    logical                                  , intent(  out) :: success
    real(wp)                                 , intent(inout) :: rho_new(:,:)

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

  end subroutine s_advector_2d_pseudo_cartesian__advect_distribution

  !-----------------------------------------------------------------------------
  subroutine s_advector_2d_pseudo_cartesian__advect_point_charges( self, h, success, point_charges )
    class(sll_t_advector_2d_pseudo_cartesian), intent(inout) :: self
    real(wp)                                 , intent(in   ) :: h
    logical                                  , intent(  out) :: success
    type(sll_t_point_charge)                 , intent(inout) :: point_charges(:)

    integer :: ic

    SLL_ASSERT( size( point_charges ) > 0 )

    associate( nc => size( point_charges ) )

      do ic = 1, nc
        call self % advect_single_coords( h, success, point_charges(ic) % location )
        if ( .not. success ) exit
      end do

    end associate

  end subroutine s_advector_2d_pseudo_cartesian__advect_point_charges

  !-----------------------------------------------------------------------------
  subroutine s_advector_2d_pseudo_cartesian__free( self )
    class(sll_t_advector_2d_pseudo_cartesian), intent(inout) :: self

    nullify( self % tau_eta1 )
    nullify( self % tau_eta2 )

    nullify( self % mapping_discrete )
    call self % jacobian_2d_pseudo_cartesian % free()

    nullify( self % sim_state )

  end subroutine s_advector_2d_pseudo_cartesian__free

  !-----------------------------------------------------------------------------
  SLL_PURE function logical_to_pseudo_cartesian( eta ) result( x )
    real(wp), intent(in) :: eta(2)
    real(wp) :: x(2)

    x(1) = eta(1) * cos( eta(2) )
    x(2) = eta(1) * sin( eta(2) )

  end function logical_to_pseudo_cartesian

  !-----------------------------------------------------------------------------
  SLL_PURE function pseudo_cartesian_to_logical( x ) result( eta )
    real(wp), intent(in) :: x(2)
    real(wp) :: eta(2)

    eta(1) = norm2( x )
    eta(2) = modulo( atan2( x(2), x(1) ), sll_p_twopi ) ! atan2 returns theta in [-pi,pi)

    ! For characteristics going out of the domain
    if ( eta(1) > 1.0_wp ) eta(1) = 1.0_wp

  end function pseudo_cartesian_to_logical

end module sll_m_advector_2d_pseudo_cartesian
