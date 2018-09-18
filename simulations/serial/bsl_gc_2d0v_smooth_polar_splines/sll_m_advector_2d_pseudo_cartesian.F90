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

  use sll_m_spline_2d, only: sll_t_spline_2d

  use sll_m_electric_field, only: sll_t_electric_field

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

    type(sll_t_polar_mapping_iga), pointer :: mapping_discrete => null()
    type(sll_t_spline_2d)        , pointer :: spline_2d_rho    => null()
    type(sll_t_electric_field)   , pointer :: electric_field   => null()

    type(sll_t_jacobian_2d_pseudo_cartesian) :: jacobian_2d_pseudo_cartesian

  contains

    procedure :: init             => s_advector_2d_pseudo_cartesian__init
    procedure :: advect           => s_advector_2d_pseudo_cartesian__advect
    procedure :: advance_position => s_advector_2d_pseudo_cartesian__advance_position
    procedure :: free             => s_advector_2d_pseudo_cartesian__free

  end type sll_t_advector_2d_pseudo_cartesian

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine s_advector_2d_pseudo_cartesian__init( &
    self            , &
    tau_eta1        , &
    tau_eta2        , &
    mapping_discrete, &
    spline_2d_rho   , &
    electric_field  , &
    abs_tol         , &
    rel_tol         , &
    maxiter )
    class(sll_t_advector_2d_pseudo_cartesian), intent(inout) :: self
    real(wp)                     , target    , intent(in   ) :: tau_eta1(:)
    real(wp)                     , target    , intent(in   ) :: tau_eta2(:)
    type(sll_t_polar_mapping_iga), target    , intent(in   ) :: mapping_discrete
    type(sll_t_spline_2d)        , target    , intent(in   ) :: spline_2d_rho
    type(sll_t_electric_field)   , target    , intent(in   ) :: electric_field
    real(wp)                                 , intent(in   ) :: abs_tol
    real(wp)                                 , intent(in   ) :: rel_tol
    integer                                  , intent(in   ) :: maxiter

    self % tau_eta1 => tau_eta1
    self % tau_eta2 => tau_eta2

    self % mapping_discrete => mapping_discrete
    self % spline_2d_rho    => spline_2d_rho
    self % electric_field   => electric_field

    call self % jacobian_2d_pseudo_cartesian % init( mapping_discrete )

    self % abs_tol = abs_tol
    self % rel_tol = rel_tol
    self % maxiter = maxiter

  end subroutine s_advector_2d_pseudo_cartesian__init

  !-----------------------------------------------------------------------------
  subroutine s_advector_2d_pseudo_cartesian__advect( self, h, success, rho_new )
    class(sll_t_advector_2d_pseudo_cartesian), intent(inout) :: self
    real(wp)                                 , intent(in   ) :: h
    logical                                  , intent(  out) :: success
    real(wp)                                 , intent(inout) :: rho_new(:,:)

    integer  :: i1, i2, ntau1, ntau2
    real(wp) :: eta(2)

    character(len=*), parameter :: this_sub_name = "sll_t_advector_2d_pseudo_cartesian % advect"
    character(len=256) :: err_msg

    ntau1 = size( self % tau_eta1 )
    ntau2 = size( self % tau_eta2 )

    success = .true.

    !$OMP PARALLEL DO PRIVATE(eta,success)
    do i2 = 1, ntau2
      do i1 = 1, ntau1

        eta(1) = self % tau_eta1(i1)
        eta(2) = self % tau_eta2(i2)

        call self % advance_position( h, success, eta )

        ! Check if integrator converged
        if ( success ) then
          rho_new(i1,i2) = self % spline_2d_rho % eval( eta(1), eta(2) )
        else
          exit
        end if

      end do
    end do
    !$OMP END PARALLEL DO

    ! Apply periodicity along theta
    if ( success ) rho_new(:,ntau2+1) = rho_new(:,1)

  end subroutine s_advector_2d_pseudo_cartesian__advect

  !-----------------------------------------------------------------------------
  recursive subroutine s_advector_2d_pseudo_cartesian__advance_position( self, h, success, eta )
    class(sll_t_advector_2d_pseudo_cartesian), intent(inout) :: self
    real(wp)                                 , intent(in   ) :: h
    logical                                  , intent(  out) :: success
    real(wp)                                 , intent(inout) :: eta(2)

    integer  :: j
    real(wp) :: tol_sqr, eta0(2), x0(2), x(2), dx(2), temp(2), k2(2), a0(2), E(2), jmat_comp(2,2)

    character(len=*), parameter :: this_sub_name = "sll_t_advector_2d_pseudo_cartesian % advance_position"
    character(len=256) :: err_msg

    success = .true.

    eta0 = eta

    x0 = logical_to_pseudo_cartesian( eta )

    tol_sqr = ( self % abs_tol + self % rel_tol * norm2( x0 ) )**2

    ! Cartesian components of electric field
    E = self % electric_field % eval( eta )

    jmat_comp = self % jacobian_2d_pseudo_cartesian % eval( eta )

    ! Pseudo-Cartesian components of advection field
    a0(1) = jmat_comp(1,1) * (-E(2)) + jmat_comp(1,2) * E(1)
    a0(2) = jmat_comp(2,1) * (-E(2)) + jmat_comp(2,2) * E(1)

    ! First iteration
    dx  = h * a0
    x   = x0 + dx
    eta = pseudo_cartesian_to_logical( x )

    ! Successive iterations if first iteration did not converge
    if ( dot_product( dx, dx ) > tol_sqr ) then

      success = .false.

      associate( k1 => a0, h_half => 0.5_wp*h, dx_old => temp, error => temp )

        do j = 2, self % maxiter

          jmat_comp = self % jacobian_2d_pseudo_cartesian % eval( eta )

          ! Cartesian components of advection field
          E = self % electric_field % eval( eta )

          ! k2 = f(t,x_{i-1})
          k2(1) = jmat_comp(1,1) * (-E(2)) + jmat_comp(1,2) * E(1)
          k2(2) = jmat_comp(2,1) * (-E(2)) + jmat_comp(2,2) * E(1)

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

!      write( err_msg, '(a,i0,a)' ) "integration of characteristics did not converge after ", self % maxiter, " iterations: trying sub-stepping"
!      SLL_WARNING( this_sub_name, err_msg )
!      return

      ! Sub-stepping
      eta = eta0
      call self % advance_position( h*0.5_wp, success, eta )
      call self % advance_position( h*0.5_wp, success, eta )
      return

    end if

  end subroutine s_advector_2d_pseudo_cartesian__advance_position

  !-----------------------------------------------------------------------------
  subroutine s_advector_2d_pseudo_cartesian__free( self )
    class(sll_t_advector_2d_pseudo_cartesian), intent(inout) :: self

    nullify( self % tau_eta1 )
    nullify( self % tau_eta2 )

    nullify( self % mapping_discrete )
    nullify( self % spline_2d_rho    )
    nullify( self % electric_field   )

    call self % jacobian_2d_pseudo_cartesian % free()

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
