!----------------------
! A note about notation
!----------------------
!  rkXY:
!       e = explicit
!       i = implicit (fully)
!       d = diagonally implicit
!       m = embedded
!       l = low-storage
!----------------------

module sll_m_rk_implicit
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_errors.h"

  use sll_m_working_precision, only: f64

  use sll_m_ode_integrator_base, only: &
    sll_c_ode, &
    sll_c_ode_integrator

  use sll_m_vector_space_base, only: &
    sll_c_vector_space

  implicit none

  public :: &
    sll_c_rk_implicit, &
    sll_t_rk1d_bwd_euler, &
    sll_t_rk1d_trapezoid

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> Working precision
  integer, parameter :: wp = f64

  !-----------------------------------------------------------------------------
  type, abstract, extends( sll_c_ode_integrator ) :: sll_c_rk_implicit

    private

    real(wp) :: rel_tol = 1.0e-12_wp
    real(wp) :: abs_tol = 1.0e-15_wp
    integer  :: maxiter = 20

  contains

    procedure :: set_params => s_rk_implicit__set_params
    procedure :: get_params => s_rk_implicit__get_params

  end type sll_c_rk_implicit

  !-----------------------------------------------------------------------------
  type, extends( sll_c_rk_implicit ) :: sll_t_rk1d_bwd_euler

  contains
    procedure :: init  => rk1d_bwd_euler__init
    procedure :: step  => rk1d_bwd_euler__step
    procedure :: clean => rk1d_bwd_euler__clean

  end type sll_t_rk1d_bwd_euler

  !-----------------------------------------------------------------------------
  type, extends( sll_c_rk_implicit ) :: sll_t_rk1d_trapezoid

  contains
    procedure :: init  => rk1d_trapezoid__init
    procedure :: step  => rk1d_trapezoid__step
    procedure :: clean => rk1d_trapezoid__clean

  end type sll_t_rk1d_trapezoid

!==============================================================================
contains
!==============================================================================

  subroutine s_rk_implicit__set_params( self, rel_tol, abs_tol, maxiter )
    class(sll_c_rk_implicit), intent(inout) :: self
    real(wp)                , intent(in   ) :: rel_tol
    real(wp)                , intent(in   ) :: abs_tol
    integer                 , intent(in   ) :: maxiter

    ! Sanity checks
    SLL_ASSERT_ALWAYS( rel_tol > 0.0_wp )
    SLL_ASSERT_ALWAYS( abs_tol > 0.0_wp )
    SLL_ASSERT_ALWAYS( maxiter > 0      )

    self % rel_tol = rel_tol
    self % abs_tol = abs_tol
    self % maxiter = maxiter

  end subroutine s_rk_implicit__set_params

  !-----------------------------------------------------------------------------
  subroutine s_rk_implicit__get_params( self, rel_tol, abs_tol, maxiter )
    class(sll_c_rk_implicit), intent(in   ) :: self
    real(wp)                , intent(  out) :: rel_tol
    real(wp)                , intent(  out) :: abs_tol
    integer                 , intent(  out) :: maxiter

    rel_tol = self % rel_tol
    abs_tol = self % abs_tol
    maxiter = self % maxiter

  end subroutine s_rk_implicit__get_params

  !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  ! Implicit Euler (backward Euler)
  !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  !-----------------------------------------------------------------------------
  subroutine rk1d_bwd_euler__init( self, ode, t0, y0 )
    class(sll_t_rk1d_bwd_euler)   , intent(  out) :: self
    class(sll_c_ode), pointer, intent(in   ) :: ode
    real(wp)                      , intent(in   ) :: t0
    class(sll_c_vector_space), intent(inout) :: y0

    self % ode => ode
    call y0 % source( self % work, 2 )
    SLL_ASSERT_ALWAYS( storage_size( y0 ) > 0 )

  end subroutine rk1d_bwd_euler__init

  !-----------------------------------------------------------------------------
  subroutine rk1d_bwd_euler__step( self, t, y, h, ynew )
    class(sll_t_rk1d_bwd_euler)   , intent(inout) :: self
    real(wp)                      , intent(in   ) :: t
    class(sll_c_vector_space), intent(in   ) :: y
    real(wp)                      , intent(in   ) :: h
    class(sll_c_vector_space), intent(inout) :: ynew

    real(wp) :: err, tol
    integer  :: iter
    logical  :: success

    character(len=*), parameter :: this_sub_name = "sll_t_rk1d_bwd_euler % step"
    character(len=256) :: err_msg

    tol = self % abs_tol + self % rel_tol * y % norm()

    success = .false.

    ! y_i = y_old
    ! work(1) = y_i
    call self % work(1) % copy( y )

    do iter = 1, self % maxiter

      ! k1 = f(t,y_i)
      ! work(2) = k1
      call self % ode % rhs( t, self % work(1), self % work(2) )

      ! y_{i+1} = y_old + h*k1
      call ynew % mult_add( h, self % work(2), y )

      ! work(1) = y_i - y_{i+1}
      call self % work(1) % incr_mult( -1.0_wp, ynew )

      err = self % work(1) % norm()

      if ( err <= tol ) then
        success = .true.
        exit
      else
        call self % work(1) % copy( ynew )
      end if

    end do

    if ( .not. success ) then
      write( err_msg, '(A,I0,A)' ) "could not converge after ", self % maxiter, " iterations"
      SLL_ERROR( this_sub_name, err_msg )
    end if

  end subroutine rk1d_bwd_euler__step

  !-----------------------------------------------------------------------------
  subroutine rk1d_bwd_euler__clean( self )
    class(sll_t_rk1d_bwd_euler), intent(inout) :: self

    deallocate( self % work )

  end subroutine rk1d_bwd_euler__clean

  !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  ! Trapezoidal
  !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  !-----------------------------------------------------------------------------
  subroutine rk1d_trapezoid__init( self, ode, t0, y0 )
    class(sll_t_rk1d_trapezoid)   , intent(  out) :: self
    class(sll_c_ode), pointer, intent(in   ) :: ode
    real(wp)                      , intent(in   ) :: t0
    class(sll_c_vector_space), intent(inout) :: y0

    self % ode => ode
    call y0 % source( self % work, 2 )
    SLL_ASSERT_ALWAYS( storage_size( y0 ) > 0 )

  end subroutine rk1d_trapezoid__init

  !-----------------------------------------------------------------------------
  subroutine rk1d_trapezoid__step( self, t, y, h, ynew )
    class(sll_t_rk1d_trapezoid)   , intent(inout) :: self
    real(wp)                      , intent(in   ) :: t
    class(sll_c_vector_space), intent(in   ) :: y
    real(wp)                      , intent(in   ) :: h
    class(sll_c_vector_space), intent(inout) :: ynew

    real(wp) :: err, tol
    integer  :: iter
    logical  :: success

    character(len=*), parameter :: this_sub_name = "sll_t_rk1d_trapezoid % step"
    character(len=256) :: err_msg

    tol = self % abs_tol + self % rel_tol * y % norm()

    success = .false.

    ! k1 = f(t,y)
    ! work(1) = k1
    call self % ode % rhs( t, y, self % work(1) )

    ! work(2) = y_i = y
    call self % work(2) % copy( y )

    do iter = 1, self % maxiter

      ! k2 = f(t+h,y_i)
      ! ynew = k2
      call self % ode % rhs( t+h, self % work(2), ynew )

      ! ynew = k1+k2
      call ynew % incr( self % work(1) )

      ! ynew = h/2*(k1+k2)
      call ynew % scal( 0.5_wp*h )

      ! ynew = y + h/2*(k1+k2)
      call ynew % incr( y )

      ! work(2) = y_i - y_{i+1}
      ! err = l2_norm( work(2) )
      call  self % work(2) % incr_mult( -1.0_wp, ynew )
      err = self % work(2) % norm()

      if ( err <= tol ) then
        success = .true.
        exit
      else
        call self % work(2) % copy( ynew )
      end if

    end do

    if ( .not. success ) then
      write( err_msg, '(A,I0,A)' ) "could not converge after ", self % maxiter, " iterations"
      SLL_ERROR( this_sub_name, err_msg )
    end if

  end subroutine rk1d_trapezoid__step

  !-----------------------------------------------------------------------------
  subroutine rk1d_trapezoid__clean( self )
    class(sll_t_rk1d_trapezoid), intent(inout) :: self

    deallocate( self % work )

  end subroutine rk1d_trapezoid__clean

end module sll_m_rk_implicit
