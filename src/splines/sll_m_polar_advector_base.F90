module sll_m_polar_advector_base
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_errors.h"

  use sll_m_working_precision, only: f64

  use sll_m_constants, only: sll_p_twopi

  use sll_m_spline_2d, only: sll_t_spline_2d

  use sll_m_polar_mapping_base, only: sll_c_polar_mapping

  use sll_m_jacobian_2d_pseudo_cartesian, only: sll_t_jacobian_2d_pseudo_cartesian

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! Working precision
  integer, parameter :: wp = f64

  type, abstract :: sll_c_polar_advector

  contains

    ! Deferred procedures
    procedure(i_sub_velocity_field), deferred :: velocity_field ! RHS of characteristic equations
    procedure(i_fun_flow_field)    , deferred :: flow_field ! analytical characteristics
    procedure(i_sub_free)          , deferred :: free

    ! Non-deferred procedures
    procedure :: advect_x1x2 => f_polar_advector__advect_x1x2 ! time integrator (Cartesian coordinates)
    procedure :: advect_y1y2 => f_polar_advector__advect_y1y2 ! time integrator (intermediate coordinates)

  end type sll_c_polar_advector

  ! Interfaces for deferred procedures
  abstract interface

    SLL_PURE subroutine i_sub_velocity_field( self, x, a )
      import sll_c_polar_advector, wp
      class(sll_c_polar_advector), intent(in   ) :: self
      real(wp)                   , intent(in   ) :: x(2)
      real(wp)                   , intent(  out) :: a(2)
    end subroutine i_sub_velocity_field

    SLL_PURE function i_fun_flow_field( self, x, h ) result( y )
      import sll_c_polar_advector, wp
      class(sll_c_polar_advector), intent(in) :: self
      real(wp)                   , intent(in) :: x(2)
      real(wp)                   , intent(in) :: h
      real(wp) :: y(2)
    end function i_fun_flow_field

    subroutine i_sub_free( self )
      import sll_c_polar_advector
      class(sll_c_polar_advector), intent(inout) :: self
    end subroutine i_sub_free

  end interface

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! Advector (time integrator RK4): works in Cartesian coordinates (x1,x2)
  function f_polar_advector__advect_x1x2( self, eta, h, mapping, spline_2d_a1, spline_2d_a2, tol, maxiter ) result( eta_new )
    class(sll_c_polar_advector), intent(in) :: self
    real(wp)                   , intent(in) :: eta(2)
    real(wp)                   , intent(in) :: h
    class(sll_c_polar_mapping) , intent(in) :: mapping
    type(sll_t_spline_2d)      , intent(in) :: spline_2d_a1
    type(sll_t_spline_2d)      , intent(in) :: spline_2d_a2
    real(wp)                   , intent(in) :: tol
    integer                    , intent(in) :: maxiter
    real(wp) :: eta_new(2)

    real(wp) :: x(2), tmp(2), k1(2), k2(2), k3(2), k4(2)

    ! Initial position
    x = mapping % eval( eta )

    !---------------------------------------------------------------------------
    ! Using interpolated advection field
    !---------------------------------------------------------------------------

    ! step #1
    k1(1) = spline_2d_a1 % eval( eta(1), eta(2) )
    k1(2) = spline_2d_a2 % eval( eta(1), eta(2) )

    ! step #2
    tmp = mapping % eval_inverse( x + 0.5_wp * h * k1, eta, tol, maxiter )

    k2(1) = spline_2d_a1 % eval( tmp(1), tmp(2) )
    k2(2) = spline_2d_a2 % eval( tmp(1), tmp(2) )

    ! step #3
    tmp = mapping % eval_inverse( x + 0.5_wp * h * k2, eta, tol, maxiter )

    k3(1) = spline_2d_a1 % eval( tmp(1), tmp(2) )
    k3(2) = spline_2d_a2 % eval( tmp(1), tmp(2) )

    ! step #4
    tmp = mapping % eval_inverse( x + h * k3, eta, tol, maxiter )

    k4(1) = spline_2d_a1 % eval( tmp(1), tmp(2) )
    k4(2) = spline_2d_a2 % eval( tmp(1), tmp(2) )

!    !---------------------------------------------------------------------------
!    ! Using analytical advection field
!    !---------------------------------------------------------------------------
!
!    call self % velocity_field( x                  , k1 )
!    call self % velocity_field( x + 0.5_wp * h * k1, k2 )
!    call self % velocity_field( x + 0.5_wp * h * k2, k3 )
!    call self % velocity_field( x +          h * k3, k4 )

    ! Final position
    x = x + h * ( k1 + 2.0_wp * k2 + 2.0_wp * k3 + k4 ) / 6.0_wp

    eta_new = mapping % eval_inverse( x, eta, tol, maxiter )

  end function f_polar_advector__advect_x1x2

  ! Advector: works in intermediate coordinates (y1,y2)
  function f_polar_advector__advect_y1y2( &
    self        , &
    eta         , &
    h           , &
    jac_2d_pcart, &
    spline_2d_a1, &
    spline_2d_a2, &
    abs_tol     , &
    rel_tol     , &
    maxiter ) result( eta_new )
    class(sll_c_polar_advector)             , intent(in) :: self
    real(wp)                                , intent(in) :: eta(2)
    real(wp)                                , intent(in) :: h
    type(sll_t_jacobian_2d_pseudo_cartesian), intent(in) :: jac_2d_pcart
    type(sll_t_spline_2d)                   , intent(in) :: spline_2d_a1
    type(sll_t_spline_2d)                   , intent(in) :: spline_2d_a2
    real(wp)                                , intent(in) :: abs_tol
    real(wp)                                , intent(in) :: rel_tol
    integer                                 , intent(in) :: maxiter
    real(wp) :: eta_new(2)

    !---------------------------------------------------------------------------
    ! Trapezoidal 2nd order
    !---------------------------------------------------------------------------

!    integer  :: j
!    real(wp) :: tol_sqr
!    real(wp) :: y0(2), y(2), dy(2), temp(2), k2(2), a0(2)
!    real(wp) :: jmat_comp(2,2)
!
!    y0 = polar_map( eta )
!
!    tol_sqr = ( abs_tol + rel_tol * norm2( y0 ) )**2
!
!    jmat_comp = jac_2d_pcart % eval( eta )
!
!    ! Pseudo-Cartesian components of advection field
!    a0(1) =   jmat_comp(1,1) * spline_2d_a1 % eval( eta(1), eta(2) ) &
!            + jmat_comp(1,2) * spline_2d_a2 % eval( eta(1), eta(2) )
!    a0(2) =   jmat_comp(2,1) * spline_2d_a1 % eval( eta(1), eta(2) ) &
!            + jmat_comp(2,2) * spline_2d_a2 % eval( eta(1), eta(2) )
!
!    ! First iteration
!    dy  = h * a0
!    y   = y0 + dy
!    eta_new = polar_map_inverse( y )
!
!    ! Successive iterations if first iteration did not converge
!    if ( dot_product( dy, dy ) > tol_sqr ) then
!
!      associate( k1 => a0, h_half => 0.5_wp*h, dy_old => temp, error => temp )
!
!        do j = 2, maxiter
!
!          jmat_comp = jac_2d_pcart % eval( eta_new )
!
!          ! k2 = f(t,x_{i-1})
!          k2(1) =   jmat_comp(1,1) * spline_2d_a1 % eval( eta_new(1), eta_new(2) ) &
!                  + jmat_comp(1,2) * spline_2d_a2 % eval( eta_new(1), eta_new(2) )
!          k2(2) =   jmat_comp(2,1) * spline_2d_a1 % eval( eta_new(1), eta_new(2) ) &
!                  + jmat_comp(2,2) * spline_2d_a2 % eval( eta_new(1), eta_new(2) )
!
!          dy_old  = dy
!          dy      = h_half*(k1+k2)
!          error   = dy_old - dy
!          y       = y0 + dy
!          eta_new = polar_map_inverse( y )
!
!          if ( dot_product( error, error ) <= tol_sqr ) exit
!
!        end do
!
!      end associate
!
!    end if

    !---------------------------------------------------------------------------
    ! Runge-Kutta 3rd order
    !---------------------------------------------------------------------------

    real(wp) :: a_x1, a_x2
    real(wp) :: y(2), tmp(2), k1(2), k2(2), k3(2)
    real(wp) :: jmat_comp(2,2)

    ! Initial position
    y = polar_map( eta )

    ! step #1
    a_x1  = spline_2d_a1 % eval( eta(1), eta(2) )
    a_x2  = spline_2d_a2 % eval( eta(1), eta(2) )

    jmat_comp = jac_2d_pcart % eval( eta )

    k1(1) = jmat_comp(1,1) * a_x1 + jmat_comp(1,2) * a_x2
    k1(2) = jmat_comp(2,1) * a_x1 + jmat_comp(2,2) * a_x2

    ! step #2
    tmp = polar_map_inverse( y + 0.5_wp * h * k1 )

    a_x1  = spline_2d_a1 % eval( tmp(1), tmp(2) )
    a_x2  = spline_2d_a2 % eval( tmp(1), tmp(2) )

    jmat_comp = jac_2d_pcart % eval( tmp )

    k2(1) = jmat_comp(1,1) * a_x1 + jmat_comp(1,2) * a_x2
    k2(2) = jmat_comp(2,1) * a_x1 + jmat_comp(2,2) * a_x2

    ! step #3
    tmp = polar_map_inverse( y + h * ( 2.0_wp * k2 - k1 ) )

    a_x1  = spline_2d_a1 % eval( tmp(1), tmp(2) )
    a_x2  = spline_2d_a2 % eval( tmp(1), tmp(2) )

    jmat_comp = jac_2d_pcart % eval( tmp )

    k3(1) = jmat_comp(1,1) * a_x1 + jmat_comp(1,2) * a_x2
    k3(2) = jmat_comp(2,1) * a_x1 + jmat_comp(2,2) * a_x2

    ! Final position
    y = y + h * ( k1 + 4.0_wp * k2 + k3 ) / 6.0_wp

    eta_new = polar_map_inverse( y )

    !---------------------------------------------------------------------------
    ! Runge-Kutta 4th order
    !---------------------------------------------------------------------------

!    real(wp) :: a_x1, a_x2
!    real(wp) :: y(2), tmp(2), k1(2), k2(2), k3(2), k4(2)
!    real(wp) :: jmat_comp(2,2)
!
!    ! Initial position
!    y = polar_map( eta )
!
!    ! step #1
!    a_x1  = spline_2d_a1 % eval( eta(1), eta(2) )
!    a_x2  = spline_2d_a2 % eval( eta(1), eta(2) )
!
!    jmat_comp = jac_2d_pcart % eval( eta )
!
!    k1(1) = jmat_comp(1,1) * a_x1 + jmat_comp(1,2) * a_x2
!    k1(2) = jmat_comp(2,1) * a_x1 + jmat_comp(2,2) * a_x2
!
!    ! step #2
!    tmp = polar_map_inverse( y + 0.5_wp * h * k1 )
!
!    a_x1  = spline_2d_a1 % eval( tmp(1), tmp(2) )
!    a_x2  = spline_2d_a2 % eval( tmp(1), tmp(2) )
!
!    jmat_comp = jac_2d_pcart % eval( tmp )
!
!    k2(1) = jmat_comp(1,1) * a_x1 + jmat_comp(1,2) * a_x2
!    k2(2) = jmat_comp(2,1) * a_x1 + jmat_comp(2,2) * a_x2
!
!    ! step #3
!    tmp = polar_map_inverse( y + 0.5_wp * h * k2 )
!
!    a_x1  = spline_2d_a1 % eval( tmp(1), tmp(2) )
!    a_x2  = spline_2d_a2 % eval( tmp(1), tmp(2) )
!
!    jmat_comp = jac_2d_pcart % eval( tmp )
!
!    k3(1) = jmat_comp(1,1) * a_x1 + jmat_comp(1,2) * a_x2
!    k3(2) = jmat_comp(2,1) * a_x1 + jmat_comp(2,2) * a_x2
!
!    ! step #4
!    tmp = polar_map_inverse( y + h * k3 )
!
!    a_x1  = spline_2d_a1 % eval( tmp(1), tmp(2) )
!    a_x2  = spline_2d_a2 % eval( tmp(1), tmp(2) )
!
!    jmat_comp = jac_2d_pcart % eval( tmp )
!
!    k4(1) = jmat_comp(1,1) * a_x1 + jmat_comp(1,2) * a_x2
!    k4(2) = jmat_comp(2,1) * a_x1 + jmat_comp(2,2) * a_x2
!
!    ! Final position
!    y = y + h * ( k1 + 2.0_wp * k2 + 2.0_wp * k3 + k4 ) / 6.0_wp
!
!    eta_new = polar_map_inverse( y )

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  contains
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    SLL_PURE function polar_map( eta ) result( y )
      real(wp), intent(in) :: eta(2)
      real(wp) :: y(2)

      y(1) = eta(1) * cos( eta(2) )
      y(2) = eta(1) * sin( eta(2) )

    end function polar_map

    SLL_PURE function polar_map_inverse( y ) result( eta )
      real(wp), intent(in) :: y(2)
      real(wp) :: eta(2)

      eta(1) = norm2( y )
      eta(2) = modulo( atan2( y(2), y(1) ), sll_p_twopi )

      if ( eta(1) > 1.0_wp ) eta(1) = 1.0_wp

    end function polar_map_inverse

  end function f_polar_advector__advect_y1y2

end module sll_m_polar_advector_base
