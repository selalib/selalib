module m_ode_collection
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_errors.h"

  use sll_m_working_precision, only: &
    f64

  use sll_m_ode_integrator_base, only: &
    sll_c_ode

  use sll_m_vector_space_base, only: &
    sll_c_vector_space

  use sll_m_vector_space_real_arrays, only: &
    sll_t_vector_space_real_1d

  implicit none

  public :: &
    harmonic_oscillator

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! Working precision
  integer, parameter :: wp = f64

  !----------------------------------------------------------------------------
  ! Harmonic oscillator
  !----------------------------------------------------------------------------
  type, extends( sll_c_ode ) :: harmonic_oscillator
    integer  :: ndof  = 2
    real(wp) :: omega = 1.0_f64
  contains
    procedure :: rhs => rhs__harmonic_oscillator
    procedure :: y_ex => y_ex__harmonic_oscillator
  end type harmonic_oscillator

!==============================================================================
contains
!==============================================================================

  !----------------------------------------------------------------------------
  ! Harmonic oscillator
  !----------------------------------------------------------------------------
  subroutine rhs__harmonic_oscillator( self, t, y, ydot )
    class( harmonic_oscillator ), intent( inout ) :: self
    real(wp)                    , intent( in    ) :: t
    class( sll_c_vector_space ) , intent( in    ) :: y
    class( sll_c_vector_space ) , intent( inout ) :: ydot

    character( len=* ), parameter :: this_sub_name = "rhs__harmonic_oscillator"
    character( len=* ), parameter :: err_msg = &
      " not of type( sll_t_vector_space_real_1d )"

    select type( y    ); type is( sll_t_vector_space_real_1d )
    select type( ydot ); type is( sll_t_vector_space_real_1d )

      SLL_ASSERT( size( y   %array ) .eq. 2 )
      SLL_ASSERT( size( ydot%array ) .eq. 2 )

      ! Actual R.H.S. computation
      ydot%array(1) = -self%omega * y%array(2)
      ydot%array(2) =  self%omega * y%array(1)

    class default; SLL_ERROR( this_sub_name, "ydot"//err_msg ); end select
    class default; SLL_ERROR( this_sub_name, "y"   //err_msg ); end select
  end subroutine rhs__harmonic_oscillator

  ! Exact solution
  subroutine y_ex__harmonic_oscillator( self, t, y0, y_ex )
    class( harmonic_oscillator ), intent( in    ) :: self
    real(wp)                    , intent( in    ) :: t
    real(wp)                    , intent( in    ) :: y0(2)
    real(wp)                    , intent(   out ) :: y_ex(2)

    real(wp) :: theta
    theta = self%omega * t

    y_ex(1) = y0(1)*cos( theta ) - y0(2)*sin( theta )
    y_ex(2) = y0(1)*sin( theta ) + y0(2)*cos( theta )
    
  end subroutine

end module m_ode_collection
