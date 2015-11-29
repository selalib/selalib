module m_ode_collection

#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_errors.h"

  use sll_m_vector_space_base       , only: sll_vector_space_base
  use sll_m_vector_space_real_arrays, only: sll_vector_space_real_1d
  use sll_m_ode_integrator_base     , only: sll_ode_base

  implicit none
!  public :: harmonic_oscillator
  private

  !----------------------------------------------------------------------------
  ! Harmonic oscillator
  !----------------------------------------------------------------------------
  type, public, extends( sll_ode_base ) :: harmonic_oscillator
    sll_int32  :: ndof  = 2
    sll_real64 :: omega = 1.0_f64
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
    class( harmonic_oscillator )  , intent( inout ) :: self
    sll_real64                    , intent( in    ) :: t
    class( sll_vector_space_base ), intent( in    ) :: y
    class( sll_vector_space_base ), intent( inout ) :: ydot

    character( len=* ), parameter :: this_sub_name = "rhs__harmonic_oscillator"
    character( len=* ), parameter :: err_msg = &
      " not of type( sll_vector_space_real_1d )"

    select type( y    ); type is( sll_vector_space_real_1d )
    select type( ydot ); type is( sll_vector_space_real_1d )

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
    class( harmonic_oscillator )  , intent( in    ) :: self
    sll_real64                    , intent( in    ) :: t
    sll_real64                    , intent( in    ) :: y0(2)
    sll_real64                    , intent(   out ) :: y_ex(2)

    sll_real64 :: theta
    theta = self%omega * t

    y_ex(1) = y0(1)*cos( theta ) - y0(2)*sin( theta )
    y_ex(2) = y0(1)*sin( theta ) + y0(2)*cos( theta )
    
  end subroutine

end module m_ode_collection
