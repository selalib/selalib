module ode_collection

#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_errors.h"

  use sll_module_vector_space_base       , only: sll_vector_space_base
  use sll_module_vector_space_real_arrays, only: sll_vector_space_real_1d
  use sll_module_ode_integrator_base     , only: sll_ode_base

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
      ydot%array(1) = -self%omega**2 * y%array(2)
      ydot%array(2) =  self%omega**2 * y%array(1)

    class default; SLL_ERROR( this_sub_name, "ydot"//err_msg ); end select
    class default; SLL_ERROR( this_sub_name, "y"   //err_msg ); end select
  end subroutine rhs__harmonic_oscillator

end module ode_collection

!==============================================================================
program test_ode_integrators

#include "sll_working_precision.h"
#include "sll_errors.h"
  use sll_working_precision

  use sll_module_vector_space_real_arrays, only: &
    sll_vector_space_real_1d

  use sll_module_rk_explicit, only: &
    sll_rk1e_fwd_euler, &
    sll_rk2e_midpoint , &
    sll_rk2e_heun     , &
    sll_rk2e_ralston  , &
    sll_rk3e_heun3    , &
    sll_rk4e_classic

  use sll_module_ode_integrator_base, only: &
    sll_ode_base, &
    sll_ode_integrator_base

  use ode_collection, only: &
    harmonic_oscillator

  implicit none
  
  !----------------------------------------------------------------------------
  ! Variable declarations
  sll_int32                                     :: i, odeint_selector
  sll_real64                                    :: h, t, tend
  sll_real64, target                            :: z(2), znew(2)
  type( sll_vector_space_real_1d )              :: y   , ynew
  class( sll_ode_base ), allocatable, target    :: ode
  class( sll_ode_base ), pointer                :: p_ode  ! pointer to ODE
  class( sll_ode_integrator_base ), allocatable :: odeint
  
  ! Numerical constant
  sll_real64, parameter :: PI = 4.0_f64 * atan( 1.0_f64 )

  !----------------------------------------------------------------------------
  ! Choose initial conditions
  t = 0.0_f64
  z = [1.0_f64, 0.0_f64]
  
  ! Choose final time and time step size
  tend = 2.0_f64*PI
  h    = tend / 20.0_f64

  ! Choose ODE integrator
  !  .0 : forward Euler
  !  .1 : explicit midpoint
  !  .2 : explicit trapezoidal rule (Heun's method)
  !  .3 : Ralston's method
  !  .4 : 3rd-order method by Heun
  !  .5 : 'classic' RK4 by Runge
  odeint_selector = 1

  !----------------------------------------------------------------------------
  ! Create state vector (current and next solution)
  call y   %attach( z    )
  call ynew%attach( znew )

  ! Create ODE system
  allocate( harmonic_oscillator::ode )
  p_ode => ode                     ! Workaround for Intel Fortran compiler

  ! Create ODE integrator
  select case( odeint_selector )
   case( 0 ); allocate( sll_rk1e_fwd_euler::odeint )
   case( 1 ); allocate( sll_rk2e_midpoint ::odeint )
   case( 2 ); allocate( sll_rk2e_heun     ::odeint )
   case( 3 ); allocate( sll_rk2e_ralston  ::odeint )
   case( 4 ); allocate( sll_rk3e_heun3    ::odeint )
   case( 5 ); allocate( sll_rk4e_classic  ::odeint )
  end select

  ! Initialize ODE integrator
!  call odeint%init( ode, t, y )   ! This does not work with IFortran yet
  call odeint%init( p_ode, t, y )

  !----------------------------------------------------------------------------
  ! Integrate in time
  do
    write(*,*) t, y%array             ! Print solution
    if ( t .ge. tend ) exit           ! Stop if final time is reached
    call odeint%step( t, y, h, ynew ) ! Advance solution by one time-step
    call y%copy( ynew )               ! Copy new solution to current
    t = t + h                         ! Update time
  end do
  

end program test_ode_integrators
