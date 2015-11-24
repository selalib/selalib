program test_ode_integrators

#include "sll_working_precision.h"
#include "sll_errors.h"
  use sll_m_working_precision
  use sll_m_errors

  use sll_m_vector_space_real_arrays, only: &
    sll_vector_space_real_1d

  use sll_m_rk_explicit, only: &
    sll_rk1e_fwd_euler, &
    sll_rk2e_midpoint , &
    sll_rk2e_heun     , &
    sll_rk2e_ralston  , &
    sll_rk3e_heun3    , &
    sll_rk4e_classic

  use sll_m_ode_integrator_base, only: &
    sll_ode_base, &
    sll_ode_integrator_base

  use m_ode_collection, only: &
    harmonic_oscillator

  implicit none
  
  !----------------------------------------------------------------------------
  ! Input arguments
  sll_int32  :: ode_selector
  sll_int32  :: odeint_selector
  sll_real64 :: h_factor

  ! ODE test-case {ode_selector}
  !  .0 : harmonic oscillator

  ! ODE integrator {odeint_selector}
  !  .0 : forward Euler
  !  .1 : explicit midpoint
  !  .2 : explicit trapezoidal rule (Heun's method)
  !  .3 : Ralston's method
  !  .4 : 3rd-order method by Heun
  !  .5 : 'classic' RK4 by Runge
  
  ! Time step selection:
  ! h = h_max / h_factor

  !----------------------------------------------------------------------------
  ! Numerical constants
  sll_real64, parameter :: PI = 4.0_f64 * atan( 1.0_f64 )

  ! Error handling
  character(len=*), parameter :: this_prog_name = "test_ode_integrators"
  character(len=128)          :: err_msg

  ! Other variables
  sll_real64                                    :: h_max, tend
  sll_real64                                    :: h, t
  sll_int32                                     :: i
  sll_real64, target                            :: z(2), znew(2)
  sll_real64                                    :: z0(2), z_ex(2)
  sll_real64                                    :: max_err
  type( sll_vector_space_real_1d )              :: y   , ynew
  class( sll_ode_base ), allocatable, target    :: ode
  class( sll_ode_base ), pointer                :: p_ode  ! pointer to ODE
  class( sll_ode_integrator_base ), allocatable :: odeint

  !----------------------------------------------------------------------------
  ! Get input arguments (from command line)
  call parse_input_arguments( ode_selector, odeint_selector, h_factor )

  ! Initial conditions (from test-case)
  t  = 0.0_f64
  z0 = [1.0_f64, 0.0_f64]
  
  ! Final time and maximum allowed time-step (from test-case)
  tend  = 2.0_f64*PI
  h_max = PI / 2.0_f64

  ! Compute actual time-step size
  h = h_max / h_factor

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

  ! Give initial conditions
  y%array = z0

  ! Maximum error during simulation
  max_err = 0.0_f64

  !----------------------------------------------------------------------------
  ! Integrate in time
  do
    write(*,*) t, y%array             ! Print solution

    select type( ode )
      class is( harmonic_oscillator )
        call ode%y_ex( t, z0, z_ex )      ! Write exact solution to z_ex
    end select

    max_err = max( max_err, norm2( y%array-z_ex ) )
    
    if ( t .ge. tend ) exit           ! Stop if final time is reached
    call odeint%step( t, y, h, ynew ) ! Advance solution by one time-step
    call y%copy( ynew )               ! Copy new solution to current
    t = t + h                         ! Update time
  end do
  
  !----------------------------------------------------------------------------
  ! Write maximum error
  write(*,*)
  write(*,*) 'Maximum error in simulation: ', max_err           ! Print error
  write(*,*)

!==============================================================================
contains
!==============================================================================

  subroutine parse_input_arguments( ode_selector, odeint_selector, h_factor )
    sll_int32 , intent( out ) :: ode_selector
    sll_int32 , intent( out ) :: odeint_selector
    sll_real64, intent( out ) :: h_factor

    character(len= 32) :: arg
    character(len=128) :: err_msg

    ! Check if 3 input arguments are found
    if (command_argument_count() /= 3) then
      err_msg = "Exactly 3 input arguments are required"
      SLL_ERROR( this_prog_name, trim( err_msg ) )
    end if

    ! Extract input arguments
    call get_command_argument( 1, arg );  read( arg, * ) ode_selector
    call get_command_argument( 2, arg );  read( arg, * ) odeint_selector
    call get_command_argument( 3, arg );  read( arg, * ) h_factor

    ! Verify acceptable values
    select case (ode_selector)
     case (0)
      ! OK
     case default
      err_msg = "ode_selector can only be = 1 for now"
      SLL_ERROR( this_prog_name, trim( err_msg ) )
    end select

    if (odeint_selector < 0 .or. 5 < odeint_selector) then
      err_msg = "odeint_selector can only assume integer values between 0 and 5"
      SLL_ERROR( this_prog_name, trim( err_msg ) )
    end if

    if (h_factor < 1.0_f64) then
      err_msg = "h_factor must be a positive number larger than 1" 
      SLL_ERROR( this_prog_name, trim( err_msg ) )
    end if

  end subroutine

end program test_ode_integrators
