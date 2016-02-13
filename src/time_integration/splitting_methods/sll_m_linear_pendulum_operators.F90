!> @ingroup sll_t_operator_splitting
!> @brief Implements split operators for linear pendulum 
!!
!> @details Solve linear pendulum problem: \f$ \frac{dx}{dt} = v \f$, 
!> \f$ \frac{dv}{dt} = - \omega^2 x \f$. The exact solution is
!> \f[ x(t)= x(0)\cos (\omega t) + \frac{v(0)}{\omega}\sin (\omega t), ~~~ 
!> v(t)= -x(0)\omega\sin (\omega t) + v(0)\cos (\omega t) \f]

module sll_m_linear_pendulum_operators
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_operator_splitting, only: &
    sll_s_do_split_steps, &
    sll_s_initialize_operator_splitting, &
    sll_t_operator_splitting

  implicit none

  public :: &
    sll_s_check_order

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> @brief 
  !> Simple operator splitting type for linear pendulum
  !> Extends operator splitting
  !> @details This should be
  !> treated as an opaque type. No access to its internals is directly allowed.
  type, extends(sll_t_operator_splitting) :: linear_pendulum_operators
     sll_real64 :: x  !< x value
     sll_real64 :: v  !< v value
   contains
     procedure, pass(this) :: operatorT => push_x !< definition of first operator of splitting
     procedure, pass(this) :: operatorV => push_v !< definition of second operator of splitting
  end type linear_pendulum_operators

  ! module global variables
  sll_real64, parameter :: omega = 2.0_f64       !< frequency
  sll_real64, parameter :: x0 = 1.0_f64          !< initial x for order checking
  sll_real64, parameter :: v0 = 2.0_f64          !< initial v for order checking 
  sll_real64, parameter :: t_final = 1.0_f64     !< final time for order checking 
contains
  
  !> Implements the first operator of splitting for linear pendulum
  subroutine push_x (this, dt)
    class(linear_pendulum_operators), intent(inout)   :: this !< object
    sll_real64, intent(in)                            :: dt   !< time step

    this%x = this%x + this%v * dt
  end subroutine push_x

  !> Implements the second operator of splitting for linear pendulum
  subroutine push_v (this, dt)
    class(linear_pendulum_operators), intent(inout)   :: this !< object
    sll_real64, intent(in)                            :: dt   !< time step            
    
    this%v = this%v -omega**2 * this%x * dt
  end subroutine push_v

  !> checks the order of a splitting method on the linear pendulum.
  !> used for unit testing.
  subroutine sll_s_check_order( method, steps_fine, expected_order, test_passed )
    sll_int32, intent(in)  :: method         !< splitting method to be chosen from those 
                                             !< implemented in sll_m_operator_splitting
    sll_real64, intent(in) :: steps_fine     !< number of steps on fine grid
    sll_int32, intent(in)  :: expected_order !< expected_order of the method 
    logical, intent(inout) :: test_passed    !< check if test successful 

    ! local variables
    sll_real64 :: dt
    sll_int32  :: number_time_steps
    sll_real64 :: x_exact, v_exact 
    sll_real64 :: error0, error1, error2, order1, order2
    type(linear_pendulum_operators), target :: my_pendulum
    class(sll_t_operator_splitting), pointer :: split
    !sll_int32 :: ierr
    
    ! initialise my_pendulum
    split => my_pendulum
    !  call sll_s_initialize_operator_splitting(split,sll_p_lie_vt)
    !  call sll_s_initialize_operator_splitting(split,sll_p_strang_tvt)
    !  call sll_s_initialize_operator_splitting(split,sll_p_triple_jump_tvt)
    call sll_s_initialize_operator_splitting(split,method)

    ! compute exact solution at final time
    x_exact = x0 * cos( omega * t_final ) + (v0 / omega) * sin( omega * t_final )
    v_exact = -x0 * omega * sin( omega* t_final ) + v0 * cos( omega * t_final )
  
    ! do iterations with smallest time step
    my_pendulum%x = x0
    my_pendulum%v = v0  
    dt = t_final/steps_fine
    number_time_steps = int(steps_fine,i32)

    call sll_s_do_split_steps(split, dt, number_time_steps)
  
    ! compute  mean square error
    error0 =sqrt( (my_pendulum%x - x_exact)**2 +  (my_pendulum%v - v_exact)**2 )

    ! do iterations with middle time step
    my_pendulum%x = x0
    my_pendulum%v = v0  
    dt = 2*dt
    number_time_steps = number_time_steps / 2
    call sll_s_do_split_steps(split, dt, number_time_steps)
  
    ! compute mean square error
    error1 = sqrt( (my_pendulum%x - x_exact)**2 +  (my_pendulum%v - v_exact)**2 )
  
    ! do iterations with largest time step
    my_pendulum%x = x0
    my_pendulum%v = v0  
    dt = 2*dt
    number_time_steps = number_time_steps / 2
    call sll_s_do_split_steps(split, dt, number_time_steps)
  
    ! compute  mean square error
    error2 =sqrt( (my_pendulum%x - x_exact)**2 +  (my_pendulum%v - v_exact)**2 )

    ! compute order
    !print*, 'error', error0, error1, error2
    order1 = log(error1/error0)/log(2.0_f64)
    order2 = log(error2/error1)/log(2.0_f64)
    if ((abs(order1-expected_order) > 5.e-2) .or. (abs(order2-expected_order) > 5.e-2)) then       
       test_passed = .false.
       print*, 'error coarse =', error2
       print*, 'error middle =', error1
       print*, '      order (coarse/middle) =', order2
       print*, 'error fine   =', error0
       print*, '      order (middle/fine) =', order1
    endif
  end subroutine sll_s_check_order

end module sll_m_linear_pendulum_operators
