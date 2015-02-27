!> @brief Unit test for operator splitting. Linear pendulum 
!> 
!> @detail Solve linear pendulum problem: \f$ \frac{dx}{dt} = v \f$, 
!> \f$ \frac{dv}{dt} = - \omega^2 x \f$. The exact solution is
!> \[ x(t)= x(0)\cos (\omega t) + \frac{v(0)}{\omega}\sin (\omega t), ~~~ 
!> v(t)= -x(0)\omega\sin (\omega t) + v(0)\cos (\omega t) \]

program test_operator_splitting_linear_pendulum
#include "sll_working_precision.h"
#include "sll_memory.h"
  
  use sll_linear_pendulum_operators
  implicit none
 

  ! variables
  sll_real64, parameter :: x0 = 1.0_f64, v0 = 2.0_f64
  sll_real64, parameter :: t_final = 4.0_f64
  sll_real64 :: steps_fine = 100
  sll_real64 :: dt
  sll_int32  :: number_time_steps
  sll_real64 :: x_exact, v_exact 
  sll_real64 :: error0, error1, error2
  type(linear_pendulum_operators), target :: my_pendulum
  class(operator_splitting), pointer :: split
  sll_int32 :: ierr
  
  ! initialise my_pendulum
  split => my_pendulum
!  call initialize_operator_splitting(split,SLL_LIE_VT)
!  call initialize_operator_splitting(split,SLL_STRANG_TVT)
!  call initialize_operator_splitting(split,SLL_TRIPLE_JUMP_TVT)
  call initialize_operator_splitting(split,SLL_ORDER6_TVT)

  ! compute exact solution at final time
  x_exact = x0 * cos( omega * t_final ) + (v0 / omega) * sin( omega * t_final )
  v_exact = -x0 * omega * sin( omega* t_final ) + v0 * cos( omega * t_final )


  ! Test 
  
  ! do iterations with smallest time step
  my_pendulum%x = x0
  my_pendulum%v = v0  
  dt = t_final/steps_fine
  number_time_steps = steps_fine
  print*, dt, number_time_steps
  call do_split_steps(split, dt, number_time_steps)
  
  ! compute error
  error0 =sqrt( (my_pendulum%x - x_exact)**2 +  (my_pendulum%v - v_exact)**2 )

  ! do iterations with middle time step
  my_pendulum%x = x0
  my_pendulum%v = v0  
  dt = 2*dt
  number_time_steps = number_time_steps / 2
  call do_split_steps(split, dt, number_time_steps)
  
  ! compute error
  error1 =sqrt( (my_pendulum%x - x_exact)**2 +  (my_pendulum%v - v_exact)**2 )
  
  ! do iterations with largest time step
  my_pendulum%x = x0
  my_pendulum%v = v0  
  dt = 2*dt
  number_time_steps = number_time_steps / 2
  call do_split_steps(split, dt, number_time_steps)
  
  ! compute error
  error2 =sqrt( (my_pendulum%x - x_exact)**2 +  (my_pendulum%v - v_exact)**2 )

  ! compute order
  print*, 'error', error0, error1, error2
  print*, 'order',  log(error1/error0)/log(2.0_f64), log(error2/error1)/log(2.0_f64)

end program test_operator_splitting_linear_pendulum

