!> \brief 
!> The time splitting module provides a generic implementation of
!> time splitting algorithms of different order 
module sll_time_splitting
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"

  implicit none

  type, abstract :: time_splitting
     sll_real64 :: current_time = 0.0_f64
   contains 
     procedure(time_splitting_operator), pass(this), deferred :: operator1
     procedure(time_splitting_operator), pass(this), deferred :: operator2
     procedure, pass(this) :: lie_splitting
     procedure, pass(this) :: strang_splitting
  end type time_splitting

  abstract interface
     subroutine time_splitting_operator (this, dt)
       use sll_working_precision
       import time_splitting
       class(time_splitting)   :: this
       sll_real64, intent(in)  :: dt
     end subroutine time_splitting_operator
  end interface

  contains
    !> first order Lie splitting
    !> \param[inout] time_splitting object
    !> \param[in] dt : time step
    !> \param[in] number_steps : number of time steps to be performed
    subroutine lie_splitting(this, dt, number_steps)
      class(time_splitting)   :: this
      sll_real64, intent(in)  :: dt
      sll_int32, intent(in)   :: number_steps
      ! local variables
      sll_int32  :: i
      
      do i = 1, number_steps
         call this%operator1(dt)
         call this%operator2(dt)
         this%current_time = this%current_time + dt
      end do
    end subroutine lie_splitting

    !> Second order Strang splitting
    !> \param[inout] this : time_splitting object
    !> \param[in] dt : time step
    !> \param[in] number_steps : number of time steps to be performed
    subroutine strang_splitting(this, dt, number_steps)
      class(time_splitting)    :: this
      sll_real64, intent(in)   :: dt
      sll_int32, intent(in)   :: number_steps
      ! local variables
      sll_int32  :: i
     
      call this%operator1(0.5_f64*dt)
      do i = 1, number_steps - 2
         call this%operator2(dt)
         call this%operator1(dt)
         ! warning this implies that operator 1 does not depend on value of current_time
         ! if this is not true the two steps implying operator one should not be combined
         this%current_time = this%current_time + dt
      end do
      call this%operator2(dt)
      call this%operator1(0.5_f64*dt)
    end subroutine strang_splitting

end module sll_time_splitting
