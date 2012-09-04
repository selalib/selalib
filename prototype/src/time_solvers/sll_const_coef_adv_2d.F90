module sll_const_coef_advection_2d
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
  use sll_time_splitting
  use sll_module_interpolators_1d_base

  implicit none

  type, extends(time_splitting) :: const_coef_advection_2d
     sll_real64, dimension(:,:), pointer :: data
     sll_int32 :: n1, n2
     sll_real64 :: a1, a2
     class(sll_interpolator_1d_base), pointer    :: interp1, interp2
   contains
     procedure, pass(this) :: operator1 => adv1
     procedure, pass(this) :: operator2 => adv2
  end type const_coef_advection_2d

contains
  subroutine initialize(this, data, n1, n2, a1, a2, interp1, interp2)
    type(const_coef_advection_2d) :: this 
    sll_real64, dimension(:,:), target :: data
    sll_int32 :: n1, n2
    sll_real64 :: a1, a2
    class(sll_interpolator_1d_base), pointer    :: interp1, interp2
    this%data => data
    this%n1 = n1
    this%n2 = n2
    this%a1 = a1
    this%a2 = a2
    this%interp1 => interp1
    this%interp2 => interp2
  end subroutine initialize

  subroutine adv1(this, dt)
    class(const_coef_advection_2d) :: this 
    sll_real64 :: dt
    ! local variables
    sll_real64, dimension(:), pointer :: f1d
    sll_real64 :: displacement
    sll_int32 :: j
    
    do j = 1, this%n2
       displacement = this%a1 * dt
       f1d => this%data(:,j)
       f1d = this%interp1%interpolate_array_disp(this%n1, f1d, displacement)
    end do
  end subroutine adv1

  subroutine adv2(this, dt)
    class(const_coef_advection_2d) :: this 
    sll_real64 :: dt
    ! local variables
    sll_real64, dimension(:), pointer :: f1d
    sll_real64 :: displacement
    sll_int32 :: i
    
    do i = 1, this%n1
       displacement = this%a2 * dt
       f1d => this%data(i,:)
       f1d = this%interp2%interpolate_array_disp(this%n2, f1d, displacement)
    end do
  end subroutine adv2
end module sll_const_coef_advection_2d
