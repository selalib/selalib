module sll_const_coef_advection_2d
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#ifdef STDF95
  use sll_cubic_spline_interpolator_1d
#else
  use sll_module_interpolators_1d_base
  use sll_time_splitting
#endif

  implicit none

#ifdef STDF95
  type :: const_coef_advection_2d
     type(cubic_spline_1d_interpolator), pointer    :: interp1, interp2
     sll_real64 :: current_time = 0.0_f64
#else
  type, extends(time_splitting) :: const_coef_advection_2d
     class(sll_interpolator_1d_base), pointer    :: interp1, interp2
#endif
     sll_real64, dimension(:,:), pointer :: data
     sll_int32 :: n1, n2
     sll_real64 :: a1, a2
#ifndef STDF95
   contains
     procedure, pass(this) :: operator1 => adv1
     procedure, pass(this) :: operator2 => adv2
#endif
  end type const_coef_advection_2d


contains
  subroutine const_coef_advection_2d_initialize(this, data, n1, n2, a1, a2, interp1, interp2)
    type(const_coef_advection_2d) :: this 
    sll_real64, dimension(:,:), target :: data
    sll_int32 :: n1, n2
    sll_real64 :: a1, a2
#ifdef STDF95
    type(cubic_spline_1d_interpolator), pointer    :: interp1, interp2
#else
    class(sll_interpolator_1d_base), pointer    :: interp1, interp2
#endif
    this%data => data
    this%n1 = n1
    this%n2 = n2
    this%a1 = a1
    this%a2 = a2
    this%interp1 => interp1
    this%interp2 => interp2
  end subroutine 

#ifdef STDF95
  subroutine const_coef_advectopn_2d_operator1(this, dt)
    type(const_coef_advection_2d) :: this 
#else
  subroutine adv1(this, dt)
    class(const_coef_advection_2d) :: this 
#endif
    sll_real64, intent(in) :: dt
    ! local variables
    sll_real64, dimension(:), pointer :: f1d
    sll_real64 :: displacement
    sll_int32 :: j
    
    do j = 1, this%n2
       displacement = this%a1 * dt
       f1d => this%data(:,j)
#ifdef STDF95
       f1d = cubic_spline_interpolate_array_at_displacement(this%interp1,this%n1, f1d, displacement)
#else
       f1d = this%interp1%interpolate_array_disp(this%n1, f1d, displacement)
#endif
    end do
  end subroutine

#ifdef STDF95
  subroutine const_coef_advection_2d_operator2(this, dt)
    type(const_coef_advection_2d) :: this 
#else
  subroutine adv2(this, dt)
    class(const_coef_advection_2d) :: this 
#endif
    sll_real64, intent(in) :: dt
    ! local variables
    sll_real64, dimension(:), pointer :: f1d
    sll_real64 :: displacement
    sll_int32 :: i
    
    do i = 1, this%n1
       displacement = this%a2 * dt
       f1d => this%data(i,:)
#ifdef STDF95
       f1d = cubic_spline_interpolate_array_at_displacement(this%interp2, this%n2, f1d, displacement)
#else
       f1d = this%interp2%interpolate_array_disp(this%n2, f1d, displacement)
#endif
    end do
  end subroutine

#ifdef STDF95
  subroutine strang_splitting(this, dt, number_steps)
    type(const_coef_advection_2d)    :: this
    sll_real64, intent(in)   :: dt
    sll_int32, intent(in)   :: number_steps
    ! local variables
    sll_int32  :: i

    call const_coef_advectopn_2d_operator1(this, 0.5_f64*dt)
    do i = 1, number_steps - 2
       call const_coef_advection_2d_operator2(this, dt)
       call const_coef_advectopn_2d_operator1(this, dt)
       ! warning this implies that operator 1 does not depend on value of current_time
       ! if this is not true the two steps implying operator one should not be combined
       this%current_time = this%current_time + dt
    end do
    call const_coef_advection_2d_operator2(this, dt)
    call const_coef_advectopn_2d_operator1(this, 0.5_f64*dt)
  end subroutine strang_splitting

  subroutine lie_splitting(this, dt, number_steps)
    type(const_coef_advection_2d)   :: this
    sll_real64, intent(in)  :: dt
    sll_int32, intent(in)   :: number_steps
    ! local variables
    sll_int32  :: i
      
    do i = 1, number_steps
       call const_coef_advectopn_2d_operator1(this, dt)
       call const_coef_advection_2d_operator2(this, dt)
       this%current_time = this%current_time + dt
    end do
  end subroutine lie_splitting
#endif

end module sll_const_coef_advection_2d
