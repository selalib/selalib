!> @ingroup operator_splitting
!> @brief Implements split operators T and V for constant coefficient advection
module sll_const_coef_advection_2d
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
  use sll_module_interpolators_1d_base
  use sll_operator_splitting

  implicit none

  !> @brief 
  !> Simple operator splitting type for 2D constant coefficient advection
  !> Extends operator splitting
  !> @details This should be
  !> treated as an opaque type. No access to its internals is directly allowed.
  type, extends(operator_splitting) :: const_coef_advection_2d
     class(sll_interpolator_1d_base), pointer :: interp1 !< x1 1d interpolator
     class(sll_interpolator_1d_base), pointer :: interp2 !< x2 1d interpolator
     sll_real64, dimension(:,:), pointer      :: data    !< 2d array
     sll_int32                                :: n1      !< first dimesnion size
     sll_int32                                :: n2      !< second dimesnion size
     sll_real64                               :: a1      !< advection first coefficient 
     sll_real64                               :: a2      !< advection second coefficient
   contains
     !> Operator for T advection
     procedure, pass(this) :: operatorT => adv1
     !> Operator for V advection
     procedure, pass(this) :: operatorV => adv2
  end type const_coef_advection_2d


contains

  !> Initialize operator splitting type
  subroutine const_coef_advection_2d_initialize(this, data, n1, n2, a1, a2, interp1, interp2)
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
  end subroutine 

  !> Advection along direction x1
  subroutine adv1(this, dt)
    class(const_coef_advection_2d), intent(inout) :: this 
    sll_real64, intent(in) :: dt
    ! local variables
    sll_real64, dimension(:), pointer :: f1d
    sll_real64 :: displacement
    sll_int32 :: j
    
    do j = 1, this%n2
       displacement = this%a1 * dt
       f1d => this%data(:,j)
       f1d = this%interp1%interpolate_array_disp(this%n1, f1d, displacement)
    end do
  end subroutine

  !> Advection along direction x2
  subroutine adv2(this, dt)
    class(const_coef_advection_2d), intent(inout) :: this 
    sll_real64, intent(in) :: dt
    ! local variables
    sll_real64, dimension(:), pointer :: f1d
    sll_real64 :: displacement
    sll_int32 :: i
    
    do i = 1, this%n1
       displacement = this%a2 * dt
       f1d => this%data(i,:)
       f1d = this%interp2%interpolate_array_disp(this%n2, f1d, displacement)
    end do
  end subroutine


end module sll_const_coef_advection_2d
