!> @ingroup operator_splitting
!> @brief Implements split operators for constant coefficient advection
module sll_m_const_coef_advection_2d
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
  use sll_m_interpolators_1d_base
  use sll_m_operator_splitting

  implicit none

  !> @brief 
  !> Simple operator splitting type for 2D constant coefficient advection
  !> Extends operator splitting
  !> @details This should be
  !> treated as an opaque type. No access to its internals is directly allowed.
  type, extends(operator_splitting) :: const_coef_advection_2d
     !> interpolator object in first direction
     class(sll_c_interpolator_1d), pointer    :: interp1
     !> interpolator object in second direction
     class(sll_c_interpolator_1d), pointer    :: interp2 
     !> function do be evolved
     sll_real64, dimension(:,:), pointer :: data
     !> dimension in first direction
     sll_int32 :: n1
      !> dimension in second direction
     sll_int32 :: n2
     !> advection coefficient first direction
     sll_real64 :: a1
     !> advection coefficient second direction
     sll_real64 :: a2
   contains
     procedure, pass(this) :: operatorT => adv1  !< advection in first direction
     procedure, pass(this) :: operatorV => adv2  !< advection in second direction
  end type const_coef_advection_2d


contains
  function new_const_coef_advection_2d( data, n1, n2, a1, a2, interp1, interp2, &
       split_case, split_step, nb_split_step, split_begin_T, dt) &
       result(this)  
    class(const_coef_advection_2d), pointer :: this  !< object to be initialised
    sll_real64, dimension(:,:), pointer, intent(in) :: data   !< initial value of function
    sll_int32, intent(in)  :: n1   !< dimension in first direction
    sll_int32, intent(in)  :: n2   !< dimension in second direction
    sll_real64, intent(in) :: a1   !< advection coeeficient in first direction
    sll_real64, intent(in) :: a2   !< advection coeeficient in first direction
    class(sll_c_interpolator_1d), pointer    :: interp1  !< interpolator for first direction
    class(sll_c_interpolator_1d), pointer    :: interp2  !< interpolator for second direction
    sll_int32, intent(in)  :: split_case  !< defines  splitting method
    sll_real64, dimension(:), intent(in), optional :: split_step  !< coefficients of split step
    sll_int32, intent(in), optional :: nb_split_step !< number of split steps
    logical, intent(in), optional :: split_begin_T   !< begin with operator T if .true.
    sll_real64, intent(in), optional :: dt  !< time step   
    ! local variable
    sll_int32 :: ierr

    SLL_ALLOCATE(this,ierr)   
    call initialize_const_coef_advection_2d( this, data, n1, n2, a1, a2, interp1, interp2, &
         split_case, split_step, nb_split_step, split_begin_T, dt )

  end function

  !> @brief Initialise const_coef_advection_2d object
  subroutine initialize_const_coef_advection_2d( this, data, n1, n2, a1, a2, interp1, interp2, &
       split_case, split_step, nb_split_step, split_begin_T, dt)
    class(const_coef_advection_2d), intent(inout)   :: this !< object 
    sll_real64, dimension(:,:), pointer, intent(in) :: data   !< initial value of function
    sll_int32, intent(in)  :: n1   !< dimension in first direction
    sll_int32, intent(in)  :: n2   !< dimension in second direction
    sll_real64, intent(in) :: a1   !< advection coefficient in first direction
    sll_real64, intent(in) :: a2   !< advection coefficient in second direction
    class(sll_c_interpolator_1d), pointer    :: interp1  !< interpolator for first direction
    class(sll_c_interpolator_1d), pointer    :: interp2  !< interpolator for second direction
    sll_int32, intent(in)  :: split_case  !< defines  splitting method
    sll_real64, dimension(:), intent(in), optional :: split_step  !< coefficients of split step
    sll_int32, intent(in), optional :: nb_split_step !< number of split steps
    logical, intent(in), optional :: split_begin_T   !< begin with operator T if .true.
    sll_real64, intent(in), optional :: dt  !< time step

    this%data => data
    this%n1 = n1
    this%n2 = n2
    this%a1 = a1
    this%a2 = a2
    this%interp1 => interp1
    this%interp2 => interp2

    call initialize_operator_splitting( &
         this, &
         split_case, &
         split_step, &
         nb_split_step, &
         split_begin_T, &
         dt)
  end subroutine 

  !> @brief Constant coefficient advection operator in first direction
  subroutine adv1(this, dt)
    class(const_coef_advection_2d), intent(inout) :: this !< object 
    sll_real64, intent(in)                        :: dt   !< time step
    ! local variables
    sll_real64, dimension(:), pointer :: f1d
    sll_real64 :: displacement
    sll_int32 :: j
    
    do j = 1, this%n2
       displacement = -this%a1 * dt
       f1d => this%data(:,j)
       call this%interp1%interpolate_array_disp(this%n1, f1d, displacement, f1d)
    end do
  end subroutine

  !> @brief Constant coefficient advection operator in second direction
  subroutine adv2(this, dt)
    class(const_coef_advection_2d), intent(inout) :: this !< object 
    sll_real64, intent(in)                        :: dt   !< time step
    ! local variables
    sll_real64, dimension(:), pointer :: f1d
    sll_real64 :: displacement
    sll_int32 :: i
    
    do i = 1, this%n1
       displacement = -this%a2 * dt
       f1d => this%data(i,:)
       call this%interp2%interpolate_array_disp(this%n2, f1d, displacement, f1d)
    end do
  end subroutine


end module sll_m_const_coef_advection_2d
