!> @ingroup operator_splitting
!> @brief Implements split operators for constant coefficient advection
module sll_advection_2d
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
  use sll_module_advection_1d_base
  use sll_operator_splitting

  implicit none

  !> @brief 
  !> Simple operator splitting type for 2D  advection
  !> Extends operator splitting
  !> @details This should be
  !> treated as an opaque type. No access to its internals is directly allowed.
  !> comment: how to change then the advection fields in time?
  type, extends(operator_splitting) :: advection_2d
     !> advector object in first direction
     class(sll_advection_1d_base), pointer    :: advector1
     !> advector object in second direction
     class(sll_advection_1d_base), pointer    :: advector2 
     !> function do be evolved
     sll_real64, dimension(:,:), pointer :: data
     !> dimension in first direction
     sll_int32 :: n1
      !> dimension in second direction
     sll_int32 :: n2
     !> advection coefficient first direction
     sll_real64, dimension(:,:), pointer  :: a1
     !> advection coefficient second direction
     sll_real64, dimension(:,:), pointer  :: a2
     !> temporary array
     sll_real64, dimension(:), pointer :: buf1d
   contains
     procedure, pass(this) :: operatorT => adv1  !< advection in first direction
     procedure, pass(this) :: operatorV => adv2  !< advection in second direction
  end type advection_2d


contains
  function new_advection_2d( &
      data, &
      n1, &
      n2, &
      a1, &
      a2, &
      advector1, &
      advector2, &
      split_case, &
      split_step, &
      nb_split_step, &
      split_begin_T, &
      dt) &
      result(this)  
    class(advection_2d), pointer :: this  !< object to be initialised
    sll_real64, dimension(:,:), pointer, intent(in) :: data   !< initial value of function
    sll_int32, intent(in)  :: n1   !< dimension in first direction
    sll_int32, intent(in)  :: n2   !< dimension in second direction
    sll_real64, dimension(:,:), pointer, intent(in) :: a1   !< advection coefficient in first direction
    sll_real64, dimension(:,:), pointer, intent(in) :: a2   !< advection coefficient in second direction
    class(sll_advection_1d_base), pointer    :: advector1  !< advector for first direction
    class(sll_advection_1d_base), pointer    :: advector2  !< advector for second direction
    sll_int32, intent(in)  :: split_case  !< defines  splitting method
    sll_real64, dimension(:), intent(in), optional :: split_step  !< coefficients of split step
    sll_int32, intent(in), optional :: nb_split_step !< number of split steps
    logical, intent(in), optional :: split_begin_T   !< begin with operator T if .true.
    sll_real64, intent(in), optional :: dt  !< time step   
    ! local variable
    sll_int32 :: ierr

    SLL_ALLOCATE(this,ierr)   
    call initialize_advection_2d( &
      this, &
      data, &
      n1, &
      n2, &
      a1, &
      a2, &
      advector1, &
      advector2, &
      split_case, &
      split_step, &
      nb_split_step, &
      split_begin_T, &
      dt )

  end function

  !> @brief Initialise advection_2d object
  subroutine initialize_advection_2d( &
      this, &
      data, &
      n1, &
      n2, &
      a1, &
      a2, &
      advector1, &
      advector2, &
      split_case, &
      split_step, &
      nb_split_step, &
      split_begin_T, &
      dt)
    class(advection_2d), intent(inout)   :: this !< object 
    sll_real64, dimension(:,:), pointer, intent(in) :: data   !< initial value of function
    sll_int32, intent(in)  :: n1   !< dimension in first direction
    sll_int32, intent(in)  :: n2   !< dimension in second direction
    sll_real64, dimension(:,:), pointer,  intent(in) :: a1   !< advection coefficient in first direction
    sll_real64, dimension(:,:), pointer,  intent(in) :: a2   !< advection coefficient in second direction
    class(sll_advection_1d_base), pointer    :: advector1  !< advector for first direction
    class(sll_advection_1d_base), pointer    :: advector2  !< advector for second direction
    sll_int32, intent(in)  :: split_case  !< defines  splitting method
    sll_real64, dimension(:), intent(in), optional :: split_step  !< coefficients of split step
    sll_int32, intent(in), optional :: nb_split_step !< number of split steps
    logical, intent(in), optional :: split_begin_T   !< begin with operator T if .true.
    sll_real64, intent(in), optional :: dt  !< time step
    sll_int32 :: ierr

    this%data => data
    this%n1 = n1
    this%n2 = n2
    this%a1 => a1
    this%a2 => a2
    this%advector1 => advector1
    this%advector2 => advector2

    SLL_ALLOCATE(this%buf1d(max(n1,n2)),ierr)

    call initialize_operator_splitting( &
      this, &
      split_case, &
      split_step, &
      nb_split_step, &
      split_begin_T, &
      dt)
  end subroutine 

  !> @brief Advection operator in first direction
  subroutine adv1(this, dt)
    class(advection_2d), intent(inout) :: this !< object 
    sll_real64, intent(in) :: dt   !< time step
    ! local variables
    sll_int32 :: j
    sll_int32 :: n1
    sll_int32 :: n2

    n1 =this%n1
    n2 =this%n2
        
    do j=1,n2
      this%buf1d(1:n1) = this%data(1:n1,j)
      call this%advector1%advect_1d( &
        this%A1(1:n1,j), &
        dt, &
        this%buf1d(1:n1), &
        this%data(1:n1,j))
    enddo
  end subroutine


  !> @brief Advection operator in second direction
  subroutine adv2(this, dt)
    class(advection_2d), intent(inout) :: this !< object 
    sll_real64, intent(in) :: dt   !< time step
    ! local variables
    sll_int32 :: i
    sll_int32 :: n1
    sll_int32 :: n2

    n1 =this%n1
    n2 =this%n2
        
    do i=1,n1
      this%buf1d(1:n2) = this%data(i,1:n2)
      call this%advector1%advect_1d( &
        this%A2(i,1:n2), &
        dt, &
        this%buf1d(1:n2), &
        this%data(i,1:n2))
    enddo
  end subroutine


end module sll_advection_2d
