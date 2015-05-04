!> @ingroup operator_splitting
!> @brief Implements split operators for constant coefficient advection
module sll_advection_2d
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
  use sll_module_advection_1d_base
  use sll_operator_splitting

  implicit none

  sll_int32, parameter :: SLL_ADVECTIVE    = 0
  sll_int32, parameter :: SLL_CONSERVATIVE = 1

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
     sll_real64, dimension(:,:), pointer :: fdata
     !> dimension in first direction
     sll_int32 :: n1
      !> dimension in second direction
     sll_int32 :: n2
     !> number of DOF in first direction
     sll_int32 :: num_dof1
     !> number of DOF in second direction
     sll_int32 :: num_dof2
     !> advection coefficient first direction
     sll_real64, dimension(:,:), pointer  :: a1
     !> advection coefficient second direction
     sll_real64, dimension(:,:), pointer  :: a2
     !> temporary array
     sll_real64, dimension(:,:), pointer :: buf1d
     !> advection form (SLL_ADVECTIVE or SLL_CONSERVATIVE)
     sll_int32 :: advection_form
   contains
     procedure, pass(this) :: operatorT => adv1  !< advection in first direction
     procedure, pass(this) :: operatorV => adv2  !< advection in second direction
  end type advection_2d


contains
  function new_advection_2d( &
      fdata, &
      n1, &
      n2, &
      num_dof1, &
      num_dof2, &
      a1, &
      a2, &
      advector1, &
      advector2, &
      advection_form, &
      split_case, &
      split_step, &
      nb_split_step, &
      split_begin_T, &
      dt) &
      result(this)  
    class(advection_2d), pointer :: this  !< object to be initialised
    sll_real64, dimension(:,:), pointer, intent(in) :: fdata   !< initial value of function
    sll_int32, intent(in)  :: n1   !< dimension in first direction
    sll_int32, intent(in)  :: n2   !< dimension in second direction
    sll_int32, intent(in)  :: num_dof1   !< number of DOF in first direction
    sll_int32, intent(in)  :: num_dof2   !< number of DOF in second direction
    sll_real64, dimension(:,:), pointer, intent(in) :: a1   !< advection coefficient in first direction
    sll_real64, dimension(:,:), pointer, intent(in) :: a2   !< advection coefficient in second direction
    class(sll_advection_1d_base), pointer    :: advector1  !< advector for first direction
    class(sll_advection_1d_base), pointer    :: advector2  !< advector for second direction
    sll_int32, intent(in) :: advection_form
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
      fdata, &
      n1, &
      n2, &
      num_dof1, &
      num_dof2, &
      a1, &
      a2, &
      advector1, &
      advector2, &
      advection_form, &
      split_case, &
      split_step, &
      nb_split_step, &
      split_begin_T, &
      dt )

  end function

  !> @brief Initialise advection_2d object
  subroutine initialize_advection_2d( &
      this, &
      fdata, &
      n1, &
      n2, &
      num_dof1, &
      num_dof2, &
      a1, &
      a2, &
      advector1, &
      advector2, &
      advection_form, &
      split_case, &
      split_step, &
      nb_split_step, &
      split_begin_T, &
      dt)
    class(advection_2d), intent(inout)   :: this !< object 
    sll_real64, dimension(:,:), pointer, intent(in) :: fdata   !< initial value of function
    sll_int32, intent(in)  :: n1   !< dimension in first direction
    sll_int32, intent(in)  :: n2   !< dimension in second direction
    sll_int32, intent(in)  :: num_dof1   !< number of DOF in first direction
    sll_int32, intent(in)  :: num_dof2   !< number of DOF in second direction
    sll_real64, dimension(:,:), pointer,  intent(in) :: a1   !< advection coefficient in first direction
    sll_real64, dimension(:,:), pointer,  intent(in) :: a2   !< advection coefficient in second direction
    class(sll_advection_1d_base), pointer    :: advector1  !< advector for first direction
    class(sll_advection_1d_base), pointer    :: advector2  !< advector for second direction
    sll_int32, intent(in) :: advection_form
    sll_int32, intent(in)  :: split_case  !< defines  splitting method
    sll_real64, dimension(:), intent(in), optional :: split_step  !< coefficients of split step
    sll_int32, intent(in), optional :: nb_split_step !< number of split steps
    logical, intent(in), optional :: split_begin_T   !< begin with operator T if .true.
    sll_real64, intent(in), optional :: dt  !< time step
    sll_int32 :: ierr

    this%fdata => fdata
    this%n1 = n1
    this%n2 = n2
    this%num_dof1 = num_dof1
    this%num_dof2 = num_dof2
    this%a1 => a1
    this%a2 => a2
    this%advector1 => advector1
    this%advector2 => advector2
    
    this%advection_form = advection_form
    SLL_ALLOCATE(this%buf1d(max(n1,n2),2),ierr)

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
    sll_int32 :: num_dof1
    sll_int32 :: num_dof2

    n1 =this%n1
    n2 =this%n2
    num_dof1 =this%num_dof1
    num_dof2 =this%num_dof2
        
    do j=1,num_dof2
      this%buf1d(1:num_dof1,1) = this%fdata(1:num_dof1,j)
      
      call this%advector1%advect_1d( &
        this%A1(1:n1,j), &
        dt, &
        this%buf1d(1:n1,1), &
        this%buf1d(1:n1,2))
      
      this%fdata(1:num_dof1,j) = this%buf1d(1:num_dof1,2)  
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
    sll_int32 :: num_dof1
    sll_int32 :: num_dof2

    n1 =this%n1
    n2 =this%n2
    num_dof1 =this%num_dof1
    num_dof2 =this%num_dof2
        
    do i=1,num_dof1
      this%buf1d(1:num_dof2,1) = this%fdata(i,1:num_dof2)
      
      call this%advector2%advect_1d( &
        this%A2(i,1:n2), &
        dt, &
        this%buf1d(1:n2,1), &
        this%buf1d(1:n2,2))

      this%fdata(i,1:num_dof2) = this%buf1d(1:num_dof2,2)
        
    enddo
  end subroutine


  subroutine function_to_primitive_adv(f,node_positions,N,M)
    sll_real64,dimension(:),intent(inout) :: f
    sll_real64,dimension(:),intent(in) :: node_positions
    sll_int32,intent(in):: N
    sll_real64,intent(out)::M
    sll_int32::i
    sll_real64::tmp,tmp2
        
    !from f compute the mean
    M=0._f64
    do i=1,N
      M=M+f(i)*(node_positions(i+1)-node_positions(i))
    enddo
    
    f(1)=(f(1)-M)*(node_positions(2)-node_positions(1))
    tmp=f(1)
    f(1)=0._f64
    do i=2,N!+1
      f(i)=(f(i)-M)*(node_positions(i+1)-node_positions(i))
      tmp2=f(i)
      f(i)=f(i-1)+tmp
      tmp=tmp2
    enddo    
    f(N+1)=f(N)+tmp
    
    
    !print *,M,f(1),f(N+1) 

  end subroutine function_to_primitive_adv


  subroutine primitive_to_function_adv( &
    f, &
    node_positions, &
    node_positions_back, &
    N, &
    M)
    sll_real64,dimension(:),intent(inout) :: f
    sll_real64,dimension(:),intent(in) :: node_positions
    sll_real64,dimension(:),intent(in) :: node_positions_back
    sll_int32,intent(in):: N
    sll_real64,intent(in)::M
    sll_int32::i
    sll_real64::tmp!,tmp2
    
    tmp=f(1)
    do i=1,N-1
      f(i)=f(i+1)-f(i)+M*(node_positions_back(i+1)-node_positions_back(i))
    enddo
    f(N)=tmp-f(N)+M*(node_positions_back(N+1)-node_positions_back(N))




    !from mean compute f
    do i=1,N
      f(i)=f(i)/(node_positions(i+1)-node_positions(i))
    enddo
    !f(N+1) = f(1)
  end subroutine primitive_to_function_adv



end module sll_advection_2d
