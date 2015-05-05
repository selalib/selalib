!> @ingroup operator_splitting
!> @brief Implements split operators for constant coefficient advection
module sll_split_advection_2d
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_errors.h"
  use sll_module_characteristics_1d_base
  use sll_module_interpolators_1d_base
  use sll_operator_splitting
  use sll_cartesian_meshes  

  implicit none

  sll_int32, parameter :: SLL_ADVECTIVE    = 0
  sll_int32, parameter :: SLL_CONSERVATIVE = 1

  !> @brief 
  !> Simple operator splitting type for 2D  advection
  !> Extends operator splitting
  !> @details This should be
  !> treated as an opaque type. No access to its internals is directly allowed.
  !> comment: how to change then the advection fields in time?
  type, extends(operator_splitting) :: split_advection_2d
     !> interpolator object in first direction
     class(sll_interpolator_1d_base), pointer  :: interp1
     !> interpolator object in second direction
     class(sll_interpolator_1d_base), pointer  :: interp2
     !> characteristics object in first direction
     class(sll_characteristics_1d_base), pointer  :: charac1
     procedure(signature_process_outside_point_1d), pointer, nopass :: process_outside_point1 !< for bdr direction 1
     !> characteristics object in second direction
     class(sll_characteristics_1d_base), pointer  :: charac2
     procedure(signature_process_outside_point_1d), pointer, nopass :: process_outside_point2 !< for bdr direction 1
     !> mesh common for charac and interp
     type(sll_cartesian_mesh_2d), pointer :: mesh_2d
     !> function do be evolved
     sll_real64, dimension(:,:), pointer :: f
     !> dimension in first direction for f
     sll_int32 :: num_dof1
     !> dimension in second direction for f
     sll_int32 :: num_dof2
     !> DOF positions in first direction for f
     sll_real64, dimension(:), pointer  :: dof_positions1
     !> DOF positions in second direction for f
     sll_real64, dimension(:), pointer  :: dof_positions2
     !> advection form (SLL_ADVECTIVE or SLL_CONSERVATIVE)
     sll_int32 :: advection_form
     !> advection coefficient in first direction
     sll_real64, dimension(:,:), pointer  :: a1
     !> advection coefficient in second direction
     sll_real64, dimension(:,:), pointer  :: a2

     !> temporary array
     sll_real64, dimension(:), pointer :: input1
     sll_real64, dimension(:), pointer :: output1
     sll_real64, dimension(:), pointer :: origin1
     sll_real64, dimension(:), pointer :: feet1
     sll_real64, dimension(:), pointer :: feet_inside1
     sll_real64, dimension(:), pointer :: input2
     sll_real64, dimension(:), pointer :: output2
     sll_real64, dimension(:), pointer :: origin2
     sll_real64, dimension(:), pointer :: feet2
     sll_real64, dimension(:), pointer :: feet_inside2
   contains
     procedure, pass(this) :: operatorT => adv1  !< advection in first direction
     procedure, pass(this) :: operatorV => adv2  !< advection in second direction
  end type split_advection_2d


contains
  function new_split_advection_2d( &
      f, &
      a1, &
      a2, &
      interp1, &
      charac1, &
      process_outside_point1, &
      interp2, &
      charac2, &
      process_outside_point2, &
      mesh_2d, &
      advection_form, &
      split_case, &
      split_step, &
      nb_split_step, &
      split_begin_T, &
      dt) &
      result(this)  
    class(split_advection_2d), pointer :: this  !< object to be initialised
    sll_real64, dimension(:,:), pointer, intent(in) :: f   !< initial value of function
    sll_real64, dimension(:,:), pointer, intent(in) :: a1   !< advection coefficient in first direction
    sll_real64, dimension(:,:), pointer, intent(in) :: a2   !< advection coefficient in second direction
    class(sll_interpolator_1d_base), pointer  :: interp1 !< interpolator direction 1
    class(sll_interpolator_1d_base), pointer  :: interp2 !< interpolator direction 1
    class(sll_characteristics_1d_base), pointer  :: charac1 !< characteristics direction 1
    procedure(signature_process_outside_point_1d), pointer :: process_outside_point1 !< for bdr direction 1
    class(sll_characteristics_1d_base), pointer  :: charac2 !< characteristics direction 2
    procedure(signature_process_outside_point_1d), pointer :: process_outside_point2 !< for bdr direction 2
    type(sll_cartesian_mesh_2d), pointer :: mesh_2d !< cartesian mesh common to interp and charac
    sll_int32, intent(in) :: advection_form
    sll_int32, intent(in)  :: split_case  !< defines  splitting method
    sll_real64, dimension(:), intent(in), optional :: split_step  !< coefficients of split step
    sll_int32, intent(in), optional :: nb_split_step !< number of split steps
    logical, intent(in), optional :: split_begin_T   !< begin with operator T if .true.
    sll_real64, intent(in), optional :: dt  !< time step   
    ! local variable
    sll_int32 :: ierr

    SLL_ALLOCATE(this,ierr)   
    call initialize_split_advection_2d( &
      this, &
      f, &
      a1, &
      a2, &
      interp1, &
      charac1, &
      process_outside_point1, &
      interp2, &
      charac2, &
      process_outside_point2, &
      mesh_2d, &
      advection_form, &
      split_case, &
      split_step, &
      nb_split_step, &
      split_begin_T, &
      dt )

  end function

  !> @brief Initialise advection_2d object
  subroutine initialize_split_advection_2d( &
      this, &
      f, &
      a1, &
      a2, &
      interp1, &
      charac1, &
      process_outside_point1, &
      interp2, &
      charac2, &
      process_outside_point2, &
      mesh_2d, &
      advection_form, &
      split_case, &
      split_step, &
      nb_split_step, &
      split_begin_T, &
      dt )
    class(split_advection_2d), intent(inout)   :: this !< object 
    sll_real64, dimension(:,:), pointer, intent(in) :: f   !< initial value of function
    sll_real64, dimension(:,:), pointer, intent(in) :: a1   !< advection coefficient in first direction
    sll_real64, dimension(:,:), pointer, intent(in) :: a2   !< advection coefficient in second direction
    class(sll_interpolator_1d_base), pointer  :: interp1 !< interpolator direction 1
    class(sll_interpolator_1d_base), pointer  :: interp2 !< interpolator direction 1
    class(sll_characteristics_1d_base), pointer  :: charac1 !< characteristics direction 1
    procedure(signature_process_outside_point_1d), pointer :: process_outside_point1 !< for bdr direction 1
    class(sll_characteristics_1d_base), pointer  :: charac2 !< characteristics direction 2
    procedure(signature_process_outside_point_1d), pointer :: process_outside_point2 !< for bdr direction 2
    type(sll_cartesian_mesh_2d), pointer :: mesh_2d !< cartesian mesh common to interp and charac
    sll_int32, intent(in) :: advection_form
    sll_int32, intent(in)  :: split_case  !< defines  splitting method
    sll_real64, dimension(:), intent(in), optional :: split_step  !< coefficients of split step
    sll_int32, intent(in), optional :: nb_split_step !< number of split steps
    logical, intent(in), optional :: split_begin_T   !< begin with operator T if .true.
    sll_real64, intent(in), optional :: dt  !< time step
    sll_int32 :: ierr
    sll_int32 :: n1
    sll_int32 :: n2
    sll_int32 :: i
    sll_real64 :: eta1_min
    sll_real64 :: eta1_max
    sll_real64 :: eta2_min
    sll_real64 :: eta2_max
    

    this%f => f
    this%a1 => a1
    this%a2 => a2
    this%interp1 => interp1
    this%interp2 => interp2
    this%charac1 => charac1
    this%process_outside_point1 => process_outside_point1
    this%process_outside_point2 => process_outside_point2
    this%charac2 => charac2
    this%mesh_2d => mesh_2d
    
    this%advection_form = advection_form
    
    n1 = mesh_2d%num_cells1+1
    n2 = mesh_2d%num_cells2+1
    
    eta1_min = mesh_2d%eta1_min
    eta1_max = mesh_2d%eta1_max
    eta2_min = mesh_2d%eta2_min
    eta2_max = mesh_2d%eta2_max
    
    select case (this%advection_form)
      case (SLL_ADVECTIVE)
        this%num_dof1 = n1
        this%num_dof2 = n2
      case (SLL_CONSERVATIVE)
        this%num_dof1 = n1-1
        this%num_dof2 = n2-1
      case default
        SLL_ERROR("initialize_split_advection_2d", "bad value of advection_form")
        print *,'advection_form',advection_form  
        stop
    end select
    
    if((size(f,1)/= this%num_dof1).or.(size(f,2)/= this%num_dof2))then
      SLL_ERROR("initialize_split_advection_2d", "bad size of f")
      print*,'size of f=',size(f)  
      stop      
    endif
    if((size(a1,1)/= n1).or.(size(a1,2)/= this%num_dof2))then
      SLL_ERROR("initialize_split_advection_2d", "bad size of a1")
      print*,'size of a1=',size(a1)  
      stop      
    endif
    if((size(a2,1)/= this%num_dof1).or.(size(a2,2)/= n2))then
      SLL_ERROR("initialize_split_advection_2d", "bad size of a2")
      print*,'size of a2=',size(a2)  
      stop      
    endif
        
    
    SLL_ALLOCATE(this%input1(n1),ierr)
    SLL_ALLOCATE(this%output1(n1),ierr)
    SLL_ALLOCATE(this%origin1(n1),ierr)
    SLL_ALLOCATE(this%feet1(n1),ierr)
    SLL_ALLOCATE(this%feet_inside1(n1),ierr)

    SLL_ALLOCATE(this%input2(n2),ierr)
    SLL_ALLOCATE(this%output2(n2),ierr)
    SLL_ALLOCATE(this%origin2(n2),ierr)
    SLL_ALLOCATE(this%feet2(n2),ierr)
    SLL_ALLOCATE(this%feet_inside2(n2),ierr)
    
    do i=1,n1
      this%origin1(i) = real(i-1,f64)/real(n1-1,f64)
      this%origin1(i) = &
        eta1_min+this%origin1(i)*(eta1_max-eta1_min)
    enddo
    do i=1,n2
      this%origin2(i) = real(i-1,f64)/real(n2-1,f64)
      this%origin2(i) = &
        eta2_min+this%origin2(i)*(eta2_max-eta2_min)
    enddo
    

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
    class(split_advection_2d), intent(inout) :: this !< object 
    sll_real64, intent(in) :: dt   !< time step
    ! local variables
    sll_int32 :: i
    sll_int32 :: j
    sll_int32 :: n1
    sll_int32 :: n2
    sll_int32 :: num_dof1
    sll_int32 :: num_dof2
    sll_real64 :: eta1_min
    sll_real64 :: eta1_max
    sll_real64 :: mean
    n1 = this%mesh_2d%num_cells1+1
    n2 = this%mesh_2d%num_cells2+1
    num_dof1 =this%num_dof1
    num_dof2 =this%num_dof2
    eta1_min = this%mesh_2d%eta1_min
    eta1_max = this%mesh_2d%eta1_max
        
    do j=1,num_dof2
      this%input1(1:num_dof1) = this%f(1:num_dof1,j)
      if(this%advection_form==SLL_CONSERVATIVE)then      
        call function_to_primitive_adv( &
          this%input1, &
          this%origin1, &
          num_dof1, &
          mean)      
      endif
      
      

      !advection
      call this%charac1%compute_characteristics( &
        this%A1(1:n1,j), &
        dt, &
        this%origin1(1:n1), &
        this%feet1(1:n1))
      do i=1,n1
        this%feet_inside1(i) = this%process_outside_point1( &
          this%feet1(i), &
          eta1_min, &
          eta1_max)
      enddo  


      this%output1(1:n1) = this%interp1%interpolate_array( &
        n1, &
        this%input1(1:n1), &
        this%feet_inside1(1:n1))      
      
      
      if(this%advection_form==SLL_CONSERVATIVE)then      
        call primitive_to_function_adv( &
          this%output1(1:n1), &
          this%origin1(1:n1), &
          this%feet1(1:n1), &
          num_dof1, &
          mean)      
      endif

      this%f(1:num_dof1,j) = this%output1(1:num_dof1)
      
      
        
    enddo


  end subroutine


  !> @brief Advection operator in first direction
  subroutine adv2(this, dt)
    class(split_advection_2d), intent(inout) :: this !< object 
    sll_real64, intent(in) :: dt   !< time step
    ! local variables
    sll_int32 :: i
    sll_int32 :: j
    sll_int32 :: n1
    sll_int32 :: n2
    sll_int32 :: num_dof1
    sll_int32 :: num_dof2
    sll_real64 :: eta2_min
    sll_real64 :: eta2_max
    sll_real64 :: mean

    n1 = this%mesh_2d%num_cells1+1
    n2 = this%mesh_2d%num_cells2+1
    num_dof1 =this%num_dof1
    num_dof2 =this%num_dof2
    eta2_min = this%mesh_2d%eta2_min
    eta2_max = this%mesh_2d%eta2_max
        
    do i=1,num_dof1
      this%input2(1:num_dof2) = this%f(i,1:num_dof2)
      
      if(this%advection_form==SLL_CONSERVATIVE)then      
        call function_to_primitive_adv( &
          this%input2, &
          this%origin2, &
          num_dof2, &
          mean)      
      endif
      
      

      !advection
      call this%charac1%compute_characteristics( &
        this%A2(i,1:n2), &
        dt, &
        this%origin2(1:n2), &
        this%feet2(1:n2))
      
      do j=1,n2
        this%feet_inside2(j) = this%process_outside_point2( &
          this%feet2(j), &
          eta2_min, &
          eta2_max)
      enddo  
      this%output2(1:n2) = this%interp2%interpolate_array( &
        n2, &
        this%input2(1:n2), &
        this%feet_inside2(1:n2))      
      
      
      if(this%advection_form==SLL_CONSERVATIVE)then      
        call primitive_to_function_adv( &
          this%output2(1:n2), &
          this%origin2(1:n2), &
          this%feet2(1:n2), &
          num_dof2, &
          mean)      
      endif

      this%f(i,1:num_dof2) = this%output2(1:num_dof2)
      
      
        
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
    !begin new
    M = M/(node_positions(N+1)-node_positions(1))
    !end new

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



end module sll_split_advection_2d
