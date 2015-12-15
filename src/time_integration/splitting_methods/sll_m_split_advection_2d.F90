!> @ingroup sll_t_operator_splitting
!> @brief Implements split operators for constant coefficient advection
module sll_m_split_advection_2d
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_errors.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_boundary_condition_descriptors, only: &
    sll_p_periodic

  use sll_m_cartesian_meshes, only: &
    sll_t_cartesian_mesh_2d

  use sll_m_characteristics_1d_base, only: &
    sll_i_signature_process_outside_point_1d, &
    sll_c_characteristics_1d_base

  use sll_m_coordinate_transformation_2d_base, only: &
    sll_c_coordinate_transformation_2d_base

  use sll_m_cubic_non_uniform_splines, only: &
    sll_t_cubic_nonunif_spline_1d, &
    sll_f_new_cubic_nonunif_spline_1d

  use sll_m_interpolators_1d_base, only: &
    sll_c_interpolator_1d

  use sll_m_operator_splitting, only: &
    sll_s_initialize_operator_splitting, &
    sll_t_operator_splitting

  implicit none

  public :: &
    sll_f_new_split_advection_2d, &
    sll_p_advective, &
    sll_p_conservative, &
    sll_t_split_advection_2d

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  sll_int32, parameter :: sll_p_advective    = 0
  sll_int32, parameter :: sll_p_conservative = 1

  !> @brief 
  !> Simple operator splitting type for 2D  advection
  !> Extends operator splitting
  !> @details This should be
  !> treated as an opaque type. No access to its internals is directly allowed.
  !> comment: how to change then the advection fields in time?
  type, extends(sll_t_operator_splitting) :: sll_t_split_advection_2d
     !> interpolator object in first direction
     class(sll_c_interpolator_1d), pointer  :: interp1
     !> interpolator object in second direction
     class(sll_c_interpolator_1d), pointer  :: interp2
     !> characteristics object in first direction
     class(sll_c_characteristics_1d_base), pointer  :: charac1
     procedure(sll_i_signature_process_outside_point_1d), pointer, nopass :: process_outside_point1 !< for bdr direction 1
     !> characteristics object in second direction
     class(sll_c_characteristics_1d_base), pointer  :: charac2
     procedure(sll_i_signature_process_outside_point_1d), pointer, nopass :: process_outside_point2 !< for bdr direction 1
     !> mesh common for charac and interp
     type(sll_t_cartesian_mesh_2d), pointer :: mesh_2d
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
     !> advection form (sll_p_advective or sll_p_conservative)
     sll_int32 :: advection_form
     !> advection coefficient in first direction
     sll_real64, dimension(:,:), pointer  :: a1
     !> advection coefficient in second direction
     sll_real64, dimension(:,:), pointer  :: a2

     class(sll_c_coordinate_transformation_2d_base), pointer :: transformation !< coordinate transformation
     type (sll_t_cubic_nonunif_spline_1d), pointer :: spl_eta1
     type (sll_t_cubic_nonunif_spline_1d), pointer :: spl_eta2

     
      
     !> temporary array
     sll_real64, dimension(:), pointer :: input1
     sll_real64, dimension(:), pointer :: output1
     sll_real64, dimension(:), pointer :: origin1
     sll_real64, dimension(:), pointer :: origin_middle1
     sll_real64, dimension(:), pointer :: feet1
     sll_real64, dimension(:), pointer :: feet_middle1
     sll_real64, dimension(:), pointer :: feet_inside1
     sll_real64, dimension(:), pointer :: input2
     sll_real64, dimension(:), pointer :: output2
     sll_real64, dimension(:), pointer :: origin2
     sll_real64, dimension(:), pointer :: origin_middle2
     sll_real64, dimension(:), pointer :: feet2
     sll_real64, dimension(:), pointer :: feet_middle2
     sll_real64, dimension(:), pointer :: feet_inside2

     logical :: csl_2012

     sll_real64, dimension(:), pointer :: primitive1
     sll_real64, dimension(:), pointer :: primitive2
     sll_real64, dimension(:), pointer :: xi1
     sll_real64, dimension(:), pointer :: xi2
     sll_real64, dimension(:), pointer :: A1jac
     sll_real64, dimension(:), pointer :: A2jac
     
     
   contains
     procedure, pass(this) :: operatorT => adv1  !< advection in first direction
     procedure, pass(this) :: operatorV => adv2  !< advection in second direction
  end type sll_t_split_advection_2d


contains
  function sll_f_new_split_advection_2d( &
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
      dt, &
      transformation, &
      csl_2012) &
      result(this)  
    class(sll_t_split_advection_2d), pointer :: this  !< object to be initialised
    sll_real64, dimension(:,:), pointer, intent(in) :: f   !< initial value of function
    sll_real64, dimension(:,:), pointer, intent(in) :: a1   !< advection coefficient in first direction
    sll_real64, dimension(:,:), pointer, intent(in) :: a2   !< advection coefficient in second direction
    class(sll_c_interpolator_1d), pointer  :: interp1 !< interpolator direction 1
    class(sll_c_interpolator_1d), pointer  :: interp2 !< interpolator direction 1
    class(sll_c_characteristics_1d_base), pointer  :: charac1 !< characteristics direction 1
    procedure(sll_i_signature_process_outside_point_1d), pointer :: process_outside_point1 !< for bdr direction 1
    class(sll_c_characteristics_1d_base), pointer  :: charac2 !< characteristics direction 2
    procedure(sll_i_signature_process_outside_point_1d), pointer :: process_outside_point2 !< for bdr direction 2
    type(sll_t_cartesian_mesh_2d), pointer :: mesh_2d !< cartesian mesh common to interp and charac
    sll_int32, intent(in) :: advection_form
    sll_int32, intent(in)  :: split_case  !< defines  splitting method
    sll_real64, dimension(:), intent(in), optional :: split_step  !< coefficients of split step
    sll_int32, intent(in), optional :: nb_split_step !< number of split steps
    logical, intent(in), optional :: split_begin_T   !< begin with operator T if .true.
    sll_real64, intent(in), optional :: dt  !< time step   
    class(sll_c_coordinate_transformation_2d_base), pointer, optional :: transformation !< coordinate transformation
    logical, intent(in), optional :: csl_2012
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
      dt, &
      transformation, &
      csl_2012)

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
      dt, &
      transformation, &
      csl_2012)
    class(sll_t_split_advection_2d), intent(inout)   :: this !< object 
    sll_real64, dimension(:,:), pointer, intent(in) :: f   !< initial value of function
    sll_real64, dimension(:,:), pointer, intent(in) :: a1   !< advection coefficient in first direction
    sll_real64, dimension(:,:), pointer, intent(in) :: a2   !< advection coefficient in second direction
    class(sll_c_interpolator_1d), pointer  :: interp1 !< interpolator direction 1
    class(sll_c_interpolator_1d), pointer  :: interp2 !< interpolator direction 1
    class(sll_c_characteristics_1d_base), pointer  :: charac1 !< characteristics direction 1
    procedure(sll_i_signature_process_outside_point_1d), pointer :: process_outside_point1 !< for bdr direction 1
    class(sll_c_characteristics_1d_base), pointer  :: charac2 !< characteristics direction 2
    procedure(sll_i_signature_process_outside_point_1d), pointer :: process_outside_point2 !< for bdr direction 2
    type(sll_t_cartesian_mesh_2d), pointer :: mesh_2d !< cartesian mesh common to interp and charac
    sll_int32, intent(in) :: advection_form
    sll_int32, intent(in)  :: split_case  !< defines  splitting method
    sll_real64, dimension(:), intent(in), optional :: split_step  !< coefficients of split step
    sll_int32, intent(in), optional :: nb_split_step !< number of split steps
    logical, intent(in), optional :: split_begin_T   !< begin with operator T if .true.
    sll_real64, intent(in), optional :: dt  !< time step
    class(sll_c_coordinate_transformation_2d_base), pointer, optional :: transformation !< coordinate transformation
    logical, intent(in), optional :: csl_2012
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
    
    if(present(transformation))then
      this%transformation => transformation
    endif
    if(present(csl_2012))then
      this%csl_2012 = csl_2012
    else
      this%csl_2012 = .false.  
    endif
    
    n1 = mesh_2d%num_cells1+1
    n2 = mesh_2d%num_cells2+1
    
    eta1_min = mesh_2d%eta1_min
    eta1_max = mesh_2d%eta1_max
    eta2_min = mesh_2d%eta2_min
    eta2_max = mesh_2d%eta2_max
    
    select case (this%advection_form)
      case (sll_p_advective)
        this%num_dof1 = n1
        this%num_dof2 = n2
      case (sll_p_conservative)
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
    SLL_ALLOCATE(this%origin_middle1(n1-1),ierr)
    SLL_ALLOCATE(this%feet_middle1(n1-1),ierr)
    SLL_ALLOCATE(this%feet1(n1),ierr)
    SLL_ALLOCATE(this%feet_inside1(n1),ierr)

    SLL_ALLOCATE(this%input2(n2),ierr)
    SLL_ALLOCATE(this%output2(n2),ierr)
    SLL_ALLOCATE(this%origin2(n2),ierr)
    SLL_ALLOCATE(this%origin_middle2(n2-1),ierr)
    SLL_ALLOCATE(this%feet_middle2(n2-1),ierr)
    SLL_ALLOCATE(this%feet2(n2),ierr)
    SLL_ALLOCATE(this%feet_inside2(n2),ierr)

    SLL_ALLOCATE(this%primitive1(n1),ierr)
    SLL_ALLOCATE(this%primitive2(n2),ierr)
    SLL_ALLOCATE(this%xi1(n1),ierr)
    SLL_ALLOCATE(this%xi2(n2),ierr)
    SLL_ALLOCATE(this%A1jac(n1),ierr)
    SLL_ALLOCATE(this%A2jac(n2),ierr)


    
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

    do i=1,n1-1
      this%origin_middle1(i) = 0.5_f64*(this%origin1(i)+this%origin1(i+1))
    enddo
    do i=1,n2-1
      this%origin_middle2(i) = 0.5_f64*(this%origin2(i)+this%origin2(i+1))
    enddo

    this%spl_eta1 => sll_f_new_cubic_nonunif_spline_1d( &
      n1-1, &
      sll_p_periodic)
    this%spl_eta2 => sll_f_new_cubic_nonunif_spline_1d( n2-1, sll_p_periodic)
    

    call sll_s_initialize_operator_splitting( &
      this, &
      split_case, &
      split_step, &
      nb_split_step, &
      split_begin_T, &
      dt)
  end subroutine 

  !> @brief Advection operator in first direction
  subroutine adv1(this, dt)
    class(sll_t_split_advection_2d), intent(inout) :: this !< object 
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
    !sll_real64 :: eta1
    !sll_real64 :: eta2
    !sll_real64 :: xi_max
    sll_real64 :: delta_eta1
    !sll_real64 :: lperiod
    !sll_real64 :: xi_new
    !sll_real64 :: mean_init
    
    n1 = this%mesh_2d%num_cells1+1
    n2 = this%mesh_2d%num_cells2+1
    num_dof1 =this%num_dof1
    num_dof2 =this%num_dof2
    eta1_min = this%mesh_2d%eta1_min
    eta1_max = this%mesh_2d%eta1_max
    delta_eta1 = (eta1_max-eta1_min)/real(n1-1,f64)
    if(this%csl_2012)then
      
      !Let Fn(eta) = \int_0^{eta} (Jf)(tn,s)ds
      ! such that Fn'(eta_j) = Jf(tn,eta_j)
      !Let F_h(eta) = Fn(H(tn;eta,t_{n+1}))
      !we know F_h(eta_k)=Fn(eta_k^*)
      !we reconstruct F_h
      !Jf(t_{n+1},eta_j) = F_h'(eta_j)
      !we suppose here periodic boundary conditions
      
      ! on [0,1] mesh
      ! eta_k = (k-1/2)/num_dof, k=1,num_dof
      !we have to know
      !  eta_k^*,k=1,num_dof
      !at first order eta_k^*=eta_k-a(eta_k)dt
      !input: a_k, k=1,num_dof
      !      Jf(tn,eta_k), k=1,num_dof
      !      dt
      !      interp1d => new_cubic_spline_interpolator_1d( &
      !        num_dof+1, &
      !        eta_1, &
      !        eta_1+eta_max-eta_min, &
      !        sll_p_periodic)

!      do j=1,num_dof2
!
!
!        eta2 = 0.5_f64*(this%origin2(j)+this%origin2(j+1))
!        do i=1,n1
!          eta1 = this%origin1(i)
!          this%A1jac(i) = &
!            this%A1(i,j)*this%transformation%jacobian(eta1,eta2)
!        enddo
!
!!        do i=1,n1-1
!!          eta1 = 0.5_f64*(this%origin1(i)+this%origin1(i+1))
!!          this%A1jac(i) = &
!!            this%A1(i,j)*this%transformation%jacobian(eta1,eta2)
!!          print *,'A1=', this%A1jac(i)
!!        enddo
!!stop
!
!!!        do i=1,n1
!!          print*,i,this%A1(i,j)
!!        enddo
!!        stop
!!        print *,'#this%f',minval(this%f(1:num_dof1,j)),maxval(this%f(1:num_dof1,j))
!!
!        this%input1(1:num_dof1) = this%f(1:num_dof1,j)
!        eta2 = 0.5_f64*(this%origin2(j)+this%origin2(j+1))
!!        do i=1,num_dof1
!!          eta1 = 0.5_f64*(this%origin1(i)+this%origin1(i+1))
!!          this%input1(i) = &
!!            this%input1(i)/this%transformation%jacobian(eta1,eta2)
!!        enddo
!!        
!        print *,'#input1=',minval(this%input1(1:num_dof1)),maxval(this%input1(1:num_dof1))
!
!        
!        this%primitive1(1) = 0._f64
!        do i=2,num_dof1+1
!          this%primitive1(i) = this%primitive1(i-1) &
!            +delta_eta1*this%input1(i-1)
!        enddo
!
!!        do i=1,n1
!!          print *,this%primitive1(i)
!!        enddo
! 
! 
! 
!        
!        this%xi1(1) = 0._f64
!        do i=2,num_dof1+1
!          eta1 = 0.5_f64*(this%origin1(i)+this%origin1(i-1))
!          this%xi1(i) = this%xi1(i-1) &
!            +delta_eta1*this%transformation%jacobian(eta1,eta2)
!        enddo
!        !jacobian(i1) = df_jac_at_i( i1-1, i2 )
!
!        mean = this%primitive1(n1)/this%xi1(n1)
!        print *,'#mean init'
!        !modify primitive so that it becomes periodic
!        do i = 2,num_dof1+1
!          this%primitive1(i) = this%primitive1(i)-mean*this%xi1(i)
!        end do
!        xi_max = this%xi1(n1)
!
!!        print *,'#this%primitive1=',minval(this%primitive1(1:n1)),maxval(this%primitive1(1:n1))
!!        do i=1,n1
!!          print *,this%primitive1(i)
!!        enddo
!!        stop
!        call implicit_ode_nonuniform( &
!          2, &
!          dt, &
!          this%xi1, &
!          num_dof1, &
!          PERIODIC_ODE, &
!          this%feet_inside1, &
!          this%A1jac(1:n1),&
!          this%A1jac)
!        
!        do i=1,n1
!          this%feet1(i) = this%xi1(i)-this%A1jac(i)*dt
!        enddo  
!           
!        call sll_s_compute_spline_nonunif( &
!          this%primitive1, &
!          this%spl_eta1, &
!          this%xi1)
!        ! interpolate primitive at origin of characteritics
!        call sll_s_interpolate_array_value_nonunif( &
!          this%feet_inside1, &
!          this%primitive1, &
!          n1, &
!          this%spl_eta1)
!        ! come back to real primitive by adding average
!        if (this%feet1(1) > 0.5_f64*(xi_max)) then
!          lperiod = -1.0_f64
!        else
!          lperiod = 0.0_f64
!        end if
!        !print *,'lperiod=',lperiod
!        xi_new = this%feet1(1) +  lperiod*xi_max
!        this%primitive1(1) = this%primitive1(1) + mean*xi_new
!        !if ((xi_new > xi_max) .or. (xi_new <xi(1))) then
!        !   print*, 1, xi_new, xi_out(1), primitive(1)
!        !end if
!        do i = 2,n1
!          ! We need here to find the points where it has been modified by periodicity
!          if (this%feet1(i) < this%feet1(i-1)) then
!             lperiod = lperiod+1.0_f64
!          end if
!          !print *,'lperiod=',i,lperiod
!          xi_new = this%feet1(i)+lperiod*xi_max
!          this%primitive1(i) = this%primitive1(i)+mean*xi_new
!          !if (i>98) then
!          !   print*, 'iii', i, xi_new, xi_out(i),xi_max, primitive(i), primitive(i-1)
!          !endif
!        end do
!        
!        do i=1,num_dof1
!          eta1 = 0.5_f64*(this%origin1(i)+this%origin1(i+1))
!          !print *,this%feet1(i+1)-this%feet1(i), &
!          !  this%xi1(i+1)-this%xi1(i), &
!          !  this%primitive1(i+1)-this%primitive1(i), &
!          !  this%transformation%jacobian(eta1,eta2)*delta_eta1
!        enddo
!        
!        do i = 1,num_dof1 
!          this%f(i,j) = (this%primitive1(i+1)-this%primitive1(i))/delta_eta1
!          eta1 = 0.5_f64*(this%origin1(i)+this%origin1(i+1))
!          !print *,i,this%f(i,j)/this%transformation%jacobian(eta1,eta2)
!          !call sll_set_df_val( dist_func_2D, i1, i2, val )
!          !if (val/df_jac_at_i(i1,i2)>1.) then
!          !   print*, 'val', i1,i2, val, primitive1(i1) , primitive1(i1+1), df_jac_at_i(i1,i2), delta_eta1
!          !end if
!       end do
!
!       !
!!     print *,'#output1=',minval(this%f(1:num_dof1,j)),maxval(this%f(1:num_dof1,j))
!!       print *,'#mean=',mean
!!       !print *,'#this%f',minval(this%f(1:num_dof1,j)),maxval(this%f(1:num_dof1,j))
!!       !stop
!
!        
!      enddo



    else
    
        
    do j=1,num_dof2
      this%input1(1:num_dof1) = this%f(1:num_dof1,j)

!      print *,'this%input1(1:num_dof1)=', &
!        minval(this%input1(1:num_dof1)), &
!        maxval(this%input1(1:num_dof1))


      if(this%advection_form==sll_p_conservative)then      
        call function_to_primitive_adv( &
          this%input1, &
          this%origin1, &
          num_dof1, &
          mean)
!        this%input1(1:num_dof1) = this%f(1:num_dof1,j)  
!        eta2 = 0.5_f64*(this%origin2(j)+this%origin2(j+1))
!        do i=1,num_dof1
!          eta1 = 0.5_f64*(this%origin1(i)+this%origin1(i+1))
!          this%input1(i) = &
!            this%input1(i)/this%transformation%jacobian(eta1,eta2)-mean
!          this%input1(i) = &
!            this%input1(i)*this%transformation%jacobian(eta1,eta2)
!        enddo        
!        mean_init = mean
!        call function_to_primitive_adv( &
!          this%input1, &
!          this%origin1, &
!          num_dof1, &
!          mean)
!        !print *,'new mean=',mean
!        !stop  
      endif
      

      !advection
!      if(this%advection_form==sll_p_conservative)then      
!        do i=1,num_dof1
!          this%A1jac(i) = 0.5_f64*(this%A1(i,j)+this%A1(i+1,j))
!        enddo
!        call this%charac1%compute_characteristics( &
!          this%A1jac(1:num_dof1), &
!          dt, &
!          this%origin_middle1(1:num_dof1), &
!          this%feet_middle1(1:num_dof1))
!      
!        this%feet1(1) = 0.5_f64*(this%feet_middle1(1)+this%feet_middle1(n1-1)) &
!         -0.5_f64*(this%origin1(n1)-this%origin1(1))
!        
!        do i=2,n1-1
!          this%feet1(i) = 0.5_f64*(this%feet_middle1(i)+this%feet_middle1(i-1))
!        enddo
!        
!        this%feet1(n1) = this%feet1(1) &
!         +(this%origin1(n1)-this%origin1(1))
!        
!      else
        call this%charac1%compute_characteristics( &
          this%A1(1:n1,j), &
          dt, &
          this%origin1(1:n1), &
          this%feet1(1:n1))
!      endif
      
        
      do i=1,n1
        this%feet_inside1(i) = this%process_outside_point1( &
          this%feet1(i), &
          eta1_min, &
          eta1_max)
      enddo  


      call this%interp1%interpolate_array( &
        n1, &
        this%input1(1:n1), &
        -this%feet_inside1(1:n1),&
        this%output1(1:n1))      
      
      
      if(this%advection_form==sll_p_conservative)then      
        call primitive_to_function_adv( &
          this%output1(1:n1), &
          this%origin1(1:n1), &
          this%feet1(1:n1), &
          num_dof1, &
          mean)
!        eta2 = 0.5_f64*(this%origin2(j)+this%origin2(j+1))
!        do i=1,num_dof1
!          eta1 = 0.5_f64*(this%origin1(i)+this%origin1(i+1))
!          this%output1(i) = &
!            this%output1(i)/this%transformation%jacobian(eta1,eta2)+mean_init
!          this%output1(i) = &
!            this%output1(i)*this%transformation%jacobian(eta1,eta2)
!        enddo        
                
      endif

      this%f(1:num_dof1,j) = this%output1(1:num_dof1)
      
!      print *,'this%output1(1:num_dof1)=', &
!        minval(this%output1(1:num_dof1)), &
!        maxval(this%output1(1:num_dof1))
!      
!      print *,'#mean_adv1=',mean
      
!      if(this%advection_form==sll_p_conservative)then      
!        call function_to_primitive_adv( &
!          this%output1, &
!          this%origin1, &
!          num_dof1, &
!          mean)      
!      endif
!      
!      print *,'#mean_adv1b=',mean
        
    enddo
    endif

  end subroutine


  !> @brief Advection operator in second direction
  subroutine adv2(this, dt)
    class(sll_t_split_advection_2d), intent(inout) :: this !< object 
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
    !sll_real64 :: mean_init
    !sll_real64 :: eta1
    !sll_real64 :: eta2

    n1 = this%mesh_2d%num_cells1+1
    n2 = this%mesh_2d%num_cells2+1
    num_dof1 =this%num_dof1
    num_dof2 =this%num_dof2
    eta2_min = this%mesh_2d%eta2_min
    eta2_max = this%mesh_2d%eta2_max
    
    !return     
    do i=1,num_dof1
      this%input2(1:num_dof2) = this%f(i,1:num_dof2)
      !print *,'#input2=',minval(this%input2(1:num_dof2)),maxval(this%input2(1:num_dof2))
     
      if(this%advection_form==sll_p_conservative)then      
        call function_to_primitive_adv( &
          this%input2, &
          this%origin2, &
          num_dof2, &
          mean)      

!        this%input2(1:num_dof1) = this%f(i,1:num_dof2)  
!        eta1 = 0.5_f64*(this%origin1(i)+this%origin1(i+1))
!        do j=1,num_dof2
!          eta2 = 0.5_f64*(this%origin2(j)+this%origin2(j+1))
!          this%input2(j) = &
!            this%input2(j)/this%transformation%jacobian(eta1,eta2)-mean
!          this%input2(j) = &
!            this%input2(j)*this%transformation%jacobian(eta1,eta2)
!        enddo        
!        mean_init = mean
!       call function_to_primitive_adv( &
!          this%input2, &
!          this%origin2, &
!          num_dof2, &
!          mean)      


      endif
      
      

      !advection
!      if(this%advection_form==sll_p_conservative)then      
!        do j=1,num_dof2
!          this%A2jac(j) = 0.5_f64*(this%A2(i,j)+this%A2(i,j+1))
!        enddo
!        call this%charac2%compute_characteristics( &
!          this%A2jac(1:num_dof2), &
!          dt, &
!          this%origin_middle2(1:num_dof2), &
!          this%feet_middle2(1:num_dof2))
!      
!        this%feet2(1) = 0.5_f64*(this%feet_middle2(1)+this%feet_middle2(n2-1)) &
!         -0.5_f64*(this%origin2(n2)-this%origin2(1))
!        
!        do j=2,n2-1
!          this%feet2(j) = 0.5_f64*(this%feet_middle2(j)+this%feet_middle2(j-1))
!        enddo
!        
!        this%feet2(n2) = this%feet2(1) &
!         +(this%origin2(n2)-this%origin2(1))
!        
!      else
        call this%charac2%compute_characteristics( &
        this%A2(i,1:n2), &
        dt, &
        this%origin2(1:n2), &
        this%feet2(1:n2))
!      endif
      !print *,'err=',maxval(this%origin2-this%feet2)
      
      do j=1,n2
        this%feet_inside2(j) = this%process_outside_point2( &
          this%feet2(j), &
          eta2_min, &
          eta2_max)
      enddo  
      call this%interp2%interpolate_array( &
           n2, &
           this%input2(1:n2), &
           -this%feet_inside2(1:n2), &
           this%output2(1:n2))      
      
      
      if(this%advection_form==sll_p_conservative)then      
        call primitive_to_function_adv( &
          this%output2(1:n2), &
          this%origin2(1:n2), &
          this%feet2(1:n2), &
          num_dof2, &
          mean)      
!        eta1 = 0.5_f64*(this%origin1(i)+this%origin1(i+1))
!        do j=1,num_dof2
!          eta2 = 0.5_f64*(this%origin2(j)+this%origin2(j+1))
!          this%output2(j) = &
!            this%output2(j)/this%transformation%jacobian(eta1,eta2)+mean_init
!          this%output2(j) = &
!            this%output2(j)*this%transformation%jacobian(eta1,eta2)
!        enddo        


      endif

      this%f(i,1:num_dof2) = this%output2(1:num_dof2)
!      print *,'#output2=',minval(this%output2(1:num_dof2)),maxval(this%output2(1:num_dof2))
!      
!      print *,'#mean_adv2=',mean
!      stop
        
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



end module sll_m_split_advection_2d
