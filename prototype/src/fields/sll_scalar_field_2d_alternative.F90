!------------------------------------------------------------------------------
! SELALIB
!------------------------------------------------------------------------------
!
! MODULE: sll_scalar_field_2d
!
!> @author
!> - Edwin
!> - Pierre
!> - Eric
!> - Aurore
!
! DESCRIPTION: 
!
!> @brief
!> Implements the geometry and mesh descriptor types
!>
!>@details
!>
!> This module depends on:
!>    - memory
!>    - precision
!>    - assert
!>    - utilities
!>    - constants
!
! REVISION HISTORY:
! DD Mmm YYYY - Initial Version
! TODO_dd_mmm_yyyy - TODO_describe_appropriate_changes - TODO_name
!------------------------------------------------------------------------------
module sll_module_scalar_field_2d_alternative
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_file_io.h"
  use sll_module_scalar_field_2d_base
  use sll_constants
  use sll_module_interpolators_2d_base
  use sll_arbitrary_degree_spline_interpolator_2d_module
  use sll_utilities
  use sll_boundary_condition_descriptors
  use sll_gnuplot
!  use sll_scalar_field_initializers_base
  implicit none

  type, extends(sll_scalar_field_2d_base) :: sll_scalar_field_2d_analytic_alt
     procedure(two_var_parametrizable_function), pointer, nopass :: func
     procedure(two_var_parametrizable_function), pointer, nopass :: first_deriv_eta1
     procedure(two_var_parametrizable_function), pointer, nopass :: first_deriv_eta2
     sll_real64, dimension(:), pointer        :: params
     character(len=64)                        :: name
     !     sll_int32                                :: plot_counter
     class(sll_coordinate_transformation_2d_base), pointer :: T
     sll_int32 :: bc_left
     sll_int32 :: bc_right
     sll_int32 :: bc_bottom
     sll_int32 :: bc_top
     ! allows to decide if the user put the derivative of the analiytic function: func
     logical :: present_deriv_eta1_int
     logical :: present_deriv_eta2_int
   contains
     procedure, pass(field) :: initialize => &
          initialize_scalar_field_2d_analytic_alt
     
     procedure, pass(field) :: get_transformation => &
          get_transformation_analytic_alt
     procedure, pass(field) :: get_logical_mesh => &
          get_logical_mesh_2d_analytic_alt
     procedure, pass(field) :: get_jacobian_matrix => &
          get_jacobian_matrix_analytic_alt
     procedure, pass(field) :: value_at_point => value_at_pt_analytic
     procedure, pass(field) :: value_at_indices => value_at_index_analytic
     procedure, pass(field) :: first_deriv_eta1_value_at_point => &
          first_deriv_eta1_value_at_pt_analytic
     procedure, pass(field) :: first_deriv_eta2_value_at_point => &
          first_deriv_eta2_value_at_pt_analytic
     procedure, pass(field) :: first_deriv_eta1_value_at_indices => &
          first_deriv_eta1_value_at_index_analytic
     procedure, pass(field) :: first_deriv_eta2_value_at_indices => &
          first_deriv_eta2_value_at_index_analytic
     procedure, pass(field) :: set_field_data => set_field_data_analytic_2d
     procedure, pass(field) :: update_interpolation_coefficients => &
          update_interpolation_coefficients_2d_analytic
     procedure, pass(field) :: write_to_file => write_to_file_analytic_2d
     procedure, pass(field) :: delete => delete_field_2d_analytic_alt
  end type sll_scalar_field_2d_analytic_alt

  type, extends(sll_scalar_field_2d_base) :: sll_scalar_field_2d_discrete_alt
     sll_real64, dimension(:,:), pointer  :: values => null()
     !sll_real64, dimension(:,:), pointer  :: coeff_spline
     !sll_int32                            :: sz_coeff1
     !sll_int32                            :: sz_coeff2
     character(len=64)                    :: name
!     sll_int32                            :: plot_counter
     class(sll_coordinate_transformation_2d_base), pointer :: T
     class(sll_interpolator_2d_base), pointer :: interp_2d
     sll_real64, dimension(:), pointer :: point1_1d
     sll_real64, dimension(:), pointer :: point2_1d
     !sll_real64, dimension(:,:), pointer :: point2d
     sll_int32 :: bc_left
     sll_int32 :: bc_right
     sll_int32 :: bc_bottom
     sll_int32 :: bc_top
   contains
     procedure, pass(field) :: initialize => &
          initialize_scalar_field_2d_discrete_alt
     procedure, pass(field) :: update_interpolation_coefficients => &
          update_interp_coeffs_2d_discrete
     procedure, pass(field) :: get_transformation => &
          get_transformation_discrete_alt
     procedure, pass(field) :: get_logical_mesh => &
          get_logical_mesh_2d_discrete_alt
     procedure, pass(field) :: get_jacobian_matrix => &
          get_jacobian_matrix_discrete_alt
     procedure, pass(field) :: value_at_point => value_at_pt_discrete
     procedure, pass(field) :: value_at_indices => value_at_index_discrete
     procedure, pass(field) :: first_deriv_eta1_value_at_point => &
          first_deriv_eta1_value_at_pt_discrete
     procedure, pass(field) :: first_deriv_eta2_value_at_point => &
          first_deriv_eta2_value_at_pt_discrete
     procedure, pass(field) :: first_deriv_eta1_value_at_indices => &
          first_deriv_eta1_value_at_index_discrete
     procedure, pass(field) :: first_deriv_eta2_value_at_indices => &
          first_deriv_eta2_value_at_index_discrete
     procedure, pass(field) :: set_field_data => set_field_data_discrete_2d
     procedure, pass(field) :: write_to_file => write_to_file_discrete_2d
     procedure, pass(field) :: delete => delete_field_2d_discrete_alt
  end type sll_scalar_field_2d_discrete_alt

!  type sll_ptr_scalar_field_2d_a
!!$ type, extends(sll_scalar_field_2d_base) :: scalar_field_2d_dis
!!$     class(sll_interpolator_1d_base), pointer :: eta1_interpolator
!!$     class(sll_interpolator_1d_base), pointer :: eta2_interpolator
!!$     sll_real64, dimension(:,:), pointer      :: data
!!$     sll_int32                                :: data_position
!!$     character(len=64)                        :: name
!!$     sll_int32                                :: plot_counter
!!$  end type scalar_field_2d_dis


  abstract interface
     function two_var_parametrizable_function( eta1, eta2, params )
       use sll_working_precision
       sll_real64 :: two_var_parametrizable_function
       sll_real64, intent(in) :: eta1
       sll_real64, intent(in) :: eta2
       sll_real64, dimension(:), intent(in) :: params
     end function two_var_parametrizable_function
  end interface

  abstract interface
     function scalar_function_2D( eta1, eta2 )
       use sll_working_precision
       sll_real64 :: scalar_function_2D
       sll_real64, intent(in)  :: eta1
       sll_real64, intent(in)  :: eta2
     end function scalar_function_2D
  end interface

  interface delete
     module procedure delete_field_2d_analytic_alt, delete_field_2d_discrete_alt
  end interface delete


contains   ! *****************************************************************




  ! **************************************************************************
  !
  !                         ANALYTIC CASE
  !
  ! **************************************************************************


  function value_at_pt_analytic( field, eta1, eta2 )
    class(sll_scalar_field_2d_analytic_alt), intent(in) :: field
    sll_real64, intent(in) :: eta1
    sll_real64, intent(in) :: eta2
    sll_real64             ::  value_at_pt_analytic
    value_at_pt_analytic = field%func(eta1,eta2,field%params)
  end function value_at_pt_analytic

  function value_at_index_analytic( field, i, j )
    class(sll_scalar_field_2d_analytic_alt), intent(in) :: field
    sll_int32, intent(in) :: i
    sll_int32, intent(in) :: j
    sll_real64            :: eta1
    sll_real64            :: eta2
    sll_real64            :: value_at_index_analytic
    eta1 = field%T%mesh%eta1_min + real(i-1,f64)*field%T%mesh%delta_eta1
    eta2 = field%T%mesh%eta2_min + real(j-1,f64)*field%T%mesh%delta_eta2
    value_at_index_analytic = field%func(eta1,eta2,field%params)
  end function value_at_index_analytic

  function first_deriv_eta1_value_at_pt_analytic( field, eta1, eta2)
    class(sll_scalar_field_2d_analytic_alt), intent(in) :: field
    sll_real64, intent(in) :: eta1
    sll_real64, intent(in) :: eta2
    sll_real64             :: first_deriv_eta1_value_at_pt_analytic
    
    if ( field%present_deriv_eta1_int ) then 
       first_deriv_eta1_value_at_pt_analytic = &
            field%first_deriv_eta1(eta1,eta2,field%params)
    else 
       print*,field%name, &
            'first_deriv_eta1_value_at_pt_analytic(), ERROR: ', &
            ': first derivative in eta1 is not given in the initialization'
    end if

  end function first_deriv_eta1_value_at_pt_analytic

  function first_deriv_eta2_value_at_pt_analytic( field, eta1, eta2)
    class(sll_scalar_field_2d_analytic_alt), intent(in) :: field
    sll_real64, intent(in) :: eta1
    sll_real64, intent(in) :: eta2
    sll_real64            :: first_deriv_eta2_value_at_pt_analytic
    
    if ( field%present_deriv_eta2_int ) then 
       first_deriv_eta2_value_at_pt_analytic = &
            field%first_deriv_eta2(eta1,eta2,field%params)
    else 
       print*, field%name, &
            'first_deriv_eta2_value_at_pt_analytic(), ERROR: ', &
            ': first derivative in eta2 is not given in the initialization'
    end if
    
  end function first_deriv_eta2_value_at_pt_analytic
  
  function first_deriv_eta1_value_at_index_analytic( field, i, j)
    class(sll_scalar_field_2d_analytic_alt), intent(in) :: field
    sll_int32, intent(in) :: i
    sll_int32, intent(in) :: j
    sll_real64            :: eta1
    sll_real64            :: eta2
    sll_real64            :: first_deriv_eta1_value_at_index_analytic
    
    eta1 = field%T%mesh%eta1_min + real(i-1,f64)*field%T%mesh%delta_eta1
    eta2 = field%T%mesh%eta2_min + real(j-1,f64)*field%T%mesh%delta_eta2
    
    if ( field%present_deriv_eta1_int ) then 
       first_deriv_eta1_value_at_index_analytic = &
            field%first_deriv_eta1(eta1,eta2,field%params)
    else 
       print*,field%name, &
            'first_deriv_eta1_value_at_index_analytic(): ERROR, ', &
            'first derivative in eta1 is not given in the initialization'
    end if
    
  end function first_deriv_eta1_value_at_index_analytic
  
  function first_deriv_eta2_value_at_index_analytic( field, i, j)
    class(sll_scalar_field_2d_analytic_alt), intent(in) :: field
    sll_int32, intent(in) :: i
    sll_int32, intent(in) :: j
    sll_real64            :: eta1
    sll_real64            :: eta2
    sll_real64            :: first_deriv_eta2_value_at_index_analytic
    
    eta1 = field%T%mesh%eta1_min + real(i-1,f64)*field%T%mesh%delta_eta1
    eta2 = field%T%mesh%eta2_min + real(j-1,f64)*field%T%mesh%delta_eta2
    
    if ( field%present_deriv_eta2_int ) then 
       first_deriv_eta2_value_at_index_analytic = &
            field%first_deriv_eta2(eta1,eta2,field%params)
    else 
       print*,' first derivative in eta2 is not given in the initialization'
    end if
    
  end function first_deriv_eta2_value_at_index_analytic
!!$  subroutine initialize_scalar_field_2d_analytic_alt( &
!!$    obj, &
!!$    func, &
!!$    field_name, &
!!$    transformation, &
!!$    bc_left, &
!!$    bc_right, &
!!$    bc_bottom, &
!!$    bc_top, &
!!$    func_params )
!!$
!!$    class(sll_scalar_field_2d_analytic_alt)         :: obj
!!$    procedure(two_var_parametrizable_function)      :: func
!!$    character(len=*), intent(in)                    :: field_name
!!$    sll_real64, dimension(:), intent(in), optional, target :: func_params
!!$    class(sll_coordinate_transformation_2d_base), target :: transformation
!!$    sll_int32, intent(in) :: bc_left
!!$    sll_int32, intent(in) :: bc_right
!!$    sll_int32, intent(in) :: bc_bottom
!!$    sll_int32, intent(in) :: bc_top
!!$    sll_int32  :: ierr
!!$ 
!!$    obj%T => transformation
!!$    !    obj%mesh%written = .false.
!!$    obj%func      => func
!!$    obj%params    => func_params   
!!$    obj%name      = trim(field_name)
!!$    obj%bc_left   = bc_left
!!$    obj%bc_right  = bc_right
!!$    obj%bc_bottom = bc_bottom
!!$    obj%bc_top    = bc_top
!!$  end subroutine initialize_scalar_field_2d_analytic_alt

  function new_scalar_field_2d_analytic_alt( &
    func, &
    field_name, &
    transformation, &
    bc_left, &
    bc_right, &
    bc_bottom, &
    bc_top, &
    func_params,&
    first_deriv_eta1,&
    first_deriv_eta2) result(obj)
    
    type(sll_scalar_field_2d_analytic_alt), pointer :: obj
    procedure(two_var_parametrizable_function)      :: func
    character(len=*), intent(in)                    :: field_name
    class(sll_coordinate_transformation_2d_base), target :: transformation
    sll_int32, intent(in) :: bc_left
    sll_int32, intent(in) :: bc_right
    sll_int32, intent(in) :: bc_bottom
    sll_int32, intent(in) :: bc_top
    sll_real64, dimension(:), intent(in) :: func_params
    procedure(two_var_parametrizable_function), optional :: first_deriv_eta1
    procedure(two_var_parametrizable_function), optional :: first_deriv_eta2
    sll_int32  :: ierr
 
    SLL_ALLOCATE(obj,ierr)
    call obj%initialize( &
    func, &
    field_name, &
    transformation, &
    bc_left, &
    bc_right, &
    bc_bottom, &
    bc_top, &
    func_params,&
    first_deriv_eta1,&
    first_deriv_eta2)
  end function new_scalar_field_2d_analytic_alt

  subroutine set_field_data_analytic_2d( field, values )
    class(sll_scalar_field_2d_analytic_alt), intent(inout) :: field
    sll_real64, dimension(:,:), intent(in) :: values
    print *, 'WARNING: set_field_data_analytic_2d(): it is useless to ', &
         'call this function on an analytic scalar field.'
  end subroutine set_field_data_analytic_2d

  subroutine update_interpolation_coefficients_2d_analytic( field )
    class(sll_scalar_field_2d_analytic_alt), intent(inout) :: field
    print *, 'WARNING: update_interpolation_coefficients_2d_analytic(): ', &
         ' it is useless to call this function on an analytic scalar field.'
  end subroutine update_interpolation_coefficients_2d_analytic

  subroutine delete_field_2d_analytic_alt( field )
    class(sll_scalar_field_2d_analytic_alt), intent(out) :: field
    ! nothing internal do deallocate, just nullify pointers. Can't call
    ! delete on them because the field does not 'own' these data.
    nullify(field%func)
    nullify(field%params)
    nullify(field%T)
  end subroutine delete_field_2d_analytic_alt

  ! For those cases in which handling pointers to field structures is not
  ! convenient, we offer the following alternative initialization.
  subroutine initialize_scalar_field_2d_analytic_alt( &
    field, &
    func, &
    field_name, &
    transformation, &
    bc_left, &
    bc_right, &
    bc_bottom, &
    bc_top, &
    func_params, &
    first_deriv_eta1,&
    first_deriv_eta2)

    class(sll_scalar_field_2d_analytic_alt), intent(out) :: field
    procedure(two_var_parametrizable_function)           :: func
    character(len=*), intent(in)                         :: field_name
    class(sll_coordinate_transformation_2d_base), target :: transformation
    sll_int32, intent(in) :: bc_left
    sll_int32, intent(in) :: bc_right
    sll_int32, intent(in) :: bc_bottom
    sll_int32, intent(in) :: bc_top
    sll_real64, dimension(:), intent(in)  :: func_params
    procedure(two_var_parametrizable_function), optional :: first_deriv_eta1
    procedure(two_var_parametrizable_function), optional :: first_deriv_eta2
    sll_int32 :: ierr
 
    field%T => transformation
    !    field%mesh%written = .false.
    field%func      => func
    SLL_ALLOCATE(field%params(size(func_params)),ierr)
    field%params(:) = func_params
    field%name      = trim(field_name)
    field%bc_left   = bc_left
    field%bc_right  = bc_right
    field%bc_bottom = bc_bottom
    field%bc_top    = bc_top
    
    if (present(first_deriv_eta1)) then
       field%first_deriv_eta1 => first_deriv_eta1
       field%present_deriv_eta1_int = .TRUE.
    end if
    if (present(first_deriv_eta2)) then
       field%first_deriv_eta2 => first_deriv_eta2
       field%present_deriv_eta2_int = .TRUE.
    end if
  end subroutine initialize_scalar_field_2d_analytic_alt

  
  ! The following pair of subroutines are tricky. We want them as general 
  ! services by the fields, hence we need this subroutine interface, yet
  ! we would also like a flexibility in how the derivatives are computed.
  ! A general interpolator interface would cover most of the cases, maybe
  ! all. It could be that a finite difference scheme would also work, if
  ! we ignore some of the interpolator services, like the ability to return
  ! values anywhere instead of at the nodes.
  ! For now, this interface would permit to have multiple implementations.
!!$  subroutine compute_eta1_derivative_on_col( field2d, ith_col, deriv_out )
!!$    type(scalar_field_2d), intent(in)    :: field2d
!!$    sll_int32, intent(in)                :: ith_col
!!$    sll_real64, dimension(:),intent(out) :: deriv_out
!!$  end subroutine compute_eta1_derivative_on_col

  function get_transformation_analytic_alt( field ) result(res)
    class(sll_scalar_field_2d_analytic_alt), intent(in) :: field
    class(sll_coordinate_transformation_2d_base), pointer :: res
    res => field%T
  end function get_transformation_analytic_alt


  function get_logical_mesh_2d_analytic_alt( field ) result(res)
    class(sll_scalar_field_2d_analytic_alt), intent(in) :: field
    type(sll_logical_mesh_2d), pointer :: res
    res => field%T%mesh
  end function get_logical_mesh_2d_analytic_alt

  function get_jacobian_matrix_analytic_alt( field, eta1, eta2 ) result(res)
    class(sll_scalar_field_2d_analytic_alt), intent(in) :: field
    sll_real64, intent(in) :: eta1
    sll_real64, intent(in) :: eta2
    sll_real64, dimension(2,2) :: res
    res = (field%T%jacobian_matrix(eta1,eta2))
  end function get_jacobian_matrix_analytic_alt


  subroutine write_to_file_analytic_2d( field, tag )
    class(sll_scalar_field_2d_analytic_alt), intent(in) :: field
    sll_int32, intent(in)                               :: tag
    sll_int32 :: nptsx1
    sll_int32 :: nptsx2
    sll_real64, dimension(:,:), allocatable :: x1coords
    sll_real64, dimension(:,:), allocatable :: x2coords
    sll_real64, dimension(:,:), allocatable :: values
    sll_real64                              :: eta1
    sll_real64                              :: eta2
    class(sll_coordinate_transformation_2d_base), pointer :: T
    type(sll_logical_mesh_2d), pointer      :: mesh
    sll_int32 :: i
    sll_int32 :: j
    sll_int32 :: ierr

    ! use the logical mesh information to find out the extent of the
    ! domain and allocate the arrays for the plotter.
    T      => field%get_transformation()
    mesh   => field%get_logical_mesh()
    nptsx1 = mesh%num_cells1 + 1
    nptsx2 = mesh%num_cells2 + 1
    SLL_ALLOCATE(x1coords(nptsx1,nptsx2),ierr)
    SLL_ALLOCATE(x2coords(nptsx1,nptsx2),ierr)
    SLL_ALLOCATE(values(nptsx1,nptsx2),ierr)

    ! Fill the arrays with the needed information.
    do j=1, nptsx2
       eta2 = mesh%eta2_min + (j-1)*mesh%delta_eta2 
       do i=1, nptsx1
          eta1 = mesh%eta1_min + (i-1)*mesh%delta_eta1
          x1coords(i,j) = T%x1(eta1,eta2)
          x2coords(i,j) = T%x2(eta1,eta2)
          values(i,j)   = field%value_at_point(eta1,eta2)
       end do
    end do

    call sll_gnuplot_curv_2d( &
         nptsx1, &
         nptsx2, &
         x1coords, &
         x2coords, &
         values, &
         trim(field%name), &
         tag, &
         ierr )

    SLL_DEALLOCATE_ARRAY(x1coords,ierr)
    SLL_DEALLOCATE_ARRAY(x2coords,ierr)
    SLL_DEALLOCATE_ARRAY(values,ierr)
  end subroutine write_to_file_analytic_2d


  ! **************************************************************************
  !
  !                         DISCRETE CASE
  !
  ! **************************************************************************

  function new_scalar_field_2d_discrete_alt( &
    field_name, &
    interpolator_2d, &
    transformation, &
    bc_left, &
    bc_right, &
    bc_bottom, &
    bc_top,&
    point1_1d, &
    sz_point1,&
    point2_1d,&
    sz_point2) result(obj)
    ! point2d)
   ! result(obj)!

    type(sll_scalar_field_2d_discrete_alt), pointer :: obj
    character(len=*), intent(in)                    :: field_name
    class(sll_interpolator_2d_base), target        :: interpolator_2d
    class(sll_coordinate_transformation_2d_base), target :: transformation
    sll_int32 :: SPLINE_DEG1
    sll_int32 :: SPLINE_DEG2
    sll_int32, intent(in) :: bc_left
    sll_int32, intent(in) :: bc_right
    sll_int32, intent(in) :: bc_bottom
    sll_int32, intent(in) :: bc_top
    sll_real64, dimension(:), optional :: point1_1d
    sll_real64, dimension(:), optional :: point2_1d
    sll_int32, optional :: sz_point1
    sll_int32, optional :: sz_point2
    ! sll_real64, dimension(:,:), optional :: point2d
    sll_int32  :: ierr
    
    SLL_ALLOCATE(obj,ierr)
    call obj%initialize( &
         field_name, &
         interpolator_2d, &
         transformation, &
         bc_left, &
         bc_right, &
         bc_bottom, &
         bc_top,&
         point1_1d,&
         sz_point1,&
         point2_1d,&
         sz_point2)
  end function new_scalar_field_2d_discrete_alt
  
  subroutine initialize_scalar_field_2d_discrete_alt( &
    field, &
    field_name, &
    interpolator_2d, &
    transformation, &
    bc_left, &
    bc_right, &
    bc_bottom, &
    bc_top,& 
    point1_1d, &
    sz_point1,&
    point2_1d,&
    sz_point2)
    
    class(sll_scalar_field_2d_discrete_alt)              :: field
    character(len=*), intent(in)                         :: field_name
    class(sll_interpolator_2d_base), target              :: interpolator_2d
    class(sll_coordinate_transformation_2d_base), target :: transformation
    type(sll_logical_mesh_2d), pointer  :: m2d    
    sll_int32, intent(in) :: bc_left
    sll_int32, intent(in) :: bc_right
    sll_int32, intent(in) :: bc_bottom
    sll_int32, intent(in) :: bc_top
    sll_real64, dimension(:), optional :: point1_1d
    sll_real64, dimension(:), optional :: point2_1d
    sll_int32,optional :: sz_point1
    sll_int32,optional :: sz_point2
    sll_int32 :: i
    sll_int32 :: ierr   

    m2d => transformation%mesh
    field%T => transformation
    field%interp_2d => interpolator_2d
    !    field%mesh%written = .false.
    field%name      = trim(field_name)
    field%bc_left   = bc_left
    field%bc_right  = bc_right
    field%bc_bottom = bc_bottom
    field%bc_top    = bc_top

    ! Allocate internal array to store locally a copy of the data.
    SLL_ALLOCATE(field%values(m2d%num_cells1+1,m2d%num_cells2+1),ierr)    
    !print*,'first line',  field%values(1,:)
    !print*, 'second line', field%values(2,:)
!!$    call field%interp_2d%compute_interpolants( &
!!$         field%values, & !array_2d, &
!!$         point1_1d, &
!!$         sz_point1, &
!!$         point2_1d, &
!!$         sz_point2 )

  end subroutine initialize_scalar_field_2d_discrete_alt
  
  ! need to do something about deallocating the field proper, when allocated
  ! in the heap...
  subroutine delete_field_2d_discrete_alt( field )
    class(sll_scalar_field_2d_discrete_alt), intent(out) :: field
    sll_int32 :: ierr
    if( associated(field%values) ) then
       SLL_DEALLOCATE(field%values,ierr)
    end if
    nullify(field%T)
    nullify(field%interp_2d)
    nullify(field%point1_1d)
    nullify(field%point2_1d)
  end subroutine delete_field_2d_discrete_alt

  subroutine set_field_data_discrete_2d( field, values )
    class(sll_scalar_field_2d_discrete_alt), intent(inout) :: field
    sll_real64, dimension(:,:), intent(in) :: values
    if( (size(field%values,1) .ne. size(values,1) ) .or. &
        (size(field%values,2) .ne. size(values,2) ) ) then
        print *, 'WARNING, set_field_data_discrete_2d(), passed array ', &
             'is not of the size originally declared for this field.'
     end if
!!$    print *, 'size(field%values) = ', size(field%values,1), &
!!$         size(field%values,2), 'size(values) = ', size(values,1), size(values,2)
    field%values(:,:) = values(:,:)
  end subroutine set_field_data_discrete_2d

  subroutine update_interp_coeffs_2d_discrete( field )
    class(sll_scalar_field_2d_discrete_alt), intent(inout) :: field
    call field%interp_2d%compute_interpolants( field%values )
  end subroutine update_interp_coeffs_2d_discrete

  function get_transformation_discrete_alt( field ) result(res)
    class(sll_scalar_field_2d_discrete_alt), intent(in) :: field
    class(sll_coordinate_transformation_2d_base), pointer :: res
    res => field%T
  end function get_transformation_discrete_alt

  function get_logical_mesh_2d_discrete_alt( field ) result(res)
    class(sll_scalar_field_2d_discrete_alt), intent(in) :: field
    type(sll_logical_mesh_2d), pointer :: res
    res => field%T%mesh
  end function get_logical_mesh_2d_discrete_alt

  function get_jacobian_matrix_discrete_alt( field, eta1, eta2 ) result(res)
    class(sll_scalar_field_2d_discrete_alt), intent(in) :: field
    sll_real64, intent(in) :: eta1
    sll_real64, intent(in) :: eta2
    sll_real64, dimension(2,2) :: res
    res(:,:) = field%T%jacobian_matrix(eta1,eta2)
  end function get_jacobian_matrix_discrete_alt

  function value_at_pt_discrete( field, eta1, eta2 )
    class(sll_scalar_field_2d_discrete_alt), intent(in) :: field
    sll_real64, intent(in) :: eta1
    sll_real64, intent(in) :: eta2
    sll_real64             :: value_at_pt_discrete

    value_at_pt_discrete = field%interp_2d%interpolate_value(eta1,eta2)
  end function value_at_pt_discrete

  function value_at_index_discrete( field, i, j )
    class(sll_scalar_field_2d_discrete_alt), intent(in) :: field
    sll_int32, intent(in) :: i
    sll_int32, intent(in) :: j
    sll_real64            :: eta1
    sll_real64            :: eta2
    sll_real64            :: value_at_index_discrete
    eta1 = field%T%mesh%eta1_min + real(i-1,f64)*field%T%mesh%delta_eta1
    eta2 = field%T%mesh%eta2_min + real(j-1,f64)*field%T%mesh%delta_eta2
    value_at_index_discrete = field%interp_2d%interpolate_value(eta1,eta2)
  end function value_at_index_discrete

  function first_deriv_eta1_value_at_pt_discrete( field, eta1, eta2 )
    class(sll_scalar_field_2d_discrete_alt), intent(in) :: field
    sll_real64, intent(in) :: eta1
    sll_real64, intent(in) :: eta2
    sll_real64             :: first_deriv_eta1_value_at_pt_discrete
    
    first_deriv_eta1_value_at_pt_discrete = &
         field%interp_2d%interpolate_derivative_eta1(eta1,eta2)
  end function first_deriv_eta1_value_at_pt_discrete

  function first_deriv_eta2_value_at_pt_discrete( field, eta1, eta2 )
    class(sll_scalar_field_2d_discrete_alt), intent(in) :: field
    sll_real64, intent(in) :: eta1
    sll_real64, intent(in) :: eta2
    sll_real64             :: first_deriv_eta2_value_at_pt_discrete
    
    first_deriv_eta2_value_at_pt_discrete = &
         field%interp_2d%interpolate_derivative_eta2(eta1,eta2)
  end function first_deriv_eta2_value_at_pt_discrete
  
  function first_deriv_eta1_value_at_index_discrete( field, i, j )
    class(sll_scalar_field_2d_discrete_alt), intent(in) :: field
    sll_int32, intent(in) :: i
    sll_int32, intent(in) :: j
    sll_real64            :: eta1
    sll_real64            :: eta2
    sll_real64            :: first_deriv_eta1_value_at_index_discrete
    eta1 = field%T%mesh%eta1_min + real(i-1,f64)*field%T%mesh%delta_eta1
    eta2 = field%T%mesh%eta2_min + real(j-1,f64)*field%T%mesh%delta_eta2
    first_deriv_eta1_value_at_index_discrete = &
         field%interp_2d%interpolate_derivative_eta1(eta1,eta2)
  end function first_deriv_eta1_value_at_index_discrete

  function first_deriv_eta2_value_at_index_discrete( field, i, j )
    class(sll_scalar_field_2d_discrete_alt), intent(in) :: field
    sll_int32, intent(in) :: i
    sll_int32, intent(in) :: j
    sll_real64            :: eta1
    sll_real64            :: eta2
    sll_real64            :: first_deriv_eta2_value_at_index_discrete
    eta1 = field%T%mesh%eta1_min + real(i-1,f64)*field%T%mesh%delta_eta1
    eta2 = field%T%mesh%eta2_min + real(j-1,f64)*field%T%mesh%delta_eta2
    first_deriv_eta2_value_at_index_discrete = &
         field%interp_2d%interpolate_derivative_eta2(eta1,eta2)
  end function first_deriv_eta2_value_at_index_discrete

  subroutine write_to_file_discrete_2d( field, tag )
    class(sll_scalar_field_2d_discrete_alt), intent(in) :: field
    sll_int32, intent(in)                               :: tag
    sll_int32 :: nptsx1
    sll_int32 :: nptsx2
    sll_real64, dimension(:,:), allocatable :: x1coords
    sll_real64, dimension(:,:), allocatable :: x2coords
    sll_real64, dimension(:,:), allocatable :: values
    sll_real64                              :: eta1
    sll_real64                              :: eta2
    class(sll_coordinate_transformation_2d_base), pointer :: T
    type(sll_logical_mesh_2d), pointer      :: mesh
    sll_int32 :: i
    sll_int32 :: j
    sll_int32 :: ierr

    ! use the logical mesh information to find out the extent of the
    ! domain and allocate the arrays for the plotter.
    T      => field%get_transformation()
    mesh   => field%get_logical_mesh()
    nptsx1 = mesh%num_cells1 + 1
    nptsx2 = mesh%num_cells2 + 1
    SLL_ALLOCATE(x1coords(nptsx1,nptsx2),ierr)
    SLL_ALLOCATE(x2coords(nptsx1,nptsx2),ierr)
    SLL_ALLOCATE(values(nptsx1,nptsx2),ierr)

    ! Fill the arrays with the needed information.
    do j=1, nptsx2
       eta2 = mesh%eta2_min + (j-1)*mesh%delta_eta2 
       do i=1, nptsx1
          eta1 = mesh%eta1_min + (i-1)*mesh%delta_eta1
          x1coords(i,j) = field%T%x1(eta1,eta2)
          x2coords(i,j) = field%T%x2(eta1,eta2)
          values(i,j)   = field%value_at_point(eta1,eta2)
       end do
    end do

    call sll_gnuplot_curv_2d( &
         nptsx1, &
         nptsx2, &
         x1coords, &
         x2coords, &
         values, &
         trim(field%name), &
         tag, &
         ierr )

    SLL_DEALLOCATE_ARRAY(x1coords,ierr)
    SLL_DEALLOCATE_ARRAY(x2coords,ierr)
    SLL_DEALLOCATE_ARRAY(values,ierr)
  end subroutine write_to_file_discrete_2d




!!$  subroutine write_scalar_field_2d( &
!!$    scalar_field, &
!!$    multiply_by_jacobian, &
!!$    output_file_name, &
!!$    output_format)
!!$
!!$    class(scalar_field_2d) :: scalar_field
!!$    class(sll_mapped_mesh_2d_base), pointer :: mesh
!!$    logical, optional      :: multiply_by_jacobian 
!!$    sll_int32, optional    :: output_format 
!!$    character(len=*), optional    :: output_file_name 
!!$    sll_int32              :: local_format 
!!$
!!$    sll_int32  :: i1
!!$    sll_int32  :: i2
!!$    sll_real64 :: eta1
!!$    sll_real64 :: eta2
!!$    !sll_real64 :: avg
!!$    sll_int32  :: ierr
!!$    sll_real64, dimension(:,:), allocatable :: val
!!$    sll_int32  :: num_pts1
!!$    sll_int32  :: num_pts2
!!$    sll_int32  :: file_id
!!$    character(len=32) :: name
!!$    character(len=4) :: counter
!!$    character(len=4) :: center
!!$
!!$
!!$    if (.not. present(output_format)) then
!!$       local_format = SLL_IO_XDMF
!!$    else
!!$       local_format = output_format
!!$    end if
!!$
!!$    mesh => scalar_field%mesh
!!$
!!$    SLL_ASSERT(associated(mesh))  
!!$    if (.not. mesh%written) then
!!$       call mesh%write_to_file(local_format)
!!$    end if
!!$
!!$    num_pts1 = mesh%nc_eta1+1
!!$    num_pts2 = mesh%nc_eta2+1
!!$    if (scalar_field%data_position == NODE_CENTERED_FIELD) then
!!$       SLL_ALLOCATE(val(num_pts1,num_pts2), ierr)
!!$    else
!!$       SLL_ALLOCATE(val(num_pts1-1,num_pts2-1), ierr)
!!$    end if
!!$
!!$    if (.not.(present(multiply_by_jacobian))) then
!!$       val =  scalar_field%data
!!$    else !if (multiply_by_jacobian) then 
!!$
!!$       if (scalar_field%data_position == CELL_CENTERED_FIELD) then
!!$          eta2 =  0.5_f64 * mesh%delta_eta2
!!$          do i2 = 1, mesh%nc_eta2
!!$             eta1 = 0.5_f64 * mesh%delta_eta1
!!$             do i1 = 1, mesh%nc_eta1
!!$                val(i1,i2) = scalar_field%data(i1,i2) / mesh%jacobian(eta1, eta2)
!!$                eta1 = eta1 + mesh%delta_eta1
!!$             end do
!!$             eta2 = eta2 + mesh%delta_eta2
!!$          end do
!!$       else
!!$          eta2 =  0.0_f64 
!!$          do i2 = 1, num_pts2
!!$             eta1 = 0.0_f64 
!!$             do i1 = 1, num_pts1
!!$                val(i1,i2) = scalar_field%data(i1,i2)
!!$                eta1 = eta1 + mesh%delta_eta1
!!$             end do
!!$             eta2 = eta2 + mesh%delta_eta2
!!$          end do
!!$       end if
!!$     
!!$    end if
!!$
!!$
!!$    select case(local_format)
!!$    case (SLL_IO_XDMF)
!!$       
!!$       if (scalar_field%data_position == NODE_CENTERED_FIELD) then
!!$          center = "Node"
!!$       else if (scalar_field%data_position == CELL_CENTERED_FIELD) then
!!$          center = "Cell"
!!$       end if
!!$
!!$       if (.not. present(output_file_name)) then
!!$          scalar_field%plot_counter = scalar_field%plot_counter+1
!!$          call int2string(scalar_field%plot_counter, counter)
!!$          name = trim(scalar_field%name)//counter
!!$       else 
!!$          name = output_file_name
!!$       end if
!!$       call sll_xdmf_open(trim(name)//".xmf", &
!!$            scalar_field%mesh%label,        &
!!$            num_pts1,num_pts2,file_id,ierr)
!!$      
!!$       call sll_xdmf_write_array(trim(name), &
!!$                                 val,&
!!$                                 scalar_field%name,ierr,file_id, &
!!$                                 "Node")
!!$       call sll_xdmf_close(file_id,ierr)
!!$
!!$    case (SLL_IO_VTK)
!!$
!!$       call sll_ascii_file_create(trim(name)//".vtr", file_id, ierr)
!!$
!!$       write(file_id,"(a)")"<VTKFile type='RectilinearGrid'>"
!!$       write(file_id,"(a,6i5,a)")"<RectilinearGrid WholeExtent='",1, num_pts1,1,num_pts2,1,1,"'>"
!!$       write(file_id,"(a,6i5,a)")"<Piece Extent='",1, num_pts1,1,num_pts2,1,1,"'>"
!!$       write(file_id,"(a)")"<PointData>"
!!$       write(file_id,"(a)")"<DataArray type='Float64' Name='"//scalar_field%name//"' format='ascii'>"
!!$       write(file_id,"(a)")"</DataArray>"
!!$       write(file_id,"(a)")"</PointData>"
!!$       write(file_id,"(a)")"<Coordinates>"
!!$       write(file_id,"(a)")"<DataArray type='Float64' Name='"//scalar_field%name//"' format='ascii'>"
!!$       write(file_id,"(a)")"</DataArray>"
!!$       write(file_id,"(a)")"</Coordinates>"
!!$       write(file_id,"(a)")"</Piece>"
!!$       write(file_id,"(a)")"</RectilinearGrid>"
!!$       write(file_id,"(a)")"</VTKFile>"
!!$
!!$       close(file_id)
!!$
!!$    case (SLL_IO_GNUPLOT)
!!$       call sll_ascii_file_create(trim(name)//".vtr", file_id, ierr)
!!$       call sll_ascii_write_array_2d(file_id, val, ierr)
!!$       close(file_id)
!!$
!!$    case default
!!$
!!$       print*, "write_scalar_field_2d: requested output format not recognized."
!!$       stop
!!$    end select
!!$
!!$
!!$    SLL_DEALLOCATE_ARRAY(val,ierr)
!!$  end subroutine write_scalar_field_2d

end module sll_module_scalar_field_2d_alternative
