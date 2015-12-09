!> @ingroup fields
!> @brief
!> Implements the field descriptor types
module sll_m_scalar_field_2d

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_errors.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_cartesian_meshes, only: &
    sll_cartesian_mesh_2d

  use sll_m_coordinate_transformation_2d_base, only: &
    sll_coordinate_transformation_2d_base

  use sll_m_gnuplot, only: &
    sll_gnuplot_2d

  use sll_m_interpolators_2d_base, only: &
    sll_c_interpolator_2d

  use sll_m_scalar_field_2d_base, only: &
    sll_scalar_field_2d_base

  use sll_m_utilities, only: &
    int2string, &
    sll_new_file_id

  use sll_m_xdmf, only: &
    sll_xdmf_curv2d_nodes

  implicit none

  public :: &
    new_scalar_field_2d_analytic, &
    new_scalar_field_2d_discrete, &
    sll_delete, &
    sll_scalar_field_2d_analytic, &
    sll_scalar_field_2d_discrete, &
    sll_scalar_field_2d_discrete_ptr, &
    two_var_parametrizable_function

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

type, extends(sll_scalar_field_2d_base) :: sll_scalar_field_2d_analytic

  private
  type(sll_cartesian_mesh_2d),                pointer         :: mesh
  procedure(two_var_parametrizable_function), pointer, nopass :: func
  procedure(two_var_parametrizable_function), pointer, nopass :: first_deriv_eta1
  procedure(two_var_parametrizable_function), pointer, nopass :: first_deriv_eta2

  sll_real64, dimension(:), pointer                     :: params
  character(len=64)                                     :: name
  class(sll_coordinate_transformation_2d_base), pointer :: t

  sll_int32 :: bc1_min
  sll_int32 :: bc1_max
  sll_int32 :: bc2_min
  sll_int32 :: bc2_max
  logical   :: present_deriv_eta1_int
  logical   :: present_deriv_eta2_int

contains

  procedure, pass(field) :: initialize => initialize_scalar_field_2d_analytic
  procedure, pass(field) :: get_transformation => get_transformation_analytic
  procedure, pass(field) :: get_cartesian_mesh => get_cartesian_mesh_2d_analytic
  procedure, pass(field) :: get_jacobian_matrix => get_jacobian_matrix_analytic
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
  procedure, pass(field) :: delete => delete_field_2d_analytic

end type sll_scalar_field_2d_analytic


type, extends(sll_scalar_field_2d_base) :: sll_scalar_field_2d_discrete

  type(sll_cartesian_mesh_2d), pointer :: mesh
  sll_real64, dimension(:,:), pointer  :: values => null()
  logical                              :: owns_memory = .true. 
  character(len=64)                    :: name

  class(sll_coordinate_transformation_2d_base), pointer :: T
  class(sll_c_interpolator_2d), pointer              :: interp_2d

  sll_real64, dimension(:), pointer :: point1_1d
  sll_real64, dimension(:), pointer :: point2_1d

  sll_int32 :: bc1_min
  sll_int32 :: bc1_max
  sll_int32 :: bc2_min
  sll_int32 :: bc2_max

contains

  procedure, pass(field) :: initialize => initialize_scalar_field_2d_discrete
  procedure, pass(field) :: update_interpolation_coefficients => &
                              update_interp_coeffs_2d_discrete
  procedure, pass(field) :: get_transformation => get_transformation_discrete
  procedure, pass(field) :: get_cartesian_mesh => get_cartesian_mesh_2d_discrete
  procedure, pass(field) :: get_jacobian_matrix => get_jacobian_matrix_discrete
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
  procedure, pass(field) :: free_internal_data_copy => free_data_discrete_2d
  procedure, pass(field) :: reset_data_pointer => reset_ptr_discrete_2d
  procedure, pass(field) :: get_data_pointer => get_data_ptr_discrete_2d
  procedure, pass(field) :: write_to_file => write_to_file_discrete_2d
  procedure, pass(field) :: delete => delete_field_2d_discrete

end type sll_scalar_field_2d_discrete

type sll_scalar_field_2d_discrete_ptr
  type(sll_scalar_field_2d_discrete), pointer :: f
end type sll_scalar_field_2d_discrete_ptr

abstract interface

  function two_var_parametrizable_function( eta1, eta2, params )
    use sll_m_working_precision
    sll_real64 :: two_var_parametrizable_function
    sll_real64, intent(in) :: eta1
    sll_real64, intent(in) :: eta2
    sll_real64, dimension(:), intent(in) :: params
  end function two_var_parametrizable_function

end interface

abstract interface

  function scalar_function_2d( eta1, eta2 )
    use sll_m_working_precision
    sll_real64 :: scalar_function_2d
    sll_real64, intent(in)  :: eta1
    sll_real64, intent(in)  :: eta2
  end function scalar_function_2d

end interface

interface sll_delete
  module procedure delete_field_2d_analytic, delete_field_2d_discrete
end interface sll_delete

contains

! **************************************************************************
!
!                         ANALYTIC CASE
!
! **************************************************************************

function value_at_pt_analytic( field, eta1, eta2 )

class(sll_scalar_field_2d_analytic), intent(in) :: field
sll_real64, intent(in) :: eta1
sll_real64, intent(in) :: eta2
sll_real64             :: value_at_pt_analytic

value_at_pt_analytic = field%func(eta1,eta2,field%params)

end function value_at_pt_analytic

function value_at_index_analytic( field, i, j )

class(sll_scalar_field_2d_analytic), intent(in) :: field
sll_int32, intent(in)                           :: i
sll_int32, intent(in)                           :: j
type(sll_cartesian_mesh_2d), pointer            :: lm
sll_real64                                      :: eta1
sll_real64                                      :: eta2
sll_real64                                      :: value_at_index_analytic

lm => field%T%mesh
eta1 = lm%eta1_min + (i-1)*lm%delta_eta1
eta2 = lm%eta2_min + (j-1)*lm%delta_eta2
value_at_index_analytic = field%func(eta1,eta2,field%params)

end function value_at_index_analytic

function first_deriv_eta1_value_at_pt_analytic( field, eta1, eta2 )

class(sll_scalar_field_2d_analytic), intent(in) :: field
sll_real64                         , intent(in) :: eta1
sll_real64                         , intent(in) :: eta2

sll_real64 :: first_deriv_eta1_value_at_pt_analytic

character(len=128)          :: err_msg
character(len=*), parameter :: this_fun_name = &
  'first_deriv_eta1_value_at_pt_analytic'

if ( field%present_deriv_eta1_int ) then 
  first_deriv_eta1_value_at_pt_analytic = &
       field%first_deriv_eta1(eta1,eta2,field%params)
else 
  first_deriv_eta1_value_at_pt_analytic = 0.0_f64
  err_msg = "In "// field%name // &
            ": first derivative in eta1 not given in the initialization."
  SLL_ERROR( this_fun_name, err_msg )
end if

end function first_deriv_eta1_value_at_pt_analytic

function first_deriv_eta2_value_at_pt_analytic( field, eta1, eta2 )

class(sll_scalar_field_2d_analytic), intent(in) :: field
sll_real64                         , intent(in) :: eta1
sll_real64                         , intent(in) :: eta2

sll_real64 :: first_deriv_eta2_value_at_pt_analytic 

character(len=128)          :: err_msg
character(len=*), parameter :: this_fun_name = &
  'first_deriv_eta2_value_at_pt_analytic' 

if ( field%present_deriv_eta2_int ) then 
  first_deriv_eta2_value_at_pt_analytic = &
       field%first_deriv_eta2(eta1,eta2,field%params)
else 
  first_deriv_eta2_value_at_pt_analytic  = 0.0_f64
  err_msg = "In "// field%name // &
            ": first derivative in eta2 &
            & not given in the initialization."
  SLL_ERROR( this_fun_name, err_msg )
end if
    
end function first_deriv_eta2_value_at_pt_analytic
  
function first_deriv_eta1_value_at_index_analytic( field, i, j)

class(sll_scalar_field_2d_analytic), intent(in) :: field
sll_int32, intent(in)                           :: i
sll_int32, intent(in)                           :: j
sll_real64                                      :: eta1
sll_real64                                      :: eta2
type(sll_cartesian_mesh_2d), pointer            :: lm

sll_real64 :: first_deriv_eta1_value_at_index_analytic

lm => field%T%mesh
eta1 = lm%eta1_min + real(i-1,f64)*lm%delta_eta1
eta2 = lm%eta2_min + real(j-1,f64)*lm%delta_eta2

if ( field%present_deriv_eta1_int ) then 
   first_deriv_eta1_value_at_index_analytic = &
        field%first_deriv_eta1(eta1,eta2,field%params)
else 
   first_deriv_eta1_value_at_index_analytic = 0.0_f64
   print*,field%name, &
        'first_deriv_eta1_value_at_index_analytic(): ERROR, ', &
        'first derivative in eta1 is not given in the initialization'
end if
  
end function first_deriv_eta1_value_at_index_analytic
  
function first_deriv_eta2_value_at_index_analytic( field, i, j)

class(sll_scalar_field_2d_analytic), intent(in) :: field
sll_int32, intent(in)                           :: i
sll_int32, intent(in)                           :: j
sll_real64                                      :: eta1
sll_real64                                      :: eta2
type(sll_cartesian_mesh_2d), pointer            :: lm

sll_real64 :: first_deriv_eta2_value_at_index_analytic

lm => field%T%mesh
eta1 = lm%eta1_min + real(i-1,f64)*lm%delta_eta1
eta2 = lm%eta2_min + real(j-1,f64)*lm%delta_eta2
  
if ( field%present_deriv_eta2_int ) then 
  first_deriv_eta2_value_at_index_analytic = &
        field%first_deriv_eta2(eta1,eta2,field%params)
else 
   print*,' first derivative in eta2 is not given in the initialization'
end if
  
end function first_deriv_eta2_value_at_index_analytic

function new_scalar_field_2d_analytic( func,             &
                                       field_name,       &
                                       transformation,   &
                                       bc1_min,          &
                                       bc1_max,          &
                                       bc2_min,          &
                                       bc2_max,          &
                                       func_params,      &
                                       first_deriv_eta1, &
                                       first_deriv_eta2) result(obj)
  
type(sll_scalar_field_2d_analytic), pointer          :: obj
procedure(two_var_parametrizable_function)           :: func
character(len=*), intent(in)                         :: field_name
class(sll_coordinate_transformation_2d_base), target :: transformation
sll_int32, intent(in)                                :: bc1_min
sll_int32, intent(in)                                :: bc1_max
sll_int32, intent(in)                                :: bc2_min
sll_int32, intent(in)                                :: bc2_max
sll_real64, dimension(:), intent(in)                 :: func_params
procedure(two_var_parametrizable_function), optional :: first_deriv_eta1
procedure(two_var_parametrizable_function), optional :: first_deriv_eta2
sll_int32                                            :: ierr

SLL_ALLOCATE(obj,ierr)

call obj%initialize( func,               &
&                    field_name,         &
&                    transformation,     &
&                    bc1_min,            &
&                    bc1_max,            &
&                    bc2_min,            &
&                    bc2_max,            &
&                    func_params,        &
&                    first_deriv_eta1,   &
&                    first_deriv_eta2)

end function new_scalar_field_2d_analytic

subroutine set_field_data_analytic_2d( field, values )

class(sll_scalar_field_2d_analytic), intent(inout) :: field
sll_real64, dimension(:,:), intent(in) :: values

print *, 'WARNING: set_field_data_analytic_2d(): it is useless to ', &
       'call this function on an analytic scalar field.'
SLL_ASSERT(associated(field%mesh))
SLL_ASSERT(size(values,1)>0)

end subroutine set_field_data_analytic_2d

subroutine update_interpolation_coefficients_2d_analytic( field )

class(sll_scalar_field_2d_analytic), intent(inout) :: field
print *, 'WARNING: update_interpolation_coefficients_2d_analytic(): ', &
     ' it is useless to call this function on an analytic scalar field.'
SLL_ASSERT(associated(field%mesh))

end subroutine update_interpolation_coefficients_2d_analytic

subroutine delete_field_2d_analytic( field )

class(sll_scalar_field_2d_analytic), intent(inout) :: field
! nothing internal do deallocate, just nullify pointers. Can't call
! delete on them because the field does not 'own' these data.
if(associated(field%func))  nullify(field%func)
if(associated(field%params))nullify(field%params)
if(associated(field%T))     nullify(field%T)

end subroutine delete_field_2d_analytic

! For those cases in which handling pointers to field structures is not
! convenient, we offer the following initialization.
subroutine initialize_scalar_field_2d_analytic( field,            &
&                                               func,             &
&                                               field_name,       &
&                                               transformation,   &
&                                               bc1_min,          &
&                                               bc1_max,          &
&                                               bc2_min,          &
&                                               bc2_max,          &
&                                               func_params,      &
&                                               first_deriv_eta1, &
&                                               first_deriv_eta2)

class(sll_scalar_field_2d_analytic), intent(out)     :: field
procedure(two_var_parametrizable_function)           :: func
character(len=*), intent(in)                         :: field_name
class(sll_coordinate_transformation_2d_base), target :: transformation
sll_int32, intent(in)                                :: bc1_min
sll_int32, intent(in)                                :: bc1_max
sll_int32, intent(in)                                :: bc2_min
sll_int32, intent(in)                                :: bc2_max
sll_real64, dimension(:), intent(in)                 :: func_params
procedure(two_var_parametrizable_function), optional :: first_deriv_eta1
procedure(two_var_parametrizable_function), optional :: first_deriv_eta2
sll_int32 :: ierr

SLL_ALLOCATE(field%params(size(func_params)),ierr)

field%params(:) = func_params
field%T         => transformation
field%func      => func
field%name      = trim(field_name)
field%bc1_min   = bc1_min
field%bc1_max   = bc1_max
field%bc2_min   = bc2_min
field%bc2_max   = bc2_max
  
if (present(first_deriv_eta1)) then
   field%first_deriv_eta1 => first_deriv_eta1
   field%present_deriv_eta1_int = .TRUE.
end if
  if (present(first_deriv_eta2)) then
     field%first_deriv_eta2 => first_deriv_eta2
     field%present_deriv_eta2_int = .TRUE.
end if

end subroutine initialize_scalar_field_2d_analytic

  
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

function get_transformation_analytic( field ) result(res)

class(sll_scalar_field_2d_analytic), intent(in) :: field
class(sll_coordinate_transformation_2d_base), pointer :: res
res => field%T

end function get_transformation_analytic

function get_cartesian_mesh_2d_analytic( field ) result(res)

class(sll_scalar_field_2d_analytic), intent(in) :: field
class(sll_cartesian_mesh_2d), pointer :: res
res => field%T%get_cartesian_mesh()

end function get_cartesian_mesh_2d_analytic

function get_jacobian_matrix_analytic( field, eta1, eta2 ) result(res)

class(sll_scalar_field_2d_analytic), intent(in) :: field
sll_real64, intent(in) :: eta1
sll_real64, intent(in) :: eta2
sll_real64, dimension(2,2) :: res
res = (field%T%jacobian_matrix(eta1,eta2))

end function get_jacobian_matrix_analytic

subroutine write_to_file_analytic_2d( field, tag )

class(sll_scalar_field_2d_analytic), intent(in) :: field
sll_int32, intent(in)                           :: tag
sll_int32 :: nptsx1
sll_int32 :: nptsx2
sll_real64, dimension(:,:), allocatable :: x1coords
sll_real64, dimension(:,:), allocatable :: x2coords
sll_real64, dimension(:,:), allocatable :: values
sll_real64                              :: eta1
sll_real64                              :: eta2
class(sll_coordinate_transformation_2d_base), pointer :: T
class(sll_cartesian_mesh_2d), pointer      :: mesh
sll_int32 :: i
sll_int32 :: j
sll_int32 :: ierr

! use the logical mesh information to find out the extent of the
! domain and allocate the arrays for the plotter.
T      => field%get_transformation()
mesh   => field%get_cartesian_mesh()
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

call sll_gnuplot_2d( nptsx1,           &
&                    nptsx2,           &
&                    x1coords,         &
&                    x2coords,         &
&                    values,           &
&                    trim(field%name), &
&                    tag,              &
&                    ierr )

SLL_DEALLOCATE_ARRAY(x1coords,ierr)
SLL_DEALLOCATE_ARRAY(x2coords,ierr)
SLL_DEALLOCATE_ARRAY(values,ierr)

end subroutine write_to_file_analytic_2d

! **************************************************************************
!
!                         DISCRETE CASE
!
! **************************************************************************

function new_scalar_field_2d_discrete( field_name,      &
                                       interpolator_2d, &
                                       transformation,  &
                                       bc1_min,         &
                                       bc1_max,         &
                                       bc2_min,         &
                                       bc2_max,         &
                                       point1_1d,       &
                                       sz_point1,       &
                                       point2_1d,       &
                                       sz_point2) result(obj)

type(sll_scalar_field_2d_discrete), pointer  :: obj
character(len=*), intent(in)                 :: field_name
class(sll_c_interpolator_2d), target      :: interpolator_2d
class(sll_coordinate_transformation_2d_base) :: transformation
sll_int32, intent(in)                        :: bc1_min
sll_int32, intent(in)                        :: bc1_max
sll_int32, intent(in)                        :: bc2_min
sll_int32, intent(in)                        :: bc2_max
sll_real64, dimension(:), optional           :: point1_1d
sll_real64, dimension(:), optional           :: point2_1d
sll_int32, optional                          :: sz_point1
sll_int32, optional                          :: sz_point2
sll_int32                                    :: ierr
class(sll_cartesian_mesh_2d), pointer        :: mesh

mesh => transformation%get_cartesian_mesh()

SLL_ALLOCATE(obj,ierr)

call obj%initialize( field_name,      &
                     interpolator_2d, &
                     transformation,  &
                     mesh%num_cells1, &
                     mesh%num_cells2, &
                     bc1_min,         &
                     bc1_max,         &
                     bc2_min,         &
                     bc2_max,         &
                     point1_1d,       &
                     sz_point1,       &
                     point2_1d,       &
                     sz_point2)

end function new_scalar_field_2d_discrete
  
subroutine initialize_scalar_field_2d_discrete( field,           &
&                                               field_name,      &
&                                               interpolator_2d, &
&                                               transformation,  &
&                                               num_cells1,      &
&                                               num_cells2,      &
&                                               bc1_min,         &
&                                               bc1_max,         &
&                                               bc2_min,         &
&                                               bc2_max,         &
&                                               point1_1d,       &
&                                               sz_point1,       &
&                                               point2_1d,       &
&                                               sz_point2)
    
class(sll_scalar_field_2d_discrete),          intent(inout) :: field
character(len=*),                             intent(in)    :: field_name
class(sll_c_interpolator_2d),              target        :: interpolator_2d
class(sll_coordinate_transformation_2d_base), target        :: transformation
sll_int32,                                    intent(in)    :: bc1_min
sll_int32,                                    intent(in)    :: bc1_max
sll_int32,                                    intent(in)    :: bc2_min
sll_int32,                                    intent(in)    :: bc2_max
sll_int32,                                    intent(in)    :: num_cells1
sll_int32,                                    intent(in)    :: num_cells2
sll_real64, dimension(:),                     optional      :: point1_1d
sll_real64, dimension(:),                     optional      :: point2_1d
sll_int32,                                    optional      :: sz_point1
sll_int32,                                    optional      :: sz_point2
sll_int32 :: ierr   

field%T => transformation

field%interp_2d => interpolator_2d

field%name    = trim(field_name)
field%bc1_min = bc1_min
field%bc1_max = bc1_max
field%bc2_min = bc2_min
field%bc2_max = bc2_max

! Allocate internal array to store locally a copy of the data.
SLL_ALLOCATE(field%values(num_cells1+1,num_cells2+1),ierr)

return
SLL_ASSERT(present(point1_1d))
SLL_ASSERT(present(point2_1d))
SLL_ASSERT(present(sz_point1))
SLL_ASSERT(present(sz_point2))
end subroutine initialize_scalar_field_2d_discrete
  

! need to do something about deallocating the field proper, when allocated
! in the heap...
subroutine delete_field_2d_discrete( field )
class(sll_scalar_field_2d_discrete), intent(inout) :: field
sll_int32 :: ierr
if(field%owns_memory) then
   if(associated(field%values)) SLL_DEALLOCATE(field%values,ierr)
end if
if(associated(field%T))         nullify(field%T)
if(associated(field%interp_2d)) nullify(field%interp_2d)
if(associated(field%point1_1d)) nullify(field%point1_1d)
if(associated(field%point2_1d)) nullify(field%point2_1d)
end subroutine delete_field_2d_discrete

subroutine set_field_data_discrete_2d( field, values )
class(sll_scalar_field_2d_discrete), intent(inout) :: field
sll_real64, dimension(:,:), intent(in) :: values
class(sll_cartesian_mesh_2d), pointer :: m

m => field%get_cartesian_mesh()
if( (size(values,1) < m%num_cells1 ) .or. &
    (size(values,2) < m%num_cells2 ) ) then
   print *, 'WARNING, set_field_data_discrete_2d(), passed array ', &
        'is smaller than the size of data originally declared for ', &
        'this field. Size of values in first dimension:', size(values,1),&
        ' Size of mesh: ', m%num_cells1, ' Size of values in second ', &
        'dimension:', size(values,2), 'Size of mesh: ', m%num_cells2
end if
field%values(:,:) = values(:,:)
end subroutine set_field_data_discrete_2d

! There is a bit of background history to the existence of the 
! free_data_discrete_2d() and reset_ptr_discrete_2d() routines. By request
! of a user, the default behavior of the field is to manage its own copy
! of the data. Thus upon creation, the object allocates the necessary
! memory to store nodal values. However, it may be desired that the object
! DOES NOT manage its own memory, but rather that it only points to some
! external block of memory. To permit the change of the behavior of the 
! field from its default, is the role of these functions. The first
! deallocates the memory and sets the flag which indicates that the field
! no longer owns its memory. The second permits to reset the data pointer
! to whatever is desired.
subroutine free_data_discrete_2d( field )
class(sll_scalar_field_2d_discrete), intent(inout) :: field
sll_int32 :: ierr

if( .not. associated(field%values) ) then
   print *, 'ERROR, free_data_discrete_2d(): the internal copy of the ', &
        'data has been already freed or never allocated.'
end if
SLL_DEALLOCATE(field%values,ierr)
field%owns_memory = .false.
end subroutine free_data_discrete_2d

subroutine reset_ptr_discrete_2d( field, values )
class(sll_scalar_field_2d_discrete), intent(inout) :: field
sll_real64, dimension(:,:), target :: values
if( field%owns_memory .eqv. .true. ) then
   print *, 'ERROR, reset_ptr_discrete_2d(): the data pointer can not ', &
        'be reset without a previous call to free_internal_data_copy().',&
        'This object is not being used properly. A memory leak has ', &
        'occurred. Continue at your peril.'
end if
field%values => values
end subroutine reset_ptr_discrete_2d

function get_data_ptr_discrete_2d( field ) result(ptr)

  sll_real64, dimension(:,:), pointer :: ptr
  class(sll_scalar_field_2d_discrete), intent(inout) :: field
  SLL_ASSERT(associated(field%values))
  
  ptr => field%values

end function get_data_ptr_discrete_2d

subroutine update_interp_coeffs_2d_discrete( field )
class(sll_scalar_field_2d_discrete), intent(inout) :: field
call field%interp_2d%compute_interpolants( field%values )
end subroutine update_interp_coeffs_2d_discrete

function get_transformation_discrete( field ) result(res)

  class(sll_scalar_field_2d_discrete), intent(in) :: field
  class(sll_coordinate_transformation_2d_base), pointer :: res
  
  res => field%T

end function get_transformation_discrete

function get_cartesian_mesh_2d_discrete( field ) result(res)

  class(sll_scalar_field_2d_discrete), intent(in) :: field
  class(sll_cartesian_mesh_2d), pointer :: res
  class(sll_coordinate_transformation_2d_base),pointer :: transf
  
  transf => field%T
  res => transf%get_cartesian_mesh()

end function get_cartesian_mesh_2d_discrete

function get_jacobian_matrix_discrete( field, eta1, eta2 ) result(res)

  class(sll_scalar_field_2d_discrete), intent(in) :: field
  sll_real64, intent(in) :: eta1
  sll_real64, intent(in) :: eta2
  sll_real64, dimension(2,2) :: res

  res(:,:) = field%T%jacobian_matrix(eta1,eta2)

end function get_jacobian_matrix_discrete

function value_at_pt_discrete( field, eta1, eta2 )

  class(sll_scalar_field_2d_discrete), intent(in) :: field
  sll_real64, intent(in) :: eta1
  sll_real64, intent(in) :: eta2
  sll_real64             :: value_at_pt_discrete

  value_at_pt_discrete = field%interp_2d%interpolate_from_interpolant_value(eta1,eta2)

end function value_at_pt_discrete

function value_at_index_discrete( field, i, j )
  class(sll_scalar_field_2d_discrete), intent(in) :: field
  sll_int32, intent(in) :: i
  sll_int32, intent(in) :: j
  sll_real64            :: eta1
  sll_real64            :: eta2
  class(sll_cartesian_mesh_2d), pointer :: lm
  sll_real64            :: value_at_index_discrete

  lm => field%get_cartesian_mesh()
  eta1 = lm%eta1_min + real(i-1,f64)*lm%delta_eta1
  eta2 = lm%eta2_min + real(j-1,f64)*lm%delta_eta2
  value_at_index_discrete = field%interp_2d%interpolate_from_interpolant_value(eta1,eta2)
end function value_at_index_discrete

function first_deriv_eta1_value_at_pt_discrete( field, eta1, eta2 )
  class(sll_scalar_field_2d_discrete), intent(in) :: field
  sll_real64, intent(in) :: eta1
  sll_real64, intent(in) :: eta2
  sll_real64             :: first_deriv_eta1_value_at_pt_discrete
    
  first_deriv_eta1_value_at_pt_discrete = &
       field%interp_2d%interpolate_from_interpolant_derivative_eta1(eta1,eta2)
end function first_deriv_eta1_value_at_pt_discrete

function first_deriv_eta2_value_at_pt_discrete( field, eta1, eta2 )
  class(sll_scalar_field_2d_discrete), intent(in) :: field
  sll_real64, intent(in) :: eta1
  sll_real64, intent(in) :: eta2
  sll_real64             :: first_deriv_eta2_value_at_pt_discrete
    
  first_deriv_eta2_value_at_pt_discrete = &
       field%interp_2d%interpolate_from_interpolant_derivative_eta2(eta1,eta2)
end function first_deriv_eta2_value_at_pt_discrete
  
function first_deriv_eta1_value_at_index_discrete( field, i, j )
  class(sll_scalar_field_2d_discrete), intent(in) :: field
  sll_int32, intent(in) :: i
  sll_int32, intent(in) :: j
  sll_real64            :: eta1
  sll_real64            :: eta2
  class(sll_cartesian_mesh_2d), pointer :: lm
  sll_real64            :: first_deriv_eta1_value_at_index_discrete

  lm => field%get_cartesian_mesh()
  eta1 = lm%eta1_min + real(i-1,f64)*lm%delta_eta1
  eta2 = lm%eta2_min + real(j-1,f64)*lm%delta_eta2
  first_deriv_eta1_value_at_index_discrete = &
        field%interp_2d%interpolate_from_interpolant_derivative_eta1(eta1,eta2)

end function first_deriv_eta1_value_at_index_discrete

function first_deriv_eta2_value_at_index_discrete( field, i, j )
  class(sll_scalar_field_2d_discrete), intent(in) :: field
  sll_int32, intent(in) :: i
  sll_int32, intent(in) :: j
  sll_real64            :: eta1
  sll_real64            :: eta2
  class(sll_cartesian_mesh_2d), pointer :: lm
  sll_real64            :: first_deriv_eta2_value_at_index_discrete

  lm => field%get_cartesian_mesh()
  eta1 = lm%eta1_min + real(i-1,f64)*lm%delta_eta1
  eta2 = lm%eta2_min + real(j-1,f64)*lm%delta_eta2
  first_deriv_eta2_value_at_index_discrete = &
       field%interp_2d%interpolate_from_interpolant_derivative_eta2(eta1,eta2)
end function first_deriv_eta2_value_at_index_discrete

subroutine write_to_file_discrete_2d( field, tag )

class(sll_scalar_field_2d_discrete), intent(in) :: field
sll_int32, intent(in)                           :: tag
sll_int32 :: nptsx1
sll_int32 :: nptsx2
sll_real64, dimension(:,:), allocatable :: x1coords
sll_real64, dimension(:,:), allocatable :: x2coords
sll_real64, dimension(:,:), allocatable :: values
sll_real64                              :: eta1
sll_real64                              :: eta2
class(sll_coordinate_transformation_2d_base), pointer :: T
class(sll_cartesian_mesh_2d), pointer      :: mesh
sll_int32 :: i
sll_int32 :: j
sll_int32 :: ierr
character(len=4) :: ctag

! use the logical mesh information to find out the extent of the
! domain and allocate the arrays for the plotter.
T      => field%get_transformation()
mesh   => field%get_cartesian_mesh()
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

call sll_gnuplot_2d( nptsx1,           &
&                    nptsx2,           &
&                    x1coords,         &
&                    x2coords,         &
&                    values,           &
&                    trim(field%name), &
&                    tag,              &
&                    ierr )

call int2string(tag, ctag)
call sll_xdmf_curv2d_nodes(trim(field%name)//ctag, &
     values, "values", x1coords, x2coords, "HDF5") 

SLL_DEALLOCATE_ARRAY(x1coords,ierr)
SLL_DEALLOCATE_ARRAY(x2coords,ierr)
SLL_DEALLOCATE_ARRAY(values,ierr)

end subroutine write_to_file_discrete_2d

end module sll_m_scalar_field_2d
