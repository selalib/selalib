module sll_module_scalar_field_2d_base
#include "sll_working_precision.h"
  use sll_coordinate_transformation_2d_base_module
  use sll_logical_meshes
  implicit none


  !> Fundamental field type
  type, abstract :: sll_scalar_field_2d_base
     ! consider eliminating this transformation from this base class,
     ! it is already in the derived classes and is confusing...
     class(sll_coordinate_transformation_2d_base), pointer :: coord_trans 
   contains
     procedure(function_get_mesh), deferred, pass :: get_logical_mesh
     procedure(function_get_transformation), deferred, pass :: &
          get_transformation
     procedure(function_get_jacobian_matrix), deferred, pass :: &
          get_jacobian_matrix
     procedure(function_evaluation_real), deferred, pass :: value_at_point
     procedure(function_evaluation_integer), deferred, pass :: value_at_indices
     procedure(first_derivative_eta1_evaluation_real), deferred, pass :: &
          first_deriv_eta1_value_at_point
     procedure(first_derivative_eta2_evaluation_real), deferred, pass :: &
          first_deriv_eta2_value_at_point
     procedure(first_derivative_eta1_evaluation_integer), deferred, pass :: &
          first_deriv_eta1_value_at_indices
     procedure(first_derivative_eta2_evaluation_integer), deferred, pass :: &
          first_deriv_eta2_value_at_indices
     procedure(set_field_data_subroutine), deferred, pass :: set_field_data
     procedure(field_2d_message_pass), deferred, pass :: &
          update_interpolation_coefficients
     procedure(field_2d_file_output), deferred, pass :: write_to_file
     procedure(field_2d_subroutine), deferred, pass :: delete
     ! here we can continue with derivatives or whatever else that might
     ! be desired.
  end type sll_scalar_field_2d_base

  type sll_scalar_field_2d_base_ptr
     class(sll_scalar_field_2d_base), pointer :: base
  end type sll_scalar_field_2d_base_ptr


  !> Function signatures
  abstract interface
     function function_get_mesh(field) result(res)
       use sll_logical_meshes
       import sll_scalar_field_2d_base
       class(sll_scalar_field_2d_base), intent(in) :: field
       type(sll_logical_mesh_2d), pointer :: res
     end function function_get_mesh
  end interface

  abstract interface
     subroutine set_field_data_subroutine( field, values )
       use sll_working_precision
       import sll_scalar_field_2d_base
       class(sll_scalar_field_2d_base), intent(inout) :: field
       sll_real64, dimension(:,:), intent(in) :: values
     end subroutine set_field_data_subroutine
  end interface

  abstract interface
     subroutine field_2d_message_pass( field )
       import sll_scalar_field_2d_base
       class(sll_scalar_field_2d_base), intent(inout) :: field
     end subroutine field_2d_message_pass
  end interface


  abstract interface
     function function_get_transformation(field) result(res)
       use sll_coordinate_transformation_2d_base_module
       import sll_scalar_field_2d_base
       class(sll_scalar_field_2d_base), intent(in) :: field
       class(sll_coordinate_transformation_2d_base), pointer :: res
     end function function_get_transformation
  end interface

  abstract interface
     function function_get_jacobian_matrix(field,eta1,eta2 ) result(res)
       use sll_working_precision
       import sll_scalar_field_2d_base
       class(sll_scalar_field_2d_base), intent(in) :: field
       sll_real64, intent(in) :: eta1
       sll_real64, intent(in) :: eta2
       sll_real64, dimension(2,2) :: res
     end function function_get_jacobian_matrix
  end interface


  abstract interface
     function function_evaluation_real( field, eta1, eta2 ) result(res)
       use sll_working_precision
       import sll_scalar_field_2d_base
       class(sll_scalar_field_2d_base), intent(in) :: field
       sll_real64, intent(in) :: eta1
       sll_real64, intent(in) :: eta2
       sll_real64             :: res
     end function function_evaluation_real
  end interface

  abstract interface
     function function_evaluation_integer( field, i, j ) result(res)
       use sll_working_precision
       import sll_scalar_field_2d_base
       class(sll_scalar_field_2d_base), intent(in) :: field
       sll_int32, intent(in)  :: i
       sll_int32, intent(in)  :: j
       sll_real64             :: res
     end function function_evaluation_integer
  end interface

  abstract interface 
     function first_derivative_eta1_evaluation_real( field, eta1, eta2 ) result(res)
       use sll_working_precision
       import sll_scalar_field_2d_base
       class(sll_scalar_field_2d_base), intent(in) :: field
       sll_real64, intent(in) :: eta1
       sll_real64, intent(in) :: eta2
       sll_real64             :: res
     end function first_derivative_eta1_evaluation_real
  end interface

  abstract interface 
     function first_derivative_eta2_evaluation_real( field, eta1, eta2 ) result(res)
       use sll_working_precision
       import sll_scalar_field_2d_base
       class(sll_scalar_field_2d_base), intent(in) :: field
       sll_real64, intent(in) :: eta1
       sll_real64, intent(in) :: eta2
       sll_real64             :: res
     end function first_derivative_eta2_evaluation_real
  end interface

  abstract interface 
     function first_derivative_eta1_evaluation_integer( field, i, j ) result(res)
       use sll_working_precision
       import sll_scalar_field_2d_base
       class(sll_scalar_field_2d_base), intent(in) :: field
       sll_int32, intent(in)  :: i
       sll_int32, intent(in)  :: j
       sll_real64             :: res
     end function first_derivative_eta1_evaluation_integer
  end interface

  abstract interface 
     function first_derivative_eta2_evaluation_integer( field, i, j ) result(res)
       use sll_working_precision
       import sll_scalar_field_2d_base
       class(sll_scalar_field_2d_base), intent(in) :: field
       sll_int32, intent(in)  :: i
       sll_int32, intent(in)  :: j
       sll_real64             :: res
     end function first_derivative_eta2_evaluation_integer
  end interface

  abstract interface
     function return_integer( field ) result(res)
       use sll_working_precision
       import sll_scalar_field_2d_base
       class(sll_scalar_field_2d_base), intent(in) :: field
       sll_int32             :: res
     end function return_integer
  end interface

  abstract interface
     subroutine field_2d_file_output( field, tag )
       use sll_working_precision
       import sll_scalar_field_2d_base
       class(sll_scalar_field_2d_base), intent(in) :: field
       sll_int32, intent(in)                       :: tag
     end subroutine field_2d_file_output
  end interface

  abstract interface
     subroutine field_2d_subroutine( field )
       import sll_scalar_field_2d_base
       class(sll_scalar_field_2d_base), intent(inout) :: field
     end subroutine field_2d_subroutine
  end interface

end module sll_module_scalar_field_2d_base
