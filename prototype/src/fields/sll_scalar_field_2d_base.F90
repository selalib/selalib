module sll_module_scalar_field_2d_base
#include "sll_working_precision.h"
  use sll_coordinate_transformation_2d_base_module
  use sll_logical_meshes
  implicit none


  ! Fundamental field type
  type, abstract :: sll_scalar_field_2d_base
     class(sll_coordinate_transformation_2d_base), pointer :: coord_trans 
   contains
     procedure(function_get_mesh), deferred, pass :: get_logical_mesh
     procedure(function_get_jacobian_matrix), deferred, pass :: &
          get_jacobian_matrix
     procedure(function_evaluation_real), deferred, pass :: value_at_point
     procedure(function_evaluation_integer), deferred, pass :: value_at_indices
     procedure(return_integer), deferred, pass :: interpolation_degree
!     procedure(file_output), deferred, pass :: write_to_file
     ! here we can continue with derivatives or whatever else that might
     ! be desired.
  end type sll_scalar_field_2d_base



  ! Function signatures
  abstract interface
     function function_get_mesh(field) result(res)
       use sll_logical_meshes
       import sll_scalar_field_2d_base
       class(sll_scalar_field_2d_base) :: field
       type(sll_logical_mesh_2d), pointer :: res
     end function function_get_mesh
  end interface

  abstract interface
     function function_get_jacobian_matrix(field ) result(res)
       use sll_working_precision
       import sll_scalar_field_2d_base
       class(sll_scalar_field_2d_base) :: field
       sll_real64, dimension(:,:), pointer :: res
     end function function_get_jacobian_matrix
  end interface


  abstract interface
     function function_evaluation_real( field, eta1, eta2 ) result(res)
       use sll_working_precision
       import sll_scalar_field_2d_base
       class(sll_scalar_field_2d_base) :: field
       sll_real64, intent(in) :: eta1
       sll_real64, intent(in) :: eta2
       sll_real64             :: res
     end function function_evaluation_real
  end interface

  abstract interface
     function function_evaluation_integer( field, i, j ) result(res)
       use sll_working_precision
       import sll_scalar_field_2d_base
       class(sll_scalar_field_2d_base) :: field
       sll_int32, intent(in)  :: i
       sll_int32, intent(in)  :: j
       sll_real64             :: res
     end function function_evaluation_integer
  end interface

  abstract interface
     function return_integer( field ) result(res)
       use sll_working_precision
       import sll_scalar_field_2d_base
       class(sll_scalar_field_2d_base) :: field
       sll_int32             :: res
     end function return_integer
  end interface

end module sll_module_scalar_field_2d_base
