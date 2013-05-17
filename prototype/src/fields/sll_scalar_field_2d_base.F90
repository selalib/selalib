module sll_module_scalar_field_2d_base
#include "sll_working_precision.h"
  use sll_logical_meshes
  implicit none

  type, abstract :: sll_scalar_field_2d_base
     type(sll_logical_mesh_2d), pointer :: mesh 
   contains
     procedure(function_evaluation_real), deferred, pass :: value_at_point
     procedure(function_evaluation_integer), deferred, pass :: value_at_indices
!     procedure(file_output), deferred, pass :: write_to_file
     ! here we can continue with derivatives or whatever else that might
     ! be desired.
  end type sll_scalar_field_2d_base

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

end module sll_module_scalar_field_2d_base
