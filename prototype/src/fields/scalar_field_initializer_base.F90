! This module aims at providing a single interface for initializing 
! scalar fields in a reasonably uniform way. As an example, in the 2D case, 
! the basic interface is the subroutine f_of_x1x2, which initializes a whole
! 2D array. The initializer object is a child of the scalar field initializer
! class and so it can have any specific data it requires. The f_of_x1x2 
! subroutine extracts all necessary information from a mesh object.

module sll_scalar_field_initializers_base
#include "sll_working_precision.h"
  use sll_module_mapped_meshes_2d_base
  implicit none

  type, abstract :: scalar_field_2d_initializer_base
     sll_int32   :: data_position
   contains
     procedure(scalar_field_initializer), deferred, pass :: f_of_x1x2
  end type scalar_field_2d_initializer_base

  abstract interface
     subroutine scalar_field_initializer( init_obj, mesh, data_out )
       use sll_working_precision
       import sll_mapped_mesh_2d_base, scalar_field_2d_initializer_base
       class(scalar_field_2d_initializer_base), intent(inout) :: init_obj
       class(sll_mapped_mesh_2d_base), intent(in)             :: mesh
       sll_real64, dimension(:,:), intent(out)                :: data_out
     end subroutine scalar_field_initializer
  end interface

  enum, bind(C)
     enumerator :: NODE_CENTERED_FIELD = 0, CELL_CENTERED_FIELD = 1
  end enum
  
end module sll_scalar_field_initializers_base
