module sll_mapped_mesh_base
#include "sll_working_precision.h"
  implicit none
  
  ! A single abstract base class is defined which will further be extended
  ! by its subclasses. The two main types of mapped meshes are those
  ! represented by an analytic transformation and those represented by a
  ! discrete transformation. 
  type, abstract :: sll_mapped_mesh_2d_base
     sll_int32  :: num_pts1
     sll_int32  :: num_pts2
     sll_real64 :: delta1
     sll_real64 :: delta2
     sll_real64, dimension(:,:), pointer :: x1_cell
     sll_real64, dimension(:,:), pointer :: x2_cell
     sll_real64, dimension(:,:), pointer :: jacobians_n
     sll_real64, dimension(:,:), pointer :: jacobians_c
   contains
     procedure(geometry_function), deferred, pass       :: x1
     procedure(geometry_function), deferred, pass       :: x2
     procedure(geometry_function_nodes), deferred, pass :: x1_at_node
     procedure(geometry_function_nodes), deferred, pass :: x2_at_node
     procedure(geometry_function), deferred, pass       :: jacobian
!     procedure(j_matrix_function_nopass), pointer, nopass :: jacobian_matrix
  end type sll_mapped_mesh_2d_base

 
   !************************************************************************
   !                       Function signatures
   !
   !************************************************************************
   abstract interface
      function geometry_function( mesh, eta1, eta2 ) result(res)
        use sll_working_precision
        import sll_mapped_mesh_2d_base
        class(sll_mapped_mesh_2d_base) :: mesh
        sll_real64, intent(in)   :: eta1
        sll_real64, intent(in)   :: eta2
        sll_real64               :: res
      end function geometry_function
   end interface

   abstract interface
      function geometry_function_nodes( mesh, i, j ) result(res)
        use sll_working_precision
        import sll_mapped_mesh_2d_base       
        class(sll_mapped_mesh_2d_base) :: mesh
        sll_int32, intent(in)   :: i
        sll_int32, intent(in)   :: j
        sll_real64              :: res
      end function geometry_function_nodes
   end interface
   
   abstract interface
      function matrix_geometry_function( mesh, eta1, eta2 ) result(res)
        use sll_working_precision       
        import sll_mapped_mesh_2d_base
        class(sll_mapped_mesh_2d_base) :: mesh
        sll_real64, intent(in)   :: eta1
        sll_real64, intent(in)   :: eta2
        sll_real64               :: res(2,2)
      end function matrix_geometry_function
   end interface
   
   abstract interface
      function transformation_func_nopass( eta1, eta2 ) result(res)
        use sll_working_precision
        sll_real64, intent(in) :: eta1
        sll_real64, intent(in) :: eta2
        sll_real64             :: res
      end function transformation_func_nopass
   end interface
   
   abstract interface
      function matrix_geometry_function_nopass( eta1, eta2 ) result(res)
        use sll_working_precision 
        sll_real64, intent(in) :: eta1
        sll_real64, intent(in) :: eta2
        sll_real64             :: res(2,2)
      end function matrix_geometry_function_nopass
   end interface
   
   ! WE SHOULD PROBABLY HAVE A SINGLE FILE WITH ALL THE SIGNATURES THAT WE
   ! GENERALLY USE.
   abstract interface
      function two_arg_scalar_function( eta1, eta2 )
        use sll_working_precision
        sll_real64             :: two_arg_scalar_function
        sll_real64, intent(in) :: eta1
        sll_real64, intent(in) :: eta2
      end function two_arg_scalar_function
   end interface
   
   abstract interface
      function two_arg_message_passing_func( map, eta1, eta2 )
        use sll_working_precision
        import     :: sll_mapped_mesh_2d_base
        sll_real64                      :: two_arg_message_passing_func
        class(sll_mapped_mesh_2d_base)  :: map
        sll_real64, intent(in)          :: eta1
        sll_real64, intent(in)          :: eta2
      end function two_arg_message_passing_func
   end interface
   
 
   
 end module sll_mapped_mesh_base
