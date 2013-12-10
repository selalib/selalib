module sll_coordinate_transformation_2d_base_module
#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_file_io.h"
  use sll_logical_meshes
  implicit none

  ! A single abstract base class is defined which will further be extended
  ! by its subclasses. The two main types of coordinate transformations are 
  ! those represented by an analytic transformation and those represented by a
  ! discrete transformation. 

  ! The coordinate transformation always transforms from a cartesian space in
  ! the logical variables eta1, eta2, ... to a physical space of variables
  ! x1, x2, ... For example, the 2D case represents the transformation:
  !
  !              x1 = x1(eta1,eta2)
  !              x2 = x2(eta1,eta2)
  !
  ! The base class contains all the services (in the form of functions) that
  ! the different 'flavors' of coordinate transformations (analytic, discrete)
  ! should implement.

  type, abstract :: sll_coordinate_transformation_2d_base
     ! The decision to include the mesh here was determined by the need to,
     ! as error checking, test the association of this pointer within the
     ! functions that receive an argument of 
     ! class(sll_coordinate_transformation_2d_base)
     type(sll_logical_mesh_2d), pointer :: mesh
     !logical to remember when the mesh has already been written to file
     character(len=64) :: label
     logical           :: written! = .false.
   contains
     ! x1 = x1(eta1,eta2)
     procedure(geometry_function_ct), deferred, pass       :: x1
     ! x2 = x2(eta1,eta2)
     procedure(geometry_function_ct), deferred, pass       :: x2
     ! jacobian = jacobian(eta1,eta2)
     procedure(geometry_function_ct), deferred, pass       :: jacobian
     ! x1_at_node = x1_at_node(i,j)
     procedure(geometry_function_indices_ct), deferred, pass :: x1_at_node
     ! x2_at_node = x2_at_node(i,j)
     procedure(geometry_function_indices_ct), deferred, pass :: x2_at_node
     !jacobian_at_node = jacobian_at_node(i,j)
     procedure(geometry_function_indices_ct), deferred, pass :: jacobian_at_node
     ! jacobian_matrix = jacobian(matrix(eta1,eta2))
     procedure(matrix_geometry_function_ct), deferred, pass   :: jacobian_matrix
!     procedure(j_matrix_function_nopass), pointer, nopass :: jacobian_matrix
     procedure(matrix_geometry_function_ct), deferred, pass :: &
          inverse_jacobian_matrix
     ! x1_at_cell = x1_at_cell(i,j)
     procedure(geometry_function_indices_ct), deferred, pass :: x1_at_cell
     ! x1_at_cell = x1_at_cell(i,j)
     procedure(geometry_function_indices_ct), deferred, pass :: x2_at_cell
     ! jacobian_at_cell = jacobian_at_cell(i,j)
     procedure(geometry_function_indices_ct), deferred, pass :: jacobian_at_cell
     procedure(write_transformation_signature), deferred, pass :: write_to_file
     procedure(transformation_subroutine), deferred, pass      :: delete 
  end type sll_coordinate_transformation_2d_base

  
  !************************************************************************
  !
  !                       Function signatures
  !
  !************************************************************************

  !************************************************************************
  ! 2D CASE:
  !************************************************************************

   abstract interface
      function geometry_function_ct( transf, eta1, eta2 ) result(res)
        use sll_working_precision
        import sll_coordinate_transformation_2d_base
        class(sll_coordinate_transformation_2d_base) :: transf
        sll_real64, intent(in)   :: eta1
        sll_real64, intent(in)   :: eta2
        sll_real64               :: res
      end function geometry_function_ct
   end interface

   abstract interface
      function geometry_function_indices_ct( transf, i, j ) result(res)
        use sll_working_precision
        import sll_coordinate_transformation_2d_base       
        class(sll_coordinate_transformation_2d_base) :: transf
        sll_int32, intent(in)   :: i
        sll_int32, intent(in)   :: j
        sll_real64              :: res
      end function geometry_function_indices_ct
   end interface
   
   abstract interface
      function matrix_geometry_function_ct( transf, eta1, eta2 )
        use sll_working_precision       
        import sll_coordinate_transformation_2d_base
        class(sll_coordinate_transformation_2d_base) :: transf
        sll_real64, intent(in)         :: eta1
        sll_real64, intent(in)         :: eta2
        sll_real64                     :: matrix_geometry_function_ct(2,2)
      end function matrix_geometry_function_ct
   end interface
   
   ! It is mandatory to pass the params array since making it optional
   ! poses problems when storing the parameters array inside the object to
   ! pass it when the function is called.
   abstract interface
      function transformation_func_nopass( eta1, eta2, params ) result(res)
        use sll_working_precision
        sll_real64, intent(in) :: eta1
        sll_real64, intent(in) :: eta2
        sll_real64, dimension(:), intent(in) :: params
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
   ! GENERALLY USE AND DO NOT DEPEND ON A BASE OR DERIVED TYPE.

   ! The following interface is meant to specify the signature of the 
   ! (possibly user-defined) functions that should be passed as arguments
   ! to initialize the analytic transformations.
   abstract interface
      function two_arg_scalar_function( eta1, eta2 )
        use sll_working_precision
        sll_real64             :: two_arg_scalar_function
        sll_real64, intent(in) :: eta1
        sll_real64, intent(in) :: eta2
      end function two_arg_scalar_function
   end interface
   
   abstract interface
      function two_arg_message_passing_func( transf, eta1, eta2 )
        use sll_working_precision
        import     :: sll_coordinate_transformation_2d_base
        sll_real64                      :: two_arg_message_passing_func
        class(sll_coordinate_transformation_2d_base)  :: transf
        sll_real64, intent(in)          :: eta1
        sll_real64, intent(in)          :: eta2
      end function two_arg_message_passing_func
   end interface
   
   abstract interface
      subroutine write_to_file_signature( transf, label )
        import     :: sll_coordinate_transformation_2d_base
        class(sll_coordinate_transformation_2d_base)  :: transf
        character(len=*), optional      :: label
      end subroutine write_to_file_signature
   end interface

   abstract interface
      subroutine write_transformation_signature( transf, output_format )
       use sll_working_precision
        import :: sll_coordinate_transformation_2d_base
         class(sll_coordinate_transformation_2d_base)  :: transf
        sll_int32, optional :: output_format
      end subroutine write_transformation_signature
   end interface

   abstract interface
      subroutine transformation_subroutine( transf )
        import :: sll_coordinate_transformation_2d_base
        class(sll_coordinate_transformation_2d_base), intent(inout) :: transf
      end subroutine transformation_subroutine
   end interface


 end module sll_coordinate_transformation_2d_base_module
