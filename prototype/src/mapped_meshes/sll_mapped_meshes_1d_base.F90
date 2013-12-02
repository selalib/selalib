module sll_module_mapped_meshes_1d_base
#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_file_io.h"

  implicit none

  ! A single abstract base class is defined which will further be extended
  ! by its subclasses. The two main types of mapped meshes are those
  ! represented by an analytic transformation and those represented by a
  ! discrete transformation. 

  type, abstract :: sll_mapped_mesh_1d_base
     sll_int32   :: nc_eta1
     sll_real64  :: delta_eta1
     character(len=64) :: label
     logical           :: written = .false.
   contains
     procedure(geometry_function_1d), deferred, pass       :: x1
     procedure(geometry_function_nodes_1d), deferred, pass :: x1_at_node
     procedure(geometry_function_nodes_1d), deferred, pass :: x1_at_cell
     procedure(geometry_function_1d), deferred, pass       :: jacobian
     procedure(geometry_function_nodes_1d), deferred, pass :: jacobian_at_node
     procedure(geometry_function_nodes_1d), deferred, pass :: jacobian_at_cell
     procedure, pass :: write_to_file => write_to_file_1d
  end type sll_mapped_mesh_1d_base
  
  !************************************************************************
  !
  !                       Function signatures
  !
  !************************************************************************

  !************************************************************************
  ! 1D CASE:
  !************************************************************************

  abstract interface
     function geometry_function_1d( mesh, eta1 ) result(res)
       use sll_working_precision
       import sll_mapped_mesh_1d_base
       class(sll_mapped_mesh_1d_base) :: mesh
       sll_real64, intent(in)    :: eta1
       sll_real64                :: res
     end function geometry_function_1d
  end interface

  abstract interface
      function geometry_function_nodes_1d( mesh, i ) result(res)
        use sll_working_precision
        import sll_mapped_mesh_1d_base       
        class(sll_mapped_mesh_1d_base) :: mesh
        sll_int32, intent(in)   :: i
        sll_real64              :: res
      end function geometry_function_nodes_1d
   end interface

   abstract interface
      function one_arg_scalar_function_nopass( eta1, params )
        use sll_working_precision
        sll_real64             :: one_arg_scalar_function_nopass
        sll_real64, intent(in) :: eta1
        sll_real64, dimension(:), intent(in) :: params
      end function one_arg_scalar_function_nopass
   end interface
   
   abstract interface
      function one_arg_message_passing_func( map, eta1 )
        use sll_working_precision
        import     :: sll_mapped_mesh_1d_base
        sll_real64                      :: one_arg_message_passing_func
        class(sll_mapped_mesh_1d_base)  :: map
        sll_real64, intent(in)          :: eta1
      end function one_arg_message_passing_func
   end interface

contains

  subroutine write_to_file_1d(mesh,output_format)
    class(sll_mapped_mesh_1d_base) :: mesh
    sll_int32, optional            :: output_format
    sll_int32                      :: error
    sll_real64, dimension(:), allocatable  :: x1_array
    sll_int32                      :: num_pts1  
    sll_int32                      :: i

    if(present(output_format))then
      print*,'There is just gnuplot format available'
    endif
    

    num_pts1 = mesh%nc_eta1+1
    SLL_ALLOCATE(x1_array(num_pts1),error)
    do i=1,num_pts1
       x1_array(i) = mesh%x1_at_node(i)
    end do
    call sll_gnuplot_write(x1_array,mesh%label,error)
    SLL_DEALLOCATE_ARRAY(x1_array,error)
  end subroutine

end module sll_module_mapped_meshes_1d_base
