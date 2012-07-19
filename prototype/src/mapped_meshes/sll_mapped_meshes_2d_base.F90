module sll_module_mapped_meshes_2d_base
#include "sll_memory.h"
#include "sll_working_precision.h"

use sll_io

  implicit none
  
  ! A single abstract base class is defined which will further be extended
  ! by its subclasses. The two main types of mapped meshes are those
  ! represented by an analytic transformation and those represented by a
  ! discrete transformation. 

  type, abstract :: sll_mapped_mesh_2d_base
     sll_int32  :: nc_eta1
     sll_int32  :: nc_eta2
     sll_real64 :: delta_eta1
     sll_real64 :: delta_eta2
     sll_real64, dimension(:,:), pointer :: x1_cell
     sll_real64, dimension(:,:), pointer :: x2_cell
     sll_real64, dimension(:,:), pointer :: jacobians_n
     sll_real64, dimension(:,:), pointer :: jacobians_c
     character(len=64) :: label
     logical           :: written! = .false.
   contains
     procedure(geometry_function), deferred, pass       :: x1
     procedure(geometry_function), deferred, pass       :: x2
     procedure(geometry_function), deferred, pass       :: jacobian
     procedure(geometry_function_nodes), deferred, pass :: x1_at_node
     procedure(geometry_function_nodes), deferred, pass :: x2_at_node
     procedure(geometry_function_nodes), deferred, pass :: jacobian_at_node
!     procedure(j_matrix_function_nopass), pointer, nopass :: jacobian_matrix
     procedure, pass :: write_to_file => write_mapped_mesh_2d_base
  end type sll_mapped_mesh_2d_base

  
  !************************************************************************
  !
  !                       Function signatures
  !
  !************************************************************************

  !************************************************************************
  ! 2D CASE:
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
   
   abstract interface
      subroutine write_to_file_signature( map, label )
        import     :: sll_mapped_mesh_2d_base
        class(sll_mapped_mesh_2d_base)  :: map
        character(len=*), optional      :: label
      end subroutine write_to_file_signature
   end interface

contains

  subroutine write_mapped_mesh_2d_base(mesh,output_format)
    class(sll_mapped_mesh_2d_base) :: mesh
    sll_int32, optional :: output_format 
    sll_int32           :: local_format 
    sll_real64, dimension(:,:), pointer :: x1mesh
    sll_real64, dimension(:,:), pointer :: x2mesh
    sll_int32  :: i1
    sll_int32  :: i2
    sll_real64 :: eta1
    sll_real64 :: eta2
    sll_int32  :: ierr
    sll_int32  :: file_id

    if (.not. present(output_format)) then
       local_format = SLL_IO_XDMF
    else
       local_format = output_format
    end if

    if ( .not. mesh%written ) then

       select case(local_format)

       case (SLL_IO_XDMF)
          SLL_ALLOCATE(x1mesh(mesh%nc_eta1+1,mesh%nc_eta2+1), ierr)
          SLL_ALLOCATE(x2mesh(mesh%nc_eta1+1,mesh%nc_eta2+1), ierr)
          eta1 = 0.0_f64
          do i1=1, mesh%nc_eta1+1
             eta2 = 0.0_f64
             do i2=1, mesh%nc_eta2+1
                x1mesh(i1,i2) = mesh%x1_at_node(i1,i2)
                x2mesh(i1,i2) = mesh%x2_at_node(i1,i2)
                eta2 = eta2 + mesh%delta_eta2 
             end do
             eta1 = eta1 + mesh%delta_eta1
          end do
       
          call sll_xdmf_open(trim(mesh%label)//".xmf",mesh%label, &
               mesh%nc_eta1+1,mesh%nc_eta2+1,file_id,ierr)
          call sll_xdmf_write_array(mesh%label,x1mesh,"x1",ierr)
          call sll_xdmf_write_array(mesh%label,x2mesh,"x2",ierr)
          call sll_xdmf_close(file_id,ierr)

       case default
          print*, 'Not recognized format to write this mesh'
          stop

       end select

    else

       print*,' Warning, you have already written the mesh '

    end if

    mesh%written = .true.

  end subroutine write_mapped_mesh_2d_base

 end module sll_module_mapped_meshes_2d_base
