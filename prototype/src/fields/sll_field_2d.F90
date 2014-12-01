!------------------------------------------------------------------------------
! SELALIB
!------------------------------------------------------------------------------
!
! MODULE: sll_field_2d
!
!> @author
!> - Edwin
!> - Pierre
!>
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
!>    - diagnostics
!
! REVISION HISTORY:
! DD Mmm YYYY - Initial Version
! TODO_dd_mmm_yyyy - TODO_describe_appropriate_changes - TODO_name
!------------------------------------------------------------------------------
module sll_field_2d
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"

  use sll_constants
  use sll_diagnostics
  !use sll_mapped_mesh_base
  use sll_coordinate_transformation_2d_base_module
  use sll_file_io
  implicit none

  enum, bind(C)
     enumerator :: NODE_FIELD = 0, CELL_CENTER_FIELD = 1
  end enum

  type vec2
     sll_real64 :: v1
     sll_real64 :: v2
  end type vec2

  type scalar_field_2d
     !class(sll_mapped_mesh_2d_base), pointer :: mesh
     class(sll_coordinate_transformation_2d_base), pointer   :: mesh
     sll_real64, dimension(:,:), pointer :: data
     sll_int32  :: data_position
     character(len=30) :: name
   !contains
   !  procedure :: get_data
  end type scalar_field_2d

  type vector_field_2d
     !class(sll_mapped_mesh_2d_base), pointer :: mesh
     class(sll_coordinate_transformation_2d_base), pointer   :: mesh
     type(vec2), dimension(:,:), pointer :: data
     sll_int32  :: data_position
     character(len=30) :: name
  end type vector_field_2d

  abstract interface
     function scalar_function_field_2d( eta1, eta2 )
       use sll_working_precision
       sll_real64 :: scalar_function_field_2d
       sll_real64, intent(in)  :: eta1
       sll_real64, intent(in)  :: eta2
     end function scalar_function_field_2d
  end interface

contains   ! *****************************************************************  
 
  subroutine new_scalar_field_2d( this, mesh, data_position, data_func)
    class(scalar_field_2d), intent(inout)   :: this
    !class(sll_mapped_mesh_2d_base), pointer, intent(in) :: mesh
    class(sll_coordinate_transformation_2d_base), pointer   :: mesh
    sll_int32 :: data_position
    procedure(scalar_function_field_2d), pointer :: data_func
    ! local variables
    sll_int32                :: ierr
    sll_int32  :: i1, i2
    sll_real64 :: eta1, eta2
    sll_real64 :: delta1, delta2

    this%mesh => mesh
    this%data_position = data_position
    if (data_position == NODE_FIELD) then
       SLL_ALLOCATE(this%data(mesh%nc_eta1+1,mesh%nc_eta2+1), ierr)
       do i2 = 1, mesh%nc_eta2+1
          do i1 = 1, mesh%nc_eta1+1
             this%data(i1,i2) = data_func(mesh%x1_at_node(i1,i2), &
                  mesh%x2_at_node(i1,i2))
          end do
       end do
    else if (data_position == CELL_CENTER_FIELD) then
       SLL_ALLOCATE(this%data(mesh%nc_eta1,mesh%nc_eta2), ierr)
       delta1 = 1.0_f64/mesh%nc_eta1
       delta2 = 1.0_f64/mesh%nc_eta2
       eta1 = 0.5_f64 * delta1
       eta2 = 0.5_f64 * delta2
       do i2 = 1, mesh%nc_eta2
          do i1 = 1, mesh%nc_eta1
             this%data(i1,i2) = data_func(mesh%x1(eta1,eta2), &
                  mesh%x2(eta1,eta2))
             eta1 = eta1 + delta1
          end do
          eta2 = eta2 + delta2
       end do
    endif
    
  end subroutine new_scalar_field_2d

!!$  function get_data(this, i1, i2)
!!$    class (scalar_field_2d) :: this
!!$    sll_int32  :: i1, i2
!!$    sll_real64 :: get_data(this%vec_length)
!!$    get_data = this%data(:,i1,i2)
!!$  end function get_data

  subroutine delete_scalar_field_2d( this )
    type(scalar_field_2d), pointer :: this
    sll_int32                    :: ierr
    nullify(this%mesh)
    SLL_DEALLOCATE(this%data, ierr)
  end subroutine delete_scalar_field_2d

 
!!$  function new_field_2D_vec2( mesh_descriptor )
!!$    type(field_2D_vec2), pointer      :: new_field_2D_vec2
!!$    type(mesh_descriptor_2D), pointer :: mesh_descriptor
!!$    sll_int32                         :: ierr
!!$    SLL_ASSERT(associated(mesh_descriptor))
!!$    SLL_ALLOCATE(new_field_2D_vec2, ierr)
!!$    new_field_2D_vec2%descriptor => mesh_descriptor
!!$    SLL_ALLOCATE(new_field_2D_vec2%data(1:mesh_descriptor%nc_eta1+1,1:mesh_descriptor%nc_eta2+1),ierr)
!!$  end function new_field_2D_vec2
!!$
!!$  subroutine delete_field_2D_vec2( f2Dv2 )
!!$    type(field_2D_vec2), pointer :: f2Dv2
!!$    sll_int32                    :: ierr
!!$    if( .not. (associated(f2Dv2))) then
!!$       write (*,'(a)') 'ERROR: delete_field_1D_vec1(), not associated argument.'
!!$       STOP
!!$    end if
!!$    nullify(f2Dv2%descriptor)
!!$    SLL_DEALLOCATE(f2Dv2%data, ierr)
!!$    SLL_DEALLOCATE(f2Dv2, ierr)
!!$  end subroutine delete_field_2D_vec2




  subroutine write_sll_mapped_mesh_2d_base(mesh,mesh_name)
    !class(sll_mapped_mesh_2d_base) :: mesh
    class(sll_coordinate_transformation_2d_base), pointer   :: mesh
    character(len=*), intent(in), optional :: mesh_name
    sll_real64, dimension(:,:), pointer :: x1mesh
    sll_real64, dimension(:,:), pointer :: x2mesh
    sll_int32  :: i1
    sll_int32  :: i2
    sll_real64 :: eta1
    sll_real64 :: eta2
    sll_int32 :: ierr

    ! create 2D mesh
    SLL_ALLOCATE(x1mesh(mesh%mesh%num_cell1+1,mesh%mesh%num_cell2+1), ierr)
    SLL_ALLOCATE(x2mesh(mesh%mesh%num_cell1+1,mesh%mesh%num_cell2+1), ierr)
    eta1 = 0.0_f64
    do i1=1, mesh%mesh%num_cell1+1
       eta2 = 0.0_f64
       do i2=1, mesh%mesh%num_cell2+1
          x1mesh(i1,i2) = mesh%x1_at_node(i1,i2)
          x2mesh(i1,i2) = mesh%x2_at_node(i1,i2)
          eta2 = eta2 + mesh%delta_eta2 
       end do
       eta1 = eta1 + mesh%delta_eta1
    end do
    if (present(mesh_name)) then
       call write_mesh(x1mesh,x2mesh,mesh%mesh%num_cell1+1,mesh%mesh%num_cell2+1,trim(mesh_name))
    else
       call write_mesh(x1mesh,x2mesh,mesh%mesh%num_cell1+1,mesh%mesh%num_cell2+1,"mesh")
    end if
  end subroutine write_sll_mapped_mesh_2d_base

  subroutine write_scalar_field_2d( f2Dv1, name, jacobian)
    class(scalar_field_2d) :: f2Dv1
    character(len=*) :: name
    logical, optional  :: jacobian ! .true. if field data multiplied by jacobian is stored
    
    !class(sll_mapped_mesh_2d_base), pointer :: mesh
    class(sll_coordinate_transformation_2d_base), pointer   :: mesh
    sll_int32  :: i1
    sll_int32  :: i2
    sll_real64 :: eta1
    sll_real64 :: eta2
    sll_real64 :: avg
    sll_int32 :: ierr

    sll_real64, dimension(:,:), pointer :: val

    ! create 2D mesh
    mesh => f2Dv1%mesh
    SLL_ASSERT(associated(mesh))    
    SLL_ALLOCATE(val(mesh%mesh%num_cell1 + 1,mesh%mesh%num_cell2 + 1), ierr)

    if (.not.(present(jacobian))) then
       call write_vec1d(f2Dv1%data,mesh%mesh%num_cell1+1,mesh%mesh%num_cell2+1,name, &
            "mesh",f2Dv1%data_position)
    else if (jacobian) then 
       ! quantity multiplied by Jacobian is stored, need to divide by jacobian for
       eta1 = 0.5_f64 * mesh%delta_eta1
       do i1 = 1, mesh%mesh%num_cell1
          eta2 =  0.5_f64 * mesh%delta_eta2
          do i2 = 1, mesh%mesh%num_cell2
             val(i1,i2) = f2Dv1%data( i1,i2) / mesh%jacobian (eta1, eta2)
             eta2 = eta2 + mesh%delta_eta2
          end do
          eta1 = eta1 + mesh%delta_eta1
       end do
       call write_vec1d(val,mesh%mesh%num_cell1+1,mesh%mesh%num_cell2+1,name,"mesh",f2Dv1%data_position)
    end if
  end subroutine write_scalar_field_2d
end module sll_field_2d
