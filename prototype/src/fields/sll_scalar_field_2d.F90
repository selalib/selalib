!------------------------------------------------------------------------------
! SELALIB
!------------------------------------------------------------------------------
!
! MODULE: sll_scalar_field_2d
!
!> @author
!> - Edwin
!> - Pierre
!> - Eric
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
!
! REVISION HISTORY:
! DD Mmm YYYY - Initial Version
! TODO_dd_mmm_yyyy - TODO_describe_appropriate_changes - TODO_name
!------------------------------------------------------------------------------
module sll_scalar_field_2d
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
  use sll_io
  use numeric_constants
  use sll_module_mapped_meshes_2d_base
  use sll_scalar_field_initializers_base
  use sll_misc_utils
  implicit none

  enum, bind(C)
     enumerator :: NODE_CENTERED_FIELD = 0, CELL_CENTERED_FIELD = 1
  end enum

  type scalar_field_2d
     class(sll_mapped_mesh_2d_base), pointer :: mesh
     sll_real64, dimension(:,:), pointer     :: data
     sll_int32                               :: data_position
     character(len=64)                       :: name
  end type scalar_field_2d

  abstract interface
     function scalar_function_2D( eta1, eta2 )
       use sll_working_precision
       sll_real64 :: scalar_function_2D
       sll_real64, intent(in)  :: eta1
       sll_real64, intent(in)  :: eta2
     end function scalar_function_2D
  end interface

contains   ! *****************************************************************  
  ! this used to be new_scalar_field_2d
  subroutine initialize_scalar_field_2d( &
    this, &
    field_name, &
    mesh, &
    data_position, &
    initializer)

    class(scalar_field_2d), intent(inout)               :: this
    character(len=*), intent(in)                        :: field_name
    class(sll_mapped_mesh_2d_base), pointer             :: mesh
    sll_int32, intent(in)                               :: data_position
    class(scalar_field_2d_initializer_base), pointer, optional    :: initializer

    sll_int32  :: ierr
    sll_int32  :: num_cells1
    sll_int32  :: num_cells2
    sll_int32  :: num_pts1
    sll_int32  :: num_pts2
    sll_int32  :: i1, i2
    sll_real64 :: eta1, eta2
    sll_real64 :: delta1, delta2

    SLL_ASSERT(associated(mesh))
    this%mesh => mesh
    this%mesh%written = .false.
    this%name  = trim(field_name)
    num_cells1 = mesh%nc_eta1
    num_cells2 = mesh%nc_eta2
    num_pts1   = mesh%nc_eta1+1
    num_pts2   = mesh%nc_eta2+1

    this%data_position = data_position
    if (data_position == NODE_CENTERED_FIELD) then
       SLL_ALLOCATE(this%data(num_pts1,num_pts2), ierr)
       if (present(initializer)) then
          call initializer%f_of_x1x2(mesh,this%data)
       else 
          this%data = 0.0_f64
       end if
    else if (data_position == CELL_CENTERED_FIELD) then
       SLL_ALLOCATE(this%data(num_cells1,num_cells2), ierr)
       delta1 = 1.0_f64/real(num_cells1,f64)
       delta2 = 1.0_f64/real(num_cells2,f64)
       eta1   = 0.5_f64 * delta1
       eta2   = 0.5_f64 * delta2
       if (present(initializer)) then
          call initializer%f_of_x1x2(mesh,this%data)
       else 
          this%data = 0.0_f64
       end if
    endif
  end subroutine initialize_scalar_field_2d

  ! need to do something about deallocating the field proper, when allocated
  ! in the heap...
  subroutine delete_scalar_field_2d( this )
    type(scalar_field_2d), pointer :: this
    sll_int32                      :: ierr
    nullify(this%mesh)
    SLL_DEALLOCATE(this%data, ierr)
  end subroutine delete_scalar_field_2d

  subroutine write_scalar_field_2d( &
    scalar_field, &
    multiply_by_jacobian, &
    output_file_name, &
    output_format)

    class(scalar_field_2d) :: scalar_field
    logical, optional      :: multiply_by_jacobian 
    sll_int32, optional    :: output_format 
    character(len=*), optional    :: output_file_name 
    class(sll_mapped_mesh_2d_base), pointer :: mesh
    sll_int32              :: local_format 

    sll_int32  :: i1
    sll_int32  :: i2
    sll_real64 :: eta1
    sll_real64 :: eta2
    sll_real64 :: avg
    sll_int32  :: ierr
    sll_real64, dimension(:,:), allocatable :: val
    sll_int32  :: num_pts1
    sll_int32  :: num_pts2
    sll_int32  :: file_id

    if (.not. present(output_format)) then
       local_format = SLL_IO_XDMF
    else
       local_format = output_format
    end if

    mesh => scalar_field%mesh

    SLL_ASSERT(associated(mesh))  
    if (.not. mesh%written) then
       call mesh%write_to_file(local_format)
    end if

    num_pts1 = mesh%nc_eta1+1
    num_pts2 = mesh%nc_eta2+1
    if (scalar_field%data_position == NODE_CENTERED_FIELD) then
       SLL_ALLOCATE(val(num_pts1,num_pts2), ierr)
    else
       SLL_ALLOCATE(val(num_pts1-1,num_pts2-1), ierr)
    end if

    if (.not.(present(multiply_by_jacobian))) then
       val =  scalar_field%data
    else !if (multiply_by_jacobian) then 

       if (scalar_field%data_position == CELL_CENTERED_FIELD) then
          eta2 =  0.5_f64 * mesh%delta_eta2
          do i2 = 1, mesh%nc_eta2
             eta1 = 0.5_f64 * mesh%delta_eta1
             do i1 = 1, mesh%nc_eta1
                val(i1,i2) = scalar_field%data(i1,i2) / mesh%jacobian(eta1, eta2)
                eta1 = eta1 + mesh%delta_eta1
             end do
             eta2 = eta2 + mesh%delta_eta2
          end do
       else
          eta2 =  0.0_f64 
          do i2 = 1, num_pts2
             eta1 = 0.0_f64 
             do i1 = 1, num_pts1
                val(i1,i2) = scalar_field%data(i1,i2) / mesh%jacobian(eta1, eta2)
                eta1 = eta1 + mesh%delta_eta1
             end do
             eta2 = eta2 + mesh%delta_eta2
          end do
       end if
     
    end if

    select case(local_format)
    case (SLL_IO_XDMF)
       
       if (.not. present(output_file_name)) then
       call sll_xdmf_open(  &
            trim(scalar_field%name)//".xmf", &
            scalar_field%mesh%label,                              &
            num_pts1,num_pts2,file_id,ierr)
       else
       call sll_xdmf_open(  &
            trim(output_file_name)//".xmf", &
            scalar_field%mesh%label,                              &
            num_pts1,num_pts2,file_id,ierr)
       end if

       if (scalar_field%data_position == NODE_CENTERED_FIELD) then
          call sll_xdmf_write_array(scalar_field%mesh%label, &
                                    val,&
                                    scalar_field%name,ierr,file_id, &
                                    "Node")
       else if (scalar_field%data_position == CELL_CENTERED_FIELD) then
          call sll_xdmf_write_array(scalar_field%mesh%label, &
                                    val,&
                                    scalar_field%name,ierr,file_id, &
                                    "Cell")
       end if
       call sll_xdmf_close(file_id,ierr)

    case default

       print*, "No recognized output format"
       stop
    end select


    SLL_DEALLOCATE_ARRAY(val,ierr)
  end subroutine write_scalar_field_2d

end module sll_scalar_field_2d
