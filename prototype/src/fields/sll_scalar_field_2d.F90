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
  use sll_xdmf
  use numeric_constants
  use sll_mapped_mesh_base
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
    mesh, &
    data_position, &
    init_function)

    class(scalar_field_2d), intent(inout)               :: this
    class(sll_mapped_mesh_2d_base), pointer, intent(in) :: mesh
    sll_int32, intent(in)                               :: data_position
    procedure(scalar_function_2D), pointer              :: init_function
    sll_int32  :: ierr
    sll_int32  :: num_cells1
    sll_int32  :: num_cells2
    sll_int32  :: num_pts1
    sll_int32  :: num_pts2
    sll_int32  :: i1, i2
    sll_real64 :: eta1, eta2
    sll_real64 :: delta1, delta2

    this%mesh => mesh
    num_cells1 = mesh%nc_eta1
    num_cells2 = mesh%nc_eta2
    num_pts1   = mesh%nc_eta1+1
    num_pts2   = mesh%nc_eta2+1

    this%data_position = data_position
    if (data_position == NODE_CENTERED_FIELD) then
       SLL_ALLOCATE(this%data(num_pts1,num_pts2), ierr)
       do i2 = 1, num_pts2
          do i1 = 1, num_pts1
             this%data(i1,i2) = init_function( mesh%x1_at_node(i1,i2), &
                                               mesh%x2_at_node(i1,i2) )
          end do
       end do
    else if (data_position == CELL_CENTERED_FIELD) then
       SLL_ALLOCATE(this%data(num_cells1,num_cells2), ierr)
       delta1 = 1.0_f64/real(num_cells1,f64)
       delta2 = 1.0_f64/real(num_cells2,f64)
       eta1   = 0.5_f64 * delta1
       eta2   = 0.5_f64 * delta2
       do i2 = 1, num_cells2
          do i1 = 1, num_cells1
             this%data(i1,i2) = init_function( mesh%x1(eta1,eta2), &
                                               mesh%x2(eta1,eta2) )
             eta1 = eta1 + delta1
          end do
          eta2 = eta2 + delta2
       end do
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
    name, &
    multiply_by_jacobian, &
    file_id ) ! por aqui: this is intent(out) in sll_xdmf.F90
    class(scalar_field_2d) :: scalar_field
    character(len=*)       :: name
    logical, optional      :: multiply_by_jacobian 
    class(sll_mapped_mesh_2d_base), pointer :: mesh
    sll_int32  :: i1
    sll_int32  :: i2
    sll_real64 :: eta1
    sll_real64 :: eta2
    sll_real64 :: avg
    sll_int32  :: ierr
    sll_real64, dimension(:,:), pointer :: val
    sll_int32  :: num_pts1
    sll_int32  :: num_pts2

    mesh => scalar_field%mesh
    num_pts1 = mesh%nc_eta1+1
    num_pts2 = mesh%nc_eta2+2

    SLL_ASSERT(associated(mesh))    
    SLL_ALLOCATE(val(num_pts1,num_pts2), ierr)

    call sll_xdmf_open( &
         scalar_field%name, &
         
    if (.not.(present(multiply_by_jacobian))) then
       call write_vec1d( &
            scalar_field%data, &
            num_pts1, &
            num_pts2, &
            name, &
            "mesh", &
            scalar_field%data_position)
    else !if (multiply_by_jacobian) then 
       eta1 = 0.5_f64 * mesh%delta_eta1
       do i1 = 1, mesh%nc_eta1
          eta2 =  0.5_f64 * mesh%delta_eta2
          do i2 = 1, mesh%nc_eta2
             val(i1,i2) = scalar_field%data(i1,i2) / mesh%jacobian(eta1, eta2)
             eta2 = eta2 + mesh%delta_eta2
          end do
          eta1 = eta1 + mesh%delta_eta1
       end do
       call write_vec1d( &
            val, &
            num_pts1, &
            num_pts2, &
            name, &
            "mesh", &
            scalar_field%data_position)
    end if
    SLL_DEALLOCATE(val,ierr)
  end subroutine write_scalar_field_2d

end module sll_scalar_field_2d
