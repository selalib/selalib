!------------------------------------------------------------------------------
! SELALIB
!------------------------------------------------------------------------------
!
! MODULE: sll_scalar_field_1d
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
module sll_scalar_field_1d
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_file_io.h"
  use sll_constants
  use sll_scalar_field_initializers_base
  use sll_module_mapped_meshes_1d_base
  use sll_utilities
  implicit none
  
  !I removed this line because, it not exists in 2d version
  !should be in base module
  !integer, parameter :: NODE_CENTERED_FIELD = 0, CELL_CENTERED_FIELD = 1

  type scalar_field_1d
     class(sll_mapped_mesh_1d_base), pointer :: mesh
     sll_real64, dimension(:), pointer       :: data
     sll_int32                               :: data_position
     character(len=64)                       :: name
  end type scalar_field_1d

  abstract interface
     function scalar_function_1D( eta1 )
       use sll_working_precision
       sll_real64 :: scalar_function_1D
       sll_real64, intent(in)  :: eta1
     end function scalar_function_1D
  end interface

contains   ! *i****************************************************************  
  ! this used to be new_scalar_field_1d
  subroutine initialize_scalar_field_1d( &
    this, &
    field_name, &
    mesh, &
    data_position, &
    init_function)

    class(scalar_field_1d), intent(inout)               :: this
    class(sll_mapped_mesh_1d_base), pointer             :: mesh
    procedure(scalar_function_1D), optional             :: init_function
    character(len=*), intent(in)                        :: field_name
    sll_int32, intent(in)                               :: data_position
    sll_int32  :: ierr
    sll_int32  :: num_cells1
    sll_int32  :: num_pts1
    sll_int32  :: i1
    sll_real64 :: eta1
    sll_real64 :: delta1

    this%mesh => mesh
    this%name  = trim(field_name)
    num_cells1 = mesh%nc_eta1
    num_pts1   = mesh%nc_eta1+1

    this%data_position = data_position
    if (data_position == NODE_CENTERED_FIELD) then
       SLL_ALLOCATE(this%data(num_pts1), ierr)
       if (present(init_function)) then
          do i1 = 1, num_pts1
             this%data(i1) = init_function( mesh%x1_at_node(i1) )
          end do
       else 
          this%data = 0.0_f64 ! initialize to zero
       end if
    else if (data_position == CELL_CENTERED_FIELD) then
       SLL_ALLOCATE(this%data(num_cells1), ierr)
       if (present(init_function)) then
          delta1 = 1.0_f64/real(num_cells1,f64)
          eta1   = 0.5_f64 * delta1
          do i1 = 1, num_cells1
             this%data(i1) = init_function( mesh%x1(eta1) )
             eta1 = eta1 + delta1
          end do
       else 
          this%data = 0.0_f64 ! initialize to zero
       end if
    endif
  end subroutine initialize_scalar_field_1d

  ! need to do something about deallocating the field proper, when allocated
  ! in the heap...
  subroutine delete_scalar_field_1d( this )
    type(scalar_field_1d), pointer :: this
    sll_int32                      :: ierr
    nullify(this%mesh)
    SLL_DEALLOCATE(this%data, ierr)
  end subroutine delete_scalar_field_1d


  subroutine write_scalar_field_1d( &
    scalar_field, &
    multiply_by_jacobian, &
    output_file_name, &
    output_format)
    class(scalar_field_1d)                  :: scalar_field
    !sll_real64, dimension(:), pointer      :: x1_array
    logical, optional                       :: multiply_by_jacobian 
    sll_int32, optional                     :: output_format 
    character(len=*), optional              :: output_file_name
    character(len=64)                       :: file_name

    sll_int32                               :: error

    if(present(multiply_by_jacobian))then
      print*,'multiply_by_jacobian option is not implemented'
    endif

    if(present(output_format))then
      print*,'There is just gnuplot format available'
    endif
    if(present(output_file_name))then
      file_name = output_file_name
    else
      file_name = scalar_field%name
    endif 

    call sll_gnuplot_write(scalar_field%data,file_name,error)
  end subroutine write_scalar_field_1d

end module sll_scalar_field_1d
