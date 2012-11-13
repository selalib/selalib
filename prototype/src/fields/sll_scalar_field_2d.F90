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
#ifdef STDF95
  use sll_module_mapped_meshes_2d
  use sll_cubic_spline_interpolator_1d
#else
  use sll_module_mapped_meshes_2d_base
  use sll_module_interpolators_1d_base
#endif
  use sll_misc_utils
  use sll_scalar_field_initializers_base

  implicit none

  type scalar_field_2d
#ifdef STDF95
     type(sll_mapped_mesh_2d_discrete), pointer  :: mesh
     type(cubic_spline_1d_interpolator), pointer :: eta1_interpolator
     type(cubic_spline_1d_interpolator), pointer :: eta2_interpolator
#else
     class(sll_mapped_mesh_2d_base), pointer  :: mesh
     class(sll_interpolator_1d_base), pointer :: eta1_interpolator
     class(sll_interpolator_1d_base), pointer :: eta2_interpolator
#endif
     sll_real64, dimension(:,:), pointer      :: data
     sll_int32                                :: data_position
     character(len=64)                        :: name
     sll_int32                                :: plot_counter
  end type scalar_field_2d

#ifndef STDF95
  abstract interface
     function scalar_function_2D( eta1, eta2 )
       use sll_working_precision
       sll_real64 :: scalar_function_2D
       sll_real64, intent(in)  :: eta1
       sll_real64, intent(in)  :: eta2
     end function scalar_function_2D
  end interface
#endif

contains   ! *****************************************************************  
  ! this used to be new_scalar_field_2d
  subroutine initialize_scalar_field_2d( &
    this, &
    field_name, &
    mesh, &
    data_position, &
    eta1_interpolator, &
    eta2_interpolator, &
    initializer )

#ifdef STDF95
    type(scalar_field_2d), intent(inout)               :: this
    type(sll_mapped_mesh_2d_discrete), pointer             :: mesh
    type(cubic_spline_1d_interpolator), pointer            :: eta1_interpolator
    type(cubic_spline_1d_interpolator), pointer            :: eta2_interpolator
    type(scalar_field_2d_initializer_base), pointer, optional :: initializer
#else
    class(scalar_field_2d), intent(inout)               :: this
    class(sll_mapped_mesh_2d_base), pointer             :: mesh
    class(sll_interpolator_1d_base), pointer            :: eta1_interpolator
    class(sll_interpolator_1d_base), pointer            :: eta2_interpolator
    class(scalar_field_2d_initializer_base), pointer, optional :: initializer
#endif
    character(len=*), intent(in)                        :: field_name
    sll_int32, intent(in)                               :: data_position

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

    ! For an initializing function, argument check should not be assertions
    ! but more permanent if-tests. There is no reason to turn these off ever.
    if( associated(eta1_interpolator) ) then
       this%eta1_interpolator => eta1_interpolator
    else
       print *, 'initialize_scalar_field_2d(): eta1_interpolator pointer ', &
            'is not associated. Exiting...'
       stop
    end if
    if( associated(eta2_interpolator) ) then
       this%eta2_interpolator => eta2_interpolator
    else
       print *, 'initialize_scalar_field_2d(): eta2_interpolator pointer ', &
            'is not associated. Exiting...'
       stop
    end if


    this%data_position = data_position
    if (data_position == NODE_CENTERED_FIELD) then
       SLL_ALLOCATE(this%data(num_pts1,num_pts2), ierr)
       if (present(initializer)) then
          call initializer%f_of_x1x2(this%data)
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
          call initializer%f_of_x1x2(this%data)
       else 
          this%data = 0.0_f64
       end if
    endif
    this%plot_counter = 0
  end subroutine initialize_scalar_field_2d

  ! The following pair of subroutines are tricky. We want them as general 
  ! services by the fields, hence we need this subroutine interface, yet
  ! we would also like a flexibility in how the derivatives are computed.
  ! A general interpolator interface would cover most of the cases, maybe
  ! all. It could be that a finite difference scheme would also work, if
  ! we ignore some of the interpolator services, like the ability to return
  ! values anywhere instead of at the nodes.
  ! For now, this interface would permit to have multiple implementations.
  subroutine compute_eta1_derivative_on_col( field2d, ith_col, deriv_out )
    type(scalar_field_2d), intent(in)    :: field2d
    sll_int32, intent(in)                :: ith_col
    sll_real64, dimension(:),intent(out) :: deriv_out
  end subroutine compute_eta1_derivative_on_col

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

#ifdef STDF95
    type(scalar_field_2d) :: scalar_field
    type(sll_mapped_mesh_2d_discrete), pointer :: mesh
#else
    class(scalar_field_2d) :: scalar_field
    class(sll_mapped_mesh_2d_base), pointer :: mesh
#endif
    logical, optional      :: multiply_by_jacobian 
    sll_int32, optional    :: output_format 
    character(len=*), optional    :: output_file_name 
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
    character(len=32) :: name
    character(len=4) :: counter
    character(len=4) :: center


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
                val(i1,i2) = scalar_field%data(i1,i2)
                eta1 = eta1 + mesh%delta_eta1
             end do
             eta2 = eta2 + mesh%delta_eta2
          end do
       end if
     
    end if


    select case(local_format)
    case (SLL_IO_XDMF)
       
       if (scalar_field%data_position == NODE_CENTERED_FIELD) then
          center = "Node"
       else if (scalar_field%data_position == CELL_CENTERED_FIELD) then
          center = "Cell"
       end if

       if (.not. present(output_file_name)) then
          scalar_field%plot_counter = scalar_field%plot_counter+1
          call int2string(scalar_field%plot_counter, counter)
          name = trim(scalar_field%name)//counter
       else 
          name = output_file_name
       end if
       call sll_xdmf_open(trim(name)//".xmf", &
            scalar_field%mesh%label,        &
            num_pts1,num_pts2,file_id,ierr)
      
       call sll_xdmf_write_array(trim(name), &
                                 val,&
                                 scalar_field%name,ierr,file_id, &
                                 "Node")
       call sll_xdmf_close(file_id,ierr)

    case default

       print*, "write_scalar_field_2d: requested output format not recognized."
       stop
    end select


    SLL_DEALLOCATE_ARRAY(val,ierr)
  end subroutine write_scalar_field_2d

end module sll_scalar_field_2d
