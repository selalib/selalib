!------------------------------------------------------------------------------
! SELALIB
!------------------------------------------------------------------------------
!
! MODULE: sll_advection_field
!
!> @author
!> - Eric
!> - Michel
!> - Pierre
!> - Edwin
!>
!
! DESCRIPTION: 
!
!> @brief
!> Implements the distribution function types
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
!>    - splines
!>    - mesh_types
!
! REVISION HISTORY:
! 21_05_2012 - Initial Version
! TODO_dd_mmm_yyyy - TODO_describe_appropriate_changes - TODO_name
!------------------------------------------------------------------------------
module sll_advection_field
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_field_2d.h"
!#include "sll_field_1d.h"
  use sll_constants
  use sll_utilities   ! for int2string
  use sll_scalar_field_1d
  implicit none
  
  type, extends(scalar_field_2d) :: hamiltonian_advection_field_2d
     sll_real64      :: pmass
     sll_real64      :: pcharge
  end type hamiltonian_advection_field_2d

contains

  subroutine initialize_advection_field_2d( &
    this, &
    mass, &
    charge, &
    field_name, &
    mesh, &
    data_position, &
    initializer, &
    eta1_interpolator, &
    eta2_interpolator )

    type(hamiltonian_advection_field_2d), intent(inout) :: this
    sll_real64, intent(in)                              :: mass
    sll_real64, intent(in)                              :: charge
    character(len=*), intent(in)                        :: field_name
    sll_int32, intent(in)                               :: data_position
    class(scalar_field_2d_initializer_base), pointer, optional :: initializer
    class(sll_interpolator_1d_base), pointer            :: eta1_interpolator
    class(sll_interpolator_1d_base), pointer            :: eta2_interpolator
    !class(sll_mapped_mesh_2d_base), pointer             :: mesh
    class(sll_coordinate_transformation_2d_base), pointer   :: mesh

    this%pmass = mass
    this%pcharge = charge
    call initialize_scalar_field_2d( &
         this, &
         field_name, &
         mesh, &
         data_position, &
         eta1_interpolator, &
         eta2_interpolator, &
         initializer )
  end subroutine initialize_advection_field_2d

  !> sets advection field to Hamiltonian function at nodes from given 
  !> self-consistent and external potentials. 
  subroutine compute_hamiltonian(this, phi_self, phi_external)
    type(hamiltonian_advection_field_2d), intent(inout)   :: this
    type(scalar_field_1d), intent(in)                  :: phi_self
    type(scalar_field_1d), intent(in), optional        :: phi_external
    type(sll_logical_mesh_2d), pointer                 :: mesh
    sll_int32 :: nc_eta1
    sll_int32 :: nc_eta2
    sll_int32 :: i1, i2
    sll_real64 :: mass

    mesh => this%transf%mesh
    nc_eta1 = mesh%num_cells1 
    nc_eta2 = mesh%num_cells2
    mass    = this%pmass    
    if (present(phi_external)) then 
       do i1 = 1, nc_eta1+1
          do i2 = 1, nc_eta2+1
             this%data(i1,i2) = 0.5_f64*mass * this%transf%x2_at_node(i1,i2)**2 &
                  + this%pcharge*phi_self%data(i1)  &
                  + phi_external%data(i1)
          end do
       end do
    else 
       do i1 = 1, nc_eta1
          do i2 = 1, nc_eta2
             this%data(i1,i2) = 0.5_f64*mass * this%transf%x2_at_node(i1,i2)**2 &
                  + this%pcharge*phi_self%data(i1)
          end do
       end do
    end if
  end subroutine compute_hamiltonian
    

#if 0
  ! The following function computes the derivative:
  !
  !                         partial H(eta1,eta2)
  !                        ---------------------
  !                           partial   eta2
  !
  ! along the ith column (constant eta2), returning the result in the array 
  ! dh_dx1. This is written with computations along rows/columns in mind.
  !
  !          p (rows, separated by long stride)
  !         (eta2)
  !          ^
  !          |
  !          +----+----+----+----+----+
  !          |    |    |    |    |    |
  !          |----+----+----+----+----|
  !          |    |    |    |    |    |
  !          |----+----+----+----+----|
  !          |    |    |    |    |    |
  !          |----+----+----+----+----|
  !          |    |    |    |    |    |
  !          +----+----+----+----+----+ ---> q (columns in contiguous memory)
  !                                          (eta1)
  subroutine compute_dh_deta2( data2d, ith_row, row_size, delta, dh_deta2 )
    sll_real64, dimension(:,:), intent(in), target  :: data2d
    sll_int32, intent(in)                           :: ith_row
    sll_int32, intent(in)                           :: row_size
    sll_real64, intent(in)                          :: delta     ! cell spacing
    sll_real64, dimension(:), intent(out)           :: dh_deta2
    sll_int32  :: i
    sll_real64 :: r_delta  ! reciprocal of delta
    sll_real64, dimension(:), pointer :: dptr

    SLL_ASSERT( size(data2d(ith_row,:)) >= row_size )
    SLL_ASSERT( row_size <= size(dh_deta2) )

    r_delta = 1.0_f64/delta
    dptr => data2d(ith_row,:)

    ! Compute derivative in first point with a forward scheme (-3/2, 2, -1/2)
    dh_deta2(1) = r_delta*(-1.5_f64*dptr(1) + 2.0_f64*dptr(2) - 0.5_f64*dptr(3))

    ! Compute derivative in the bulk of the array with a centered scheme.
    do i=2,row_size-1
       dh_deta2(i) = 0.5_f64*r_delta*(dptr(i+1)-dptr(i-1))
    end do

    ! Compute derivative in the last point with a backward scheme (1/2, -2, 3/2)
    dh_deta2(row_size) = r_delta*(0.5_f64*dptr(row_size-2) - &
                         2.0_f64*dptr(row_size-1) + 1.5_f64*dptr(row_size))
  end subroutine compute_dh_deta2


  ! The following function computes the derivative:
  !
  !                         partial H(eta1,eta2)
  !                        ---------------------
  !                           partial   eta1
  !
  ! along the ith row (constant eta2), returning the result in the array 
  ! dh_deta1. This is written with computations along lines in mind.

  subroutine compute_dh_deta1( data2d, ith_col, col_size, delta, dh_deta1 )
    sll_real64, dimension(:,:), intent(in), target  :: data2d
    sll_int32, intent(in)                           :: ith_col
    sll_int32, intent(in)                           :: col_size
    sll_real64, intent(in)                          :: delta     ! cell spacing
    sll_real64, dimension(:), intent(out)           :: dh_deta1
    sll_int32  :: i
    sll_real64 :: r_delta  ! reciprocal of delta
    sll_real64, dimension(:), pointer :: dptr

    SLL_ASSERT( size(data2d(ith_col,:)) >= col_size )
    SLL_ASSERT( col_size <= size(dh_deta1) )

    r_delta = 1.0_f64/delta
    dptr => data2d(ith_col,:)

    ! Compute derivative in first point with a forward scheme (-3/2, 2, -1/2)
    dh_deta1(1) = r_delta*(-1.5_f64*dptr(1) + 2.0_f64*dptr(2) - 0.5_f64*dptr(3))

    ! Compute derivative in the bulk of the array with a centered scheme.
    do i=2,col_size-1
       dh_deta1(i) = 0.5_f64*r_delta*(dptr(i+1)-dptr(i-1))
    end do

    ! Compute derivative in the last point with a backward scheme (1/2, -2, 3/2)
    dh_deta1(col_size) = r_delta*(0.5_f64*dptr(col_size-2) - &
                       2.0_f64*dptr(col_size-1) + 1.5_f64*dptr(col_size))
  end subroutine compute_dh_deta1
#endif

end module sll_advection_field
