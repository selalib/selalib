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
#include "sll_field_1d.h"
  use numeric_constants
  use sll_misc_utils   ! for int2string
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
    data_position )

    type(hamiltonian_advection_field_2d), intent(inout) :: this
    sll_real64, intent(in)                              :: mass
    sll_real64, intent(in)                              :: charge
    character(len=*), intent(in)                        :: field_name
    class(sll_mapped_mesh_2d_base), target              :: mesh
    sll_int32, intent(in)                               :: data_position

    this%pmass = mass
    this%pcharge = charge
    call initialize_scalar_field_2d( &
         this, &
         field_name, &
         mesh, &
         data_position)
  end subroutine initialize_advection_field_2d

  !> sets advection field to Hamiltonian function at nodes from given 
  !> self-consistent and external potentials
  subroutine compute_hamiltonian(this, phi_self, phi_external)
    type(hamiltonian_advection_field_2d), intent(inout)   :: this
    type(scalar_field_1d), intent(in)                  :: phi_self
    type(scalar_field_1d), intent(in), optional        :: phi_external

    sll_int32 :: nc_eta1
    sll_int32 :: nc_eta2
    sll_int32 :: i1, i2
    sll_real64 :: mass

    nc_eta1 = GET_FIELD_NC_ETA1( this ) 
    nc_eta2 = GET_FIELD_NC_ETA2( this )     
    mass    = this%pmass    
    if (present(phi_external)) then 
       do i1 = 1, nc_eta1+1
          do i2 = 1, nc_eta2+1 
             this%data(i1,i2) = 0.5_f64*mass * this%mesh%x2_at_node(i1,i2)**2 &
                  + this%pcharge*(FIELD_1D_AT_I(phi_self,i1)  &
                  + FIELD_1D_AT_I(phi_external,i1))
          end do
       end do
    else 
       do i1 = 1, nc_eta1
          do i2 = 1, nc_eta2
             this%data(i1,i2) = 0.5_f64*mass * this%mesh%x2_at_node(i1,i2)**2 &
                  + this%pcharge*FIELD_1D_AT_I(phi_self,i1)
          end do
       end do
    end if
  end subroutine compute_hamiltonian
    


  ! The following function computes the derivative:
  !
  !                         partial H(q,p)
  !                        ---------------
  !                           partial p
  !
  ! Which in our case is:
  !                         partial H(x1,x2)
  !                        ---------------
  !                           partial x2
  !
  ! along the ith column (constant q), returning the result in the array dh_dx1.
  ! This is written with computations along rows/columns in mind.
  !
  !          p (rows, separated by long stride)
  !         (x2)
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
  !                                          (x1)
  subroutine compute_dh_dx2( data2d, ith_row, row_size, delta, dh_dx2 )
    sll_real64, dimension(:,:), intent(in), target  :: data2d
    sll_int32, intent(in)                           :: ith_row
    sll_int32, intent(in)                           :: row_size
    sll_real64, intent(in)                          :: delta     ! cell spacing
    sll_real64, dimension(:), intent(out)           :: dh_dx2
    sll_int32  :: i
    sll_real64 :: r_delta  ! reciprocal of delta
    sll_real64, dimension(:), pointer :: dptr

    SLL_ASSERT( size(data2d(ith_row,:)) >= row_size )
    SLL_ASSERT( row_size <= size(dh_dx2) )

    r_delta = 1.0_f64/delta
    dptr => data2d(ith_row,:)

    ! Compute derivative in first point with a forward scheme (-3/2, 2, -1/2)
    dh_dx2(1) = r_delta*(-1.5_f64*dptr(1) + 2.0_f64*dptr(2) - 0.5_f64*dptr(3))

    ! Compute derivative in the bulk of the array with a centered scheme.
    do i=2,row_size-1
       dh_dx2(i) = 0.5_f64*r_delta*(dptr(i+1)-dptr(i-1))
    end do

    ! Compute derivative in the last point with a backward scheme (1/2, -2, 3/2)
    dh_dx2(row_size) = r_delta*(0.5_f64*dptr(row_size-2) - &
                       2.0_f64*dptr(row_size-1) + 1.5_f64*dptr(row_size))
  end subroutine compute_dh_dx2


  ! The following function computes the derivative:
  !
  !                         partial H(q,p)
  !                        ---------------
  !                           partial q
  !
  ! which in this case is:
  !
  !                         partial H(x1,x2)
  !                        ---------------
  !                           partial x1
  !
  ! along the ith row (constant p), returning the result in the array dh_dx1.
  ! This is written with computations along lines in mind.

  subroutine compute_dh_dx1( data2d, ith_col, col_size, delta, dh_dx1 )
    sll_real64, dimension(:,:), intent(in), target  :: data2d
    sll_int32, intent(in)                           :: ith_col
    sll_int32, intent(in)                           :: col_size
    sll_real64, intent(in)                          :: delta     ! cell spacing
    sll_real64, dimension(:), intent(out)           :: dh_dx1
    sll_int32  :: i
    sll_real64 :: r_delta  ! reciprocal of delta
    sll_real64, dimension(:), pointer :: dptr

    SLL_ASSERT( size(data2d(ith_col,:)) >= col_size )
    SLL_ASSERT( col_size <= size(dh_dx1) )

    r_delta = 1.0_f64/delta
    dptr => data2d(ith_col,:)

    ! Compute derivative in first point with a forward scheme (-3/2, 2, -1/2)
    dh_dx1(1) = r_delta*(-1.5_f64*dptr(1) + 2.0_f64*dptr(2) - 0.5_f64*dptr(3))

    ! Compute derivative in the bulk of the array with a centered scheme.
    do i=2,col_size-1
       dh_dx1(i) = 0.5_f64*r_delta*(dptr(i+1)-dptr(i-1))
    end do

    ! Compute derivative in the last point with a backward scheme (1/2, -2, 3/2)
    dh_dx1(col_size) = r_delta*(0.5_f64*dptr(col_size-2) - &
                       2.0_f64*dptr(col_size-1) + 1.5_f64*dptr(col_size))
  end subroutine compute_dh_dx1

end module sll_advection_field
