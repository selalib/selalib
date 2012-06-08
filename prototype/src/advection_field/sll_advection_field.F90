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

    nc_eta1 = GET_FIELD_NC_ETA1( this ) 
    nc_eta2 = GET_FIELD_NC_ETA2( this )     
    
    if (present(phi_external)) then 
       do i1 = 1, nc_eta1+1
          do i2 = 1, nc_eta2+1 
             this%data(i1,i2) = 0.5_f64*this%pmass * this%mesh%x2_at_node(i1,i2)**2 &
                  + this%pcharge*(FIELD_1D_AT_I(phi_self,i1)  &
                  + FIELD_1D_AT_I(phi_external,i1))
          end do
       end do
    else 
       do i1 = 1, nc_eta1
          do i2 = 1, nc_eta2
             this%data(i1,i2) = 0.5_f64*this%pmass * this%mesh%x2_at_node(i1,i2)**2 &
                  + this%pcharge*FIELD_1D_AT_I(phi_self,i1)
          end do
       end do
    end if
  end subroutine compute_hamiltonian
    
end module sll_advection_field
