!------------------------------------------------------------------------------
! SELALIB
!------------------------------------------------------------------------------
!
! MODULE: sll_maxwell_2d
!
!> @author
!> Pierre Navaro Philippe Helluy
!>
!
! DESCRIPTION: 
!
!> @brief
!> Implements the Maxwell solver in 2D
!>
!>@details
!>This module depends on:
!> - memory
!> - precision
!> - assert 
!> - numerical_utilities
!> - constants
!> - mesh_types
!> - diagnostics
!> - sll_utilities
!>
! REVISION HISTORY:
! 03 02 2012 - Initial Version  (fevrier)
! TODO_dd_mmm_yyyy - TODO_describe_appropriate_changes - TODO_name
!------------------------------------------------------------------------------

module sll_maxwell_2d

#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_mesh_types.h"

use numeric_constants

implicit none
private

public :: new_maxwell_2d
public :: solve_maxwell_2d
public :: delete_maxwell_2d

!> Object with data to solve Maxwell equation on 2d domain
!> Maxwell in TE mode: (Ex,Ey,Hz)
type, public :: maxwell_2d
  type(field_2d_vec3), pointer        :: TEfield,dt_TEfield  
  type(mesh_descriptor_2d), pointer   :: descriptor
end type maxwell_2d

contains

function new_maxwell_2d(TEfield)

   type(maxwell_2d),pointer :: new_maxwell_2d
   type(field_2D_vec3),      pointer :: TEfield
   type(mesh_descriptor_2d), pointer :: mesh
   sll_int32                         :: error
   sll_int32                         :: ncx
   sll_int32                         :: ncy

   mesh => TEfield%descriptor
   ncx = GET_MESH_NC_ETA1(mesh)
   ncy = GET_MESH_NC_ETA2(mesh)

   SLL_ALLOCATE(new_maxwell_2d,                   error)
   new_maxwell_2d%dt_TEfield => new_field_2D_vec3(mesh)
!>  SLL_ALLOCATE(new_maxwell_2d%descriptor,        error)  ! useless ??

   new_maxwell_2d%descriptor => mesh


 end function new_maxwell_2d

subroutine solve_maxwell_2d(this,error)

   type(maxwell_2d), intent(inout) :: this
   sll_int32 , intent(out)                  :: error
   sll_int32                                :: ncx,ncy
   sll_int32                                :: i, j


   ncx = GET_FIELD_NC_ETA1(this)
   ncy = GET_FIELD_NC_ETA2(this)

 end subroutine solve_maxwell_2d

!> Delete the Maxwell object
subroutine delete_maxwell_2d(this)
  
  type(maxwell_2d) :: this
  !deallocate(this%rhst)

end subroutine delete_maxwell_2d

end module sll_maxwell_2d
