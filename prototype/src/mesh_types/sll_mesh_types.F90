module sll_mesh_types
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
  use numeric_constants
  use sll_misc_utils                  ! for is_power_of_2()
  implicit none

  type mesh_cylindrical_3D
     sll_real64 :: rmin
     sll_real64 :: rmax
     sll_real64 :: zmin
     sll_real64 :: zmax
     sll_int32  :: ncr
     sll_int32  :: nctheta
     sll_int32  :: ncz
     sll_real64 :: delta_r
     sll_real64 :: delta_theta
     sll_real64 :: delta_z
     sll_real64, dimension(:,:,:), pointer :: data
  end type mesh_cylindrical_3D

contains

  function new_mesh_cylindrical_3D(rmin, rmax, zmin, zmax, ncr, nctheta, ncz )
    intrinsic :: real, int
    type(mesh_cylindrical_3D), pointer :: new_mesh_cylindrical_3D
    sll_real64, intent(in) :: rmin
    sll_real64, intent(in) :: rmax
    sll_real64, intent(in) :: zmin
    sll_real64, intent(in) :: zmax
    sll_int32,  intent(in) :: ncr
    sll_int32,  intent(in) :: nctheta
    sll_int32,  intent(in) :: ncz
    sll_int32              :: ierr

    SLL_ASSERT(rmax .ge. rmin)
    SLL_ASSERT(zmax .ge. zmin)
    SLL_ASSERT(is_power_of_2(int(ncr,i64)))
    SLL_ASSERT(is_power_of_2(int(nctheta,i64)))
    ! The following permits to have a number of cells in z equal to 0, for
    ! 2D tests. Something similar could be done for the other directions...
    SLL_ASSERT((is_power_of_2(int(ncz,i64))) .or. (ncz .eq. 0))

    SLL_ALLOCATE(new_mesh_cylindrical_3D, ierr)
    new_mesh_cylindrical_3D%rmin        = rmin
    new_mesh_cylindrical_3D%rmax        = rmax
    new_mesh_cylindrical_3D%zmin        = zmin
    new_mesh_cylindrical_3D%zmax        = zmax
    new_mesh_cylindrical_3D%ncr         = ncr
    new_mesh_cylindrical_3D%nctheta     = nctheta
    new_mesh_cylindrical_3D%ncz         = ncz
    new_mesh_cylindrical_3D%delta_r     = (rmax-rmin)/real(ncr,f64)
    new_mesh_cylindrical_3D%delta_theta = 2.0_f64*sll_pi/real(nctheta,f64)
    new_mesh_cylindrical_3D%delta_z     = (zmax-zmin)/real(ncz,f64)
    SLL_ALLOCATE(new_mesh_cylindrical_3D%data(1:ncr+1,1:nctheta+1,1:ncz+1),ierr)
  end function new_mesh_cylindrical_3D

  subroutine delete_mesh_cylindrical_3D( mesh )
    intrinsic :: associated
    type(mesh_cylindrical_3D), pointer :: mesh
    sll_int32 :: ierr

    if(.not. associated(mesh)) then
       write (*,'(a)') 'ERROR: non-associated mesh pointer passed to DELETE'
       STOP 'delete_mesh_cylindrical_3D'
    end if
    SLL_DEALLOCATE( mesh%data, ierr )
    SLL_DEALLOCATE( mesh, ierr )
  end subroutine delete_mesh_cylindrical_3D

end module sll_mesh_types
