module sll_mesh_types
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
  use numeric_constants
  use sll_misc_utils                  ! for is_power_of_2()
  implicit none

  type mesh_descriptor_1D
     sll_real64 :: x1_min
     sll_real64 :: x1_max
     sll_int32  :: ncx1
     sll_real64 :: delta_x1
  end type mesh_descriptor_1D

  type mesh_descriptor_2D
     sll_real64 :: x1_min
     sll_real64 :: x1_max
     sll_real64 :: x2_min
     sll_real64 :: x2_max
     sll_int32  :: ncx1
     sll_int32  :: ncx2
     sll_real64 :: delta_x1
     sll_real64 :: delta_x2
  end type mesh_descriptor_2D

  type mesh_descriptor_3D
     sll_real64 :: x1_min
     sll_real64 :: x1_max
     sll_real64 :: x2_min
     sll_real64 :: x2_max
     sll_real64 :: x3_min
     sll_real64 :: x3_max
     sll_int32  :: ncx1
     sll_int32  :: ncx2
     sll_int32  :: ncx3
     sll_real64 :: delta_x1
     sll_real64 :: delta_x2
     sll_real64 :: delta_x3
  end type mesh_descriptor_3D

  ! It would have been nice to explore the following type, for uniformity
  ! purposes, but presumably, the functions that use this data type, like
  ! a Poisson solver, would operate on simple arrays. Hence the following
  ! type should not be used:
  !
  !  type vec1
  !     sll_real64 :: v1
  !  end type vec1
  !
  ! Similar considerations could prevent the use of the following types.
  ! But it is reasonable to expect that the tuple that describes the 
  ! N-dimensional information at a point would be stored in a structure,
  ! as it would be used together.

  type vec2
     sll_real64 :: v1
     sll_real64 :: v2
  end type vec2

  type vec3
     sll_real64 :: v1
     sll_real64 :: v2
     sll_real64 :: v3
  end type vec3

  type field_1D_vec1
     type(mesh_descriptor_1D), pointer :: descriptor
     sll_real64, dimension(:), pointer :: data ! different from other cases
  end type field_1D_vec1

  type field_1D_vec2
     type(mesh_descriptor_1D), pointer :: descriptor
     type(vec2), dimension(:), pointer :: data
  end type field_1D_vec2

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

contains   ! *****************************************************************

  function new_mesh_descriptor_1D( x1_min, x1_max, ncells )
    intrinsic                         :: real
    type(mesh_descriptor_1D), pointer :: new_mesh_descriptor_1D
    sll_int32                         :: ierr
    sll_real64, intent(in)            :: x1_min
    sll_real64, intent(in)            :: x1_max
    sll_int32, intent(in)             :: ncells
    SLL_ALLOCATE(new_mesh_descriptor_1D, ierr)
    new_mesh_descriptor_1D%x1_min   = x1_min
    new_mesh_descriptor_1D%x1_min   = x1_max
    new_mesh_descriptor_1D%ncx1     = ncells
    new_mesh_descriptor_1D%delta_x1 = (x1_max-x1_min)/real(ncells)
  end function new_mesh_descriptor_1D

  function new_field_1D_vec1( mesh_descriptor )
    type(field_1D_vec1), pointer      :: new_field_1D_vec1
    type(mesh_descriptor_1D), pointer :: mesh_descriptor
    sll_int32                         :: ierr
    SLL_ASSERT(associated(mesh_descriptor))
    SLL_ALLOCATE(new_field_1D_vec1, ierr)
    new_field_1D_vec1%descriptor => mesh_descriptor
    SLL_ALLOCATE(new_field_1D_vec1%data(1:mesh_descriptor%ncx1+1),ierr)
  end function new_field_1D_vec1

  subroutine delete_field_1D_vec1( f1Dv1 )
    type(field_1D_vec1), pointer :: f1Dv1
    sll_int32                    :: ierr
    if( .not. (associated(f1Dv1))) then
       write (*,'(a)') 'ERROR: delete_field_1D_vec1(), not associated argument.'
       STOP
    end if
    nullify(f1Dv1%descriptor)
    SLL_DEALLOCATE(f1Dv1%data, ierr)
    SLL_DEALLOCATE(f1Dv1, ierr)
  end subroutine delete_field_1D_vec1

!#define NEW_FIELD_CONSTRUCTOR_FUNCTION( func_name, descriptor_t, field_t )
!  function func_name( mesh_descriptor )
!    type(descriptor_t)

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
