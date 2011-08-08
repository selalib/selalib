module sll_mesh_types
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"

  use numeric_constants
  use sll_misc_utils                  ! for is_power_of_2()
  use sll_diagnostics
  use geometry_functions 
  implicit none

  enum, bind(C)
     enumerator :: PERIODIC = 0, COMPACT = 1
  end enum

  abstract interface
     function scalar_function_2D( eta1, eta2 )
       use sll_working_precision
       sll_real64 :: scalar_function_2D
       sll_real64, intent(in)  :: eta1
       sll_real64, intent(in)  :: eta2
     end function scalar_function_2D
  end interface

  type mesh_descriptor_1D
     sll_real64 :: eta1_min
     sll_real64 :: eta1_max
     sll_int32  :: nc_eta1
     sll_real64 :: delta_eta1
     sll_int32  :: boundary_type
  end type mesh_descriptor_1D

  type geometry_2D
     procedure(scalar_function_2D), pointer, nopass :: x1
     procedure(scalar_function_2D), pointer, nopass :: x2
     procedure(scalar_function_2D), pointer, nopass :: Jacobian11
     procedure(scalar_function_2D), pointer, nopass :: Jacobian12
     procedure(scalar_function_2D), pointer, nopass :: Jacobian21
     procedure(scalar_function_2D), pointer, nopass :: Jacobian22
  end type geometry_2D

  type mesh_descriptor_2D
     sll_real64 :: eta1_min
     sll_real64 :: eta1_max
     sll_real64 :: eta2_min
     sll_real64 :: eta2_max
     sll_int32  :: nc_eta1
     sll_int32  :: nc_eta2
     sll_real64 :: delta_eta1
     sll_real64 :: delta_eta2
     sll_int32  :: boundary1_type
     sll_int32  :: boundary2_type
     type (geometry_2D), pointer :: geom 
  end type mesh_descriptor_2D

  type mesh_descriptor_3D
     sll_real64 :: eta1_min
     sll_real64 :: eta1_max
     sll_real64 :: eta2_min
     sll_real64 :: eta2_max
     sll_real64 :: eta3_min
     sll_real64 :: eta3_max
     sll_int32  :: nc_eta1
     sll_int32  :: nc_eta2
     sll_int32  :: nc_eta3
     sll_real64 :: delta_eta1
     sll_real64 :: delta_eta2
     sll_real64 :: delta_eta3
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

  type field_2D_vec1
     type(mesh_descriptor_2D), pointer :: descriptor
     sll_real64, dimension(:,:), pointer :: data
  end type field_2D_vec1

  type field_2D_vec2
     type(mesh_descriptor_2D), pointer :: descriptor
     type(vec2), dimension(:,:), pointer :: data
  end type field_2D_vec2

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
  

  function new_geometry_2D ( x1, x2, jac11, jac12, jac21, jac22 )
    intrinsic  :: present 
    type(geometry_2D), pointer  ::  new_geometry_2D
    procedure(scalar_function_2D), pointer, optional :: x1
    procedure(scalar_function_2D), pointer, optional :: x2
    procedure(scalar_function_2D), pointer, optional :: jac11
    procedure(scalar_function_2D), pointer, optional :: jac12
    procedure(scalar_function_2D), pointer, optional :: jac21
    procedure(scalar_function_2D), pointer, optional :: jac22

    sll_int32  :: ierr
    
    SLL_ALLOCATE(new_geometry_2D,ierr)
    if (present(x1) .and. present(x2) .and. present(jac11) .and. present(jac12) .and. present(jac21) &
         .and. present(jac22)) then
       new_geometry_2D%x1 => x1
       new_geometry_2D%x2 => x2
       new_geometry_2D%Jacobian11 => jac11
       new_geometry_2D%Jacobian12 => jac12
       new_geometry_2D%Jacobian21 => jac21
       new_geometry_2D%Jacobian22 => jac22
    else
       new_geometry_2D%x1 => default_x1
       new_geometry_2D%x2 => default_x2
       new_geometry_2D%Jacobian11 => default_jac11
       new_geometry_2D%Jacobian12 => default_jac12
       new_geometry_2D%Jacobian21 => default_jac21
       new_geometry_2D%Jacobian22 => default_jac22
    end if
  end function new_geometry_2D

  function new_mesh_descriptor_1D( eta1_min, eta1_max, ncells, bt )
    intrinsic                         :: real
    type(mesh_descriptor_1D), pointer :: new_mesh_descriptor_1D
    sll_int32                         :: ierr
    sll_real64, intent(in)            :: eta1_min
    sll_real64, intent(in)            :: eta1_max
    sll_int32, intent(in)             :: ncells
    sll_int32, intent(in)             :: bt
    SLL_ALLOCATE(new_mesh_descriptor_1D, ierr)
    new_mesh_descriptor_1D%eta1_min   = eta1_min
    new_mesh_descriptor_1D%eta1_max   = eta1_max
    new_mesh_descriptor_1D%nc_eta1     = ncells
    new_mesh_descriptor_1D%delta_eta1 = (eta1_max-eta1_min)/real(ncells)
    new_mesh_descriptor_1D%boundary_type = bt
  end function new_mesh_descriptor_1D

  function new_mesh_descriptor_2D( eta1_min, eta1_max, ncells1, bt1, eta2_min, eta2_max, ncells2, bt2, geom )
    intrinsic                         :: real
    type(mesh_descriptor_2D), pointer :: new_mesh_descriptor_2D
    sll_int32                         :: ierr
    sll_real64, intent(in)            :: eta1_min, eta2_min
    sll_real64, intent(in)            :: eta1_max, eta2_max
    sll_int32, intent(in)             :: ncells1, ncells2
    sll_int32, intent(in)             :: bt1, bt2
    type(geometry_2D), pointer        :: geom

    if (.not.(associated(geom))) then
       stop 'new_mesh_descriptor_2D: pointer geom not associated'
    end if
    SLL_ALLOCATE(new_mesh_descriptor_2D, ierr)
    new_mesh_descriptor_2D%eta1_min   = eta1_min
    new_mesh_descriptor_2D%eta1_max   = eta1_max
    new_mesh_descriptor_2D%nc_eta1     = ncells1
    new_mesh_descriptor_2D%delta_eta1 = (eta1_max-eta1_min)/real(ncells1)
    new_mesh_descriptor_2D%boundary1_type = bt1
    new_mesh_descriptor_2D%eta2_min   = eta2_min
    new_mesh_descriptor_2D%eta2_max   = eta2_max
    new_mesh_descriptor_2D%nc_eta2     = ncells2
    new_mesh_descriptor_2D%delta_eta2 = (eta2_max-eta2_min)/real(ncells2)
    new_mesh_descriptor_2D%boundary2_type = bt2
    new_mesh_descriptor_2D%geom => geom
  end function new_mesh_descriptor_2D


  function new_field_1D_vec1( mesh_descriptor )
    type(field_1D_vec1), pointer      :: new_field_1D_vec1
    type(mesh_descriptor_1D), pointer :: mesh_descriptor
    sll_int32                         :: ierr
    SLL_ASSERT(associated(mesh_descriptor))
    SLL_ALLOCATE(new_field_1D_vec1, ierr)
    new_field_1D_vec1%descriptor => mesh_descriptor
    SLL_ALLOCATE(new_field_1D_vec1%data(1:mesh_descriptor%nc_eta1+1),ierr)
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

  function new_field_2D_vec1( mesh_descriptor )
    type(field_2D_vec1), pointer      :: new_field_2D_vec1
    type(mesh_descriptor_2D), pointer :: mesh_descriptor
    sll_int32                         :: ierr
    SLL_ASSERT(associated(mesh_descriptor))
    SLL_ALLOCATE(new_field_2D_vec1, ierr)
    new_field_2D_vec1%descriptor => mesh_descriptor
    SLL_ALLOCATE(new_field_2D_vec1%data(1:mesh_descriptor%nc_eta1+1,1:mesh_descriptor%nc_eta2+1),ierr)
  end function new_field_2D_vec1

  subroutine delete_field_2D_vec1( f2Dv1 )
    type(field_2D_vec1), pointer :: f2Dv1
    sll_int32                    :: ierr
    if( .not. (associated(f2Dv1))) then
       write (*,'(a)') 'ERROR: delete_field_1D_vec1(), not associated argument.'
       STOP
    end if
    nullify(f2Dv1%descriptor)
    SLL_DEALLOCATE(f2Dv1%data, ierr)
    SLL_DEALLOCATE(f2Dv1, ierr)
  end subroutine delete_field_2D_vec1

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

  subroutine write_mesh_2D (mesh)
    type(mesh_descriptor_2D), pointer :: mesh
    sll_real64, dimension(:,:), pointer :: x1mesh
    sll_real64, dimension(:,:), pointer :: x2mesh
    sll_int32  :: i1
    sll_int32  :: i2
    sll_real64 :: eta1
    sll_real64 :: eta2
    sll_int32 :: ierr

    ! create 2D mesh
    SLL_ALLOCATE(x1mesh(mesh%nc_eta1+1,mesh%nc_eta2+1), ierr)
    SLL_ALLOCATE(x2mesh(mesh%nc_eta1+1,mesh%nc_eta2+1), ierr)
    eta1 = mesh%eta1_min
    do i1=1, mesh%nc_eta1+1
       eta2 = mesh%eta2_min
       do i2=1, mesh%nc_eta2+1
          x1mesh(i1,i2) = mesh%geom%x1(eta1,eta2)
          x2mesh(i1,i2) = mesh%geom%x2(eta1,eta2)
          eta2 = eta2 + mesh%delta_eta2 
       end do
       eta1 = eta1 + mesh%delta_eta1
    end do
    ! missing the actual function to write the mesh in hdf5. Should come from diagnostics
  end subroutine write_mesh_2D

  ! writes field along with mesh for the moment. The mesh part can be removed when it can be written separately
  subroutine write_field_2D_vec1( f2Dv1 )
    type(field_2D_vec1), pointer :: f2Dv1

    type(mesh_descriptor_2D), pointer :: mesh
    sll_real64, dimension(:,:), pointer :: x1mesh
    sll_real64, dimension(:,:), pointer :: x2mesh
    sll_int32  :: i1
    sll_int32  :: i2
    sll_real64 :: eta1
    sll_real64 :: eta2
    sll_int32 :: ierr
    
    ! create 2D mesh
    mesh => f2Dv1%descriptor
    SLL_ASSERT(associated(mesh))
    SLL_ALLOCATE(x1mesh(mesh%nc_eta1+1,mesh%nc_eta2+1), ierr)
    SLL_ALLOCATE(x2mesh(mesh%nc_eta1+1,mesh%nc_eta2+1), ierr)
    eta1 = mesh%eta1_min
    do i1=1, mesh%nc_eta1+1
       eta2 = mesh%eta2_min
       do i2=1, mesh%nc_eta2+1
          x1mesh(i1,i2) = mesh%geom%x1(eta1,eta2)
          x2mesh(i1,i2) = mesh%geom%x2(eta1,eta2)
          eta2 = eta2 + mesh%delta_eta2 
       end do
       eta1 = eta1 + mesh%delta_eta1
    end do
    call write_fxv(f2Dv1%data,x1mesh,x2mesh)
  end subroutine write_field_2D_vec1
end module sll_mesh_types
