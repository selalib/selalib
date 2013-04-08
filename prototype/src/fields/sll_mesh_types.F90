!------------------------------------------------------------------------------
! SELALIB
!------------------------------------------------------------------------------
!
! MODULE: sll_mesh_types
!
!> @author
!> - Edwin
!> - Pierre
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
!>    - diagnostics
!
! REVISION HISTORY:
! DD Mmm YYYY - Initial Version
! TODO_dd_mmm_yyyy - TODO_describe_appropriate_changes - TODO_name
!------------------------------------------------------------------------------
module sll_mesh_types
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"

  use sll_constants
  use sll_utilities                  ! for is_power_of_two()
  use sll_diagnostics
  use geometry_functions 
  implicit none

  integer, parameter :: PERIODIC = 0, COMPACT = 1
  !enum, bind(C)
  !   enumerator :: PERIODIC = 0, COMPACT = 1
  !end enum

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
     procedure(scalar_function_2D), pointer, nopass :: eta1
     procedure(scalar_function_2D), pointer, nopass :: eta2
     procedure(scalar_function_2D), pointer, nopass :: Jacobian11
     procedure(scalar_function_2D), pointer, nopass :: Jacobian12
     procedure(scalar_function_2D), pointer, nopass :: Jacobian21
     procedure(scalar_function_2D), pointer, nopass :: Jacobian22
     procedure(scalar_function_2D), pointer, nopass :: Jacobian
     sll_real64, dimension(:,:), pointer :: x1_at_i 
     sll_real64, dimension(:,:), pointer :: x2_at_i
     sll_real64, dimension(:,:), pointer :: x1c_at_i 
     sll_real64, dimension(:,:), pointer :: x2c_at_i
     sll_real64, dimension(:,:), pointer :: eta1_at_i
     sll_real64, dimension(:,:), pointer :: eta2_at_i
     sll_real64, dimension(:,:), pointer :: Jacobian11_at_i
     sll_real64, dimension(:,:), pointer :: Jacobian12_at_i
     sll_real64, dimension(:,:), pointer :: Jacobian21_at_i
     sll_real64, dimension(:,:), pointer :: Jacobian22_at_i
     sll_real64, dimension(:,:), pointer :: Jacobian_at_i
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
     !type (geometry_2D), pointer :: geom
   contains
     procedure(geom_func_int), pass :: x1_int
     procedure(geom_func_int), pass :: x2_int
     procedure(geom_func_int), pass :: jac_int
     procedure(geom_func_real), pass :: x1_real
     procedure(geom_func_real), pass :: x2_real
     procedure(geom_func_real), pass :: jac_real
     generic :: x1 => x1_int, x1_real
     generic :: x2 => x2_int, x2_real
     generic :: Jacobian => jac_int, jac_real
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

  !This type is defined only in splitted meshes case
  type field_4D_vec1
     !type(mesh_descriptor_4D), pointer :: descriptor
     type(mesh_descriptor_2D), pointer :: descriptor_x
     type(mesh_descriptor_2D), pointer :: descriptor_v
     sll_real64, dimension(:,:,:,:), pointer :: data
  end type field_4D_vec1

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
  

  function new_geometry_2D ( name )
    type(geometry_2D), pointer  ::  new_geometry_2D
    character(len=*)               :: name
    sll_int32 :: nc1
    sll_int32 :: nc2
    sll_real64, dimension(:,:), optional :: x1_in
    sll_real64, dimension(:,:), optional :: x2_in
    sll_real64, dimension(:,:), optional :: x1c_in
    sll_real64, dimension(:,:), optional :: x2c_in
    sll_real64, dimension(:,:), optional :: jac_in
    sll_int32, optional :: bt1
    sll_int32, optional :: bt2
    
    ! local variables
    sll_real64 :: eta1, eta2, delta1, delta2
    sll_int32  :: i1, i2, ierr
    
    SLL_ALLOCATE(new_geometry_2D,ierr)
    SLL_ALLOCATE(new_geometry_2D%x1_at_i(nc1+1,nc2+1),ierr)
    SLL_ALLOCATE(new_geometry_2D%x2_at_i(nc1+1,nc2+1),ierr)
    SLL_ALLOCATE(new_geometry_2D%x1c_at_i(nc1+1,nc2+1),ierr)
    SLL_ALLOCATE(new_geometry_2D%x2c_at_i(nc1+1,nc2+1),ierr)
    SLL_ALLOCATE(new_geometry_2D%Jacobian_at_i(nc1+1,nc2+1),ierr)

    ! cartesian coordinates correspond to identity mapping
    if ((name(1:8)=='identity').or.(name(1:9)=='cartesian')) then
       new_geometry_2D%x1         => identity_x1
       new_geometry_2D%x2         => identity_x2
       new_geometry_2D%eta1       => identity_eta1
       new_geometry_2D%eta2       => identity_eta2
       new_geometry_2D%Jacobian11 => identity_jac11
       new_geometry_2D%Jacobian12 => identity_jac12
       new_geometry_2D%Jacobian21 => identity_jac21
       new_geometry_2D%Jacobian22 => identity_jac22
       new_geometry_2D%Jacobian => identity_jac
    ! polar coordinates
    else if (name(1:5) == 'polar') then
       new_geometry_2D%x1         => polar_x1
       new_geometry_2D%x2         => polar_x2
       new_geometry_2D%eta1       => polar_eta1
       new_geometry_2D%eta2       => polar_eta2
       new_geometry_2D%Jacobian11 => polar_jac11
       new_geometry_2D%Jacobian12 => polar_jac12
       new_geometry_2D%Jacobian21 => polar_jac21
       new_geometry_2D%Jacobian22 => polar_jac22
       new_geometry_2D%Jacobian => polar_jac
    else if (name(1:7) == 'sinprod') then
       new_geometry_2D%x1         => sinprod_x1
       new_geometry_2D%x2         => sinprod_x2
       new_geometry_2D%eta1       => sinprod_eta1
       new_geometry_2D%eta2       => sinprod_eta2
       new_geometry_2D%Jacobian11 => sinprod_jac11
       new_geometry_2D%Jacobian12 => sinprod_jac12
       new_geometry_2D%Jacobian21 => sinprod_jac21
       new_geometry_2D%Jacobian22 => sinprod_jac22
       new_geometry_2D%Jacobian => sinprod_jac
    else if (name(1:4) == 'test') then
       new_geometry_2D%x1         => test_x1
       new_geometry_2D%x2         => test_x2
       new_geometry_2D%eta1       => test_eta1
       new_geometry_2D%eta2       => test_eta2
       new_geometry_2D%Jacobian11 => test_jac11
       new_geometry_2D%Jacobian12 => test_jac12
       new_geometry_2D%Jacobian21 => test_jac21
       new_geometry_2D%Jacobian22 => test_jac22
       new_geometry_2D%Jacobian => test_jac
    else if (name(1:10) == 'from_array') then
       if (present(x1_in) .and. &
            present(x2_in) .and. present(x1c_in) .and. &
            present(x2c_in) .and. present( jac_in)) then

          new_geometry_2D%x1_at_i = x1_in
          new_geometry_2D%x2_at_i = x2_in
          new_geometry_2D%x1c_at_i = x1c_in
          new_geometry_2D%x2c_at_i = x2c_in
          new_geometry_2D%Jacobian_at_i = jac_in
          ! initialize cubic spline mesh interpolation
          call init_cubic_spline(x1_in, x2_in, bt1, bt2)

          new_geometry_2D%x1         => cubic_spline_x1
          new_geometry_2D%x2         => cubic_spline_x2
          new_geometry_2D%eta1       => cubic_spline_eta1
          new_geometry_2D%eta2       => cubic_spline_eta2
          new_geometry_2D%Jacobian11 => cubic_spline_jac11
          new_geometry_2D%Jacobian12 => cubic_spline_jac12
          new_geometry_2D%Jacobian21 => cubic_spline_jac21
          new_geometry_2D%Jacobian22 => cubic_spline_jac22
          new_geometry_2D%Jacobian => cubic_spline_jac
       else 
          print*, 'new_geometry_2D: from_array. Need input'
          stop
       end if
       if (.not.(name(1:10) == 'from_array')) then
          delta1 = 1.0_f64 / nc1
          delta1 = 1.0_f64 / nc2
          eta1 = 0.0_f64
          eta2 = 0.0_f64
          do i2 = 1, nc2
             do i1 = 1, nc1
                new_geometry_2D%x1_at_i = new_geometry_2D%x1(eta1,eta2)
                new_geometry_2D%x2_at_i = new_geometry_2D%x2(eta1,eta2)
                new_geometry_2D%x1c_at_i = new_geometry_2D%x1(eta1+0.5_f64*delta1,eta2+0.5_f64*delta2)
                new_geometry_2D%x2c_at_i = new_geometry_2D%x2(eta1+0.5_f64*delta1,eta2+0.5_f64*delta2)
                new_geometry_2D%Jacobian_at_i = new_geometry_2D%Jacobian(eta1+0.5_f64*delta1,eta2+0.5_f64*delta2)
                eta1 = eta1 + delta1
             end do
             eta2 = eta2 + delta2
          end do
       end if
    else
       print*, 'new_geometry_2D: mapping ',name, ' is not implemented'
       stop
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

  function new_field_4D_vec1( mesh_descriptor_x, mesh_descriptor_v )
    type(field_4D_vec1), pointer      :: new_field_4D_vec1
    type(mesh_descriptor_2D), pointer :: mesh_descriptor_x
    type(mesh_descriptor_2D), pointer :: mesh_descriptor_v
    sll_int32                         :: ierr
    SLL_ASSERT(associated(mesh_descriptor_x))
    SLL_ASSERT(associated(mesh_descriptor_v))
    SLL_ALLOCATE(new_field_4D_vec1, ierr)
    new_field_4D_vec1%descriptor_x => mesh_descriptor_x
    new_field_4D_vec1%descriptor_v => mesh_descriptor_v
    SLL_ALLOCATE(new_field_4D_vec1%data(mesh_descriptor_x%nc_eta1+1, mesh_descriptor_x%nc_eta2+1,mesh_descriptor_v%nc_eta1+1, mesh_descriptor_v%nc_eta2+1),ierr)
  end function new_field_4D_vec1

  subroutine delete_field_4D_vec1( f4Dv1 )
    type(field_4D_vec1), pointer :: f4Dv1
    sll_int32                    :: ierr
    if( .not. (associated(f4Dv1))) then
       write (*,'(a)') 'ERROR: delete_field_4D_vec1(), not associated argument.'
       STOP
    end if
    nullify(f4Dv1%descriptor_x)
    nullify(f4Dv1%descriptor_v)
    SLL_DEALLOCATE(f4Dv1%data, ierr)
    SLL_DEALLOCATE(f4Dv1, ierr)
  end subroutine delete_field_4D_vec1

  function new_field_2D_vec2( mesh_descriptor )
    type(field_2D_vec2), pointer      :: new_field_2D_vec2
    type(mesh_descriptor_2D), pointer :: mesh_descriptor
    sll_int32                         :: ierr
    SLL_ASSERT(associated(mesh_descriptor))
    SLL_ALLOCATE(new_field_2D_vec2, ierr)
    new_field_2D_vec2%descriptor => mesh_descriptor
    SLL_ALLOCATE(new_field_2D_vec2%data(1:mesh_descriptor%nc_eta1+1,1:mesh_descriptor%nc_eta2+1),ierr)
  end function new_field_2D_vec2

  subroutine delete_field_2D_vec2( f2Dv2 )
    type(field_2D_vec2), pointer :: f2Dv2
    sll_int32                    :: ierr
    if( .not. (associated(f2Dv2))) then
       write (*,'(a)') 'ERROR: delete_field_1D_vec1(), not associated argument.'
       STOP
    end if
    nullify(f2Dv2%descriptor)
    SLL_DEALLOCATE(f2Dv2%data, ierr)
    SLL_DEALLOCATE(f2Dv2, ierr)
  end subroutine delete_field_2D_vec2


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
    SLL_ASSERT(is_power_of_two(int(ncr,i64)))
    SLL_ASSERT(is_power_of_two(int(nctheta,i64)))
    ! The following permits to have a number of cells in z equal to 0, for
    ! 2D tests. Something similar could be done for the other directions...
    SLL_ASSERT((is_power_of_two(int(ncz,i64))) .or. (ncz .eq. 0))

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

  subroutine write_mesh_2D(mesh,mesh_name)
    type(mesh_descriptor_2D), pointer :: mesh
    character(len=*), intent(in), optional :: mesh_name
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
          x1mesh(i1,i2) = mesh%geom%x1_at_i(i1,i2)
          x2mesh(i1,i2) = mesh%geom%x2_at_i(i1,i2)
          eta2 = eta2 + mesh%delta_eta2 
       end do
       eta1 = eta1 + mesh%delta_eta1
    end do
    if (present(mesh_name)) then
       call write_mesh(x1mesh,x2mesh,mesh%nc_eta1+1,mesh%nc_eta2+1,trim(mesh_name))
    else
       call write_mesh(x1mesh,x2mesh,mesh%nc_eta1+1,mesh%nc_eta2+1,"mesh")
    end if
  end subroutine write_mesh_2D

  subroutine write_field_2D_vec1( f2Dv1, name, jacobian, average, center )
    type(field_2D_vec1), pointer :: f2Dv1
    character(len=*) :: name
    logical, optional       :: jacobian   ! .true. if field data multiplied by jacobian is stored
    sll_real64, optional    :: average    ! average value to add to field
    sll_int32, intent(in), optional     :: center     ! node or cell centered values

    type(mesh_descriptor_2D), pointer :: mesh
    sll_real64, dimension(:,:), pointer :: x1mesh
    sll_real64, dimension(:,:), pointer :: x2mesh
    sll_int32  :: i1
    sll_int32  :: i2
    sll_real64 :: eta1
    sll_real64 :: eta2
    sll_real64 :: avg
    sll_int32 :: ierr

    sll_real64, dimension(:,:), pointer :: val

    ! create 2D mesh
    mesh => f2Dv1%descriptor
    SLL_ASSERT(associated(mesh))

    SLL_ALLOCATE(val(mesh%nc_eta1+1,mesh%nc_eta2 + 1), ierr)

    if (present(average)) then
       avg = average
    else
       avg = 0.0_f64
    end if

    if (.not.(present(jacobian))) then
       call write_vec1d(f2Dv1%data,mesh%nc_eta1+1,mesh%nc_eta2+1,name,"mesh",center)
    else if (jacobian) then 
       ! quantity multiplied by Jacobian is stored, need to divide by jacobian for
       eta1 = mesh%eta1_min + 0.5_f64 * mesh%delta_eta1
       do i1 = 1, mesh%nc_eta1
          eta2 = mesh%eta2_min + 0.5_f64 * mesh%delta_eta2
          do i2 = 1, mesh%nc_eta2
!             SLL_ASSERT( mesh%geom%Jacobian (eta1, eta2) > 0.0_f64 )
!             if (mesh%geom%Jacobian (eta1, eta2) > 1.0D-14) then 
!                val(i1,i2) = f2Dv1%data( i1,i2) / mesh%geom%Jacobian (eta1, eta2) + avg
!             else 
!                val(i1,i2) = 0.0
!             end if
             val(i1,i2) = f2Dv1%data( i1,i2) / mesh%geom%Jacobian_at_i (i1, i2) + avg
             eta2 = eta2 + mesh%delta_eta2
          end do
          eta1 = eta1 + mesh%delta_eta1
       end do
       call write_vec1d(val,mesh%nc_eta1+1,mesh%nc_eta2+1,name,"mesh",center)
    end if
  end subroutine write_field_2D_vec1
end module sll_mesh_types
