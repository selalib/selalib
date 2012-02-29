program unit_test
#include "sll_mesh_types.h"
#include "sll_working_precision.h"
#include "sll_memory.h"

  use numeric_constants
  implicit none
  
  call test_1d()

  call test_1d_product()

  call test_cylindrical_3d()
  print *, 'Successful, exiting program.'

  contains

  subroutine test_1d()

  type(mesh_descriptor_1D), pointer :: m1D_descriptor
  type(field_1D_vec1), pointer :: field1D_scalar

  m1D_descriptor => new_mesh_descriptor_1D(0.0_f64, 1.0_f64, 1024, PERIODIC)

  field1D_scalar => new_field_1D_vec1(m1D_descriptor)
  print *, 'Allocated a field1d_vec1 type successfully'
  
  print *, 'Proceeding to delete the field_1D_vec1 instance.'
  call delete_field_1D_vec1( field1D_scalar )

  end subroutine test_1d

  subroutine test_1d_product()
  type(mesh_descriptor_1d), pointer :: mesh_x
  type(mesh_descriptor_1d), pointer :: mesh_v
  type(mesh_descriptor_2d), pointer :: mesh_xv
  sll_real64 :: x_min, x_max
  sll_int32  :: nc_x
  sll_real64 :: v_min, v_max
  sll_int32  :: nc_v

  x_min =  0.0_f64 ; x_max = 4.0_f64 * sll_pi
  v_min = -6.0_f64 ; v_max = 6.0_f64

  nc_x = 64
  nc_v = 64

  mesh_x => new_mesh_descriptor_1d(x_min,x_max,nc_x, PERIODIC)
  mesh_v => new_mesh_descriptor_1d(v_min,v_max,nc_v, PERIODIC)

  mesh_xv =>  mesh_x * mesh_v

  call mesh_x%dump()

  call write_mesh_2D(mesh_xv)

  end subroutine test_1d_product

  subroutine test_cylindrical_3d()

  type(mesh_cylindrical_3D), pointer :: mesh
  integer :: i,j,k
  sll_real64 :: val

  print *, 'Proceeding to allocate and initialize a 3D mesh...'
  mesh => new_mesh_cylindrical_3D(0.2_f64, 0.8_f64, 0.0_f64, 10.0_f64, &
       32, 32, 256)
  print *, 'Allocation successful'
  
  do i=1, GET_MESH_NCR(mesh)
     do j=1, GET_MESH_NCTHETA(mesh)
        do k=1, GET_MESH_NCR(mesh)
           val = sin(int(i)*2.0*sll_pi/real(GET_MESH_NCR(mesh),f64))*&
                sin(int(j)*2.0*sll_pi/real(GET_MESH_NCTHETA(mesh),f64))*&
                sin(int(k)*2.0*sll_pi/real(GET_MESH_NCZ(mesh),f64))
           SET_MESH_VALUE( mesh, i, j, k, val )
        end do
     end do
  end do
  print *, GET_MESH_DELTA_R(mesh)
  print *, GET_MESH_DELTA_THETA(mesh)
  print *, GET_MESH_DELTA_Z(mesh)
  
  print *, 'Proceeding to delete the mesh...'
  call delete_mesh_cylindrical_3D( mesh )

  end subroutine test_cylindrical_3d

  
end program unit_test
