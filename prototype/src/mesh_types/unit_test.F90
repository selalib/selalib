program unit_test
#include "sll_mesh_types.h"
#include "sll_working_precision.h"
  use numeric_constants
  implicit none
  
  type(mesh_cylindrical_3D), pointer :: mesh
  integer :: i,j,k
  sll_real64 :: val
  type(mesh_descriptor_1D), pointer :: m1D_descriptor
  type(field_1D_vec1), pointer :: field1D_scalar
  m1D_descriptor => new_mesh_descriptor_1D(0.0_f64, 1.0_f64, 1024)
  field1D_scalar => new_field_1D_vec1(m1D_descriptor)
  print *, 'Allocated a field1d_vec1 type successfully'
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
  print *, 'Proceeding to delete the field_1D_vec1 instance.'
  call delete_field_1D_vec1( field1D_scalar )
  print *, 'Successful, exiting program.'
  
end program unit_test
