!> @internal [example]
program unit_test_triangular_mesh
#include "sll_working_precision.h"
  use sll_triangular_meshes
  implicit none

  type(sll_triangular_mesh_2d) :: mesh
  sll_int32  :: nc_eta1=4
  sll_real64 :: eta1_min=0.0_f64
  sll_real64 :: eta1_max=1.0_f64
  sll_int32  :: nc_eta2=4
  sll_real64 :: eta2_min=0.0_f64
  sll_real64 :: eta2_max=1.0_f64

  call sll_create(mesh,     &
               nc_eta1,  &
               eta1_min, &
               eta1_max, &
               nc_eta2,  &
               eta2_min, &
               eta2_max)

  call sll_display(mesh)
  call write_triangular_mesh_mtv(mesh, "tri_mesh.mtv")
  call sll_delete(mesh)

  print *, 'PASSED'

end program unit_test_triangular_mesh
!> @internal [example]
