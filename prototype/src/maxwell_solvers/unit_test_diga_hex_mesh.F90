program test_box_splines

#include "sll_working_precision.h"
#include "sll_memory.h"
use hex_mesh
use sll_constants
use sll_maxwell_diga_hex_mesh

implicit none

type(maxwell_dg_hex_mesh)   :: maxwell
type(hex_mesh_2d), pointer  :: mesh
sll_int32                   :: num_cells
sll_real64, pointer         :: field(:)
sll_int32                   :: error
sll_real64                  :: x1
sll_real64                  :: x2
sll_int32                   :: i
sll_int32                   :: degree

num_cells = 40

print *, ""
print *, "Creating a mesh with 40 cells, mesh coordinates written in ./hex_mesh_coo.txt"
mesh => new_hex_mesh_2d(num_cells)
call sll_display(mesh)
call write_hex_mesh_2d(mesh,"hex_mesh_coo.txt")
print *, ""

degree = 1
call initialize(maxwell, mesh, degree)

end program test_box_splines

