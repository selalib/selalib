program test_box_splines

#include "sll_working_precision.h"
#include "sll_memory.h"
use hex_mesh
implicit none

type(hex_mesh_2d), pointer  :: mesh
sll_int32                   :: num_cells

num_cells = 40

print *, ""
print *, "Creating a mesh with 40 cells, mesh coordinates written in ./hex_mesh_coo.txt"
mesh => new_hex_mesh_2d(num_cells)
call sll_display(mesh)
call write_hex_mesh_2d(mesh,"hex_mesh_coo.txt")
print *, ""

end program test_box_splines
