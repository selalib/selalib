program test_hex_meshes

#include "sll_working_precision.h"
#include "sll_memory.h"
use sll_hex_meshes
use sll_constants

implicit none

type(sll_hex_mesh_2d), pointer  :: mesh
sll_int32                   :: num_cells
sll_real64, pointer         :: field(:)
sll_int32                   :: error
sll_real64                  :: x1
sll_real64                  :: x2
sll_int32                   :: i
sll_int32                   :: nei1, nei2, nei3

num_cells = 2

print *, ""
print *, "Creating a mesh with 40 cells, mesh coordinates written in ./hex_mesh_coo.txt"
mesh => new_hex_mesh_2d(num_cells)
call sll_display(mesh)
call write_hex_mesh_2d(mesh,"hex_mesh_coo.txt")
call write_hex_mesh_mtv(mesh,"hex_mesh_coo.mtv")
print *, ""

SLL_ALLOCATE(field(mesh%num_pts_tot), error)

do i = 1, mesh%num_pts_tot
   x1 = mesh%global_to_x1(i)
   x2 = mesh%global_to_x2(i)
   field(i) = cos(2*sll_pi*x1)*sin(2*sll_pi*x2)
end do

call get_neighbours(mesh, 1, nei1, nei2, nei3)
print *, "neis =",nei1, nei2, nei3

call write_field_hex_mesh_xmf(mesh, field, 'field')

call write_caid_files(mesh)

call delete(mesh)

end program test_hex_meshes
