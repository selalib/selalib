program unit_test_positions
#include "sll_working_precision.h"
#include "sll_memory.h"
use sll_hex_meshes
use sll_triangular_meshes
use sll_advection_2d_tri_mesh

implicit none

type(sll_triangular_mesh_2d), pointer :: t_mesh
type(sll_hex_mesh_2d), pointer        :: h_mesh
type(sll_advection_tri_mesh), pointer :: t_adv

sll_real64, dimension(:), allocatable :: df
sll_real64, dimension(:), allocatable :: ex 
sll_real64, dimension(:), allocatable :: ey

sll_int32 :: num_cells
sll_int32 :: ierr

sll_real64 :: dt = 0.1

!Create a triangular mesh from an hex mesh
!Reference on the boundary is set to "one"

num_cells = 2
h_mesh => new_hex_mesh_2d( num_cells, 0._f64, 0._f64) 
t_mesh => new_triangular_mesh_2d(h_mesh) 

SLL_CLEAR_ALLOCATE(df(1:t_mesh%num_nodes), ierr)
SLL_CLEAR_ALLOCATE(ex(1:t_mesh%num_nodes), ierr)
SLL_CLEAR_ALLOCATE(ey(1:t_mesh%num_nodes), ierr)

call map_to_circle(t_mesh, num_cells)

call write_triangular_mesh_mtv(t_mesh, "positions_mesh.mtv")

ex = 1.0_f64

t_adv => new_advection_2d_tri_mesh(t_mesh)
call positions(t_adv, df, ex, ey, dt)

call sll_delete(t_mesh)

end program unit_test_positions
