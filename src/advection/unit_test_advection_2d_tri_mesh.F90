program unit_test_positions
#include "sll_working_precision.h"
#include "sll_memory.h"
use sll_hex_meshes
use sll_triangular_meshes
use sll_advection_2d_tri_mesh
use sll_gnuplot

implicit none

type(sll_triangular_mesh_2d), pointer :: t_mesh
type(sll_hex_mesh_2d), pointer        :: h_mesh
type(sll_advection_tri_mesh), pointer :: t_adv

sll_real64, dimension(:), allocatable :: df
sll_real64, dimension(:), allocatable :: ex 
sll_real64, dimension(:), allocatable :: ey

sll_real64, dimension(:), pointer :: x1
sll_real64, dimension(:), pointer :: x2

sll_int32 :: num_cells
sll_int32 :: ierr
sll_int32 :: istep

sll_real64 :: dt = 0.01

!Create a triangular mesh from an hex mesh
!Reference on the boundary is set to "one"

num_cells = 50
h_mesh => new_hex_mesh_2d( num_cells, 0._f64, 0._f64) 
t_mesh => new_triangular_mesh_2d(h_mesh) 

SLL_CLEAR_ALLOCATE(df(1:t_mesh%num_nodes), ierr)
SLL_CLEAR_ALLOCATE(ex(1:t_mesh%num_nodes), ierr)
SLL_CLEAR_ALLOCATE(ey(1:t_mesh%num_nodes), ierr)

call map_to_circle(t_mesh, num_cells)

call write_triangular_mesh_mtv(t_mesh, "positions_mesh.mtv")

x1 => t_mesh%coord(1,:)
x2 => t_mesh%coord(2,:)

df = exp(-((x1-0.5)**2+x2*x2)/0.04_f64)
ex = - 1.0_f64 !- x2
ey =   0.0_f64 !+ x1

t_adv => new_advection_2d_tri_mesh(t_mesh)

call sll_gnuplot_2d( df, "f_tri", t_mesh%coord, t_mesh%nodes, 1)
do istep = 1, 50
  call positions(t_adv, df, ex, ey, dt)
end do
call sll_gnuplot_2d( df, "f_tri", t_mesh%coord, t_mesh%nodes, 2)

print*, 'error =', sum(abs(df-exp(-(x1*x1+x2*x2)/0.04_f64)))/t_mesh%num_nodes


call sll_delete(t_mesh)

end program unit_test_positions
