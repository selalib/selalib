!
!  Contact : Pierre Navaro http://wwww-irma.u-strasbg.fr/~navaro
!
program test_maxwell_fdtd_hex_mesh

#include "sll_working_precision.h"
#include "sll_memory.h"
use hex_mesh
use sll_constants
use sll_maxwell_diga_hex_mesh

implicit none

type(maxwell_dg_hex_mesh)   :: maxwell
type(hex_mesh_2d), pointer  :: mesh
sll_int32                   :: num_cells
sll_int32                   :: error
sll_int32                   :: istep
sll_int32                   :: nstep = 1000
sll_real64                  :: cfl, dt 
sll_real64                  :: time

num_cells = 20

print *, ""
print *, "Creating a mesh with 40 cells, mesh coordinates written in ./hex_mesh_coo.txt"
mesh => new_hex_mesh_2d(num_cells,             &
                        0.0_f64,               &
                        0.0_f64,               &
                        0.5_f64*sqrt(3.0_f64), &
                        0.5_f64,               &
                       -0.5_f64*sqrt(3.0_f64), &
                        0.5_f64,               &
                        0.0_f64,               &
                        1.0_f64,               &
                        10.0_f64 )

call sll_display(mesh)
call write_hex_mesh_2d(mesh, "hex_mesh_coo.txt")
call write_hex_mesh_mtv(mesh, "hex_mesh.mtv")

end program test_maxwell_fdtd_hex_mesh
