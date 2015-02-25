program test_tri_poisson
#include "sll_working_precision.h"
                        
use sll_hex_meshes
use tri_poisson
use sll_triangular_meshes
use sll_mesh_calculus_2d_module
use sll_gnuplot

!----------------------------------------------------------------------

implicit none

real(8) :: tcpu
real(8), dimension(:), allocatable :: rho, phi
real(8), dimension(:), allocatable :: ex, ey
type(sll_triangular_poisson_2d)       :: solver
type(sll_triangular_mesh_2d), pointer :: t_mesh
type(sll_hex_mesh_2d), pointer        :: h_mesh

sll_int32  :: num_cells
sll_int32  :: ntypfr(5)
sll_real64 :: potfr(5)

!mesh => new_triangular_mesh_2d("diode.maa") 

num_cells = 2

h_mesh => new_hex_mesh_2d( num_cells, 0._f64, 0._f64) 
  
t_mesh => new_triangular_mesh_2d(h_mesh) 
call write_triangular_mesh_mtv(t_mesh, "test_tri_poisson.mtv")

call analyze_triangular_mesh(t_mesh) 

allocate(ex(t_mesh%num_nodes));  ex  = 0.0_f64 
allocate(ey(t_mesh%num_nodes));  ey  = 0.0_f64
allocate(rho(t_mesh%num_nodes)); rho = 0.0_f64
allocate(phi(t_mesh%num_nodes)); phi = 0.0_f64

call cpu_time(tcpu)  !Initialisation du temps CPU

ntypfr = 1
potfr  = 0.0_f64
call sll_create(solver, t_mesh, ntypfr, potfr)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!   Equation de POISSON - elements finis !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call poissn(solver, ex, ey, rho, phi)
call poliss(solver, phi, ex, ey)
call poifrc(solver, ex, ey)

call sll_gnuplot_2d( phi, "phi", t_mesh%coord, t_mesh%nodes, 1)

end program test_tri_poisson
