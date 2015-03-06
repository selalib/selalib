program test_tri_poisson
#include "sll_working_precision.h"
#include "sll_memory.h"
                        
use sll_hex_meshes
use tri_poisson
use sll_triangular_meshes
use sll_mesh_calculus_2d_module
use sll_gnuplot

!----------------------------------------------------------------------

implicit none

real(8) :: tcpu
real(8), dimension(:), allocatable    :: x1
real(8), dimension(:), allocatable    :: x2
real(8), dimension(:), allocatable    :: rho
real(8), dimension(:), allocatable    :: phi
real(8), dimension(:), allocatable    :: e_x
real(8), dimension(:), allocatable    :: e_y
type(sll_triangular_poisson_2d)       :: solver
type(sll_triangular_mesh_2d), pointer :: t_mesh
type(sll_hex_mesh_2d), pointer        :: h_mesh

sll_int32  :: num_cells
sll_int32  :: ntypfr(5)
sll_real64 :: potfr(5)

sll_int32  :: nc_x1 = 64
sll_real64 :: x1_min = -1.0_f64, x1_max = 1.0_f64
sll_int32  :: nc_x2 = 64
sll_real64 :: x2_min = -1.0_f64, x2_max = 1.0_f64

sll_int32  :: i, ierr

!mesh => new_triangular_mesh_2d("diode.maa") 

num_cells = 10

h_mesh => new_hex_mesh_2d( num_cells, 0._f64, 0._f64) 
  
!t_mesh => new_triangular_mesh_2d(h_mesh) 

t_mesh => new_triangular_mesh_2d(nc_x1, x1_min, x1_max, nc_x2, x2_min, x2_max) 

call analyze_triangular_mesh(t_mesh) 
call write_triangular_mesh_mtv(t_mesh, "test_tri_poisson.mtv")

SLL_CLEAR_ALLOCATE(x1(1:t_mesh%num_nodes),ierr)
SLL_CLEAR_ALLOCATE(x2(1:t_mesh%num_nodes),ierr)

SLL_CLEAR_ALLOCATE(e_x(1:t_mesh%num_nodes),ierr)
SLL_CLEAR_ALLOCATE(e_y(1:t_mesh%num_nodes),ierr)
SLL_CLEAR_ALLOCATE(rho(1:t_mesh%num_nodes),ierr)
SLL_CLEAR_ALLOCATE(phi(1:t_mesh%num_nodes),ierr)

call cpu_time(tcpu)  !Initialisation du temps CPU

ntypfr(1) = 1; potfr(1) = 0.0_f64
ntypfr(2) = 1; potfr(2) = 0.0_f64
ntypfr(3) = 1; potfr(3) = 0.0_f64
ntypfr(4) = 1; potfr(4) = 0.0_f64
ntypfr(5) = 1; potfr(5) = 0.0_f64

call sll_create(solver, t_mesh, ntypfr, potfr)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!   Equation de POISSON - elements finis !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

x1 = t_mesh%coord(1,:)
x2 = t_mesh%coord(2,:)

rho = 8*sll_pi*sll_pi*sin(2*sll_pi*x1)*sin(2*sll_pi*x2)

call sll_compute_phi_from_rho(solver, rho, phi)
call sll_gnuplot_2d( phi, "phi", t_mesh%coord, t_mesh%nodes, 1)
!call sll_compute_e_from_phi(solver, phi, e_x, e_y)
!call sll_gnuplot_2d( e_x, "e_x", t_mesh%coord, t_mesh%nodes, 1)
!call sll_gnuplot_2d( e_y, "e_y", t_mesh%coord, t_mesh%nodes, 1)

print*,'error=', maxval(abs(phi-sin(2*sll_pi*x1)*sin(2*sll_pi*x2)))

end program test_tri_poisson
