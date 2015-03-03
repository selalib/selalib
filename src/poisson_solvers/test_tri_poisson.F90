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

sll_int32  :: nc_x1 = 10
sll_real64 :: x1_min = 0.0_f64, x1_max = 1.0_f64
sll_int32  :: nc_x2 = 10
sll_real64 :: x2_min = 0.0_f64, x2_max = 1.0_f64

sll_int32  :: i
sll_real64 :: x1, x2

!mesh => new_triangular_mesh_2d("diode.maa") 

num_cells = 2

h_mesh => new_hex_mesh_2d( num_cells, 0._f64, 0._f64) 
  
t_mesh => new_triangular_mesh_2d(h_mesh) 

!t_mesh => new_triangular_mesh_2d(nc_x1, x1_min, x1_max, nc_x2, x2_min, x2_max) 

call write_triangular_mesh_mtv(t_mesh, "test_tri_poisson.mtv")

call analyze_triangular_mesh(t_mesh) 

allocate(e_x(t_mesh%num_nodes)); e_x = 0.0_f64 
allocate(e_y(t_mesh%num_nodes)); e_y = 0.0_f64
allocate(rho(t_mesh%num_nodes)); rho = 0.0_f64
allocate(phi(t_mesh%num_nodes)); phi = 0.0_f64

call cpu_time(tcpu)  !Initialisation du temps CPU

ntypfr(1) = 1; potfr(1) = 0.0_f64
ntypfr(2) = 1; potfr(2) = 0.0_f64
ntypfr(3) = 1; potfr(3) = 0.0_f64
ntypfr(4) = 1; potfr(4) = 0.0_f64
call sll_create(solver, t_mesh, ntypfr, potfr)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!   Equation de POISSON - elements finis !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


do i = 1, h_mesh%num_pts_tot
   x1 = h_mesh%global_to_x1(i)
   x2 = h_mesh%global_to_x2(i)
   rho(i) = exp(-(x1*x1+x2*x2)/0.01)
end do

call write_field_hex_mesh_xmf(h_mesh, rho, 'rho')

call poissn(solver, e_x, e_y, rho, phi)
call poliss(solver, phi, e_x, e_y)
call poifrc(solver, e_x, e_y)

call sll_gnuplot_2d( phi, "phi", t_mesh%coord, t_mesh%nodes, 1)
call write_field_hex_mesh_xmf(h_mesh, phi, 'phi')

end program test_tri_poisson
