program test_poisson_2d_tri
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

use sll_m_constants, only: &
  sll_p_pi

use sll_m_gnuplot, only: &
  sll_o_gnuplot_2d

use sll_m_hexagonal_meshes, only: &
  sll_f_new_hex_mesh_2d, &
  sll_t_hex_mesh_2d

use sll_m_timer, only: &
  sll_s_set_time_mark, &
  sll_f_time_elapsed_between, &
  sll_t_time_mark

use sll_m_poisson_2d_tri, only: &
  sll_f_new_triangular_poisson_2d, &
  sll_s_compute_e_from_phi, &
  sll_s_compute_e_from_rho, &
  sll_s_compute_phi_from_rho, &
  sll_t_triangular_poisson_2d, &
  sll_o_delete

use sll_m_triangular_meshes, only: &
  sll_s_analyze_triangular_mesh, &
  sll_s_map_to_circle, &
  sll_o_new_triangular_mesh_2d, &
  sll_o_delete, &
  sll_t_triangular_mesh_2d, &
  sll_s_write_triangular_mesh_mtv

implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

sll_real64,                        allocatable :: rho(:)
sll_real64,                        allocatable :: phi(:)
sll_real64,                        allocatable :: e_x(:)
sll_real64,                        allocatable :: e_y(:)

type(sll_t_triangular_poisson_2d), pointer     :: poisson
type(sll_t_triangular_mesh_2d),    pointer     :: square

sll_int32  :: ntypfr(5)
sll_real64 :: potfr(5)
sll_int32  :: nc_x1 = 40
sll_real64 :: x1_min = -1.0_f64
sll_real64 :: x1_max = 1.0_f64
sll_int32  :: nc_x2 = 40
sll_real64 :: x2_min = -1.0_f64
sll_real64 :: x2_max = 1.0_f64
sll_int32  :: ierr

!type(sll_t_triangular_poisson_2d), pointer :: solver
!type(sll_t_triangular_mesh_2d),    pointer :: t_mesh
!type(sll_t_hex_mesh_2d),           pointer :: h_mesh
!sll_int32                                  :: num_cells
!type(sll_t_time_mark)                      :: t0
!type(sll_t_time_mark)                      :: t1
!sll_real64,                        pointer :: x1(:)
!sll_real64,                        pointer :: x2(:)
!sll_int32                                  :: i
!sll_real64                                 :: r
!sll_real64,                    allocatable :: sol(:)

!First test, the unstructured mesh is created 
!by meshing a square with triangles.

square => sll_o_new_triangular_mesh_2d(nc_x1, x1_min, x1_max, &
                                       nc_x2, x2_min, x2_max) 

call sll_s_analyze_triangular_mesh(square) 
call sll_s_write_triangular_mesh_mtv(square, "tri_square.mtv")

!ref 1 is the south boundary (Neumann)
!ref 2 is the east  boundary (Dirichlet phi = +1)
!ref 3 is the north boundary (Neumann)
!ref 4 is the west  boundary (Dirichlet phi = -1)

ntypfr(1:4) = [3,1,3,1]
potfr(1:4)  = [0.0_f64,1.0_f64,0.0_f64,-1.0_f64]

SLL_CLEAR_ALLOCATE(e_x(1:square%num_nodes),ierr)
SLL_CLEAR_ALLOCATE(e_y(1:square%num_nodes),ierr)
SLL_CLEAR_ALLOCATE(rho(1:square%num_nodes),ierr)
SLL_CLEAR_ALLOCATE(phi(1:square%num_nodes),ierr)

rho = 0.0_f64

!Create the Poisson solver on unstructured mesh
poisson => sll_f_new_triangular_poisson_2d(square, ntypfr, potfr)
!We compute phi
call sll_s_compute_phi_from_rho(poisson, rho, phi)

!Check result
print*,'error phi=', maxval(abs(phi-square%coord(1,:)))

!We compute ex and ey from phi
call sll_s_compute_e_from_phi(poisson, phi, e_x, e_y)

!Check result
print*,'error e_x=', maxval(abs(e_x+1.0_f64))
print*,'error e_y=', maxval(abs(e_y))

deallocate(e_x,e_y,phi,rho)

call sll_o_delete(square)
call sll_o_delete(poisson)

print*, 'PASSED'

!!Second test, the unstructured mesh is created from an hexagonal mesh.
!
!num_cells = 50
!h_mesh => sll_f_new_hex_mesh_2d( num_cells, 0._f64, 0._f64) 
!t_mesh => sll_o_new_triangular_mesh_2d(h_mesh) 
!
!!The hexagone is mapped to a disk
!call sll_s_map_to_circle(t_mesh, num_cells)
!call sll_s_analyze_triangular_mesh(t_mesh) 
!call sll_s_write_triangular_mesh_mtv(t_mesh, "hex_circle.mtv")
!
!SLL_CLEAR_ALLOCATE(e_x(1:t_mesh%num_nodes),ierr)
!SLL_CLEAR_ALLOCATE(e_y(1:t_mesh%num_nodes),ierr)
!SLL_CLEAR_ALLOCATE(rho(1:t_mesh%num_nodes),ierr)
!SLL_CLEAR_ALLOCATE(phi(1:t_mesh%num_nodes),ierr)
!SLL_CLEAR_ALLOCATE(sol(1:t_mesh%num_nodes),ierr)
!
!!The boundary reference is set to "1"
!ntypfr(1) = 1
!potfr(1)  = 0.0_f64
!
!!Create the Poisson solver on unstructured mesh
!solver => sll_f_new_triangular_poisson_2d(t_mesh, ntypfr, potfr)
!
!!We set the RHS and analytic solution (see functions below)
!!Positions of nodes
!x1 => t_mesh%coord(1,:)
!x2 => t_mesh%coord(2,:)
!do i = 1, t_mesh%num_nodes
!  r = sqrt(x1(i)*x1(i)+x2(i)*x2(i))
!  rho(i) = 4._f64 * sll_p_pi * f(r)
!  sol(i) = u(r)
!end do
!
!!Plot fields drawable by gnuplot
!call sll_o_gnuplot_2d( rho, "rho", t_mesh%coord, t_mesh%nodes, 1)
!call sll_o_gnuplot_2d( sol, "sol", t_mesh%coord, t_mesh%nodes, 1)
!
!call sll_s_set_time_mark(t0)
!call sll_s_compute_e_from_rho(solver, rho, phi, e_x, e_y)
!call sll_s_set_time_mark(t1)
!print *, 'Time elapsed to solve Poisson ', sll_f_time_elapsed_between(t0,t1)
!
!call sll_o_gnuplot_2d( phi, "phi", t_mesh%coord, t_mesh%nodes, 1)
!call sll_o_gnuplot_2d( sol, "sol", t_mesh%coord, t_mesh%nodes, 1)
!
!print*,'error phi=', maxval(abs(phi-sol))
!
!call sll_o_delete(t_mesh)
!call sll_o_delete(solver)

contains

!Charge density is a solid cylinder of radius 0.2
function f(r)

  sll_real64 :: f
  sll_real64 :: r

  if ( 0._f64 <= r .and. r <= 0.2_f64 ) then
    f = 1.0_f64
  else 
    f = 0.0_f64
  end if

end function f

!We have the equation :
! -4 pi f(r) = Laplacian(u(r))
!If
!  f(r) = rho                    for 0   <= r <= 0.2
!  f(r) = 0.0                    for 0.2 <  r <= 1
!Then
!  u(r) = -pi * f(r) * r^2 + a_1  for 0   <= r <= 0.2
!  u(r) = a_2 * ln(r)            for 0.2 <  r <= 1
!
!  a_1 =  0.04 * pi * f(r) * (-2*ln(0.2)+1)
!  a_2 = -0.08 * pi * f(r)

function u(r)

  sll_real64 :: u
  sll_real64 :: r
  sll_real64 :: pi
  sll_real64 :: a_0, a_1, a_2, a_3

  pi  =  4.0_f64 * atan(1._f64)
  a_0 =  0.0_f64
  a_1 =  0.04_f64 * pi * f(r) * (-2.0_f64*log(0.2_f64)+1.0_f64)
  a_2 = -0.08_f64 * pi !* f(r)
  a_3 =  0.0_f64

  if (0._f64 < r .and. r <= 0.2_f64) then
    u = -pi * f(r) * r*r + a_0*log(r) + a_1 
  else if ( 0.2_f64 < r .and. r <= 1._f64) then
    u = a_2 * log(r) + a_3
  else 
    u = 0._f64
  end if

end function u

end program test_poisson_2d_tri
