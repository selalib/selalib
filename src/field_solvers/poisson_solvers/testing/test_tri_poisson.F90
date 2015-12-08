program test_tri_poisson
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_constants, only: &
    sll_pi

  use sll_m_gnuplot, only: &
    sll_gnuplot_2d

  use sll_m_hexagonal_meshes, only: &
    new_hex_mesh_2d, &
    sll_hex_mesh_2d

  use sll_m_timer, only: &
    sll_set_time_mark, &
    sll_time_elapsed_between, &
    sll_time_mark

  use sll_m_tri_poisson, only: &
    new_triangular_poisson_2d, &
    sll_compute_e_from_phi, &
    sll_compute_e_from_rho, &
    sll_compute_phi_from_rho, &
    sll_triangular_poisson_2d, &
    sll_delete

  use sll_m_triangular_meshes, only: &
    analyze_triangular_mesh, &
    map_to_circle, &
    new_triangular_mesh_2d, &
    sll_delete, &
    sll_triangular_mesh_2d, &
    write_triangular_mesh_mtv

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

type(sll_time_mark)     ::  t0
type(sll_time_mark)     ::  t1
sll_real64, pointer     ::  x1(:)
sll_real64, pointer     ::  x2(:)
sll_real64, allocatable :: rho(:)
sll_real64, allocatable :: phi(:)
sll_real64, allocatable :: sol(:)
sll_real64, allocatable :: e_x(:)
sll_real64, allocatable :: e_y(:)
type(sll_triangular_poisson_2d), pointer :: solver
type(sll_triangular_poisson_2d), pointer :: poisson
type(sll_triangular_mesh_2d), pointer    :: square
type(sll_triangular_mesh_2d), pointer    :: t_mesh
type(sll_hex_mesh_2d), pointer           :: h_mesh

sll_int32  :: num_cells
sll_int32  :: ntypfr(5)
sll_real64 :: potfr(5)

sll_int32  :: nc_x1 = 40
sll_real64 :: x1_min = -1.0_f64, x1_max = 1.0_f64
sll_int32  :: nc_x2 = 40
sll_real64 :: x2_min = -1.0_f64, x2_max = 1.0_f64

sll_int32  :: i, ierr
sll_real64 :: r

!First test, the unstructured mesh is created by meshing a square with triangles.

square => new_triangular_mesh_2d(nc_x1, x1_min, x1_max, nc_x2, x2_min, x2_max) 
call analyze_triangular_mesh(square) 
call write_triangular_mesh_mtv(square, "tri_square.mtv")

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
!RHS is ste to zero
rho = 0.0_f64

!Create the Poisson solver on unstructured mesh
poisson => new_triangular_poisson_2d(square, ntypfr, potfr)
!We compute phi
call sll_compute_phi_from_rho(poisson, rho, phi)

!Check result
print*,'error phi=', maxval(abs(phi-square%coord(1,:)))

!We compute ex and ey from phi
call sll_compute_e_from_phi(poisson, phi, e_x, e_y)

!Check result
print*,'error e_x=', maxval(abs(e_x+1.0_f64))
print*,'error e_y=', maxval(abs(e_y))

deallocate(e_x,e_y,phi,rho)

call sll_delete(square)
call sll_delete(poisson)

!Second test, the unstructured mesh is created from an hexagonal mesh.

num_cells = 20
h_mesh => new_hex_mesh_2d( num_cells, 0._f64, 0._f64) 
t_mesh => new_triangular_mesh_2d(h_mesh) 

!The hexagone is mapped to a disk
call map_to_circle(t_mesh, num_cells)
call analyze_triangular_mesh(t_mesh) 
call write_triangular_mesh_mtv(t_mesh, "hex_circle.mtv")

SLL_CLEAR_ALLOCATE(e_x(1:t_mesh%num_nodes),ierr)
SLL_CLEAR_ALLOCATE(e_y(1:t_mesh%num_nodes),ierr)
SLL_CLEAR_ALLOCATE(rho(1:t_mesh%num_nodes),ierr)
SLL_CLEAR_ALLOCATE(phi(1:t_mesh%num_nodes),ierr)
SLL_CLEAR_ALLOCATE(sol(1:t_mesh%num_nodes),ierr)

!The boundary reference is set to "1"
ntypfr(1) = 1
potfr(1)  = 0.0_f64

!Create the Poisson solver on unstructured mesh
solver => new_triangular_poisson_2d(t_mesh, ntypfr, potfr)

!We set the RHS and analytic solution (see functions below)
!Positions of nodes
x1 => t_mesh%coord(1,:)
x2 => t_mesh%coord(2,:)
do i = 1, t_mesh%num_nodes
  r = sqrt(x1(i)*x1(i)+x2(i)*x2(i))
  rho(i) = 4._f64 * sll_pi * f(r)
  sol(i) = u(r)
end do

!Plot fields drawable by gnuplot
call sll_gnuplot_2d( rho, "rho", t_mesh%coord, t_mesh%nodes, 1)
call sll_gnuplot_2d( sol, "sol", t_mesh%coord, t_mesh%nodes, 1)

call sll_set_time_mark(t0)
call sll_compute_e_from_rho(solver, rho, phi, e_x, e_y)
call sll_set_time_mark(t1)
print *, 'Time elapsed to solve Poisson ', sll_time_elapsed_between(t0,t1)

call sll_gnuplot_2d( phi, "phi", t_mesh%coord, t_mesh%nodes, 1)

print*,'error phi=', maxval(abs(phi-sol))

call sll_delete(t_mesh)
call sll_delete(solver)

contains

!We have the equation :
! -4 pi f(r) = Laplacian(u(r))
!If
!  f(r) = rho                    for 0   <= r <= 0.2
!  f(r) = 0.0                    for 0.2 <  r <= 1
!Then
!  u(r) = -pi * rho * r^2 + a_1  for 0   <= r <= 0.2
!  u(r) = a_2 * ln(r)            for 0.2 <  r <= 1
!
!  a_1 =  0.04 * pi * rho * (-2*ln(0.2)+1)
!  a_2 = -0.08 * pi * rho

function u(r)

  sll_real64 :: u
  sll_real64 :: r
  sll_real64 :: pi
  sll_real64, parameter :: one = 1._f64
  sll_real64 :: a_1, a_2

  pi  =  4.0_f64 * atan(1._f64)
  a_1 =  0.04_f64 * pi * one * (-2._f64*log(0.2_f64)+1._f64)
  a_2 = -0.08_f64 * pi * one

  if (0._f64 <= r .and. r <= 0.2_f64) then
    u = -pi * r*r + a_1 
  else if ( 0.2_f64 < r .and. r <= 1._f64) then
    u = a_2 * log(r) 
  else 
    u = 0._f64
  end if

end function u

function f(r)

  sll_real64 :: f
  sll_real64 :: r
  sll_real64, parameter :: one = 1._f64

  if ( 0._f64 <= r .and. r <= 0.2_f64 ) then
    f = one
  else 
    f = 0._f64
  end if

end function f

end program test_tri_poisson
