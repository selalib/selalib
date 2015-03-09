program test_tri_poisson
#include "sll_working_precision.h"
#include "sll_memory.h"
                        
use sll_hex_meshes
use tri_poisson
use sll_triangular_meshes
use sll_mesh_calculus_2d_module
use sll_gnuplot
use sll_timer

!----------------------------------------------------------------------

implicit none

type(sll_time_mark)  :: t0 
type(sll_time_mark)  :: t1 
real(8), dimension(:), pointer        :: x1
real(8), dimension(:), pointer        :: x2
real(8), dimension(:), allocatable    :: rho
real(8), dimension(:), allocatable    :: phi
real(8), dimension(:), allocatable    :: sol
real(8), dimension(:), allocatable    :: e_x
real(8), dimension(:), allocatable    :: e_y
type(sll_triangular_poisson_2d)       :: solver
type(sll_triangular_mesh_2d), pointer :: t_mesh
type(sll_hex_mesh_2d), pointer        :: h_mesh

sll_int32  :: num_cells
sll_int32  :: ntypfr(5)
sll_real64 :: potfr(5)

sll_int32  :: nc_x1 = 100
sll_real64 :: x1_min = -1.0_f64, x1_max = 1.0_f64
sll_int32  :: nc_x2 = 100
sll_real64 :: x2_min = -1.0_f64, x2_max = 1.0_f64

sll_int32  :: i, ierr
sll_real64 :: r

!mesh => new_triangular_mesh_2d("diode.maa") 
!t_mesh => new_triangular_mesh_2d(nc_x1, x1_min, x1_max, nc_x2, x2_min, x2_max) 

num_cells = 15
h_mesh => new_hex_mesh_2d( num_cells, 0._f64, 0._f64) 
t_mesh => new_triangular_mesh_2d(h_mesh) 
call map_to_circle(t_mesh, num_cells)
call analyze_triangular_mesh(t_mesh) 
call write_triangular_mesh_mtv(t_mesh, "test_tri_poisson.mtv")

SLL_CLEAR_ALLOCATE(e_x(1:t_mesh%num_nodes),ierr)
SLL_CLEAR_ALLOCATE(e_y(1:t_mesh%num_nodes),ierr)
SLL_CLEAR_ALLOCATE(rho(1:t_mesh%num_nodes),ierr)
SLL_CLEAR_ALLOCATE(phi(1:t_mesh%num_nodes),ierr)
SLL_CLEAR_ALLOCATE(sol(1:t_mesh%num_nodes),ierr)

ntypfr(1) = 1
potfr(1)  = 0.0_f64

call sll_create(solver, t_mesh, ntypfr, potfr)

x1 => t_mesh%coord(1,:)
x2 => t_mesh%coord(2,:)

do i = 1, t_mesh%num_nodes
  r = sqrt(x1(i)*x1(i)+x2(i)*x2(i))
  rho(i) = 4 * sll_pi * f(r)
  sol(i) = u(r)
end do

call sll_gnuplot_2d( rho, "rho", t_mesh%coord, t_mesh%nodes, 1)
call sll_gnuplot_2d( sol, "sol", t_mesh%coord, t_mesh%nodes, 1)

call sll_compute_phi_from_rho(solver, rho, phi)

call sll_gnuplot_2d( phi, "phi", t_mesh%coord, t_mesh%nodes, 1)

print*,'error phi=', maxval(abs(phi-sol))

!call sll_compute_e_from_phi(solver, phi, e_x, e_y)
!
!print*,'error phi=', maxval(abs(phi-x1))
!print*,'error e_x=', maxval(abs(e_x+1.0_f64))
!print*,'error e_y=', maxval(abs(e_y))
!
!call sll_set_time_mark(t0)
!do i = 1, 100
!  call sll_compute_e_from_rho(solver, rho, phi, e_x, e_y)
!end do
!call sll_set_time_mark(t1)
!print *, 'Time elapsed to solve Poisson ', sll_time_elapsed_between(t0,t1)
!
!
!print*,'error phi=', maxval(abs(phi-x1))
!print*,'error e_x=', maxval(abs(e_x+1.0_f64))
!print*,'error e_y=', maxval(abs(e_y))
!ntypfr(2) = 1; potfr(2) = +1.0_f64
!ntypfr(3) = 3; potfr(3) =  0.0_f64
!ntypfr(4) = 1; potfr(4) = -1.0_f64
!call sll_gnuplot_2d( phi, "phi", t_mesh%coord, t_mesh%nodes, 1)
!call sll_gnuplot_2d( e_x, "e_x", t_mesh%coord, t_mesh%nodes, 1)
!call sll_gnuplot_2d( e_y, "e_y", t_mesh%coord, t_mesh%nodes, 1)

!rho = 8*sll_pi*sll_pi*sin(2*sll_pi*x1)*sin(2*sll_pi*x2)
!
!call sll_set_time_mark(t0)
!call sll_compute_phi_from_rho(solver, rho, phi)
!call sll_set_time_mark(t1)
!print*,'error phi=', maxval(abs(phi-sin(2*sll_pi*x1)*sin(2*sll_pi*x2)))
!
!print *, 'Time elapsed to solve Poisson ', sll_time_elapsed_between(t0,t1)
!
!rho = 0.0_f64
!call sll_compute_e_from_phi(solver, phi, e_x, e_y)
!
!print*,'error e_x=', maxval(abs(e_x+2*sll_pi*cos(2*sll_pi*x1)*sin(2*sll_pi*x2)))
!print*,'error e_y=', maxval(abs(e_y+2*sll_pi*sin(2*sll_pi*x1)*cos(2*sll_pi*x2)))
!
!call sll_gnuplot_2d( phi, "phi", t_mesh%coord, t_mesh%nodes, 1)
!call sll_gnuplot_2d( e_x, "e_x", t_mesh%coord, t_mesh%nodes, 1)
!call sll_gnuplot_2d( e_y, "e_y", t_mesh%coord, t_mesh%nodes, 1)

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

  real(8) :: u
  real(8) :: r
  real(8) :: pi 
  real(8), parameter :: one = 1d0 
  real(8) :: a_0, a_1, a_2, a_3

  pi  =  4d0 * atan(1d0)
  a_1 =  0.04 * pi * one * (-2*log(0.2d0)+1)
  a_2 = -0.08 * pi * one

  if ( 0d0 <= r .and. r <= 0.2d0 ) then
    u = -pi * r*r + a_1 
  else if ( 0.2d0 < r .and. r <= 1d0) then
    u = a_2 * log(r) 
  else 
    u = 0.0d0
  end if

end function u

function f(r)

  real(8) :: f
  real(8) :: r
  real(8), parameter :: one = 1d0 

  if ( 0d0 <= r .and. r <= 0.2d0 ) then
    f = one
  else 
    f = 0.0d0
  end if

end function f

end program test_tri_poisson
