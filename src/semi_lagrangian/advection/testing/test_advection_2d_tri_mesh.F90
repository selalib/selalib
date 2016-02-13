program test_advection_2d_tri_mesh
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_advection_2d_tri_mesh, only: &
    sll_s_advection_2d, &
    sll_f_new_advection_2d_tri_mesh, &
    sll_t_advection_tri_mesh

  use sll_m_gnuplot, only: &
    sll_o_gnuplot_2d

  use sll_m_hexagonal_meshes, only: &
    sll_f_new_hex_mesh_2d, &
    sll_t_hex_mesh_2d

  use sll_m_triangular_meshes, only: &
    sll_s_analyze_triangular_mesh, &
    sll_s_map_to_circle, &
    sll_o_new_triangular_mesh_2d, &
    sll_o_delete, &
    sll_t_triangular_mesh_2d, &
    sll_s_write_triangular_mesh_mtv

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

type(sll_t_triangular_mesh_2d), pointer :: t_mesh
type(sll_t_hex_mesh_2d), pointer        :: h_mesh
type(sll_t_advection_tri_mesh), pointer :: t_adv

sll_real64, dimension(:), allocatable :: df
sll_real64, dimension(:), allocatable :: ex 
sll_real64, dimension(:), allocatable :: ey

sll_real64, dimension(:), pointer :: x1
sll_real64, dimension(:), pointer :: x2

sll_int32  :: num_cells
sll_int32  :: ierr
sll_int32  :: istep

sll_real64 :: dt = 0.1_f64

sll_int32, allocatable  :: mitchell_corners(:,:)
sll_int32  :: is1, is2, ic, iv, it, iac, nbc
!sll_int32  :: ic1, ic2, ic3, ic4

!Create a triangular mesh from an hex mesh
!Reference on the boundary is set to "one"

num_cells = 2
h_mesh => sll_f_new_hex_mesh_2d( num_cells, 0._f64, 0._f64) 
t_mesh => sll_o_new_triangular_mesh_2d(h_mesh) 

SLL_CLEAR_ALLOCATE(df(1:t_mesh%num_nodes), ierr)
SLL_CLEAR_ALLOCATE(ex(1:t_mesh%num_nodes), ierr)
SLL_CLEAR_ALLOCATE(ey(1:t_mesh%num_nodes), ierr)

call sll_s_map_to_circle(t_mesh, num_cells)

call sll_s_write_triangular_mesh_mtv(t_mesh, "positions_mesh.mtv")

x1 => t_mesh%coord(1,:)
x2 => t_mesh%coord(2,:)

df = exp(-((x1-0.5)**2+x2*x2)/0.04_f64)
ex = -1.0_f64 !- x2
ey =  0.0_f64 !+ x1

t_adv => sll_f_new_advection_2d_tri_mesh(t_mesh)

do istep = 1, 1
  call sll_s_advection_2d(t_adv, df, ex, ey, dt)
  call sll_o_gnuplot_2d( df, "f_tri", t_mesh%coord, t_mesh%nodes, istep)
end do

print*, 'error =', sum(abs(df-exp(-(x1*x1+x2*x2)/0.04_f64)))/real(t_mesh%num_nodes,f64)

call sll_s_analyze_triangular_mesh(t_mesh) 

!Compute  Mitchell corners
allocate(mitchell_corners(3,t_mesh%num_triangles))

do it = 1, t_mesh%num_triangles
  do iv = 1, 3
    is1 = t_mesh%nodes(mod(iv-1,3)+1,it)
    is2 = t_mesh%nodes(mod(iv  ,3)+1,it) 
    !Nombre de triangles commun avec le noeud is1
    nbc=t_mesh%npoel1(is1+1)-t_mesh%npoel1(is1)
    do ic = 1, nbc
      !Numero du triangle 
      iac=t_mesh%npoel2(t_mesh%npoel1(is1)+ic)
      if ( iac == it) then
        mitchell_corners(iv,it) = ic
        exit
      end if
    end do
  end do
end do

do it = 1, t_mesh%num_triangles
  write(*,"(i3,2x,3i4,2x,3i4)") it, t_mesh%nodes(:,it), mitchell_corners(:,it)
end do
call sll_o_delete(t_mesh)

end program test_advection_2d_tri_mesh
