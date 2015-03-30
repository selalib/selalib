program unit_test_advection_2d
#include "sll_working_precision.h"
#include "sll_memory.h"
use sll_hex_meshes
use sll_triangular_meshes
use sll_advection_2d_tri_mesh
use sll_gnuplot
use sll_mesh_calculus_2d_module

implicit none

type(sll_triangular_mesh_2d), pointer :: t_mesh
type(sll_hex_mesh_2d), pointer        :: h_mesh
type(sll_advection_tri_mesh), pointer :: t_adv

sll_real64, dimension(:), allocatable :: df
sll_real64, dimension(:), allocatable :: ex 
sll_real64, dimension(:), allocatable :: ey

sll_real64, dimension(:), pointer :: x1
sll_real64, dimension(:), pointer :: x2

sll_int32  :: num_cells
sll_int32  :: ierr
sll_int32  :: istep

sll_real64 :: dt = 0.1

sll_int32, allocatable  :: mitchell_corners(:,:)
sll_int32  :: is1, is2, is3, ic, it, iac, nbc, ic1, ic2

!Create a triangular mesh from an hex mesh
!Reference on the boundary is set to "one"

num_cells = 3
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
ex = -1. !- x2
ey =  0. !+ x1

t_adv => new_advection_2d_tri_mesh(t_mesh)

do istep = 1, 5
  call advection_2d(t_adv, df, ex, ey, dt)
  call sll_gnuplot_2d( df, "f_tri", t_mesh%coord, t_mesh%nodes, istep)
end do

print*, 'error =', sum(abs(df-exp(-(x1*x1+x2*x2)/0.04_f64)))/t_mesh%num_nodes

call analyze_triangular_mesh(t_mesh) 

!Compute  Mitchell corners
allocate(mitchell_corners(3,t_mesh%num_triangles))

do it = 1, t_mesh%num_triangles !Loop  over triangles

  is1 = t_mesh%nodes(1,it)
  is2 = t_mesh%nodes(2,it)
  is3 = t_mesh%nodes(3,it)

  print*, '##############'
  print*, is1 !, is2, is3
  !Nombre de cote avec le noeud is1
  nbc=t_mesh%nbcov(is1+1)-t_mesh%nbcov(is1)
  !Numero du premier cote
  !Numero des sommets a l'autre extermite
  do ic = 1, nbc
    !Numero du cote 
    iac=t_mesh%nbcov(is1)+1
    !Numero globaux des extremites
    ic1 = t_mesh%nuvac(1,t_mesh%nugcv(iac+ic-1))
    ic2 = t_mesh%nuvac(2,t_mesh%nugcv(iac+ic-1))
    if(ic1 == is1 .and. ic2 == is2) then
      mitchell_corners(1,it) = ic
    else if(ic1 == is1 .and. ic1 == is3) then
      mitchell_corners(3,it) = ic
    else if(ic1 == is1 .and. ic1 == is3) then
      mitchell_corners(2,it) = ic
    end if
  end do

end do



call sll_delete(t_mesh)

end program unit_test_advection_2d
