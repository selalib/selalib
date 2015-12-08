program test_hexagonal_meshes

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_constants, only: &
    sll_pi

  use sll_m_hexagonal_meshes, only: &
    delete, &
    new_hex_mesh_2d, &
    sll_hex_mesh_2d

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type(sll_hex_mesh_2d), pointer  :: mesh
  sll_int32                   :: num_cells
  sll_real64, pointer         :: field(:)
  sll_int32                   :: error
  sll_real64                  :: x1
  sll_real64                  :: x2
  sll_int32                   :: i
  sll_int32                   :: nei1
  sll_int32                   :: nei2
  sll_int32                   :: nei3
  sll_int32                   :: type
  sll_int32                   :: spline_degree

  num_cells = 4
  spline_degree = 1

  print *, ""
  print *, "Creating a mesh with", num_cells, &
       "cells, mesh coordinates written in ./hex_mesh_coo.txt"
  mesh => new_hex_mesh_2d(num_cells)
  call mesh%display()
  call mesh%write_hex_mesh_2d( "hex_mesh_coo.txt")
  call mesh%write_hex_mesh_mtv("hex_mesh_coo.mtv")
  print *, ""

  SLL_ALLOCATE(field(mesh%num_pts_tot), error)

  do i = 1, mesh%num_pts_tot
     x1 = mesh%global_to_x1(i)
     x2 = mesh%global_to_x2(i)
     field(i) = cos(2*sll_pi*x1)*sin(2*sll_pi*x2)
  end do

  call mesh%write_field_hex_mesh_xmf(field, 'field')

  call delete(mesh)

  ! TESTING NEIGHBOURS :
  num_cells = 2
  mesh => new_hex_mesh_2d(num_cells)

  do i = 1, mesh%num_triangles
     call mesh%get_neighbours(i, nei1, nei2, nei3)
     call mesh%cell_type(i, type)
     print *, "i =", i, "type:", type, "neighbourcells =", nei1, nei2, nei3
  end do

  call delete(mesh)

  print *, 'PASSED'

end program test_hexagonal_meshes
