program test_hexagonal_meshes

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_constants, only: &
    sll_p_pi

  use sll_m_hexagonal_meshes, only: &
    sll_o_delete, &
    sll_f_new_hex_mesh_2d, &
    sll_t_hex_mesh_2d

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type(sll_t_hex_mesh_2d), pointer  :: mesh
  sll_int32                   :: num_cells
  sll_real64, pointer         :: field(:)
  sll_int32                   :: error
  sll_real64                  :: x1
  sll_real64                  :: x2
  sll_real64                  :: x1_new
  sll_real64                  :: x2_new
  sll_int32                   :: i
  sll_int32                   :: nei1
  sll_int32                   :: nei2
  sll_int32                   :: nei3
  sll_int32                   :: type
  sll_int32                   :: spline_degree
  sll_real64,  dimension(2,2) :: transf_matA
  sll_real64,  dimension(2)   :: transf_vecB

  num_cells = 4
  spline_degree = 1

  print *, ""
  print *, "Creating a mesh with", num_cells, &
       "cells, mesh coordinates written in ./hex_mesh_coo.txt"
  mesh => sll_f_new_hex_mesh_2d(num_cells)
  call mesh%display()
  call mesh%sll_s_write_hex_mesh_2d( "hex_mesh_coo.txt")
  call mesh%sll_s_write_hex_mesh_mtv("hex_mesh_coo.mtv")
  print *, ""

  SLL_ALLOCATE(field(mesh%num_pts_tot), error)

  do i = 1, mesh%num_pts_tot
     x1 = mesh%global_to_x1(i)
     x2 = mesh%global_to_x2(i)
     field(i) = cos(2*sll_p_pi*x1)*sin(2*sll_p_pi*x2)
  end do

  call mesh%sll_s_write_field_hex_mesh_xmf(field, 'field')

  call sll_o_delete(mesh)

  ! TESTING NEIGHBOURS :
  num_cells = 2
  mesh => sll_f_new_hex_mesh_2d(num_cells)

  do i = 1, mesh%num_triangles
     call mesh%get_neighbours(i, nei1, nei2, nei3)
     call mesh%cell_type(i, type)
     print *, "i =", i, "type:", type, "neighbourcells =", nei1, nei2, nei3
     print *, "element =", i
     call mesh%hex_to_aligned_elmt(i, "CIRCLE", transf_matA, transf_vecB)
     print *, transf_matA(1, :), transf_matA(2, :)
     print *, transf_vecB(:)
     print *, "---------------------------------"
  end do


  print *, transf_matA(:, :)
  print *, transf_vecB(:)

  call sll_o_delete(mesh)

  print *, 'PASSED'

end program test_hexagonal_meshes
