!> @internal [example]
program test_triangular_meshes
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"

  use sll_m_hexagonal_meshes, only: &
    sll_f_new_hex_mesh_2d, &
    sll_t_hex_mesh_2d

  use sll_m_triangular_meshes, only: &
    sll_s_map_to_circle, &
    sll_o_new_triangular_mesh_2d, &
    sll_s_read_from_file, &
    sll_o_delete, &
    sll_o_display, &
    sll_t_triangular_mesh_2d, &
    sll_s_write_triangular_mesh_mtv

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type(sll_t_triangular_mesh_2d), pointer :: t_mesh
  type(sll_t_hex_mesh_2d), pointer        :: h_mesh

  sll_int32    :: nc_x1  = 4
  sll_real64   :: x1_min = 0.0_f64
  sll_real64   :: x1_max = 1.0_f64
  sll_int32    :: nc_x2  = 4
  sll_real64   :: x2_min = 0.0_f64
  sll_real64   :: x2_max = 1.0_f64
  sll_int32    :: num_cells
  logical                          :: file_exists
  character(len=256)               :: reference_filename
  !Create a triangular mesh from square
  t_mesh => sll_o_new_triangular_mesh_2d(nc_x1, x1_min, x1_max, nc_x2, x2_min, x2_max)

  call sll_o_display(t_mesh)
  call sll_s_write_triangular_mesh_mtv(t_mesh, "tri_mesh.mtv")
  call sll_o_delete(t_mesh)

  ! Read name of reference file from input argument
  !------------------------------------------------
  call get_command_argument( 1, reference_filename )

  ! Check that file exists
  !-----------------------
  inquire( file=trim( reference_filename ), exist=file_exists )
  if (.not. file_exists) then
     write(*,*)  "ERROR: reference file '"&
          //trim( reference_filename )//"' does not exist"
     stop
  end if

  !Create a triangular mesh from a file
  call sll_s_read_from_file(t_mesh, reference_filename)
  call sll_o_display(t_mesh)
  call sll_s_write_triangular_mesh_mtv(t_mesh, "diode_mesh.mtv")
  call sll_o_delete(t_mesh)

  !Create a triangular mesh from an hex mesh
  !Reference on the boundary is set to "one"
  num_cells = 3
  h_mesh => sll_f_new_hex_mesh_2d( num_cells, 0._f64, 0._f64)
  t_mesh => sll_o_new_triangular_mesh_2d(h_mesh)

  call sll_s_map_to_circle(t_mesh, num_cells, 1)
  call sll_s_write_triangular_mesh_mtv(t_mesh, "circle_hex_mesh.mtv")


  call sll_o_delete(t_mesh)

  print *, 'PASSED'

end program test_triangular_meshes
!> @internal [example]
