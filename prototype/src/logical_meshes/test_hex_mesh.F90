program test_box_splines

use hex_mesh
implicit none

type(hex_mesh_2d), pointer  :: mesh
!sll_int32    :: num_cells


!num_cells = 4
mesh => new_hex_mesh_2d(4)
!call sll_display(mesh)
call write_hex_mesh_2d(mesh,"hex_mesh_coo.txt")

end program test_box_splines
