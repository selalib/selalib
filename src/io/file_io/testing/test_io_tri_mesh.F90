program test_tri_mesh

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  use sll_generate_tri_mesh, only: &
    plaqx, &
    plaqy

  use sll_m_constants, only: &
    sll_p_pi

  use sll_m_tri_mesh_xmf, only: &
    sll_s_write_tri_mesh_xmf

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

integer :: nbox, nboy
real(8), dimension(:,:), pointer :: coor1, coor2
integer, dimension(:,:), pointer :: ntri1, ntri2
integer :: nbt1, nbt2
integer :: nbs1, nbs2
real(8) :: x_min, x_max, y_min, y_max
real(8), dimension(:), pointer :: field1
real(8), dimension(:), pointer :: field2

write(6,*) '-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-'
write(6,*)
write(6,*) 'Mesh Generation on a square '
write(6,*)
write(6,*) '-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-'

x_min = 0.0_8; x_max = 1.0_8
y_min = 0.0_8; y_max = 1.0_8

nbox = 33; nboy = 33

call plaqx( x_min, x_max, nbox, &
            y_min, y_max, nboy, &
            coor1, ntri1, nbs1, nbt1)

allocate(field1(nbs2))
field1  = cos(2*sll_p_pi*coor1(1,:))*sin(2*sll_p_pi*coor1(2,:))
call sll_s_write_tri_mesh_xmf("tri_mesh_1", coor1, ntri1, nbs1, nbt1, field1, 'field_1')

call plaqy( x_min, x_max, nbox, &
            y_min, y_max, nboy, &
            coor2, ntri2, nbs2, nbt2)

allocate(field2(nbs2))
field2  = cos(2*sll_p_pi*coor2(1,:))*sin(2*sll_p_pi*coor2(2,:))
call sll_s_write_tri_mesh_xmf("tri_mesh_2", coor2, ntri2, nbs2, nbt2, field2, 'field_2')

end program test_tri_mesh
