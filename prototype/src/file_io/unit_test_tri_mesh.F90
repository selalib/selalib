program test_tri_mesh

#include "sll_constants.h"
use sll_tri_mesh_xmf
use sll_generate_tri_mesh

implicit none

integer :: nbox, nboy
real(8), dimension(:,:), pointer :: coor1, coor2
integer, dimension(:,:), pointer :: ntri1, ntri2
integer :: nbt1, nbt2
integer :: nbs1, nbs2
real(8) :: x_min, x_max, y_min, y_max
real(8), dimension(:), pointer :: field
integer :: i


write(6,*) '-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-'
write(6,*)
write(6,*) 'Mesh Generation on a square '
write(6,*)
write(6,*) '-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-'

x_min = 0.; x_max = 1.0
y_min = 0.; y_max = 1.0

nbox = 33; nboy = 33

!call plaqx( x_min, x_max, nbox, &
!            y_min, y_max, nboy, &
!            coor1, ntri1, nbs1, nbt1)

call plaqy( x_min, x_max, nbox, &
            y_min, y_max, nboy, &
            coor2, ntri2, nbs2, nbt2)

allocate(field(nbs2))
field  = cos(2*sll_pi*coor2(1,:))*sin(2*sll_pi*coor2(2,:))
call write_tri_mesh_xmf("plaqy", coor2, ntri2, nbs2, nbt2, field)

end program test_tri_mesh
