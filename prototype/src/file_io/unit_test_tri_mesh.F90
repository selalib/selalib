program test_tri_mesh

use sll_generate_tri_mesh

implicit none

integer :: nbox, nboy
real(8) :: alx, aly
real(8), dimension(:,:), allocatable :: coor
integer, dimension(:,:), allocatable :: ntri
integer :: nbt
integer :: nbs
real(8) :: x_min, x_max, y_min, y_max


write(6,*) '-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-'
write(6,*)
write(6,*) 'Mesh Generation on a square '
write(6,*)
write(6,*) '-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-'

x_min = 0.; x_max = 1.0
y_min = 0.; y_max = 1.0

nbox = 10; nboy = 10
alx = 1.;  aly = 1.

call plaqx( x_min, x_max, nbox, &
            y_min, y_max, nboy, &
            coor, ntri, nbs, nbt)

call plaqy( x_min, x_max, nbox, &
            y_min, y_max, nboy, &
            coor, ntri, nbs, nbt)


end program test_tri_mesh
