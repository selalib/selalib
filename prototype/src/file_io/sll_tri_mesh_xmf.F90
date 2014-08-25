module sll_tri_mesh_xmf
implicit none

!*****************************************
!   XDMF output for VisIt and paraview
!*****************************************

contains


subroutine write_tri_mesh_xmf(filename, coor, ntri, nbs, nbt, field)

character(len=*), intent(in) :: filename
integer, intent(in)          :: nbs
integer, intent(in)          :: nbt
real(8), intent(in)          :: coor(2,nbs)
integer, intent(in)          :: ntri(3,nbt)
integer, parameter           :: xmf = 99
real(8), intent(in)          :: field(:)
integer                      :: i, j

open(xmf, file=filename//".xmf")

write(xmf,"(a)") "<?xml version='1.0'?>"
write(xmf,"(a)") "<!DOCTYPE Xdmf SYSTEM 'Xdmf.dtd' []> "
write(xmf,"(a)") "<Xdmf Version='2.0'>"
write(xmf,"(a)") "<Domain>"
write(xmf,"(a)") "<Grid Name='Mesh' GridType='Uniform' >"
write(xmf,"(a,i6,a)") "<Topology Type='Triangle' NumberOfElements='",nbt,"'>"
write(xmf,"(a,i6,a)") "<DataItem Name='Connections' Format='XML' DataType='Int' Dimensions='", &
                      nbt, " 3'>" 
do i=1, nbt
   write(xmf,"(3i6)") (ntri(j,i)-1,j=1,3)
end do
write(xmf,"(a)") "</DataItem>"
write(xmf,"(a)") "</Topology>"
write(xmf,"(a)") "<Geometry GeometryType='XY'>"
write(xmf,"(a,i6,a)") "<DataItem Format='XML' Dimensions='", nbs, " 2'>"
do  i=1, nbs
   write(xmf,"(2f12.6)") coor(:,i)
end do
write(xmf,"(a)") "</DataItem>"
write(xmf,"(a)") "</Geometry>"
write(xmf,"(a)") "<Attribute Name='Node_centered_field' Center='Node'>"
write(xmf,"(a,i6,a)") "<DataItem Format='XML' Datatype='Float' Dimensions='", nbs, "'>"
do  i=1, nbs
   write(xmf,"(f12.6)") field(i)
end do
write(xmf,"(a)") "</DataItem>"
write(xmf,"(a)") "</Attribute>"
write(xmf,"(a)") "</Grid>"
write(xmf,"(a)") "</Domain>"
write(xmf,"(a)") "</Xdmf>" 

end subroutine write_tri_mesh_xmf

end module sll_tri_mesh_xmf
