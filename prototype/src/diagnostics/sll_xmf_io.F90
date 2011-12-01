module sll_xmf_io
#include "sll_working_precision.h"
#include "sll_assert.h"

contains

subroutine sll_xmf_file_create(filename,file_id,error)
character(len=*) , intent(in)  :: filename   
sll_int32        , intent(out) :: file_id   
sll_int32        , intent(out) :: error
logical                        :: lopen

error=0

do 100 file_id=20,99

   inquire(unit=file_id,opened=lopen)
   if(lopen) then
      cycle
   else
      open(file_id,status='SCRATCH',err=100)
      close(file_id,status='DELETE',err=100)
      goto 200
   end if
 
100 continue
error=1
200 continue
error=0

! Create a new file using default properties
inquire(file=filename,opened=lopen)
SLL_ASSERT(.not. lopen)

open(file_id,FILE=filename,FORM='FORMATTED',IOSTAT=error)
rewind(file_id)

write(file_id,"(a)")"<?xml version='1.0' ?>"
write(file_id,"(a)")"<!DOCTYPE Xdmf SYSTEM 'Xdmf.dtd' []>"
write(file_id,"(a)")"<Xdmf Version='2.0'>"
write(file_id,"(a)")"<Domain>"

end subroutine sll_xmf_file_create

subroutine sll_xmf_file_close(file_id,error)
sll_int32, intent(in)  :: file_id
sll_int32, intent(out) :: error
write(file_id,"(a)")"</Domain>"
write(file_id,"(a)")"</Xdmf>"
close(file_id)
error = 0
end subroutine sll_xmf_file_close

subroutine sll_xmf_grid_geometry_2d(file_id, nnodes_x1, nnodes_x2)
sll_int32, intent(in) :: file_id
sll_int32, intent(in) :: nnodes_x1
sll_int32, intent(in) :: nnodes_x2
write(file_id,"(a)")"<Grid Name='mesh' GridType='Uniform'>"
write(file_id,"(a,2i5,a)")"<Topology TopologyType='2DSMesh' NumberOfElements='", &
                     nnodes_x2, nnodes_x1,"'/>"
write(file_id,"(a)")"<Geometry GeometryType='X_Y'>"
end subroutine sll_xmf_grid_geometry_2d

subroutine sll_xmf_dataitem_2d(file_id, filename, nnodes_x1, nnodes_x2, filetype)
sll_int32, intent(in) :: file_id
character(len=*), intent(in) :: filename
character(len=*), intent(in) :: filetype
sll_int32, intent(in) :: nnodes_x1
sll_int32, intent(in) :: nnodes_x2
write(file_id,"(a,2i5,a)")"<DataItem Dimensions='",nnodes_x2,nnodes_x1, &
"' NumberType='Float' Precision='8' Format='"//trim(filetype)//"'>"
write(file_id,"(a)")trim(filename)
write(file_id,"(a)")"</DataItem>"
end subroutine sll_xmf_dataitem_2d


subroutine sll_xmf_open_grid_2d(file_id, filename, nnodes_x1, nnodes_x2)
sll_int32, intent(in) :: file_id
character(len=*), intent(in) :: filename
sll_int32, intent(in) :: nnodes_x1
sll_int32, intent(in) :: nnodes_x2

write(file_id,"(a)")"<Grid Name='mesh' GridType='Uniform'>"
write(file_id,"(a,2i5,a)")"<Topology TopologyType='2DSMesh' NumberOfElements='", &
                          nnodes_x2,nnodes_x1,"'/>"
write(file_id,"(a)")"<Geometry GeometryType='X_Y'>"

#ifdef NOHDF5

call sll_xmf_dataitem_2d(file_id,trim(filename)//"-x1.bin",nnodes_x1,nnodes_x2,'Binary')
call sll_xmf_dataitem_2d(file_id,trim(filename)//"-x2.bin",nnodes_x1,nnodes_x2,'Binary')

#else

call sll_xmf_dataitem_2d(file_id,trim(filename)//".h5:/x1",nnodes_x1,nnodes_x2,'HDF')
call sll_xmf_dataitem_2d(file_id,trim(filename)//".h5:/x2",nnodes_x1,nnodes_x2,'HDF')

#endif

write(file_id,"(a)")"</Geometry>"

end subroutine sll_xmf_open_grid_2d



end module sll_xmf_io
