!------------------------------------------------------------------------------
! SELALIB
!------------------------------------------------------------------------------
!
! MODULE: sll_xmf_io
!
!> @author
!> Pierre Navaro
!>
!
! DESCRIPTION: 
!
!> @brief
!> Implements the functions to write xmf file to store light data
!>
!>@details
!> With XDMF file you can describe data to plot them with VisIt
!>
!> <h2>How to use this module: </h2>
!>
!> \code use sll_xmf_io \endcode
!>
!> External links:
!> - https://wci.llnl.gov/codes/visit/
!> - http://www.xdmf.org/index.php/Main_Page
!>
!><h2> Example </h2>
!>The call sequence to create an xdmf file is
!>
!><pre>
!>
!>call sll_xmf_file_create("example.xmf",xmffile_id,error)
!>call sll_xmf_grid_geometry_2d(xmffile_id, "mesh", nnodes_x1, nnodes_x2)
!>call sll_xmf_field_2d(xmffile_id,'NodeValues',"field",nnodes_x1,nnodes_x2,'Node')
!>call sll_xmf_file_close(xmffile_id,error)
!>
!></pre>
!>Don't forget to write heavy data in hdf5 file or binary file !
!
! REVISION HISTORY:
! 05 12 2011 - Initial Version
! TODO_dd_mmm_yyyy - TODO_describe_appropriate_changes - TODO_name
!------------------------------------------------------------------------------
module sll_xmf_io
#include "sll_working_precision.h"
#include "sll_assert.h"

contains

!> Create the XML file and begin to write first lines.
!> You get the file unit number.
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

!> Close the XML file and finish to write last lines.
!> You give the file unit number.
!> \param[in] file_id is the unit number or your xmf file
subroutine sll_xmf_file_close(file_id,error)
sll_int32, intent(in)  :: file_id
sll_int32, intent(out) :: error

write(file_id,"(a)")"</Grid>"
write(file_id,"(a)")"</Domain>"
write(file_id,"(a)")"</Xdmf>"
close(file_id)
error = 0
end subroutine sll_xmf_file_close

!> Write the description of a scalar field on a 2D mesh.
!> \param[in] file_id is the unit number or your xmf file
!> \param[in] filename is the file name where the heavy data are (bin or h5)
!> \param[in] nnodes_x1 - nodes number along direction 1
!> \param[in] nnodes_x2 - nodes number along direction 2
!> \param[in] filetype  - heavy data format 'HDF' or 'Binary'
!>
!> The file named filename must exist.
!>
subroutine sll_xmf_dataitem_2d(file_id, filename, nnodes_x1, nnodes_x2, filetype)
sll_int32, intent(in) :: file_id
character(len=*), intent(in) :: filename
character(len=*), intent(in) :: filetype
sll_int32, intent(in) :: nnodes_x1
sll_int32, intent(in) :: nnodes_x2

SLL_ASSERT(filetype == 'HDF' .or. filetype == 'Binary')
write(file_id,"(a,2i5,a)")"<DataItem Dimensions='",nnodes_x2,nnodes_x1, &
"' NumberType='Float' Precision='8' Format='"//trim(filetype)//"'>"
write(file_id,"(a)")trim(filename)
write(file_id,"(a)")"</DataItem>"
end subroutine sll_xmf_dataitem_2d

!> Write the description of a scalar field on a 2D mesh.
!> \param[in] file_id is the unit number or your xmf file
!> \param[in] fieldname is the dataset name where the heavy data are (hdf5 case)
!> \param[in] filename is the file name where the heavy data are (bin or h5)
!> \param[in] npoints_1 - nodes or cells number along direction 1
!> \param[in] npoints_2 - nodes or cells number along direction 2
!>
!> The file named filename-fieldname.bin must exist in case of binary output.
!> The file named filename.h5 with dataset fieldname must exist in case of hdf5 output.
!>
subroutine sll_xmf_field_2d(file_id,fieldname,filename,npoints_1,npoints_2,center)
sll_int32, intent(in) :: file_id
character(len=*), intent(in) :: filename
character(len=*), intent(in) :: fieldname
character(len=*), intent(in) :: center
sll_int32, intent(in) :: npoints_1
sll_int32, intent(in) :: npoints_2

write(file_id,"(a)")"<Attribute Name='"//fieldname//"' AttributeType='Scalar' Center='"//center//"'>"
#ifdef NOHDF5
call sll_xmf_dataitem_2d(file_id,trim(filename)//"-"//fieldname//".bin",npoints_1,npoints_2,'Binary')
#else
call sll_xmf_dataitem_2d(file_id,trim(filename)//".h5:/"//fieldname,npoints_1,npoints_2,'HDF')
#endif
write(file_id,"(a)")"</Attribute>"

end subroutine sll_xmf_field_2d


!> Write the description of a 2D strutured grid
!> mesh with its nodes coordinates contains in filename-x1 and filename-x2.
!> \param[in] file_id is the unit number or your xmf file
!> \param[in] filename is the file name where the coordinates data are (bin or h5)
!> \param[in] nnodes_x1 - nodes number along direction 1
!> \param[in] nnodes_x2 - nodes number along direction 2
!>
!> The file named filename-x1.bin and filename-x2.bin must exist in case of binary output.
!> The file named filename.h5 with dataset x1 and x2 must exist in case of hdf5 output.
!>
subroutine sll_xmf_grid_geometry_2d(file_id, filename, nnodes_x1, nnodes_x2)
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

end subroutine sll_xmf_grid_geometry_2d


end module sll_xmf_io
