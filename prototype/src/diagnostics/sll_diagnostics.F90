!> \brief Implements the functions to write xdmf file plotable by VisIt
!>
!> https://wci.llnl.gov/codes/visit/
!> http://www.xdmf.org/index.php/Main_Page
!> 
!> Two versions of output are available:
!> Data are stored in HDF5 file (http://www.hdfgroup.org/HDF5/) or in classic
!> binary file. This behaviour is control by the CPPDEFINE NOHDF5
!> HDF5 is the best choice but il you prefer binary just add 
!> CPPDEFINES = ['NOHDF5'] in your SCons environment
!> You can also do this with :
!> env.Append(CPPDEFINES=['NOHDF5']
!>
module sll_diagnostics
use sll_ascii_io
use sll_binary_io
use sll_xmf_io
#ifndef NOHDF5
use sll_hdf5_io
#endif

#include "sll_working_precision.h"
#include "sll_assert.h"

implicit none

private
enum, bind(C)
   enumerator :: NODE_CENTERED_DF = 0, CELL_CENTERED_DF = 1
end enum

public write_vec1d, write_vec2d, write_mesh

contains  ! ****************************************************************

#ifdef NOHDF5

!> Write XDMF file for 2D mesh plot
!> Data description is in the mesh_name.xmf file
!> Coordinates of mesh nodes are in mesh_name-x1.bin and mesh_name-x1.bin
!>\param[in] mesh_name filename prefix for the mesh
!>\param[in] x1 nodes coordinates  along direction 1
!>\param[in] x2 nodes coordinates  along direction 2
!>\param[in] nnodes_x1 nodes number along direction 1
!>\param[in] nnodes_x2 nodes number along direction 2
!> file_id file unit number
!> error  parameter that should be 0
subroutine write_mesh_data_2d(mesh_name,x1,x2,nnodes_x1,nnodes_x2,error)
character(len=*), intent(in) :: mesh_name 
sll_real64, dimension(:,:), intent(in) :: x1
sll_real64, dimension(:,:), intent(in) :: x2
sll_int32, intent(in) :: nnodes_x1
sll_int32, intent(in) :: nnodes_x2
sll_int32 :: xmffile_id, binfile_id
sll_int32, intent(out) :: error

SLL_ASSERT(size(x1,1) == nnodes_x1)
SLL_ASSERT(size(x2,1) == nnodes_x1)
SLL_ASSERT(size(x1,2) == nnodes_x2)
SLL_ASSERT(size(x2,2) == nnodes_x2)

call sll_xmf_file_create(trim(mesh_name)//".xmf",xmffile_id,error)
write(xmffile_id,"(a)")"<Grid Name='mesh' GridType='Uniform'>"
write(xmffile_id,"(a,2i5,a)")"<Topology TopologyType='2DSMesh' NumberOfElements='", &
                     nnodes_x2,nnodes_x1,"'/>"
write(xmffile_id,"(a)")"<Geometry GeometryType='X_Y'>"

#define WRITE_2D_ARRAY(array) \
call sll_binary_file_create(mesh_name//"-array.bin",binfile_id,error);\
SLL_ASSERT(error==0); \
call sll_binary_write_array_2d(binfile_id,array,error); \
SLL_ASSERT(error==0); \
call sll_binary_file_close(binfile_id,error); \
write(xmffile_id,"(a,2i5,a)")"<DataItem Dimensions='",nnodes_x2,nnodes_x1 \
                    ,"' NumberType='Float' Precision='8' Format='Binary'>"; \
write(xmffile_id,"(a)")mesh_name//"-array.bin";\
write(xmffile_id,"(a)")"</DataItem>"

WRITE_2D_ARRAY(x1)
WRITE_2D_ARRAY(x2)

write(xmffile_id,"(a)")"</Geometry>"
write(xmffile_id,"(a)")"</Grid>"
call sll_xmf_file_close(xmffile_id,error)

end subroutine write_mesh_data_2d

!> Write XDMF file for 3D mesh plot
!> Data description is in the mesh_name.xmf file
!> Coordinates of mesh nodes are in mesh_name-x1.bin, mesh_name-x1.bin and
!> mesh_name-x3.bin
!>\param[in] mesh_name filename prefix for the mesh
!>\param[in] x1 nodes coordinates  along direction 1
!>\param[in] x2 nodes coordinates  along direction 2
!>\param[in] x3 nodes coordinates  along direction 3
!>\param[in] nnodes_x1 nodes number along direction 1
!>\param[in] nnodes_x2 nodes number along direction 2
!>\param[in] nnodes_x3 nodes number along direction 3
!> file_id file unit number
!> error  parameter that should be 0
subroutine write_mesh_data_3d(mesh_name,x1,x2,x3,nnodes_x1,nnodes_x2,nnodes_x3,error)
character(len=*), intent(in) :: mesh_name 
sll_real64, dimension(:,:,:), intent(in) :: x1
sll_real64, dimension(:,:,:), intent(in) :: x2
sll_real64, dimension(:,:,:), intent(in) :: x3
sll_int32, intent(in) :: nnodes_x1
sll_int32, intent(in) :: nnodes_x2
sll_int32, intent(in) :: nnodes_x3
sll_int32, intent(out):: error
sll_int32 :: file_id

SLL_ASSERT(size(x1,1) == nnodes_x1)
SLL_ASSERT(size(x2,1) == nnodes_x1)
SLL_ASSERT(size(x3,1) == nnodes_x3)
SLL_ASSERT(size(x1,2) == nnodes_x2)
SLL_ASSERT(size(x2,2) == nnodes_x2)
SLL_ASSERT(size(x3,2) == nnodes_x2)
SLL_ASSERT(size(x1,3) == nnodes_x3)
SLL_ASSERT(size(x2,3) == nnodes_x3)
SLL_ASSERT(size(x3,3) == nnodes_x3)

call sll_xmf_file_create(trim(mesh_name)//".xmf",file_id,error)
write(file_id,"(a)")"<Grid Name='mesh' GridType='Uniform'>"
write(file_id,"(a,3i5,a)")"<Topology TopologyType='3DSMesh' NumberOfElements='", &
                     nnodes_x3,nnodes_x2,nnodes_x1,"'/>"
write(file_id,"(a)")"<Geometry GeometryType='X_Y'>"

#define WRITE_ARRAY(array) \
call sll_binary_file_create(mesh_name//"-array.bin",file_id,error); \
SLL_ASSERT(error==0); \
call sll_binary_write_array_3d(file_id,array,error); \
SLL_ASSERT(error==0); \
call sll_binary_file_close(file_id,error); \
SLL_ASSERT(error==0); \
write(file_id,"(a,3i5,a)")"<DataItem Dimensions='",nnodes_x3,nnodes_x2,nnodes_x1 \
                  ,"' NumberType='Float' Precision='8' Format='Binary'>"; \
write(file_id,"(a)")mesh_name//"-array.bin"; \
write(file_id,"(a)")"</DataItem>"

WRITE_ARRAY(x1)
WRITE_ARRAY(x2)
WRITE_ARRAY(x3)

write(file_id,"(a)")"</Geometry>"
write(file_id,"(a)")"</Grid>"
call sll_xmf_file_close(file_id,error)

end subroutine write_mesh_data_3d

!** In order to keep old behaviour
subroutine write_mesh(x1,x2,nnodes_x1,nnodes_x2,mesh_name)
character(len=*), intent(in) :: mesh_name 
sll_real64, dimension(:,:), intent(in) :: x1
sll_real64, dimension(:,:), intent(in) :: x2
sll_int32, intent(in) :: nnodes_x1
sll_int32, intent(in) :: nnodes_x2
sll_int32 :: error
    
call write_mesh_data_2d(mesh_name,x1,x2,nnodes_x1,nnodes_x2,error)

end subroutine write_mesh
    
!> Write XDMF file for vector 1d on 2d mesh with data in binary files
!> Data description is in the field_prefix.xmf file
!> Coordinates of mesh nodes are in mesh_name-x1.bin, mesh_name-x1.bin 
!> Vector field data are in field_prefix-vec1_nodes.bin or field_prefix-vec1_cells.bin
!>\param[in] Array with vector values
!>\param[in] nnodes_x1 nodes number along direction 1
!>\param[in] nnodes_x2 nodes number along direction 2
!>\param[in] field_prefix filename prefix for the vector
!>\param[in] mesh_name filename prefix for the mesh
!>\param[in] icenter parameter to distinguish cells values or nodes values
subroutine write_vec1d( vec_values, nnodes_x1, nnodes_x2, field_prefix, mesh_name, icenter)
character(len=*), intent(in) :: field_prefix !prefix for field file
character(len=*), intent(in) :: mesh_name !prefix for mesh file
sll_real64, dimension(:,:), intent(in) :: vec_values
sll_int32, intent(in) :: nnodes_x1
sll_int32, intent(in) :: nnodes_x2
sll_int32             :: ncells_x1
sll_int32             :: ncells_x2

sll_int32 :: file_id
sll_int32 :: data_dims(2)
integer :: error
logical :: flag
integer, intent(in), optional :: icenter   ! centering ('node' or 'cell')

SLL_ASSERT(size(vec_values,1) == nnodes_x1)
SLL_ASSERT(size(vec_values,2) == nnodes_x2)
ncells_x1 = nnodes_x1 - 1
ncells_x2 = nnodes_x2 - 1

inquire(file=trim(mesh_name)//".xmf", exist=flag) 

if (.not. flag) then
   call errout(6,"W","sll_diagnostics:write_vec1d","Mesh file does not exist")
end if

if (.not. present(icenter)) then
   call errout(6,"W","sll_diagnostics:write_vec1d","Default: cell centered value")
end if
    
call sll_xmf_file_create(trim(field_prefix)//".xmf",file_id,error)

write(file_id,"(a)")"<Grid Name='mesh' GridType='Uniform'>"
write(file_id,"(a,2i5,a)")"<Topology TopologyType='2DSMesh' NumberOfElements='", &
                          nnodes_x2,nnodes_x1,"'/>"
write(file_id,"(a)")"<Geometry GeometryType='X_Y'>"
write(file_id,"(a,2i5,a)")"<DataItem Dimensions='",nnodes_x2,nnodes_x1 &
                   ,"' NumberType='Float' Precision='8' Format='Binary'>"
write(file_id,"(a)")mesh_name//"-x1.bin"
write(file_id,"(a)")"</DataItem>"
write(file_id,"(a,2i5,a)")"<DataItem Dimensions='",nnodes_x2,nnodes_x1, &
                   "' NumberType='Float' Precision='8' Format='Binary'>"
write(file_id,"(a)")mesh_name//"-x2.bin"
write(file_id,"(a)")"</DataItem>"
write(file_id,"(a)")"</Geometry>"


if (present(icenter) .and. icenter == NODE_CENTERED_DF) then

   data_dims(1) = nnodes_x1
   data_dims(2) = nnodes_x2

   write(file_id,"(a)")"<Attribute Name='NodesValues' AttributeType='Scalar' Center='Node'>"
   write(file_id,"(a,2i5,a)")"<DataItem Dimensions='",data_dims(2),data_dims(1),              &
                        "' NumberType='Float' Precision='8' Format='Binary'>"
   write(file_id,"(a)")field_prefix//".bin"
   write(file_id,"(a)")"</DataItem>"
   write(file_id,"(a)")"</Attribute>"

else

   data_dims(1) = nnodes_x1-1
   data_dims(2) = nnodes_x2-1

   write(file_id,"(a)")"<Attribute Name='CellsValues' AttributeType='Scalar' Center='Cell'>"
   write(file_id,"(a,2i5,a)")"<DataItem Dimensions='",data_dims(2),data_dims(1),             &
                        "' NumberType='Float' Precision='8' Format='Binary'>"
   write(file_id,"(a)")field_prefix//".bin"
   write(file_id,"(a)")"</DataItem>"
   write(file_id,"(a)")"</Attribute>"

end if

write(file_id,"(a)")"</Grid>"
call sll_xmf_file_close(file_id,error)

if (present(icenter) .and. icenter == NODE_CENTERED_DF) then

   call sll_binary_file_create(trim(field_prefix)//".bin",file_id,error)
   SLL_ASSERT(error==0)
   call sll_binary_write_array_2d(file_id,vec_values,error)
   SLL_ASSERT(error==0)
   call sll_binary_file_close(file_id,error)
   SLL_ASSERT(error==0)

else

   call sll_binary_file_create(trim(field_prefix)//".bin",file_id,error)
   SLL_ASSERT(error==0)
   call sll_binary_write_array_2d(file_id,vec_values(1:ncells_x1,1:ncells_x2),error)
   SLL_ASSERT(error==0)
   call sll_binary_file_close(file_id,error)
   SLL_ASSERT(error==0)

end if
    
end subroutine write_vec1d
  
!> Write XDMF file for vector 2d on 2d mesh with data in binary files
!> Data description is in the field_prefix.xmf file
!> Coordinates of mesh nodes are in mesh_name-x1.bin, mesh_name-x1.bin
!> Vector field data are in field_prefix-vec1_nodes.bin or field_prefix-vec1_cells.bin
!>\param[in] vec_values_x1 Array with vector values
!>\param[in] vec_values_x2 Array with vector values
!>\param[in] nnodes_x1 nodes number along direction 1
!>\param[in] nnodes_x2 nodes number along direction 2
!>\param[in] field_prefix filename prefix for the vector
!>\param[in] mesh_name filename prefix for the mesh
!>\param[in] icenter parameter to distinguish cells values or nodes values
subroutine write_vec2d( vec_values_x1, vec_values_x2,       &
                        nnodes_x1, nnodes_x2,               &
                        field_prefix, mesh_prefix, icenter)

character(len=*), intent(in) :: field_prefix !prefix for field file
character(len=*), intent(in) :: mesh_prefix  !prefix for mesh file
sll_real64, dimension(:,:), intent(in) :: vec_values_x1
sll_real64, dimension(:,:), intent(in) :: vec_values_x2
sll_int32, intent(in) :: nnodes_x1
sll_int32, intent(in) :: nnodes_x2
sll_int32             :: ncells_x1
sll_int32             :: ncells_x2
integer, intent(in), optional :: icenter   ! centering of field, one of ('node' or 'cell')

sll_int32        :: file_id
sll_int32        :: data_dims(2)
sll_int32        :: error
logical          :: flag

SLL_ASSERT(size(vec_values_x1,1) == nnodes_x1)
SLL_ASSERT(size(vec_values_x1,2) == nnodes_x2)
SLL_ASSERT(size(vec_values_x2,1) == nnodes_x1)
SLL_ASSERT(size(vec_values_x2,2) == nnodes_x2)

ncells_x1 = nnodes_x1-1
ncells_x2 = nnodes_x2-1

inquire(file=trim(mesh_prefix)//".xmf", exist=flag) 

if (.not. flag) then
   call errout(6,"W","sll_diagnostics:write_vec2d", "Mesh file does not exist" )
end if
if (.not. present(icenter)) then
   call errout(6,"W","sll_diagnostics:write_vec2d","Default: cell centered value")
end if
    
call sll_xmf_file_create(trim(field_prefix)//".xmf",file_id,error)
write(file_id,"(a)")"<Grid Name='mesh' GridType='Uniform'>"
write(file_id,"(a,2i5,a)")"<Topology TopologyType='2DSMesh' NumberOfElements='",nnodes_x2,nnodes_x1,"'/>"
write(file_id,"(a)")"<Geometry GeometryType='X_Y'>"
call sll_xmf_dataitem_2d(file_id, trim(mesh_prefix)//"-x1.bin", nnodes_x1, nnodes_x2)
call sll_xmf_dataitem_2d(file_id, trim(mesh_prefix)//"-x2.bin", nnodes_x1, nnodes_x2)
write(file_id,"(a)")"</Geometry>"

if (present(icenter) .and. icenter == NODE_CENTERED_DF) then

   data_dims(1) = nnodes_x1
   data_dims(2) = nnodes_x2

   write(file_id,"(a)")"<Attribute Name='NodesValues1' AttributeType='Scalar' Center='Node'>"
   call sll_xmf_dataitem_2d(file_id,trim(field_prefix)//"-x1.bin",nnodes_x1, nnodes_x2)
   write(file_id,"(a)")"</Attribute>"
   write(file_id,"(a)")"<Attribute Name='NodesValues2' AttributeType='Scalar' Center='Node'>"
   call sll_xmf_dataitem_2d(file_id,trim(field_prefix)//"-x2.bin",nnodes_x1, nnodes_x2)
   write(file_id,"(a)")"</Attribute>"

else
    
   data_dims(1) = ncells_x1
   data_dims(2) = ncells_x2

   write(file_id,"(a)")"<Attribute Name='CellsValues1' AttributeType='Scalar' Center='Cell'>"
   call sll_xmf_dataitem_2d(file_id,trim(field_prefix)//"-x1.bin",ncells_x1, ncells_x2)
   write(file_id,"(a)")"</Attribute>"
   write(file_id,"(a)")"<Attribute Name='CellsValues2' AttributeType='Scalar' Center='Cell'>"
   call sll_xmf_dataitem_2d(file_id,trim(field_prefix)//"-x2.bin",ncells_x1, ncells_x2)
   write(file_id,"(a)")"</Attribute>"

end if
    
write(file_id,"(a)")"</Grid>"
call sll_xmf_file_close(file_id,error)

if (present(icenter) .and. icenter == NODE_CENTERED_DF) then

   call sll_binary_file_create(trim(field_prefix)//"-x1.bin",file_id,error)
   SLL_ASSERT(error==0)
   call sll_binary_write_array_2d(file_id,vec_values_x1,error)
   SLL_ASSERT(error==0)
   call sll_binary_file_close(file_id,error)
   SLL_ASSERT(error==0)

   call sll_binary_file_create(trim(field_prefix)//"-x2.bin",file_id,error)
   SLL_ASSERT(error==0)
   call sll_binary_write_array_2d(file_id,vec_values_x2,error)
   SLL_ASSERT(error==0)
   call sll_binary_file_close(file_id,error)
   SLL_ASSERT(error==0)

else

   call sll_binary_file_create(trim(field_prefix)//"-x1.bin",file_id,error)
   SLL_ASSERT(error==0)
   call sll_binary_write_array_2d(file_id,vec_values_x1(1:ncells_x1,1:ncells_x2),error)
   SLL_ASSERT(error==0)
   call sll_binary_file_close(file_id,error)
   SLL_ASSERT(error==0)

   call sll_binary_file_create(trim(field_prefix)//"-x2.bin",file_id,error)
   SLL_ASSERT(error==0)
   call sll_binary_write_array_2d(file_id,vec_values_x2(1:ncells_x1,1:ncells_x2),error)
   SLL_ASSERT(error==0)
   call sll_binary_file_close(file_id,error)
   SLL_ASSERT(error==0)

end if

end subroutine write_vec2d

#else

subroutine write_mesh( &
   x1,                 &
   x2,                 &
   nnodes_x1,          &
   nnodes_x2,          &
   mesh_name)

character(len=*), intent(in) :: mesh_name 
sll_real64, dimension(:,:), intent(in) :: x1
sll_real64, dimension(:,:), intent(in) :: x2

sll_int32, intent(in) :: nnodes_x1
sll_int32, intent(in) :: nnodes_x2
   
integer(hid_t)   :: file_id
integer(hsize_t) :: data_dims(2)
character(len=3) :: coord_names(2)
integer :: error

    
SLL_ASSERT(size(x1,1) == nnodes_x1)
SLL_ASSERT(size(x2,1) == nnodes_x1)
SLL_ASSERT(size(x1,2) == nnodes_x2)
SLL_ASSERT(size(x2,2) == nnodes_x2)
    
!Initialize FORTRAN interface.
!Create a new file using default properties.
call sll_hdf5_file_create(trim(mesh_name)//".h5",file_id,error)
SLL_ASSERT(error==0)

!Write the data file.
coord_names(1) = "/x1"
coord_names(2) = "/x2"
    
!Write separate coordinate arrays for the x and v coordinates.
data_dims(1) = nnodes_x1
data_dims(2) = nnodes_x2

call sll_hdf5_write_array_2d(file_id,x1,coord_names(1),error)
SLL_ASSERT(error==0)
call sll_hdf5_write_array_2d(file_id,x2,coord_names(2),error)
SLL_ASSERT(error==0)

!Terminate access to the file.
call sll_hdf5_file_close(file_id, error)
SLL_ASSERT(error==0)

!write the xmf file readable by VisIt
call sll_xmf_file_create(trim(mesh_name)//".xmf",file_id,error)
write(file_id,"(a)")"<Grid Name='mesh' GridType='Uniform'>"
write(file_id,"(a,2i5,a)")"<Topology TopologyType='2DSMesh' NumberOfElements='", &
                     data_dims(2),data_dims(1),"'/>"
write(file_id,"(a)")"<Geometry GeometryType='X_Y'>"
write(file_id,"(a,2i5,a)")"<DataItem Dimensions='",data_dims(2),data_dims(1) &
                    ,"' NumberType='Float' Precision='8' Format='HDF'>"
write(file_id,"(a)")trim(mesh_name)//".h5:"//coord_names(1)
write(file_id,"(a)")"</DataItem>"
write(file_id,"(a,2i5,a)")"<DataItem Dimensions='",data_dims(2),data_dims(1) &
                    ,"' NumberType='Float' Precision='8' Format='HDF'>"
write(file_id,"(a)")trim(mesh_name)//".h5:"//coord_names(2)
write(file_id,"(a)")"</DataItem>"
write(file_id,"(a)")"</Geometry>"
write(file_id,"(a)")"</Grid>"
call sll_xmf_file_close(file_id,error)

end subroutine write_mesh
    
subroutine write_vec1d( vec_values, nnodes_x1, nnodes_x2, field_prefix, mesh_name, icenter)
character(len=*), intent(in) :: field_prefix !prefix for field file
character(len=*), intent(in) :: mesh_name !prefix for mesh file
sll_real64, dimension(:,:), intent(in) :: vec_values
sll_int32, intent(in) :: nnodes_x1
sll_int32, intent(in) :: nnodes_x2
sll_int32             :: ncells_x1
sll_int32             :: ncells_x2

integer(hid_t)   :: file_id
integer(hsize_t) :: data_dims(2)
character(len=3) :: coord_names(2)
integer :: error
logical :: flag
integer, intent(in), optional :: icenter   ! centering ('node' or 'cell')

SLL_ASSERT(size(vec_values,1) == nnodes_x1)
SLL_ASSERT(size(vec_values,2) == nnodes_x2)
ncells_x1 = nnodes_x1 - 1
ncells_x2 = nnodes_x2 - 1

inquire(file=trim(mesh_name)//".h5", exist=flag) 

if (.not. flag) then
   call errout(6,"W","sll_diagnostics:write_vec1d","Mesh file does not exist")
end if

if (.not. present(icenter)) then
   call errout(6,"W","sll_diagnostics:write_vec1d","Default: cell centered value")
end if

coord_names(1) = "/x1"
coord_names(2) = "/x2"

call sll_xmf_file_create(trim(field_prefix)//".xmf",file_id,error)
write(file_id,"(a)")"<Grid Name='mesh' GridType='Uniform'>"
write(file_id,"(a,2i5,a)")"<Topology TopologyType='2DSMesh' NumberOfElements='",nnodes_x2,nnodes_x1,"'/>"
write(file_id,"(a)")"<Geometry GeometryType='X_Y'>"
write(file_id,"(a,2i5,a)")"<DataItem Dimensions='",nnodes_x2,nnodes_x1,&
"' NumberType='Float' Precision='8' Format='HDF'>"
write(file_id,"(a)")trim(mesh_name)//".h5:"//coord_names(1)
write(file_id,"(a)")"</DataItem>"
write(file_id,"(a,2i5,a)")"<DataItem Dimensions='",nnodes_x2,nnodes_x1,&
"' NumberType='Float' Precision='8' Format='HDF'>"
write(file_id,"(a)")trim(mesh_name)//".h5:"//coord_names(2)
write(file_id,"(a)")"</DataItem>"
write(file_id,"(a)")"</Geometry>"

if (present(icenter) .and. icenter == NODE_CENTERED_DF) then

   data_dims(1) = nnodes_x1
   data_dims(2) = nnodes_x2
   
   write(file_id,"(a)")"<Attribute Name='NodesValues' AttributeType='Scalar' Center='Node'>"
   write(file_id,"(a,2i5,a)")"<DataItem Dimensions='",data_dims(2),data_dims(1),              &
                        "' NumberType='Float' Precision='8' Format='HDF'>"
   write(file_id,"(a)")trim(field_prefix)//".h5:/vec1_nodes"
   write(file_id,"(a)")"</DataItem>"
   write(file_id,"(a)")"</Attribute>"
   
else

   data_dims(1) = nnodes_x1-1
   data_dims(2) = nnodes_x2-1

   write(file_id,"(a)")"<Attribute Name='CellsValues' AttributeType='Scalar' Center='Cell'>"
   write(file_id,"(a,2i5,a)")"<DataItem Dimensions='",data_dims(2),data_dims(1),              &
                        "' NumberType='Float' Precision='8' Format='HDF'>"
   write(file_id,"(a)")trim(field_prefix)//".h5:/vec1_cells"
   write(file_id,"(a)")"</DataItem>"
   write(file_id,"(a)")"</Attribute>"

end if
write(file_id,"(a)")"</Grid>"
call sll_xmf_file_close(file_id,error)

!write data
call sll_hdf5_file_create(trim(field_prefix)//".h5",file_id,error)
SLL_ASSERT(error==0)
if (present(icenter) .and. icenter == NODE_CENTERED_DF) then
   call sll_hdf5_write_array_2d(file_id,vec_values,"vec1_nodes",error)
   SLL_ASSERT(error==0)
else
   call sll_hdf5_write_array_2d(file_id,vec_values(1:ncells_x1,1:ncells_x2),"vec1_cells",error)
   SLL_ASSERT(error==0)
end if
call sll_hdf5_file_close(file_id,error)
SLL_ASSERT(error==0)

    
end subroutine write_vec1d
  
subroutine write_vec2d( &
   vec_values_x1,       &
   vec_values_x2,       &
   nnodes_x1,           &
   nnodes_x2,           &
   field_prefix,     &
   mesh_prefix,     &
   icenter)

character(len=*), intent(in) :: field_prefix !prefix for field file
character(len=*), intent(in) :: mesh_prefix  !prefix for mesh file
sll_real64, dimension(:,:), intent(in) :: vec_values_x1
sll_real64, dimension(:,:), intent(in) :: vec_values_x2
sll_int32, intent(in) :: nnodes_x1
sll_int32, intent(in) :: nnodes_x2
sll_int32             :: ncells_x1
sll_int32             :: ncells_x2
integer, intent(in), optional :: icenter   ! centering of field, one of ('node' or 'cell')

integer(hid_t)   :: file_id
integer(hsize_t) :: data_dims(2)
character(len=3) :: coord_names(2)
sll_int32        :: error
logical          :: flag

SLL_ASSERT(size(vec_values_x1,1) == nnodes_x1)
SLL_ASSERT(size(vec_values_x1,2) == nnodes_x2)
SLL_ASSERT(size(vec_values_x2,1) == nnodes_x1)
SLL_ASSERT(size(vec_values_x2,2) == nnodes_x2)

ncells_x1 = nnodes_x1-1
ncells_x2 = nnodes_x2-1

inquire(file=trim(mesh_prefix)//".h5", exist=flag) 

if (.not. flag) then
   call errout(6,"W","sll_diagnostics:write_vec2d", "Mesh file does not exist" )
end if
if (.not. present(icenter)) then
   call errout(6,"W","sll_diagnostics:write_vec2d","Default: cell centered value")
end if
    
coord_names(1) = "/x1"
coord_names(2) = "/x2"

call sll_xmf_file_create(trim(field_prefix)//".xmf",file_id,error)
write(file_id,"(a)")"<Grid Name='mesh' GridType='Uniform'>"
write(file_id,"(a,2i5,a)")"<Topology TopologyType='2DSMesh' NumberOfElements='",nnodes_x2,nnodes_x1,"'/>"
write(file_id,"(a)")"<Geometry GeometryType='X_Y'>"
write(file_id,"(a,2i5,a)")"<DataItem Dimensions='",nnodes_x2,nnodes_x1, &
                          "' NumberType='Float' Precision='8' Format='HDF'>"
write(file_id,"(a)")trim(mesh_prefix)//".h5:"//coord_names(1)
write(file_id,"(a)")"</DataItem>"
write(file_id,"(a,2i5,a)")"<DataItem Dimensions='",nnodes_x2,nnodes_x1, &
                          "' NumberType='Float' Precision='8' Format='HDF'>"
write(file_id,"(a)")trim(mesh_prefix)//".h5:"//coord_names(2)
write(file_id,"(a)")"</DataItem>"
write(file_id,"(a)")"</Geometry>"

if (present(icenter) .and. icenter == NODE_CENTERED_DF) then

   write(file_id,"(a)")"<Attribute Name='x_NodesValues' AttributeType='Scalar' Center='Node'>"
   write(file_id,"(a,2i5,a)")"<DataItem Dimensions='",nnodes_x2,nnodes_x1, &
                        "' NumberType='Float' Precision='8' Format='HDF'>"
   write(file_id,"(a)")trim(field_prefix)//".h5:/vec2_x1_nodes"
   write(file_id,"(a)")"</DataItem>"
   write(file_id,"(a)")"</Attribute>"
   write(file_id,"(a)")"<Attribute Name='v_NodesValues' AttributeType='Scalar' Center='Node'>"
   write(file_id,"(a,2i5,a)")"<DataItem Dimensions='",nnodes_x2,nnodes_x1, &
                        "' NumberType='Float' Precision='8' Format='HDF'>"
   write(file_id,"(a)")trim(field_prefix)//".h5:/vec2_x2_nodes"
   write(file_id,"(a)")"</DataItem>"
   write(file_id,"(a)")"</Attribute>"

else

   write(file_id,"(a)")"<Attribute Name='x_CellsValues' AttributeType='Scalar' Center='Cell'>"
   write(file_id,"(a,2i5,a)")"<DataItem Dimensions='",ncells_x2,ncells_x1,"' NumberType='Float' Precision='8' Format='HDF'>"
   write(file_id,"(a)")trim(field_prefix)//".h5:/vec2_x1_cells"
   write(file_id,"(a)")"</DataItem>"
   write(file_id,"(a)")"</Attribute>"
   write(file_id,"(a)")"<Attribute Name='v_CellsValues' AttributeType='Scalar' Center='Cell'>"
   write(file_id,"(a,2i5,a)")"<DataItem Dimensions='",ncells_x2,ncells_x1,"' NumberType='Float' Precision='8' Format='HDF'>"
   write(file_id,"(a)")trim(field_prefix)//".h5:/vec2_x2_cells"
   write(file_id,"(a)")"</DataItem>"
   write(file_id,"(a)")"</Attribute>"
    
end if

write(file_id,"(a)")"</Grid>"
call sll_xmf_file_close(file_id,error)
    
call sll_hdf5_file_create(trim(field_prefix)//".h5",file_id,error)
SLL_ASSERT(error==0)

if (present(icenter) .and. icenter == NODE_CENTERED_DF) then

   call sll_hdf5_write_array_2d(file_id,vec_values_x1,"vec2_x1_nodes",error)
   SLL_ASSERT(error==0)

   call sll_hdf5_write_array_2d(file_id,vec_values_x2,"vec2_x2_nodes",error)
   SLL_ASSERT(error==0)

else

   call sll_hdf5_write_array_2d(file_id,vec_values_x1(1:ncells_x1,1:ncells_x2),"vec2_x1_cells", error)
   SLL_ASSERT(error==0)

   call sll_hdf5_write_array_2d(file_id,vec_values_x2(1:ncells_x1,1:ncells_x2),"vec2_x2_cells", error)
   SLL_ASSERT(error==0)
end if

call sll_hdf5_file_close(file_id,error)
SLL_ASSERT(error==0)

end subroutine write_vec2d

#endif
  

!Title: Subroutine Errout
!Subroutine: errout
!     Outputs an error message
!     PRTFIL - unit number for print-out
!     SEVRTY - 'W' - Warning 'F' - Fatal
!     WHERE  - in which program or subroutine
!     ErrMsg - error message
subroutine errout( prtfil, sevrty, lwhere, ErrMsg )
implicit none
integer, intent(in) ::  prtfil
character(len=1) :: sevrty 
character(len=*) :: lwhere , ErrMsg

write( prtfil, * )
select case ( sevrty )  !     *** Severity ***
case ( 'W' )
  write(prtfil,"(/10x,a)") '*** WARNING ***'
case ( 'F' )
  write(prtfil,"(/10x,a)") '*** FATAL ERROR ***'
case default
  write(prtfil,"(/10x,a)") '*** FATAL ERROR ***'
  write(prtfil,"(/10x,a)") 'Error handler (ERROUT) called with unknown severity level: ', SEVRTY
end select
write( prtfil,"(/10x,a)") 'Generated by program or subroutine: ', trim(lwhere)
write( prtfil,"(/10x,a)") trim(ErrMsg)
write( prtfil,"(/10x,a)")

! return or stop depending on severity

if ( sevrty == 'W' ) then
   return
else
   stop 'Fatal Error: See print file for details'
end if

end subroutine errout

end module sll_diagnostics
