!------------------------------------------------------------------------------
! SELALIB
!------------------------------------------------------------------------------
!
! MODULE: sll_diagnostics
!
!> @author
!> Pierre Navaro
!>
!
! DESCRIPTION: 
!
!> @brief
!> Implements the functions to write xdmf file plotable by VisIt
!>
!>@details
!> In <b> XDMF </b> (eXtensible Data Model and Format) the description of the data 
!> is separate from the values themselves. Light data is stored using XML, Heavy data 
!> is stored using :
!> - HDF5 file 
!> - Binary file. 
!>
!> This is control by the variable <code>NOHDF5</code>.
!> HDF5 is set by default but il you prefer binary just add 
!>
!> <code> env.Append(CPPDEFINES=['NOHDF5']) </code>
!>
!> in your SCons script
!>
!> <h2>How to use this module: </h2>
!>
!> <p> Just use the module \a sll_diagnostics 
!> \code use sll_diagnostics \endcode
!>
!> External links:
!> - https://wci.llnl.gov/codes/visit/
!> - http://www.xdmf.org/index.php/Main_Page
!> - HDF5 file (http://www.hdfgroup.org/HDF5/)
!
! REVISION HISTORY:
! DD Mmm YYYY - Initial Version
! TODO_dd_mmm_yyyy - TODO_describe_appropriate_changes - TODO_name
!------------------------------------------------------------------------------
module sll_diagnostics
#ifndef NOHDF5
use hdf5
use sll_hdf5_io
#else
use sll_binary_io
#endif
use sll_ascii_io
use sll_xmf_io

#include "sll_working_precision.h"
#include "sll_assert.h"

implicit none

private
enum, bind(C)
   enumerator :: NODE_CENTERED_DF = 0, CELL_CENTERED_DF = 1
end enum

public write_vec1d, write_vec2d, write_mesh

contains  
!>Writes data for mesh plotting:
!> \param[in] mesh_prefix is the filename prefix for the mesh
!> \param[in] x1 is array with nodes coordinates  along direction 1
!> \param[in] x2 is array with nodes coordinates  along direction 2
!> \param[in] nnodes_x1 is the number of nodes number along direction 1
!> \param[in] nnodes_x2 is the number of nodes number along direction 2
!> \param[out] error is a parameter that should be 0
subroutine write_mesh_data_2d(mesh_prefix,x1,x2,nnodes_x1,nnodes_x2,error)
character(len=*), intent(in) :: mesh_prefix  
sll_real64, dimension(:,:), intent(in) :: x1 
sll_real64, dimension(:,:), intent(in) :: x2 
sll_int32, intent(in)  :: nnodes_x1          
sll_int32, intent(in)  :: nnodes_x2          
sll_int32, intent(out) :: error              

sll_int32 :: xmffile_id

#ifndef NOHDF5
integer(hid_t)   :: hdffile_id
integer(hsize_t) :: data_dims(2)
character(len=3) :: coord_names(2)
#else
sll_int32 :: binfile_id 
#endif

SLL_ASSERT(size(x1,1) == nnodes_x1)
SLL_ASSERT(size(x2,1) == nnodes_x1)
SLL_ASSERT(size(x1,2) == nnodes_x2)
SLL_ASSERT(size(x2,2) == nnodes_x2)

call sll_xmf_file_create(trim(mesh_prefix)//".xmf",xmffile_id,error)
call sll_xmf_grid_geometry_2d(xmffile_id, trim(mesh_prefix), nnodes_x1, nnodes_x2)

#ifdef NOHDF5

call sll_binary_file_create(mesh_prefix//"-x1.bin",binfile_id,error); SLL_ASSERT(error==0)
call sll_binary_write_array_2d(binfile_id,x1,error); SLL_ASSERT(error==0)
call sll_binary_file_close(binfile_id,error)
call sll_binary_file_create(mesh_prefix//"-x2.bin",binfile_id,error); SLL_ASSERT(error==0)
call sll_binary_write_array_2d(binfile_id,x2,error); SLL_ASSERT(error==0)
call sll_binary_file_close(binfile_id,error); SLL_ASSERT(error==0)

#else

call sll_hdf5_file_create(trim(mesh_prefix)//".h5",hdffile_id,error)
SLL_ASSERT(error==0)

coord_names(1) = "/x1"
coord_names(2) = "/x2"

data_dims(1) = nnodes_x1
data_dims(2) = nnodes_x2

call sll_hdf5_write_array_2d(hdffile_id,x1,coord_names(1),error); SLL_ASSERT(error==0)
call sll_hdf5_write_array_2d(hdffile_id,x2,coord_names(2),error); SLL_ASSERT(error==0)
call sll_hdf5_file_close(hdffile_id, error);                      SLL_ASSERT(error==0)

#endif

call sll_xmf_file_close(xmffile_id,error)

end subroutine write_mesh_data_2d

subroutine write_mesh_data_3d(mesh_prefix,x1,x2,x3,nnodes_x1,nnodes_x2,nnodes_x3,error)
character(len=*), intent(in) :: mesh_prefix 
sll_real64, dimension(:,:,:), intent(in) :: x1
sll_real64, dimension(:,:,:), intent(in) :: x2
sll_real64, dimension(:,:,:), intent(in) :: x3
sll_int32, intent(in) :: nnodes_x1
sll_int32, intent(in) :: nnodes_x2
sll_int32, intent(in) :: nnodes_x3
sll_int32, intent(out):: error
sll_int32 :: xmffile_id

#ifndef NOHDF5
sll_int32 :: i
character(len=3), dimension(3) :: coord_names
integer(hsize_t), dimension(3) :: data_dims
integer(hid_t)   :: hdffile_id
#else
sll_int32 :: binfile_id
#endif

SLL_ASSERT(size(x1,1) == nnodes_x1)
SLL_ASSERT(size(x2,1) == nnodes_x2)
SLL_ASSERT(size(x3,1) == nnodes_x3)
SLL_ASSERT(size(x1,2) == nnodes_x1)
SLL_ASSERT(size(x2,2) == nnodes_x2)
SLL_ASSERT(size(x3,2) == nnodes_x3)
SLL_ASSERT(size(x1,3) == nnodes_x1)
SLL_ASSERT(size(x2,3) == nnodes_x2)
SLL_ASSERT(size(x3,3) == nnodes_x3)

call sll_xmf_file_create(trim(mesh_prefix)//".xmf",xmffile_id,error)
write(xmffile_id,"(a)")"<Grid Name='mesh' GridType='Uniform'>"
write(xmffile_id,"(a,3i5,a)")"<Topology TopologyType='3DSMesh' NumberOfElements='", &
                     nnodes_x3,nnodes_x2,nnodes_x1,"'/>"
write(xmffile_id,"(a)")"<Geometry GeometryType='X_Y_Z'>"

#define WRITE_BIN_ARRAY(array) \
call sll_binary_file_create(mesh_prefix//"-array.bin",binfile_id,error); \
SLL_ASSERT(error==0); \
call sll_binary_write_array_3d(binfile_id,array,error); \
SLL_ASSERT(error==0); \
call sll_binary_file_close(binfile_id,error); \
SLL_ASSERT(error==0); \
write(xmffile_id,"(a,3i5,a)")"<DataItem Dimensions='",nnodes_x3,nnodes_x2,nnodes_x1 \
                  ,"' NumberType='Float' Precision='8' Format='Binary'>"; \
write(xmffile_id,"(a)")mesh_prefix//"-array.bin"; \
write(xmffile_id,"(a)")"</DataItem>"

#ifdef NOHDF5

WRITE_BIN_ARRAY(x1)
WRITE_BIN_ARRAY(x2)
WRITE_BIN_ARRAY(x3)

#else

coord_names(1) = "/x1"
coord_names(2) = "/x2"
coord_names(3) = "/x3"
    
data_dims(1) = nnodes_x1
data_dims(2) = nnodes_x2
data_dims(3) = nnodes_x3

call sll_hdf5_file_create(trim(mesh_prefix)//".h5",hdffile_id,error); SLL_ASSERT(error==0)
call sll_hdf5_write_array_3d(hdffile_id,x1,coord_names(1),error);   SLL_ASSERT(error==0)
call sll_hdf5_write_array_3d(hdffile_id,x2,coord_names(2),error);   SLL_ASSERT(error==0)
call sll_hdf5_write_array_3d(hdffile_id,x3,coord_names(3),error);   SLL_ASSERT(error==0)
call sll_hdf5_file_close(hdffile_id, error);                        SLL_ASSERT(error==0)

do i = 1, 3
   write(xmffile_id,"(a,3i5,a)")"<DataItem Dimensions='",nnodes_x3,nnodes_x2,nnodes_x1 &
                     ,"' NumberType='Float' Precision='8' Format='HDF'>"
   write(xmffile_id,"(a)")mesh_prefix//".h5:"//coord_names(i)
   write(xmffile_id,"(a)")"</DataItem>"
end do

#endif

write(xmffile_id,"(a)")"</Geometry>"
call sll_xmf_file_close(xmffile_id,error)

end subroutine write_mesh_data_3d

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> Write XDMF file for 2D mesh plot
!!
!! Light data are in the file:
!!  - <c>mesh_prefix.xmf</c>.
!!
!! Heavy data (mesh nodes coordinates) are in :
!!  - <c>mesh_prefix.h5</c>
!!  - <c>mesh_prefix-x1.bin</c> and <c>mesh_prefix-x2.bin</c>
!!
!! \param[in] x1 array with first coordinates
!! \param[in] x2 array with second coordinates
!! \param[in] nnodes_x1 nodes number along direction 1
!! \param[in] nnodes_x2 nodes number along direction 2
!! \param[in] mesh_prefix filename prefix for the mesh
!!
!! Example of a cartesian grid (xmin,xmax) x (vmin,vmax):
!!\verbatim
!!program grid
!!#include "sll_memory.h"
!!use sll_diagnostics
!!sll_int64, dimension(:,:), allocatable :: x, v
!!sll_int32 :: nx, nv
!!
!!nx = 32; xmin = 0.; xmax = 1.
!!nv = 64; vmin = 0.; vmax = 1.
!!SLL_ALLOCATE(x(nx,nv))
!!SLL_ALLOCATE(v(nx,nv)) 
!!do j = 1, nv                                             
!!   do i = 1, nx                                         
!!      x(i,j) = float(i-1) * (xmax-xmin) / (nx-1)       
!!      v(i,j) = float(j-1) * (vmax-vmin) / (nv-1)      
!!   end do                                            
!!end do                                              
!!call write_mesh(x,v,nx,nv,"mesh")
!!end program
!!\endverbatim
!! A file <i>mesh.xmf</i> is created, readable by VisIt or Paraview. For non
!!moving mesh, call it only one time, heavy data produced will be reused.
!!
subroutine write_mesh(x1,x2,nnodes_x1,nnodes_x2,mesh_prefix)
character(len=*), intent(in) :: mesh_prefix 
sll_real64, dimension(:,:), intent(in) :: x1
sll_real64, dimension(:,:), intent(in) :: x2
sll_int32, intent(in) :: nnodes_x1
sll_int32, intent(in) :: nnodes_x2
sll_int32 :: error
    
call write_mesh_data_2d(mesh_prefix,x1,x2,nnodes_x1,nnodes_x2,error)

end subroutine write_mesh
    
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> Write XDMF file for vector 1d on 2d mesh with data in binary or HDF5 files
!!
!! Light data are in :
!!  - <c>field_prefix.xmf</c>.
!! Heavy data :
!!  - mesh nodes coordinates are in :
!!     - <c>mesh_prefix.h5</c>
!!     - <c>mesh_prefix-x1.bin</c> and <c>mesh_prefix-x2.bin</c>
!!  - scalar values are in :
!!     - <c>field_prefix.h5</c>
!!     - <c>field_prefix.bin</c>
!!
!! \param[in] vec_values_x1 Array with vector values
!! \param[in] nnodes_x1 nodes number along direction 1
!! \param[in] nnodes_x2 nodes number along direction 2
!! \param[in] field_prefix filename prefix for the vector
!! \param[in] mesh_prefix filename prefix for the mesh
!! \param[in] icenter parameter to distinguish cells (0) or nodes (1) values
!!
!! Example of f(x,v) = xv on a cartesian grid (xmin,xmax) x (vmin,vmax):
!!\verbatim
!!program grid
!!#include "sll_memory.h"
!!#include "sll_working_precision.h"
!!use sll_diagnostics
!!sll_int64, dimension(:,:), allocatable :: x, v
!!sll_int32 :: nx, nv, error
!!
!!nx = 32; xmin = 0.; xmax = 1.
!!nv = 64; vmin = 0.; vmax = 1.
!!SLL_ALLOCATE(x(nx,nv),error)
!!SLL_ALLOCATE(v(nx,nv),error) 
!!SLL_ALLOCATE(f(nx,nv),error) 
!!do j = 1, nv                                             
!!   do i = 1, nx                                         
!!      x(i,j) = float(i-1) * (xmax-xmin) / (nx-1)       
!!      v(i,j) = float(j-1) * (vmax-vmin) / (nv-1)      
!!      f(i,j) = x(i,j) * v(i,j)
!!   end do                                            
!!end do                                              
!!call write_mesh(x,v,nx,nv,"mesh")
!!call write_vec1d(f,nx,nv,"f","mesh",0) !write cell centered values
!!end program
!!\endverbatim
!! A file <i>f.xmf</i> is created, readable by VisIt or Paraview. 
!!
subroutine write_vec1d(vec_values,nnodes_x1,nnodes_x2,field_prefix,mesh_prefix,icenter)
character(len=*), intent(in) :: field_prefix
character(len=*), intent(in) :: mesh_prefix
sll_real64, dimension(:,:), intent(in) :: vec_values
sll_int32, intent(in) :: nnodes_x1
sll_int32, intent(in) :: nnodes_x2
sll_int32             :: ncells_x1
sll_int32             :: ncells_x2

#ifndef NOHDF5
integer(hid_t)   :: hdffile_id
#else
sll_int32 :: binfile_id
#endif

sll_int32 :: xmffile_id

integer :: error
logical :: flag
integer, intent(in), optional :: icenter   ! centering ('node' or 'cell')

SLL_ASSERT(size(vec_values,1) == nnodes_x1)
SLL_ASSERT(size(vec_values,2) == nnodes_x2)
ncells_x1 = nnodes_x1 - 1
ncells_x2 = nnodes_x2 - 1

inquire(file=trim(mesh_prefix)//".xmf", exist=flag) 

if (.not. flag) then
   call errout(6,"W","sll_diagnostics:write_vec1d","Mesh file does not exist")
end if

if (.not. present(icenter)) then
   call errout(6,"W","sll_diagnostics:write_vec1d","Default: cell centered value")
end if

call sll_xmf_file_create(trim(field_prefix)//".xmf",xmffile_id,error)

call sll_xmf_grid_geometry_2d(xmffile_id, trim(mesh_prefix), nnodes_x1, nnodes_x2)

if (present(icenter) .and. icenter == NODE_CENTERED_DF) then
 call sll_xmf_field_2d(xmffile_id,'NodeValues',trim(field_prefix),nnodes_x1,nnodes_x2,'Node')
else
 call sll_xmf_field_2d(xmffile_id,'CellValues',trim(field_prefix),ncells_x1,ncells_x2,'Cell')
end if

call sll_xmf_file_close(xmffile_id,error)


#ifndef NOHDF5

call sll_hdf5_file_create(trim(field_prefix)//".h5",hdffile_id,error)
SLL_ASSERT(error==0)
if (present(icenter) .and. icenter == NODE_CENTERED_DF) then
   call sll_hdf5_write_array_2d(hdffile_id,vec_values,"NodeValues",error)
   SLL_ASSERT(error==0)
else
   call sll_hdf5_write_array_2d(hdffile_id,vec_values(1:ncells_x1,1:ncells_x2),"CellValues",error)
   SLL_ASSERT(error==0)
end if
call sll_hdf5_file_close(hdffile_id,error)
SLL_ASSERT(error==0)

#else

if (present(icenter) .and. icenter == NODE_CENTERED_DF) then
   call sll_binary_file_create(trim(field_prefix)//"-NodeValues.bin",binfile_id,error); SLL_ASSERT(error==0)
   call sll_binary_write_array_2d(binfile_id,vec_values,error); SLL_ASSERT(error==0)
   call sll_binary_file_close(binfile_id,error); SLL_ASSERT(error==0)
else
   call sll_binary_file_create(trim(field_prefix)//"-CellValues.bin",binfile_id,error); SLL_ASSERT(error==0)
   call sll_binary_write_array_2d(binfile_id,vec_values(1:ncells_x1,1:ncells_x2),error); SLL_ASSERT(error==0)
   call sll_binary_file_close(binfile_id,error); SLL_ASSERT(error==0)
end if

#endif

end subroutine write_vec1d
  
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> Write XDMF file for vector 2d on 2d mesh with data in binary files
!!
!! Light data are in :
!!  - <c>field_prefix.xmf</c>.
!! Heavy data :
!!  - mesh nodes coordinates created by write_mesh call are in :
!!     - <c>mesh_prefix.h5</c>
!!     - <c>mesh_prefix-x1.bin</c> and <c>mesh_prefix-x2.bin</c>
!!  - scalar values are in :
!!     - <c>field_prefix.h5</c>
!!     - <c>field_prefix-x1.bin</c> and <c>field_prefix-x2.bin</c>
!!
!! \param[in] vec_values_x1 Array with vector values
!! \param[in] vec_values_x2 Array with vector values
!! \param[in] nnodes_x1 nodes number along direction 1
!! \param[in] nnodes_x2 nodes number along direction 2
!! \param[in] field_prefix filename prefix for the vector
!! \param[in] mesh_prefix filename prefix for the mesh
!! \param[in] icenter parameter to distinguish cells values or nodes values
!!
!! Example of f1(x,v) = xv , f2(x,v) = x+v on a cartesian grid (xmin,xmax) x (vmin,vmax):
!!\verbatim
!!program grid
!!#include "sll_memory.h"
!!#include "sll_working_precision.h"
!!use sll_diagnostics
!!sll_int64, dimension(:,:), allocatable :: x, v
!!sll_int32 :: nx, nv, error
!!
!!nx = 32; xmin = 0.; xmax = 1.
!!nv = 64; vmin = 0.; vmax = 1.
!!SLL_ALLOCATE(x(nx,nv),error)
!!SLL_ALLOCATE(v(nx,nv),error) 
!!SLL_ALLOCATE(f1(nx,nv),error) 
!!SLL_ALLOCATE(f2(nx,nv),error) 
!!do j = 1, nv                                             
!!   do i = 1, nx                                         
!!      x(i,j) = float(i-1) * (xmax-xmin) / (nx-1)       
!!      v(i,j) = float(j-1) * (vmax-vmin) / (nv-1)      
!!      f1(i,j) = x(i,j) * v(i,j)
!!      f2(i,j) = x(i,j) + v(i,j)
!!   end do                                            
!!end do                                              
!!call write_mesh(x,v,nx,nv,"mesh")
!!call write_vec2d(f1,f2,nx,nv,"f","mesh",0) !write cell centered values
!!end program
!!\endverbatim
!! A file <i>f.xmf</i> is created, readable by VisIt or Paraview. 
!!
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

sll_int32        :: error
logical          :: flag
sll_int32        :: xmffile_id

#ifndef NOHDF5
integer(hid_t)   :: hdffile_id
#else
sll_int32        :: binfile_id
#endif

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
    
call sll_xmf_file_create(trim(field_prefix)//".xmf",xmffile_id,error)

call sll_xmf_grid_geometry_2d(xmffile_id, trim(mesh_prefix), nnodes_x1, nnodes_x2)

if (present(icenter) .and. icenter == NODE_CENTERED_DF) then
   call sll_xmf_field_2d(xmffile_id,'NodeValues1',trim(field_prefix),nnodes_x1,nnodes_x2,'Node')
   call sll_xmf_field_2d(xmffile_id,'NodeValues2',trim(field_prefix),nnodes_x1,nnodes_x2,'Node')
else
   call sll_xmf_field_2d(xmffile_id,'CellValues1',trim(field_prefix),ncells_x1,ncells_x2,'Cell')
   call sll_xmf_field_2d(xmffile_id,'CellValues2',trim(field_prefix),ncells_x1,ncells_x2,'Cell')
end if
    
call sll_xmf_file_close(xmffile_id,error)

#ifdef NOHDF5

if (present(icenter) .and. icenter == NODE_CENTERED_DF) then

   call sll_binary_file_create(trim(field_prefix)//"-NodeValues1.bin",binfile_id,error); SLL_ASSERT(error==0)
   call sll_binary_write_array_2d(binfile_id,vec_values_x1,error); SLL_ASSERT(error==0)
   call sll_binary_file_close(binfile_id,error); SLL_ASSERT(error==0)

   call sll_binary_file_create(trim(field_prefix)//"-NodeValues2.bin",binfile_id,error); SLL_ASSERT(error==0)
   call sll_binary_write_array_2d(binfile_id,vec_values_x2,error); SLL_ASSERT(error==0)
   call sll_binary_file_close(binfile_id,error); SLL_ASSERT(error==0)

else

   call sll_binary_file_create(trim(field_prefix)//"-CellValues1.bin",binfile_id,error); SLL_ASSERT(error==0)
   call sll_binary_write_array_2d(binfile_id,vec_values_x1(1:ncells_x1,1:ncells_x2),error); SLL_ASSERT(error==0)
   call sll_binary_file_close(binfile_id,error); SLL_ASSERT(error==0)

   call sll_binary_file_create(trim(field_prefix)//"-CellValues2.bin",binfile_id,error); SLL_ASSERT(error==0)
   call sll_binary_write_array_2d(binfile_id,vec_values_x2(1:ncells_x1,1:ncells_x2),error); SLL_ASSERT(error==0)
   call sll_binary_file_close(binfile_id,error); SLL_ASSERT(error==0)

end if

#else

call sll_hdf5_file_create(trim(field_prefix)//".h5",hdffile_id,error)
SLL_ASSERT(error==0)

if (present(icenter) .and. icenter == NODE_CENTERED_DF) then

   call sll_hdf5_write_array_2d(hdffile_id,vec_values_x1,"/NodeValues1",error); SLL_ASSERT(error==0)
   call sll_hdf5_write_array_2d(hdffile_id,vec_values_x2,"/NodeValues2",error); SLL_ASSERT(error==0)

else

   call sll_hdf5_write_array_2d(hdffile_id,vec_values_x1(1:ncells_x1,1:ncells_x2),"/CellValues1",error); SLL_ASSERT(error==0)
   call sll_hdf5_write_array_2d(hdffile_id,vec_values_x2(1:ncells_x1,1:ncells_x2),"/CellValues2",error); SLL_ASSERT(error==0)

end if
call sll_hdf5_file_close(hdffile_id,error); SLL_ASSERT(error==0)

#endif

end subroutine write_vec2d

!------------------------------------------------------------------------------
!> Outputs an error message:
!>   - PRTFIL : unit number for print-out
!>   - SEVRTY : 'W' - Warning 'F' - Fatal
!>   - WHERE  : in which program or subroutine
!>   - ErrMsg : error message
subroutine errout( prtfil, sevrty, lwhere, ErrMsg )
integer, intent(in) ::  prtfil
character(len=1),intent(in) :: sevrty 
character(len=*),intent(in) :: lwhere , ErrMsg

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
