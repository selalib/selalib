program diagnostics_tester
#ifndef NODHF5
use hdf5
#endif

#include "sll_working_precision.h"
use sll_binary_io
use sll_hdf5_io
use sll_xmf_io
use sll_diagnostics

implicit none

sll_int32  :: nx, nv, i, j
sll_real64 :: angle, xt, vt, pi, R

sll_real64, allocatable, dimension(:,:) :: x
sll_real64, allocatable, dimension(:,:) :: v
sll_real64, allocatable, dimension(:,:) :: f

sll_int32 :: error

#ifdef NOHDF5
sll_int32 :: file_id
#else
integer(hid_t)   :: file_id
integer(hsize_t) :: data_dims(2)
character(len=2) :: coord_names(2)
#endif
  
pi = 4.*atan(1.)

nx = 32
nv = 128

! Create the coordinate data.
allocate(x(nx,nv))
allocate(v(nx,nv))

do j = 1, nv
   vt = real(j-1)/(nv-1)
   angle = vt * 2. * pi
   do i = 1, nx
      xt = real(i-1) / float(nx-1)
      R = (1.-xt)*2. + xt*5.
      x(i,j) = R * cos(angle)
      v(i,j) = R * sin(angle)
   end do
end do 

! Create the scalar data.
allocate(f(nx,nv))
do j = 1, nv
   do i = 1, nx
      f(i,j) = i*sin(float(j-1))
      f(i,j) = i*sin(float(j-1))
   end do
end do
 
call write_mesh(x,v,nx,nv,"mesh")

!cells values
call write_vec1d(f,nx,nv,"vec1d_on_cells","mesh")
call write_vec2d(x,v,nx,nv,"vec2d_on_cells","mesh")

!nodes values
call write_vec1d(f,nx,nv,"vec1d_on_nodes","mesh", 0)
call write_vec2d(x,v,nx,nv,"vec2d_on_nodes","mesh", 0)

!Binary version

call sll_binary_file_create("test-f.bin",file_id,error)
call sll_binary_write_array_2d(file_id,f,error)
call sll_binary_file_close(file_id,error)

call sll_binary_file_create("test-x.bin",file_id,error)
call sll_binary_write_array_2d(file_id,f,error)
call sll_binary_file_close(file_id,error)

call sll_binary_file_create("test-v.bin",file_id,error)
call sll_binary_write_array_2d(file_id,f,error)
call sll_binary_file_close(file_id,error)

call sll_xmf_file_create("test.xmf",file_id,error)
write(file_id,"(a)")"<Grid Name='mesh' GridType='Uniform'>"
write(file_id,"(a,2i4,a)")"<Topology TopologyType='2DSMesh' NumberOfElements='",nv,nx,"'/>"
write(file_id,"(a)")"<Geometry GeometryType='X_Y'>"
write(file_id,"(a,2i4,a)")"<DataItem Dimensions='",nv,nx,&
                     "' NumberType='Float' Precision='8' Endian='Little' Format='Binary'>"
write(file_id,"(a)")"test-x.bin"
write(file_id,"(a)")"</DataItem>"
write(file_id,"(a,2i4,a)")"<DataItem Dimensions='",nv,nx,&
                     "' NumberType='Float' Precision='8' Endian='Little' Format='Binary'>"
write(file_id,"(a)")"test-v.bin"
write(file_id,"(a)")"</DataItem>"
write(file_id,"(a)")"</Geometry>"
write(file_id,"(a)")"<Attribute Name='NodesValues' AttributeType='Scalar' Center='Node'>"
write(file_id,"(a,2i4,a)")"<DataItem Dimensions='",nv,nx,&
                     "' NumberType='Float' Precision='8' Endian='Little' Format='Binary'>"
write(file_id,"(a)")"test-f.bin"
write(file_id,"(a)")"</DataItem>"
write(file_id,"(a)")"</Attribute>"
write(file_id,"(a)")"</Grid>"
call sll_xmf_file_close(file_id,error)

!HDF5 version

#ifndef NOHDF5

coord_names(1) = "/x"
coord_names(2) = "/v"
    
data_dims(1) = nx
data_dims(2) = nv

call sll_hdf5_file_create("test.h5",file_id,error)
call sll_hdf5_write_array_2d(file_id,x,coord_names(1),error)
call sll_hdf5_write_array_2d(file_id,v,coord_names(2),error)
call sll_hdf5_file_close(file_id, error)

!write the xmf file readable by VisIt
call sll_xmf_file_create("test_h5.xmf",file_id,error)
write(file_id,"(a)")"<Grid Name='mesh' GridType='Uniform'>"
write(file_id,"(a,2i5,a)")"<Topology TopologyType='2DSMesh' NumberOfElements='", &
                     data_dims(2),data_dims(1),"'/>"
write(file_id,"(a)")"<Geometry GeometryType='X_Y'>"
write(file_id,"(a,2i5,a)")"<DataItem Dimensions='",data_dims(2),data_dims(1) &
                    ,"' NumberType='Float' Precision='8' Format='HDF'>"
write(file_id,"(a)")"test.h5:"//coord_names(1)
write(file_id,"(a)")"</DataItem>"
write(file_id,"(a,2i5,a)")"<DataItem Dimensions='",data_dims(2),data_dims(1) &
                    ,"' NumberType='Float' Precision='8' Format='HDF'>"
write(file_id,"(a)")"test.h5:"//coord_names(2)
write(file_id,"(a)")"</DataItem>"
write(file_id,"(a)")"</Geometry>"
write(file_id,"(a)")"</Grid>"
call sll_xmf_file_close(file_id,error)

#endif

!ASCII version

call sll_xmf_file_create("test_ascii.xmf",file_id,error)
call sll_xmf_grid_geometry_2d(file_id, nx, nv)
write(file_id,"(a,2i5,a)")"<DataItem Dimensions='", nx, nv &
                    ,"' NumberType='Float' Precision='8' Format='XML'>"
do i = 1, nx
   write(file_id,*)(x(i,j),j=1,nv)
end do
write(file_id,"(a)")"</DataItem>"
write(file_id,"(a,2i5,a)")"<DataItem Dimensions='", nx, nv &
                    ,"' NumberType='Float' Precision='8' Format='XML'>"
do i = 1, nx
   write(file_id,*)(v(i,j),j=1,nv)
end do
write(file_id,"(a)")"</DataItem>"
call sll_xmf_file_close(file_id,error)

end program diagnostics_tester
