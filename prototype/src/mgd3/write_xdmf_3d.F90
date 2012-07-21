subroutine write_xdmf_3d(my_id,f,sx,ex,sy,ey,sz,ez,hx,hy,hz,error)

use hdf5 
use sll_misc_utils

implicit none

integer, intent(in) :: my_id
integer, intent(in) :: sx,ex,sy,ey,sz,ez
real(8) :: f(sx-1:ex+1,sy-1:ey+1,sz-1:ez+1)
real    :: x(sx-1:ex+1,sy-1:ey+1,sz-1:ez+1)
real    :: y(sx-1:ex+1,sy-1:ey+1,sz-1:ez+1)
real    :: z(sx-1:ex+1,sy-1:ey+1,sz-1:ez+1)
real(8) :: hx, hy, hz
integer :: nx, ny, nz
character(len=2), dimension(3) :: coordNames
integer, parameter :: xmf = 77
character(len=72) :: field_name, mesh_name
character(len=4)  :: my_proc

integer          :: error, i, j, k
integer(hid_t)   :: file_id, dataset_id, dataspace_id
integer(hsize_t) :: data_dims(3)
logical, save    :: mesh_writed = .false.

call int2string(my_id,my_proc)
field_name = "xdmf3d"
mesh_name = "mesh3d"

!Write the data file.
coordNames(1) = "/X"
coordNames(2) = "/Y"
coordNames(3) = "/Z"
    
!Write separate coordinate arrays for the x and y coordinates.
nx = size(f,1)
ny = size(f,2)
nz = size(f,3)
data_dims = (/nx,ny,nz/)

!Initialize FORTRAN interface.
call H5open_f (error)

if (.not. mesh_writed) then

   print*,'Writing mesh data...'
   do k=sz-1,ez+1
      do j=sy-1,ey+1
         do i=sx-1,ex+1
            x(i,j,k)=(float(i))/sngl(hx)
            y(i,j,k)=(float(j))/sngl(hy)
            z(i,j,k)=(float(k))/sngl(hz)
         end do
      end do
   end do

   !Create a new file using default properties.
   call H5Fcreate_f(trim(mesh_name)//my_proc//".h5", H5F_ACC_TRUNC_F, file_id, error);
   
   
   call H5Screate_simple_f(3, data_dims, dataspace_id, error)
   call H5Dcreate_f(file_id, coordnames(1), H5T_NATIVE_REAL, &
                    dataspace_id, dataset_id, error)
   call H5Dwrite_f(dataset_id, H5T_NATIVE_REAL, x(sx-1,sy-1,sz-1), data_dims, error)
   call H5Dclose_f(dataset_id,error);
   call H5Sclose_f(dataspace_id,error);

   call H5Screate_simple_f(3, data_dims, dataspace_id, error)
   call H5Dcreate_f(file_id, coordnames(2), H5T_NATIVE_REAL, &
                    dataspace_id, dataset_id, error)
   call H5Dwrite_f(dataset_id, H5T_NATIVE_REAL, y(sx-1,sy-1,sz-1), data_dims, error)
   call H5Dclose_f(dataset_id,error);
   call H5Sclose_f(dataspace_id,error);

   call H5Screate_simple_f(3, data_dims, dataspace_id, error)
   call H5Dcreate_f(file_id, coordnames(3), H5T_NATIVE_REAL, &
                    dataspace_id, dataset_id, error)
   call H5Dwrite_f(dataset_id, H5T_NATIVE_REAL, z(sx-1,sy-1,sz-1), data_dims, error)
   call H5Dclose_f(dataset_id,error);
   call H5Sclose_f(dataspace_id,error);

   !Terminate access to the file.
   call H5fclose_f(file_id, error)

   mesh_writed = .true.

end if

!Create a new file using default properties.
call H5Fcreate_f(trim(field_name)//my_proc//".h5", H5F_ACC_TRUNC_F, file_id, error);
 
!Write the scalar data.
call H5Screate_simple_f(3, data_dims, dataspace_id, error);
call H5Dcreate_f(file_id, "/F", H5T_NATIVE_DOUBLE, &
                 dataspace_id, dataset_id, error);
call H5Dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, f(sx-1,sy-1,sz-1), data_dims, error);
call H5Dclose_f(dataset_id,error);
call H5Sclose_f(dataspace_id,error);
!Terminate access to the file.
call H5fclose_f(file_id, error)

!Close FORTRAN interface.
call H5close_f(error) 



!<Grid Name="Car Wheel" GridType="Tree">
!<Grid Name="Tire" GridType="Uniform">
!<Topology ...
!<Geometry ...
!</Grid>
!<Grid Name="Lug Nuts" GridType="Collection">
!<Grid Name="Lug Nut 0" GridType="Uniform"
!<Topology ....
!<Geometry ...
!</Grid>
!<Grid Name="Lug Nut 1" GridType="Uniform"
!<Topology ...
!<Geometry ...
!</Grid>
!<Grid Name="Lug Nut 2" GridType="Uniform"
!<Topology ...
!<Geometry ...
!</Grid>
!</Grid>

!Open the file and write the XML description of the mesh..
open(xmf,file=trim(field_name)//my_proc//".xmf")
write(xmf,'(a)')"<?xml version=""1.0"" ?>"
write(xmf,'(a)')"<!DOCTYPE Xdmf SYSTEM ""Xdmf.dtd"" []>"
write(xmf,'(a)')"<Xdmf Version=""2.0"">"
write(xmf,'(a)')"<Domain>"
write(xmf,'(a)')"<Grid Name=""mesh3D"" GridType=""Uniform"">"
write(xmf,'(a,3i5,a)')"<Topology TopologyType=""3DSMesh"" NumberOfElements='",nz,ny,nx,"'/>"
write(xmf,'(a)')"<Geometry GeometryType=""X_Y_Z"">"
write(xmf,'(a,3i5,a)')"<DataItem Dimensions='",nz,ny,nx,"' NumberType=""Float"" Precision=""4"" Format=""HDF"">"
write(xmf,'(a)')trim(mesh_name)//my_proc//".h5:"//coordnames(1)
write(xmf,'(a)')"</DataItem>"
write(xmf,'(a,3i5,a)')"<DataItem Dimensions='",nz,ny,nx,"' NumberType=""Float"" Precision=""4"" Format=""HDF"">"
write(xmf,'(a)')trim(mesh_name)//my_proc//".h5:"//coordnames(2)
write(xmf,'(a)')"</DataItem>"
write(xmf,'(a,3i5,a)')"<DataItem Dimensions='",nz,ny,nx,"' NumberType=""Float"" Precision=""4"" Format=""HDF"">"
write(xmf,'(a)')trim(mesh_name)//my_proc//".h5:"//coordnames(3)
write(xmf,'(a)')"</DataItem>"
write(xmf,'(a)')"</Geometry>"
write(xmf,'(a)')"<Attribute Name=""U"" AttributeType=""Scalar"" Center=""Node"">"
write(xmf,'(a,3i5,a)')"<DataItem Dimensions='",nz,ny,nx,"' NumberType=""Float"" Precision=""8"" Format=""HDF"">"
write(xmf,'(a)')trim(field_name)//my_proc//".h5:/F"
write(xmf,'(a)')"</DataItem>"
write(xmf,'(a)')"</Attribute>"
write(xmf,'(a)')"</Grid>"
write(xmf,'(a)')"</Domain>"
write(xmf,'(a)')"</Xdmf>"
close(xmf)

end subroutine write_xdmf_3d 
