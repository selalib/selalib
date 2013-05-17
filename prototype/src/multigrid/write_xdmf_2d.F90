subroutine write_xdmf_2d(my_id,nproc,f,sx,ex,sy,ey,hx,hy,error)

use hdf5 
use sll_utilities
!use mgd2

implicit none

integer, intent(in) :: my_id, nproc
integer, intent(in) :: sx,ex,sy,ey
real(8)             :: f(sx-1:ex+1,sy-1:ey+1)
real                :: x(sx-1:ex+1,sy-1:ey+1)
real                :: y(sx-1:ex+1,sy-1:ey+1)
real(8)             :: hx, hy
integer             :: nx, ny
integer, parameter  :: xmf = 77
character(len=72)   :: field_label, mesh_label
character(len=72)   :: field_name, mesh_name
character(len=4)    :: my_proc, cproc

character(len=2), dimension(2) :: coordNames

integer          :: error, i, j, iproc
integer(hid_t)   :: file_id, dataset_id, dataspace_id
integer(hsize_t) :: data_dims(2)
logical, save    :: mesh_writed = .false.
integer, save    :: iplot = 0
character(len=4) :: cplot

call int2string(my_id,my_proc)
field_label = "xdmf2d"
mesh_label = "mesh2d"

!Write the data file.
coordNames(1) = "/X"
coordNames(2) = "/Y"
    
!Write separate coordinate arrays for the x and y coordinates.
nx = size(f,1)
ny = size(f,2)
data_dims = (/nx,ny/)

iplot = iplot+1
call int2string(iplot,cplot)
mesh_name  = trim(mesh_label)//my_proc//".h5"
field_name = trim(field_label)//my_proc//"-"//cplot//".h5"

!Open the file and write the XML description of the mesh..
open(xmf,file=trim(field_label)//cplot//"-"//my_proc//".xmf")
write(xmf,'(a)')"<?xml version=""1.0"" ?>"
write(xmf,'(a)')"<!DOCTYPE Xdmf SYSTEM ""Xdmf.dtd"" []>"
write(xmf,'(a)')"<Xdmf Version=""2.0"">"
write(xmf,'(a)')"<Domain>"
write(xmf,'(a)')"<Grid Name=""mesh2D"" GridType=""Uniform"">"
write(xmf,'(a,2i5,a)')"<Topology TopologyType=""2DSMesh"" NumberOfElements='",ny,nx,"'/>"
write(xmf,'(a)')"<Geometry GeometryType=""X_Y"">"
write(xmf,'(a,2i5,a)')"<DataItem Dimensions='",ny,nx,"' NumberType=""Float"" Precision=""4"" Format=""HDF"">"
write(xmf,'(a)')trim(mesh_name)//":"//coordnames(1)
write(xmf,'(a)')"</DataItem>"
write(xmf,'(a,2i5,a)')"<DataItem Dimensions='",ny,nx,"' NumberType=""Float"" Precision=""4"" Format=""HDF"">"
write(xmf,'(a)')trim(mesh_name)//":"//coordnames(2)
write(xmf,'(a)')"</DataItem>"
write(xmf,'(a)')"</Geometry>"
write(xmf,'(a)')"<Attribute Name=""U"" AttributeType=""Scalar"" Center=""Node"">"
write(xmf,'(a,2i5,a)')"<DataItem Dimensions='",ny,nx,"' NumberType=""Float"" Precision=""8"" Format=""HDF"">"
write(xmf,'(a)')trim(field_name)//":/F"
write(xmf,'(a)')"</DataItem>"
write(xmf,'(a)')"</Attribute>"
write(xmf,'(a)')"</Grid>"
write(xmf,'(a)')"</Domain>"
write(xmf,'(a)')"</Xdmf>"
close(xmf)

!Initialize FORTRAN interface.
call H5open_f (error)
if (.not. mesh_writed) then

   do j=sy-1,ey+1
      do i=sx-1,ex+1
         x(i,j)=(float(i))/sngl(hx)
         y(i,j)=(float(j))/sngl(hy)
      end do
   end do

   !Create a new file using default properties.
   call H5Fcreate_f(trim(mesh_name), H5F_ACC_TRUNC_F, file_id, error);
   
   call H5Screate_simple_f(2, data_dims, dataspace_id, error)
   call H5Dcreate_f(file_id, coordnames(1), H5T_NATIVE_REAL, &
                    dataspace_id, dataset_id, error)
   call H5Dwrite_f(dataset_id, H5T_NATIVE_REAL, x(sx-1,sy-1), data_dims, error)
   call H5Dclose_f(dataset_id,error)
   call H5Sclose_f(dataspace_id,error)

   call H5Screate_simple_f(2, data_dims, dataspace_id, error)
   call H5Dcreate_f(file_id, coordnames(2), H5T_NATIVE_REAL, &
                    dataspace_id, dataset_id, error)
   call H5Dwrite_f(dataset_id, H5T_NATIVE_REAL, y(sx-1,sy-1), data_dims, error)
   call H5Dclose_f(dataset_id,error)
   call H5Sclose_f(dataspace_id,error)

   call H5fclose_f(file_id, error)

   mesh_writed = .true.

end if

call H5Fcreate_f(trim(field_name),H5F_ACC_TRUNC_F,file_id,error)
call H5Screate_simple_f(2,data_dims,dataspace_id,error)
call H5Dcreate_f(file_id,"/F",H5T_NATIVE_DOUBLE,dataspace_id,dataset_id,error)
call H5Dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,f(sx-1,sy-1),data_dims,error)
call H5Dclose_f(dataset_id,error)
call H5Sclose_f(dataspace_id,error)
call H5fclose_f(file_id,error)
call H5close_f(error) 

if (my_id == 0) then
   open(xmf,file="domain"//cplot//".xmf")
   write(xmf,'(a)')"<?xml version=""1.0"" ?>"
   write(xmf,'(a)')"<!DOCTYPE Xdmf SYSTEM ""Xdmf.dtd"" []>"
   write(xmf,'(a)')"<Xdmf Version=""2.0"">"
   write(xmf,'(a)')"<Domain Name=""mesh2d"">"
   write(xmf,'(a)')"<Grid Name=""Domain"" GridType=""Collection"" CollectionType=""Spatial"">"
   do iproc = 0, nproc-1
      call int2string(iproc,cproc)
      write(xmf,'(a)')"<Grid Name=""SubDomain"" GridType=""Uniform"">"
      write(xmf,'(a,2i5,a)')"<Topology TopologyType=""2DSMesh"" NumberOfElements='", &
                            ny,nx,"'/>"
      write(xmf,'(a)')"<Geometry GeometryType=""X_Y"">"
      write(xmf,'(a,2i5,a)')"<DataItem Dimensions='",ny,nx, &
      "' NumberType=""Float"" Precision=""4"" Format=""HDF"">"
      write(xmf,'(a)')trim(mesh_label)//cproc//".h5:"//coordnames(1)
      write(xmf,'(a)')"</DataItem>"
      write(xmf,'(a,2i5,a)')"<DataItem Dimensions='",ny,nx, &
      "' NumberType=""Float"" Precision=""4"" Format=""HDF"">"
      write(xmf,'(a)')trim(mesh_label)//cproc//".h5:"//coordnames(2)
      write(xmf,'(a)')"</DataItem>"
      write(xmf,'(a)')"</Geometry>"
      write(xmf,'(a)')"<Attribute Name=""U"" AttributeType=""Scalar"" Center=""Node"">"
      write(xmf,'(a,2i5,a)')"<DataItem Dimensions='",ny,nx, &
      "' NumberType=""Float"" Precision=""8"" Format=""HDF"">"
      write(xmf,'(a)')trim(field_label)//cproc//"-"//cplot//".h5:/F"
      write(xmf,'(a)')"</DataItem>"
      write(xmf,'(a)')"</Attribute>"
      write(xmf,'(a)')"</Grid>"
   end do
   write(xmf,'(a)')"</Grid>"
   write(xmf,'(a)')"</Domain>"
   write(xmf,'(a)')"</Xdmf>"
   close(xmf)
end if

end subroutine write_xdmf_2d 
