program diagnostics_tester

#include "sll_memory.h"
#include "sll_working_precision.h"
use sll_binary_io
use sll_xmf_io
use sll_diagnostics
use numeric_constants

#ifndef NOHDF5
use hdf5
use sll_hdf5_io
#endif

implicit none
sll_int32  :: i, j, k
sll_int32 :: error

call test_2d()

call test_3d()

contains

subroutine test_2d()

sll_int32  :: nx, nv
sll_real64 :: angle, xt, vt, R

sll_real64, allocatable, dimension(:,:) :: x
sll_real64, allocatable, dimension(:,:) :: v
sll_real64, allocatable, dimension(:,:) :: f

nx = 32
nv = 64

! Create the coordinate data.
SLL_ALLOCATE(x(nx,nv),error)
SLL_ALLOCATE(v(nx,nv),error)

do j = 1, nv
   vt = real(j-1)/(nv-1)
   angle = vt * 2. * sll_pi
   do i = 1, nx
      xt = real(i-1) / float(nx-1)
      R = (1.-xt)*2. + xt*5.
      x(i,j) = R * cos(angle)
      v(i,j) = R * sin(angle)
   end do
end do 

! Create the scalar data.
SLL_ALLOCATE(f(nx,nv),error)
do j = 1, nv
   do i = 1, nx
      f(i,j) = i*sin(float(j-1))
   end do
end do
 
call write_mesh(x,v,nx,nv,"mesh")

!cells values
call write_vec1d(f,nx,nv,"vec1d_on_cells","mesh",1)
!nodes values
call write_vec1d(f,nx,nv,"vec1d_on_nodes","mesh",0)

!cells values
call write_vec2d(x,v,nx,nv,"vec2d_on_cells","mesh",1)
!nodes values
call write_vec2d(x,v,nx,nv,"vec2d_on_nodes","mesh",0)

end subroutine test_2d

#ifdef NOHDF5

subroutine test_binary(x, v, f, nx, nv)
sll_real64, dimension(:,:) :: x, v, f
sll_int32 :: nx, nv
sll_int32 :: file_id
sll_int32 :: error
!Binary version
call sll_binary_file_create("test-f.bin",file_id,error)
call sll_binary_write_array_2d(file_id,f,error)
call sll_binary_file_close(file_id,error)

call sll_binary_file_create("test-x.bin",file_id,error)
call sll_binary_write_array_2d(file_id,x,error)
call sll_binary_file_close(file_id,error)

call sll_binary_file_create("test-v.bin",file_id,error)
call sll_binary_write_array_2d(file_id,v,error)
call sll_binary_file_close(file_id,error)

call sll_xmf_file_create("test.xmf",file_id,error)
write(file_id,"(a)"      ) "<Grid Name='mesh' GridType='Uniform'>"
write(file_id,"(a,2i5,a)") "<Topology TopologyType='2DSMesh' NumberOfElements='",nv,nx,"'/>"
write(file_id,"(a)"      ) "<Geometry GeometryType='X_Y'>"
write(file_id,"(a,2i5,a)") "<DataItem Dimensions='",nv,nx,"' NumberType='Float' Precision='8' Format='Binary'>"
write(file_id,"(a)"      ) "test-x.bin"
write(file_id,"(a)"      ) "</DataItem>"
write(file_id,"(a,2i5,a)") "<DataItem Dimensions='",nv,nx,"' NumberType='Float' Precision='8' Format='Binary'>"
write(file_id,"(a)"      ) "test-v.bin"
write(file_id,"(a)"      ) "</DataItem>"
write(file_id,"(a)"      ) "</Geometry>"
write(file_id,"(a)"      ) "<Attribute Name='NodesValues' AttributeType='Scalar' Center='Node'>"
write(file_id,"(a,2i5,a)") "<DataItem Dimensions='",nv,nx,"' NumberType='Float' Precision='8' Format='Binary'>"
write(file_id,"(a)"      ) "test-f.bin"
write(file_id,"(a)"      ) "</DataItem>"
write(file_id,"(a)"      ) "</Attribute>"
call sll_xmf_file_close(file_id,error)

end subroutine test_binary

#else

subroutine test_hdf5(x, v, nx, nv)
sll_real64, dimension(:,:) :: x, v
sll_int32 :: nx, nv

integer(hid_t)   :: file_id
integer(hsize_t) :: data_dims(2)
character(len=2) :: coord_names(2)

data_dims(1) = nx
data_dims(2) = nv

coord_names(1) = "/x"
coord_names(2) = "/v"

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
call sll_xmf_file_close(file_id,error)

end subroutine test_hdf5

#endif

!ASCII version
subroutine test_ascii( x, v, f, nx, nv )
sll_real64, dimension(:,:) :: x, v, f
sll_int32 :: nx, nv
sll_int32 :: file_id

call sll_xmf_file_create("test_ascii.xmf",file_id,error)
write(file_id,"(a)")"<Grid Name='mesh' GridType='Uniform'>"
write(file_id,"(a,2i5,a)")"<Topology TopologyType='2DSMesh' NumberOfElements='",nx,nv,"'/>"
write(file_id,"(a)")"<Geometry GeometryType='X_Y'>"
write(file_id,"(a,2i5,a)")"<DataItem Dimensions='",nx,nv,"' NumberType='Float' Precision='8' Format='XML'>"
do i = 1, nx
write(file_id,*)(x(i,j),j=1,nv)
end do
write(file_id,"(a)")"</DataItem>"
write(file_id,"(a,2i5,a)")"<DataItem Dimensions='",nx,nv,"' NumberType='Float' Precision='8' Format='XML'>"
do i = 1, nx
write(file_id,*)(v(i,j),j=1,nv)
end do
write(file_id,"(a)")"</DataItem>"
write(file_id,"(a)")"</Geometry>"
write(file_id,"(a)")"<Attribute Name='NodesValues' AttributeType='Scalar' Center='Node'>"
write(file_id,"(a,2i5,a)")"<DataItem Dimensions='",nx,nv,"' NumberType='Float' Precision='8' Format='XML'>"
do i = 1, nx
write(file_id,*)(f(i,j),j=1,nv)
end do
write(file_id,"(a)")"</DataItem>"
write(file_id,"(a)")"</Attribute>"
call sll_xmf_file_close(file_id,error)

end subroutine test_ascii

subroutine test_3d()
sll_int32  :: n1, n2, n3
sll_real64 :: theta, a, b, phi
sll_real64, allocatable, dimension(:,:,:) :: p1, p2, p3, pf

n1 = 32
n2 = 32
n3 = 32
SLL_ALLOCATE(p1(n1,n2,n3),error)
SLL_ALLOCATE(p2(n1,n2,n3),error)
SLL_ALLOCATE(p3(n1,n2,n3),error)
SLL_ALLOCATE(pf(n1,n2,n3),error)

a = 3
phi = 0
do k = 1, n3
   theta = 0
   do j = 1, n2
      b = 0
      do i = 1, n1
         p1(i,j,k) =  (a + b*cos(phi))*cos(theta)
         p2(i,j,k) =  (a + b*cos(phi))*sin(theta)
         p3(i,j,k) =  b*sin(phi)
         pf(i,j,k) =  sin(phi)*cos(theta)
         b = b + 1._f64/(n1-1)
      end do
      theta = theta + 2._f64*sll_pi / (n2-1)
   end do 
   phi = phi + 2._f64*sll_pi / (n3-1)
end do 

call write_mesh(p1,p2,p3,n1,n2,n3,"tore3d")
call write_vec1d(pf,n1,n2,n3,"vec3d_on_cells","tore3d",1)
call write_vec1d(pf,n1,n2,n3,"vec3d_on_nodes","tore3d",0)

do k = 1, n3
   do j = 1, n2
      do i = 1, n1
         p1(i,j,k) =  i
         p2(i,j,k) =  j
         p3(i,j,k) =  k
         pf(i,j,k) =  (i+j)*k
      end do
   end do 
end do 

call write_mesh(p1,p2,p3,n1,n2,n3,"cube3d")
call write_vec1d(pf,n1,n2,n3,"vec3d_on_cells","cube3d",1)
call write_vec1d(pf,n1,n2,n3,"vec3d_on_nodes","cube3d",0)

end subroutine test_3d

end program diagnostics_tester
