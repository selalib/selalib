module test_io
#include "sll_memory.h"
#include "sll_working_precision.h"
use sll_xdmf
#ifndef NOHDF5
use sll_hdf5_io
#endif
use sll_constants
use sll_utilities, only: int2string

sll_int32, private :: i, j, k !< indices
sll_int32 :: error            !< error code
sll_int32 :: file_id          !< file unit number
sll_int32 :: iplot            !< plot counter
character(len=4) :: cplot     !< plot counter

#ifndef NOHDF5
integer(hid_t) :: hfile_id    !< file unit number
#endif

contains

!>Unit test program for xdmf outputs in 2d
subroutine test_io_2d()

sll_int32 :: nnodes_x1, nnodes_x2
sll_int32 :: ncells_x1, ncells_x2

sll_real64 :: angle, xt, vt, R
sll_real64, allocatable, dimension(:)   :: theta
sll_real64, allocatable, dimension(:)   :: ray
sll_real64, allocatable, dimension(:,:) :: x1
sll_real64, allocatable, dimension(:,:) :: x2
sll_real64, allocatable, dimension(:,:) :: df

character(6)  :: mesh_name = "grid2d"
character(10) :: file_name = "test2d.xmf"

nnodes_x1 = 32
nnodes_x2 = 64
ncells_x1 = nnodes_x1 - 1
ncells_x2 = nnodes_x2 - 1

SLL_ALLOCATE(theta(nnodes_x2),error)
SLL_ALLOCATE(ray(nnodes_x1),error)
SLL_ALLOCATE(x1(nnodes_x1,nnodes_x2),error)
SLL_ALLOCATE(x2(nnodes_x1,nnodes_x2),error)

do j = 1, nnodes_x2
   vt = real(j-1)/(nnodes_x2-1)
   angle = vt * 2. * sll_pi
   theta(j) = angle
   do i = 1, nnodes_x1
      xt = real(i-1) / float(nnodes_x1-1)
      R =  1 + xt
      ray(i) = R
      x1(i,j) = R * cos(angle)
      x2(i,j) = R * sin(angle)
   end do
end do 

SLL_ALLOCATE(df(nnodes_x1,nnodes_x2),error)

df = cos(2.*x1)*exp(-x2*x2)
 
call sll_xdmf_open(file_name,mesh_name,nnodes_x1,nnodes_x2,file_id,error)
call sll_xdmf_write_array(mesh_name,x1,'x1',error)
call sll_xdmf_write_array(mesh_name,x2,'x2',error)
call sll_xdmf_write_array("test2d",df,"NodeVal",error,file_id,"Node")
call sll_xdmf_write_array("test2d",df(1:ncells_x1,1:ncells_x2),"CellVal",error,file_id,"Cell")
call sll_xdmf_close(file_id,error)



!ASCII version just in case of problem with binary format
call sll_xml_file_create("test_ascii.xmf",file_id,error)
write(file_id,"(a)")"<Grid Name='mesh' GridType='Uniform'>"
write(file_id,"(a,2i5,a)")"<Topology TopologyType='2DSMesh' NumberOfElements='", &
                          nnodes_x2,nnodes_x1,"'/>"
write(file_id,"(a)")"<Geometry GeometryType='X_Y'>"
write(file_id,"(a,2i5,a)")"<DataItem Dimensions='",nnodes_x2,nnodes_x1, &
                          "' NumberType='Float' Precision='8' Format='XML'>"
call sll_ascii_write_array(file_id,x1,error)
write(file_id,"(a)")"</DataItem>"
write(file_id,"(a,2i5,a)")"<DataItem Dimensions='",nnodes_x2,nnodes_x1, &
                          "' NumberType='Float' Precision='8' Format='XML'>"
call sll_ascii_write_array(file_id,x2,error)
write(file_id,"(a)")"</DataItem>"
write(file_id,"(a)")"</Geometry>"
write(file_id,"(a)")"<Attribute Name='NodesVal' AttributeType='Scalar' Center='Node'>"
write(file_id,"(a,2i5,a)")"<DataItem Dimensions='",nnodes_x2,nnodes_x1, &
                          "' NumberType='Float' Precision='8' Format='XML'>"
call sll_ascii_write_array(file_id,df,error)
write(file_id,"(a)")"</DataItem>"
write(file_id,"(a)")"</Attribute>"
call sll_xml_file_close(file_id,error)

!The field is describe on a cartesian mesh
!Axis are perpendicular and spacing is constant
call sll_xdmf_corect2d_nodes( "test_corect2d", df, "f1_2d", &
                              ray(1), (ray(nnodes_x1)-ray(1))/(nnodes_x1-1), &
                              theta(1), (theta(nnodes_x2)-theta(1))/(nnodes_x2-1), &
                              "HDF5") 
!The field is describe on a cartesian mesh
!Axis are perpendicular and spacing is define by eta1 and eta2 arrays
call sll_xdmf_rect2d_nodes( "test_rect2d", df, "f2_2d", ray, theta, "HDF5") 

!The field is describe on a structured curvilinear mesh.
!Nodes coordinates are defined by eta1 and eta2 that are 2d arrays.
call sll_xdmf_curv2d_nodes( "test_curv2d", df, "f3_2d", x1, x2, "HDF5") 


#ifndef NOHDF5
!Init step, create h5 files with mesh coordinates
call sll_hdf5_file_create("polar_mesh-x1.h5",hfile_id,error)
call sll_hdf5_write_array(hfile_id,x1,"/x1",error)
call sll_hdf5_file_close(hfile_id, error)
call sll_hdf5_file_create("polar_mesh-x2.h5",hfile_id,error)
call sll_hdf5_write_array(hfile_id,x2,"/x2",error)
call sll_hdf5_file_close(hfile_id, error)
#endif

!plot 10 fields using mesh coordinates written before
do iplot = 1, 10
   call int2string(iplot,cplot)
   call sll_xdmf_open("f"//cplot//".xmf","polar_mesh", &
                      nnodes_x1,nnodes_x2,file_id,error)
   df = cos(2.*x1)*exp(-x2*x2)*sin(2*sll_pi*iplot/10.)
   call sll_xdmf_write_array("f"//cplot,df,"f_values",error,file_id,"Node")
   call sll_xdmf_close(file_id,error)
end do


end subroutine test_io_2d

!> Unit test program for xdmf outputs in 3d
subroutine test_io_3d()
implicit none
sll_int32  :: nnodes_x1, nnodes_x2, nnodes_x3
sll_int32  :: ncells_x1, ncells_x2, ncells_x3
sll_real64 :: a, b, phi , theta
sll_real64, allocatable, dimension(:) :: vx, vy, vz
sll_real64, allocatable, dimension(:,:,:) :: x1, x2, x3, df

character(6)  :: mesh_name = "mesh3d"
character(10) :: file_name = "test3d.xmf"

nnodes_x1 = 32
nnodes_x2 = 64
nnodes_x3 = 128

ncells_x1 = nnodes_x1-1
ncells_x2 = nnodes_x2-1
ncells_x3 = nnodes_x3-1

SLL_ALLOCATE(x1(nnodes_x1,nnodes_x2,nnodes_x3),error)
SLL_ALLOCATE(x2(nnodes_x1,nnodes_x2,nnodes_x3),error)
SLL_ALLOCATE(x3(nnodes_x1,nnodes_x2,nnodes_x3),error)
SLL_ALLOCATE(df(nnodes_x1,nnodes_x2,nnodes_x3),error)

SLL_ALLOCATE(vx(nnodes_x1),error)
SLL_ALLOCATE(vy(nnodes_x2),error)
SLL_ALLOCATE(vz(nnodes_x3),error)

a = 3
phi = 0
do k = 1, nnodes_x3
   theta = 0
   do j = 1, nnodes_x2
      b = 0
      do i = 1, nnodes_x1
         x1(i,j,k) =  (a + b*cos(phi))*cos(theta)
         x2(i,j,k) =  (a + b*cos(phi))*sin(theta)
         x3(i,j,k) =  b*sin(phi)
         df(i,j,k) =  phi*theta
         vx(i) = b
         b = b + 1._f64/(nnodes_x1-1)
      end do
      vy(j) = theta
      theta = theta + 2._f64*sll_pi / (nnodes_x2-1)
   end do 
   vz(k) = phi
   phi = phi + 2._f64*sll_pi / (nnodes_x3-1)
end do 

call sll_xdmf_open(file_name,mesh_name,nnodes_x1,nnodes_x2,nnodes_x3,file_id,error)
call sll_xdmf_write_array(mesh_name,x1,'x1',error)
call sll_xdmf_write_array(mesh_name,x2,'x2',error)
call sll_xdmf_write_array(mesh_name,x3,'x3',error)
call sll_xdmf_write_array("field3d",df,"NodeVal",error,file_id,"Node")
call sll_xdmf_write_array("field3d",df(1:ncells_x1,1:ncells_x2,1:ncells_x3), &
                          "CellVal",error,file_id,"Cell")
call sll_xdmf_close(file_id,error)

!The field is describe on a cartesian mesh
!Axis are perpendicular and spacing is constant
call sll_xdmf_corect3d_nodes( "test_corect3d", df, "f1_3d", &
                              vx(1), (vx(nnodes_x1)-vx(1))/(nnodes_x1-1), &
                              vy(1), (vy(nnodes_x2)-vy(1))/(nnodes_x2-1), &
                              vz(1), (vz(nnodes_x2)-vz(1))/(nnodes_x3-1), &
                              "HDF5") 

!The field is describe on a cartesian mesh
!Axis are perpendicular and spacing is define by eta1 and eta2 arrays
call sll_xdmf_rect3d_nodes( "test_rect3d", df, "f2_3d", vx, vy, vz, "HDF5") 

!The field is describe on a structured curvilinear mesh.
!Nodes coordinates are defined by eta1 and eta2 that are 2d arrays.
call sll_xdmf_curv3d_nodes( "test_curv3d", df, "f3_3d", x1, x2, x3, "HDF5") 

end subroutine test_io_3d

end module test_io

!>Unit test program for xdmf outputs
program test_io_module

use test_io

implicit none


call test_io_2d()

call test_io_3d()

print*,"PASSED"

end program test_io_module
