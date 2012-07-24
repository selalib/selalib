program unit_test

#include "sll_memory.h"
#include "sll_working_precision.h"
use sll_xdmf
use numeric_constants

implicit none

sll_int32 :: i, j, k
sll_int32 :: error
sll_int32 :: file_id

call test_2d()

call test_3d()

print*,"PASSED"

contains

subroutine test_2d()

sll_int32 :: nnodes_x1, nnodes_x2
sll_int32 :: ncells_x1, ncells_x2

sll_real64 :: angle, xt, vt, R
sll_real64, allocatable, dimension(:,:) :: x1
sll_real64, allocatable, dimension(:,:) :: x2
sll_real64, allocatable, dimension(:,:) :: df

character(6)  :: mesh_name = "grid2d"
character(10) :: file_name = "test2d.xmf"

nnodes_x1 = 32
nnodes_x2 = 64
ncells_x1 = nnodes_x1 - 1
ncells_x2 = nnodes_x2 - 1

SLL_ALLOCATE(x1(nnodes_x1,nnodes_x2),error)
SLL_ALLOCATE(x2(nnodes_x1,nnodes_x2),error)

do j = 1, nnodes_x2
   vt = real(j-1)/(nnodes_x2-1)
   angle = vt * 2. * sll_pi
   do i = 1, nnodes_x1
      xt = real(i-1) / float(nnodes_x1-1)
      R =  1 + xt
      x1(i,j) = R * cos(angle)
      x2(i,j) = R * sin(angle)
   end do
end do 

! Create the scalar data.
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

end subroutine test_2d

subroutine test_3d()
sll_int32  :: nnodes_x1, nnodes_x2, nnodes_x3
sll_int32  :: ncells_x1, ncells_x2, ncells_x3
sll_real64 :: theta, a, b, phi
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
         b = b + 1._f64/(nnodes_x1-1)
      end do
      theta = theta + 2._f64*sll_pi / (nnodes_x2-1)
   end do 
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

end subroutine test_3d

end program unit_test
