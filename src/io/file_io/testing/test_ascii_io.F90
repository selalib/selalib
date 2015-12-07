!>@internal [example] 
program test_ascii_io
#include "sll_memory.h"
#include "sll_working_precision.h"
use sll_m_ascii_io
use sll_m_xml_io
use sll_m_constants

sll_int32 :: i, j             !< indices
sll_int32 :: error            !< error code
sll_int32 :: file_id          !< file unit number

sll_int32 :: nnodes_x1, nnodes_x2

sll_real64 :: angle, xt, vt, R
sll_real64, allocatable, dimension(:)   :: theta
sll_real64, allocatable, dimension(:)   :: ray
sll_real64, allocatable, dimension(:,:) :: x1
sll_real64, allocatable, dimension(:,:) :: x2
sll_real64, allocatable, dimension(:,:) :: df


nnodes_x1 = 32
nnodes_x2 = 64

SLL_ALLOCATE(theta(nnodes_x2),error)
SLL_ALLOCATE(ray(nnodes_x1),error)
SLL_ALLOCATE(x1(nnodes_x1,nnodes_x2),error)
SLL_ALLOCATE(x2(nnodes_x1,nnodes_x2),error)

do j = 1, nnodes_x2
   vt = real(j-1,f64)/real(nnodes_x2-1,f64)
   angle = vt * 2. * sll_pi
   theta(j) = angle
   do i = 1, nnodes_x1
      xt = real(i-1,f64) / real(nnodes_x1-1,f64)
      R =  1 + xt
      ray(i) = R
      x1(i,j) = R * cos(angle)
      x2(i,j) = R * sin(angle)
   end do
end do 

SLL_ALLOCATE(df(nnodes_x1,nnodes_x2),error)

df = cos(2.*x1)*exp(-x2*x2)
 
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

print*, 'PASSED'

end program test_ascii_io
!>@internal [example] 
