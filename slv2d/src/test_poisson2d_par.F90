program test_io_parallel

use mpi
use sll_collective
use sll_hdf5_io_parallel
use sll_xml_io

use geometry_module
use poisson2dpp_seq
use sll_xdmf_parallel

#include "sll_remap.h"
#include "sll_memory.h"
#include "sll_working_precision.h"
#include "misc_utils.h"

implicit none

sll_int32   :: myrank
sll_int64   :: colsz 
sll_int32   :: comm, info
sll_int32   :: error
sll_int32   :: i, j, k
sll_real64  :: tcpu1
sll_real64  :: tcpu2

character(len=8), parameter :: xfile = "xdata.h5" ! File name
character(len=8), parameter :: yfile = "ydata.h5" ! File name
character(len=8), parameter :: zfile = "zdata.h5" ! File name
character(len=8), parameter :: xdset = "xdataset" ! Dataset name
character(len=8), parameter :: ydset = "ydataset" ! Dataset name
character(len=8), parameter :: zdset = "zdataset" ! Dataset name


type(geometry)    :: geom     
type(poisson2dpp) :: poisson 

sll_int32       :: nx, ny      
sll_real64      :: x_min, x_max
sll_real64      :: y_min, y_max

sll_real64, dimension(:,:),pointer :: x
sll_real64, dimension(:,:),pointer :: y
sll_real64, dimension(:,:),pointer :: rho
sll_real64, dimension(:,:),pointer :: ex
sll_real64, dimension(:,:),pointer :: ey

sll_int32  :: my_num
sll_int32  :: num_threads

character(len=4)  :: prefix = "mesh"
sll_int32         :: file_id
integer(HSSIZE_T) :: offset(2)
integer(HSIZE_T)  :: global_dims(2)

! Boot parallel environment
call sll_boot_collective()

colsz  = sll_get_collective_size(sll_world_collective)
myrank = sll_get_collective_rank(sll_world_collective)
comm   = sll_world_collective%comm
info   = MPI_INFO_NULL

tcpu1 = MPI_WTIME()

x_min = 0.0
x_max = 2.*sll_pi
y_min = 0.0
y_max = 2.*sll_pi
nx  = 64
ny  = 64

SLL_ALLOCATE(rho(nx,ny),error)
SLL_ALLOCATE(ex(nx,ny),error)
SLL_ALLOCATE(ey(nx,ny),error)

call new_geometry2(geom,x_min,y_min,x_max,y_max,nx,ny,error,"perxy")

call meshgrid(geom%xgrid, geom%ygrid, x, y)
rho = -2_f64 * sin(x) * sin(y)
global_dims(1) = nx
global_dims(2) = ny
offset = 0

call sll_xdmf_open("fields.xmf",prefix,nx,ny,file_id,error)
call sll_xdmf_write_array(prefix,global_dims,offset,x,'x1',error)
call sll_xdmf_write_array(prefix,global_dims,offset,y,'x2',error)
call sll_xdmf_write_array(prefix,global_dims,offset,rho,"rho",error,file_id,"Node")

call new(poisson, rho, geom, error)

write(*,*) 'physical space: nx, ny, x0, x1, y0, y1, dx, dy'
write(*,"(2(i3,1x),6(g13.3,1x))") geom%nx, geom%ny, geom%x0, &
                                  geom%x0+(geom%nx)*geom%dx, &
                                  geom%y0, geom%y0+(geom%ny)*geom%dy, &
                                  geom%dx, geom%dy   

call solve(poisson,ex,ey,rho)

call sll_xdmf_write_array(prefix,global_dims,offset,ex,"ex",error,file_id,"Node")

call sll_xdmf_write_array(prefix,global_dims, offset,&
                          ey,"ey",error,file_id,"Node")

call sll_xdmf_close(file_id,error)

print*, " error ex : ", sum(abs(ex - cos(x)*sin(y)))/(nx*ny)
print*, " error ey : ", sum(abs(ey - sin(x)*cos(y)))/(nx*ny)

tcpu2 = MPI_WTIME()
write(*,"(//10x,' Wall time = ', G15.3, ' sec' )") (tcpu2-tcpu1)*num_threads


tcpu2 = MPI_WTIME()
if (myrank == 0) &
   write(*,"(//10x,' Temps CPU = ', G15.3, ' sec' )") (tcpu2-tcpu1)*colsz


if( myrank .eq. 0) print *, 'PASSED'

call sll_halt_collective()
  
contains

subroutine meshgrid(eta1, eta2, x, y)

   sll_real64, dimension(:), intent(in) :: eta1
   sll_real64, dimension(:), intent(in) :: eta2
   sll_real64, dimension(:,:), intent(out), pointer :: x
   sll_real64, dimension(:,:), intent(out), pointer :: y
   sll_int32 :: i, j, error
   
   SLL_ALLOCATE(x(size(eta1),size(eta2)),error)
   SLL_ALLOCATE(y(size(eta1),size(eta2)),error)
   
   do j = 1, size(eta2)
      do i = 1, size(eta2)
         x(i,j) = eta1(i)
         y(i,j) = eta2(j)
      end do
   end do

end subroutine meshgrid


end program test_io_parallel
