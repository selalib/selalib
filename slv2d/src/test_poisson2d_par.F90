program test_2d_poisson_par

#include "selalib.h"
use hdf5
use sll_xdmf_parallel
use geometry_module
use poisson2dpp_seq

implicit none

sll_int32   :: error
sll_real64  :: tcpu1
sll_real64  :: tcpu2

character(len=8), parameter :: xfile = "xdata.h5" ! File name
character(len=8), parameter :: yfile = "ydata.h5" ! File name
character(len=8), parameter :: zfile = "zdata.h5" ! File name
character(len=8), parameter :: xdset = "xdataset" ! Dataset name
character(len=8), parameter :: ydset = "ydataset" ! Dataset name
character(len=8), parameter :: zdset = "zdataset" ! Dataset name

type(geometry)    :: geom     

sll_int32       :: nx, ny      
sll_real64      :: x_min, x_max
sll_real64      :: y_min, y_max

sll_real64, dimension(:,:),pointer :: x
sll_real64, dimension(:,:),pointer :: y
sll_real64, dimension(:,:),pointer :: rho
sll_real64, dimension(:,:),pointer :: phi
sll_real64, dimension(:,:),pointer :: ex
sll_real64, dimension(:,:),pointer :: ey

sll_int32  :: my_num
sll_int32  :: num_threads
sll_int32  :: comm

character(len=4)  :: prefix = "mesh"
sll_int32         :: file_id
integer(HSSIZE_T) :: offset(2)
integer(HSIZE_T)  :: global_dims(2)

! Boot parallel environment
call sll_boot_collective()

num_threads  = sll_get_collective_size(sll_world_collective)
my_num       = sll_get_collective_rank(sll_world_collective)
comm         = sll_world_collective%comm

x_min = 0.0
x_max = 2.*sll_pi
y_min = 0.0
y_max = 2.*sll_pi
nx  = 1024
ny  = 1024

SLL_ALLOCATE(rho(nx,ny),error)
SLL_ALLOCATE(phi(nx,ny),error)
SLL_ALLOCATE(ex(nx,ny),error)
SLL_ALLOCATE(ey(nx,ny),error)

call new_geometry2(geom,x_min,y_min,x_max,y_max,nx,ny,error,"perxy")

write(*,*) 'physical space: nx, ny, x0, x1, y0, y1, dx, dy'
write(*,"(2(i5,1x),6(g13.3,1x))") geom%nx, geom%ny, geom%x0, &
                                  geom%x0+(geom%nx)*geom%dx, &
                                  geom%y0, geom%y0+(geom%ny)*geom%dy, &
                                  geom%dx, geom%dy   

call meshgrid(geom%xgrid, geom%ygrid, x, y)
rho = -2_f64 * sin(x) * sin(y)
global_dims(1) = nx
global_dims(2) = ny
offset = 0

call sll_xdmf_open(my_num,"fields.xmf",prefix,nx,ny,file_id,error)
call sll_xdmf_write_array(prefix,global_dims,offset,x,'x1',error)
call sll_xdmf_write_array(prefix,global_dims,offset,y,'x2',error)
call sll_xdmf_write_array(prefix,global_dims,offset,rho,"rho",error,file_id,"Node")

tcpu1 = MPI_WTIME()
call solver_with_fftpack()
tcpu2 = MPI_WTIME()
if (my_num == 0) then
   write(*,"(//10x,' Wall time = ', G15.3, ' sec' )") (tcpu2-tcpu1)*num_threads
   print *, 'PASSED'
end if

#ifdef _FFTW

tcpu1 = MPI_WTIME()
call solver_with_fftw3()
tcpu2 = MPI_WTIME()
if (my_num == 0) then
   write(*,"(//10x,' Wall time = ', G15.3, ' sec' )") (tcpu2-tcpu1)*num_threads
   print *, 'PASSED'
end if

#endif

call sll_xdmf_write_array(prefix,global_dims,offset,phi,"phi",error,file_id,"Node")
call sll_xdmf_write_array(prefix,global_dims,offset,ex ,"ex" ,error,file_id,"Node")
call sll_xdmf_write_array(prefix,global_dims,offset,ey ,"ey" ,error,file_id,"Node")
call sll_xdmf_close(file_id,error)


call sll_halt_collective()

contains

subroutine solver_with_fftpack()
use poisson2dpp_seq
type(poisson2dpp) :: poisson 
sll_int32 :: i

call new(poisson, rho, geom, error)

do i = 1, 2
   rho = -2_f64 * sin(x) * sin(y)
   phi = 0_f64
   call solve(poisson,ex,ey,rho)
   print*, " error ex : ", sum(abs(ex-cos(x)*sin(y)))/(nx*ny)
   print*, " error ey : ", sum(abs(ey-sin(x)*cos(y)))/(nx*ny)
end do


end subroutine solver_with_fftpack

#ifdef _FFTW
subroutine solver_with_fftw3()
use poisson2d_periodic
type(poisson2dpp) :: poisson 
sll_int32 :: i, error

call new(poisson,ex,ey,geom,error)
do i = 1, 50
   rho = -2_f64 * sin(x) * sin(y)
   phi = 0_f64
   call solve(poisson, ex,ey,rho)
end do

print*, " error ex : ", sum(abs(ex-cos(x)*sin(y)))/(nx*ny)
print*, " error ey : ", sum(abs(ey-sin(x)*cos(y)))/(nx*ny)

end subroutine solver_with_fftw3
#endif

subroutine meshgrid(vec_x, vec_y, mat_x, mat_y)

   sll_real64, dimension(:), intent(in) :: vec_x
   sll_real64, dimension(:), intent(in) :: vec_y
   sll_real64, dimension(:,:), intent(out), pointer :: mat_x
   sll_real64, dimension(:,:), intent(out), pointer :: mat_y
   sll_int32 :: i, j, error
   sll_int32 :: nx, ny

   nx = size(vec_x)
   ny = size(vec_y)
   
   SLL_ALLOCATE(x(nx,ny),error)
   SLL_ALLOCATE(y(nx,ny),error)
   
   do j=1,ny
   do i=1,nx
      mat_x(i,j) = vec_x(i)
      mat_y(i,j) = vec_y(j)
   end do
   end do

end subroutine meshgrid

end program test_2d_poisson_par
