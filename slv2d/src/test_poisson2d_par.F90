program test_io_parallel

use mpi
use sll_collective
use sll_hdf5_io_parallel
use sll_xml_io

use geometry_module
use poisson2dpp_seq
!use sll_xdmf_parallel

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

!call sll_xdmf_open("fields.xmf",prefix,nx,ny,file_id,error)
!call sll_xdmf_write_array(prefix,global_dims,offset,x,'x1',error)
!call sll_xdmf_write_array(prefix,global_dims,offset,y,'x2',error)
!call sll_xdmf_write_array(prefix,global_dims,offset,rho,"rho",error,file_id,"Node")

call new(poisson, rho, geom, error)

write(*,*) 'physical space: nx, ny, x0, x1, y0, y1, dx, dy'
write(*,"(2(i3,1x),6(g13.3,1x))") geom%nx, geom%ny, geom%x0, &
                                  geom%x0+(geom%nx)*geom%dx, &
                                  geom%y0, geom%y0+(geom%ny)*geom%dy, &
                                  geom%dx, geom%dy   

call solve(poisson,ex,ey,rho)

!call sll_xdmf_write_array(prefix,global_dims,offset,ex,"ex",error,file_id,"Node")

!call sll_xdmf_write_array(prefix,global_dims, offset,&
!                          ey,"ey",error,file_id,"Node")

!call sll_xdmf_close(file_id,error)

print*, " error ex : ", sum(abs(ex - cos(x)*sin(y)))/(nx*ny)
print*, " error ey : ", sum(abs(ey - sin(x)*cos(y)))/(nx*ny)

tcpu2 = MPI_WTIME()
write(*,"(//10x,' Wall time = ', G15.3, ' sec' )") (tcpu2-tcpu1)*num_threads


call plot_layout2d()

tcpu2 = MPI_WTIME()
if (myrank == 0) &
   write(*,"(//10x,' Temps CPU = ', G15.3, ' sec' )") (tcpu2-tcpu1)*colsz


if( myrank .eq. 0) print *, 'PASSED'

call sll_halt_collective()
  
contains

! Take a 2D array of dimensions ni*nj where ni, nj are the dimensions of
! the full array.
 subroutine plot_layout2d()

  integer , parameter       :: nx = 512
  integer , parameter       :: ny = 256
  integer                   :: mx, my    ! Local sizes
  integer                   :: npi, npj
  sll_int32                 :: gi, gj
  
  sll_int32, dimension(2)   :: global_indices
  type(layout_2D), pointer  :: layout
  
  real(8), dimension(:,:), allocatable :: xdata, ydata, zdata
  integer(HID_T) :: file_id
  integer(HSIZE_T), dimension(2) :: datadims = (/nx,ny/)
  integer(HSSIZE_T), dimension(2) :: offset 
  
  layout => new_layout_2D( sll_world_collective )        

  call two_power_rand_factorization(colsz, npi, npj)
  
  if( myrank .eq. 0 ) then
     print *, '2d layout configuration: ', npi, npj
  end if
  
  call initialize_layout_with_distributed_2D_array( &
       nx, ny, npi, npj, layout )
       
  call sll_collective_barrier(sll_world_collective)
  
  call compute_local_sizes_2d( layout, mx, my)        
  
  SLL_ALLOCATE(xdata(mx,my),error)
  SLL_ALLOCATE(ydata(mx,my),error)
  SLL_ALLOCATE(zdata(mx,my),error)
  
  do j = 1, my
     do i = 1, mx
        global_indices =  local_to_global_2D( layout, (/i, j/) )
        gi = global_indices(1)
        gj = global_indices(2)
        xdata(i,j) = float(gi-1)/(nx-1)
        ydata(i,j) = float(gj-1)/(ny-1)
        zdata(i,j) = myrank !* xdata(i,j) * ydata(i,j)
     end do
  end do
  
  offset(1) =  get_layout_2D_i_min( layout, myrank ) - 1
  offset(2) =  get_layout_2D_j_min( layout, myrank ) - 1
  
  call sll_hdf5_file_create(xfile, file_id, error)
  call sll_hdf5_write_array(file_id, datadims,offset,xdata,xdset,error)
  call sll_hdf5_file_close(file_id,error)
  
  call sll_hdf5_file_create(yfile, file_id, error)
  call sll_hdf5_write_array(file_id, datadims,offset,ydata,ydset,error)
  call sll_hdf5_file_close(file_id,error)
  
  call sll_hdf5_file_create(zfile, file_id, error)
  call sll_hdf5_write_array(file_id, datadims,offset,zdata,zdset,error)
  call sll_hdf5_file_close(file_id,error)

  if (myrank == 0) then
  
     call sll_xml_file_create("layout2d.xmf",file_id,error)
     call sll_xml_grid_geometry(file_id, xfile, nx, yfile, ny, xdset, ydset )
     call sll_xml_field(file_id,'values', "zdata.h5:/zdataset",nx,ny,'HDF','Node')
     call sll_xml_file_close(file_id,error)
     print *, 'Printing 2D layout: '
     call sll_view_lims_2D( layout )
     print *, '--------------------'

  end if
  
  call delete_layout_2D( layout )
  
 end subroutine plot_layout2d

  subroutine two_power_rand_factorization(n, n1, n2, n3)
    sll_int64, intent(in) :: n
    integer, intent(out) ::n1, n2
    integer, intent(out), optional :: n3
    integer   :: expo, expo1, expo2, expo3
    sll_real64                :: rand_real
    if (.not.is_power_of_two(colsz)) then   
       print*, 'The number of processors must be a power of 2'
       call sll_halt_collective()
       stop
    endif 
    expo = int(log(real(n))/log(2.))  
    call random_number(rand_real)
    expo1 = int(rand_real*expo)
    if (present(n3)) then
       call random_number(rand_real)
       expo2 = int(rand_real*(expo-expo1))
       expo3 = expo - (expo1+expo2)
       n3 = 2**expo3
    else
       expo2 = expo - expo1
    end if

    n1 = 2**expo1
    n2 = 2**expo2

  end subroutine two_power_rand_factorization

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
