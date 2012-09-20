program test_io_parallel

use mpi
use sll_collective
use sll_hdf5_io_parallel
use sll_xml_io

#include "sll_remap.h"
#include "sll_memory.h"
#include "sll_working_precision.h"
#include "misc_utils.h"

implicit none

! Take a 2D array of dimensions ni*nj where ni, nj are the dimensions of
! the full array.
integer , parameter                       :: nx = 512
integer , parameter                       :: ny = 256
! Local sizes
integer                                   :: loc_sz_i_init
integer                                   :: loc_sz_j_init

! the process mesh
integer                                   :: npi
integer                                   :: npj
sll_int32                                 :: gi, gj
integer                                   :: ierr
integer                                   :: myrank
sll_int64                                 :: colsz        ! collective size
type(layout_2D), pointer                  :: layout

sll_real64                                :: rand_real
integer                                   :: i, j
sll_int32, dimension(2)                   :: global_indices
logical                                   :: test_passed

integer                  :: mx, my

!Parameters
real(8)                  :: tcpu1
real(8)                  :: tcpu2

real(8), dimension(:,:), allocatable :: xdata, ydata, zdata
character(len=8), parameter :: xfile = "xdata.h5" ! File name
character(len=8), parameter :: yfile = "ydata.h5" ! File name
character(len=8), parameter :: zfile = "zdata.h5" ! File name
character(len=8), parameter :: xdset = "xdataset" ! Dataset name
character(len=8), parameter :: ydset = "ydataset" ! Dataset name
character(len=8), parameter :: zdset = "zdataset" ! Dataset name
type(phdf5_file) :: my_file
integer(HSIZE_T), dimension(2) :: datadims = (/nx,ny/)
integer(HSSIZE_T), dimension(2) :: offset 

sll_int32 :: file_id, error

test_passed = .true.

! Boot parallel environment
call sll_boot_collective()

colsz  = sll_get_collective_size(sll_world_collective)
myrank = sll_get_collective_rank(sll_world_collective)

if( myrank .eq. 0) then
   print *, ' '
   print *, '--------------- HDF5 parallel test ---------------------'
   print *, ' '
   print"('Running a test on ',i4,' processes')", colsz
   call flush()
end if

if (.not. is_power_of_two(colsz)) then     
   print *, 'This test needs to run in a number of processes which is ',&
        'a power of 2.'
   call sll_halt_collective()
   stop
end if

layout => new_layout_2D( sll_world_collective )        
call factorize_in_random_2powers_2d(colsz, npi, npj)

if( myrank .eq. 0 ) then
   print *, 'source configuration: ', npi, npj
end if

call initialize_layout_with_distributed_2D_array( &
     nx, ny, npi, npj, layout )
     
call sll_collective_barrier(sll_world_collective)

tcpu1 = MPI_WTIME()

call compute_local_sizes_2d( layout, mx, my)        

SLL_ALLOCATE(xdata(mx,my),ierr)
SLL_ALLOCATE(ydata(mx,my),ierr)
SLL_ALLOCATE(zdata(mx,my),ierr)

do j = 1, my
   do i = 1, mx
      global_indices =  local_to_global_2D( layout, (/i, j/) )
      gi = global_indices(1)
      gj = global_indices(2)
      xdata(i,j) = float(gi-1)
      ydata(i,j) = float(gj-1)
      zdata(i,j) = myrank * xdata(i,j) * ydata(i,j)
   end do
end do

offset(1) =  get_layout_2D_i_min( layout, myrank ) - 1;          \
offset(2) =  get_layout_2D_j_min( layout, myrank ) - 1;          \

call my_file%create(xfile, error)
call my_file%write_array(datadims,offset,xdata,xdset,error)
call my_file%close(error)

call my_file%create(yfile, error)
call my_file%write_array(datadims,offset,ydata,ydset,error)
call my_file%close(error)

call my_file%create(zfile, error)
call my_file%write_array(datadims,offset,zdata,zdset,error)
call my_file%close(error)

call sll_xml_file_create("parallel.xmf",file_id,error)
call sll_xml_grid_geometry(file_id, xfile, xdset, nx, yfile, ydset, ny )
call sll_xml_field(file_id,'Z', "zdata.h5:/zdataset",nx,ny,'HDF','Node')
call sll_xml_file_close(file_id,error)

call my_file%create('layout.h5', error)
call my_file%write_array(datadims,offset,xdata,'x1',error)
call my_file%write_array(datadims,offset,ydata,'x2',error)
call my_file%close(error)
call sll_xml_file_create("layout.xmf",file_id,error)
call sll_xml_grid_geometry(file_id, 'layout', nx, ny )
call sll_xml_field(file_id,'Z', "zdata.h5:/zdataset",nx,ny,'HDF','Node')
call sll_xml_file_close(file_id,error)

tcpu2 = MPI_WTIME()
if (myrank == 0) &
   write(*,"(//10x,' Temps CPU = ', G15.3, ' sec' )") (tcpu2-tcpu1)*colsz

call delete_layout_2D( layout )
call sll_halt_collective()

if( myrank .eq. 0) print *, 'PASSED'
  
contains

  subroutine factorize_in_random_2powers_2d(n, n1, n2)
    sll_int64, intent(in) :: n
    integer, intent(out)  ::n1, n2
    integer   :: expo, expo1, expo2
    if (.not.is_power_of_two(n)) then   
       print*, 'The number of processors must be a power of 2'
       stop
    endif 
    expo = int(log(real(n))/log(2.))  
    call random_number(rand_real)
    expo1 = int(rand_real*expo)
    expo2 = expo - expo1
    n1 = 2**expo1
    n2 = 2**expo2
  end subroutine factorize_in_random_2powers_2d

  subroutine compute_local_sizes_2d( layout, loc_sz_i, loc_sz_j )
    type(layout_2D), pointer :: layout
    sll_int32, intent(out) :: loc_sz_i
    sll_int32, intent(out) :: loc_sz_j
    sll_int32 :: i_min
    sll_int32 :: i_max
    sll_int32 :: j_min
    sll_int32 :: j_max
    sll_int32 :: my_rank
    if( .not. associated(layout) ) then
       print *, 'not-associated layout passed to new_distributed_mesh_2D'
       print *, 'Exiting...'
       STOP
    end if
    my_rank = sll_get_collective_rank(get_layout_2D_collective(layout))
    i_min = get_layout_2D_i_min( layout, my_rank )
    i_max = get_layout_2D_i_max( layout, my_rank )
    j_min = get_layout_2D_j_min( layout, my_rank )
    j_max = get_layout_2D_j_max( layout, my_rank )
    loc_sz_i = i_max - i_min + 1
    loc_sz_j = j_max - j_min + 1
  end subroutine compute_local_sizes_2d

        
    

end program test_io_parallel
