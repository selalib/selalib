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

sll_int32   :: myrank
sll_int64   :: colsz 
sll_int32   :: comm, info
sll_int32   :: file_id, error
sll_int32   :: i, j, k
sll_real64  :: tcpu1
sll_real64  :: tcpu2

character(len=8), parameter :: xfile = "xdata.h5" ! File name
character(len=8), parameter :: yfile = "ydata.h5" ! File name
character(len=8), parameter :: zfile = "zdata.h5" ! File name
character(len=8), parameter :: xdset = "xdataset" ! Dataset name
character(len=8), parameter :: ydset = "ydataset" ! Dataset name
character(len=8), parameter :: zdset = "zdataset" ! Dataset name


! Boot parallel environment
call sll_boot_collective()

colsz  = sll_get_collective_size(sll_world_collective)
myrank = sll_get_collective_rank(sll_world_collective)
comm   = sll_world_collective%comm
info   = MPI_INFO_NULL

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

tcpu1 = MPI_WTIME()

call plot_layout2d()
call plot_layout3d()

tcpu2 = MPI_WTIME()
if (myrank == 0) &
   write(*,"(//10x,' Temps CPU = ', G15.3, ' sec' )") (tcpu2-tcpu1)*colsz

call sll_halt_collective()

if( myrank .eq. 0) print *, 'PASSED'
  
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
  type(phdf5_file) :: my_file
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
  
  call my_file%create(xfile, error)
  call my_file%write_array(datadims,offset,xdata,xdset,error)
  call my_file%close(error)
  
  call my_file%create(yfile, error)
  call my_file%write_array(datadims,offset,ydata,ydset,error)
  call my_file%close(error)
  
  call my_file%create(zfile, error)
  call my_file%write_array(datadims,offset,zdata,zdset,error)
  call my_file%close(error)

  if (myrank == 0) then
  
     call sll_xml_file_create("layout2d.xmf",file_id,error)
     call sll_xml_grid_geometry(file_id, xfile, xdset, nx, yfile, ydset, ny )
     call sll_xml_field(file_id,'values', "zdata.h5:/zdataset",nx,ny,'HDF','Node')
     call sll_xml_file_close(file_id,error)
     print *, 'Printing 2D layout: '
     call sll_view_lims_2D( layout )
     print *, '--------------------'

  end if
  
  call delete_layout_2D( layout )
  
 end subroutine plot_layout2d

 subroutine plot_layout3d()

  ! Test of the 3D remapper takes a 3D array whose global size Nx*Ny*Nz,
  ! distributed among NPi*NPj*NPk processors.
  sll_real64, dimension(:,:,:), allocatable :: local_array
  ! Take a 3D array of dimensions ni*nj*nk
  ! ni, nj, nk: global sizes
  integer , parameter                       :: ni = 32
  integer , parameter                       :: nj = 64
  integer , parameter                       :: nk = 128
  ! Local sizes
  integer                                   :: loc_sz_i_init
  integer                                   :: loc_sz_j_init
  integer                                   :: loc_sz_k_init

  ! the process mesh
  integer                                   :: npi
  integer                                   :: npj
  integer                                   :: npk
  sll_int32                                 :: gi, gj, gk

  type(layout_3D), pointer                  :: layout

  sll_int32, dimension(3)                   :: global_indices
  type(phdf5_file)                          :: hdf_file

  integer(HID_T) :: file_id       ! File identifier 

  integer(HSIZE_T), dimension(3) :: datadims = (/ni,nj,nk/) ! Dataset dimensions.

  integer, PARAMETER :: rank = 3

  integer(HSSIZE_T), dimension(rank) :: offset 

  real(8), dimension(:,:,:), allocatable :: xdata, ydata, zdata

  layout  => new_layout_3D( sll_world_collective )        
  call two_power_rand_factorization(colsz, npi, npj, npk)

  if( myrank .eq. 0 ) &
     print *, '3D layout configuration: ', npi, npj, npk

  call initialize_layout_with_distributed_3D_array( &
                     ni, nj, nk, npi, npj, npk, layout)
     
  call compute_local_sizes( layout,       &
                            loc_sz_i_init, &
                            loc_sz_j_init, &
                            loc_sz_k_init )        

  SLL_ALLOCATE(local_array(loc_sz_i_init,loc_sz_j_init,loc_sz_k_init),error)
  SLL_ALLOCATE(xdata(loc_sz_i_init,loc_sz_j_init,loc_sz_k_init),error)
  SLL_ALLOCATE(ydata(loc_sz_i_init,loc_sz_j_init,loc_sz_k_init),error)
  SLL_ALLOCATE(zdata(loc_sz_i_init,loc_sz_j_init,loc_sz_k_init),error)
 
  ! initialize the local data    
  do k=1,loc_sz_k_init
     do j=1,loc_sz_j_init 
        do i=1,loc_sz_i_init
           global_indices =  local_to_global_3D( layout, (/i, j, k/) )
           gi = global_indices(1)
           gj = global_indices(2)
           gk = global_indices(3)
           local_array(i,j,k) = myrank !gi + (gj-1)*ni + (gk-1)*ni*nj
           xdata(i,j,k) = float(gi-1) / (ni-1)
           ydata(i,j,k) = float(gj-1) / (nj-1)
           zdata(i,j,k) = float(gk-1) / (nk-1)
        enddo
     enddo
  enddo

  offset(1) = get_layout_3D_i_min( layout, myrank ) - 1
  offset(2) = get_layout_3D_j_min( layout, myrank ) - 1
  offset(3) = get_layout_3D_k_min( layout, myrank ) - 1

  call hdf_file%create('layout3d-x.h5',error)
  call hdf_file%write_array(datadims,offset,xdata,'x',error)
  call hdf_file%close(error)

  call hdf_file%create('layout3d-y.h5',error)
  call hdf_file%write_array(datadims,offset,ydata,'y',error)
  call hdf_file%close(error)

  call hdf_file%create('layout3d-z.h5',error)
  call hdf_file%write_array(datadims,offset,zdata,'z',error)
  call hdf_file%close(error)

  call hdf_file%create('layout3d.h5',error)
  call hdf_file%write_array(datadims,offset,local_array,'array',error)
  call hdf_file%close(error)

  if (myrank == 0) then
     call sll_xml_file_create("layout3d.xmf",file_id,error)
     call sll_xml_grid_geometry(file_id, 'layout3d-x.h5', 'x', ni, &
                                         'layout3d-y.h5', 'y', nj, &
                                         'layout3d-z.h5', 'z', nk )
     call sll_xml_field(file_id,'values', "layout3d.h5:/array", &
                        ni,nj,nk,'HDF','Node')
     call sll_xml_file_close(file_id,error)

     print *, 'Printing 3D layout: '
     call sll_view_lims_3D( layout )
     print *, '--------------------'
  end if


  call sll_collective_barrier(sll_world_collective)
  
  call delete_layout_3D( layout )
  SLL_DEALLOCATE_ARRAY(local_array, error)

  end subroutine plot_layout3d

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

end program test_io_parallel
