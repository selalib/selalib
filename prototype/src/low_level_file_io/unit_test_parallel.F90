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

!call plot_layout2d()
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
  call factorize_in_random_2powers_2d(colsz, npi, npj)
  
  if( myrank .eq. 0 ) then
     print *, 'source configuration: ', npi, npj
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
        zdata(i,j) = myrank * xdata(i,j) * ydata(i,j)
     end do
  end do
  
  offset(1) =  get_layout_2D_i_min( layout, myrank ) - 1
  offset(2) =  get_layout_2D_j_min( layout, myrank ) - 1
  
  !data stored in separate file
  !convenient to use xml interface
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
  
     call sll_xml_file_create("layout_2d_0.xmf",file_id,error)
     call sll_xml_grid_geometry(file_id, xfile, xdset, nx, yfile, ydset, ny )
     call sll_xml_field(file_id,'Z', "zdata.h5:/zdataset",nx,ny,'HDF','Node')
     call sll_xml_file_close(file_id,error)

  end if
  
  !data stored in the same file
  call my_file%create('layout_2d.h5', error)
  call my_file%write_array(datadims,offset,xdata,'x',error)
  call my_file%write_array(datadims,offset,ydata,'y',error)
  call my_file%write_array(datadims,offset,zdata,'z',error)
  call my_file%close(error)

  !write xml file
  if (myrank == 0) then

     call sll_xml_file_create("layout_2d_1.xmf",file_id,error)
     write(file_id,"(a)")"<Grid Name='mesh' GridType='Uniform'>"
     write(file_id,"(a,2i5,a)") &
              "<Topology TopologyType='2DSMesh' NumberOfElements='", ny,nx,"'/>"
     write(file_id,"(a)")"<Geometry GeometryType='X_Y'>"
     write(file_id,"(a,2i5,a)")"<DataItem Dimensions='",ny,nx, &
             "' NumberType='Float' Precision='8' Format='HDF'>"
     write(file_id,"(a)")"layout.h5:/x"
     write(file_id,"(a)")"</DataItem>"
     write(file_id,"(a,2i5,a)")"<DataItem Dimensions='",ny,nx, &
             "' NumberType='Float' Precision='8' Format='HDF'>"
     write(file_id,"(a)")"layout.h5:/y"
     write(file_id,"(a)")"</DataItem>"
     write(file_id,"(a)")"</Geometry>"
     write(file_id,"(a)") &
             "<Attribute Name='z_values' AttributeType='Scalar' Center='Node'>"
     write(file_id,"(a,2i5,a)")"<DataItem Dimensions='",ny,nx, &
             "' NumberType='Float' Precision='8' Format='HDF'>"
     write(file_id,"(a)")"layout.h5:/z"
     write(file_id,"(a)")"</DataItem>"
     write(file_id,"(a)")"</Attribute>"
     call sll_xml_file_close(file_id,error)
  
  end if

  call delete_layout_2D( layout )
  
 end subroutine plot_layout2d

 subroutine plot_layout3d()

  ! Test of the 3D remapper takes a 3D array whose global size Nx*Ny*Nz,
  ! distributed among NPi*NPj*NPk processors.
  sll_real64, dimension(:,:,:), allocatable :: local_array1
  sll_real64, dimension(:,:,:), allocatable :: local_array2
  ! Take a 3D array of dimensions ni*nj*nk
  ! ni, nj, nk: global sizes
  integer , parameter                       :: ni = 64
  integer , parameter                       :: nj = 64
  integer , parameter                       :: nk = 64
  ! Local sizes
  integer                                   :: loc_sz_i_init
  integer                                   :: loc_sz_j_init
  integer                                   :: loc_sz_k_init
  integer                                   :: loc_sz_i_final
  integer                                   :: loc_sz_j_final
  integer                                   :: loc_sz_k_final

  ! the process mesh
  integer                                   :: npi
  integer                                   :: npj
  integer                                   :: npk
  sll_int32                                 :: gi, gj, gk
  ! Remap stuff
  type(layout_3D), pointer                  :: layout1
  type(layout_3D), pointer                  :: layout2
  type(remap_plan_3D), pointer              :: rmp3

  sll_int32, dimension(3)                   :: global_indices
  type(phdf5_file)                          :: file_layout1
  type(phdf5_file)                          :: file_layout2
  integer                                   :: istart, jstart, kstart

  integer(HID_T) :: file_id       ! File identifier 
  integer(HID_T) :: dset_id       ! Dataset identifier 
  integer(HID_T) :: filespace     ! Dataspace identifier in file 
  integer(HID_T) :: plist_id      ! Property list identifier 

  integer(HSIZE_T), dimension(3) :: dimsf = (/ni,nj,nk/) ! Dataset dimensions.
  integer(HSIZE_T), dimension(3) :: dimsfi
  integer(HSIZE_T), dimension(3) :: chunk_dims

  integer, PARAMETER :: rank = 3

  integer(HSIZE_T),  dimension(rank) :: count  
  integer(HSSIZE_T), dimension(rank) :: offset 
  integer(HSIZE_T),  dimension(rank) :: stride
  integer(HSIZE_T),  dimension(rank) :: block
  integer(HID_T) :: memspace 

  real(8), dimension(:,:,:), allocatable :: xdata, ydata, zdata


  layout1  => new_layout_3D( sll_world_collective )        
  call two_power_rand_factorization(colsz, npi, npj, npk)

  if( myrank .eq. 0 ) &
     print *, '3D layout configuration: ', npi, npj, npk

  call initialize_layout_with_distributed_3D_array( &
                     ni, nj, nk, npi, npj, npk, layout1 )
     
  call compute_local_sizes( layout1,       &
                            loc_sz_i_init, &
                            loc_sz_j_init, &
                            loc_sz_k_init )        

  SLL_ALLOCATE(local_array1(loc_sz_i_init,loc_sz_j_init,loc_sz_k_init),error)
  SLL_ALLOCATE(xdata(loc_sz_i_init,loc_sz_j_init,loc_sz_k_init),error)
  SLL_ALLOCATE(ydata(loc_sz_i_init,loc_sz_j_init,loc_sz_k_init),error)
  SLL_ALLOCATE(zdata(loc_sz_i_init,loc_sz_j_init,loc_sz_k_init),error)
 
  ! initialize the local data    
  do k=1,loc_sz_k_init
     do j=1,loc_sz_j_init 
        do i=1,loc_sz_i_init
           global_indices =  local_to_global_3D( layout1, (/i, j, k/) )
           gi = global_indices(1)
           gj = global_indices(2)
           gk = global_indices(3)
           local_array1(i,j,k) = myrank !gi + (gj-1)*ni + (gk-1)*ni*nj
           xdata(i,j,k) = float(gi-1) / (ni-1)
           ydata(i,j,k) = float(gj-1) / (nj-1)
           zdata(i,j,k) = float(gk-1) / (nk-1)
        enddo
     enddo
  enddo

  chunk_dims(1) = loc_sz_i_init
  chunk_dims(2) = loc_sz_j_init
  chunk_dims(3) = loc_sz_k_init

  stride = 1 
  count =  1 
  block = chunk_dims

  istart = get_layout_3D_i_min( layout1, myrank ) - 1
  jstart = get_layout_3D_j_min( layout1, myrank ) - 1
  kstart = get_layout_3D_k_min( layout1, myrank ) - 1

  offset(1) = istart
  offset(2) = jstart
  offset(3) = kstart

  call file_layout1%create('mesh3d-x.h5',error)
  call file_layout1%write_array(dimsf,offset,xdata,'x',error)
  call file_layout1%close(error)
  call file_layout1%create('mesh3d-y.h5',error)
  call file_layout1%write_array(dimsf,offset,ydata,'y',error)
  call file_layout1%close(error)
  call file_layout1%create('mesh3d-z.h5',error)
  call file_layout1%write_array(dimsf,offset,zdata,'z',error)
  call file_layout1%close(error)

  if (myrank == 0) then
  
     call sll_xml_file_create("layout_3d_0.xmf",file_id,error)
     call sll_xml_grid_geometry(file_id, 'mesh3d-x.h5', 'x', ni, &
                                         'mesh3d-y.h5', 'y', nj, &
                                         'mesh3d-z.h5', 'z', nk )
     call sll_xml_field(file_id,'array_values', "layout_0.h5:/array", &
                        ni,nj,nk,'HDF','Node')
     call sll_xml_file_close(file_id,error)

  end if

  !write(*,"(i2,a16,3i3)") myrank, " offset = ", offset
  !write(*,"(i2,a16,3i3)") myrank, " blocks = ", block

  !write hdf5 file with low with native API
  call h5open_f(error) 
  call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
  call h5pset_fapl_mpio_f(plist_id, comm, info, error)
  call h5fcreate_f("layout_0.h5",H5F_ACC_TRUNC_F,file_id, &
                   error,access_prp = plist_id)
  call h5pclose_f(plist_id, error)
  call h5screate_simple_f(rank, dimsf, filespace, error)
  call h5screate_simple_f(rank, chunk_dims, memspace, error)
  call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)
  call h5pset_chunk_f(plist_id, rank, chunk_dims, error)
  call h5dcreate_f(file_id, "array", H5T_NATIVE_DOUBLE, filespace, &
                   dset_id, error, plist_id)
  call h5sclose_f(filespace, error)
  call h5dget_space_f(dset_id, filespace, error)
  call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, &
                              count, error, stride, block)
  call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error) 
  call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
  call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, local_array1, dimsfi, error, &
                  file_space_id = filespace, mem_space_id = memspace, &
                  xfer_prp = plist_id)
  call h5sclose_f(filespace, error)
  call h5sclose_f(memspace, error)
  call h5dclose_f(dset_id, error)
  call h5pclose_f(plist_id, error)
  call h5fclose_f(file_id, error)
  call h5close_f(error)

  layout2  => new_layout_3D( sll_world_collective )
  call two_power_rand_factorization(colsz, npi, npj, npk)
  if( myrank .eq. 0 ) &
     print *, 'target configuration: ', npi, npj, npk

  call initialize_layout_with_distributed_3D_array( &
       ni, nj, nk, npi, npj, npk, layout2 )

  call reorganize_randomly(layout2)
     
  call compute_local_sizes(layout2,        &
                           loc_sz_i_final, &
                           loc_sz_j_final, &
                           loc_sz_k_final )

  SLL_ALLOCATE(local_array2(loc_sz_i_final, loc_sz_j_final, loc_sz_k_final), error )
    
  rmp3 => NEW_REMAP_PLAN_3D( layout1, layout2, local_array1)     

  call apply_remap_3D( rmp3, local_array1, local_array2 ) 

  if( myrank .eq. 0 ) then
     print *, 'Printing layout1: '
     call sll_view_lims_3D( layout1 )
     print *, 'Printing layout2: '
     call sll_view_lims_3D( layout2 )
  end if

  !write with slelalib interface sll_hdf5_io_parallel
  call file_layout1%create('layout_1.h5',error)
  call file_layout1%write_array(dimsf,offset,local_array1,'array',error)
  call file_layout1%close(error)

  call file_layout2%create('layout_2.h5',error)
  call file_layout2%write_array(dimsf,offset,local_array2,'array',error)
  call file_layout2%close(error)

  call sll_collective_barrier(sll_world_collective)
  
  call delete_layout_3D( layout1 )
  call delete_layout_3D( layout2 )
  SLL_DEALLOCATE_ARRAY(local_array1, error)
  SLL_DEALLOCATE_ARRAY(local_array2, error)

  end subroutine plot_layout3d

  subroutine reorganize_randomly(layout)
    implicit none
    type(layout_3D), pointer   :: layout
    integer                    :: i, colsz, proc_n, proc_p
    real                       :: rand_real
    colsz = sll_get_num_nodes(layout)
    do i=0, colsz-1
       call random_number(rand_real)
       proc_n = int(rand_real*(colsz-1))
       call random_number(rand_real)
       proc_p = int(rand_real*(colsz-1))
       call swap_box_3D(proc_n, proc_p, layout)
    enddo    
  end subroutine reorganize_randomly

  subroutine swap_box_3D(proc_n, proc_p, layout)
    implicit none
    integer                  :: proc_n, proc_p
    type(layout_3D), pointer :: layout
    integer                  :: i_min_n, i_max_n, j_min_n, j_max_n, &
    k_min_n, k_max_n, i_min_p, i_max_p, j_min_p, j_max_p, k_min_p, k_max_p    
    ! Get proc_n contents from layout
    i_min_n = get_layout_3D_i_min( layout, proc_n )
    i_max_n = get_layout_3D_i_max( layout, proc_n )
    j_min_n = get_layout_3D_j_min( layout, proc_n )
    j_max_n = get_layout_3D_j_max( layout, proc_n )
    k_min_n = get_layout_3D_k_min( layout, proc_n )
    k_max_n = get_layout_3D_k_max( layout, proc_n )
    ! Get proc_p contents from layout
    i_min_p = get_layout_3D_i_min( layout, proc_p )
    i_max_p = get_layout_3D_i_max( layout, proc_p )
    j_min_p = get_layout_3D_j_min( layout, proc_p )
    j_max_p = get_layout_3D_j_max( layout, proc_p )
    k_min_p = get_layout_3D_k_min( layout, proc_p )
    k_max_p = get_layout_3D_k_max( layout, proc_p )
    ! Set proc_n contents in layout
    call set_layout_3D_i_min( layout, proc_n, i_min_p )
    call set_layout_3D_i_max( layout, proc_n, i_max_p)
    call set_layout_3D_j_min( layout, proc_n, j_min_p )
    call set_layout_3D_j_max( layout, proc_n, j_max_p )
    call set_layout_3D_k_min( layout, proc_n, k_min_p )
    call set_layout_3D_k_max( layout, proc_n, k_max_p )
    ! Set proc_p contents in layout
    call set_layout_3D_i_min( layout, proc_p, i_min_n )
    call set_layout_3D_i_max( layout, proc_p, i_max_n )
    call set_layout_3D_j_min( layout, proc_p, j_min_n )
    call set_layout_3D_j_max( layout, proc_p, j_max_n )
    call set_layout_3D_k_min( layout, proc_p, k_min_n )
    call set_layout_3D_k_max( layout, proc_p, k_max_n )   
  end subroutine swap_box_3D
  
  
  subroutine split_interval_randomly(n, num_halvings, ans)
    integer, intent(in) :: n
    integer, intent(in) :: num_halvings
    integer, dimension(:), pointer :: ans
    integer :: load_i = 1
    SLL_ALLOCATE(ans(2*(2**num_halvings)),error)
    call split_interval_randomly_aux( 1, n, 0, num_halvings, load_i, ans )
  end subroutine split_interval_randomly
      
  recursive subroutine split_interval_randomly_aux( lo, hi, gen, lim, loadi, ans )
    intrinsic :: random_number
    integer, intent(in)                :: lo
    integer, intent(in)                :: hi
    integer, intent(in)                :: gen ! splitting generation
    integer, intent(in)                :: lim ! maximum number of splittings
    integer :: mid
    integer, parameter                 :: spread = 25
    ! dangerous... the following should be found in the enclosing scope...
    integer, intent(inout) :: loadi
    integer, dimension(:), pointer :: ans
    if( (hi-lo).eq. 1 ) then
       ans(loadi) = lo
       ans(loadi+1) = hi
       loadi = loadi + 2
    else if (gen .eq. lim) then
       ans(loadi) = lo
       ans(loadi+1) = hi
       loadi = loadi + 2
    else
       ! call random_number(rand)
       ! decided to use the gaussian because the uniform deviate can give
       ! a number very close to the 0 or one and such a small interval is
       ! hard to subdivide. In such case one may end up with less intervals
       ! and this is a problem since the number of processors is fixed from
       ! the beginning.
       mid = int((hi+lo)/2 + gaussian_dev()*spread)
       !       mid = int(rand*(hi-lo))+lo
       call split_interval_randomly_aux( lo,   mid, gen+1, lim, loadi, ans )
       call split_interval_randomly_aux( mid+1, hi, gen+1, lim, loadi, ans )
    end if
  end subroutine split_interval_randomly_aux
      
  function gaussian_dev()
    intrinsic :: random_number
    real :: gaussian_dev
    real :: v1, v2
    real :: ran1, ran2
    real :: rsq
    real :: fac
    do
       call random_number(ran1)
       call random_number(ran2)
       v1 = 2.0*ran1-1.0
       v2 = 2.0*ran2-1.0
       rsq = v1*v1 + v2*v2
       if( (rsq .lt. 1.0) .and. (rsq .gt. 0.0) ) exit
    end do
    fac = sqrt(-2.0*log(rsq)/rsq)
    gaussian_dev = v1*fac
  end function gaussian_dev

  subroutine two_power_rand_factorization(n, n1, n2, n3)
    sll_int64, intent(in) :: n
    integer, intent(out) ::n1, n2, n3
    integer   :: expo, expo1, expo2, expo3
    sll_real64                :: rand_real
    if (.not.is_power_of_two(colsz)) then   
       print*, 'The number of processors must be a power of 2'
       stop
    endif 
    expo = int(log(real(n))/log(2.))  
    call random_number(rand_real)
    expo1 = int(rand_real*expo)
    call random_number(rand_real)
    expo2 = int(rand_real*(expo-expo1))
    expo3 = expo - (expo1+expo2)
    n1 = 2**expo1
    n2 = 2**expo2
    n3 = 2**expo3
  end subroutine two_power_rand_factorization

  subroutine factorize_in_random_2powers_2d(n, n1, n2)
    sll_int64, intent(in) :: n
    integer, intent(out)  ::n1, n2
    integer   :: expo, expo1, expo2
    sll_real64                :: rand_real
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
  
end program test_io_parallel
