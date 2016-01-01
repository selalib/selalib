!> Unit test for parallel output
program test_io_parallel

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

  use iso_fortran_env, only: &
    output_unit

  use sll_m_collective, only: &
    sll_s_boot_collective, &
    sll_s_collective_barrier, &
    sll_f_get_collective_rank, &
    sll_f_get_collective_size, &
    sll_s_halt_collective, &
    sll_v_world_collective

  use sll_m_gnuplot_parallel, only: &
    sll_s_gnuplot_curv_2d_parallel, &
    sll_s_gnuplot_rect_2d_parallel

  use sll_m_hdf5_io_parallel, only: &
    sll_o_hdf5_file_create, &
    sll_o_hdf5_write_array

  use sll_m_hdf5_io_serial, only: &
    sll_o_hdf5_file_close

  use sll_m_remapper, only: &
    sll_o_compute_local_sizes, &
    sll_o_get_layout_i_min, &
    sll_o_get_layout_j_min, &
    sll_o_get_layout_k_min, &
    sll_o_initialize_layout_with_distributed_array, &
    sll_t_layout_2d, &
    sll_t_layout_3d, &
    sll_o_local_to_global, &
    sll_f_new_layout_2d, &
    sll_f_new_layout_3d, &
    sll_o_delete, &
    sll_o_view_lims

  use sll_m_utilities, only: &
    sll_f_is_power_of_two

  use sll_m_xdmf_parallel, only: &
    sll_o_xdmf_close, &
    sll_o_xdmf_open, &
    sll_o_xdmf_write_array

  use sll_m_xml_io, only: &
    sll_o_xml_field, &
    sll_s_xml_file_close, &
    sll_s_xml_file_create, &
    sll_o_xml_grid_geometry

  use sll_mpi, only: &
    mpi_info_null, &
    mpi_thread_single, &
    mpi_wtime

#ifndef NOHDF5
  use hdf5, only: &
    hid_t, &
    hsize_t, &
    hssize_t

#endif
  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

sll_int32   :: myrank
sll_int64   :: colsz 
sll_int32   :: comm, info
sll_int32   :: error
sll_int32   :: i, j, k
sll_real64  :: tcpu1
sll_real64  :: tcpu2

character(len=8), parameter :: xfile = "xdata.h5" !< File x coordinates
character(len=8), parameter :: yfile = "ydata.h5" !< File y coordinates
character(len=8), parameter :: zfile = "zdata.h5" !< File z coordinates
character(len=8), parameter :: xdset = "xdataset" !< x dataset name
character(len=8), parameter :: ydset = "ydataset" !< y dataset name
character(len=8), parameter :: zdset = "zdataset" !< z dataset name


! Boot parallel environment
call sll_s_boot_collective(MPI_THREAD_SINGLE)

colsz  = int(sll_f_get_collective_size(sll_v_world_collective),i64)
myrank = sll_f_get_collective_rank(sll_v_world_collective)
comm   = sll_v_world_collective%comm
info   = MPI_INFO_NULL

if( myrank .eq. 0) then
   print *, ' '
   print *, '--------------- HDF5 parallel test ---------------------'
   print *, ' '
   print"('Running a test on ',i4,' processes')", colsz
   flush( output_unit )
end if

if (.not. sll_f_is_power_of_two(colsz)) then     
   print *, 'This test needs to run in a number of processes which is ',&
        'a power of 2.'
   call sll_s_halt_collective()
   stop
end if

tcpu1 = MPI_WTIME()

call plot_layout2d()
call plot_layout3d()

tcpu2 = MPI_WTIME()
if (myrank == 0) &
   write(*,"(//10x,' Temps CPU = ', G15.3, ' sec' )") (tcpu2-tcpu1)*colsz


print *, myrank, 'PASSED'

call sll_s_halt_collective()
  
contains

!> Take a 2D array of dimensions ni*nj where ni, nj are the dimensions of
!> the full array.
!> @internal [example]
 subroutine plot_layout2d()

  sll_int32                :: mx, my    ! Local sizes
  sll_int32                :: npi, npj
  sll_int32                :: gi, gj
 
  sll_int32, dimension(2)  :: global_indices
  type(sll_t_layout_2d), pointer :: layout
  
  real(8), dimension(:,:), allocatable :: xdata, ydata, zdata
  sll_int32      :: xml_id
  sll_int32, parameter    :: nx = 64
  sll_int32, parameter    :: ny = 32
#ifndef NOHDF5
  integer(HID_T) :: file_id
  integer(HSIZE_T), dimension(2) :: datadims = (/int(nx,HSIZE_T),int(ny,HSIZE_T)/)
  integer(HSSIZE_T), dimension(2) :: offset 
#else
  sll_int32, dimension(2) :: offset 
#endif
  character(len=4) :: prefix = "mesh"
  
  layout => sll_f_new_layout_2d( sll_v_world_collective )        

  call two_power_rand_factorization(colsz, npi, npj)
  
  if( myrank .eq. 0 ) then
     print *, '2d layout configuration: ', npi, npj
  end if
  
  call sll_o_initialize_layout_with_distributed_array( &
       nx, ny, npi, npj, layout )
       
  call sll_s_collective_barrier(sll_v_world_collective)
  
  call sll_o_compute_local_sizes( layout, mx, my)        
  
  SLL_ALLOCATE(xdata(mx,my),error)
  SLL_ALLOCATE(ydata(mx,my),error)
  SLL_ALLOCATE(zdata(mx,my),error)
  
  do j = 1, my
     do i = 1, mx
        global_indices =  sll_o_local_to_global( layout, (/i, j/) )
        gi = global_indices(1)
        gj = global_indices(2)
        xdata(i,j) = real(myrank,f64) !float(gi-1)!/(nx-1)
        ydata(i,j) = real(gj-1,f64)!/(ny-1)
        zdata(i,j) = (myrank+1) * xdata(i,j) * ydata(i,j)
     end do
  end do
  
#ifdef NOHDF5
#define HSSIZE_T i32
#endif
  offset(1) =  int(sll_o_get_layout_i_min( layout, myrank ) - 1, HSSIZE_T)
  offset(2) =  int(sll_o_get_layout_j_min( layout, myrank ) - 1, HSSIZE_T)

  !Gnuplot output
  call sll_s_gnuplot_rect_2d_parallel(dble(offset(1)), dble(1), &
                                    dble(offset(2)), dble(1), &
                                    mx, my, &
                                    zdata, "rect_mesh", 1, error)  

  call sll_s_gnuplot_curv_2d_parallel(xdata, ydata, zdata, "curv_mesh", 1, error)  
  
  
#ifndef NOHDF5
  !Begin high level version

  call sll_o_xdmf_open(myrank,"zdata.xmf",prefix,nx,ny,xml_id,error)
  call sll_o_xdmf_write_array(prefix,datadims,offset,xdata,'x1',error)
  call sll_o_xdmf_write_array(prefix,datadims,offset,ydata,'x2',error)
  call sll_o_xdmf_write_array(prefix,datadims,offset,zdata,"x3",error,xml_id,"Node")
  call sll_o_xdmf_close(xml_id,error)

  !End high level version

!---------------------------------------------------------------------------------!

  !Begin low level version

  call sll_o_hdf5_file_create(xfile,comm,file_id, error)
  call sll_o_hdf5_write_array(file_id,datadims,offset,xdata,xdset,error)
  call sll_o_hdf5_file_close(file_id,error)

  
  call sll_o_hdf5_file_create(yfile,comm,file_id,error)
  call sll_o_hdf5_write_array(file_id,datadims,offset,ydata,ydset,error)
  call sll_o_hdf5_file_close(file_id,error)
  
  call sll_o_hdf5_file_create(zfile,comm,file_id,error)
  call sll_o_hdf5_write_array(file_id,datadims,offset,zdata,zdset,error)
  call sll_o_hdf5_file_close(file_id,error)

  if (myrank == 0) then
  
     call sll_s_xml_file_create("layout2d.xmf",xml_id,error)
     call sll_o_xml_grid_geometry(xml_id, xfile, nx, yfile, ny, &
                                xdset, ydset, 'Uniform' )
     call sll_o_xml_field(xml_id,'values', "zdata.h5:/zdataset",nx,ny,'HDF','Node')
     call sll_s_xml_file_close(xml_id,error)
     print *, 'Printing 2D layout: '
     call sll_o_view_lims( layout )
     print *, '--------------------'

  end if

#endif

  !End low level version

  call sll_o_delete( layout )
  
 end subroutine plot_layout2d
!>@internal [example]

 subroutine plot_layout3d()

  ! Test of the 3D remapper takes a 3D array whose global size Nx*Ny*Nz,
  ! distributed among NPi*NPj*NPk processors.
  sll_real64, dimension(:,:,:), allocatable :: local_array
  ! Take a 3D array of dimensions ni*nj*nk
  ! ni, nj, nk: global sizes
  ! Local sizes
  sll_int32                                   :: loc_sz_i_init
  sll_int32                                   :: loc_sz_j_init
  sll_int32                                   :: loc_sz_k_init

  ! the process mesh
  sll_int32                                   :: npi
  sll_int32                                   :: npj
  sll_int32                                   :: npk
  sll_int32                                 :: gi, gj, gk

  type(sll_t_layout_3d), pointer                  :: layout

  sll_int32, dimension(3)                   :: global_indices

  sll_int32      :: xml_id
  sll_int32, PARAMETER :: rank = 3
  sll_int32 , parameter                       :: ni = 32
  sll_int32 , parameter                       :: nj = 64
  sll_int32 , parameter                       :: nk = 128
#ifndef NOHDF5
  integer(HID_T) :: file_id       ! File identifier 
  integer(HSIZE_T), dimension(3) :: datadims = (/int(ni,HSIZE_T),int(nj,HSIZE_T),int(nk,HSIZE_T)/) ! Dataset dimensions.
  integer(HSSIZE_T), dimension(rank) :: offset 
#else
  sll_int32, dimension(rank) :: offset 
#endif



  real(8), dimension(:,:,:), allocatable :: xdata, ydata, zdata

  layout  => sll_f_new_layout_3d( sll_v_world_collective )        
  call two_power_rand_factorization(colsz, npi, npj, npk)

  if( myrank .eq. 0 ) &
     print *, '3D layout configuration: ', npi, npj, npk

  call sll_o_initialize_layout_with_distributed_array( &
                     ni, nj, nk, npi, npj, npk, layout)
     
  call sll_o_compute_local_sizes( layout,       &
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
           global_indices =  sll_o_local_to_global( layout, (/i, j, k/) )
           gi = global_indices(1)
           gj = global_indices(2)
           gk = global_indices(3)
           local_array(i,j,k) = real(myrank,f64) !gi + (gj-1)*ni + (gk-1)*ni*nj
           xdata(i,j,k) = dble(gi-1) / dble(ni-1)
           ydata(i,j,k) = dble(gj-1) / dble(nj-1)
           zdata(i,j,k) = dble(gk-1) / dble(nk-1)
        enddo
     enddo
  enddo

#ifdef NOHDF5
#define HSSIZE_T i32
#endif

  offset(1) = int(sll_o_get_layout_i_min( layout, myrank ) - 1, HSSIZE_T)
  offset(2) = int(sll_o_get_layout_j_min( layout, myrank ) - 1, HSSIZE_T)
  offset(3) = int(sll_o_get_layout_k_min( layout, myrank ) - 1, HSSIZE_T)

#ifndef NOHDF5

  !Begin high level version

  call sll_o_xdmf_open(myrank,"3d-data.xmf","mesh3d",ni,nj,nk,xml_id,error)
  call sll_o_xdmf_write_array("mesh3d",datadims,offset,xdata,'x1',error)
  call sll_o_xdmf_write_array("mesh3d",datadims,offset,ydata,'x2',error)
  call sll_o_xdmf_write_array("mesh3d",datadims,offset,zdata,"x3",error,xml_id,"Node")
  call sll_o_xdmf_close(xml_id,error)

  !End high level version

  call sll_o_hdf5_file_create('layout3d-x.h5',comm,file_id,error)
  call sll_o_hdf5_write_array(file_id,datadims,offset,xdata,'x',error)
  call sll_o_hdf5_file_close(file_id, error)

  call sll_o_hdf5_file_create('layout3d-y.h5',comm,file_id, error)
  call sll_o_hdf5_write_array(file_id,datadims,offset,ydata,'y',error)
  call sll_o_hdf5_file_close(file_id, error)

  call sll_o_hdf5_file_create('layout3d-z.h5',comm,file_id, error)
  call sll_o_hdf5_write_array(file_id,datadims,offset,zdata,'z',error)
  call sll_o_hdf5_file_close(file_id, error)

  call sll_o_hdf5_file_create('layout3d.h5',comm,file_id, error)
  call sll_o_hdf5_write_array(file_id,datadims,offset,local_array,'array',error)
  call sll_o_hdf5_file_close(file_id,error)

  if (myrank == 0) then
     call sll_s_xml_file_create("layout3d.xmf",xml_id,error)
     call sll_o_xml_grid_geometry(xml_id, 'layout3d-x.h5', ni, &
                                        'layout3d-y.h5', nj, &
                                        'layout3d-z.h5', nk, &
                                        'x', 'y', 'z', &
                                        'Uniform' )
     call sll_o_xml_field(xml_id,'values', "layout3d.h5:/array", &
                        ni,nj,nk,'HDF','Node')
     call sll_s_xml_file_close(xml_id,error)

     print *, 'Printing 3D layout: '
     call sll_o_view_lims( layout )
     print *, '--------------------'
  end if

#endif

  call sll_s_collective_barrier(sll_v_world_collective)
  
  call sll_o_delete( layout )
  SLL_DEALLOCATE_ARRAY(local_array, error)

  end subroutine plot_layout3d

  subroutine two_power_rand_factorization(n, n1, n2, n3)
    sll_int64, intent(in) :: n
    integer, intent(out) ::n1, n2
    integer, intent(out), optional :: n3
    integer   :: expo, expo1, expo2, expo3
    sll_real64                :: rand_real
    if (.not.sll_f_is_power_of_two(colsz)) then   
       print*, 'The number of processors must be a power of 2'
       call sll_s_halt_collective()
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
