!===========================================================================
!> 3D scalar field in (x,y,eta3) for
!>  4D drift-kinetic hybrid simulation
!>
!> \date 2014-08-20
!> \author V. Grandgirard
!---------------------------------------------------------------------------
module field2d_VP4D_module
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_utilities.h"

  use input_VP4D_module
  use mesh_VP4D_module
  use sll_m_collective, only : sll_v_world_collective, &
      sll_f_get_collective_size, sll_f_get_collective_rank
  use sll_m_remapper

  implicit none

  !===========================================================================
  !> 
  !---------------------------------------------------------------------------
  type, public :: field2d_VP4D_t

    !> Name of the field
    character(len=30) :: name

    !> Boundary conditions
    type(boundary_conditions_2d_t):: bound_cond

    !> Parallel in (x1,x2)
    type(sll_t_layout_2d), pointer :: layout2d_parx1x2
    sll_real64, dimension(:,:), pointer :: val2d_parx1x2 

    !> Sequential in (x1,x2)
    sll_real64, dimension(:,:), pointer :: val2d_seqx1x2 

    !> For interpolations
    type(spline_degree_2d_t) :: spline_degree

  end type field2d_VP4D_t
  !---------------------------------------------------------------------------

contains

  !===========================================================================
  !> Field2d: Allocation 
  !---------------------------------------------------------------------------
  subroutine new_field2d_VP4D( &
      field2d, &
      field2d_name, &
      mesh4d, &
      bound_cond, &
      spline_degree )

    use sll_m_utilities, only : sll_f_is_even
    
    type(field2d_VP4D_t)          , intent(inout) :: field2d
    character(len=*)              , intent(in)    :: field2d_name
    type(mesh_VP4D_t)             , intent(in)    :: mesh4d
    type(boundary_conditions_2d_t), intent(in)    :: bound_cond
    type(spline_degree_2d_t)      , intent(in)    :: spline_degree

    !-> Local variables 
    sll_int32 :: ierr
    !--> For MPI parallelization
    sll_int32 :: world_size
    sll_int32 :: my_rank
    sll_int32 :: power2
    sll_int32 :: nproc_x1, nproc_x2
    !--> For parallel distribution
    sll_int32 :: npoints_x1, npoints_x2
    sll_int32 :: loc2d_sz_x1, loc2d_sz_x2

    !**********************************************************
    !*** Initialization of the name of the field            ***
    !**********************************************************
    field2d%name = trim(field2d_name)

    !**********************************************************
    !*** Initialization of the parallel distribution        ***
    !**********************************************************
    world_size = sll_f_get_collective_size(sll_v_world_collective)
    my_rank    = sll_f_get_collective_rank(sll_v_world_collective)
    power2 = int(log(real(world_size))/log(2.0))

    !**********************************************************************
    !*** Initialization of parallel layout of field2d in (x1,x2)        ***
    !***  directions                                                    *** 
    !***  (x1,x2) : parallelized layout                                 ***
    !**********************************************************************
    !--> special case N = 1, so power2 = 0
    if ( power2 == 0 ) then
       nproc_x1 = 1
       nproc_x2 = 1
    end if
    
    if ( sll_f_is_even(power2) ) then
       nproc_x1 = 2**(power2/2)
       nproc_x2 = 2**(power2/2)
    else 
       nproc_x1 = 2**((power2-1)/2)
       nproc_x2 = 2**((power2+1)/2)
    end if

    field2d%layout2d_parx1x2  => sll_f_new_layout_2d( sll_v_world_collective )
    npoints_x1 = mesh4d%eta1_eta2_mesh2d%num_cells1 + 1
    npoints_x2 = mesh4d%eta1_eta2_mesh2d%num_cells2 + 1

    call sll_o_initialize_layout_with_distributed_array( &
      npoints_x1, &
      npoints_x2, &
      nproc_x1, &
      nproc_x2, &
      field2d%layout2d_parx1x2 )
    call sll_o_compute_local_sizes( field2d%layout2d_parx1x2, &
      loc2d_sz_x1, &
      loc2d_sz_x2 )

    SLL_ALLOCATE(field2d%val2d_parx1x2(loc2d_sz_x1,loc2d_sz_x2), ierr )

    !*************************************************************
    !*** Initialization of the complete field2d                ***
    !***   (i.e sequential in (x1,x2))                         ***
    !*************************************************************
    SLL_ALLOCATE( field2d%val2d_seqx1x2(npoints_x1,npoints_x2), ierr )

    !*************************************************
    !*** Initialization of the boundary conditions ***
    !*************************************************
    field2d%bound_cond%left_eta1  = bound_cond%left_eta1
    field2d%bound_cond%right_eta1 = bound_cond%right_eta1
    field2d%bound_cond%left_eta2  = bound_cond%left_eta2
    field2d%bound_cond%right_eta2 = bound_cond%right_eta2

    !*************************************************
    !*** Initialization of the interpolators       ***
    !*************************************************
    !--> Initialize the spline degree
    field2d%spline_degree%eta1 = spline_degree%eta1
    field2d%spline_degree%eta2 = spline_degree%eta2

  end subroutine new_field2d_VP4D
  !---------------------------------------------------------------------------


  !===========================================================================
  !> Field2d: Deallocation 
  !---------------------------------------------------------------------------
  subroutine delete_field2d_VP4D( field2d )

    type(field2d_VP4D_t), intent(inout) :: field2d

    sll_int32 :: ierr

    call sll_o_delete( field2d%layout2d_parx1x2 )
    SLL_DEALLOCATE( field2d%val2d_seqx1x2, ierr )
    SLL_DEALLOCATE( field2d%val2d_parx1x2, ierr )

  end subroutine delete_field2d_VP4D
  !---------------------------------------------------------------------------


  !===========================================================================
  !> Field2d: remapping from (x1,x2) parallel to (x1,x2) sequential
  !---------------------------------------------------------------------------
  subroutine remap_parx1x2_to_seqx1x2_field2d_VP4D( field2d )

    use sll_m_collective
    use utils_VP4D_module

    type(field2d_VP4D_t), intent(inout) :: field2d

    sll_int32 :: loc_sz_x1, loc_sz_x2
    sll_int32 :: Nx1, Nx2
    sll_int32 :: world_size
    sll_int32 :: send_size
    sll_int32 :: ierr
    sll_real64, dimension(:), pointer :: send_buf
    sll_real64, dimension(:), pointer :: recv_buf
    sll_int32 , dimension(:), pointer :: recv_sz
    sll_int32 , dimension(:), pointer :: disps 

    world_size = sll_f_get_collective_size( sll_v_world_collective )
    SLL_ALLOCATE( recv_sz(world_size), ierr )
    SLL_ALLOCATE( disps(world_size), ierr )

    call sll_o_compute_local_sizes( field2d%layout2d_parx1x2, loc_sz_x1, loc_sz_x2 )
    SLL_ALLOCATE( send_buf(loc_sz_x1*loc_sz_x2), ierr )

    Nx1 = size(field2d%val2d_seqx1x2,1)
    Nx2 = size(field2d%val2d_seqx1x2,2)
    SLL_ALLOCATE( recv_buf(Nx1*Nx2), ierr )

    call load_buffer( field2d%layout2d_parx1x2, &
        field2d%val2d_parx1x2, send_buf )
    recv_sz(:) = receive_counts_array( field2d%layout2d_parx1x2, &
        world_size )
    send_size = size(send_buf)
    call compute_displacements_array( &
        field2d%layout2d_parx1x2, &
        world_size, &
        disps )
    call sll_s_collective_allgatherv_real64( &
         sll_v_world_collective, &
         send_buf, &
         send_size, &
         recv_sz, &
         disps, &
         recv_buf )
    call unload_buffer( field2d%layout2d_parx1x2, &
        recv_buf, field2d%val2d_seqx1x2 )
    
    SLL_DEALLOCATE( recv_sz, ierr )
    SLL_DEALLOCATE( disps, ierr )
    SLL_DEALLOCATE( send_buf, ierr )
    SLL_DEALLOCATE( recv_buf, ierr )

  end subroutine remap_parx1x2_to_seqx1x2_field2d_VP4D
  !---------------------------------------------------------------------------


  !===========================================================================
  !> Field2d: remapping from (x1,x2) parallel to (x1,x2) sequential
  !---------------------------------------------------------------------------
  subroutine remap_seqx1x2_to_parx1x2_field2d_VP4D( field2d )
    
    type(field2d_VP4D_t), intent(inout) :: field2d

    sll_int32 :: i1, i2
    sll_int32 :: iloc1, iloc2
    sll_int32 :: loc4d_sz_x1, loc4d_sz_x2
    sll_int32, dimension(1:2) :: glob_ind2d

    call sll_o_compute_local_sizes( field2d%layout2d_parx1x2, &
      loc4d_sz_x1, &
      loc4d_sz_x2 )
    do iloc2 = 1,loc4d_sz_x2
      do iloc1 = 1,loc4d_sz_x1
        glob_ind2d(:) = sll_o_local_to_global( &
            field2d%layout2d_parx1x2, &
            (/iloc1,iloc2/) )
        i1 = glob_ind2d(1)
        i2 = glob_ind2d(2)
      end do
    end do
    
  end subroutine remap_seqx1x2_to_parx1x2_field2d_VP4D
  !---------------------------------------------------------------------------


  !===========================================================================
  !> Field2d: Printing of different cross-sections in HDF5 format
  !---------------------------------------------------------------------------
  subroutine print_field2d_VP4D( &
      field2d, &
      idiag_num )

    use sll_m_collective
    use sll_m_hdf5_io_serial, only: sll_s_hdf5_ser_file_create, &
         sll_o_hdf5_ser_write_array, &
         sll_s_hdf5_ser_file_close, &
         sll_t_hdf5_ser_handle

    type(field2d_VP4D_t), intent(in) :: field2d
    sll_int32 , optional, intent(in) :: idiag_num 

    !-> Local variables
    !--> For parallelisation
    sll_int32 :: my_rank
    !--> For  HDF5 saving
    integer   :: file_err
    type(sll_t_hdf5_ser_handle) :: handle    !< file handle
    character(len=80)    :: filename_HDF5
    character(20) , save :: numfmt = "'_d',i5.5"
    character(len=30)    :: var_name

    !*** Construction of the name for the output HDF5 file ***
    if ( present(idiag_num) ) then
      write(filename_HDF5,'(A,'//numfmt//',A)') &
          trim(field2d%name), idiag_num, ".h5"
    else
      write(filename_HDF5,'(A,A)') &
          trim(field2d%name), ".h5"      
    end if

    !*** HDF5 saving ***
    my_rank = sll_f_get_collective_rank(sll_v_world_collective)

    !--> Saving of the cross-section in (x1,x2)
    if ( my_rank==0 ) then
      print*,'===> File ',trim(filename_HDF5),' saved by'
      call sll_s_hdf5_ser_file_create( trim(filename_HDF5), handle, file_err )
      !--> Saving of val2d_seqx1x2 of the 2d field
      write(var_name,'(A,A)') trim(field2d%name), "_x1x2"
      call sll_o_hdf5_ser_write_array( &
          handle, &
          field2d%val2d_seqx1x2(:,:), &
          trim(var_name), &
          file_err )
      call sll_s_hdf5_ser_file_close( handle, file_err )
    end if

  end subroutine print_field2d_VP4D
  !---------------------------------------------------------------------------

end module field2d_VP4D_module
!---------------------------------------------------------------------------
