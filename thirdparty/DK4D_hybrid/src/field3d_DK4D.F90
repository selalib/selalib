!===========================================================================
!> 3D scalar field in (x,y,eta3) for
!>  4D drift-kinetic hybrid simulation
!>
!> \date 2014-08-20
!> \author V. Grandgirard
!---------------------------------------------------------------------------
module field3d_DK4D_module
#include "sll_working_precision.h"
#include "sll_memory.h"

  use input_DK4D_module
  use mesh_DK4D_module
  use sll_m_collective, only : sll_v_world_collective, &
      sll_f_get_collective_size, sll_f_get_collective_rank
  use sll_m_remapper
  use sll_m_utilities, only : sll_f_is_even
  use utils_DK4D_module

  implicit none

  !===========================================================================
  !> 
  !---------------------------------------------------------------------------
  type, public :: field3d_DK4D_t

    !> Name of the field
    character(len=30) :: name

    !> Boundary conditions
    type(boundary_conditions_3d_t):: bound_cond

    !> Parallel in x3 and sequential in (x1,x2)
    type(sll_t_layout_3d), pointer :: layout3d_seqx1x2
    sll_real64, dimension(:,:,:), pointer :: val3d_seqx1x2 

    !> Parallel in (x1,x2) and sequential in x3
    type(sll_t_layout_3d), pointer :: layout3d_seqx3
    sll_real64, dimension(:,:,:), pointer :: val3d_seqx3 

    !> For remapping
    type(sll_t_remap_plan_3d_real64), pointer :: seqx1x2_to_seqx3
    type(sll_t_remap_plan_3d_real64), pointer :: seqx3_to_seqx1x2

    !> For interpolations
    type(spline_degree_3d_t) :: spline_degree

  end type field3d_DK4D_t
  !---------------------------------------------------------------------------

contains

  !===========================================================================
  !> Field3D: Allocation 
  !---------------------------------------------------------------------------
  subroutine new_field3d_DK4D( &
      field3d, &
      field3d_name, &
      mesh4d, &
      bound_cond, &
      spline_degree )

    type(field3d_DK4D_t)          , intent(inout) :: field3d
    character(len=*)              , intent(in)    :: field3d_name
    type(mesh_DK4D_t)             , intent(in)    :: mesh4d
    type(boundary_conditions_3d_t), intent(in)    :: bound_cond
    type(spline_degree_3d_t)      , intent(in)    :: spline_degree

    !-> Local variables 
    sll_int32 :: ierr
    !--> For MPI parallelization
    sll_int32 :: world_size
    sll_int32 :: my_rank
    sll_int32 :: power2
    sll_int32 :: nproc_x1, nproc_x2, nproc_x3
    !--> For parallel distribution
    sll_int32 :: npoints_x1, npoints_x2, npoints_x3
    sll_int32 :: loc3d_sz_x1, loc3d_sz_x2, loc3d_sz_x3

    !**********************************************************
    !*** Initialization of the name of the field            ***
    !**********************************************************
    field3d%name = trim(field3d_name)

    !**********************************************************
    !*** Initialization of the parallel distribution        ***
    !**********************************************************
    world_size = sll_f_get_collective_size(sll_v_world_collective)
    my_rank    = sll_f_get_collective_rank(sll_v_world_collective)
    power2 = int(log(real(world_size))/log(2.0))

    !**********************************************************************
    !*** Initialization of parallel layout of field3d in (x1,x2)        ***
    !***  directions                                                    *** 
    !***  (x1,x2) : parallelized layout                                 ***
    !***  x3      : sequential                                          ***
    !**********************************************************************
    !--> special case N = 1, so power2 = 0
    if ( power2 == 0 ) then
       nproc_x1 = 1
       nproc_x2 = 1
       nproc_x3 = 1
    end if
    
    if ( sll_f_is_even(power2) ) then
       nproc_x1 = 2**(power2/2)
       nproc_x2 = 2**(power2/2)
       nproc_x3 = 1
    else 
       nproc_x1 = 2**((power2-1)/2)
       nproc_x2 = 2**((power2+1)/2)
       nproc_x3 = 1
    end if

    field3d%layout3d_seqx3  => sll_f_new_layout_3d( sll_v_world_collective )
    npoints_x1 = mesh4d%eta1_eta2_mesh2d%num_cells1 + 1
    npoints_x2 = mesh4d%eta1_eta2_mesh2d%num_cells2 + 1
    npoints_x3 = mesh4d%eta3_mesh1d%num_cells + 1

    call sll_o_initialize_layout_with_distributed_array( &
      npoints_x1, &
      npoints_x2, &
      npoints_x3, &
      nproc_x1, &
      nproc_x2, &
      nproc_x3, &
      field3d%layout3d_seqx3 )
    call sll_o_compute_local_sizes( field3d%layout3d_seqx3, &
      loc3d_sz_x1, &
      loc3d_sz_x2, &
      loc3d_sz_x3)

    SLL_ALLOCATE(field3d%val3d_seqx3(loc3d_sz_x1,loc3d_sz_x2,loc3d_sz_x3), ierr )

    !**********************************************************************
    !*** Initialization of parallel layout of field3d in x3             ***
    !***  directions                                                    *** 
    !***  (x1,x2) : sequential                                          *** 
    !***  x3      : parallelized layout                                 ***
    !**********************************************************************
    nproc_x3 = nproc_x1*nproc_x2
    nproc_x1 = 1
    nproc_x2 = 1

    field3d%layout3d_seqx1x2  => sll_f_new_layout_3d( sll_v_world_collective )
    call sll_o_initialize_layout_with_distributed_array( &
      npoints_x1, & 
      npoints_x2, & 
      npoints_x3, &
      nproc_x1, &
      nproc_x2, &
      nproc_x3, &
      field3d%layout3d_seqx1x2 )
    call sll_o_compute_local_sizes( field3d%layout3d_seqx1x2, &
      loc3d_sz_x1, &
      loc3d_sz_x2, &
      loc3d_sz_x3)
    SLL_ALLOCATE( field3d%val3d_seqx1x2(loc3d_sz_x1,loc3d_sz_x2,loc3d_sz_x3), ierr )

    !****************************************
    !*** Initialization of the remaping   ***
    !****************************************
    field3d%seqx3_to_seqx1x2 => &
      sll_o_new_remap_plan( field3d%layout3d_seqx3, &
      field3d%layout3d_seqx1x2, field3d%val3d_seqx3 )
    
    field3d%seqx1x2_to_seqx3 => &
      sll_o_new_remap_plan( field3d%layout3d_seqx1x2, &
      field3d%layout3d_seqx3,field3d%val3d_seqx1x2 )

    !*************************************************
    !*** Initialization of the boundary conditions ***
    !*************************************************
    field3d%bound_cond%left_eta1  = bound_cond%left_eta1
    field3d%bound_cond%right_eta1 = bound_cond%right_eta1
    field3d%bound_cond%left_eta2  = bound_cond%left_eta2
    field3d%bound_cond%right_eta2 = bound_cond%right_eta2
    field3d%bound_cond%left_eta3  = bound_cond%left_eta3
    field3d%bound_cond%right_eta3 = bound_cond%right_eta3

    !*************************************************
    !*** Initialization of the interpolators       ***
    !*************************************************
    !--> Initialize the spline degree
    field3d%spline_degree%eta1 = spline_degree%eta1
    field3d%spline_degree%eta2 = spline_degree%eta2
    field3d%spline_degree%eta3 = spline_degree%eta3

  end subroutine new_field3d_DK4D
  !---------------------------------------------------------------------------


  !===========================================================================
  !> Field3D: Deallocation 
  !---------------------------------------------------------------------------
  subroutine delete_field3d_DK4D( field3d )

    type(field3d_DK4D_t), intent(inout) :: field3d

    sll_int32 :: ierr

    call sll_o_delete( field3d%layout3d_seqx1x2 )
    call sll_o_delete( field3d%layout3d_seqx3 )
    SLL_DEALLOCATE( field3d%val3d_seqx1x2, ierr )
    SLL_DEALLOCATE( field3d%val3d_seqx3, ierr )

  end subroutine delete_field3d_DK4D
  !---------------------------------------------------------------------------


  !===========================================================================
  !> Field3D: Printing of different cross-sections in HDF5 format
  !---------------------------------------------------------------------------
  subroutine print_field3d_DK4D( &
      field3d, &
      ix1_diag, &
      ix2_diag, &
      ix3_diag, &
      idiag_num )

    use sll_m_collective
    use sll_m_hdf5_io_serial, only: sll_s_hdf5_ser_file_create, &
        sll_o_hdf5_ser_write_array, sll_o_hdf5_ser_write_array, &
        sll_s_hdf5_ser_file_close, sll_t_hdf5_ser_handle

    type(field3d_DK4D_t), intent(in) :: field3d
    sll_int32           , intent(in) :: ix1_diag
    sll_int32           , intent(in) :: ix2_diag
    sll_int32           , intent(in) :: ix3_diag
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
    !--> For cross-section saving
    sll_int32 :: ieta3_min_loc, ieta3_max_loc
    sll_int32, dimension(1:3) :: indx_loc
    sll_int32 :: ix3_diag_loc
    sll_real64, dimension(1:3) :: indx_diag

    !*** Construction of the name for the output HDF5 file ***
    if ( present(idiag_num) ) then
      write(filename_HDF5,'(A,'//numfmt//',A)') &
          trim(field3d%name), idiag_num, ".h5"
    else
      write(filename_HDF5,'(A,A)') &
          trim(field3d%name), ".h5"      
    end if

    !*** HDF5 saving ***
    my_rank = sll_f_get_collective_rank(sll_v_world_collective)

    !--> Saving of the cross-section in (x1,x2)
    ieta3_min_loc = sll_o_get_layout_k_min( field3d%layout3d_seqx1x2, my_rank )
    ieta3_max_loc = sll_o_get_layout_k_max( field3d%layout3d_seqx1x2, my_rank )

    if ( (ieta3_min_loc.le.ix3_diag) .and. (ix3_diag.le.ieta3_max_loc) ) then

      print*,'===> File ',trim(filename_HDF5),' saved by proc ',my_rank
      call sll_s_hdf5_ser_file_create( trim(filename_HDF5), handle, file_err )
      indx_loc = sll_o_global_to_local( field3d%layout3d_seqx1x2, &
          (/1,1,ix3_diag/) )
      ix3_diag_loc = indx_loc(3)
      !--> Saving of the array of indexes ix1_diag, ix2_diag and ix3_diag
      write(var_name,'(A,A)') trim(field3d%name), "_indx_diag"
      indx_diag(1) = ix1_diag
      indx_diag(2) = ix2_diag
      indx_diag(3) = ix3_diag
      call sll_o_hdf5_ser_write_array( &
          handle,  &
          indx_diag, &
          var_name, file_err )
      !--> Saving of val3d_seqx1x2 of the 3d field
      write(var_name,'(A,A)') trim(field3d%name), "_x1x2"
      call sll_o_hdf5_ser_write_array( &
          handle, &
          field3d%val3d_seqx1x2(:,:,ix3_diag_loc), &
          trim(var_name), &
          file_err )
      call sll_s_hdf5_ser_file_close( handle, file_err )

    end if

  end subroutine print_field3d_DK4D
  !---------------------------------------------------------------------------


end module field3d_DK4D_module
!---------------------------------------------------------------------------
