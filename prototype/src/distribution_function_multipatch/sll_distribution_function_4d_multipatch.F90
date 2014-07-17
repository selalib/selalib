!------------------------------------------------------------------------------
! SELALIB
!------------------------------------------------------------------------------
!
! MODULE: sll_distribution_function_4d_multipatch
!
!> @author
!> - Edwin
!
! DESCRIPTION: 
!
!> @brief
!> Encapsulates a group of distribution functions, each associated with a 
!> patch.
!>
!>@details
!>
!
! REVISION HISTORY:
! 08 jul 2014 - Initial Version
! TODO_dd_mmm_yyyy - TODO_describe_appropriate_changes - TODO_name
!------------------------------------------------------------------------------
module sll_distribution_function_4d_multipatch_module
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_file_io.h"
  use sll_coordinate_transformation_multipatch_module
  use sll_constants
  use sll_remapper
  use sll_collective
  use sll_module_interpolators_2d_base
  use sll_arbitrary_degree_spline_interpolator_2d_module
  use sll_utilities
  use sll_boundary_condition_descriptors
  use sll_gnuplot
  use sll_parallel_array_initializer_module
  use sll_module_scalar_field_2d_multipatch
  implicit none


  type ::  sll_distribution_function_4d_multipatch
     type(sll_collective_t), pointer :: collective
     sll_int32 :: num_patches
     sll_int32 :: nproc_factor1
     sll_int32 :: nproc_factor2
     logical   :: ready_for_sequential_ops_in_x1x2 = .false.
     type(sll_logical_mesh_2d), pointer :: mesh_v ! same for all patches
     type(sll_coordinate_transformation_multipatch_2d), pointer :: transf
     type(layout_4d_ptr), dimension(:), pointer :: layouts_x1x2
     type(layout_4d_ptr), dimension(:), pointer :: layouts_x3x4
     type(layout_2d_ptr), dimension(:), pointer :: layouts_split
     type(layout_2d_ptr), dimension(:), pointer :: layouts_full
     type(data_4d_ptr), dimension(:), pointer :: f_x1x2
     type(data_4d_ptr), dimension(:), pointer :: f_x3x4
     type(multipatch_data_2d_real), dimension(:), pointer :: rho_split
     type(remap_plan_2d_real64_ptr), dimension(:), pointer :: remap_split2full
     type(remap_plan_4d_real64_ptr), dimension(:), pointer :: remap_x1x2tox3x4
     type(remap_plan_4d_real64_ptr), dimension(:), pointer :: remap_x3x4tox1x2

   contains
     procedure, pass(df) :: allocate_memory => allocate_memory_df_4d_mp
     procedure, pass(df) :: initialize => initialize_df_4d_mp 
     procedure, pass(df) :: set_to_sequential_x1x2 => x3x4_to_x1x2
     procedure, pass(df) :: set_to_sequential_x3x4 => x1x2_to_x3x4
     procedure, pass(df) :: delete => delete_df_4d_mp
  end type sll_distribution_function_4d_multipatch

  type :: data_4d_ptr
     sll_real64, dimension(:,:,:,:), pointer :: f
  end type data_4d_ptr

  interface sll_delete
     module procedure delete_df_4d_mp_ptr
  end interface sll_delete

contains

  ! Note how all the dimensions  in the velocity space are the same for all
  ! patches.
  ! nproc_factor1*nproc_factor2 = N, where N is the total number of processors
  ! in the collective. The distribution function multipatch object will
  ! use this information to create the different layouts it needs.
  function sll_new_distribution_function_4d_multipatch( &
       collective, &
       transf_mp, &
       mesh_v, & ! same for all patches
       nproc_factor1, &
       nproc_factor2 ) result(df)

    type(sll_distribution_function_4d_multipatch), pointer :: df
    type(sll_collective_t), pointer :: collective
    type(sll_coordinate_transformation_multipatch_2d), intent(in), target:: &
         transf_mp
    type(sll_logical_mesh_2d), pointer :: mesh_v
    sll_int32, intent(in) :: nproc_factor1
    sll_int32, intent(in) :: nproc_factor2
    sll_int32 :: ierr
    sll_int32 :: num_proc_total
    sll_int32 :: myrank

    num_proc_total = sll_get_collective_size(collective)
    myrank = sll_get_collective_rank(collective)

    SLL_ALLOCATE(df,ierr)
    df%collective => collective
    df%mesh_v => mesh_v
    df%transf => transf_mp
    df%nproc_factor1 = nproc_factor1
    df%nproc_factor2 = nproc_factor2
    df%num_patches = transf_mp%get_number_patches()
    SLL_ALLOCATE( df%f_x1x2(df%num_patches), ierr )
    SLL_ALLOCATE( df%f_x3x4(df%num_patches), ierr )
    SLL_ALLOCATE( df%layouts_x1x2(df%num_patches), ierr )
    SLL_ALLOCATE( df%layouts_x3x4(df%num_patches), ierr )
    SLL_ALLOCATE( df%layouts_split(df%num_patches), ierr )
    SLL_ALLOCATE( df%layouts_full(df%num_patches), ierr )
    SLL_ALLOCATE( df%rho_split(df%num_patches), ierr )
    SLL_ALLOCATE( df%remap_split2full(df%num_patches), ierr )
    SLL_ALLOCATE( df%remap_x1x2tox3x4(df%num_patches), ierr )
    SLL_ALLOCATE( df%remap_x3x4tox1x2(df%num_patches), ierr )

    call df%allocate_memory()

  end function sll_new_distribution_function_4d_multipatch

  subroutine allocate_memory_df_4d_mp( df )
    class(sll_distribution_function_4d_multipatch), intent(inout) :: df
    sll_int32 :: ierr
    sll_int32 :: i
    sll_int32 :: j
    sll_int32 :: ip
    sll_int32 :: np_x1
    sll_int32 :: np_x2
    sll_int32 :: np_x3
    sll_int32 :: np_x4
    sll_int32 :: locsz1
    sll_int32 :: locsz2
    sll_int32 :: locsz3
    sll_int32 :: locsz4
    sll_int32 :: myrank
    sll_int32 :: col_size
    sll_int32 :: imax
    sll_int32 :: jmax
    type(sll_logical_mesh_2d), pointer :: lm

    np_x3 = df%mesh_v%num_cells1+1
    np_x4 = df%mesh_v%num_cells2+1
    myrank = sll_get_collective_rank(df%collective)
    col_size = sll_get_collective_size(df%collective)

    do ip=0,df%num_patches-1
       ! obtain global dimensions for each array associated with each patch.
       np_x1 = df%transf%get_num_cells_eta1(ip) + 1
       np_x2 = df%transf%get_num_cells_eta2(ip) + 1
       df%layouts_x1x2(ip+1)%l => new_layout_4D(df%collective)
       call initialize_layout_with_distributed_4D_array( &
            np_x1, &
            np_x2, &
            np_x3, &
            np_x4, &
            1, &
            1, &
            df%nproc_factor1, &
            df%nproc_factor2, &
            df%layouts_x1x2(ip+1)%l )
       call compute_local_sizes_4d( &
            df%layouts_x1x2(ip+1)%l, &
            locsz1, &
            locsz2, &
            locsz3, &
            locsz4)
       SLL_ALLOCATE(df%f_x1x2(ip+1)%f(locsz1,locsz2,locsz3,locsz4),ierr)

       df%layouts_x3x4(ip+1)%l => new_layout_4D(df%collective)
       call initialize_layout_with_distributed_4D_array( &
            np_x1, &
            np_x2, &
            np_x3, &
            np_x4, &
            df%nproc_factor1, &
            df%nproc_factor2, &
            1, &
            1, &
            df%layouts_x3x4(ip+1)%l )
       call compute_local_sizes_4d( &
            df%layouts_x3x4(ip+1)%l, &
            locsz1, &
            locsz2, &
            locsz3, &
            locsz4)
       SLL_ALLOCATE(df%f_x3x4(ip+1)%f(locsz1,locsz2,locsz3,locsz4),ierr)
       SLL_ALLOCATE(df%rho_split(ip+1)%array(locsz1,locsz2),ierr)

       df%layouts_split(ip+1)%l => &
            new_layout_2D_from_layout_4D( df%layouts_x3x4(ip+1)%l )

       df%layouts_full(ip+1)%l => new_layout_2d(df%collective)
       ! The layouts_full are used to execute an all to all operation. The rho
       ! data which after the reduction are split in the x1 and x2 dimensions
       ! are to be collected into each processor redundantly. Solving for the 
       ! potential will therefore be redundantly done. For lack of something
       ! better, we fill the layout manually.
       imax = df%transf%get_num_cells_eta1(ip) + 1
       jmax = df%transf%get_num_cells_eta2(ip) + 1

       do j=0,col_size-1
          call set_layout_i_min( df%layouts_full(ip+1)%l, j, 1)
          call set_layout_j_min( df%layouts_full(ip+1)%l, j, 1)
          call set_layout_i_max( df%layouts_full(ip+1)%l, j, imax)
          call set_layout_j_max( df%layouts_full(ip+1)%l, j, jmax)
       end do

       ! call sll_view_lims_2d(df%layouts_split(ip+1)%l)
       ! call sll_view_lims_2d(df%layouts_full(ip+1)%l)

       df%remap_split2full(ip+1)%r => &
            new_remap_plan( df%layouts_split(ip+1)%l, &
                            df%layouts_full(ip+1)%l, &
                            df%rho_split(ip+1)%array )

       df%remap_x1x2tox3x4(ip+1)%r => &
            new_remap_plan( df%layouts_x1x2(ip+1)%l, &
                            df%layouts_x3x4(ip+1)%l, &
                            df%f_x1x2(ip+1)%f )

       df%remap_x3x4tox1x2(ip+1)%r => &
            new_remap_plan( df%layouts_x3x4(ip+1)%l, &
                            df%layouts_x1x2(ip+1)%l, &
                            df%f_x3x4(ip+1)%f )

    end do
  end subroutine allocate_memory_df_4d_mp

  subroutine initialize_df_4d_mp( &
       df, &
       init_func, &
       init_func_params )
    class(sll_distribution_function_4d_multipatch), intent(inout) :: df
    procedure(sll_scalar_initializer_4d)        :: init_func
    sll_real64, dimension(:), optional          :: init_func_params
    sll_int32 :: num_patches
    sll_int32 :: i
    type(layout_4d), pointer :: layout
    type(sll_logical_mesh_2d), pointer :: lm
    class(sll_coordinate_transformation_2d_base), pointer :: t
    num_patches = df%num_patches

    do i=0,num_patches-1
       layout => df%layouts_x3x4(i+1)%l
       ! please correct this to:
       ! lm => df%transf%get_logical_mesh(i)
       ! whenever gfortan 4.6 is no longer supported by Selalib.
       lm => df%transf%transfs(i+1)%t%mesh
       t => df%transf%transfs(i+1)%t
       call sll_4d_parallel_array_initializer( &
            layout, &
            lm, &
            df%mesh_v, &
            df%f_x3x4(i+1)%f, &
            init_func, &
            init_func_params, &
            t )
    end do
  end subroutine initialize_df_4d_mp

  ! Note that to carry out multiple remap operations is expensive. The
  ! latency would be multiplied by the number of patches... Other means
  ! of parallelization are more interesting, like setting each patch in its
  ! own process, if the advections can be properly made to work.
  subroutine x3x4_to_x1x2( df )
    class(sll_distribution_function_4d_multipatch), intent(inout) :: df
    sll_int32 :: num_patches
    sll_int32 :: i

    if(df%ready_for_sequential_ops_in_x1x2 .eqv. .true.) then
       print *, 'ERROR, x3x4_to_x1x2(): the distribution function multipatch ',&
            'is already configured for x1x2 operations.'
    end if

    num_patches = df%num_patches
    do i=0, num_patches-1
       call apply_remap_4D_double( &
            df%remap_x3x4tox1x2(i+1)%r, &
            df%f_x3x4(i+1)%f, &
            df%f_x1x2(i+1)%f )
    end do

    df%ready_for_sequential_ops_in_x1x2 = .true.

  end subroutine x3x4_to_x1x2

  subroutine x1x2_to_x3x4( df )
    class(sll_distribution_function_4d_multipatch), intent(inout) :: df
    sll_int32 :: num_patches
    sll_int32 :: i

    if(df%ready_for_sequential_ops_in_x1x2 .eqv. .false.) then
       print *, 'ERROR, x1x2_to_x3x4(): the distribution function multipatch ',&
            'is already configured for x3x4 operations.'
    end if

    num_patches = df%num_patches
    do i=0, num_patches-1
       call apply_remap_4D_double( &
            df%remap_x1x2tox3x4(i+1)%r, &
            df%f_x1x2(i+1)%f, &
            df%f_x3x4(i+1)%f )
    end do

    df%ready_for_sequential_ops_in_x1x2 = .false.

  end subroutine x1x2_to_x3x4


  ! This assumes that df is configured for sequential operations in x3 and x4.
  subroutine compute_charge_density_multipatch( df, rho )
    type(sll_distribution_function_4d_multipatch), intent(in) :: df
    class(sll_scalar_field_multipatch_2d), intent(inout)      :: rho
    sll_real64 :: delta3
    sll_real64 :: delta4
    sll_int32  :: ipatch
    sll_int32  :: i
    sll_int32  :: j
    sll_int32  :: num_patches
    sll_int32  :: locsz1
    sll_int32  :: locsz2
    sll_real64, dimension(:,:), pointer :: remap_out

    num_patches =  df%num_patches
    delta3  = df%mesh_v%delta_eta1
    delta4  = df%mesh_v%delta_eta2

    do ipatch=0, num_patches - 1
       call compute_local_sizes( df%layouts_split(ipatch+1)%l, locsz1, locsz2 )
       do j=1, locsz2
          do i=1, locsz1
             df%rho_split(ipatch+1)%array(i,j) = &
                  delta3*delta4*sum(df%f_x3x4(ipatch+1)%f(i,j,:,:))
          end do
       end do

       ! Reconfigure the data to store redundantly all the values of rho.
       remap_out => rho%get_patch_data_pointer(ipatch)
       call apply_remap_2D_double( &
            df%remap_split2full(ipatch+1)%r, &
            df%rho_split(ipatch+1)%array, &
            remap_out )
    end do

  end subroutine compute_charge_density_multipatch


  subroutine delete_df_4d_mp( df )
    class(sll_distribution_function_4d_multipatch), intent(inout) :: df
    sll_int32 :: i
    sll_int32 :: num_patches
    sll_int32 :: ierr

    num_patches = df%num_patches
    do i = 0,num_patches-1
       SLL_DEALLOCATE( df%f_x1x2(i+1)%f, ierr )
       SLL_DEALLOCATE( df%f_x3x4(i+1)%f, ierr )
       SLL_DEALLOCATE( df%rho_split(i+1)%array, ierr )
       call sll_delete( df%layouts_x1x2(i+1)%l )
       call sll_delete( df%layouts_x3x4(i+1)%l )
       call sll_delete( df%layouts_split(i+1)%l )
       call sll_delete( df%layouts_full(i+1)%l )
       call sll_delete( df%remap_split2full(i+1)%r )
       call sll_delete( df%remap_x1x2tox3x4(i+1)%r )
       call sll_delete( df%remap_x3x4tox1x2(i+1)%r )
    end do
    
    SLL_DEALLOCATE( df%f_x1x2, ierr )
    SLL_DEALLOCATE( df%f_x3x4, ierr )
    SLL_DEALLOCATE( df%layouts_x1x2, ierr )
    SLL_DEALLOCATE( df%layouts_x3x4, ierr )
    SLL_DEALLOCATE( df%layouts_split, ierr )
    SLL_DEALLOCATE( df%layouts_full, ierr )
    SLL_DEALLOCATE( df%rho_split, ierr )
    SLL_DEALLOCATE( df%remap_split2full, ierr )
    SLL_DEALLOCATE( df%remap_x1x2tox3x4, ierr )
    SLL_DEALLOCATE( df%remap_x3x4tox1x2, ierr )

  end subroutine delete_df_4d_mp

  subroutine delete_df_4d_mp_ptr( df )
    type(sll_distribution_function_4d_multipatch), pointer :: df
    sll_int32 :: ierr
    call df%delete()
    SLL_DEALLOCATE(df, ierr)
  end subroutine delete_df_4d_mp_ptr

end module sll_distribution_function_4d_multipatch_module
