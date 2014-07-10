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
  implicit none


  type ::  sll_distribution_function_4d_multipatch
     type(sll_collective_t), pointer :: collective
     sll_int32 :: num_patches
     sll_int32 :: nproc_factor1
     sll_int32 :: nproc_factor2
     type(sll_logical_mesh_2d), pointer :: mesh_v ! same for all patches
     type(sll_coordinate_transformation_multipatch_2d), pointer :: transf
     type(layout_4d_ptr), dimension(:), pointer :: layouts_x1x2
     type(layout_4d_ptr), dimension(:), pointer :: layouts_x3x4
     type(data_4d_ptr), dimension(:), pointer :: f_x1x2
     type(data_4d_ptr), dimension(:), pointer :: f_x3x4
   contains
     procedure, pass(df) :: allocate_memory => allocate_memory_df_4d_mp
     procedure, pass(df) :: initialize => initialize_df_4d_mp 
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

    call df%allocate_memory()

  end function sll_new_distribution_function_4d_multipatch

  subroutine allocate_memory_df_4d_mp( df )
    class(sll_distribution_function_4d_multipatch), intent(inout) :: df
    sll_int32 :: ierr
    sll_int32 :: i
    sll_int32 :: np_x1
    sll_int32 :: np_x2
    sll_int32 :: np_x3
    sll_int32 :: np_x4
    sll_int32 :: locsz1
    sll_int32 :: locsz2
    sll_int32 :: locsz3
    sll_int32 :: locsz4
    type(sll_logical_mesh_2d), pointer :: lm

    np_x3 = df%mesh_v%num_cells1+1
    np_x4 = df%mesh_v%num_cells2+1

    do i=0,df%num_patches-1
       ! obtain global dimensions for each array associated with each patch.
       np_x1 = df%transf%get_num_cells_eta1(i) + 1
       np_x2 = df%transf%get_num_cells_eta2(i) + 1
       df%layouts_x1x2(i+1)%l => new_layout_4D(df%collective)
       call initialize_layout_with_distributed_4D_array( &
            np_x1, &
            np_x2, &
            np_x3, &
            np_x4, &
            1, &
            1, &
            df%nproc_factor1, &
            df%nproc_factor2, &
            df%layouts_x1x2(i+1)%l )
       call compute_local_sizes_4d( &
            df%layouts_x1x2(i+1)%l, &
            locsz1, &
            locsz2, &
            locsz3, &
            locsz4)
       SLL_ALLOCATE(df%f_x1x2(i+1)%f(locsz1,locsz2,locsz3,locsz4),ierr)

       df%layouts_x3x4(i+1)%l => new_layout_4D(df%collective)
       call initialize_layout_with_distributed_4D_array( &
            np_x1, &
            np_x2, &
            np_x3, &
            np_x4, &
            df%nproc_factor1, &
            df%nproc_factor2, &
            1, &
            1, &
            df%layouts_x3x4(i+1)%l )
       call compute_local_sizes_4d( &
            df%layouts_x3x4(i+1)%l, &
            locsz1, &
            locsz2, &
            locsz3, &
            locsz4)
       SLL_ALLOCATE(df%f_x3x4(i+1)%f(locsz1,locsz2,locsz3,locsz4),ierr)


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
!       call sll_4d_parallel_array_initializer( &
call sll_2d_times_2d_parallel_array_initializer(&
            layout, &
            lm, &
            df%mesh_v, &
            df%f_x3x4(i+1)%f, &
            init_func, &
            init_func_params, &
            t )
    end do
  end subroutine initialize_df_4d_mp


  subroutine delete_df_4d_mp( df )
    class(sll_distribution_function_4d_multipatch), intent(inout) :: df
    sll_int32 :: i
    sll_int32 :: num_patches
    sll_int32 :: ierr

    num_patches = df%num_patches
    do i = 0,num_patches-1
       SLL_DEALLOCATE(df%f_x1x2(i+1)%f, ierr)
       SLL_DEALLOCATE(df%f_x3x4(i+1)%f, ierr)
    end do
    
    SLL_DEALLOCATE(df%f_x1x2, ierr)
    SLL_DEALLOCATE(df%f_x3x4, ierr)
    SLL_DEALLOCATE(df%layouts_x1x2, ierr)
    SLL_DEALLOCATE(df%layouts_x3x4, ierr)
    
  end subroutine delete_df_4d_mp

  subroutine delete_df_4d_mp_ptr( df )
    type(sll_distribution_function_4d_multipatch), pointer :: df
    sll_int32 :: ierr
    call df%delete()
    SLL_DEALLOCATE(df, ierr)
  end subroutine delete_df_4d_mp_ptr

end module sll_distribution_function_4d_multipatch_module
