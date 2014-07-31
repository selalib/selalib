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
  use sll_timer
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
     procedure, pass(df) :: get_x1x2_data_slice_pointer => get_x1x2_slice_4d
     procedure, pass(df) :: get_x3_line_pointer => get_x3_line_df4d
     procedure, pass(df) :: get_x4_line_pointer => get_x4_line_df4d
     procedure, pass(df) :: get_full_patch_data_pointer => get_full_patch_df4d
     procedure, pass(df) :: set_to_sequential_x1x2 => x3x4_to_x1x2
     procedure, pass(df) :: set_to_sequential_x3x4 => x1x2_to_x3x4
     procedure, pass(df) :: get_local_data_sizes => get_locsz_df4d
     procedure, pass(df) :: get_eta_coordinates => get_eta_coords_df4d
     procedure, pass(df) :: compute_moment => moments_df4d
     procedure, pass(df) :: compute_Lp_norms => Lp_norms_df4d
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
    !type(sll_time_mark)  :: t0 ! delete this when done timing/debugging
    !sll_real64 :: time ! delete

    num_patches = df%num_patches

    do i=0,num_patches-1
       layout => df%layouts_x3x4(i+1)%l
       ! please correct this to:
       ! lm => df%transf%get_logical_mesh(i)
       ! whenever gfortan 4.6 is no longer supported by Selalib.
       lm => df%transf%transfs(i+1)%t%mesh
       t => df%transf%transfs(i+1)%t
       !call sll_set_time_mark(t0) ! delete this
       call sll_4d_parallel_array_initializer( &
            layout, &
            lm, &
            df%mesh_v, &
            df%f_x3x4(i+1)%f, &
            init_func, &
            init_func_params, &
            transf_x1_x2=t )
       !time =  sll_time_elapsed_since(t0)
       !print *, 'inside initialize_df_4d_mp: time to initialize array:', time
    end do
  end subroutine initialize_df_4d_mp


  function get_x1x2_slice_4d( df, patch, k, l ) result(ptr)
    sll_real64, dimension(:,:), pointer :: ptr
    class(sll_distribution_function_4d_multipatch), intent(in) :: df
    sll_int32, intent(in) :: patch
    sll_int32, intent(in) :: k  ! third index in 4D array
    sll_int32, intent(in) :: l  ! fourth index in 4D array

    SLL_ASSERT( (patch >= 0) .and. (patch <= df%num_patches - 1) )

    ptr => df%f_x1x2(patch+1)%f(:,:,k,l)
  end function get_x1x2_slice_4d

  function get_x3_line_df4d( df, patch, i, j, l ) result(ptr)
    sll_real64, dimension(:), pointer :: ptr
    class(sll_distribution_function_4d_multipatch), intent(in) :: df
    sll_int32, intent(in) :: patch
    sll_int32, intent(in) :: i  ! first index in 4D array
    sll_int32, intent(in) :: j  ! second index in 4D array
    sll_int32, intent(in) :: l  ! fourth index in 4D array

    SLL_ASSERT( (patch >= 0) .and. (patch <= df%num_patches - 1) )

    if(df%ready_for_sequential_ops_in_x1x2 .eqv. .true.) then
       print *, 'ERROR, get_x3_line_df4d(): data is configured for ', &
            'sequential operations in x1 and x2. Need to reconfigure the ', &
            'distribution function before making this call.'
    end if
    ptr => df%f_x3x4(patch+1)%f(i,j,:,l)
  end function get_x3_line_df4d

  function get_x4_line_df4d( df, patch, i, j, k ) result(ptr)
    sll_real64, dimension(:), pointer :: ptr
    class(sll_distribution_function_4d_multipatch), intent(in) :: df
    sll_int32, intent(in) :: patch
    sll_int32, intent(in) :: i  ! first index in 4D array
    sll_int32, intent(in) :: j  ! second index in 4D array
    sll_int32, intent(in) :: k  ! third index in 4D array

    SLL_ASSERT( (patch >= 0) .and. (patch <= df%num_patches - 1) )

    if(df%ready_for_sequential_ops_in_x1x2 .eqv. .true.) then
       print *, 'ERROR, get_x3_line_df4d(): data is configured for ', &
            'sequential operations in x1 and x2. Need to reconfigure the ', &
            'distribution function before making this call.'
    end if
    ptr => df%f_x3x4(patch+1)%f(i,j,k,:)
  end function get_x4_line_df4d

  function get_full_patch_df4d( df, patch ) result(ptr)
    sll_real64, dimension(:,:,:,:), pointer :: ptr
    class(sll_distribution_function_4d_multipatch), intent(in) :: df
    sll_int32, intent(in) :: patch

    SLL_ASSERT( (patch >= 0) .and. (patch <= df%num_patches - 1) )

    ptr => df%f_x3x4(patch+1)%f
  end function get_full_patch_df4d



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

  subroutine get_locsz_df4d( &
       df, &
       patch, &
       loc_sz_x1, &
       loc_sz_x2, &
       loc_sz_x3, &
       loc_sz_x4 )
 
    class(sll_distribution_function_4d_multipatch), intent(in) :: df
    sll_int32, intent(in)  :: patch
    sll_int32, intent(out) :: loc_sz_x1
    sll_int32, intent(out) :: loc_sz_x2
    sll_int32, intent(out) :: loc_sz_x3
    sll_int32, intent(out) :: loc_sz_x4

    SLL_ASSERT( (patch >=0) .and. (patch <= df%num_patches - 1) )

    if(df%ready_for_sequential_ops_in_x1x2 .eqv. .true.) then
       call compute_local_sizes_4d( &
            df%layouts_x1x2(patch+1)%l, &
            loc_sz_x1, &
            loc_sz_x2, &
            loc_sz_x3, &
            loc_sz_x4)
    end if

    if(df%ready_for_sequential_ops_in_x1x2 .eqv. .false.) then
       call compute_local_sizes_4d( &
            df%layouts_x3x4(patch+1)%l, &
            loc_sz_x1, &
            loc_sz_x2, &
            loc_sz_x3, &
            loc_sz_x4)
    end if
  end subroutine get_locsz_df4d

  ! The purpose of this function is to return the eta-coordinates of a given
  ! 4-tuple of local indices (i,j,k,l). Since the distribution function has
  ! multiple arrays corresponding to the different patches and furthermore it
  ! can be parallel, this work should be centralized.
  function get_eta_coords_df4d( df, patch, ijkl ) result(etas)
    sll_real64, dimension(4) :: etas
    class(sll_distribution_function_4d_multipatch), intent(in) :: df
    sll_int32,  intent(in)   :: patch
    sll_int32, dimension(4), intent(in) :: ijkl
    sll_int32, dimension(4)  :: gi  ! global indices
    sll_real64 :: eta1_min
    sll_real64 :: eta2_min
    sll_real64 :: eta3_min
    sll_real64 :: eta4_min
    sll_real64 :: delta1
    sll_real64 :: delta2
    sll_real64 :: delta3
    sll_real64 :: delta4

    SLL_ASSERT( (patch >=0) .and. (patch <= df%num_patches - 1) )

    eta1_min = df%transf%get_eta1_min( patch )
    eta2_min = df%transf%get_eta2_min( patch )
    eta3_min = df%mesh_v%eta1_min
    eta4_min = df%mesh_v%eta2_min

    delta1 = df%transf%get_delta_eta1( patch )
    delta2 = df%transf%get_delta_eta2( patch )
    delta3 = df%mesh_v%delta_eta1
    delta4 = df%mesh_v%delta_eta2

    if(df%ready_for_sequential_ops_in_x1x2 .eqv. .true.) then
       gi = local_to_global_4D(df%layouts_x1x2(patch+1)%l, ijkl)
       etas(1) = eta1_min + (gi(1)-1)*delta1
       etas(2) = eta2_min + (gi(2)-1)*delta2
       etas(3) = eta3_min + (gi(3)-1)*delta3
       etas(4) = eta4_min + (gi(4)-1)*delta4
    end if

    if(df%ready_for_sequential_ops_in_x1x2 .eqv. .false.) then
       gi = local_to_global_4D(df%layouts_x3x4(patch+1)%l, ijkl)
       etas(1) = eta1_min + (gi(1)-1)*delta1
       etas(2) = eta2_min + (gi(2)-1)*delta2
       etas(3) = eta3_min + (gi(3)-1)*delta3
       etas(4) = eta4_min + (gi(4)-1)*delta4
    end if
  end function get_eta_coords_df4d

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

  function moments_df4d( df, moment ) result(mom)
    sll_real64 :: mom
    class(sll_distribution_function_4d_multipatch), intent(inout) :: df
    sll_int32, intent(in) :: moment

    select case ( moment )
    case (0)
       mom = mass_df4d(df)
    case (1)
       mom = momentum_df4d(df)
    case (2)
       mom = kinetic_energy_df4d(df)
    end select
  end function moments_df4d

  function mass_df4d( df ) result(mass)
    sll_real64 :: mass
    class(sll_distribution_function_4d_multipatch), intent(inout) :: df
    sll_int32 :: num_patches
    sll_int32 :: ip
    sll_int32 :: i
    sll_int32 :: j
    sll_real64 :: accumulator
    sll_real64, dimension(:,:,:,:), pointer :: f
    type(sll_coordinate_transformation_multipatch_2d), pointer :: t
    type(sll_logical_mesh_2d), pointer :: mv
    sll_real64 :: delta1
    sll_real64 :: delta2
    sll_real64 :: delta3
    sll_real64 :: delta4
    sll_real64, dimension(1) :: send_b  ! send buffer
    sll_real64, dimension(1) :: recv_b  ! receive buffer
    sll_int32 :: loc_sz1
    sll_int32 :: loc_sz2
    sll_int32 :: loc_sz3
    sll_int32 :: loc_sz4
    sll_real64 :: rho_ij
    sll_real64, dimension(4) :: eta
    sll_real64 :: jacobian
    sll_int32 :: ierr

    num_patches =  df%num_patches
    accumulator =  0.0_f64
    t           => df%transf
    mv          => df%mesh_v

    delta3 = mv%delta_eta1
    delta4 = mv%delta_eta2

    do ip=0,num_patches - 1
       f => df%get_full_patch_data_pointer(ip)
       call df%get_local_data_sizes(ip, loc_sz1, loc_sz2, loc_sz3, loc_sz4)
       delta1 = t%get_delta_eta1(ip)
       delta2 = t%get_delta_eta2(ip)

       do j=1,loc_sz2
          do i=1, loc_sz1
             rho_ij      = sum(f(i,j,:,:))*delta3*delta4
             eta(:)      = df%get_eta_coordinates(ip, (/i,j,1,1/) )
             jacobian    = t%jacobian(eta(1),eta(2),ip)
             accumulator = accumulator + rho_ij*jacobian*delta1*delta2
          end do
       end do
    end do
    ! Finally, reduce the result so that all processes have the same answer.
    send_b(1) = accumulator
    call sll_collective_allreduce(df%collective, send_b, 1, MPI_SUM, recv_b)
    mass = recv_b(1)
  end function mass_df4d


  function momentum_df4d( df ) result(mom)
    sll_real64 :: mom
    class(sll_distribution_function_4d_multipatch), intent(inout) :: df
    sll_int32 :: num_patches
    sll_int32 :: ip
    sll_int32 :: i
    sll_int32 :: j
    sll_int32 :: k
    sll_int32 :: l
    sll_real64 :: accumulator
    sll_real64, dimension(:,:,:,:), pointer :: f
    type(sll_coordinate_transformation_multipatch_2d), pointer :: t
    type(sll_logical_mesh_2d), pointer :: mv
    sll_real64 :: delta1
    sll_real64 :: delta2
    sll_real64 :: delta3
    sll_real64 :: delta4
    sll_real64, dimension(1) :: send_b  ! send buffer
    sll_real64, dimension(1) :: recv_b  ! receive buffer
    sll_int32 :: loc_sz1
    sll_int32 :: loc_sz2
    sll_int32 :: loc_sz3
    sll_int32 :: loc_sz4
    sll_real64, dimension(4) :: eta
    sll_real64 :: jacobian
    sll_int32 :: ierr

    num_patches =  df%num_patches
    accumulator =  0.0_f64
    t           => df%transf
    mv          => df%mesh_v

    delta3 = mv%delta_eta1
    delta4 = mv%delta_eta2
    print *, 'momentum of the distribution function multipatch: ', &
         ' this function has not been implemented.'
    do ip=0,num_patches - 1
       f => df%get_full_patch_data_pointer(ip)
       call df%get_local_data_sizes(ip, loc_sz1, loc_sz2, loc_sz3, loc_sz4)
       delta1 = t%get_delta_eta1(ip)
       delta2 = t%get_delta_eta2(ip)
       do l=1, loc_sz4
          do k=1, loc_sz3
             do j=1, loc_sz2
                do i=1, loc_sz1
                   eta(:)      = df%get_eta_coordinates(ip, (/i,j,k,l/) )
                   accumulator = accumulator + sum(f(i,j,:,:))*delta3*delta4
                   jacobian    = t%jacobian(eta(1),eta(2),ip)
                !   accumulator = accumulator + rho_ij*jacobian*delta1*delta2
                end do
             end do
          end do
       end do
    end do
    send_b(1) = accumulator
    ! Finally, reduce the result so that all processes have the same answer.
    call sll_collective_allreduce(df%collective, send_b, 1, MPI_SUM, recv_b)
    mom = recv_b(1)
  end function momentum_df4d


  function kinetic_energy_df4d( df ) result(ke)
    sll_real64 :: ke
    class(sll_distribution_function_4d_multipatch), intent(inout) :: df
    sll_int32 :: num_patches
    sll_int32 :: ip
    sll_int32 :: i
    sll_int32 :: j
    sll_int32 :: k
    sll_int32 :: l
    sll_real64 :: accumulator
    sll_real64, dimension(:,:,:,:), pointer :: f
    type(sll_coordinate_transformation_multipatch_2d), pointer :: t
    type(sll_logical_mesh_2d), pointer :: mv
    sll_real64 :: delta1
    sll_real64 :: delta2
    sll_real64 :: delta3
    sll_real64 :: delta4
    sll_real64, dimension(1) :: send_b  ! send buffer
    sll_real64, dimension(1) :: recv_b  ! receive buffer
    sll_int32 :: loc_sz1
    sll_int32 :: loc_sz2
    sll_int32 :: loc_sz3
    sll_int32 :: loc_sz4
    sll_real64 :: rho_ij
    sll_real64, dimension(4) :: eta
    sll_real64 :: jacobian
    sll_real64 :: v2
    sll_real64 :: metric
    sll_int32 :: ierr

    num_patches =  df%num_patches
    accumulator =  0.0_f64
    t           => df%transf
    mv          => df%mesh_v

    delta3 = mv%delta_eta1
    delta4 = mv%delta_eta2

    do ip=0,num_patches - 1
       f => df%get_full_patch_data_pointer(ip)
       call df%get_local_data_sizes(ip, loc_sz1, loc_sz2, loc_sz3, loc_sz4)
       delta1 = t%get_delta_eta1(ip)
       delta2 = t%get_delta_eta2(ip)
       do l=1, loc_sz4
          do k=1, loc_sz3
             do j=1,loc_sz2
                do i=1, loc_sz1
                   eta(:)      = df%get_eta_coordinates(ip, (/i,j,k,l/) )
                   v2          = eta(3)*eta(3) + eta(4)*eta(4)
                   jacobian    = t%jacobian(eta(1),eta(2),ip)
                   metric      = jacobian*delta1*delta2*delta3*delta4
                   accumulator = accumulator + f(i,j,k,l)*v2*metric
                end do
             end do
          end do
       end do
    end do
    ! Finally, reduce the result so that all processes have the same answer.
    send_b(1) = accumulator
    call sll_collective_allreduce(df%collective, send_b, 1, MPI_SUM, recv_b)
    ke = recv_b(1)
  end function kinetic_energy_df4d

  function Lp_norms_df4d( df ) result(Lp_norms)
    sll_real64, dimension(3) :: Lp_norms
    class(sll_distribution_function_4d_multipatch), intent(inout) :: df    
    sll_int32  :: num_patches
    sll_int32  :: ip
    sll_int32  :: i
    sll_int32  :: j
    sll_int32  :: k
    sll_int32  :: l
    sll_real64 :: accumulator1
    sll_real64 :: accumulator2
    sll_real64 :: accumulator3
    sll_real64, dimension(:,:,:,:), pointer :: f
    type(sll_coordinate_transformation_multipatch_2d), pointer :: t
    type(sll_logical_mesh_2d), pointer :: mv
    sll_real64 :: delta1
    sll_real64 :: delta2
    sll_real64 :: delta3
    sll_real64 :: delta4
    sll_real64, dimension(3) :: send_b  ! send buffer
    sll_real64, dimension(3) :: recv_b  ! receive buffer
    sll_int32  :: loc_sz1
    sll_int32  :: loc_sz2
    sll_int32  :: loc_sz3
    sll_int32  :: loc_sz4
    sll_real64, dimension(4) :: eta
    sll_real64 :: jacobian
    sll_real64 :: metric
    sll_int32  :: ierr

    num_patches  =  df%num_patches
    accumulator1 =  0.0_f64
    accumulator2 =  0.0_f64
    accumulator3 =  0.0_f64
    t            => df%transf
    mv           => df%mesh_v

    delta3 = mv%delta_eta1
    delta4 = mv%delta_eta2

    do ip=0,num_patches - 1
       f => df%get_full_patch_data_pointer(ip)
       call df%get_local_data_sizes(ip, loc_sz1, loc_sz2, loc_sz3, loc_sz4)
       delta1 = t%get_delta_eta1(ip)
       delta2 = t%get_delta_eta2(ip)
       do l=1, loc_sz4
          do k=1, loc_sz3
             do j=1,loc_sz2
                do i=1, loc_sz1
                   eta(:)       = df%get_eta_coordinates(ip, (/i,j,k,l/) )
                   jacobian     = t%jacobian(eta(1),eta(2),ip)
                   metric       = jacobian*delta1*delta2*delta3*delta4
                   accumulator1 = accumulator1 + abs(f(i,j,k,l))*metric
                   accumulator2 = accumulator2  +     f(i,j,k,l)**2*metric
                   accumulator3 = max(abs(f(i,j,k,l)),accumulator3)
                end do
             end do
          end do
       end do
    end do
    ! Finally, reduce the result so that all processes have the same answer.
    send_b(1) = accumulator1
    send_b(2) = accumulator2
    send_b(3) = accumulator3
    call sll_collective_allreduce(df%collective, send_b, 3, MPI_SUM, recv_b)
    Lp_norms(:) = recv_b(:)
  end function Lp_norms_df4d



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
