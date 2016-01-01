program parallel_advection

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

  use iso_fortran_env, only: &
    output_unit

  use sll_m_boundary_condition_descriptors, only: &
    sll_p_periodic

  use sll_m_collective, only: &
    sll_s_boot_collective, &
    sll_f_get_collective_rank, &
    sll_f_get_collective_size, &
    sll_s_halt_collective, &
    sll_v_world_collective

  use sll_m_cubic_spline_interpolator_1d, only: &
    sll_t_cubic_spline_interpolator_1d, &
    sll_o_delete

  use sll_m_gnuplot_parallel, only: &
    sll_o_gnuplot_2d_parallel

  use sll_m_interpolators_1d_base, only: &
    sll_c_interpolator_1d

  use sll_m_remapper, only: &
    sll_o_apply_remap_2d, &
    sll_o_compute_local_sizes, &
    sll_o_initialize_layout_with_distributed_array, &
    sll_t_layout_2d, &
    sll_o_local_to_global, &
    sll_f_new_layout_2d, &
    sll_o_new_remap_plan, &
    sll_t_remap_plan_2d_real64, &
    sll_o_view_lims, &
    sll_o_delete

  use sll_m_utilities, only: &
    sll_f_is_power_of_two

  use sll_mpi, only: &
    mpi_wtime

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#define MPI_MASTER 0

  class(sll_c_interpolator_1d), pointer   :: interp_eta1
  class(sll_c_interpolator_1d), pointer   :: interp_eta2
  type(sll_t_cubic_spline_interpolator_1d), target :: spl_eta1
  type(sll_t_cubic_spline_interpolator_1d), target :: spl_eta2
  sll_real64, dimension(:,:),  pointer       :: f_eta1
  sll_real64, dimension(:,:),  pointer       :: f_eta2
  type(sll_t_layout_2d), pointer                   :: layout_eta1
  type(sll_t_layout_2d), pointer                   :: layout_eta2
  type(sll_t_remap_plan_2d_real64), pointer        :: eta1_to_eta2 
  type(sll_t_remap_plan_2d_real64), pointer        :: eta2_to_eta1

  sll_int32  :: i_step
  sll_int32  :: prank, comm
  sll_int32  :: loc_sz_i, loc_sz_j
  sll_int64  :: psize
  sll_real64 :: tcpu1
  sll_real64 :: tcpu2
  sll_int32  :: gi
  sll_int32  :: gj
  sll_int32  :: global_indices(2)
  sll_int32  :: error
  sll_real64 :: eta1
  sll_real64 :: eta2
  sll_int32  :: offset(2)
  sll_real64 :: offset_eta1, offset_eta2

!###########################################################
  !Simulation parameters and geometry sizes                !
                                                           !
  sll_int32,  parameter :: nc_eta1 = 64                    !
  sll_int32,  parameter :: nc_eta2 = 64                    !
  sll_real64, parameter :: eta1_min = - 5.0_f64            !
  sll_real64, parameter :: eta1_max =   5.0_f64            !
  sll_real64, parameter :: eta2_min = - 5.0_f64            !
  sll_real64, parameter :: eta2_max = + 5.0_f64            !
  sll_real64 :: delta_eta1 = (eta1_max-eta1_min)/nc_eta1   !
  sll_real64 :: delta_eta2 = (eta2_max-eta2_min)/nc_eta2   !
  sll_real64, parameter :: delta_t = 0.05_f64              !
  sll_int32,  parameter :: n_step  = 200                  !
                                                           !
!###########################################################

  sll_int32 :: i, j

  call sll_s_boot_collective()

  prank = sll_f_get_collective_rank(sll_v_world_collective)
  psize = int(sll_f_get_collective_size(sll_v_world_collective),kind=i64)
  comm  = sll_v_world_collective%comm

  tcpu1 = MPI_WTIME()
  if (.not. sll_f_is_power_of_two(psize)) then     
     print *, 'This test needs to run in a number of processes which is ',&
          'a power of 2.'
     stop
  end if

  call spl_eta1%initialize(nc_eta1+1,eta1_min,eta1_max,sll_p_periodic)

  call spl_eta2%initialize(nc_eta2+1,eta2_min,eta2_max,sll_p_periodic)

  interp_eta1 => spl_eta1
  interp_eta2 => spl_eta2

  layout_eta1 => sll_f_new_layout_2d( sll_v_world_collective )        

  call sll_o_initialize_layout_with_distributed_array( &
             nc_eta1+1, nc_eta2+1, 1,int(psize,4),layout_eta1)

  if ( prank == MPI_MASTER ) call sll_o_view_lims( layout_eta1 )
  flush( output_unit )

  call sll_o_compute_local_sizes(layout_eta1,loc_sz_i,loc_sz_j)        
  SLL_CLEAR_ALLOCATE(f_eta1(1:loc_sz_i,1:loc_sz_j),error)

  layout_eta2 => sll_f_new_layout_2d( sll_v_world_collective )

  call sll_o_initialize_layout_with_distributed_array( &
              nc_eta1+1, nc_eta2+1, int(psize,4),1,layout_eta2)

  if ( prank == MPI_MASTER ) call sll_o_view_lims( layout_eta2 )
  flush( output_unit )

  call sll_o_compute_local_sizes(layout_eta2,loc_sz_i,loc_sz_j)        
  SLL_CLEAR_ALLOCATE(f_eta2(1:loc_sz_i,1:loc_sz_j),error)

  eta1_to_eta2 => sll_o_new_remap_plan( layout_eta1, layout_eta2, f_eta1)     
  eta2_to_eta1 => sll_o_new_remap_plan( layout_eta2, layout_eta1, f_eta2)     
  
  do j=1,loc_sz_j
  do i=1,loc_sz_i

     global_indices = sll_o_local_to_global(layout_eta2,(/i,j/)) 
     gi = global_indices(1)
     gj = global_indices(2)

     eta1  = eta1_min+(gi-1)*delta_eta1
     eta2  = eta2_min+(gj-1)*delta_eta2

     f_eta2(i,j)=exp(-.5*(eta1*eta1+eta2*eta2))

  end do
  end do

  call sll_o_apply_remap_2d( eta2_to_eta1, f_eta2, f_eta1 )
  call advection_eta1(0.5*delta_t)

  offset = sll_o_local_to_global(layout_eta1,(/1,1/)) 
  offset_eta1 = (offset(1)-1)*delta_eta1
  offset_eta2 = (offset(2)-1)*delta_eta2

  do i_step=1, n_step

      call sll_o_apply_remap_2d( eta1_to_eta2, f_eta1, f_eta2 )

      call advection_eta2(delta_t)

      call sll_o_apply_remap_2d( eta2_to_eta1, f_eta2, f_eta1 )

      call advection_eta1(delta_t)

      call sll_o_gnuplot_2d_parallel( offset_eta1, delta_eta1, &
                                    offset_eta2, delta_eta2, &
                                    size(f_eta1,1), size(f_eta1,2), &
                                    f_eta1, 'f_parallel', &
                                    i_step, error)

  end do

  tcpu2 = MPI_WTIME()
  if (prank == MPI_MASTER) &
       write(*,"(//10x,' Wall time = ', G15.3, ' sec' )") (tcpu2-tcpu1)*psize

  call sll_o_delete(layout_eta1)
  call sll_o_delete(layout_eta2)
  SLL_DEALLOCATE_ARRAY(f_eta1, error)
  SLL_DEALLOCATE_ARRAY(f_eta2, error)

  call sll_s_halt_collective()

contains

 subroutine advection_eta1(dt)

  sll_real64, intent(in) :: dt
  sll_real64 :: alpha

  call sll_o_compute_local_sizes(layout_eta1,loc_sz_i,loc_sz_j)

  do j=1,loc_sz_j

     global_indices = sll_o_local_to_global(layout_eta1,(/1,j/)) 
     gj = global_indices(2)
     alpha = -dt
     call interp_eta1%interpolate_array_disp_inplace(loc_sz_i,f_eta1(:,j),alpha)

  end do

 end subroutine advection_eta1

 subroutine advection_eta2(dt)

  sll_real64, intent(in) :: dt
  sll_real64 :: alpha

  call sll_o_compute_local_sizes(layout_eta2,loc_sz_i,loc_sz_j)        

  do i=1,loc_sz_i

     global_indices = sll_o_local_to_global(layout_eta2,(/i,1/)) 
     gi = global_indices(1)
     alpha = -dt
     call interp_eta2%interpolate_array_disp_inplace(loc_sz_j,f_eta2(i,:),alpha)

  end do

 end subroutine advection_eta2

end program parallel_advection
