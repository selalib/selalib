! Test the test_particle visualization interface module
! author: Martin Campos Pinto, CNRS

program test_particle_visualization_interface
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_particle_visualization_interface, only : &
    sll_t_plotting_params_2d, &
    sll_s_visualize_particle_group

  use sll_m_particle_group_base, only: &
    sll_c_particle_group_base

  use sll_m_particle_group_2d2v_lbf, only: &
    sll_s_new_particle_group_2d2v_lbf_ptr

  use sll_m_initial_distribution, only : &
     sll_c_distribution_params, &
     sll_s_initial_distribution_new

  use sll_m_particle_sampling_interface, only : &
    sll_s_sample_particle_group, &
    sll_s_resample_particle_group

  use sll_m_io_utilities, only : &
     sll_s_concatenate_filename_and_path


  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  class(sll_c_particle_group_base), pointer  :: particle_group

  class(sll_c_distribution_params), allocatable :: distribution_params
  character(len=256) :: distribution_str
  character(len=256) :: filename
  sll_int32  :: file_id

  type(sll_t_plotting_params_2d)          :: plotting_params_2d

  ! parameters of the lbf particle group
  logical :: domain_is_x_periodic
  logical :: domain_is_y_periodic
  sll_real64, dimension(4)  :: remapping_grid_eta_min
  sll_real64, dimension(4)  :: remapping_grid_eta_max
  sll_int32,  dimension(4)  :: remapping_sparse_grid_max_levels
  sll_int32   :: remap_degree
  sll_int32   :: n_particles_x
  sll_int32   :: n_particles_y
  sll_int32   :: n_particles_vx
  sll_int32   :: n_particles_vy

  sll_real64 :: charge
  sll_real64 :: mass

  sll_int32  :: dim_x
  sll_int32  :: dim_v

  ! plotting parameters
  sll_int32 :: plot_np_x    !< nb of points in the x  plotting grid for a (x,vx) plot
  sll_int32 :: plot_np_vx   !< nb of points in the vx plotting grid for a (x,vx) plot
  sll_real64 :: slice_y
  sll_real64 :: slice_vy

  sll_int32 :: plot_count

  sll_int32 :: io_stat
  logical    :: fail

  ! Note: the initialization and sampling of the lbf particle group here is copied from the test test_particle_sampling_interface
  ! in the pic_sampling library

  charge = -1.0_f64
  mass = 1.0_f64

  domain_is_x_periodic = .true.
  domain_is_y_periodic = .true.
  remapping_grid_eta_min = 0.0_f64
  remapping_grid_eta_max = 1.0_f64
  remapping_sparse_grid_max_levels = 5
  remap_degree = 3
  n_particles_x = 6
  n_particles_y = 6
  n_particles_vx = 6
  n_particles_vy = 6

  ! initialize the particle group, particles not sampled yet
  call sll_s_new_particle_group_2d2v_lbf_ptr( &
      particle_group, &
      charge,    &
      mass,      &
      domain_is_x_periodic,    &
      domain_is_y_periodic,    &
      remap_degree,    &
      remapping_grid_eta_min, &
      remapping_grid_eta_max, &
      remapping_sparse_grid_max_levels, &
      n_particles_x,  &
      n_particles_y,  &
      n_particles_vx, &
      n_particles_vy &
  )

  ! initialize the distribution
  distribution_str = "cossum_onegaussian"
  dim_x = 2
  dim_v = 2
  ! file must have a namelist of the form
  !   /cossum_onegaussian/ kx, alpha, v_thermal, v_mean
  ! with
  !   sll_real64 :: kx(dim_x)
  !   sll_real64 :: alpha
  !   sll_real64 :: v_thermal(dim_v)
  !   sll_real64 :: v_mean(dim_v)
  call sll_s_concatenate_filename_and_path( "distribution_for_the_test.nml", __FILE__,&
       filename)
  open(newunit = file_id, file=trim(filename), IOStat=io_stat)
  if (io_stat /= 0) then
    print*, 'test_particle_sampling_interface failed to open file ', filename
    STOP
  end if

  call sll_s_initial_distribution_new( trim(distribution_str), [dim_x, dim_v], file_id, distribution_params )
  close(file_id)

  ! deterministic sampling for lbf particles -- few arguments needed, as some are already encoded in the particle group
  call sll_s_sample_particle_group( &
      particle_group,  &
      distribution_params = distribution_params &
  )

  ! ok now we test the visualization interface
  plot_np_x  = 10
  plot_np_vx = 10
  slice_y  = 0.0_f64
  slice_vy  = 0.0_f64

  call plotting_params_2d%reset_params( &
      'f_slice', &
      plot_np_x, &
      plot_np_vx, &
      slice_y, &
      slice_vy &
      )
  plot_count = 1
  call sll_s_visualize_particle_group( particle_group, plotting_params_2d, plot_count)

  fail = .FALSE.

  ! maybe we should test something -- else :)

  if (fail .EQV. .FALSE.) then
     print*, 'PASSED'
  else
     print*, 'FAILED'
     stop
  end if

  call particle_group%free()
  deallocate(particle_group)
  call distribution_params%free()
  deallocate(distribution_params)

end program test_particle_visualization_interface
