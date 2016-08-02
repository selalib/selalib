program test_particle_group_2d2v_lbf
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_particle_group_2d2v_lbf, only: &
    sll_s_new_particle_group_2d2v_lbf_ptr

    !    , &
    !    sll_t_particle_group_2d2v_lbf

  use sll_m_particle_group_base, only: &
    sll_c_particle_group_base

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


  !  type(sll_t_particle_group_2d2v_lbf) :: particle_group
  !   class(sll_c_particle_group_base), pointer :: pgp

  class(sll_c_particle_group_base), pointer :: particle_group

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

  sll_int32  :: n_particles
  ! sll_int32  :: n_total_particles
  sll_real64 :: charge
  sll_real64 :: mass
  ! sll_int32  :: n_weights

  sll_int32  :: i_part
  sll_real64 :: x(3)
  logical    :: fail

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

  ! initialize the particle group, no particle are sampled yet
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

  ! artificial sampling of particles, just for the test
  n_particles = particle_group%n_particles
  do i_part = 1, n_particles
     call particle_group%set_x(i_part, [real(i_part,f64), 0.0_f64, 0.0_f64])
     call particle_group%set_v(i_part, [real(i_part,f64)**2, 0.0_f64, 0.0_f64])
     call particle_group%set_weights(i_part, [real(i_part/n_particles, f64)])
  end do

  fail = .FALSE.
  i_part = 4
  x =  particle_group%get_x(i_part)
  if ( abs(x(1)- real(i_part,f64))> 1E-15) then
     fail = .TRUE.
  end if
  x = particle_group%get_v(i_part)
  if ( abs(x(1)- real(i_part,f64)**2)> 1E-15) then
     fail = .TRUE.
  end if
  x(1:1) = particle_group%get_charge(i_part)
  if ( abs(x(1)- charge*real(i_part/n_particles, f64))> 1E-15) then
     fail = .TRUE.
  end if
  x(1:1) = particle_group%get_mass(i_part)
  if ( abs(x(1)- mass*real(i_part/n_particles, f64))> 1E-15) then
     fail = .TRUE.
  end if

  if( particle_group%n_particles .ne. n_particles_x * n_particles_y * n_particles_vx * n_particles_vy )then
    fail = .TRUE.
  end if

  if (fail .EQV. .FALSE.) then
     print*, 'PASSED'
  else
     print*, 'FAILED'
     stop
  end if

  call particle_group%free()

end program test_particle_group_2d2v_lbf
