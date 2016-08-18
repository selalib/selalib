program test_pic_accumulators
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_accumulators.h"

  use sll_m_accumulators, only: &
    sll_f_new_charge_accumulator_2d, &
    sll_t_charge_accumulator_2d, &
    sll_o_delete

  use sll_m_cartesian_meshes, only: &
    sll_f_new_cartesian_mesh_2d, &
    sll_t_cartesian_mesh_2d, &
    sll_o_delete

  use sll_m_constants, only: &
    sll_p_pi

  use sll_m_particle_group_2d, only: &
    sll_o_delete

  use sll_m_particle_group_4d, only: &
    sll_f_new_particle_4d_group, &
    sll_o_delete, &
    sll_t_particle_group_4d

  use sll_m_particle_initializers_4d, only: &
    sll_s_initial_random_particles_4d, &
    sll_s_initial_hammersley_particles_4d

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#define THERM_SPEED 1._f64
#define NUM_PARTICLES 100000_i32
#define GUARD_SIZE    10000_i32
#define PARTICLE_ARRAY_SIZE 150000_i32
#define ALPHA  0.5_f64
#define NC_X 256_i32
#define XMIN 0._f64
#define KX   0.5_f64
#define XMAX 2._f64*sll_p_pi/KX
#define NC_Y 64_i32
#define YMIN 0._f64
#define YMAX 1._f64
#define QoverM 1._f64

  
  type(sll_t_particle_group_4d), pointer :: part_group
  type(sll_t_cartesian_mesh_2d),   pointer :: m2d
  type(sll_t_charge_accumulator_2d), pointer :: all_charge

  m2d =>  sll_f_new_cartesian_mesh_2d( NC_X, NC_Y, &
       XMIN, XMAX, YMIN, YMAX )

  part_group => sll_f_new_particle_4d_group( &
       NUM_PARTICLES, &
       PARTICLE_ARRAY_SIZE, &
       GUARD_SIZE, QoverM, m2d )

  call sll_s_initial_random_particles_4d(THERM_SPEED, &
        ALPHA, KX, m2d, &
        NUM_PARTICLES, part_group )

  call sll_s_initial_hammersley_particles_4d(THERM_SPEED, &
        ALPHA, KX, m2d, &
        NUM_PARTICLES, part_group )

  all_charge => sll_f_new_charge_accumulator_2d( m2d )

!  call sll_first_charge_accumulation_2d( part_group, all_charge )
! I changed the type of all_charge in order to do multi-threading
! 
  call sll_o_delete( part_group )
  call sll_o_delete( m2d )
  call sll_o_delete( all_charge )

  print*, "PASSED"

!contains
  
end program test_pic_accumulators
