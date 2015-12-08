program test_pic_accumulators
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_accumulators.h"

  use sll_m_accumulators, only: &
    new_charge_accumulator_2d, &
    sll_charge_accumulator_2d, &
    sll_delete

  use sll_m_cartesian_meshes, only: &
    new_cartesian_mesh_2d, &
    sll_cartesian_mesh_2d, &
    sll_delete

  use sll_m_constants, only: &
    sll_pi

  use sll_m_particle_group_2d, only: &
    sll_delete

  use sll_m_particle_group_4d, only: &
    new_particle_4d_group, &
    sll_delete, &
    sll_particle_group_4d

  use sll_m_particle_initializers_4d, only: &
    sll_initial_particles_4d

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
#define XMAX 2._f64*sll_pi/KX
#define NC_Y 64_i32
#define YMIN 0._f64
#define YMAX 1._f64
#define QoverM 1._f64

  
  type(sll_particle_group_4d), pointer :: part_group
  type(sll_cartesian_mesh_2d),   pointer :: m2d
  type(sll_charge_accumulator_2d), pointer :: all_charge

  m2d =>  new_cartesian_mesh_2d( NC_X, NC_Y, &
       XMIN, XMAX, YMIN, YMAX )

  part_group => new_particle_4d_group( &
       NUM_PARTICLES, &
       PARTICLE_ARRAY_SIZE, &
       GUARD_SIZE, QoverM, m2d )

  call sll_initial_particles_4d(THERM_SPEED, &
        ALPHA, KX, m2d, &
        NUM_PARTICLES, part_group )

  all_charge => new_charge_accumulator_2d( m2d )

!  call sll_first_charge_accumulation_2d( part_group, all_charge )
! I changed the type of all_charge in order to do multi-threading
! 
  call sll_delete( part_group )
  call sll_delete( m2d )
  call sll_delete( all_charge )

  print*, "PASSED"

!contains
  
end program test_pic_accumulators
