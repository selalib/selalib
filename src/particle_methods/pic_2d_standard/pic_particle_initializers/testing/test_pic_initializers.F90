program test_pic_initializers
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"

  use sll_m_constants, only: sll_pi
  use sll_m_particle_group_4d
  use sll_m_particle_group_2d
  use sll_m_particle_initializers_2d
  use sll_m_particle_initializers_4d
  use sll_m_cartesian_meshes

  implicit none

#define THERM_SPEED 1._f64
#define NUM_PARTICLES 100000_i32
#define GUARD_SIZE 10000_i32
#define PARTICLE_ARRAY_SIZE 150000_i32
#define ALPHA 0.5_f64
#define NC_X 128_i32
#define XMIN 0._f64
#define KX   0.5_f64
#define XMAX 2._f64*sll_pi/KX
#define NC_Y 32_i32
#define YMIN 0._f64
#define YMAX 1._f64
#define QoverM 1._f64
  
  type(sll_particle_group_4d), pointer :: init_group
  type(sll_particle_group_2d), pointer :: init_group_GC
  type(sll_cartesian_mesh_2d),   pointer :: m2d

  m2d =>  new_cartesian_mesh_2d( NC_X, NC_Y, &
       XMIN, XMAX, YMIN, YMAX )

  init_group => new_particle_4d_group( &
       NUM_PARTICLES, &
       PARTICLE_ARRAY_SIZE, &
       GUARD_SIZE, QoverM, m2d )

  call sll_initial_particles_4d(THERM_SPEED, &
        ALPHA, KX, m2d, NUM_PARTICLES, init_group )

  init_group_GC => new_particle_2d_group( &
       NUM_PARTICLES, &
       PARTICLE_ARRAY_SIZE, &
       GUARD_SIZE, QoverM, m2d )
  call sll_initial_particles_2d( ALPHA, &
       KX, m2d, NUM_PARTICLES, init_group_GC )
  
  call sll_delete( init_group )
  call sll_delete( m2d )
  print*, "PASSED"

!contains
  
end program test_pic_initializers
