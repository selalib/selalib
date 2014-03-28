program accumulate_tester
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"

  use sll_accumulators
  use sll_constants, only: sll_pi
  use sll_particle_group_2d_module
  use sll_particle_initializers
  use sll_logical_meshes


#define THERM_SPEED 1._f64
#define NUM_PARTICLES 100000_i64
#define GUARD_SIZE    10000_i64
#define PARTICLE_ARRAY_SIZE 150000_i64
#define ALPHA  0.5_f64
#define NC_X 256_i32
#define XMIN 0._f64
#define KX   0.5_f64
#define XMAX 2._f64*sll_pi/KX
#define NC_Y 64_i32
#define YMIN 0._f64
#define YMAX 1._f64

  
  implicit none
  type(sll_particle_group_2d), pointer :: part_group
  type(sll_logical_mesh_2d),   pointer :: m2d
  type(charge_accumulator_cell), dimension(:), pointer :: all_charge
  sll_int64 :: j

  part_group => new_particle_2d_group( &
       NUM_PARTICLES, &
       PARTICLE_ARRAY_SIZE, &
       GUARD_SIZE )
  
  m2d =>  new_logical_mesh_2d( NC_X, NC_Y, &
       XMIN, XMAX, YMIN, YMAX )

  call sll_initialize_some4Dfunction( THERM_SPEED, &
       ALPHA, KX, m2d, &
       NUM_PARTICLES, part_group )

  all_charge => new_accumulate_charge( NC_X*NC_Y )

  do j=1, part_group%number_particles
     call sll_accumulate_charge( part_group%p_list(j), &
                     all_charge( part_group%p_list(j)%ic ) )
  enddo

  call sll_delete( part_group )
  call delete( m2d )
  call sll_delete( all_charge )

  print*, "PASSED"

!contains
  
end program accumulate_tester
