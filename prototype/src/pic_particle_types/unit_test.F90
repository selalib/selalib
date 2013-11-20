program particles_tester
  use sll_particle_group_2d_module
  implicit none

#define NUM_PARTICLES 1000000_i64
#define GUARD_SIZE    100000_i64
#define PARTICLE_ARRAY_SIZE 1500000_i64

  type(sll_particle_group_2d), pointer :: p_group

  p_group => new_particle_2d_group( &
       NUM_PARTICLES, &
       PARTICLE_ARRAY_SIZE, &
       GUARD_SIZE )

  call delete( p_group )

  print *, 'PASSED'
end program particles_tester
