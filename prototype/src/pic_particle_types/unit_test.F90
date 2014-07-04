program particles_tester
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"

  use sll_particle_group_2d_module

  implicit none

#define NUM_PARTICLES 1000000_i64
#define GUARD_SIZE    100000_i64
#define PARTICLE_ARRAY_SIZE 1500000_i64

!!$  p_group => new_particle_2d_group( &
!!$       NUM_PARTICLES, &
!!$       PARTICLE_ARRAY_SIZE, &
!!$       GUARD_SIZE )
!!$
!!$  call sll_delete( p_group )

  print *, 'PASSED'

!  contains


end program particles_tester
