program test_pic_particles
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#define NUM_PARTICLES 1000000_i64
#define GUARD_SIZE    100000_i64
#define PARTICLE_ARRAY_SIZE 1500000_i64

!!$  p_group => new_particle_4d_group( &
!!$       NUM_PARTICLES, &
!!$       PARTICLE_ARRAY_SIZE, &
!!$       GUARD_SIZE )
!!$
!!$  call sll_o_delete( p_group )

  print *, 'PASSED'

!  contains


end program test_pic_particles
