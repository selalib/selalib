program test_pic_initializers
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_cartesian_meshes, only: &
    sll_f_new_cartesian_mesh_2d, &
    sll_t_cartesian_mesh_2d, &
    sll_o_delete

  use sll_m_constants, only: &
    sll_p_pi

  use sll_m_particle_group_2d, only: &
    sll_f_new_particle_2d_group, &
    sll_o_delete, &
    sll_t_particle_group_2d

  use sll_m_particle_group_4d, only: &
    sll_f_new_particle_4d_group, &
    sll_o_delete, &
    sll_t_particle_group_4d

  use sll_m_particle_initializers_2d, only: &
    sll_s_initial_particles_2d

  use sll_m_particle_initializers_4d, only: &
    sll_s_initial_random_particles_4d, &
    sll_s_initial_hammersley_particles_4d

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#define THERM_SPEED 1._f64
#define NUM_PARTICLES 50000_i32
#define GUARD_SIZE 10000_i32
#define PARTICLE_ARRAY_SIZE 15000000_i32
#define ALPHA 0.5_f64
#define NC_X 128_i32
#define XMIN 0._f64
#define KX   0.5_f64
#define XMAX 2._f64*sll_p_pi/KX
#define NC_Y 32_i32
#define YMIN 0._f64
#define YMAX 1._f64
#define QoverM 1._f64
  
  type(sll_t_particle_group_4d), pointer :: init_group_random
  type(sll_t_particle_group_4d), pointer :: init_group_hamm
  type(sll_t_particle_group_2d), pointer :: init_group_GC
  type(sll_t_cartesian_mesh_2d), pointer :: m2d
  sll_real64 :: mean_ref(1:2), stddev_ref(1:2)
    
  m2d =>  sll_f_new_cartesian_mesh_2d( NC_X, NC_Y, &
       XMIN, XMAX, YMIN, YMAX )

  init_group_random => sll_f_new_particle_4d_group( &
       NUM_PARTICLES, &
       PARTICLE_ARRAY_SIZE, &
       GUARD_SIZE, QoverM, m2d )

  init_group_hamm => sll_f_new_particle_4d_group( &
       NUM_PARTICLES, &
       PARTICLE_ARRAY_SIZE, &
       GUARD_SIZE, QoverM, m2d )

  mean_ref   = [0._f64, 0._f64]! the 2d Gaussian for velocity
  stddev_ref = [1._f64, 1._f64]!  

  print*, 'the Random initialization' 
  call sll_s_initial_random_particles_4d(THERM_SPEED, &
        ALPHA, KX, m2d, NUM_PARTICLES, init_group_random )
  call test_mean_stddev(init_group_random, 3._f64/sqrt(real(NUM_PARTICLES,f64)))

  print*, 'the HAMMERSLEY initialization' 
  call sll_s_initial_hammersley_particles_4d(THERM_SPEED, &
        ALPHA, KX, m2d, NUM_PARTICLES, init_group_hamm )
  call test_mean_stddev(init_group_hamm, 3._f64/sqrt(real(NUM_PARTICLES,f64)))

  init_group_GC => sll_f_new_particle_2d_group( &
       NUM_PARTICLES, &
       PARTICLE_ARRAY_SIZE, &
       GUARD_SIZE, QoverM, m2d )

  call sll_s_initial_particles_2d( ALPHA, &
       KX, m2d, NUM_PARTICLES, init_group_GC )

  call sll_o_delete( init_group_random )
  call sll_o_delete( init_group_hamm )
  call sll_o_delete( init_group_GC )
  call sll_o_delete( m2d )

  print*, "PASSED"

contains

  subroutine test_mean_stddev (part_group, tolerance)
    type(sll_t_particle_group_4d), pointer, intent(in) :: part_group
    sll_real64, intent(in)                             :: tolerance
    sll_real64 :: mean(1:2), stddev(1:2)
    sll_int32  :: i
    sll_real64 :: out(1:2)

    mean = 0._f64
    do i = 1, NUM_PARTICLES
       mean(1) = mean(1) + part_group%p_list(i)%vx
       mean(2) = mean(2) + part_group%p_list(i)%vy
    end do
    mean = mean/real(NUM_PARTICLES,f64)

    out = mean - mean_ref
    if ( (abs(out(1)) > tolerance) .or. (abs(out(2)) > tolerance) ) then
       print*, 'FAILED'
    endif

   stddev = 0._f64
    do i = 1, NUM_PARTICLES
       stddev(1) = stddev(1) + (part_group%p_list(i)%vx - mean(1))**2
       stddev(2) = stddev(2) + (part_group%p_list(i)%vy - mean(2))**2
    end do
    stddev = stddev/real(NUM_PARTICLES-1,f64)
    print*, 'Std deviation error=', stddev - stddev_ref

  end subroutine test_mean_stddev

end program test_pic_initializers
