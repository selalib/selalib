program test_pic_initializers_2d
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_particle_representation.h"

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

  use sll_m_particle_initializers_2d, only: &
    sll_s_initial_random_particles_2d

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
  
  type(sll_t_particle_group_2d), pointer :: init_group_GC
  type(sll_t_cartesian_mesh_2d), pointer :: m2d
  sll_real64 :: mean_ref, variance_ref, len
  logical    :: passed

  passed = .true.

  m2d =>  sll_f_new_cartesian_mesh_2d( NC_X, NC_Y, &
       XMIN, XMAX, YMIN, YMAX )

  init_group_GC => sll_f_new_particle_2d_group( &
       NUM_PARTICLES, &
       PARTICLE_ARRAY_SIZE, &
       GUARD_SIZE, QoverM, m2d )

  mean_ref     = sll_p_pi/KX! the Landau 1d perturbation for positions
  variance_ref = (sll_p_pi**2/3._f64 + 2._f64*ALPHA)/(KX**2)
! the confidence interval for the mean
  len = 3._f64*sqrt(variance_ref)/sqrt(real(NUM_PARTICLES,f64))

  call sll_s_initial_random_particles_2d( ALPHA, KX, &
       m2d, NUM_PARTICLES, init_group_GC )

  call test_mean_variance_x (init_group_GC, len)

  call sll_o_delete( init_group_GC )
  call sll_o_delete( m2d )

  if ( passed .eqv. .true.) then
     print*, "PASSED"
  else 
     print*, "FAILED"
  endif

contains

  subroutine test_mean_variance_x (part_group, tolerance)
    type(sll_t_particle_group_2d), pointer, intent(in) :: part_group
    sll_real64, intent(in)                             :: tolerance
    sll_real64 :: mean, variance! estimators of the mean and the variance
    sll_int32  :: i
    sll_real64 :: out
    sll_real64 :: x, y
    
    mean = 0._f64
    do i = 1, NUM_PARTICLES
       GET_PARTICLE_POSITION( part_group%p_list(i), m2d, x, y)
       mean = mean + x
    end do
    mean = mean/real(NUM_PARTICLES,f64)
    out = mean - mean_ref

    if ( abs(out) > tolerance ) then
       passed = .false.
       print*, 'Error in the expected value'
    endif
    
   variance = 0._f64
    do i = 1, NUM_PARTICLES
       GET_PARTICLE_POSITION( part_group%p_list(i), m2d, x, y)
       variance = variance + (x - mean)**2
    end do
    variance = variance/real(NUM_PARTICLES-1,f64)
    print*, 'the error of the variance approximation in X=', variance - variance_ref

  end subroutine test_mean_variance_x

end program test_pic_initializers_2d
