program test_pic_initializers_4d
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
  type(sll_t_cartesian_mesh_2d), pointer :: m2d
  sll_real64 :: mean_ref_v(1:2), variance_ref_v(1:2)
  sll_real64 :: mean_ref_x, variance_ref_x
  sll_real64 :: len_x, len_v
  logical    :: passed

  passed = .true.

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

  mean_ref_v     = [0._f64, 0._f64]! the 2d Gaussian for velocity
  variance_ref_v = [1._f64, 1._f64]!  
  mean_ref_x     = sll_p_pi/KX! the Landau 1d perturbation for positions
  variance_ref_x = (sll_p_pi**2/3._f64 + 2._f64*ALPHA)/(KX**2)
! the confidence interval for the mean
  len_x = 3._f64*sqrt(variance_ref_x)/sqrt(real(NUM_PARTICLES,f64))
  len_v = 3._f64*sqrt(variance_ref_v(1))/sqrt(real(NUM_PARTICLES,f64))

  print*, 'the Random initialization for the Landau damping' 
  call sll_s_initial_random_particles_4d(THERM_SPEED, &
        ALPHA, KX, m2d, NUM_PARTICLES, init_group_random )
  call test_mean_variance_xv (init_group_random, len_x, len_v)

  print*, 'the Hammersley initialization for the Landau damping' 
  call sll_s_initial_hammersley_particles_4d(THERM_SPEED, &
        ALPHA, KX, m2d, NUM_PARTICLES, init_group_hamm )
  call test_mean_variance_xv (init_group_hamm, len_x, len_v)

  call sll_o_delete( init_group_random )
  call sll_o_delete( init_group_hamm )
  call sll_o_delete( m2d )

  if ( passed .eqv. .true.) then
     print*, "PASSED"
  else 
     print*, "FAILED"
  endif

contains

  subroutine test_mean_variance_xv (part_group, tol_x, tol_v)
    type(sll_t_particle_group_4d), pointer, intent(in) :: part_group
    sll_real64, intent(in)                             :: tol_x, tol_v
    sll_real64 :: mean_v(1:2), variance_v(1:2)! estimators in velocity
    sll_real64 :: mean_x, variance_x!           estimators in position
    sll_int32  :: i
    sll_real64 :: out_v(1:2), out_x
    sll_real64 :: x, y

    mean_x = 0._f64
    mean_v = 0._f64
    do i = 1, NUM_PARTICLES
       GET_PARTICLE_POSITION( part_group%p_list(i), m2d, x, y)
       mean_x = mean_x + x
       mean_v(1) = mean_v(1) + part_group%p_list(i)%vx
       mean_v(2) = mean_v(2) + part_group%p_list(i)%vy
    end do
    mean_v = mean_v/real(NUM_PARTICLES,f64)
    mean_x = mean_x/real(NUM_PARTICLES,f64)

    out_v = mean_v - mean_ref_v
    out_x = mean_x - mean_ref_x

    if ( (abs(out_v(1)) > tol_v) .or. (abs(out_v(2)) > tol_v) ) then
       passed = .false.
       print*, 'Error in the expected value in velocity'
    endif
    if ( abs(out_x) > tol_x ) then
       passed = .false.
       print*, 'Error in the expected value in position'
    endif

    variance_v = 0._f64
    variance_x = 0._f64
    do i = 1, NUM_PARTICLES
       variance_v(1) = variance_v(1) + (part_group%p_list(i)%vx - mean_v(1))**2
       variance_v(2) = variance_v(2) + (part_group%p_list(i)%vy - mean_v(2))**2
       GET_PARTICLE_POSITION( part_group%p_list(i), m2d, x, y)
       variance_x = variance_x + (x - mean_x)**2
    end do
    variance_v = variance_v/real(NUM_PARTICLES-1,f64)
    variance_x = variance_x/real(NUM_PARTICLES-1,f64)

    print*, 'the error of the variance approximation in X=', &
         variance_x - variance_ref_x
    print*, 'the error of the variance approximation in V=', &
         variance_v - variance_ref_v

  end subroutine test_mean_variance_xv

end program test_pic_initializers_4d
