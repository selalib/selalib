program test_particle_initializer
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_constants, only: &
    sll_p_pi

  use sll_m_particle_group_1d2v, only: &
    sll_s_new_particle_group_1d2v, &
    sll_t_particle_group_1d2v

  use sll_m_particle_group_base, only: &
    sll_c_particle_group_base

  use sll_m_particle_initializer

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  class(sll_c_particle_group_base), allocatable:: particle_group
  sll_real64                                   :: landau_param(2) !< parameter defining the perturbation: landau_param(1)*cos(landau_param(2)*x1)
  sll_real64                                   :: xmin !< lower bound of the domain
  sll_real64                                   :: Lx !< length of the domain.
  sll_real64                                   :: thermal_velocity(2) !< Value of the thermal velocity along each dimension.
  sll_int64                                    :: rnd_seed !< Random seed.
  sll_int32, allocatable                       :: rnd_seeds(:) 
  sll_int32                                    :: rnd_seed_size
  sll_int32                                    :: n_particles      
  sll_int32                                    :: i_part, j
  sll_real64                                   :: xi(3), mean(3), sigma(3), mean_ref(3), sigma_ref(3)

  n_particles = 80000

  allocate( sll_t_particle_group_1d2v :: particle_group)
  select type (particle_group)
  type is (sll_t_particle_group_1d2v)
     call particle_group%init( n_particles, n_particles, 1.0_F64, 1.0_F64, 1)
  end select
  landau_param = [0.01_f64, 0.5_f64]
  xmin = 1.0_f64 
  Lx = sll_p_pi*4.0_f64
  thermal_velocity = [ 0.1_f64, 2.0_f64 ]
  rnd_seed = int(10,8)
  call random_seed(size=rnd_seed_size)
  SLL_ALLOCATE(rnd_seeds(rnd_seed_size), j)
  do j=1, rnd_seed_size
     rnd_seeds(j) = (-1)**j*(100 + 15*j)
  end do
  mean_ref = [Lx*0.5_f64+xmin, 0.0_f64, 0.0_f64]
  sigma_ref = [Lx**2/12.0_f64, thermal_velocity(1)**2, thermal_velocity(2)**2 ]

  ! Sobol

  call sll_s_particle_initialize_sobol_landau_1d2v(particle_group, &
            landau_param, xmin, Lx, &
            thermal_velocity, rnd_seed)

  mean = 0.0_f64
  sigma = 0.0_f64
  do i_part = 1, n_particles
     xi = particle_group%get_x(i_part)
     mean(1) = mean(1) + xi(1)
     xi = particle_group%get_v(i_part)
     mean(2) = mean(2) + xi(1)
     mean(3) = mean(3) + xi(2)
  end do
  mean = mean/real(n_particles, f64)
  do i_part = 1, n_particles
     xi = particle_group%get_x(i_part)
     sigma(1) = sigma(1) + (xi(1)-mean(1))**2
     xi = particle_group%get_v(i_part)
     sigma(2) = sigma(2) + (xi(1)-mean(2))**2
     sigma(3) = sigma(3) + (xi(2)-mean(3))**2
  end do

  sigma = sigma/real(n_particles-1, f64)

  mean = mean - mean_ref
  sigma = sigma - sigma_ref

  print*, 'Mean error:',  mean
  print*, 'Variance error:', sigma 

  if ( maxval(abs(mean)) > 1d2/sqrt(real(n_particles,f64))) then
     print*, 'FAILED'
     stop
  end if

  ! Sobol symmetric

  call sll_s_particle_initialize_sobol_landau_symmetric_1d2v(particle_group, &
            landau_param, xmin, Lx, &
            thermal_velocity, rnd_seed)

  mean = 0.0_f64
  sigma = 0.0_f64
  do i_part = 1, n_particles
     xi = particle_group%get_x(i_part)
     mean(1) = mean(1) + xi(1)
     xi = particle_group%get_v(i_part)
     mean(2) = mean(2) + xi(1)
     mean(3) = mean(3) + xi(2)
  end do
  mean = mean/real(n_particles, f64)
  do i_part = 1, n_particles
     xi = particle_group%get_x(i_part)
     sigma(1) = sigma(1) + (xi(1)-mean(1))**2
     xi = particle_group%get_v(i_part)
     sigma(2) = sigma(2) + (xi(1)-mean(2))**2
     sigma(3) = sigma(3) + (xi(2)-mean(3))**2
  end do

  sigma = sigma/real(n_particles-1, f64)


  mean = mean - mean_ref
  sigma = sigma - sigma_ref
  print*, 'Mean error:',  mean
  print*, 'Variance error:', sigma

  if ( maxval(abs(mean)) > 1d-12) then
     print*, 'FAILED'
     stop
  end if

  ! Random

  call sll_s_particle_initialize_random_landau_1d2v(particle_group, &
       landau_param, xmin, Lx, &
       thermal_velocity, rnd_seeds)

  mean = 0.0_f64
  sigma = 0.0_f64
  do i_part = 1, n_particles
     xi = particle_group%get_x(i_part)
     mean(1) = mean(1) + xi(1)
     xi = particle_group%get_v(i_part)
     mean(2) = mean(2) + xi(1)
     mean(3) = mean(3) + xi(2)
  end do
  mean = mean/real(n_particles, f64)
  do i_part = 1, n_particles
     xi = particle_group%get_x(i_part)
     sigma(1) = sigma(1) + (xi(1)-mean(1))**2
     xi = particle_group%get_v(i_part)
     sigma(2) = sigma(2) + (xi(1)-mean(2))**2
     sigma(3) = sigma(3) + (xi(2)-mean(3))**2
  end do

  sigma = sigma/real(n_particles-1, f64)


  mean = mean - mean_ref
  sigma = sigma - sigma_ref
  print*, 'Mean error:',  mean
  print*, 'Variance error:', sigma

  if ( maxval(abs(mean)) > 1d2/sqrt(real(n_particles,f64))) then
     print*, 'FAILED'
     stop
  end if

  ! Random symmetric

  call sll_s_particle_initialize_random_landau_symmetric_1d2v(particle_group, &
       landau_param, xmin, Lx, &
       thermal_velocity, rnd_seeds)

  mean = 0.0_f64
  sigma = 0.0_f64
  do i_part = 1, n_particles
     xi = particle_group%get_x(i_part)
     mean(1) = mean(1) + xi(1)
     xi = particle_group%get_v(i_part)
     mean(2) = mean(2) + xi(1)
     mean(3) = mean(3) + xi(2)
  end do
  mean = mean/real(n_particles, f64)
  do i_part = 1, n_particles
     xi = particle_group%get_x(i_part)
     sigma(1) = sigma(1) + (xi(1)-mean(1))**2
     xi = particle_group%get_v(i_part)
     sigma(2) = sigma(2) + (xi(1)-mean(2))**2
     sigma(3) = sigma(3) + (xi(2)-mean(3))**2
  end do

  sigma = sigma/real(n_particles-1, f64)

  mean = mean - mean_ref
  sigma = sigma - sigma_ref


  print*, 'Mean error:',  mean
  print*, 'Variance error:', sigma
  
  if ( maxval(abs(mean)) > 1d-12) then
     print*, 'FAILED'
     stop
  end if

  ! Expected mean:
  ! 2Pi+1, 0, 0

  ! Expected variance:
  ! 16/12 pi**2 , 0.01, 4


  ! If we never stop due to tolerance not met, the test passed.
  print*, 'PASSED.'

  call particle_group%free()
  deallocate(particle_group)
  deallocate(rnd_seeds)


end program test_particle_initializer
