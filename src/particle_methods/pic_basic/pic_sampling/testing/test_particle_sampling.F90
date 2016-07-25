program test_particle_sampling
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_constants, only: &
       sll_p_pi, &
       sll_p_twopi

  use sll_m_initial_distribution, only: &
       sll_t_cos_gaussian
  
  use sll_m_particle_group_1d2v, only: &
    sll_s_new_particle_group_1d2v, &
    sll_t_particle_group_1d2v

  use sll_m_particle_group_base, only: &
    sll_c_particle_group_base

  use sll_m_particle_sampling
  
  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  class(sll_c_particle_group_base), allocatable:: particle_group
  type(sll_t_cos_gaussian)                     :: params
  sll_real64                                   :: xmin !< lower bound of the domain
  sll_real64                                   :: Lx !< length of the domain.
  sll_int32                                    :: n_particles      
  sll_int32                                    :: j, k
  sll_real64                                   :: mean_ref(3), sigma_ref(3), delta(2)

  n_particles = 80000

  allocate( sll_t_particle_group_1d2v :: particle_group)
  select type (particle_group)
  type is (sll_t_particle_group_1d2v)
     call particle_group%init( n_particles, n_particles, 1.0_F64, 1.0_F64, 1)
  end select
  xmin = 1.0_f64 
  Lx = sll_p_pi*4.0_f64
  mean_ref = [Lx*0.5_f64+xmin, 0.0_f64, 0.0_f64]
  


  ! Set initial parameters
  ! One gaussian parameters
  params%dims = [1,2]
  allocate( params%kx(params%dims(1),1) )
  allocate( params%alpha(params%dims(1)) )
  allocate( params%v_thermal(params%dims(2),1) )
  allocate( params%v_mean(params%dims(2),1) )
  allocate( params%normal(1) )
  allocate( params%delta(1) )
  params%n_cos = 1
  params%n_gaussians = 1
  params%kx = 0.5_f64
  params%alpha = 0.01_f64
  params%v_thermal(:,1) = [0.1_f64, 2.0_f64]
  params%v_mean = 0.0_f64
  params%delta = 0.0_f64
  params%normal = 1.0_f64/(sll_p_twopi**(0.5_f64*real(params%dims(2),f64))*&
       product(params%v_thermal(:,1)))

  sigma_ref = [Lx**2/12.0_f64, params%v_thermal(1,1)**2, params%v_thermal(2,1)**2 ]

  ! Landau
  print*, 'Sobol'
  call test( sll_p_particle_sampling_sobol, [1,2], particle_group,  &
       1d2/sqrt(real(n_particles,f64)) )
  print*, 'Sobol symmetric'
  call test( sll_p_particle_sampling_sobol_symmetric, [1,2], particle_group, &
       1d-12)
  print*, 'Random'
  call test( sll_p_particle_sampling_random, [1,2], particle_group,  &
       1d2/sqrt(real(n_particles,f64)) )
  print*, 'Random symmetric'
  call test( sll_p_particle_sampling_random_symmetric, [1,2], particle_group, &
       1d-12)
   
  
  ! Expected mean:
  ! 2Pi+1, 0, 0

  ! Expected variance:
  ! 16/12 pi**2 , 0.01, 4


  ! Next we consider an example with two Gaussians
  call params%free()

  ! Two gaussian parameters
  params%dims = [1,2]
  allocate( params%kx(params%dims(1),1) )
  allocate( params%alpha(params%dims(1)) )
  allocate( params%v_thermal(params%dims(2),2) )
  allocate( params%v_mean(params%dims(2),2) )
  allocate( params%normal(2) )
  allocate( params%delta(2) )
  params%n_cos = 1
  params%n_gaussians = 2
  params%kx = 0.5_f64
  params%alpha = 0.01_f64
  params%v_thermal(:,1) = [0.1_f64, 2.0_f64]
  params%v_thermal(:,2) = [2.0_f64, 2.0_f64]
  params%v_mean(:,1) = 0.0_f64
  params%v_mean(:,2) = 1.0_f64
  params%delta(1)  = 0.7_f64
  params%delta(2) = 0.3_f64
  params%normal(1) = 1.0_f64/(sll_p_twopi**(0.5_f64*real(params%dims(2),f64))*&
       product(params%v_thermal(:,1)))
  params%normal(2) = 1.0_f64/(sll_p_twopi**(0.5_f64*real(params%dims(2),f64))*&
       product(params%v_thermal(:,2)))

  
  mean_ref = [Lx*0.5_f64+xmin, 0.3_f64, 0.3_f64]
  sigma_ref(1) = Lx**2/12.0_f64
  
  do j=1,2
     sigma_ref(j+1) = - (params%delta(1)*params%v_mean(j,1)+&
          (params%delta(2)-params%delta(1))*params%v_mean(j,2))**2
     do k=1,2
        sigma_ref(j+1) = sigma_ref(j+1) + &
             params%delta(k) * (params%v_thermal(j,k)**2+params%v_mean(j,k)**2)
     end do
          
  end do
  
  ! TSI
  print*, 'Sobol'
  call test( sll_p_particle_sampling_sobol, [1,2], particle_group,  &
       1d2/sqrt(real(n_particles,f64)) )
  print*, 'Sobol symmetric'
  call test( sll_p_particle_sampling_sobol_symmetric, [1,2], particle_group, &
       1d2/sqrt(real(n_particles,f64)) )

  ! If we never stop due to tolerance not met, the test passed.
  print*, 'PASSED.'

  call particle_group%free()
  call params%free()
  deallocate(particle_group)

contains


 subroutine test( sampling_type, dims, particle_group, tolerance )
   sll_int32, intent( in    ) :: sampling_type
   sll_int32, intent( in    ) :: dims(2)
   class( sll_c_particle_group_base ), intent( inout ) :: particle_group
   sll_real64, intent( in   ) :: tolerance

   !local
   sll_real64 :: mean(3), sigma(3), xi(3)
   sll_int32  :: i_part
   type(sll_t_particle_sampling)                :: sampling
   sll_int32  :: n_particles

   n_particles = particle_group%n_particles
   
   ! Sobol symmetric
   ! Set sampling type: Sobol symmetric
   call sampling%init( sampling_type, dims, n_particles )
   ! Sample particles
   call sampling%sample( particle_group, params, [xmin], [Lx] )
   call sampling%free()
   
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
   
   if ( maxval(abs(mean)) > tolerance) then
      print*, 'FAILED'
      stop
   end if
   
 end subroutine test

end program test_particle_sampling
