!> @ingroup pic_sampling
!> @author Katharina Kormann
!> @brief Particle initializer class with various functions to initialize a particle.
!> @details ...
module sll_m_particle_sampling
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_errors.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_working_precision.h"

  use sll_m_control_variate, only : &
       sll_t_control_variate

  use sll_m_initial_distribution, only: &
       sll_c_distribution_params, &
       sll_t_cos_gaussian
  
  use sll_m_particle_group_base, only: &
    sll_c_particle_group_base

  use sll_m_prob, only: &
    sll_s_normal_cdf_inv

  use sll_m_sobol, only: &
    sll_s_i8_sobol

  implicit none

  public :: &
       sll_t_particle_sampling, &
       sll_p_particle_sampling_random, &
       sll_p_particle_sampling_sobol, &
       sll_p_particle_sampling_random_symmetric, &
       sll_p_particle_sampling_sobol_symmetric

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  sll_int32, parameter :: sll_p_particle_sampling_random = 0
  sll_int32, parameter :: sll_p_particle_sampling_sobol = 1
  sll_int32, parameter :: sll_p_particle_sampling_random_symmetric = 2
  sll_int32, parameter :: sll_p_particle_sampling_sobol_symmetric = 3

  sll_int32, parameter :: sll_p_random_numbers = 0
  sll_int32, parameter :: sll_p_sobol_numbers = 1

!!$  type, abstract :: sll_c_number_generator
!!$     logical :: symmetric
!!$
!!$   contains
!!$     procedure(signature_get), deferred :: get_numbers
!!$
!!$     
!!$  end type sll_c_number_generator
!!$
!!$  type, extends(sll_c_number_generator) :: sll_t_number_generator_random
!!$     sll_int32, allocatable :: seed(:)
!!$
!!$   contains
!!$     procedure :: init => init_random
!!$     procedure :: get_numbers => get_random
!!$
!!$  end type sll_t_number_generator_random

  type :: sll_t_particle_sampling
     logical :: symmetric
     sll_int32               :: random_numbers
     sll_int32,  allocatable :: random_seed(:)
     sll_int64               :: sobol_seed

   contains
     procedure :: init => init_particle_sampling
     procedure :: sample => sample_particle_sampling
     procedure :: sample_cv => sample_cv_particle_sampling
     procedure :: free => free_particle_sampling

  end type sll_t_particle_sampling
  

contains

  subroutine free_particle_sampling( self )
    class( sll_t_particle_sampling ), intent( inout ) :: self !< particle sampling object

    if (allocated(self%random_seed)) deallocate( self%random_seed )
    
  end subroutine free_particle_sampling
    
  subroutine init_particle_sampling( self, sampling_type, dims, n_particles_local, rank )
    class( sll_t_particle_sampling ), intent( out ) :: self !< particle sampling object
    sll_int32,                         intent(in )              :: sampling_type !< sampling_type
    sll_int32,                         intent( in )             :: dims(:)
    sll_int32,                         intent(inout )           :: n_particles_local !< number of particles on processor
    sll_int32, optional,               intent(in )              :: rank !< optional argument to set random seed dependent on processor rank

    ! local
    sll_int32 :: prank
    sll_int32 :: ncopies
    sll_int32 :: np, j, rnd_seed_size

    prank = 0
    if( present(rank)) prank = rank

    select case( sampling_type )
    case( sll_p_particle_sampling_random )
       self%symmetric = .false.
       self%random_numbers = sll_p_random_numbers
    case( sll_p_particle_sampling_sobol )
       self%symmetric = .false.
       self%random_numbers = sll_p_sobol_numbers       
    case( sll_p_particle_sampling_random_symmetric )
       self%symmetric = .true.
       self%random_numbers = sll_p_random_numbers
    case( sll_p_particle_sampling_sobol_symmetric )
       self%symmetric = .true.
       self%random_numbers = sll_p_sobol_numbers
    case default
       SLL_ERROR("init_particle_sampling","Sampling type not implemented")
    end select

    ! Make sure that the particle number is conforming with symmetric sampling (if necessary)
    if (self%symmetric) then
       ncopies = 2**(sum(dims))
       np = modulo(n_particles_local,ncopies)
       if ( np .ne. 0 ) then
          n_particles_local = n_particles_local + np
          
       end if
    else
       ncopies = 1
    end if

    ! Set random numbers to default values (overwrite manually if other values are required)
    select case ( self%random_numbers )
    case( sll_p_random_numbers )
       call random_seed(size=rnd_seed_size)
       allocate(self%random_seed(rnd_seed_size))
       do j=1, rnd_seed_size
          self%random_seed(j) = (-1)**j*(100 + 15*j)*(2*prank + 1)
       end do
    case( sll_p_sobol_numbers )      
       self%sobol_seed = int(10 + prank*n_particles_local/ncopies, 8)
    end select
    
    
  end subroutine init_particle_sampling

  subroutine sample_cv_particle_sampling( self, particle_group, params, xmin, Lx, control_variate, time )
    class( sll_t_particle_sampling ), intent( inout ) :: self !< particle sampling object
    class(sll_c_particle_group_base), target, intent(inout)        :: particle_group
    class( sll_c_distribution_params ),  target,     intent( in )      :: params
    sll_real64,                        intent(in)               :: xmin(:) !< lower bound of the domain
    sll_real64,                        intent(in)               :: Lx(:) !< length of the domain.
    class(sll_t_control_variate),      intent(in)               :: control_variate
    sll_real64, optional,              intent(in)               :: time

    sll_real64 :: vi(3), xi(3), wi(particle_group%n_weights)
    sll_int32  :: i_part
    sll_real64 :: time0

    time0 = 0.0_f64
    if( present(time) ) time0 = time

    
    ! First sample the particles and set weights for full f
    call self%sample( particle_group, params, xmin, Lx )

    ! Fill wi(2) with value of initial distribution at initial positions (g0)
    ! and wi(3) with the delta f weights
    do i_part = 1, particle_group%n_particles
       xi = particle_group%get_x( i_part )
       vi = particle_group%get_v( i_part )
       wi = particle_group%get_weights( i_part )
       ! TODO: Distinguish here between different initial sampling distributions
       wi(2) = params%evalv( vi(1:params%dims(2)) )/product(Lx)
       wi(3) = control_variate%update_df_weight( xi, vi, time0, wi(1), wi(2) )
       call particle_group%set_weights( i_part, wi )
    end do
    

  end subroutine sample_cv_particle_sampling
  

  subroutine sample_particle_sampling( self, particle_group, params, xmin, Lx )
    class( sll_t_particle_sampling ), intent( inout ) :: self !< particle sampling object
    class(sll_c_particle_group_base), target, intent(inout)        :: particle_group
    class( sll_c_distribution_params ),  target,     intent( in )      :: params
    sll_real64,                        intent(in)               :: xmin(:) !< lower bound of the domain
    sll_real64,                        intent(in)               :: Lx(:) !< length of the domain.

    select type( params )
    type is( sll_t_cos_gaussian)
    
       if( self%symmetric .eqv. .false. ) then
          call sample_particle_sampling_all( self, particle_group, params, xmin, Lx )
       else
          if ( params%dims(1) == 1 .and. params%dims(2) == 2 ) then
             call sample_particle_sampling_sym_1d2v( self, particle_group, params, xmin, Lx )
          else
             SLL_ERROR("sample_particle_sampling", "symmetric sampling not implemented for given dimension")
          end if
       end if
    end select
    
  end subroutine sample_particle_sampling

  subroutine sample_particle_sampling_all( self,  particle_group, params, xmin, Lx )
    type( sll_t_particle_sampling ), intent( inout ) :: self !< particle sampling object
    class(sll_c_particle_group_base), target, intent(inout)        :: particle_group
    class( sll_t_cos_gaussian ),  target,     intent( in )      :: params
    sll_real64,                        intent(in)               :: xmin(:) !< lower bound of the domain
    sll_real64,                        intent(in)               :: Lx(:) !< length of the domain.

    sll_int32 :: n_rnds
    sll_real64                                         :: x(3),v(3)
    sll_int32                                          :: i_part
    sll_int32                                          :: i_v, i_gauss
    sll_real64, allocatable                            :: rdn(:)
    sll_real64                                         :: wi(particle_group%n_weights)
    sll_real64                                         :: rnd_no
    sll_real64                                         :: delta(params%n_gaussians)

    n_rnds = 0
    if ( params%n_gaussians > 1 ) then
       n_rnds = 1
    end if

    do i_v=1,params%n_gaussians
       delta(i_v) = sum(params%delta(1:i_v))
    end do
    
    
    n_rnds = n_rnds+params%dims(1)+params%dims(2)
    allocate( rdn(params%dims(1)+params%dims(2)+1) )
    rdn = 0.0_f64
    
    ! 1/Np in common weight
    call particle_group%set_common_weight &
         (1.0_f64/real(particle_group%n_total_particles, f64))

    if ( self%random_numbers == sll_p_random_numbers ) then
       call random_seed(put=self%random_seed)
    end if
   
    do i_part = 1, particle_group%n_particles
       ! Generate Random or Sobol numbers on [0,1]
       select case( self%random_numbers )
       case( sll_p_sobol_numbers )
          call sll_s_i8_sobol( int(n_rnds,8), self%sobol_seed, rdn(1:n_rnds))
       case( sll_p_random_numbers )
          call random_number( rdn(1:n_rnds) )
       end select
   
       ! Transform rdn to the interval
       x(1:params%dims(1)) = xmin + Lx * rdn(1:params%dims(1))

       ! Set weight according to value of perturbation
       wi(1) = params%evalx(x(1:params%dims(1)))*product(Lx)

       ! Maxwellian distribution of the temperature
       do i_v = 1,params%dims(2)
          call sll_s_normal_cdf_inv( rdn(i_v+params%dims(1)), 0.0_f64, 1.0_f64, &
               v(i_v))
       end do
       ! For multiple Gaussian, draw which one to take
       rnd_no = rdn(params%dims(1)+params%dims(2)+1)
       i_gauss = 1
       do while( rnd_no > delta(i_gauss) )
          i_gauss = i_gauss+1
       end do
       v(1:params%dims(2)) = v(1:params%dims(2)) * params%v_thermal(:,i_gauss) + params%v_mean(:,i_gauss)
       
       ! Copy the generated numbers to the particle
       call particle_group%set_x(i_part, x)
       call particle_group%set_v(i_part, v)
       ! Set weights.
       call particle_group%set_weights(i_part, &
            wi)
       
    end do
       
    

  end subroutine sample_particle_sampling_all

subroutine sample_particle_sampling_sym_1d2v( self, particle_group, params, xmin, Lx )
    type( sll_t_particle_sampling ), intent( inout ) :: self !< particle sampling object
    class(sll_c_particle_group_base), target, intent(inout)        :: particle_group
    class( sll_t_cos_gaussian ),  target,     intent( in )      :: params
    sll_real64,                        intent(in)               :: xmin(:) !< lower bound of the domain
    sll_real64,                        intent(in)               :: Lx(:) !< length of the domain.

    sll_int32 :: n_rnds
    sll_real64                                         :: x(3),v(3)
    sll_int32                                          :: i_part
    sll_int32                                          :: i_v
    sll_real64, allocatable                            :: rdn(:)
    sll_real64                                         :: wi(1)
    sll_real64                                         :: rnd_no
    sll_int32                                          :: ip, i_gauss
    sll_real64                                         :: delta(params%n_gaussians)

    n_rnds = 0
    if ( params%n_gaussians > 1 ) then
       n_rnds = 1
    end if

    do i_v=1,params%n_gaussians
       delta(i_v) = sum(params%delta(1:i_v))
    end do
    
    n_rnds = n_rnds+params%dims(1)+params%dims(2)
    allocate( rdn(params%dims(1)+params%dims(2)+1) )
    rdn = 0.0_f64
    
    ! 1/Np in common weight
    call particle_group%set_common_weight &
         (1.0_f64/real(particle_group%n_total_particles, f64))

    if ( self%random_numbers == sll_p_random_numbers ) then
       call random_seed(put=self%random_seed)
    end if

    do i_part = 1, particle_group%n_particles
       ip = modulo(i_part, 8 )
       if ( ip == 1) then
          ! Generate Random or Sobol numbers on [0,1]
          select case( self%random_numbers )
          case( sll_p_sobol_numbers )
             call sll_s_i8_sobol( int(n_rnds,8), self%sobol_seed, rdn(1:n_rnds))
          case( sll_p_random_numbers )
             call random_number( rdn(1:n_rnds) )
          end select
          
          ! Transform rdn to the interval
          x(1:params%dims(1)) = xmin + Lx * rdn(1:params%dims(1))
          
          ! Set weight according to value of perturbation
          wi(1) = params%evalx(x(1:params%dims(1)))*product(Lx)
          
          ! Maxwellian distribution of the temperature
          do i_v = 1,params%dims(2)
             call sll_s_normal_cdf_inv( rdn(i_v+params%dims(1)), 0.0_f64, 1.0_f64, &
                  v(i_v))
          end do
          ! For multiple Gaussian, draw which one to take
          rnd_no = rdn(params%dims(1)+params%dims(2)+1)
          i_gauss = 1
          do while( rnd_no > delta(i_gauss) )
             i_gauss = i_gauss+1
          end do
          v(1:params%dims(2)) = v(1:params%dims(2)) * params%v_thermal(:,i_gauss) + params%v_mean(:,i_gauss)
       elseif ( ip == 5 ) then
          x(1) = Lx(1) - x(1) + 2.0_f64*xmin(1)
       elseif ( modulo(ip,2) == 0 ) then
          v(1) = -v(1) + 2.0_f64*params%v_mean(1,i_gauss)
       else          
          v(2) = -v(2) + 2.0_f64*params%v_mean(2,i_gauss)
       end if
          
       ! Copy the generated numbers to the particle
       call particle_group%set_x(i_part, x)
       call particle_group%set_v(i_part, v)
       ! Set weights.
       call particle_group%set_weights(i_part, &
            wi)
       
    end do    

  end subroutine sample_particle_sampling_sym_1d2v
  

  ! From here old

  subroutine sll_s_particle_sampling_sobol_gaussian_cos_weights( &
       particle_group, &
       params, &
       xmin, &
       Lx, &
       seed )
    class(sll_c_particle_group_base),  intent(inout)            :: particle_group
    class( sll_t_cos_gaussian ),       intent( in )             :: params
    sll_real64,                        intent(in)               :: xmin(:) !< lower bound of the domain
    sll_real64,                        intent(in)               :: Lx(:) !< length of the domain.
    sll_int64,                         intent(inout)            :: seed !< Random seed.


    sll_int32 :: n_rnds
    sll_real64                                         :: x(3),v(3)
    sll_int32                                          :: i_part
    sll_int32                                          :: i_v
    sll_real64, allocatable                            :: rdn(:)
    sll_real64                                         :: wi
    sll_real64                                         :: rnd_no

    n_rnds = 0
    if ( params%n_gaussians > 1 ) then
       n_rnds = 1
    end if
    
    n_rnds = n_rnds+params%dims(1)+params%dims(2)
    allocate( rdn(params%dims(1)+params%dims(2)+1) )
    rdn = 0.0_f64
    
    ! 1/Np in common weight
    call particle_group%set_common_weight &
         (1.0_f64/real(particle_group%n_total_particles, f64))
    

    do i_part = 1, particle_group%n_particles
       ! Generate Sobol numbers on [0,1]
       call sll_s_i8_sobol( int(n_rnds,8), seed, rdn(1:n_rnds))

       ! Transform rdn to the interval
       x(1:params%dims(1)) = xmin + Lx * rdn(1:params%dims(1))

       ! Set weight according to value of perturbation
       wi = params%evalx(x(1:params%dims(1)))*product(Lx)

       ! Maxwellian distribution of the temperature
       do i_v = 1,params%dims(2)
          call sll_s_normal_cdf_inv( rdn(i_v+params%dims(1)), 0.0_f64, 1.0_f64, &
               v(i_v))
       end do
       ! For multiple Gaussian, draw which one to take
       rnd_no = rdn(params%dims(1)+params%dims(2)+1)
       i_v = 1
       do while( rnd_no > params%delta(i_v) )
          i_v = i_v+1
       end do
       v(1:params%dims(2)) = v(1:params%dims(2)) * params%v_thermal(:,i_v) + params%v_mean(:,i_v)
       
       ! Copy the generated numbers to the particle
       call particle_group%set_x(i_part, x)
       call particle_group%set_v(i_part, v)
       ! Set weights.
       call particle_group%set_weights(i_part, &
            [wi])
       
    end do
    
  end subroutine sll_s_particle_sampling_sobol_gaussian_cos_weights

  subroutine sll_s_particle_sampling_random_gaussian_cos_weights( &
       particle_group, &
       params, &
       xmin, &
       Lx, &
       seed )
    class(sll_c_particle_group_base),  intent(inout)            :: particle_group
    class( sll_t_cos_gaussian ),       intent( in )             :: params
    sll_real64,                        intent(in)               :: xmin(:) !< lower bound of the domain
    sll_real64,                        intent(in)               :: Lx(:) !< length of the domain.
    sll_int32,                         intent(inout)            :: seed(:) !< Random seed.


    sll_int32 :: n_rnds
    sll_real64                                         :: x(3),v(3)
    sll_int32                                          :: i_part
    sll_int32                                          :: i_v
    sll_real64, allocatable                            :: rdn(:)
    sll_real64                                         :: wi
    sll_real64                                         :: rnd_no


    ! Set random seed
    call random_seed(put=seed)
    
    n_rnds = 0
    if ( params%n_gaussians > 1 ) then
       n_rnds = 1
    end if
    
    n_rnds = n_rnds+params%dims(1)+params%dims(2)
    allocate( rdn(params%dims(1)+params%dims(2)+1) )
    rdn = 0.0_f64
    
    ! 1/Np in common weight
    call particle_group%set_common_weight &
         (1.0_f64/real(particle_group%n_total_particles, f64))
    

    do i_part = 1, particle_group%n_particles
       ! Generate Sobol numbers on [0,1]
       call random_number( rdn(1:n_rnds))

       ! Transform rdn to the interval
       x(1:params%dims(1)) = xmin + Lx * rdn(1:params%dims(1))

       ! Set weight according to value of perturbation
       wi = params%evalx(x(1:params%dims(1)))*product(Lx)

       ! Maxwellian distribution of the temperature
       do i_v = 1,params%dims(2)
          call sll_s_normal_cdf_inv( rdn(i_v+params%dims(1)), 0.0_f64, 1.0_f64, &
               v(i_v))
       end do
       ! For multiple Gaussian, draw which one to take
       rnd_no = rdn(params%dims(1)+params%dims(2)+1)
       i_v = 1
       do while( rnd_no > params%delta(i_v) )
          i_v = i_v+1
       end do
       v(1:params%dims(2)) = v(1:params%dims(2)) * params%v_thermal(:,i_v) + params%v_mean(:,i_v)
       
       ! Copy the generated numbers to the particle
       call particle_group%set_x(i_part, x)
       call particle_group%set_v(i_part, v)
       ! Set weights.
       call particle_group%set_weights(i_part, &
            [wi])
       
    end do
    
  end subroutine sll_s_particle_sampling_random_gaussian_cos_weights

end module sll_m_particle_sampling
