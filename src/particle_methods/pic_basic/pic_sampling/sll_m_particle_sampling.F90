!> @ingroup pic_sampling
!> @author Katharina Kormann, Benedikt Perse
!> @brief Particle initializer class with various functions to initialize a particle.
!> @details ...
module sll_m_particle_sampling
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_errors.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_working_precision.h"

  use sll_m_collective, only: &
       sll_o_collective_allreduce, &
       sll_v_world_collective

  use sll_m_control_variate, only : &
       sll_t_control_variate

  use sll_m_initial_distribution, only: &
       sll_c_distribution_params, &
       sll_t_params_cos_gaussian, &
       sll_t_params_cos_gaussian_screwpinch, &
       sll_t_params_noise_gaussian

  use sll_m_particle_group_base, only: &
       sll_c_particle_group_base

  use sll_m_prob, only: &
       sll_s_normal_cdf_inv
  
  use sll_m_sobol, only: &
       sll_s_i8_sobol

  use sll_m_mapping_3d, only: &
       sll_t_mapping_3d

  use sll_mpi, only: &
       mpi_sum

  implicit none

  public :: &
       sll_t_particle_sampling, &
       sll_s_particle_sampling_randomized_weights

  private
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! Internal parameter to distinguish how to draw
  sll_int32, parameter :: sll_p_random_numbers = 0 !< draw random numbers
  sll_int32, parameter :: sll_p_sobol_numbers = 1 !< draw sobol numbers

  sll_int32, parameter :: sll_p_standard = 0 !< each particle is drawn separately
  sll_int32, parameter :: sll_p_symmetric_all = 1 !<  2^(dim_x+dim_v) points constructed from each drawn particle, all possible combinations of reflections along each direction
  sll_int32, parameter :: sll_p_symmetric_negative = 2 !< one additional point for each drawd point, the one reflected along all directions
  sll_int32, parameter :: sll_p_uniformx_negative = 3

  !> Data type for particle sampling
  type :: sll_t_particle_sampling
     sll_int32               :: symmetric !< Various cases of symmetric loading. 0: non-symmetric; symmetry is obtained by reflecting at the mean value of the Gaussian (for v) or the domain mid point (for x); 1: 2^(dim_x+dim_v) points with all possible combinations;  2: only one additional point per randomly drawn particle, the one reflected in each direction
     sll_int32               :: random_numbers !< How to draw (currently random or Sobol) defined by descriptors
     logical                 :: uniform = .false.!< Uniform loading.
     sll_int32,  allocatable :: random_seed(:) !< seed for random numbers
     sll_int32,  allocatable :: random_seed_start(:) !< save starting seed
     sll_int64               :: sobol_seed     !< seed for Sobol numbers
     sll_int64               :: sobol_seed_start !< save starting seed
     logical                 :: inverse = .false. !< true if the mapping has an analytical inverse
     logical                 :: xiprofile = .false. !< logical temperature profile

     logical                 :: delta_perturb = .false. !< true if delta perturbation for the velocity distribution
     sll_real64              :: eps(3) = [0.6, 0.6, 0.6] !< sigma quantile in which the velocity is perturbed
     sll_real64              :: a(3) = [0.0, 0.0, 0.0] !< value of the perturbation

   contains
     procedure :: init => init_particle_sampling !> Initializer
     procedure :: sample => sample_particle_sampling !> Sample particles
     procedure :: sample_cv => sample_cv_particle_sampling !> Sample particles and set delta f weights
     procedure :: free => free_particle_sampling !> Destructor
     procedure :: reset_seed_jump !< reset seed jump

  end type sll_t_particle_sampling


contains

  !> Descructor
  subroutine free_particle_sampling( self )
    class( sll_t_particle_sampling ), intent( inout ) :: self !< particle sampling object

    if (allocated(self%random_seed)) deallocate( self%random_seed )

  end subroutine free_particle_sampling

  !> Initializer
  subroutine init_particle_sampling( self, sampling_type, dims, n_particles_local, rank, delta_perturb, delta_eps )
    class( sll_t_particle_sampling ), intent(   out ) :: self !< particle sampling object
    character(len=*),                 intent( in    ) :: sampling_type !< sampling_type
    sll_int32,                        intent( in    ) :: dims(:) !< \a dims(1) number of spatial dimensions, \a dims(2) number of velocity dimensions
    sll_int32,                        intent( inout ) :: n_particles_local !< number of particles on processor
    sll_int32, optional,              intent( in    ) :: rank !< optional argument to set random seed dependent on processor rank
    logical, optional,                intent( in    ) :: delta_perturb !< delta perturbation of the velocity distribution
    sll_real64, optional,             intent( in    ) :: delta_eps(6) !< values for delta perturbation
    ! local variables
    sll_int32 :: prank
    sll_int32 :: ncopies
    sll_int32 :: np, j, rnd_seed_size

    if( present(delta_perturb) )then
       self%delta_perturb = delta_perturb
       if( delta_perturb ) then
          self%eps = delta_eps(1:3)
          self%a = delta_eps(4:6)
       end if
    end if

    prank = 0
    if( present(rank)) prank = rank

    select case( trim(sampling_type) )
    case( "particle_sampling_random" )
       self%symmetric = 0
       self%random_numbers = sll_p_random_numbers
    case( "particle_sampling_random_xiprofile" )
       self%symmetric = 0
       self%random_numbers = sll_p_random_numbers
       self%xiprofile = .true.
    case( "particle_sampling_random_inverse" )
       self%symmetric = 0
       self%random_numbers = sll_p_random_numbers
       self%inverse = .true.
    case( "particle_sampling_random_inverse_xiprofile" )
       self%symmetric = 0
       self%random_numbers = sll_p_random_numbers
       self%inverse = .true.
       self%xiprofile = .true.
    case( "particle_sampling_sobol" )
       self%symmetric = 0
       self%random_numbers = sll_p_sobol_numbers
    case( "particle_sampling_sobol_xiprofile" )
       self%symmetric = 0
       self%random_numbers = sll_p_sobol_numbers
       self%xiprofile = .true.
    case( "particle_sampling_sobol_inverse" )
       self%symmetric = 0
       self%random_numbers = sll_p_sobol_numbers
       self%inverse = .true.
    case( "particle_sampling_sobol_inverse_xiprofile" )
       self%symmetric = 0
       self%random_numbers = sll_p_sobol_numbers
       self%inverse = .true.
       self%xiprofile = .true.
    case( "particle_sampling_random_symmetric" )
       self%symmetric = 1
       self%random_numbers = sll_p_random_numbers
    case( "particle_sampling_random_symmetric_xiprofile" )
       self%symmetric = 1
       self%random_numbers = sll_p_random_numbers
       self%xiprofile = .true.
    case( "particle_sampling_random_symmetric_inverse" )
       self%symmetric = 1
       self%random_numbers = sll_p_random_numbers
       self%inverse = .true.
    case( "particle_sampling_random_symmetric_inverse_xiprofile" )
       self%symmetric = 1
       self%random_numbers = sll_p_random_numbers
       self%inverse = .true.
       self%xiprofile = .true.
    case( "particle_sampling_sobol_symmetric" )
       self%symmetric = 1
       self%random_numbers = sll_p_sobol_numbers
    case( "particle_sampling_sobol_symmetric_xiprofile" )
       self%symmetric = 1
       self%random_numbers = sll_p_sobol_numbers
       self%xiprofile = .true.
    case( "particle_sampling_sobol_symmetric_inverse" )
       self%symmetric = 1
       self%random_numbers = sll_p_sobol_numbers
       self%inverse = .true.
    case( "particle_sampling_sobol_symmetric_inverse_xiprofile" )
       self%symmetric = 1
       self%random_numbers = sll_p_sobol_numbers
       self%inverse = .true.
       self%xiprofile = .true.
    case( "particle_sampling_sobol_symmetric_uniform" )
       self%symmetric = 1
       self%uniform = .true.
       self%random_numbers = sll_p_sobol_numbers  
    case( "particle_sampling_random_negative" )
       self%symmetric = 2
       self%random_numbers = sll_p_random_numbers
    case( "particle_sampling_sobol_negative" )
       self%symmetric = 2
       self%random_numbers = sll_p_sobol_numbers  
    case( "particle_sampling_random_unix" )
       self%symmetric = 3
       self%random_numbers = sll_p_random_numbers
    case( "particle_sampling_sobol_unix" )
       self%symmetric = 3
       self%random_numbers = sll_p_sobol_numbers
    case default
       SLL_ERROR("init_particle_sampling","Sampling type not implemented")
    end select

    ! Make sure that the particle number is conforming with symmetric sampling (if necessary)
    if (self%symmetric == 1 .or. self%symmetric == 3 ) then
       ncopies = 2**(sum(dims))
       np = modulo(n_particles_local,ncopies)
       if ( np .ne. 0 ) then
          n_particles_local = n_particles_local + ncopies - np
          !print*,  n_particles_local, np
          print*, 'Warning: Number of particles has been increased so that the particle number is compatable with symmetric sampling.'
       end if
    elseif (self%symmetric == 2 ) then
       ncopies = 2
       np = modulo(n_particles_local,ncopies)
       if ( np .ne. 0 ) then
          n_particles_local = n_particles_local + ncopies - np
          !print*,  n_particles_local, np
          print*, 'Warning: Number of particles has been increased so that the particle number is compatable with symmetric sampling.'
       end if
    else
       ncopies = 1
    end if
    ! Set random numbers to default values (overwrite manually if other values are required)
    select case ( self%random_numbers )
    case( sll_p_random_numbers )
       call random_seed(size=rnd_seed_size)
       allocate(self%random_seed(1:rnd_seed_size))
       allocate(self%random_seed_start(1:rnd_seed_size))
       do j=1, rnd_seed_size
          self%random_seed(j) = (-1)**j*(100 + 15*j)*(2*prank + 1)
       end do
       self%random_seed_start = self%random_seed
    case( sll_p_sobol_numbers )      
       self%sobol_seed = int(10 + prank*n_particles_local/ncopies, 8)
       self%sobol_seed_start = self%sobol_seed
    end select


  end subroutine init_particle_sampling

  subroutine reset_seed_jump( self, jump )
    class( sll_t_particle_sampling ), intent( inout ) :: self !< particle sampling object
    sll_int32, intent( in ) :: jump !< jump

    ! Set random numbers to default values (overwrite manually if other values are required)
    select case ( self%random_numbers )
    case( sll_p_random_numbers )
       self%random_seed = self%random_seed_start+jump
    case( sll_p_sobol_numbers )      
       self%sobol_seed = self%sobol_seed_start+int(jump,8)
    end select

  end subroutine reset_seed_jump

  !> Sample with control variate (we assume that the particle weights are given in the following order:
  !> (full f weight, value of initial distribution at time 0, delta f weights)
  subroutine sample_cv_particle_sampling( self, particle_group, params, xmin, Lx, control_variate, time, map, lindf )
    class( sll_t_particle_sampling ), intent( inout ) :: self !< particle sampling object
    class(sll_c_particle_group_base), target, intent(inout)        :: particle_group !< particle group
    class( sll_c_distribution_params ),  target,     intent( inout )      :: params !< parameters for initial distribution
    sll_real64,                        intent(in)               :: xmin(:) !< lower bound of the domain
    sll_real64,                        intent(in)               :: Lx(:) !< length of the domain.
    class(sll_t_control_variate),      intent(in)               :: control_variate !< PIC control variate
    sll_real64, optional,              intent(in)               :: time !< initial time (default: 0)
    type(sll_t_mapping_3d), optional                :: map !< coordinate transformation
    logical,  optional                              :: lindf                                    
    !local variables
    sll_real64 :: vi(3), xi(3), x(3), wi(particle_group%n_weights)
    sll_int32  :: i_part, counter
    sll_real64 :: time0, Rx(3), q, m, g0, df_weight
    logical :: delta_perturb, lindf0

    time0 = 0.0_f64
    lindf0 = .false.
    if( present(time) ) time0 = time
    if( present(lindf) ) lindf0 = lindf

    delta_perturb = self%delta_perturb
    self%delta_perturb = .false.

    if( present(map)) then
       ! First sample the particles and set weights for full f
       call self%sample( particle_group, params, xmin, Lx, map )
    else
       ! First sample the particles and set weights for full f
       call self%sample( particle_group, params, xmin, Lx )
    end if

    counter = 0
    ! Fill g0 with value of initial distribution at initial positions (g0)
    ! and df_weight with the delta f weights
    do i_part = 1, particle_group%n_particles
       q = particle_group%species%q
       m = particle_group%species%m
       xi = particle_group%get_x( i_part )
       vi = particle_group%get_v( i_part )
       wi = particle_group%get_weights( i_part )
       ! TODO: Distinguish here between different initial sampling distributions
       if( self%xiprofile ) then
          if( present(map)) then
             if( self%inverse ) then
                g0 = params%eval_v_density( vi(1:params%dims(2)), xi(1:params%dims(1)), m=particle_group%species%m )/map%volume
             else
                g0 = params%eval_v_density( vi(1:params%dims(2)), xi(1:params%dims(1)), m=particle_group%species%m )/map%jacobian(xi)
             end if
             Rx = xi
             Rx(1) = map%get_x1(xi) + m*vi(2)/q
             Rx(1) = (Rx(1) - xmin(1))/Lx(1)
             df_weight = control_variate%update_df_weight( Rx, vi, time0, wi(1), g0 )
          else
             x = (xi - xmin)/Lx
             g0 = params%eval_v_density( vi(1:params%dims(2)), x(1:params%dims(1)), m=particle_group%species%m )/product(Lx)
             Rx = x
             Rx(1) = xi(1) + m*vi(2)/q
             Rx(1) = (Rx(1) - xmin(1))/Lx(1)
             df_weight = control_variate%update_df_weight( Rx, vi, time0, wi(1), g0 )
          end if
       else
          if( present(map)) then
             x = xi
             g0 = params%eval_v_density( vi(1:params%dims(2)), xi(1:params%dims(1)), m=particle_group%species%m )/map%jacobian(xi)
          else
             x = (xi - xmin)/Lx
             g0 = params%eval_v_density( vi(1:params%dims(2)), xi(1:params%dims(1)), m=particle_group%species%m )/product(Lx)
          end if
          df_weight = control_variate%update_df_weight( xi, vi, time0, wi(1), g0 )
          if(delta_perturb) then
             if( particle_group%species%q < 0._f64) then
                if( x(1) > 0.5_f64-self%eps(1) .and. x(1) < 0.5_f64+self%eps(1) )then
                   if( x(2) > 0.5_f64-self%eps(2) .and. x(2) < 0.5_f64+self%eps(2) ) then
                      if(xi(3) > 0.5_f64-self%eps(3) .and. x(3) < 0.5_f64+self%eps(3) )then
                         vi = vi + self%a
                         counter = counter + 1
                      end if
                   end if
                end if
             end if
          end if
          call particle_group%set_v(i_part, vi)
       end if
       if(lindf0) then
          wi(1) = df_weight
       else
          wi(2) = g0
          wi(3) = df_weight
       end if
       call particle_group%set_weights( i_part, wi )
    end do
    if(delta_perturb) then
       print*, 'delta_perturb_cv: ', self%eps, self%a, 'shifted particles:', counter
    end if
  end subroutine sample_cv_particle_sampling


!!$   !> Sample with control variate (we assume that the particle weights are given in the following order:
!!$  !> (full f weight, value of initial distribution at time 0, delta f weights)
!!$  subroutine sample_lindf_particle_sampling( self, particle_group, params, xmin, Lx, control_variate, time, map )
!!$    class( sll_t_particle_sampling ), intent( inout ) :: self !< particle sampling object
!!$    class(sll_c_particle_group_base), target, intent(inout)        :: particle_group !< particle group
!!$    class( sll_c_distribution_params ),  target,     intent( inout )      :: params !< parameters for initial distribution
!!$    sll_real64,                        intent(in)               :: xmin(:) !< lower bound of the domain
!!$    sll_real64,                        intent(in)               :: Lx(:) !< length of the domain.
!!$    class(sll_t_control_variate),      intent(in)               :: control_variate !< PIC control variate
!!$    sll_real64, optional,              intent(in)               :: time !< initial time (default: 0)
!!$    type(sll_t_mapping_3d), optional                :: map !< coordinate transformation 
!!$    !local variables
!!$    sll_real64 :: vi(3), xi(3), x(3), wi(1)
!!$    sll_int32  :: i_part, counter
!!$    sll_real64 :: time0, Rx(3), q, m
!!$    logical :: delta_perturb
!!$
!!$    time0 = 0.0_f64
!!$    if( present(time) ) time0 = time
!!$
!!$    delta_perturb = self%delta_perturb
!!$    self%delta_perturb = .false.
!!$
!!$    if( present(map)) then
!!$       ! First sample the particles and set weights for full f
!!$       call self%sample( particle_group, params, xmin, Lx, map )
!!$    else
!!$       ! First sample the particles and set weights for full f
!!$       call self%sample( particle_group, params, xmin, Lx )
!!$    end if
!!$
!!$    counter = 0
!!$    ! Fill wi(2) with value of initial distribution at initial positions (g0)
!!$    ! and wi(3) with the delta f weights
!!$    do i_part = 1, particle_group%n_particles
!!$       q = particle_group%species%q
!!$       m = particle_group%species%m
!!$       xi = particle_group%get_x( i_part )
!!$       vi = particle_group%get_v( i_part )
!!$       wi = particle_group%get_weights( i_part )
!!$       ! TODO: Distinguish here between different initial sampling distributions
!!$       if( self%xiprofile ) then
!!$          if( present(map)) then
!!$             Rx = xi
!!$             Rx(1) = map%get_x1(xi) + m*vi(2)/q
!!$             Rx(1) = (Rx(1) - xmin(1))/Lx(1)
!!$             if( self%inverse ) then
!!$                wi(1) = wi(1) -  params%eval_v_density( vi(1:params%dims(2)), Rx(1:params%dims(1)), m=particle_group%species%m )/params%eval_v_density( vi(1:params%dims(2)), xi(1:params%dims(1)), m=particle_group%species%m )*map%volume!control_variate%update_df_weight( Rx, vi, time0, wi(1), wi(2) )
!!$             else
!!$                wi(1) = wi(1) -  params%eval_v_density( vi(1:params%dims(2)), Rx(1:params%dims(1)), m=particle_group%species%m )/params%eval_v_density( vi(1:params%dims(2)), xi(1:params%dims(1)), m=particle_group%species%m )*map%jacobian(xi)
!!$             end if
!!$          else
!!$             x = (xi - xmin)/Lx
!!$             Rx = x
!!$             Rx(1) = xi(1) + m*vi(2)/q
!!$             Rx(1) = (Rx(1) - xmin(1))/Lx(1)
!!$             wi(1) = wi(1) - params%eval_v_density( vi(1:params%dims(2)), Rx(1:params%dims(1)), m=particle_group%species%m )/params%eval_v_density( vi(1:params%dims(2)), x(1:params%dims(1)), m=particle_group%species%m )*product(Lx)!control_variate%update_df_weight( Rx, vi, time0, wi(1), wi(2) )
!!$          end if
!!$       else
!!$          if(delta_perturb) then
!!$             if( particle_group%species%q < 0._f64) then
!!$                if( x(1) > 0.5_f64-self%eps(1) .and. x(1) < 0.5_f64+self%eps(1) )then
!!$                   if( x(2) > 0.5_f64-self%eps(2) .and. x(2) < 0.5_f64+self%eps(2) ) then
!!$                      if(xi(3) > 0.5_f64-self%eps(3) .and. x(3) < 0.5_f64+self%eps(3) )then
!!$                         vi = vi + self%a
!!$                         counter = counter + 1
!!$                      end if
!!$                   end if
!!$                end if
!!$             end if
!!$          end if
!!$          if( present(map)) then
!!$             wi(1) = wi(1) - map%jacobian(xi)!control_variate%update_df_weight( xi, vi, time0, wi(1), wi(2) )
!!$          else
!!$             wi(1) = wi(1) - product(Lx)
!!$          end if
!!$          call particle_group%set_v(i_part, vi)
!!$       end if
!!$
!!$       call particle_group%set_weights( i_part, wi )
!!$    end do
!!$    if(delta_perturb) then
!!$       print*, 'delta_perturb_cv: ', self%eps, self%a, 'shifted particles:', counter
!!$    end if
!!$  end subroutine sample_lindf_particle_sampling

  !> Sample from distribution defined by \a params
  subroutine sample_particle_sampling( self, particle_group, params, xmin, Lx, map )
    class( sll_t_particle_sampling ), intent( inout ) :: self !< particle sampling object
    class(sll_c_particle_group_base), target, intent(inout)        :: particle_group !< particle group
    class( sll_c_distribution_params ),  target,     intent( inout )      :: params !< parameters for initial distribution
    sll_real64,                        intent(in)               :: xmin(:) !< lower bound of the domain
    sll_real64,                        intent(in)               :: Lx(:) !< length of the domain.
    type(sll_t_mapping_3d), optional                :: map  !< coordinate transformation

    !select type( params )
    !type is( sll_t_params_cos_gaussian)

    if( present(map)) then

       if( self%symmetric == 0 ) then
          call sample_particle_sampling_all_trafo( self, particle_group, params, xmin, Lx, map )
       elseif ( self%symmetric == 1 ) then
          if ( params%dims(1) == 1 .and. params%dims(2) == 2 ) then
             call sample_particle_sampling_sym_1d2v_trafo( self, particle_group, params, Lx, map )
          elseif ( params%dims(1) == 3 .and. params%dims(2) == 3 ) then
             call sample_particle_sampling_sym_3d3v_trafo( self, particle_group, params, xmin, Lx, map )
          else
             SLL_ERROR("sample_particle_sampling", "symmetric sampling not implemented for given dimension")
          end if
       else
          SLL_ERROR("sample_particle_sampling", "this symmetric sampling is not implemented")
       end if

    else
       if( self%symmetric == 0 ) then
          call sample_particle_sampling_all( self, particle_group, params, xmin, Lx )
       elseif ( self%symmetric == 1 ) then
          if ( params%dims(1) == 1 .and. params%dims(2) == 2 ) then
             call sample_particle_sampling_sym_1d2v( self, particle_group, params, xmin, Lx )
          elseif ( params%dims(1) == 3 .and. params%dims(2) == 3 ) then
             if( self%uniform ) then
                call sample_particle_sampling_sym_uni_3d3v( self, particle_group, params, xmin, Lx )
             else
                call sample_particle_sampling_sym_3d3v( self, particle_group, params, xmin, Lx )
             end if
          else
             SLL_ERROR("sample_particle_sampling", "symmetric sampling not implemented for given dimension")
          end if
       elseif ( self%symmetric == 2) then
          if ( params%dims(1) == 1 .and. params%dims(2) == 2 ) then
             call sample_particle_sampling_sym2_1d2v( self, particle_group, params, xmin, Lx )
          else
             SLL_ERROR("sample_particle_sampling", "symmetric sampling not implemented for given dimension")
          end if
       elseif( self%symmetric == 3) then
          if ( params%dims(1) == 1 .and. params%dims(2) == 2 ) then
             call sample_particle_sampling_unix_1d2v( self, particle_group, params, xmin, Lx )
          end if
       end if
    end if
    !end select

  end subroutine sample_particle_sampling


  !> Helper function for pure sampling
  subroutine sample_particle_sampling_all( self,  particle_group, params, xmin, Lx )
    type( sll_t_particle_sampling ),            intent( inout ) :: self !< particle sampling object
    class(sll_c_particle_group_base), target,   intent( inout ) :: particle_group !< particle group
    class( sll_c_distribution_params ), target, intent( inout ) :: params !< distribution parameter
    sll_real64,                                 intent( in    ) :: xmin(:) !< lower bound of the domain
    sll_real64,                                 intent( in    ) :: Lx(:) !< length of the domain.
    !local variables
    sll_int32 :: n_rnds
    sll_real64                                         :: x(3), xi(3), v(3)
    sll_int32                                          :: i_part
    sll_int32                                          :: i_v, i_gauss
    sll_real64, allocatable                            :: rdn(:)
    sll_real64                                         :: wi(particle_group%n_weights)
    sll_real64                                         :: rnd_no
    sll_real64                                         :: delta(params%n_gaussians)
    sll_real64                                         :: common_weight
    sll_real64                                         :: r, Rx(3)

    n_rnds = 0
    if ( params%n_gaussians > 1 ) then
       n_rnds = 1
    end if

    do i_v=1,params%n_gaussians
       delta(i_v) = sum(params%delta(1:i_v))
    end do


    n_rnds = n_rnds+params%dims(1)+params%dims(2)
    allocate( rdn(1:params%dims(1)+params%dims(2)+1) )
    rdn = 0.0_f64

    ! 1/Np in common weight
    common_weight =  particle_group%get_common_weight()
    call particle_group%set_common_weight &
         (common_weight/real(particle_group%n_total_particles, f64))

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
       xi(1:params%dims(1)) = rdn(1:params%dims(1))
       x(1:params%dims(1)) = xmin + Lx * xi(1:params%dims(1))

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
       if( self%xiprofile ) then
          select type(p=>params)
          type is(sll_t_params_cos_gaussian_screwpinch)
             if( particle_group%species%q > 0._f64) then
                v(1:p%dims(2)) = v(1:p%dims(2)) * sqrt(p%profile%T_i(xi(1))/particle_group%species%m) + p%v_mean(:,i_gauss)
             else
                v(1:p%dims(2)) = v(1:p%dims(2)) * sqrt(p%profile%T_e(xi(1))/particle_group%species%m) + p%v_mean(:,i_gauss)
             end if

             Rx = xi
             Rx(1) = x(1) + particle_group%species%m*v(2)/particle_group%species%q
             Rx(1) = (Rx(1) - xmin(1))/Lx(1)

             !special perturbation in gyrocoordinates
             wi(1) = p%eval_x_density(xi(1:p%dims(1)), v(1:p%dims(2)))*&
                  p%eval_v_density(v(1:p%dims(2)), Rx(1:p%dims(1)), particle_group%species%m)/p%eval_v_density(v(1:p%dims(2)), xi(1:p%dims(1)), particle_group%species%m)  *product(Lx)
          type is(sll_t_params_noise_gaussian)
             if( particle_group%species%q > 0._f64) then
                v(1:p%dims(2)) = v(1:p%dims(2)) * sqrt(p%profile%T_i(xi(1))/particle_group%species%m) + p%v_mean(:,i_gauss)
             else
                v(1:p%dims(2)) = v(1:p%dims(2)) * sqrt(p%profile%T_e(xi(1))/particle_group%species%m) + p%v_mean(:,i_gauss)
             end if
             wi(1) = p%eval_x_density(xi(1:p%dims(1))) *product(Lx)
          end select
       else
          v(1:params%dims(2)) = v(1:params%dims(2)) * params%v_thermal(:,i_gauss) + params%v_mean(:,i_gauss)
          ! Set weight according to value of perturbation
          wi(1) = params%eval_x_density(x(1:params%dims(1)),v(1:params%dims(2)))*product(Lx) 
       end if
       ! Copy the generated numbers to the particle
       call particle_group%set_x(i_part, x)
       call particle_group%set_v(i_part, v)
       ! Set weights.
       call particle_group%set_weights(i_part, &
            wi)
    end do

    if( particle_group%species%q < 0._f64) then
       if( self%delta_perturb ) then
          call delta_function_perturbation( self, particle_group, xmin, Lx )
       end if
    end if

  end subroutine sample_particle_sampling_all


  !> Helper function for antithetic sampling in 1d2v
  subroutine sample_particle_sampling_sym_1d2v( self, particle_group, params, xmin, Lx )
    type( sll_t_particle_sampling ),          intent( inout ) :: self !< particle sampling object
    class(sll_c_particle_group_base), target, intent( inout ) :: particle_group !< particle group
    class( sll_c_distribution_params ), target, intent( in )  :: params !< distribution parameter
    sll_real64,                        intent(in)             :: xmin(:) !< lower bound of the domain
    sll_real64,                        intent(in)             :: Lx(:) !< length of the domain.
    !local variables
    sll_int32 :: n_rnds
    sll_real64                                         :: x(3),v(3)
    sll_int32                                          :: i_part
    sll_int32                                          :: i_v
    sll_real64, allocatable                            :: rdn(:)
    sll_real64                                         :: wi(1)
    sll_real64                                         :: rnd_no
    sll_int32                                          :: ip, i_gauss, counter
    sll_real64                                         :: delta(params%n_gaussians)

    n_rnds = 0
    if ( params%n_gaussians > 1 ) then
       n_rnds = 1
    end if

    do i_v=1,params%n_gaussians
       delta(i_v) = sum(params%delta(1:i_v))
    end do

    n_rnds = n_rnds+params%dims(1)+params%dims(2)
    allocate( rdn(1:params%dims(1)+params%dims(2)+1) )
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
          wi(1) = params%eval_x_density(x(1:params%dims(1)))*product(Lx)

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
          ! Set weight according to value of perturbation
          wi(1) = params%eval_x_density(x(1:params%dims(1)))*product(Lx)

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

    if( particle_group%species%q < 0._f64) then
       if( self%delta_perturb ) then
          call delta_function_perturbation( self, particle_group, xmin, Lx )
       end if
    end if

  end subroutine sample_particle_sampling_sym_1d2v



  !> Helper function for antithetic sampling in 1d2v
  subroutine sample_particle_sampling_sym2_1d2v( self, particle_group, params, xmin, Lx )
    type( sll_t_particle_sampling ), intent( inout ) :: self !< particle sampling object
    class(sll_c_particle_group_base), target, intent(inout)        :: particle_group !< particle group
    class( sll_c_distribution_params ), target, intent( inout )      :: params !< distribution parameter
    sll_real64,                        intent(in)               :: xmin(:) !< lower bound of the domain
    sll_real64,                        intent(in)               :: Lx(:) !< length of the domain.
    !local variables
    sll_int32 :: n_rnds
    sll_real64                                         :: x(3),v(3), xi(3)
    sll_int32                                          :: i_part
    sll_int32                                          :: i_v
    sll_real64, allocatable                            :: rdn(:)
    sll_real64                                         :: wi(1)
    sll_real64                                         :: rnd_no
    sll_int32                                          :: ip, i_gauss
    sll_real64                                         :: delta(params%n_gaussians)
    sll_real64                                         :: common_weight

    print*, 'new sampling'

    n_rnds = 0
    if ( params%n_gaussians > 1 ) then
       n_rnds = 1
    end if

    do i_v=1,params%n_gaussians
       delta(i_v) = sum(params%delta(1:i_v))
    end do

    n_rnds = n_rnds+params%dims(1)+params%dims(2)
    allocate( rdn(1:params%dims(1)+params%dims(2)+1) )
    rdn = 0.0_f64

    ! 1/Np in common weight
    common_weight =  particle_group%get_common_weight()
    call particle_group%set_common_weight &
         (common_weight/real(particle_group%n_total_particles, f64))

    if ( self%random_numbers == sll_p_random_numbers ) then
       call random_seed(put=self%random_seed)
    end if

    do i_part = 1, particle_group%n_particles
       ip = modulo(i_part, 2 )
       if ( ip == 1) then
          ! Generate Random or Sobol numbers on [0,1]
          select case( self%random_numbers )
          case( sll_p_sobol_numbers )
             call sll_s_i8_sobol( int(n_rnds,8), self%sobol_seed, rdn(1:n_rnds))
          case( sll_p_random_numbers )
             call random_number( rdn(1:n_rnds) )
          end select

          ! Transform rdn to the interval
          xi(1:params%dims(1)) = rdn(1:params%dims(1))
          x(1:params%dims(1)) = xmin + Lx * xi(1:params%dims(1))

          ! Set weight according to value of perturbation
          wi(1) = params%eval_x_density(x(1:params%dims(1)))*product(Lx)

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
       else
          x(1) = Lx(1) - x(1) + 2.0_f64*xmin(1)
          v(1) = -v(1) + 2.0_f64*params%v_mean(1,i_gauss)
          v(2) = -v(2) + 2.0_f64*params%v_mean(2,i_gauss)
       end if

       ! Copy the generated numbers to the particle
       call particle_group%set_x(i_part, x)
       call particle_group%set_v(i_part, v)
       ! Set weights.
       call particle_group%set_weights(i_part, &
            wi)

    end do

  end subroutine sample_particle_sampling_sym2_1d2v

  !> Helper function for antithetic sampling in 1d2v
  subroutine sample_particle_sampling_unix( self, particle_group, params, xmin, Lx )
    type( sll_t_particle_sampling ), intent( inout ) :: self !< particle sampling object
    class(sll_c_particle_group_base), target, intent(inout)        :: particle_group !< particle group
    class( sll_c_distribution_params ),  target,     intent( inout )      :: params !< distribution parameter
    sll_real64,                        intent(in)               :: xmin(:) !< lower bound of the domain
    sll_real64,                        intent(in)               :: Lx(:) !< length of the domain.
    !local variables
    sll_int32 :: n_rnds
    sll_real64                                         :: x(3),v(3), xi(3)
    sll_int32                                          :: i_part
    sll_int32                                          :: i_v
    sll_real64, allocatable                            :: rdn(:)
    sll_real64                                         :: wi(1)
    sll_real64                                         :: rnd_no
    sll_int32                                          :: i_gauss
    sll_real64                                         :: delta(params%n_gaussians)
    sll_real64                                         :: dx(params%dims(1))
    sll_real64                                         :: common_weight

    print*, 'new sampling 2'

    n_rnds = 0
    if ( params%n_gaussians > 1 ) then
       n_rnds = 1
    end if

    do i_v=1,params%n_gaussians
       delta(i_v) = sum(params%delta(1:i_v))
    end do

    n_rnds = n_rnds+params%dims(2)
    allocate( rdn(1:params%dims(1)+params%dims(2)+1) )
    rdn = 0.0_f64

    ! 1/Np in common weight
    common_weight =  particle_group%get_common_weight()
    call particle_group%set_common_weight &
         (common_weight/real(particle_group%n_total_particles, f64))

    if ( self%random_numbers == sll_p_random_numbers ) then
       call random_seed(put=self%random_seed)
    end if

    dx = Lx/real(particle_group%n_particles/2,f64)

    do i_part = 1, particle_group%n_particles/4

       ! Set x value
       x(1:params%dims(1)) = xmin + real(i_part-1, f64)*dx
       ! Set weight according to value of perturbation
       wi(1) = params%eval_x_density(x(1:params%dims(1)))*product(Lx)

       ! Generate first draw for v
       ! Generate Random or Sobol numbers on [0,1]
       select case( self%random_numbers )
       case( sll_p_sobol_numbers )
          call sll_s_i8_sobol( int(n_rnds,8), self%sobol_seed, rdn(1:n_rnds))
       case( sll_p_random_numbers )
          call random_number( rdn(1:n_rnds) )
       end select
       ! Maxwellian distribution of the temperature
       do i_v = 1,params%dims(2)
          call sll_s_normal_cdf_inv( rdn(i_v), 0.0_f64, 1.0_f64, &
               v(i_v))
       end do
       ! For multiple Gaussian, draw which one to take
       rnd_no = rdn(params%dims(2)+1)
       i_gauss = 1
       do while( rnd_no > delta(i_gauss) )
          i_gauss = i_gauss+1
       end do
       v(1:params%dims(2)) = v(1:params%dims(2)) * params%v_thermal(:,i_gauss) + params%v_mean(:,i_gauss)

       ! Copy the generated numbers to the particle
       call particle_group%set_x(2*i_part-1, x)
       call particle_group%set_v(2*i_part-1, v)
       ! Set weights.
       call particle_group%set_weights(2*i_part-1, &
            wi)

       ! Generate second random number
       ! Generate Random or Sobol numbers on [0,1]
       select case( self%random_numbers )
       case( sll_p_sobol_numbers )
          call sll_s_i8_sobol( int(n_rnds,8), self%sobol_seed, rdn(1:n_rnds))
       case( sll_p_random_numbers )
          call random_number( rdn(1:n_rnds) )
       end select
       ! Maxwellian distribution of the temperature
       do i_v = 1,params%dims(2)
          call sll_s_normal_cdf_inv( rdn(i_v), 0.0_f64, 1.0_f64, &
               v(i_v))
       end do
       ! For multiple Gaussian, draw which one to take
       rnd_no = rdn(params%dims(2)+1)
       i_gauss = 1
       do while( rnd_no > delta(i_gauss) )
          i_gauss = i_gauss+1
       end do
       v(1:params%dims(2)) = v(1:params%dims(2)) * params%v_thermal(:,i_gauss) + params%v_mean(:,i_gauss)

       ! Copy the generated numbers to the particle
       call particle_group%set_x(2*i_part, x)
       call particle_group%set_v(2*i_part, v)
       ! Set weights.
       call particle_group%set_weights(2*i_part, &
            wi)
    end do

    x = xmin + Lx*0.5_f64
    call particle_group%set_x(1, x)

    ! Add the negative particles
    do i_part = 1,particle_group%n_particles/2

       x = particle_group%get_x( i_part )
       v = -particle_group%get_v( i_part )

       x(1:params%dims(1)) = Lx - x(1:params%dims(1)) + 2.0_f64*xmin
       ! xi(1:params%dims(1)) = 1._f64 - xi(1:params%dims(1))
       wi(1) = params%eval_x_density(x(1:params%dims(1)))*product(Lx)

       ! Copy the generated numbers to the particle
       call particle_group%set_x(i_part+ particle_group%n_particles/2, x)
       call particle_group%set_v(i_part+ particle_group%n_particles/2, v)
       ! Set weights.
       call particle_group%set_weights(i_part+ particle_group%n_particles/2, &
            wi)

    end do

  end subroutine sample_particle_sampling_unix


  !> Helper function for antithetic sampling in 1d2v
  subroutine sample_particle_sampling_unix_1d2v( self, particle_group, params, xmin, Lx )
    type( sll_t_particle_sampling ), intent( inout ) :: self !< particle sampling object
    class(sll_c_particle_group_base), target, intent(inout)        :: particle_group !< particle group
    class( sll_c_distribution_params ),  target, intent( inout )      :: params !< distribution parameter
    sll_real64,                        intent(in)               :: xmin(:) !< lower bound of the domain
    sll_real64,                        intent(in)               :: Lx(:) !< length of the domain.
    !local variables
    sll_int32 :: n_rnds
    sll_real64                                         :: x(3),v(3), xi(3)
    sll_int32                                          :: i_part
    sll_int32                                          :: j
    sll_int32                                          :: i_v
    sll_real64, allocatable                            :: rdn(:)
    sll_real64                                         :: wi(1)
    sll_real64                                         :: rnd_no
    sll_int32                                          :: i_gauss
    sll_real64                                         :: delta(params%n_gaussians)
    sll_real64                                         :: dx(params%dims(1))
    sll_real64                                         :: common_weight

    print*, 'new sampling 2'

    n_rnds = 0
    if ( params%n_gaussians > 1 ) then
       n_rnds = 1
    end if

    do i_v=1,params%n_gaussians
       delta(i_v) = sum(params%delta(1:i_v))
    end do

    n_rnds = n_rnds+params%dims(2)
    allocate( rdn(1:params%dims(1)+params%dims(2)+1) )
    rdn = 0.0_f64

    ! 1/Np in common weight
    common_weight =  particle_group%get_common_weight()
    call particle_group%set_common_weight &
         (common_weight/real(particle_group%n_total_particles, f64))

    if ( self%random_numbers == sll_p_random_numbers ) then
       call random_seed(put=self%random_seed)
    end if

    dx = Lx/real(particle_group%n_particles/4,f64)

    do i_part = 1, particle_group%n_particles/8+1

       ! Set x value
       x(1:params%dims(1)) = xmin + real(i_part-1,f64)*dx
       ! Set weight according to value of perturbation
       wi(1) = params%eval_x_density(x(1:params%dims(1)))*product(Lx)

       ! Generate first draw for v
       ! Generate Random or Sobol numbers on [0,1]
       select case( self%random_numbers )
       case( sll_p_sobol_numbers )
          call sll_s_i8_sobol( int(n_rnds,8), self%sobol_seed, rdn(1:n_rnds))
       case( sll_p_random_numbers )
          call random_number( rdn(1:n_rnds) )
       end select
       ! Maxwellian distribution of the temperature
       do i_v = 1,params%dims(2)
          call sll_s_normal_cdf_inv( rdn(i_v), 0.0_f64, 1.0_f64, &
               v(i_v))
       end do
       ! For multiple Gaussian, draw which one to take
       rnd_no = rdn(params%dims(2)+1)
       i_gauss = 1
       do while( rnd_no > delta(i_gauss) )
          i_gauss = i_gauss+1
       end do
       v(1:params%dims(2)) = v(1:params%dims(2)) * params%v_thermal(:,i_gauss) + params%v_mean(:,i_gauss)

       ! Copy the generated numbers to the particle
       do j = 0,3
          call particle_group%set_x(4*i_part-j, x)
          ! Set weights.
          call particle_group%set_weights(4*i_part-j, &
               wi)
       end do

       call particle_group%set_v(4*i_part-3, v)
       call particle_group%set_v(4*i_part-2, [-v(1),v(2),v(3)])
       call particle_group%set_v(4*i_part-1, [v(1),-v(2),v(3)])
       call particle_group%set_v(4*i_part, [-v(1),-v(2),v(3)])


    end do

    ! Add the negative particles
    do i_part = 5,particle_group%n_particles/2

       x = particle_group%get_x( i_part )
       v = -particle_group%get_v( i_part )

       x(1:params%dims(1)) = Lx - x(1:params%dims(1)) + 2.0_f64*xmin
       !xi(1:params%dims(1)) = 1._f64 - xi(1:params%dims(1))
       wi(1) = params%eval_x_density(x(1:params%dims(1)))*product(Lx)

       ! Copy the generated numbers to the particle
       call particle_group%set_x(i_part+ particle_group%n_particles/2, x)
       call particle_group%set_v(i_part+ particle_group%n_particles/2, v)
       ! Set weights.
       call particle_group%set_weights(i_part+ particle_group%n_particles/2, &
            wi)

    end do

  end subroutine sample_particle_sampling_unix_1d2v



  !> Helper function for antithetic sampling in 3d3v
  subroutine sample_particle_sampling_sym_3d3v( self, particle_group, params, xmin, Lx )
    type( sll_t_particle_sampling ), intent( inout ) :: self !< particle sampling object
    class(sll_c_particle_group_base), target, intent(inout)        :: particle_group !< particle group
    class( sll_c_distribution_params ), target, intent( in )      :: params !< distribution parameter
    sll_real64,                        intent(in)               :: xmin(:) !< lower bound of the domain
    sll_real64,                        intent(in)               :: Lx(:) !< length of the domain.
    !local variables
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
    allocate( rdn(1:params%dims(1)+params%dims(2)+1) )
    rdn = 0.0_f64

    ! 1/Np in common weight
    call particle_group%set_common_weight &
         (1.0_f64/real(particle_group%n_total_particles, f64))

    if ( self%random_numbers == sll_p_random_numbers ) then
       call random_seed(put=self%random_seed)
    end if

    do i_part = 1, particle_group%n_particles
       ip = modulo(i_part, 64 )
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
          wi(1) = params%eval_x_density(x(1:params%dims(1)))*product(Lx)

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
       elseif ( modulo(i_part, 2) == 0 ) then
          v(3) = - v(3) + 2.0_f64*params%v_mean(3,i_gauss)
       elseif ( modulo(i_part, 4) == 3 ) then
          v(2) = - v(2) +  2.0_f64*params%v_mean(2,i_gauss)
       elseif ( modulo(i_part, 8) == 5 ) then
          v(1) = - v(1) +  2.0_f64*params%v_mean(1,i_gauss)
       elseif ( modulo(i_part, 16) == 9 ) then
          x(3) = Lx(3) - x(3) + 2.0_f64*xmin(3)
          ! Set weight according to value of perturbation
          wi(1) = params%eval_x_density(x(1:params%dims(1)))*product(Lx)
       elseif ( modulo(i_part, 32) == 17 ) then
          x(2) = Lx(2) - x(2) + 2.0_f64*xmin(2)
          ! Set weight according to value of perturbation
          wi(1) = params%eval_x_density(x(1:params%dims(1)))*product(Lx)
       elseif ( modulo(i_part, 64) == 33 ) then
          x(1) = Lx(1) - x(1) + 2.0_f64*xmin(1)
          ! Set weight according to value of perturbation
          wi(1) = params%eval_x_density(x(1:params%dims(1)))*product(Lx)
       end if

       ! Copy the generated numbers to the particle
       call particle_group%set_x(i_part, x)
       call particle_group%set_v(i_part, v)
       ! Set weights.
       call particle_group%set_weights(i_part, &
            wi)

    end do

    if( particle_group%species%q < 0._f64) then
       if( self%delta_perturb ) then
          call delta_function_perturbation( self, particle_group, xmin, Lx )
       end if
    end if

  end subroutine sample_particle_sampling_sym_3d3v


  !> Helper function for antithetic sampling in 1d2v
  subroutine sample_particle_sampling_sym_uni_3d3v( self, particle_group, params, xmin, Lx )
    type( sll_t_particle_sampling ), intent( inout ) :: self !< particle sampling object
    class(sll_c_particle_group_base), target, intent(inout)        :: particle_group !< particle group
    class( sll_c_distribution_params ),  target,     intent( inout )      :: params !< distribution parameter
    sll_real64,                        intent(in)               :: xmin(:) !< lower bound of the domain
    sll_real64,                        intent(in)               :: Lx(:) !< length of the domain.
    !local variables
    sll_int32 :: n_rnds
    sll_real64                                         :: x(3), v(3), xi(3)
    sll_int32                                          :: i_part, i_v
    sll_real64, allocatable                            :: rdn(:)
    sll_real64                                         :: wi(1)
    sll_real64                                         :: rnd_no
    sll_int32                                          :: ip, i_gauss
    sll_real64                                         :: v_min(3), Lv(3)

    n_rnds = 0
    if ( params%n_gaussians > 1 ) then
       n_rnds = 1
    end if


    n_rnds = n_rnds+params%dims(1)+params%dims(2)
    allocate( rdn(1:params%dims(1)+params%dims(2)+1) )
    rdn = 0.0_f64

    do i_v = 1, 3
       v_min(i_v) = minval(params%v_mean(i_v,:)) -3._f64 * maxval(params%v_thermal(i_v,:))
       Lv(i_v) = v_min(i_v) - maxval(params%v_mean(i_v,:)) +3._f64 * maxval(params%v_thermal(i_v,:))
    end do

    ! 1/Np in common weight
    call particle_group%set_common_weight &
         (1.0_f64/real(particle_group%n_total_particles, f64))

    if ( self%random_numbers == sll_p_random_numbers ) then
       call random_seed(put=self%random_seed)
    end if

    do i_part = 1, particle_group%n_particles
       ip = modulo(i_part, 64 )
       if ( ip == 1) then
          ! Generate Random or Sobol numbers on [0,1]
          select case( self%random_numbers )
          case( sll_p_sobol_numbers )
             call sll_s_i8_sobol( int(n_rnds,8), self%sobol_seed, rdn(1:n_rnds))
          case( sll_p_random_numbers )
             call random_number( rdn(1:n_rnds) )
          end select

          ! Transform rdn to the interval, sample x in [x_min, x_min+L_x]
          xi(1:params%dims(1)) = rdn(1:params%dims(1))
          x(1:params%dims(1)) = xmin + Lx * xi(1:params%dims(1))

          ! Sample v uniform in [v_min, v_min+Lv]
          v(1:params%dims(2)) = v_min + Lv * rdn(params%dims(1)+1:params%dims(1)+params%dims(2))

          ! Set weight according to value of perturbation
          wi(1) = params%eval_xv_density(x(1:params%dims(1)), v(1:params%dims(2)))*product(Lx)*product(Lv)
       elseif ( modulo(i_part, 2) == 0 ) then
          v(3) = - v(3) + 2.0_f64*params%v_mean(3,i_gauss)
       elseif ( modulo(i_part, 4) == 3 ) then
          v(2) = - v(2) +  2.0_f64*params%v_mean(2,i_gauss)
       elseif ( modulo(i_part, 8) == 5 ) then
          v(1) = - v(1) +  2.0_f64*params%v_mean(1,i_gauss)
       elseif ( modulo(i_part, 16) == 9 ) then
          x(3) = Lx(3) - x(3) + 2.0_f64*xmin(3)

          ! Set weight according to value of perturbation
          wi(1) = params%eval_x_density(x(1:params%dims(1)))*product(Lx)
       elseif ( modulo(i_part, 32) == 17 ) then
          x(2) = Lx(2) - x(2) + 2.0_f64*xmin(2)

          ! Set weight according to value of perturbation
          wi(1) = params%eval_x_density(x(1:params%dims(1)))*product(Lx)
       elseif ( modulo(i_part, 64) == 33 ) then
          x(1) = Lx(1) - x(1) + 2.0_f64*xmin(1)

          ! Set weight according to value of perturbation
          wi(1) = params%eval_x_density(x(1:params%dims(1)))*product(Lx)
       end if

       ! Copy the generated numbers to the particle
       call particle_group%set_x(i_part, x)
       call particle_group%set_v(i_part, v)
       ! Set weights.
       call particle_group%set_weights(i_part, &
            wi)

    end do

  end subroutine sample_particle_sampling_sym_uni_3d3v


  !> Helper function for pure sampling
  subroutine sample_particle_sampling_all_trafo( self,  particle_group, params, xmin, Lx, map )
    type(sll_t_particle_sampling),            intent( inout ) :: self !< particle sampling object
    class(sll_c_particle_group_base), target, intent( inout ) :: particle_group !< particle group
    class(sll_c_distribution_params), target, intent( inout ) :: params !< distribution parameter
    sll_real64,                               intent( in    ) :: xmin(:) !< lower bound of the domain
    sll_real64,                               intent( in    ) :: Lx(:) !< domain length           
    type(sll_t_mapping_3d),                   intent( inout ) :: map !< coordinate transformation
    !local variables
    sll_int32 :: n_rnds
    sll_real64                                         :: x(3), xi(3), v(3)
    sll_int32                                          :: i_part
    sll_int32                                          :: i_v, i_gauss
    sll_real64, allocatable                            :: rdn(:)
    sll_real64                                         :: wi(particle_group%n_weights)
    sll_real64                                         :: rnd_no
    sll_real64                                         :: delta(params%n_gaussians)
    sll_real64                                         :: w_total(1), w_local(1)
    sll_real64                                         :: common_weight
    sll_real64                                         :: r, Rx(3)

    w_local=0._f64
    w_total=0._f64
    n_rnds = 0
    if ( params%n_gaussians > 1 ) then
       n_rnds = 1
    end if

    do i_v=1,params%n_gaussians
       delta(i_v) = sum(params%delta(1:i_v))
    end do

    n_rnds = n_rnds+params%dims(1)+params%dims(2)
    allocate( rdn(1:params%dims(1)+params%dims(2)+1) )
    rdn = 0.0_f64

    ! 1/Np in common weight
    common_weight =  particle_group%get_common_weight()
    call particle_group%set_common_weight &
         (common_weight/real(particle_group%n_total_particles, f64))

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

       xi = 0._f64
       x = 0._f64
       v = 0._f64

       if( self%inverse ) then
          ! Transform rdn to the interval
          x(1:2) =  2._f64*map%params(2) * rdn(1:2) - map%params(2)
          x(3) = map%params(3) * rdn(3)
          r = sqrt(x(1)**2 +x(2)**2)
          do while( r > map%params(2) .or. r < map%params(1)  )
             select case( self%random_numbers )
             case( sll_p_sobol_numbers )
                call sll_s_i8_sobol( int(n_rnds,8), self%sobol_seed, rdn(1:n_rnds))
             case( sll_p_random_numbers )
                call random_number( rdn(1:n_rnds) )
             end select
             ! Transform rdn to the interval
             x(1:2) =  2._f64*map%params(2) * rdn(1:2) - map%params(2)
             x(3) = map%params(3) * rdn(3)
             ! inverse transformation of the uniformly sampled particles in physical coordinates
             r = sqrt(x(1)**2 +x(2)**2) 
          end do
          xi = map%get_xi(x)
          ! inverse transformation of the uniformly sampled particles in physical coordinates
          ! Set weight according to value of perturbation
          wi(1) = params%eval_x_density(x) * map%volume
       else
          xi(1:params%dims(1)) = rdn(1:params%dims(1))
          x = map%get_x(xi)
       end if
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
       if( self%xiprofile ) then
          select type(p=>params)
          type is(sll_t_params_cos_gaussian_screwpinch)
             if( particle_group%species%q > 0._f64) then
                v(1:p%dims(2)) = v(1:p%dims(2)) * sqrt(p%profile%T_i(xi(1))/particle_group%species%m) + p%v_mean(:,i_gauss)
             else
                v(1:p%dims(2)) = v(1:p%dims(2)) * sqrt(p%profile%T_e(xi(1))/particle_group%species%m) + p%v_mean(:,i_gauss)
             end if
             Rx = xi
             Rx(1) = x(1) + particle_group%species%m*v(2)/particle_group%species%q
             Rx(1) = (Rx(1) - xmin(1))/Lx(1)

             !special perturbation in gyrocoordinates
             wi(1) = p%eval_x_density(xi(1:p%dims(1)), v(1:p%dims(2)))*&
                  p%eval_v_density(v(1:p%dims(2)), Rx(1:p%dims(1)), particle_group%species%m)/p%eval_v_density(v(1:p%dims(2)), xi(1:p%dims(1)), particle_group%species%m)  *map%jacobian(xi)
          class default
             wi(1) = params%eval_x_density( xi(1:params%dims(1)) )*map%jacobian(xi)
          end select
       else
          v(1:params%dims(2)) = v(1:params%dims(2)) * params%v_thermal(:,i_gauss) + params%v_mean(:,i_gauss)
          ! Set weight according to value of perturbation
          wi(1) = params%eval_x_density( x(1:params%dims(1)) )*map%jacobian(xi)
       end if
       ! Copy the generated numbers to the particle
       call particle_group%set_x(i_part, xi)
       call particle_group%set_v(i_part, v)
       ! Set weights.
       call particle_group%set_weights(i_part, &
            wi)
       !w_local=w_local+wi(1)
    end do

    if( particle_group%species%q < 0._f64) then
       if( self%delta_perturb ) then
          call delta_function_perturbation( self, particle_group, [0._f64,0._f64,0._f64], [1._f64,1._f64,1._f64] )
       end if
    end if

  end subroutine sample_particle_sampling_all_trafo

  subroutine delta_function_perturbation( self,  particle_group, xmin, Lx )
    type(sll_t_particle_sampling), intent( inout )            :: self !< particle sampling object
    class(sll_c_particle_group_base), target, intent( inout ) :: particle_group
    sll_real64,                        intent(in)               :: xmin(:) !< lower bound of the domain
    sll_real64,                        intent(in)               :: Lx(:) !< length of the domain.
    !local variables
    sll_int32 :: n_rnds
    sll_real64                                         :: xi(3), vi(3)
    sll_int32                                          :: i_part, counter

    counter = 0
    do i_part = 1, particle_group%n_particles
       xi = particle_group%get_x(i_part)
       vi  = particle_group%get_v(i_part)

       xi = (xi - xmin)/Lx
       if( xi(1) > 0.5_f64-self%eps(1) .and. xi(1) < 0.5_f64+self%eps(1) )then
          if( xi(2) > 0.5_f64-self%eps(2) .and. xi(2) < 0.5_f64+self%eps(2) ) then
             if(xi(3) > 0.5_f64-self%eps(3) .and. xi(3) < 0.5_f64+self%eps(3) )then
                vi = vi + self%a
                counter = counter + 1
             end if
          end if
       end if
       call particle_group%set_v(i_part, vi)
    end do
    print*, 'delta_perturb, ', self%eps, self%a, 'shifted particles:', counter
  end subroutine delta_function_perturbation

  
  !> Helper function for antithetic sampling in 1d2v
  subroutine sample_particle_sampling_sym_1d2v_trafo( self, particle_group, params, Lx, map )
    type(sll_t_particle_sampling), intent( inout )            :: self !< particle sampling object
    class(sll_c_particle_group_base), target, intent( inout ) :: particle_group !< particle group
    class(sll_c_distribution_params), target, intent( inout ) :: params !< distribution parameter
    sll_real64, intent( in    )                               :: Lx(:) !< domain length   
    type(sll_t_mapping_3d),                   intent( inout ) :: map  !< coordinate transformation
    !local variables
    sll_int32                                          :: n_rnds
    sll_real64                                         :: x(3),v(3), xi(3)
    sll_int32                                          :: i_part
    sll_int32                                          :: i_v
    sll_real64, allocatable                            :: rdn(:)
    sll_real64                                         :: wi(1)
    sll_real64                                         :: rnd_no
    sll_int32                                          :: ip, i_gauss
    sll_real64                                         :: delta(params%n_gaussians)
    sll_real64                                         :: common_weight

    n_rnds = 0
    if ( params%n_gaussians > 1 ) then
       n_rnds = 1
    end if

    do i_v=1,params%n_gaussians
       delta(i_v) = sum(params%delta(1:i_v))
    end do

    n_rnds = n_rnds+params%dims(1)+params%dims(2)
    allocate( rdn(1:params%dims(1)+params%dims(2)+1) )
    rdn = 0.0_f64

    ! 1/Np in common weight
    common_weight =  particle_group%get_common_weight()
    call particle_group%set_common_weight &
         (common_weight/real(particle_group%n_total_particles, f64))

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

          xi = 0._f64
          x = 0._f64
          v = 0._f64

          ! Transform rdn to the interval
          xi(1:params%dims(1)) = rdn(1:params%dims(1))
          x(1) = map%get_x1(xi)
          ! Set weight according to value of perturbation
          wi(1) = params%eval_x_density(x(1:params%dims(1)))*abs(map%jacobian(xi))

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
          xi(1) = 1._f64 - xi(1)
          x = 0._f64
          x(1) = map%get_x1(xi)
          ! Set weight according to value of perturbation
          wi(1) = params%eval_x_density(x(1:params%dims(1)))*abs(map%jacobian(xi))
       elseif ( modulo(ip,2) == 0 ) then
          v(1) = -v(1) + 2.0_f64*params%v_mean(1,i_gauss)
       else          
          v(2) = -v(2) + 2.0_f64*params%v_mean(2,i_gauss)
       end if

       ! Copy the generated numbers to the particle
       call particle_group%set_x(i_part, xi)
       call particle_group%set_v(i_part, v)
       ! Set weights.
       call particle_group%set_weights(i_part, &
            wi)
    end do

  end subroutine sample_particle_sampling_sym_1d2v_trafo


  !> Helper function for antithetic sampling in 1d2v
  subroutine sample_particle_sampling_sym_3d3v_trafo( self, particle_group, params, xmin, Lx, map )
    type(sll_t_particle_sampling), intent( inout )            :: self !< particle sampling object
    class(sll_c_particle_group_base), target, intent( inout ) :: particle_group !< particle group
    class(sll_c_distribution_params), target, intent( inout ) :: params !< distribution parameter
    sll_real64,                               intent( in    ) :: xmin(:) !< lower bound of the domain
    sll_real64,                               intent( in    ) :: Lx(:) !< domain length  
    type(sll_t_mapping_3d),                   intent( inout ) :: map !< coordinate transformation 
    !local variables
    sll_int32 :: n_rnds
    sll_real64                                         :: x(3), xi(3), v(3)
    sll_int32                                          :: i_part
    sll_int32                                          :: i_v
    sll_real64, allocatable                            :: rdn(:)
    sll_real64                                         :: wi(particle_group%n_weights)
    sll_real64                                         :: rnd_no
    sll_int32                                          :: ip, i_gauss
    sll_real64                                         :: delta(params%n_gaussians)
    sll_real64                                         :: w_total(1), w_local(1)
    sll_real64                                         :: common_weight, weight
    sll_real64                                         :: r

    w_local=0._f64
    w_total=0._f64
    n_rnds = 0
    if ( params%n_gaussians > 1 ) then
       n_rnds = 1
    end if

    do i_v=1,params%n_gaussians
       delta(i_v) = sum(params%delta(1:i_v))
    end do

    n_rnds = n_rnds+params%dims(1)+params%dims(2)
    allocate( rdn(1:params%dims(1)+params%dims(2)+1) )
    rdn = 0.0_f64

    ! 1/Np in common weight
    common_weight =  particle_group%get_common_weight()
    call particle_group%set_common_weight &
         (common_weight/real(particle_group%n_total_particles, f64))

    if ( self%random_numbers == sll_p_random_numbers ) then
       call random_seed(put=self%random_seed)
    end if

    do i_part = 1, particle_group%n_particles
       ip = modulo(i_part, 4 )
       if ( ip == 1) then
          ! Generate Random or Sobol numbers on [0,1]
          select case( self%random_numbers )
          case( sll_p_sobol_numbers )
             call sll_s_i8_sobol( int(n_rnds,8), self%sobol_seed, rdn(1:n_rnds))
          case( sll_p_random_numbers )
             call random_number( rdn(1:n_rnds) )
          end select

          if( self%inverse ) then
             ! Transform rdn to the interval
             x(1:params%dims(1)-1) =  2._f64*map%params(2) * rdn(1:params%dims(1)-1) - map%params(2)
             x(params%dims(1)) = map%params(3) * rdn(params%dims(1))
             ! inverse transformation of the uniformly sampled particles in physical coordinates
             r = sqrt(x(1)**2 +x(2)**2)
             do while( r > map%params(2) .or. r < map%params(1)  )
                select case( self%random_numbers )
                case( sll_p_sobol_numbers )
                   call sll_s_i8_sobol( int(n_rnds,8), self%sobol_seed, rdn(1:n_rnds))
                case( sll_p_random_numbers )
                   call random_number( rdn(1:n_rnds) )
                end select
                ! Transform rdn to the interval
                x(1:params%dims(1)-1) =  2._f64*map%params(2) * rdn(1:params%dims(1)-1) - map%params(2)
                x(params%dims(1)) = map%params(3) * rdn(params%dims(1))
                ! inverse transformation of the uniformly sampled particles in physical coordinates
                r = sqrt(x(1)**2 +x(2)**2) 
             end do
             xi = map%get_xi(x)
             ! Set weight according to value of perturbation
             wi(1) =  params%eval_x_density(x) * map%volume
          else
             xi(1:params%dims(1)) = rdn(1:params%dims(1))
             x = map%get_x(xi)
          end if
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
          if( self%xiprofile ) then
             select type(p=>params)
             type is(sll_t_params_cos_gaussian_screwpinch)
                if( particle_group%species%q > 0._f64) then
                   v(1:p%dims(2)) = v(1:p%dims(2)) * sqrt(p%profile%T_i(xi(1))/particle_group%species%m) + p%v_mean(:,i_gauss)
                else
                   v(1:p%dims(2)) = v(1:p%dims(2)) * sqrt(p%profile%T_e(xi(1))/particle_group%species%m) + p%v_mean(:,i_gauss)
                end if
                class default
                   ! Set weight according to value of perturbation
                wi(1) = params%eval_x_density(xi)*map%jacobian(xi)
             end select
          else
             ! Set weight according to value of perturbation
             wi(1) = params%eval_x_density(x)*map%jacobian(xi)
             v(1:params%dims(2)) = v(1:params%dims(2)) * params%v_thermal(:,i_gauss) + params%v_mean(:,i_gauss)
          end if
       elseif ( modulo(i_part, 2) == 0 ) then
          v(3) = - v(3) + 2.0_f64*params%v_mean(3,i_gauss)
       elseif ( modulo(i_part, 4) == 3 ) then
          v(2) = - v(2) +  2.0_f64*params%v_mean(2,i_gauss)
          v(1) = - v(1) +  2.0_f64*params%v_mean(1,i_gauss)
       end if
       ! Copy the generated numbers to the particle
       call particle_group%set_x(i_part, xi)
       call particle_group%set_v(i_part, v)
       ! Set weights.
       call particle_group%set_weights(i_part, &
            wi)
       ! w_local=w_local+wi(1)
    end do

    if( particle_group%species%q < 0._f64) then
       if( self%delta_perturb ) then
          call delta_function_perturbation( self, particle_group, [0._f64,0._f64,0._f64], [1._f64,1._f64,1._f64] )
       end if
    end if

  end subroutine sample_particle_sampling_sym_3d3v_trafo

  subroutine sll_s_particle_sampling_randomized_weights(particle_group, alpha, seed)
    class(sll_c_particle_group_base), target, intent(inout)        :: particle_group
    sll_real64, intent(in) :: alpha
    sll_int32, intent(in), optional :: seed(:)

    sll_int32 :: i_part
    sll_real64 :: wi(particle_group%N_weights)
    sll_real64 :: pert, rnd(1)

    if (present(seed)) call random_seed(put=seed)


    do i_part = 1, particle_group%n_particles
       call random_number( rnd )
       wi = particle_group%get_weights(i_part)
       call sll_s_normal_cdf_inv( rnd(1), 0.0_f64, 1.0_f64, &
            pert)
       if (particle_group%n_weights > 2) then
          wi(3) = wi(3)-wi(1)
       end if
       wi(1) = wi(1) + alpha * pert
       if (particle_group%n_weights > 2) then
          wi(3) = wi(3)+wi(1)
       end if
       call particle_group%set_weights(i_part, wi)
    end do


  end subroutine sll_s_particle_sampling_randomized_weights


end module sll_m_particle_sampling
