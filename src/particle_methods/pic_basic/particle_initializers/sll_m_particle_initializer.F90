!> @ingroup particle_initializers
!> @author Katharina Kormann
!> @brief Particle initializer class with various functions to initialize a particle.
!> @details ...
module sll_m_particle_initializer

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_gaussian, only: &
    gaussian_deviate_2d

  use sll_m_particle_group_base, only: &
    sll_particle_group_base

  use sll_m_prob, only: &
    normal_cdf_inv

  use sll_m_sobol, only: &
    i8_sobol

  implicit none

  public :: &
    sll_particle_initialize_random_landau_1d2v, &
    sll_particle_initialize_random_landau_2d2v, &
    sll_particle_initialize_sobol_landau_1d2v, &
    sll_particle_initialize_sobol_landau_2d2v

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! Different ansatz. Use or remove
!!$  type :: sll_particle_initializer_2d2v
!!$     
!!$
!!$   contains
!!$     procedure :: draw_product_density => draw_product_density_2x2v
!!$     
!!$  end type sll_particle_initializer_2x2v
 

contains
! Alternative versions. Decide if use or remove.
!!$  !---------------------------------------------------------------------------!
!!$  
!!$  subroutine draw_product_density_2d2v (this, particle_group, random_seed)
!!$    class(sll_particle_initializer_2d2v), intent(inout)         :: this
!!$    class(sll_particle_group_base),  intent(inout)              :: particle_group
!!$    !sll_real64, intent(in)                                      :: thermal_velocity
!!$    sll_int32,  intent(in)                                      :: random_seed
!!$
!!$
!!$    !print*, this%densityx1(0.0)
!!$    
!!$
!!$  end subroutine draw_product_density_2x2v

  !> Initialization to initialize a 2x2v particle group. Maxwellian distribution of V, cosine perturbation along x1, equal weights.
!!$  subroutine sll_particle_initialize_landaux1_2d2v(&
!!$       particle_group, &
!!$       landau_param, &
!!$       xmin, &
!!$       Lx, &
!!$       thermal_velocity, &
!!$       rnd_seed)
!!$    class(sll_particle_group_base),  intent(inout)              :: particle_group
!!$    sll_real64, intent(in)                                      :: landau_param(2) !< parameter defining the perturbation: landau_param(1)*cos(landau_param(2)*x1)
!!$    sll_real64, intent(in)                                      :: xmin(2) !< lower bound of the domain
!!$    sll_real64, intent(in)                                      :: Lx(2) !< length of the domain.
!!$    sll_real64, intent(in)                                      :: thermal_velocity(2) !< Value of the thermal velocity along each dimension.
!!$    sll_int32,  intent(in)                                      :: rnd_seed(:) !< Random seed.
!!$
!!$    sll_real64                                                  :: x(3),v(3)
!!$    sll_real64                                                  :: rnd_no
!!$    sll_int32                                                   :: i_part
!!$    
!!$
!!$    x = 0.0_f64
!!$    v = 0.0_f64
!!$
!!$    ! Set random seed.
!!$    call random_seed (put=rnd_seed)
!!$
!!$    ! Loop over particles to set coordinates and weights.
!!$    i_part = 0
!!$    do while( i_part < particle_group%n_particles)
!!$       ! Rejection sampling for x(1)
!!$       call random_number(rnd_no)
!!$       x(1) = Lx(1) * rnd_no + xmin(1)
!!$       call random_number(rnd_no)
!!$       rnd_no = (1.0_f64 + landau_param(1)) * rnd_no
!!$       if ( (1.0_f64 + landau_param(1) * cos(landau_param(2)*x(1))) >= rnd_no ) then  
!!$          ! Draw x(2) from uniform distribution
!!$          i_part = i_part + 1
!!$          call random_number(rnd_no)
!!$          x(2) = Lx(2) * rnd_no + xmin(2)
!!$
!!$          ! Draw velocities from 2D Gauss distribution
!!$          call gaussian_deviate_2D(v(1:2))
!!$          v(1:2) = v(1:2) * thermal_velocity
!!$
!!$          call particle_group%set_x(i_part, x)
!!$          call particle_group%set_v(i_part, v)
!!$          ! Set weights.
!!$          call particle_group%set_weights(i_part, &
!!$               [Lx(1)*Lx(2)/real(particle_group%n_total_particles, f64)])
!!$       end if
!!$    end do
!!$
!!$    
!!$
!!$  end subroutine sll_particle_initialize_landaux1_2d2v


!!$!> Initialization to initialize a 1d2v particle group. Maxwellian distribution of V, cosine perturbation along x, equal weights.
!!$  subroutine sll_particle_initialize_landau_1d2v(&
!!$       particle_group, &
!!$       landau_param, &
!!$       xmin, &
!!$       Lx, &
!!$       thermal_velocity, &
!!$       rnd_seed)
!!$    class(sll_particle_group_base),  intent(inout)              :: particle_group
!!$    sll_real64, intent(in)                                      :: landau_param(2) !< parameter defining the perturbation: landau_param(1)*cos(landau_param(2)*x1)
!!$    sll_real64, intent(in)                                      :: xmin !< lower bound of the domain
!!$    sll_real64, intent(in)                                      :: Lx !< length of the domain.
!!$    sll_real64, intent(in)                                      :: thermal_velocity(2) !< Value of the thermal velocity along each dimension.
!!$    sll_int32,  intent(in)                                      :: rnd_seed(:) !< Random seed.
!!$
!!$    sll_real64                                                  :: x(3),v(3)
!!$    sll_real64                                                  :: rnd_no
!!$    sll_int32                                                   :: i_part
!!$    
!!$
!!$    x = 0.0_f64
!!$    v = 0.0_f64
!!$
!!$    ! Set random seed.
!!$    call random_seed (put=rnd_seed)
!!$
!!$    ! Loop over particles to set coordinates and weights.
!!$    i_part = 0
!!$    do while( i_part < particle_group%n_particles)
!!$       ! Rejection sampling for x(1)
!!$       call random_number(rnd_no)
!!$       x(1) = Lx * rnd_no + xmin
!!$       call random_number(rnd_no)
!!$       rnd_no = (1.0_f64 + landau_param(1)) * rnd_no
!!$       if ( (1.0_f64 + landau_param(1) * cos(landau_param(2)*x(1))) >= rnd_no ) then  
!!$          ! Draw x(2) from uniform distribution
!!$          i_part = i_part + 1
!!$
!!$          ! Draw velocities from 2D Gauss distribution
!!$          call gaussian_deviate_2D(v(1:2))
!!$          v(1:2) = v(1:2) * thermal_velocity
!!$
!!$          call particle_group%set_x(i_part, x)
!!$          call particle_group%set_v(i_part, v)
!!$          ! Set weights.
!!$          call particle_group%set_weights(i_part, &
!!$               [Lx/real(particle_group%n_total_particles, f64)])
!!$       end if
!!$    end do
!!$
!!$    
!!$
!!$  end subroutine sll_particle_initialize_landau_1d2v


!> Initialize of a 1d2v particle group with Sobol pseudorandom numbers. Maxwellian distribution of V, cosine perturbation along x, equal weights.
  subroutine sll_particle_initialize_sobol_landau_1d2v(&
       particle_group, &
       landau_param, &
       xmin, &
       Lx, &
       thermal_velocity, &
       rnd_seed)
    class(sll_particle_group_base),  intent(inout)              :: particle_group
    sll_real64, intent(in)                                      :: landau_param(2) !< parameter defining the perturbation: landau_param(1)*cos(landau_param(2)*x1)
    sll_real64, intent(in)                                      :: xmin !< lower bound of the domain
    sll_real64, intent(in)                                      :: Lx !< length of the domain.
    sll_real64, intent(in)                                      :: thermal_velocity(2) !< Value of the thermal velocity along each dimension.
    sll_int64,  intent(inout)                                   :: rnd_seed !< Random seed.

    sll_real64                                                  :: x(3),v(3)
    sll_int32                                                   :: i_part
    sll_int32                                                   :: i_v
    sll_real64                                                  :: rdn(3)
    sll_real64                                                  :: wi
    

    x = 0.0_f64
    v = 0.0_f64

    !rdn(1) = 0.1_f64
    !call normal_cdf_inv(rdn(1), 0.0_f64, 1.0_f64, v(1))
    !print*, 'cdf', v(1)


    do i_part = 1, particle_group%n_particles
       ! Generate Sobol numbers on [0,1]
       call i8_sobol( int(3,8), rnd_seed, rdn)

       ! Transform rdn to the interval
       x(1) = xmin + Lx * rdn(1)
       ! Landau perturbation for position 
       wi = (1.0_f64 + landau_param(1)*cos(landau_param(2)*x(1)))*Lx

       ! Maxwellian distribution of the temperature
       do i_v = 1,2
          call normal_cdf_inv( rdn(i_v+1), 0.0_f64, 1.0_f64, &
               v(i_v))
          v(i_v) = v(i_v)*(thermal_velocity(i_v))
       end do

       ! Copy the generated numbers to the particle
       call particle_group%set_x(i_part, x)
       call particle_group%set_v(i_part, v)
       ! Set weights.
       call particle_group%set_weights(i_part, &
            [wi/real(particle_group%n_total_particles, f64)])
    end do


  end subroutine sll_particle_initialize_sobol_landau_1d2v

!> Initialize of a 1d2v particle group with Sobol pseudorandom numbers. Maxwellian distribution of V, cosine perturbation along x, equal weights.
  subroutine sll_particle_initialize_random_landau_1d2v(&
       particle_group, &
       landau_param, &
       xmin, &
       Lx, &
       thermal_velocity, &
       rnd_seed)
    class(sll_particle_group_base),  intent(inout)              :: particle_group
    sll_real64, intent(in)                                      :: landau_param(2) !< parameter defining the perturbation: landau_param(1)*cos(landau_param(2)*x1)
    sll_real64, intent(in)                                      :: xmin !< lower bound of the domain
    sll_real64, intent(in)                                      :: Lx !< length of the domain.
    sll_real64, intent(in)                                      :: thermal_velocity(2) !< Value of the thermal velocity along each dimension.
    sll_int32,  intent(inout)                                   :: rnd_seed(:) !< Random seed.

    sll_real64                                                  :: x(3),v(3)
    sll_real64                                                  :: rnd_no
    sll_int32                                                   :: i_part
    sll_int32                                                   :: i_v
    sll_real64                                                  :: wi
    
    ! Set random seed.
    call random_seed (put=rnd_seed)

    x = 0.0_f64
    v = 0.0_f64


    do i_part = 1, particle_group%n_particles
       ! Generate Sobol numbers on [0,1]
       call random_number(rnd_no)
       ! Transform rdn to the interval
       x(1) = xmin + Lx * rnd_no
       ! Landau perturbation for position 
       wi = (1.0_f64 + landau_param(1)*cos(landau_param(2)*x(1)))*Lx

       ! Maxwellian distribution of the temperature
       do i_v = 1,2
          call random_number(rnd_no)
          call normal_cdf_inv( rnd_no, 0.0_f64, 1.0_f64, &
               v(i_v))
          v(i_v) = v(i_v)*(thermal_velocity(i_v))
       end do

       ! Copy the generated numbers to the particle
       call particle_group%set_x(i_part, x)
       call particle_group%set_v(i_part, v)
       ! Set weights.
       call particle_group%set_weights(i_part, &
            [wi/real(particle_group%n_total_particles, f64)])
    end do


  end subroutine sll_particle_initialize_random_landau_1d2v



!> Initialize of a 2d2v particle group with Sobol pseudorandom numbers. Maxwellian distribution of V, cosine perturbation along x, equal weights.
  subroutine sll_particle_initialize_sobol_landau_2d2v(&
       particle_group, &
       landau_param, &
       xmin, &
       Lx, &
       thermal_velocity, &
       rnd_seed)
    class(sll_particle_group_base),  intent(inout)              :: particle_group
    sll_real64, intent(in)                                      :: landau_param(2) !< parameter defining the perturbation: landau_param(1)*cos(landau_param(2)*x1)
    sll_real64, intent(in)                                      :: xmin(2) !< lower bound of the domain
    sll_real64, intent(in)                                      :: Lx(2) !< length of the domain.
    sll_real64, intent(in)                                      :: thermal_velocity(2) !< Value of the thermal velocity along each dimension.
    sll_int64,  intent(inout)                                   :: rnd_seed !< Random seed.

    sll_real64                                                  :: x(3),v(3)
    sll_int32                                                   :: i_part
    sll_int32                                                   :: i_v
    sll_real64                                                  :: rdn(4)
    sll_real64                                                  :: wi
    

    x = 0.0_f64
    v = 0.0_f64


    do i_part = 1, particle_group%n_particles
       ! Generate Sobol numbers on [0,1]
       call i8_sobol( int(4,8), rnd_seed, rdn)

       ! Transform rdn to the interval
       x(1:2) = xmin + Lx * rdn(1:2)
       ! Landau perturbation for position 
       wi = (1.0_f64 + landau_param(1)*cos(landau_param(2)*x(1)))*Lx(1)*Lx(2)


       ! Maxwellian distribution of the temperature
       do i_v = 1,2
          call normal_cdf_inv( rdn(i_v+2), 0.0_f64, 1.0_f64, &
               v(i_v))
          v(i_v) = v(i_v)*(thermal_velocity(i_v))
       end do

       ! Copy the generated numbers to the particle
       call particle_group%set_x(i_part, x)
       call particle_group%set_v(i_part, v)
       ! Set weights.
       call particle_group%set_weights(i_part, &
            [wi/real(particle_group%n_total_particles, f64)])
    end do


  end subroutine sll_particle_initialize_sobol_landau_2d2v


!> Initialize of a 1d2v particle group with Sobol pseudorandom numbers. Maxwellian distribution of V, cosine perturbation along x, equal weights.
  subroutine sll_particle_initialize_random_landau_2d2v(&
       particle_group, &
       landau_param, &
       xmin, &
       Lx, &
       thermal_velocity, &
       rnd_seed)
    class(sll_particle_group_base),  intent(inout)              :: particle_group
    sll_real64, intent(in)                                      :: landau_param(2) !< parameter defining the perturbation: landau_param(1)*cos(landau_param(2)*x1)
    sll_real64, intent(in)                                      :: xmin(2) !< lower bound of the domain
    sll_real64, intent(in)                                      :: Lx(2) !< length of the domain.
    sll_real64, intent(in)                                      :: thermal_velocity(2) !< Value of the thermal velocity along each dimension.
    sll_int32,  intent(inout)                                   :: rnd_seed(:) !< Random seed.

    sll_real64                                                  :: x(3),v(3)
    sll_real64                                                  :: rnd_no
    sll_int32                                                   :: i_part
    sll_int32                                                   :: i_v
    sll_real64                                                  :: wi
    
    ! Set random seed.
    call random_seed (put=rnd_seed)

    x = 0.0_f64
    v = 0.0_f64


    do i_part = 1, particle_group%n_particles
       ! Generate Sobol numbers on [0,1]
       call random_number(rnd_no)
       ! Transform rdn to the interval
       x(1) = xmin(1) + Lx(1) * rnd_no

       call random_number(rnd_no)
       x(2) = xmin(2) + Lx(2) * rnd_no
       ! Landau perturbation for position 
       wi = (1.0_f64 + landau_param(1)*cos(landau_param(2)*x(1)))*Lx(1)*Lx(2)

       ! Maxwellian distribution of the temperature
       do i_v = 1,2
          call random_number(rnd_no)
          call normal_cdf_inv( rnd_no, 0.0_f64, 1.0_f64, &
               v(i_v))
          v(i_v) = v(i_v)*(thermal_velocity(i_v))
       end do

       ! Copy the generated numbers to the particle
       call particle_group%set_x(i_part, x)
       call particle_group%set_v(i_part, v)
       ! Set weights.
       call particle_group%set_weights(i_part, &
            [wi/real(particle_group%n_total_particles, f64)])
    end do


  end subroutine sll_particle_initialize_random_landau_2d2v

end module sll_m_particle_initializer
