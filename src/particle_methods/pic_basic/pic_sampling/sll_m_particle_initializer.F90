!> @ingroup pic_sampling
!> @author Katharina Kormann
!> @brief Particle initializer class with various functions to initialize a particle.
!> @details ...
module sll_m_particle_initializer

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_errors.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_particle_group_base, only: &
    sll_c_particle_group_base

  use sll_m_prob, only: &
    sll_s_normal_cdf_inv

  use sll_m_sobol, only: &
    sll_s_i8_sobol

  implicit none

  public :: &
    sll_s_particle_initialize_random_landau_1d2v, &
    sll_s_particle_initialize_random_landau_2d2v, &
    sll_s_particle_initialize_random_landau_symmetric_1d2v, &
    sll_s_particle_initialize_sobol_landau_1d2v, &
    sll_s_particle_initialize_sobol_landau_2d2v, &
    sll_s_particle_initialize_sobol_landau_symmetric_1d2v

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



contains


!> Initialize of a 1d2v particle group with Sobol pseudorandom numbers. Maxwellian distribution of V, cosine perturbation along x, equal weights.
  subroutine sll_s_particle_initialize_sobol_landau_1d2v(&
       particle_group, &
       landau_param, &
       xmin, &
       Lx, &
       thermal_velocity, &
       rnd_seed)
    class(sll_c_particle_group_base),  intent(inout)   :: particle_group
    sll_real64,                        intent(in)      :: landau_param(2) !< parameter defining the perturbation: landau_param(1)*cos(landau_param(2)*x1)
    sll_real64,                        intent(in)      :: xmin !< lower bound of the domain
    sll_real64,                        intent(in)      :: Lx !< length of the domain.
    sll_real64,                        intent(in)      :: thermal_velocity(2) !< Value of the thermal velocity along each dimension.
    sll_int64,                         intent(inout)   :: rnd_seed !< Random seed.

    sll_real64                                         :: x(3),v(3)
    sll_int32                                          :: i_part
    sll_int32                                          :: i_v
    sll_real64                                         :: rdn(3)
    sll_real64                                         :: wi
    

    x = 0.0_f64
    v = 0.0_f64

    ! 1/Np in common weight
    call particle_group%set_common_weight &
         (1.0_f64/real(particle_group%n_total_particles, f64))

    do i_part = 1, particle_group%n_particles
       ! Generate Sobol numbers on [0,1]
       call sll_s_i8_sobol( int(3,8), rnd_seed, rdn)

       ! Transform rdn to the interval
       x(1) = xmin + Lx * rdn(1)
       ! Landau perturbation for position 
       wi = (1.0_f64 + landau_param(1)*cos(landau_param(2)*x(1)))*Lx

       ! Maxwellian distribution of the temperature
       do i_v = 1,2
          call sll_s_normal_cdf_inv( rdn(i_v+1), 0.0_f64, 1.0_f64, &
               v(i_v))
          v(i_v) = v(i_v)*(thermal_velocity(i_v))
       end do

       ! Copy the generated numbers to the particle
       call particle_group%set_x(i_part, x)
       call particle_group%set_v(i_part, v)
       ! Set weights.
       call particle_group%set_weights(i_part, &
            [wi])
    end do


  end subroutine sll_s_particle_initialize_sobol_landau_1d2v

!> Initialize of a 1d2v particle group with Sobol pseudorandom numbers. Maxwellian distribution of V, cosine perturbation along x, equal weights. Symmetry around 0 for v and mirrored around center of domain in x.
  subroutine sll_s_particle_initialize_sobol_landau_symmetric_1d2v(&
       particle_group, &
       landau_param, &
       xmin, &
       Lx, &
       thermal_velocity, &
       rnd_seed)
    class(sll_c_particle_group_base),  intent(inout)            :: particle_group
    sll_real64,                        intent(in)               :: landau_param(2) !< parameter defining the perturbation: landau_param(1)*cos(landau_param(2)*x1)
    sll_real64,                        intent(in)               :: xmin !< lower bound of the domain
    sll_real64,                        intent(in)               :: Lx !< length of the domain.
    sll_real64,                        intent(in)               :: thermal_velocity(2) !< Value of the thermal velocity along each dimension.
    sll_int64,                         intent(inout)            :: rnd_seed !< Random seed.

    sll_real64                                                  :: x(3),v(3)
    sll_int32                                                   :: i_part
    sll_int32                                                   :: i_v
    sll_real64                                                  :: rdn(3)
    sll_real64                                                  :: wi
    sll_int32                                                   :: ip
    

    x = 0.0_f64
    v = 0.0_f64

    if (modulo(particle_group%n_particles,8) .NE. 0) then
       SLL_ERROR('sll_s_particle_initialize_random_landau_symmetric_1d2v', 'particle number not multiple of 8')
    end if

    ! 1/Np in common weight
    call particle_group%set_common_weight &
         (1.0_f64/real(particle_group%n_total_particles, f64))

    do i_part = 1, particle_group%n_particles
       ip = modulo(i_part, 8)
       if ( ip == 1) then
          ! Generate Sobol numbers on [0,1]
          call sll_s_i8_sobol( int(3,8), rnd_seed, rdn)

          ! Transform rdn to the interval
          x(1) = xmin + Lx * rdn(1)
          ! Landau perturbation for position 
          wi = (1.0_f64 + landau_param(1)*cos(landau_param(2)*x(1)))*Lx

          ! Maxwellian distribution of the temperature
          do i_v = 1,2
             call sll_s_normal_cdf_inv( rdn(i_v+1), 0.0_f64, 1.0_f64, &
                  v(i_v))
             v(i_v) = v(i_v)*(thermal_velocity(i_v))
          end do
       elseif ( ip == 5) then
          x(1) = Lx - x(1) + 2.0_f64*xmin
       elseif ( modulo(ip,2) == 0 ) then
          v(2) = - v(2)
       else
          v(1) = -v(1)
       end if

       ! Copy the generated numbers to the particle
       call particle_group%set_x(i_part, x)
       call particle_group%set_v(i_part, v)
       ! Set weights.
       call particle_group%set_weights(i_part, &
            [wi])
    end do


  end subroutine sll_s_particle_initialize_sobol_landau_symmetric_1d2v

!> Initialize of a 1d2v particle group with Sobol pseudorandom numbers. Maxwellian distribution of V, cosine perturbation along x, equal weights.
  subroutine sll_s_particle_initialize_random_landau_1d2v(&
       particle_group, &
       landau_param, &
       xmin, &
       Lx, &
       thermal_velocity, &
       rnd_seed)
    class(sll_c_particle_group_base),  intent(inout)            :: particle_group
    sll_real64,                        intent(in)               :: landau_param(2) !< parameter defining the perturbation: landau_param(1)*cos(landau_param(2)*x1)
    sll_real64,                        intent(in)               :: xmin !< lower bound of the domain
    sll_real64,                        intent(in)               :: Lx !< length of the domain.
    sll_real64,                        intent(in)               :: thermal_velocity(2) !< Value of the thermal velocity along each dimension.
    sll_int32,                         intent(inout)            :: rnd_seed(:) !< Random seed.

    sll_real64                                                  :: x(3),v(3)
    sll_real64                                                  :: rnd_no
    sll_int32                                                   :: i_part
    sll_int32                                                   :: i_v
    sll_real64                                                  :: wi
    
    ! Set random seed.
    call random_seed (put=rnd_seed)

    x = 0.0_f64
    v = 0.0_f64

    ! 1/Np in common weight
    call particle_group%set_common_weight &
         (1.0_f64/real(particle_group%n_total_particles, f64))

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
          call sll_s_normal_cdf_inv( rnd_no, 0.0_f64, 1.0_f64, &
               v(i_v))
          v(i_v) = v(i_v)*(thermal_velocity(i_v))
       end do

       ! Copy the generated numbers to the particle
       call particle_group%set_x(i_part, x)
       call particle_group%set_v(i_part, v)
       ! Set weights.
       call particle_group%set_weights(i_part, &
            [wi])
    end do


  end subroutine sll_s_particle_initialize_random_landau_1d2v


!> Initialize of a 1d2v particle group with Sobol pseudorandom numbers. Maxwellian distribution of V, cosine perturbation along x, equal weights.  Symmetry around 0 for v and mirrored around center of domain in x.
  subroutine sll_s_particle_initialize_random_landau_symmetric_1d2v(&
       particle_group, &
       landau_param, &
       xmin, &
       Lx, &
       thermal_velocity, &
       rnd_seed)
    class(sll_c_particle_group_base),  intent(inout)            :: particle_group
    sll_real64,                        intent(in)               :: landau_param(2) !< parameter defining the perturbation: landau_param(1)*cos(landau_param(2)*x1)
    sll_real64,                        intent(in)               :: xmin !< lower bound of the domain
    sll_real64,                        intent(in)               :: Lx !< length of the domain.
    sll_real64,                        intent(in)               :: thermal_velocity(2) !< Value of the thermal velocity along each dimension.
    sll_int32,                         intent(inout)            :: rnd_seed(:) !< Random seed.

    sll_real64                                                  :: x(3),v(3), wi(1)
    sll_real64                                                  :: rnd_no
    sll_int32                                                   :: i_part
    sll_int32                                                   :: i_v
    sll_int32                                                   :: ip
    
    ! Set random seed.
    call random_seed (put=rnd_seed)

    x = 0.0_f64
    v = 0.0_f64

    if (modulo(particle_group%n_particles,8) .NE. 0) then
       SLL_ERROR('sll_s_particle_initialize_random_landau_symmetric_1d2v', 'particle number not multiple of 8.')
    end if


    ! 1/Np in common weight
    call particle_group%set_common_weight &
         (1.0_f64/real(particle_group%n_total_particles, f64))

    do i_part = 1, particle_group%n_particles
       ip = modulo(i_part, 8)
       if ( ip == 1) then
          ! Generate random number on [0,1]
          call random_number(rnd_no)
          ! Transform rdn to the interval
          x(1) = xmin + Lx * rnd_no
          ! Landau perturbation for position 
          wi = (1.0_f64 + landau_param(1)*cos(landau_param(2)*x(1)))*Lx

          ! Maxwellian distribution of the temperature
          do i_v = 1,2
             call random_number(rnd_no)
             call sll_s_normal_cdf_inv( rnd_no, 0.0_f64, 1.0_f64, &
                  v(i_v))
             v(i_v) = v(i_v)*(thermal_velocity(i_v))
          end do
       elseif ( ip == 5) then
          x(1) = Lx - x(1) + 2.0_f64*xmin
       elseif ( modulo(ip,2) == 0 ) then
          v(2) = - v(2)
       else
          v(1) = -v(1)
       end if

       ! Copy the generated numbers to the particle
       call particle_group%set_x(i_part, x)
       call particle_group%set_v(i_part, v)
       ! Set weights.
       call particle_group%set_weights(i_part, &
            [wi])
    end do



  end subroutine sll_s_particle_initialize_random_landau_symmetric_1d2v



  !> Initialize of a 2d2v particle group with Sobol pseudorandom numbers. Maxwellian distribution of V, cosine perturbation along x, equal weights.
  subroutine sll_s_particle_initialize_sobol_landau_2d2v(&
       particle_group, &
       landau_param, &
       xmin, &
       Lx, &
       thermal_velocity, &
       rnd_seed)
    class(sll_c_particle_group_base),  intent(inout)            :: particle_group
    sll_real64,                        intent(in)               :: landau_param(2) !< parameter defining the perturbation: landau_param(1)*cos(landau_param(2)*x1)
    sll_real64,                        intent(in)               :: xmin(2) !< lower bound of the domain
    sll_real64,                        intent(in)               :: Lx(2) !< length of the domain.
    sll_real64,                        intent(in)               :: thermal_velocity(2) !< Value of the thermal velocity along each dimension.
    sll_int64,                         intent(inout)            :: rnd_seed !< Random seed.
    

    sll_real64                                                  :: x(3),v(3)
    sll_int32                                                   :: i_part
    sll_int32                                                   :: i_v
    sll_real64                                                  :: rdn(4)
    sll_real64                                                  :: wi(3)

    x = 0.0_f64
    v = 0.0_f64

    ! 1/Np in common weight
    call particle_group%set_common_weight &
         (1.0_f64/real(particle_group%n_total_particles, f64))


    do i_part = 1, particle_group%n_particles
       ! Generate Sobol numbers on [0,1]
       call sll_s_i8_sobol( int(4,8), rnd_seed, rdn)

       ! Transform rdn to the interval
       x(1:2) = xmin + Lx * rdn(1:2)
       ! Landau perturbation for position 
       wi(1) = (1.0_f64 + landau_param(1)*cos(landau_param(2)*x(1)))*Lx(1)*Lx(2)
          

       ! Maxwellian distribution of the temperature
       do i_v = 1,2
          call sll_s_normal_cdf_inv( rdn(i_v+2), 0.0_f64, 1.0_f64, &
               v(i_v))
          v(i_v) = v(i_v)*(thermal_velocity(i_v))
       end do
       
       ! If control variate, set also g0 and df weight
       if (particle_group%n_weights == 3) then
          wi(2) = 1.0_f64/product(Lx)*&
               exp(-0.5_f64*((v(1)/thermal_velocity(1))**2+&
               (v(2)/thermal_velocity(2))**2))
          wi(3) = wi(1) - exp(-0.5_f64*((v(1)/thermal_velocity(1))**2+&
               (v(2)/thermal_velocity(2))**2))/wi(2)
       end if

       ! Copy the generated numbers to the particle
       call particle_group%set_x(i_part, x)
       call particle_group%set_v(i_part, v)
       ! Set weights.
       call particle_group%set_weights(i_part, &
            wi(1:particle_group%n_weights))
    end do


  end subroutine sll_s_particle_initialize_sobol_landau_2d2v


!> Initialize of a 1d2v particle group with Sobol pseudorandom numbers. Maxwellian distribution of V, cosine perturbation along x, equal weights.
  subroutine sll_s_particle_initialize_random_landau_2d2v(&
       particle_group, &
       landau_param, &
       xmin, &
       Lx, &
       thermal_velocity, &
       rnd_seed)
    class(sll_c_particle_group_base),  intent(inout)            :: particle_group
    sll_real64,                        intent(in)               :: landau_param(2) !< parameter defining the perturbation: landau_param(1)*cos(landau_param(2)*x1)
    sll_real64,                        intent(in)               :: xmin(2) !< lower bound of the domain
    sll_real64,                        intent(in)               :: Lx(2) !< length of the domain.
    sll_real64,                        intent(in)               :: thermal_velocity(2) !< Value of the thermal velocity along each dimension.
    sll_int32,                         intent(inout)            :: rnd_seed(:) !< Random seed.

    sll_real64                                                  :: x(3),v(3)
    sll_real64                                                  :: rnd_no
    sll_int32                                                   :: i_part
    sll_int32                                                   :: i_v
    sll_real64                                                  :: wi
    
    ! Set random seed.
    call random_seed (put=rnd_seed)

    x = 0.0_f64
    v = 0.0_f64


    ! 1/Np in common weight
    call particle_group%set_common_weight &
         (1.0_f64/real(particle_group%n_total_particles, f64))

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
          call sll_s_normal_cdf_inv( rnd_no, 0.0_f64, 1.0_f64, &
               v(i_v))
          v(i_v) = v(i_v)*(thermal_velocity(i_v))
       end do

       ! Copy the generated numbers to the particle
       call particle_group%set_x(i_part, x)
       call particle_group%set_v(i_part, v)
       ! Set weights.
       call particle_group%set_weights(i_part, &
            [wi])
    end do


  end subroutine sll_s_particle_initialize_random_landau_2d2v

end module sll_m_particle_initializer
