program test_control_variate
#include "sll_working_precision.h"
  use sll_m_constants, only : &
       sll_p_pi

  use sll_m_control_variate, only : &
       sll_t_control_variate, &
       sll_f_control_variate
  use sll_m_particle_group_base, only : &
       sll_c_particle_group_base
  use sll_m_particle_group_2d2v, only : &
       sll_s_new_particle_group_2d2v

  class(sll_t_control_variate), pointer :: control_variate
  sll_real64, pointer :: control_variate_parameter(:)
  class(sll_c_particle_group_base), pointer :: particle_group
  sll_int32 :: n_particles
  sll_real64 :: x_vec(4,2)
  sll_real64 :: v_vec(4,2)
  sll_real64 :: w_vec(4,3)
  sll_real64 :: xi(3), wi(3)

  logical :: passed 
  passed = .TRUE.
  
  n_particles = 4

  ! Initialize the control variate
  allocate(control_variate_parameter(2))
  control_variate_parameter = [1.0_f64, 1.0_f64]
  allocate(control_variate)
  call control_variate%initialize(control_variate_equi, &
       control_variate_parameter)

  ! Initialize a particle group
  x_vec(:,1) = [0.1_f64, 0.65_f64, 0.7_f64, 1.5_f64] ! Particle positions
  x_vec(:,2) = [0.0_f64, 0.0_f64, 0.0_f64, 0.0_f64] ! Particle positions
  v_vec(:,1) = [1.5_f64, 0.0_f64, 0.0_f64, 0.0_f64]
  v_vec(:,2) = [0.0_f64, 0.5_f64, 0.0_f64, 0.0_f64]

  w_vec(:,1) = 1.0_f64
  w_vec(:,2) = [0.1_f64, 0.01_f64, 0.01_f64, 1.0_f64]
  w_vec(:,3) = w_vec(:,1) - &
       [5.1670044967061561D-002, 0.14045374430962521D0, &
       0.15915494309189535D0, 0.15915494309189535D0]/w_vec(:,2)

  ! We need to initialize the particle group
  call sll_s_new_particle_group_2d2v(particle_group, n_particles, &
       n_particles ,1.0_f64, 1.0_f64, 3)
  
  call particle_group%set_common_weight( 1.0_f64/real(n_particles,f64))
  do i_part = 1,n_particles
     xi(1:2) = x_vec(i_part,:)
     call particle_group%set_x(i_part, xi)
     wi(1:2) = w_vec(i_part,1:2)
     wi(3) = control_variate%update_df_weight(xi,  v_vec(i_part,:), 0.0_f64, wi(1),wi(2))
     call particle_group%set_weights(i_part, wi)
     xi(1:2) = v_vec(i_part,:)
     call particle_group%set_v(i_part, xi)
  end do

  do i_part = 1, n_particles
     wi = particle_group%get_weights(i_part)
     if (abs(wi(3)-w_vec(i_part,3))> 1D-15) then
        passed = .FALSE.
     end if
  end do

  
  if (passed .EQV. .TRUE.) then
     print*, 'PASSED'
  else
     print*, 'FAILED'
     stop
  end if

  deallocate(control_variate_parameter)

contains

!------------------------------------------------------------------------------!

  function control_variate_equi( this, xi, vi, time) result(sll_f_control_variate)
    class(sll_t_control_variate) :: this
    sll_real64, optional,  intent( in ) :: xi(:) !< particle position
    sll_real64, optional, intent( in ) :: vi(:) !< particle velocity
    sll_real64, optional, intent( in ) :: time  !< current time
    sll_real64               :: sll_f_control_variate


    sll_f_control_variate = exp(-0.5_f64*&
         ((vi(1)/this%control_variate_parameters(1))**2+&
         (vi(2)/this%control_variate_parameters(2))**2))/&
         (2.0_f64*sll_p_pi*product(this%control_variate_parameters))

  end function control_variate_equi


end program test_control_variate
