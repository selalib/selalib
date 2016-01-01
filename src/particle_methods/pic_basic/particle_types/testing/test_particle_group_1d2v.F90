program test_particle_group_1d2v
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_particle_group_1d2v, only: &
    sll_f_new_particle_group_1d2v, &
    sll_t_particle_group_1d2v

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  
  class(sll_t_particle_group_1d2v), pointer :: particle_group

  sll_int32  :: n_particles
  sll_int32  :: n_total_particles
  sll_real64 :: charge
  sll_real64 :: mass
  sll_int32  :: n_weights

  sll_int32  :: i_part 
  sll_real64 :: x(3)
  logical    :: fail


  n_particles = 10
  n_total_particles = n_particles
  charge = -1.0_f64
  mass = 1.0_f64
  n_weights = 1

  particle_group =>  sll_f_new_particle_group_1d2v(n_particles, n_total_particles, charge, mass, n_weights)

  do i_part = 1, n_particles
     call particle_group%set_x(i_part, [real(i_part,f64), 0.0_f64, 0.0_f64])
     call particle_group%set_v(i_part, [real(i_part,f64)**2, 0.0_f64, 0.0_f64])
     call particle_group%set_weights(i_part, [real(i_part/n_particles, f64)])
  end do

  fail = .FALSE.
  i_part = 4
  x =  particle_group%get_x(i_part)
  if ( abs(x(1)- real(i_part,f64))> 1E-15) then
     fail = .TRUE.
  end if
  x = particle_group%get_v(i_part)
  if ( abs(x(1)- real(i_part,f64)**2)> 1E-15) then
     fail = .TRUE.
  end if
  x(1:1) = particle_group%get_charge(i_part)
  if ( abs(x(1)- charge*real(i_part/n_particles, f64))> 1E-15) then
     fail = .TRUE.
  end if
  x(1:1) = particle_group%get_mass(i_part)
  if ( abs(x(1)- mass*real(i_part/n_particles, f64))> 1E-15) then
     fail = .TRUE.
  end if

  if (fail .EQV. .FALSE.) then
     print*, 'PASSED'
  else
     print*, 'FAILED'
     stop
  end if

end program test_particle_group_1d2v
