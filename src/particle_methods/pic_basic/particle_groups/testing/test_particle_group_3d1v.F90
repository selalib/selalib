program test_particle_group_3d1v
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_particle_group_3d1v, only: &
       sll_t_particle_group_3d1v, &
       sll_s_new_particle_group_3d1v, &
       sll_s_new_particle_group_3d1v_ptr

  use sll_m_particle_group_base, only : &
       sll_c_particle_group_base

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  type(sll_t_particle_group_3d1v) :: particle_group
  class(sll_c_particle_group_base), pointer :: pgp
  
  sll_int32 :: n_particles
  sll_int32 :: n_total_particles
  sll_int32 :: n_weights
  sll_real64 :: charge
  sll_real64 :: mass
  
  sll_int32 :: i_part
  sll_real64 :: x(3)
  logical :: fail
  
  
  n_particles = 15
  n_total_particles = n_particles
  charge= -1.0_f64
  mass= 1.0_f64
  n_weights = 1
  fail=.false.
  
  call sll_s_new_particle_group_3d1v_ptr(pgp, n_particles, n_total_particles, charge, mass, n_weights)

  call particle_group%init(n_particles, n_total_particles, charge, mass, n_weights)

  do i_part = 1,n_particles
     call particle_group%set_x(i_part,[real(i_part, f64),real(i_part, f64)**0.5_f64,real(i_part, f64)**0.2_f64])
     call Particle_group%set_v(i_part,[real(i_part,f64)**2,0.0_f64,0.0_f64])
     call Particle_group%set_weights(i_part,[real(i_part/n_particles,f64)])
  end do

  i_part=4
  !do i_part = 1,n_particles,2
     
     x = particle_group%get_x(i_part)
     if(abs(x(1)-real(i_part, f64))>1d-15) fail=.true.
     if(abs(x(2)-real(i_part, f64)**0.5_f64)>1d-15) fail=.true.
     if(abs(x(3)-real(i_part, f64)**0.2_f64)>1d-15) fail=.true.
     
     x = particle_group%get_v(i_part)
     if(abs(x(1)-real(i_part,f64)**2)>1d-15) fail=.true.
     
     x(1) = particle_group%get_mass(i_part)
     if(abs(x(1)-mass*real(i_part/n_particles, f64))>1d-15) fail=.true.
     
     x(1) = particle_group%get_charge(i_part)
     if(abs(x(1)-charge*real(i_part/n_particles, f64))>1d-15) fail=.true.
     
  !end do

  call particle_group%free()

  if(fail .eqv. .false.) then
     write(*,*) 'PASSED'
  else  
     write(*,*)'FAILED'
     stop
  end if


end program  test_particle_group_3d1v
