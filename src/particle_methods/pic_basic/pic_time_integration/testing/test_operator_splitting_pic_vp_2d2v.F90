! TODO: Use input from file to sll_o_initialize and compare

program test_operator_splitting_pic_vp_2d2v
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_constants, only: &
    sll_p_pi

  use sll_m_particle_mesh_coupling_base_1d, only: &
    sll_p_collocation, &
    sll_c_particle_mesh_coupling_1d

  use sll_m_particle_mesh_coupling_spline_2d, only: &
    sll_t_particle_mesh_coupling_spline_2d, &
    sll_s_new_particle_mesh_coupling_spline_2d_ptr

  use sll_m_operator_splitting_pic_vp_2d2v, only: &
    sll_t_operator_splitting_pic_vp_2d2v

  use sll_m_particle_group_2d2v, only: &
    sll_s_new_particle_group_2d2v_ptr, &
    sll_t_particle_group_2d2v

  use sll_m_particle_group_base, only: &
    sll_c_particle_group_base

  use sll_m_poisson_2d_periodic, only: &
    sll_f_new_poisson_2d_periodic

  use sll_m_poisson_2d_base, only: &
    sll_c_poisson_2d_base

 use sll_m_pic_poisson_base, only : &
       sll_c_pic_poisson

  use sll_m_pic_poisson_2d, only : &
       sll_s_new_pic_poisson_2d

  use sll_m_collective, only : &
       sll_s_boot_collective, sll_s_halt_collective

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! Tolerance for comparison of real numbers: set it here!
  sll_real64, parameter :: EQV_TOL = 1.0e-14_f64

  ! Abstract particle group
  class(sll_c_particle_group_base), pointer :: particle_group

  ! Array for efield
  !sll_real64, pointer :: efield(:,:)

  ! Abstract kernel smoother
  class(sll_c_particle_mesh_coupling_1d), pointer :: kernel_smoother
  
  ! Poisson solver
  class(sll_c_poisson_2d_base), pointer :: poisson_solver 

  ! PIC Poisson solver
  class(sll_c_pic_poisson), pointer  :: solver
  
  ! Specific operator splitting
  type(sll_t_operator_splitting_pic_vp_2d2v) :: propagator
  
  ! Parameters
  sll_int32  :: n_particles
  sll_real64 :: eta_min(2), eta_max(2), domain(2,2)
  sll_int32  :: num_cells(2)
  sll_real64 :: delta_t
  sll_int32  :: degree_smoother
  sll_int64  :: rnd_seed

  ! Helper 
  sll_int32  :: i_part
  sll_real64 :: xi(3)
  logical    :: passed

  ! References
  sll_real64, allocatable :: particle_info_ref(:,:)
   
  call sll_s_boot_collective()

  ! Set parameters
  n_particles = 2!10
  eta_min = 0.0_f64
  eta_max = 4.0_f64*sll_p_pi
  num_cells = 10
  delta_t = 0.1_f64
  degree_smoother = 3
  passed = .TRUE.
  rnd_seed = 10
  
  ! Initialize
  call sll_s_new_particle_group_2d2v_ptr(particle_group, n_particles, &
       n_particles ,1.0_f64, 1.0_f64, 1)

  ! Initial particle information   
  ! Data produce with following call
  !call sll_s_particle_initialize_sobol_landau_2d2v(&
  !     particle_group, &
  !     [0.1_f64, 0.5_f64], &
  !     eta_min, &
  !     eta_max-eta_min, &
  !     [1.0_f64, 1.0_f64], &
  !     rnd_seed)
  SLL_ALLOCATE(particle_info_ref(n_particles,5), i_part)
  particle_info_ref = 0.0_f64
  particle_info_ref = reshape([   11.780972450961723_f64 ,       &
5.4977871437821380_f64,       0.78539816339744828_f64,        &
7.0685834705770345_f64,       0.15731068461017067_f64,       &
-1.5341205443525459_f64,        1.5341205443525459_f64,      &
-0.15731068461017067_f64,        86.251495608834688_f64,     &
   71.662174808595040_f64], [n_particles, 5])

  ! sll_o_initialize particles from particle_info_ref
  xi = 0.0_f64
  do i_part =1, n_particles
     xi(1:2) = particle_info_ref(i_part, 1:2) 
     call particle_group%set_x(i_part, xi)
     xi(1:2) = particle_info_ref(i_part, 3:4)
     call particle_group%set_v(i_part, xi)
     xi(1) = particle_info_ref(i_part, 5)
     call particle_group%set_weights(i_part, xi(1))
  end do

  domain(:,1) = eta_min
  domain(:,2) = eta_max
  call sll_s_new_particle_mesh_coupling_spline_2d_ptr(kernel_smoother, &
         domain, num_cells, n_particles, &
         degree_smoother, sll_p_collocation)
  
  poisson_solver => sll_f_new_poisson_2d_periodic( &
       eta_min(1), eta_max(1), num_cells(1), &
       eta_min(2), eta_max(2), num_cells(2))

  call sll_s_new_pic_poisson_2d(solver, num_cells, poisson_solver, kernel_smoother)

  call propagator%init(solver, particle_group)

  call propagator%operatorT(delta_t)

  ! Compare to reference
  ! Particle information after operatorV application 
  particle_info_ref = reshape([   11.796703519422740_f64,      &
  5.3443750893468831_f64,       0.93881021783270291_f64,       &
 7.0528524021160175_f64,       0.15731068461017067_f64,       &
-1.5341205443525459_f64,        1.5341205443525459_f64,      &
-0.15731068461017067_f64,        86.251495608834688_f64,      &
  71.662174808595040_f64], [n_particles, 5])
  ! Compare computed values to reference values
  do i_part=1,n_particles
     xi = particle_group%get_x(i_part)
     if (abs(xi(1)-particle_info_ref(i_part,1))> EQV_TOL) then
        passed = .FALSE.
     elseif (abs(xi(2)-particle_info_ref(i_part,2))> EQV_TOL) then
        passed = .FALSE.
     end if
     xi = particle_group%get_v(i_part)
     if (abs(xi(1)-particle_info_ref(i_part,3))> EQV_TOL) then
        passed = .FALSE.
     elseif (abs(xi(2)-particle_info_ref(i_part,4))> EQV_TOL) then
        passed = .FALSE.
     end if
     xi(1:1) = particle_group%get_charge(i_part)
     if (abs(xi(1)-particle_info_ref(i_part,5))> EQV_TOL) then
        passed = .FALSE.
     end if
  end do
  if (passed .EQV. .FALSE.) then
     print*, 'Error in operatorT.'
  end if

  call propagator%charge_deposition()
  call propagator%field_solver()
  call propagator%operatorV(delta_t)
  ! Compare to reference
  ! Particle information after operatorV application 
  particle_info_ref = reshape([   11.796703519422740_f64,    &
    5.3443750893468831_f64,       0.93881021783270291_f64,   &
     7.0528524021160175_f64,       0.15346588344558551_f64,  &
     -1.5294930005799796_f64,        1.5302579166423358_f64, &
     -0.15266168508042444_f64,        86.251495608834688_f64,&
        71.662174808595040_f64  ], [n_particles, 5])
  ! Compare computed values to reference values
  do i_part=1,n_particles
     xi = particle_group%get_x(i_part)
     if (abs(xi(1)-particle_info_ref(i_part,1))> EQV_TOL) then
        passed = .FALSE.
     elseif (abs(xi(2)-particle_info_ref(i_part,2))> EQV_TOL) then
        passed = .FALSE.
     end if
     xi = particle_group%get_v(i_part)
     if (abs(xi(1)-particle_info_ref(i_part,3))> EQV_TOL) then
        passed = .FALSE.
     elseif (abs(xi(2)-particle_info_ref(i_part,4))> EQV_TOL) then
        passed = .FALSE.
     end if
     xi(1:1) = particle_group%get_charge(i_part)
     if (abs(xi(1)-particle_info_ref(i_part,5))> EQV_TOL) then
        passed = .FALSE.
     end if
  end do

  if (passed .EQV. .FALSE.) then
     print*, 'Error in operatorV.'
  end if

  if (passed .EQV. .TRUE.) then
     print*, 'PASSED'
  else
     print*, 'FAILED'
     stop
  end if

  call particle_group%free()
  deallocate(particle_group)
  call kernel_smoother%free()
  deallocate(kernel_smoother)
  call poisson_solver%free()
  deallocate(poisson_solver)
  call solver%free()
  deallocate(solver)
  
  call sll_s_halt_collective()

end program test_operator_splitting_pic_vp_2d2v
