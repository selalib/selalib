! TODO: Use input from file to sll_o_initialize and compare

program test_operator_splitting_pic_vp_2d2v
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_constants, only: &
    sll_p_pi

  use sll_m_kernel_smoother_base, only: &
    sll_p_collocation, &
    sll_c_kernel_smoother_base

  use sll_m_kernel_smoother_spline_2d, only: &
    sll_t_kernel_smoother_spline_2d, &
    sll_f_new_smoother_spline_2d

  use sll_m_operator_splitting_pic_vp_2d2v, only: &
    sll_f_new_hamiltonian_splitting_pic_vp_2d2v, &
    sll_t_operator_splitting_pic_vp_2d2v

  use sll_m_particle_group_2d2v, only: &
    sll_f_new_particle_group_2d2v, &
    sll_t_particle_group_2d2v

  use sll_m_particle_group_base, only: &
    sll_c_particle_group_base

  use sll_m_poisson_2d_periodic_fft, only: &
    sll_f_new_poisson_2d_periodic_fft, &
    sll_t_poisson_2d_periodic_fft

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! Tolerance for comparison of real numbers: set it here!
  sll_real64, parameter :: EQV_TOL = 1.0e-14_f64

  ! Abstract particle group
  class(sll_c_particle_group_base), pointer :: particle_group
  ! Specific particle group
  class(sll_t_particle_group_2d2v), pointer :: specific_particle_group 

  ! Array for efield
  sll_real64, pointer :: efield(:,:)

  ! Abstract kernel smoother
  class(sll_c_kernel_smoother_base), pointer :: kernel_smoother
  ! Specific kernel smoother
  class(sll_t_kernel_smoother_spline_2d), pointer :: specific_kernel_smoother
  
  ! Poisson solver
  class(sll_t_poisson_2d_periodic_fft), pointer :: poisson_solver 
  
  ! Specific operator splitting
  class(sll_t_operator_splitting_pic_vp_2d2v), pointer :: propagator
  
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
  sll_int32  :: ierr   ! error code for SLL_ALLOCATE

  ! Reference
  sll_real64, allocatable :: particle_info_ref(:,:)
   
  !call sll_s_boot_collective()

  ! Set parameters
  n_particles = 2!10
  eta_min = 0.0_f64
  eta_max = 4.0_f64*sll_p_pi
  num_cells = 10
  delta_t = 0.1_f64
  degree_smoother = 3
  passed = .TRUE.
  rnd_seed = 10
  
  ! sll_o_initialize
  specific_particle_group => sll_f_new_particle_group_2d2v(n_particles, &
       n_particles ,1.0_f64, 1.0_f64, 1)
  particle_group => specific_particle_group

  ! Initial particle information   
  ! Data produce with following call
  !call sll_particle_initialize_sobol_landau_2d2v(&
  !     particle_group, &
  !     [0.1_f64, 0.5_f64], &
  !     eta_min, &
  !     eta_max-eta_min, &
  !     [1.0_f64, 1.0_f64], &
  !     rnd_seed)
  SLL_ALLOCATE(particle_info_ref(n_particles,5), i_part)
  particle_info_ref = 0.0_f64
  particle_info_ref = reshape([   11.780972450961723_f64 ,       5.4977871437821380_f64,       0.78539816339744828_f64,        7.0685834705770345_f64,       0.15731068461017067_f64,       -1.5341205443525459_f64,        1.5341205443525459_f64,      -0.15731068461017067_f64,        86.251495608834688_f64,        71.662174808595040_f64], [n_particles, 5])

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
  specific_kernel_smoother => sll_f_new_smoother_spline_2d(&
         domain, num_cells, n_particles, &
         degree_smoother, sll_p_collocation)
  kernel_smoother => specific_kernel_smoother
  
  poisson_solver => sll_f_new_poisson_2d_periodic_fft( &
       eta_min(1), eta_max(1), num_cells(1), &
       eta_min(2), eta_max(2), num_cells(2))
  SLL_ALLOCATE(efield(kernel_smoother%n_dofs,2),ierr)
  propagator => sll_f_new_hamiltonian_splitting_pic_vp_2d2v(poisson_solver, kernel_smoother, &
       particle_group, efield)


  call propagator%operatorT(delta_t)

  ! Compare to reference
  ! Particle information after operatorV application 
  particle_info_ref = reshape([   11.796703519422740_f64,        5.3443750893468831_f64,       0.93881021783270291_f64,        7.0528524021160175_f64,       0.15731068461017067_f64,       -1.5341205443525459_f64,        1.5341205443525459_f64,      -0.15731068461017067_f64,        86.251495608834688_f64,        71.662174808595040_f64], [n_particles, 5])
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

  call propagator%operatorV(delta_t)
  ! Compare to reference
  ! Particle information after operatorV application 
  particle_info_ref = reshape([   11.796703519422740_f64,        5.3443750893468831_f64,       0.93881021783270291_f64,        7.0528524021160175_f64,       0.15346588344558551_f64,       -1.5294930005799796_f64,        1.5302579166423358_f64,      -0.15266168508042444_f64,        86.251495608834688_f64,        71.662174808595040_f64  ], [n_particles, 5])
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

  if (passed .EQV. .TRUE.) then
     print*, 'PASSED'
  else
     print*, 'FAILED'
     stop
  end if
  
  !call sll_s_halt_collective()
end program test_operator_splitting_pic_vp_2d2v
