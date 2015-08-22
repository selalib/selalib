! This unit test does not test the whole pusher but just the accumulation and integration of j.

program unit_test_operator_splitting_pic_1d2v_vm

#include "sll_working_precision.h"
  
  use sll_collective

  use sll_m_kernel_smoother_base
  use sll_m_kernel_smoother_spline_1d  
  use sll_module_pic_base  
  use sll_m_maxwell_1d_base
  use sll_m_maxwell_1d_fem
  use sll_m_particle_group_1d2v

  use sll_m_operator_splitting_pic_1d2v_vm

  ! Split operator
  class(sll_operator_splitting_pic_1d2v_vm), pointer :: propagator

  ! Abstract kernel smoother
  class(sll_kernel_smoother_base), pointer :: kernel_0
  class(sll_kernel_smoother_base), pointer :: kernel_1
  ! Specific kernel smoother
  class(sll_kernel_smoother_spline_1d),pointer :: specific_kernel_0
  class(sll_kernel_smoother_spline_1d),pointer :: specific_kernel_1
  ! Abstract particle group
  class(sll_particle_group_base), pointer :: particle_group
  ! Specific particle group
  class(sll_particle_group_1d2v), pointer :: specific_particle_group 
  ! Abstract maxwell solver
  class(sll_maxwell_1d_base), pointer :: maxwell_solver
  ! Specific Maxwell solver
  class(sll_maxwell_1d_fem), pointer :: specific_maxwell_solver
  ! Parameters for the test
  sll_int32 :: n_cells
  sll_int32 :: n_particles
  sll_int32 :: spline_degree
  sll_real64 :: domain(2)
  sll_real64 :: x_vec(4)
  sll_real64 :: v_vec(4,2)
  sll_real64 :: delta_t

  ! helper variables
  sll_int32  :: i_part 
  sll_real64 :: xi(3)
  logical :: passed
  sll_real64 :: error

  ! J dofs
  sll_real64 :: j_dofs_ref(10,2)


  ! Note: Test is only serial.
  call sll_boot_collective()
  
  ! 
  passed = .TRUE.

  ! This tests the kernel smoother for a fixed particle and grid and spline degree 3.
  ! Test parameters
  n_cells = 10; ! Number of cells
  n_particles = 4 ! Number of particles
  spline_degree = 3 ! Spline degree
  domain = [0.0_f64, 2.0_f64] ! x_min, x_max
  x_vec = [0.1_f64, 0.65_f64, 0.7_f64, 1.5_f64] ! Particle positions
  v_vec(:,1) = [1.5_f64, 0.0_f64, 0.0_f64, 0.0_f64]
  v_vec(:,2) = v_vec(:,1)
  !v_vec(:,2) = [0.0_f64, 0.0_f64, 0.1_f64, 0.0_f64]

  ! We need to initialize the particle group
  specific_particle_group => sll_new_particle_group_1d2v(n_particles, &
       n_particles ,1.0_f64, 1.0_f64)
  
  do i_part = 1,n_particles
     xi(1) = x_vec(i_part)
     call specific_particle_group%set_x(i_part, xi)
     call specific_particle_group%set_weight(i_part, 1/real(n_particles,f64))
     xi(1:2) = v_vec(i_part,:)
     call specific_particle_group%set_v(i_part, xi)
  end do
  
  particle_group => specific_particle_group

  ! Initialize the field solver 
  specific_maxwell_solver => sll_new_maxwell_1d_fem(domain, n_cells, spline_degree)
  maxwell_solver => specific_maxwell_solver

  ! Initialize the kernel
  specific_kernel_0 => sll_new_smoother_spline_1d(&
       domain, [n_cells], n_particles, spline_degree-1)
  kernel_0 => specific_kernel_0
  specific_kernel_1 => sll_new_smoother_spline_1d(&
       domain, [n_cells], n_particles, spline_degree)
  kernel_1 => specific_kernel_1

  ! Initialize propagator
  propagator => sll_new_splitting_pic_1d2v_vm(&
       maxwell_solver, kernel_0, kernel_1,&
       particle_group, domain(1), domain(2)-domain(1))

  ! Tests over one interval
  delta_t = 0.03_f64
  call propagator%operatorHf(delta_t)
  j_dofs_ref(:,1) = [ 1.6004882812500010E-002_f64,   0.0000000000000000_f64, &
       0.0000000000000000_f64,        0.0000000000000000_f64, &
       0.0000000000000000_f64,        0.0000000000000000_f64, &
       0.0000000000000000_f64,        0.0000000000000000_f64, &
       6.5126953124999970E-003_f64,   6.1857421875000040E-002_f64];
  j_dofs_ref(:,2) =   [3.3403381347656275E-003_f64,   0.0000000000000000_f64, &
       0.0000000000000000_f64,        0.0000000000000000_f64, &
       0.0000000000000000_f64,        0.0000000000000000_f64, &
       0.0000000000000000_f64,        8.8720092773437315E-004_f64, &
       3.4261248779296871E-002_f64,   4.5886212158203141E-002_f64];

  error = maxval(abs(j_dofs_ref-propagator%j_dofs))
  if (error > 1.e-14) then
     passed = .FALSE.
  end if


  ! Reset x for second test
  do i_part = 1,n_particles
     xi(1) = x_vec(i_part)
     call specific_particle_group%set_x(i_part, xi)
  end do

  ! Test over three intervals
  delta_t = 0.25_f64
  call propagator%operatorHf(delta_t)
  j_dofs_ref(:,1) = [0.35192871093750000_f64,       0.15258789062499989_f64, &
       3.2958984374999883E-003_f64,   0.0000000000000000_f64, &
       0.0000000000000000_f64,        0.0000000000000000_f64, &
       0.0000000000000000_f64,        0.0000000000000000_f64, &
       7.8124999999999931E-003_f64,  0.18750000000000000_f64];
  j_dofs_ref(:,2) =   [0.27460861206054676_f64, &
       5.4615020751953063E-002_f64,   3.0899047851562354E-004_f64,   &
       0.0000000000000000_f64,         0.0000000000000000_f64, &
       0.0000000000000000_f64,        0.0000000000000000_f64, &
       9.7656249999999913E-004_f64,   7.5195312499999986E-002_f64, &
       0.29742050170898438_f64 ];

  error = maxval(abs(j_dofs_ref-propagator%j_dofs))
  if (error > 1.e-14) then
     passed = .FALSE.
  end if

  if (passed .EQV. .TRUE.) then
     print*, 'PASSED'
  else
     print*, 'FAILED'
     stop
  end if

  call sll_halt_collective()

end program unit_test_operator_splitting_pic_1d2v_vm
