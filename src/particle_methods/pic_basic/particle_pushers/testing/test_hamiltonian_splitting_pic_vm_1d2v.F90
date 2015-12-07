! TODO: Use input from file to initialize and compare

program test_hamiltonian_splitting_pic_1d2v_vm
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_constants, only: &
    sll_pi

  use sll_m_hamiltonian_splitting_pic_vm_1d2v, only: &
    sll_new_hamiltonian_splitting_pic_vm_1d2v, &
    sll_t_hamiltonian_splitting_pic_vm_1d2v

  use sll_m_kernel_smoother_base, only: &
    sll_galerkin, &
    sll_kernel_smoother_base

  use sll_m_kernel_smoother_spline_1d, only: &
    sll_kernel_smoother_spline_1d, &
    sll_new_smoother_spline_1d

  use sll_m_maxwell_1d_base, only: &
    sll_maxwell_1d_base

  use sll_m_maxwell_1d_fem, only: &
    sll_maxwell_1d_fem, &
    sll_new_maxwell_1d_fem

  use sll_m_particle_group_1d2v, only: &
    sll_new_particle_group_1d2v, &
    sll_particle_group_1d2v

  use sll_m_particle_group_base, only: &
    sll_particle_group_base

  use sll_m_particle_initializer, only: &
    sll_particle_initialize_sobol_landau_1d2v

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! Tolerance for comparison of real numbers: set it here!
  sll_real64, parameter :: EQV_TOL = 1.0e-14_f64

  ! Abstract particle group
  class(sll_particle_group_base), pointer :: particle_group
  ! Specific particle group
  class(sll_particle_group_1d2v), pointer :: specific_particle_group 

  ! Arrays for the fields
  sll_real64, pointer :: efield(:,:), efield_ref(:,:)
  sll_real64, pointer :: bfield(:), bfield_ref(:)

  ! Abstract kernel smoothers
  class(sll_kernel_smoother_base), pointer :: kernel_smoother_0     
  class(sll_kernel_smoother_base), pointer :: kernel_smoother_1
  ! Specific kernel smoother
  class(sll_kernel_smoother_spline_1d), pointer :: specific_kernel_smoother_0
  class(sll_kernel_smoother_spline_1d), pointer :: specific_kernel_smoother_1
  
  ! Maxwell solver 
  ! Abstract 
  class(sll_maxwell_1d_base), pointer :: maxwell_solver
  ! Specific
  class(sll_maxwell_1d_fem), pointer :: specific_maxwell_solver
  
  ! Specific Hamiltonian splitting
  class(sll_t_hamiltonian_splitting_pic_vm_1d2v), pointer :: propagator
  
  ! Parameters
  sll_int32  :: n_particles
  sll_real64 :: eta_min, eta_max, domain(3)
  sll_int32  :: num_cells
  sll_real64 :: delta_t
  sll_int32  :: degree_smoother
  sll_int64  :: rnd_seed

  ! Helper 
  sll_int32  :: i_part
  sll_real64 :: xi(3)
  logical    :: passed
  sll_real64 :: error
  sll_int32  :: ierr   ! error code for SLL_ALLOCATE

  ! Reference
  sll_real64, allocatable :: particle_info_ref(:,:)
   
!  call sll_boot_collective()

  ! Set parameters
  n_particles = 2!10
  eta_min = 0.0_f64
  eta_max = 4.0_f64*sll_pi
  num_cells = 10
  delta_t = 0.1_f64
  degree_smoother = 3
  passed = .TRUE.
  rnd_seed = 10

  domain = [eta_min, eta_max, eta_max - eta_min]
  
  ! Initialize
  specific_particle_group => sll_new_particle_group_1d2v(n_particles, &
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
  call sll_particle_initialize_sobol_landau_1d2v(particle_group, &
       [0.1_f64, 0.5_f64], domain(1),domain(3), &
       [1.0_f64, 1.0_f64] , rnd_seed)
  

  SLL_ALLOCATE(particle_info_ref(n_particles,4), i_part)
  particle_info_ref = 0.0_f64
  particle_info_ref = reshape([11.780972450961723_f64,        5.4977871437821380_f64,       -1.5341205443525459_f64,       0.15731068461017067_f64,       0.15731068461017067_f64,       -1.5341205443525459_f64,        6.8636759376074723_f64,        5.7026946767517002_f64   ], [n_particles, 4])

  ! Initialize particles from particle_info_ref
  xi = 0.0_f64
  do i_part =1, n_particles
     xi(1) = particle_info_ref(i_part, 1) 
     call particle_group%set_x(i_part, xi)
     xi(1:2) = particle_info_ref(i_part, 2:3)
     call particle_group%set_v(i_part, xi)
     xi(1) = particle_info_ref(i_part, 4)
     call particle_group%set_weights(i_part, xi(1))
  end do

  ! Initialize kernel smoother    
  specific_kernel_smoother_1 => sll_new_smoother_spline_1d(&
       domain(1:2), [num_cells], &
       n_particles, degree_smoother-1, SLL_GALERKIN) 
  kernel_smoother_1 => specific_kernel_smoother_1
  specific_kernel_smoother_0 => &
       sll_new_smoother_spline_1d(domain(1:2), [num_cells], &
       n_particles, degree_smoother, SLL_GALERKIN) 
  kernel_smoother_0 => specific_kernel_smoother_0
  
  ! Initialize Maxwell solver
  specific_maxwell_solver => sll_new_maxwell_1d_fem([eta_min, eta_max], num_cells, &
       degree_smoother)
  maxwell_solver => specific_maxwell_solver

  SLL_ALLOCATE(efield(kernel_smoother_0%n_dofs,2),ierr)
  SLL_ALLOCATE(bfield(kernel_smoother_0%n_dofs),ierr)
  SLL_ALLOCATE(efield_ref(kernel_smoother_0%n_dofs,2),ierr)
  SLL_ALLOCATE(bfield_ref(kernel_smoother_0%n_dofs),ierr)

  efield = 1.0_f64
  bfield = 1.0_f64

  propagator => sll_new_hamiltonian_splitting_pic_vm_1d2v(maxwell_solver, &
            kernel_smoother_0, kernel_smoother_1, particle_group, &
            efield, bfield, &
            eta_min, eta_max-eta_min)

  call propagator%operatorHp1(delta_t)

  ! Compare to reference
  ! Particle information after operatorV application 
  particle_info_ref = reshape([  11.627560396526469_f64,        5.5135182122431550_f64,       -1.5341205443525459_f64,       0.15731068461017067_f64,       0.31072273904542569_f64,       -1.5498516128135633_f64,        6.8636759376074723_f64,        5.7026946767517002_f64    ], [n_particles, 4])
  ! Compare computed values to reference values
  do i_part=1,n_particles
     xi = particle_group%get_x(i_part)
     if (abs(xi(1)-particle_info_ref(i_part,1))> EQV_TOL) then
        passed = .FALSE.
     end if
     xi = particle_group%get_v(i_part)
     if (abs(xi(1)-particle_info_ref(i_part,2))> EQV_TOL) then
        passed = .FALSE.
     elseif (abs(xi(2)-particle_info_ref(i_part,3))> EQV_TOL) then
        passed = .FALSE.
     end if
     xi(1:1) = particle_group%get_charge(i_part)
     if (abs(xi(1)-particle_info_ref(i_part,4))> EQV_TOL) then
        passed = .FALSE.
     end if
  end do

  
  call propagator%operatorHp2(delta_t)

  ! Compare to reference
  ! Particle information after operatorV application 
  particle_info_ref = reshape([ 11.627560396526469_f64,        5.5135182122431550_f64,       -1.5030482704480033_f64,        2.3255233288143329D-003,  0.31072273904542569_f64,       -1.5498516128135633_f64,        6.8636759376074723_f64,        5.7026946767517002_f64   ], [n_particles, 4])
  ! Compare computed values to reference values
    do i_part=1,n_particles
     xi = particle_group%get_x(i_part)
     if (abs(xi(1)-particle_info_ref(i_part,1))> EQV_TOL) then
        passed = .FALSE.
     end if
     xi = particle_group%get_v(i_part)
     if (abs(xi(1)-particle_info_ref(i_part,2))> EQV_TOL) then
        passed = .FALSE.
     elseif (abs(xi(2)-particle_info_ref(i_part,3))> EQV_TOL) then
        passed = .FALSE.
     end if
     xi(1:1) = particle_group%get_charge(i_part)
     if (abs(xi(1)-particle_info_ref(i_part,4))> EQV_TOL) then
        passed = .FALSE.
     end if
  end do

  call propagator%operatorHE(delta_t)

  ! Compare to reference
  ! Particle information after operatorV application 
  particle_info_ref = reshape([    11.627560396526469_f64,        5.5135182122431550_f64,       -1.3696420089578860_f64,       9.9430804248216542e-2_f64,  0.40491306194184096_f64,       -1.4259284606804052_f64,        6.8636759376074723_f64,        5.7026946767517002_f64    ], [n_particles, 4])
  ! Compare computed values to reference values
    do i_part=1,n_particles
     xi = particle_group%get_x(i_part)
     if (abs(xi(1)-particle_info_ref(i_part,1))> EQV_TOL) then
        passed = .FALSE.
     end if
     xi = particle_group%get_v(i_part)
     if (abs(xi(1)-particle_info_ref(i_part,2))> EQV_TOL) then
        passed = .FALSE.
     elseif (abs(xi(2)-particle_info_ref(i_part,3))> EQV_TOL) then
        passed = .FALSE.
     end if
     xi(1:1) = particle_group%get_charge(i_part)
     if (abs(xi(1)-particle_info_ref(i_part,4))> EQV_TOL) then
        passed = .FALSE.
     end if
  end do


  call propagator%operatorHB(delta_t)

  ! Compare to reference
  ! Particle information after operatorV application 
  particle_info_ref = reshape([   11.627560396526469_f64,        5.5135182122431550_f64,       -1.3696420089578860_f64,        9.9430804248216542D-002,  0.40491306194184096_f64,       -1.4259284606804052_f64,        6.8636759376074723_f64,        5.7026946767517002_f64      ], [n_particles, 4])
  ! Compare computed values to reference values
    do i_part=1,n_particles
     xi = particle_group%get_x(i_part)
     if (abs(xi(1)-particle_info_ref(i_part,1))> EQV_TOL) then
        passed = .FALSE.
     end if
     xi = particle_group%get_v(i_part)
     if (abs(xi(1)-particle_info_ref(i_part,2))> EQV_TOL) then
        passed = .FALSE.
     elseif (abs(xi(2)-particle_info_ref(i_part,3))> EQV_TOL) then
        passed = .FALSE.
     end if
     xi(1:1) = particle_group%get_charge(i_part)
     if (abs(xi(1)-particle_info_ref(i_part,4))> EQV_TOL) then
        passed = .FALSE.
     end if
  end do


  bfield_ref = [  0.99770541216286657_f64,       0.99234240623721359_f64,       0.98825065652414523_f64,        1.0022336143723096_f64,        1.0119601928906674_f64,        1.0060821139662497_f64,        1.0027406207578249_f64,        1.0027048521768871_f64,       0.99884366018700488_f64,       0.99713647072483091_f64]

  efield_ref = reshape([1.0139333322572233_f64,       0.99694980433046210_f64,       0.98105828960960473_f64,       0.96702117187072789_f64,       0.98564985298325469_f64,        1.0000855014667056_f64,        1.0477798381501124_f64,        1.2387410309599800_f64,        1.3810186780785210_f64,        1.1543013648190512_f64,        1.0145964909357603_f64,        1.1109252342944447_f64,        1.2566590132766358_f64,        1.2290902864314954_f64,        1.0803226658682838_f64,        1.0034722970882306_f64,       0.96871172275213535_f64,       0.93518373011620015_f64,       0.94946045559169057_f64,       0.98519521240929375_f64 ],[num_cells,2])

  error = maxval(bfield-bfield_ref)

  if (error> EQV_TOL) then
     passed = .FALSE.
  end if

  error = maxval(efield-efield_ref)

  if (error> EQV_TOL) then
     passed = .FALSE.
  end if



  if (passed .EQV. .TRUE.) then
     print*, 'PASSED'
  else
     print*, 'FAILED'
     stop
  end if
  
 ! call sll_halt_collective()
end program test_hamiltonian_splitting_pic_1d2v_vm
