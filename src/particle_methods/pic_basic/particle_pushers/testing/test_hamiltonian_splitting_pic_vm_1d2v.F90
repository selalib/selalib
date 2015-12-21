! TODO: Use input from file to initialize and compare

program test_hamiltonian_splitting_pic_1d2v_vm
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_mpi, only: &
    MPI_SUM

  use sll_m_collective, only: &
    sll_s_boot_collective, &
    sll_o_collective_allreduce, &
    sll_s_halt_collective, &
    sll_v_world_collective

  use sll_m_constants, only: &
    sll_p_pi

  use sll_m_hamiltonian_splitting_pic_vm_1d2v, only: &
    sll_f_new_hamiltonian_splitting_pic_vm_1d2v, &
    sll_t_hamiltonian_splitting_pic_vm_1d2v

  use sll_m_kernel_smoother_base, only: &
    sll_p_galerkin, &
    sll_c_kernel_smoother

  use sll_m_kernel_smoother_spline_1d, only: &
    sll_t_kernel_smoother_spline_1d, &
    sll_f_new_smoother_spline_1d

  use sll_m_maxwell_1d_base, only: &
    sll_c_maxwell_1d_base

  use sll_m_maxwell_1d_fem, only: &
    sll_t_maxwell_1d_fem, &
    sll_f_new_maxwell_1d_fem

  use sll_m_particle_group_1d2v, only: &
    sll_f_new_particle_group_1d2v, &
    sll_t_particle_group_1d2v

  use sll_m_particle_group_base, only: &
    sll_c_particle_group_base

  use sll_m_particle_initializer, only: &
    sll_s_particle_initialize_sobol_landau_1d2v

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! Tolerance for comparison of real numbers: set it here!
  sll_real64, parameter :: EQV_TOL = 1.0e-14_f64

  ! Abstract particle group
  class(sll_c_particle_group_base), pointer :: particle_group
  class(sll_t_particle_group_1d2v), pointer :: pg

  ! Arrays for the fields
  sll_real64, pointer :: efield(:,:), efield_ref(:,:)
  sll_real64, pointer :: bfield(:), bfield_ref(:)
  sll_real64, pointer :: rho(:), rho_local(:)

  ! Abstract kernel smoothers
  class(sll_c_kernel_smoother), pointer :: kernel_smoother_0     
  class(sll_c_kernel_smoother), pointer :: kernel_smoother_1
  
  ! Maxwell solver 
  ! Abstract 
  class(sll_c_maxwell_1d_base), pointer :: maxwell_solver
  
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
  sll_real64 :: xi(3), wi(1)
  logical    :: passed
  sll_real64 :: error
  sll_int32  :: ierr   ! error code for SLL_ALLOCATE

  ! Reference
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

  domain = [eta_min, eta_max, eta_max - eta_min]

  ! Initialize
  pg => sll_f_new_particle_group_1d2v(n_particles, &
       n_particles ,1.0_f64, 1.0_f64, 1)
  particle_group => pg

  call particle_group%set_common_weight (1.0_f64)


  ! Initial particle information   
  ! Data produce with following call
  !call sll_s_particle_initialize_sobol_landau_2d2v(&
  !     particle_group, &
  !     [0.1_f64, 0.5_f64], &
  !     eta_min, &
  !     eta_max-eta_min, &
  !     [1.0_f64, 1.0_f64], &
  !     rnd_seed)
  !call sll_s_particle_initialize_sobol_landau_1d2v(particle_group, &
  !     [0.1_f64, 0.5_f64], domain(1),domain(3), &
  !     [1.0_f64, 1.0_f64] , rnd_seed)
  

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

  call particle_group%set_common_weight (1.0_f64)

  ! Initialize kernel smoother    
  kernel_smoother_1 => sll_f_new_smoother_spline_1d(&
       domain(1:2), [num_cells], &
       n_particles, degree_smoother-1, sll_p_galerkin) 
  kernel_smoother_0 => &
       sll_f_new_smoother_spline_1d(domain(1:2), [num_cells], &
       n_particles, degree_smoother, sll_p_galerkin) 
  
  ! Initialize Maxwell solver
  maxwell_solver => sll_f_new_maxwell_1d_fem([eta_min, eta_max], num_cells, &
       degree_smoother)
  
  SLL_ALLOCATE(efield(kernel_smoother_0%n_dofs,2),ierr)
  SLL_ALLOCATE(bfield(kernel_smoother_0%n_dofs),ierr)
  SLL_ALLOCATE(efield_ref(kernel_smoother_0%n_dofs,2),ierr)
  SLL_ALLOCATE(bfield_ref(kernel_smoother_0%n_dofs),ierr)
  SLL_ALLOCATE(rho(kernel_smoother_0%n_dofs),ierr)
  SLL_ALLOCATE(rho_local(kernel_smoother_0%n_dofs),ierr)

  efield = 1.0_f64

  rho_local = 0.0_f64
  do i_part = 1,n_particles
     xi = particle_group%get_x(i_part)
     wi(1) = particle_group%get_charge( i_part)
     call kernel_smoother_0%add_charge(xi(1), wi(1), rho_local)
  end do
  ! MPI to sum up contributions from each processor
  rho = 0.0_f64
  call sll_o_collective_allreduce( sll_v_world_collective, &
       rho_local, &
       kernel_smoother_0%n_dofs, MPI_SUM, rho)
  ! Solve Poisson problem
  call maxwell_solver%compute_E_from_rho(efield(:,1),&
       rho)
  bfield = 1.0_f64

  propagator => sll_f_new_hamiltonian_splitting_pic_vm_1d2v(maxwell_solver, &
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
  
  if (passed .EQV. .FALSE.) then
     print*, 'Error in operatorHp1.'
  end if
  
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

  if (passed .EQV. .FALSE.) then
     print*, 'Error in operatorHp2.'
  end if

  !print*, 'v(1)', size(pg%particle_array,1), size(pg%particle_array,2),pg%particle_array(2,:)

  call propagator%operatorHE(delta_t)

  ! Compare to reference
  ! Particle information after operatorV application 
  particle_info_ref = reshape([    11.627560396526469_f64,        5.5135182122431550_f64,       -1.4311461960364338_f64,       2.7611337138493959E-002_f64,  0.39511488788223620_f64,       -1.3903131848268702_f64,        6.8636759376074723_f64,        5.7026946767517002_f64    ], [n_particles, 4])
  ! Compare computed values to reference values
    do i_part=1,n_particles
     xi = particle_group%get_x(i_part)
     if (abs(xi(1)-particle_info_ref(i_part,1))> EQV_TOL) then
        passed = .FALSE.
        print*, 'x(1) of particle ', i_part,' wrong'
     end if
     xi = particle_group%get_v(i_part)
     if (abs(xi(1)-particle_info_ref(i_part,2))> EQV_TOL) then
        passed = .FALSE.
        print*, 'v(1) of particle ', i_part,' wrong'
        print*, i_part, xi(1), particle_info_ref(i_part,2)
     elseif (abs(xi(2)-particle_info_ref(i_part,3))> EQV_TOL) then
        passed = .FALSE.
        print*, 'v(2) of particle ', i_part,' wrong'
        print*,   i_part, xi(2), particle_info_ref(i_part,3)
     end if
     xi(1:1) = particle_group%get_charge(i_part)
     if (abs(xi(1)-particle_info_ref(i_part,4))> EQV_TOL) then
        passed = .FALSE.
        print*, 'weight of particle ', i_part,' wrong'
     end if
  end do

  if (passed .EQV. .FALSE.) then
     print*, 'Error in operatorE'
  end if

  call propagator%operatorHB(delta_t)

  ! Compare to reference
  ! Particle information after operatorV application 
  particle_info_ref = reshape([   11.627560396526469_f64,        5.5135182122431550_f64,       -1.4311461960364338_f64,       2.7611337138493959E-002_f64,  0.39511488788223620_f64,       -1.3903131848268702_f64,       6.8636759376074723_f64,        5.7026946767517002_f64      ], [n_particles, 4])
  ! Compare computed values to reference values
    do i_part=1,n_particles
     xi = particle_group%get_x(i_part)
     if (abs(xi(1)-particle_info_ref(i_part,1))> EQV_TOL) then
        passed = .FALSE.
        print*, 'x(1) of particle ', i_part,' wrong'
     end if
     xi = particle_group%get_v(i_part)
     if (abs(xi(1)-particle_info_ref(i_part,2))> EQV_TOL) then
        passed = .FALSE.
        print*, 'v(1) of particle ', i_part,' wrong'
        print*, i_part, xi(1), particle_info_ref(i_part,2)
     elseif (abs(xi(2)-particle_info_ref(i_part,3))> EQV_TOL) then
        passed = .FALSE.
        print*, 'v(2) of particle ', i_part,' wrong'
        print*,   i_part, xi(2), particle_info_ref(i_part,3)
     end if
     xi(1:1) = particle_group%get_charge(i_part)
     if (abs(xi(1)-particle_info_ref(i_part,4))> EQV_TOL) then
        passed = .FALSE.
        print*, 'weight of particle ', i_part,' wrong'
     end if
  end do

  if (passed .EQV. .FALSE.) then
     print*, 'Error in operatorHB.'
  end if


  bfield_ref = [ 0.96241425319343565_f64,        1.0704842857730226_f64, &
       0.87749470807080665_f64,        1.0585080035598995_f64, &
       1.0312067899733146_f64,       0.98524486372490983_f64, &
       0.99163866889590213_f64,        1.0359265946162064_f64, &
       0.96576229888372844_f64,        1.0213195333087741_f64 ];

  efield_ref = reshape([ 0.32668967300827889_f64,       0.21111816665325256 _f64, &
       -3.2797353014133206_f64,       0.93182059869704992_f64, &
       2.4099520529926513_f64,      -0.54074935301675420_f64, &
       -6.6508143734355582E-002_f64,  -4.4526750013638647_f64, &
       2.7314012140039408_f64,        2.4952249586987651_f64, &
       1.2896538770270727_f64,        0.45428428470806997_f64, &
       1.9372031866314203_f64,        1.2380946765497305 _f64, &
       0.83732795515227731_f64,        1.0244885213955153_f64, &
       1.1187987978241094_f64,       0.68789307296286761_f64, &
       1.0938576576751671_f64,       0.85201507883794003_f64 ], [num_cells,2]);


  error = maxval(bfield-bfield_ref)

  if (error> EQV_TOL) then
     passed = .FALSE.
     print*, 'b field wrong.'
  end if

  error = maxval(efield-efield_ref)

  if (error> EQV_TOL) then
     passed = .FALSE.
     print*, 'e field wrong.'
  end if



  if (passed .EQV. .TRUE.) then
     print*, 'PASSED'
  else
     print*, 'FAILED'
     stop
  end if
  
  call sll_s_halt_collective()

end program test_hamiltonian_splitting_pic_1d2v_vm
