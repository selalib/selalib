! Unit test for Boris-Yee scheme
! author: Katharina Kormann, IPP

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

  use sll_m_hamiltonian_splitting_pic_vm_1d2v_boris, only: &
    sll_t_hamiltonian_splitting_pic_vm_1d2v_boris

  use sll_m_particle_mesh_coupling_base, only: &
    sll_p_galerkin, &
    sll_c_particle_mesh_coupling

  use sll_m_particle_mesh_coupling_spline_1d, only: &
    sll_t_particle_mesh_coupling_spline_1d, &
    sll_s_new_particle_mesh_coupling_spline_1d_ptr

  use sll_m_maxwell_1d_base, only: &
    sll_c_maxwell_1d_base

  use sll_m_maxwell_1d_fem, only: &
    sll_t_maxwell_1d_fem

  use sll_m_particle_group_1d2v, only: &
    sll_t_particle_group_1d2v

  use sll_m_particle_group_base, only: &
    sll_c_particle_group_base

  use sll_m_particle_initializer, only: &
    sll_s_particle_initialize_sobol_landau_1d2v

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! Tolerance for comparison of real numbers: set it here!
  sll_real64, parameter :: eqv_tol = 1.0e-14_f64

  ! Abstract particle group
  class(sll_c_particle_group_base), pointer :: particle_group
  class(sll_t_particle_group_1d2v), pointer :: pg

  ! Arrays for the fields
  sll_real64, pointer :: efield(:,:), efield_ref(:,:)
  sll_real64, pointer :: bfield(:), bfield_ref(:)
  sll_real64, pointer :: rho(:), rho_local(:)

  ! Abstract kernel smoothers
  class(sll_c_particle_mesh_coupling), pointer :: kernel_smoother_0     
  class(sll_c_particle_mesh_coupling), pointer :: kernel_smoother_1
  
  ! Maxwell solver 
  ! Abstract 
  class(sll_c_maxwell_1d_base), pointer :: maxwell_solver
  
  ! Specific Hamiltonian splitting
  type(sll_t_hamiltonian_splitting_pic_vm_1d2v_boris) :: propagator
  
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
  num_cells = 16
  delta_t = 0.1_f64
  degree_smoother = 3
  passed = .TRUE.
  rnd_seed = 10

  domain = [eta_min, eta_max, eta_max - eta_min]

  ! Initialize
  allocate(pg)
  call pg%init(n_particles, &
       n_particles ,1.0_f64, 1.0_f64, 1)
  particle_group => pg

  call particle_group%set_common_weight (1.0_f64)


  SLL_ALLOCATE(particle_info_ref(4,n_particles), i_part)
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
  call sll_s_new_particle_mesh_coupling_spline_1d_ptr(kernel_smoother_1, &
       domain(1:2), [num_cells], &
       n_particles, degree_smoother-1, sll_p_galerkin) 
  call sll_s_new_particle_mesh_coupling_spline_1d_ptr(kernel_smoother_0, &
       domain(1:2), [num_cells], &
       n_particles, degree_smoother, sll_p_galerkin) 
  
  ! Initialize Maxwell solver
  allocate( sll_t_maxwell_1d_fem :: maxwell_solver )
  select type ( maxwell_solver )
  type is ( sll_t_maxwell_1d_fem )
     call maxwell_solver%init( [eta_min, eta_max], num_cells, &
          degree_smoother)
  end select
  
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
  
  call propagator%init( maxwell_solver, &
       kernel_smoother_0, kernel_smoother_1, particle_group, &
       efield, bfield, &
       eta_min, eta_max-eta_min)
  call propagator%staggering( 0.5_f64*delta_t )

  call propagator%strang_splitting( delta_t, 1 )
  
  ! Compare to reference
  ! Particle information after operatorV application 
  particle_info_ref = reshape([11.590250917552364_f64,       -1.5236851980054489_f64,       0.40632305631363136_f64,        6.8636759376074723_f64,        5.5028810915315152_f64,        1.1611806341228641E-002_f64,  -1.4090781780562809_f64,        5.7026946767517002_f64 ], [4,n_particles])
  ! Compare computed values to reference values
    do i_part=1,n_particles
       xi = particle_group%get_x(i_part)
     if (abs(xi(1)-particle_info_ref(1,i_part))> EQV_TOL) then
        passed = .FALSE.
        print*, 'x(1) of particle ', i_part,' wrong'
     end if
     xi = particle_group%get_v(i_part)
     if (abs(xi(1)-particle_info_ref(2,i_part))> EQV_TOL) then
        passed = .FALSE.
        print*, 'v(1) of particle ', i_part,' wrong'
        print*, i_part, xi(1), particle_info_ref(i_part,2)
     elseif (abs(xi(2)-particle_info_ref(3,i_part))> EQV_TOL) then
        passed = .FALSE.
        print*, 'v(2) of particle ', i_part,' wrong'
        print*,   i_part, xi(2), particle_info_ref(3,i_part)
     end if
     xi(1:1) = particle_group%get_charge(i_part)
     if (abs(xi(1)-particle_info_ref(4,i_part))> EQV_TOL) then
        passed = .FALSE.
        print*, 'weight of particle ', i_part,' wrong'
     end if
  end do


  bfield_ref = [1.0007152502482244_f64,       0.99445857477120148_f64,        1.0126078645236389_f64,       0.97492340832308344_f64,        1.0502493982082513_f64,       0.88560153553873056_f64,        1.1141170903091009_f64,       0.95020604146500376_f64,        1.0248316646998006_f64,       0.98745124326634226_f64,        1.0056520876592949_f64,       0.99896954592982978_f64,       0.99645963114519265_f64,        1.0122347396384064_f64,       0.98742382899460557_f64,        1.0040980952792931_f64   ]
  error = maxval(abs(bfield-bfield_ref))  
 if (error > eqv_tol ) then
     print*, 'bfield error too large.'
     passed = .false.
  end if
  
  efield_ref = reshape( [1.7801632633330595_f64,       0.51717913695931361_f64,        4.2412959490364575E-002_f64,  -1.1212332464258261_f64,       -1.1684637217416221_f64,       -3.7806966401897091_f64,        3.7363585338991756_f64,        1.1932927327940017_f64,        1.1033167281181586_f64,       -0.018180082909540971_f64, -0.56595990082616909_f64,       -1.6702327921823792_f64,       -1.8045561507212267_f64,       -4.3141455300350842_f64,        4.8236042130256775_f64,        1.5537560432220632_f64,       0.98879315272141322_f64,        1.0323154046944396_f64,       0.93329346828214133_f64,        1.1302445587553216_f64,       0.73558670810944515_f64,        1.6340701469431378_f64,       0.73779661553285170_f64,        1.1288774513495987_f64,       0.93385001285635261_f64,       1.0324077177717137_f64,       0.98801632510199811_f64,       0.99610949244376890 _f64,       1.0239154844069209_f64,       0.92782406399041817_f64,        1.0265970800925086_f64,       0.99441071503466349_f64   ], [num_cells, 2])
  error = maxval(abs(efield-efield_ref))  
  if (error > eqv_tol ) then
     print*, 'efield error too large.'
     passed = .false.
  end if
  
  if (passed .eqv. .true.) then
     print*, 'PASSED'
  else
     print*, 'FAILED'
     stop
  end if

  particle_group => null()
  call propagator%free()
  call pg%free()
  deallocate(pg)
  deallocate(efield)
  deallocate(efield_ref)
  deallocate(bfield)
  deallocate(bfield_ref)
  deallocate(rho)
  deallocate(rho_local)
  call kernel_smoother_0%free()
  deallocate(kernel_smoother_0)
  call kernel_smoother_1%free()
  deallocate(kernel_smoother_1)
  call maxwell_solver%free()
  deallocate(maxwell_solver)
  
  call sll_s_halt_collective()

end program test_hamiltonian_splitting_pic_1d2v_vm
