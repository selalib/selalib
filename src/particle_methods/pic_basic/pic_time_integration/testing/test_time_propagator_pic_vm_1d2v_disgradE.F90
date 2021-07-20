! TODO: Use input from file to initialize and compare
! Unit test for symplectic splitting
! author: Katharina Kormann, IPP
program test_time_propagator_pic_1d2v_vm_disgradE
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_mpi, only: &
       MPI_SUM

  use sll_m_binomial_filter, only : &
       sll_t_binomial_filter

  use sll_m_collective, only: &
    sll_s_boot_collective, &
    sll_o_collective_allreduce, &
    sll_s_halt_collective, &
    sll_v_world_collective

  use sll_m_constants, only: &
    sll_p_pi

  use sll_m_time_propagator_pic_vm_1d2v_disgradE, only: &
    sll_t_time_propagator_pic_vm_1d2v_disgradE

  use sll_m_particle_mesh_coupling_base_1d, only: &
    sll_p_galerkin, &
    sll_c_particle_mesh_coupling_1d

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
    sll_t_particle_array

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! Tolerance for comparison of real numbers: set it here!
  sll_real64, parameter :: EQV_TOL = 1.5e-12_f64

  ! Abstract particle group
  class(sll_t_particle_array), pointer :: particle_group
  class(sll_t_particle_group_1d2v), pointer :: pg

  ! Arrays for the fields
  sll_real64, pointer :: phi(:), efield(:,:), efield_ref(:,:)
  sll_real64, pointer :: bfield(:), bfield_ref(:)
  sll_real64, pointer :: rho(:), rho_local(:)

  ! Abstract kernel smoothers
  class(sll_c_particle_mesh_coupling_1d), pointer :: kernel_smoother_0     
  class(sll_c_particle_mesh_coupling_1d), pointer :: kernel_smoother_1
  
  ! Maxwell solver 
  ! Abstract 
  class(sll_c_maxwell_1d_base), pointer :: maxwell_solver
  
  ! Specific Hamiltonian splitting
  type(sll_t_time_propagator_pic_vm_1d2v_disgradE) :: propagator

  ! Binomial filter object
  type( sll_t_binomial_filter), target :: filter
  
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
  sll_real64, allocatable :: particle_info_check(:,:)
   
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
  allocate(pg)
  call pg%init(n_particles, &
       n_particles ,1.0_f64, 1.0_f64, 1)
  allocate( particle_group)
  particle_group%n_species = 1
  allocate( particle_group%group(1), source=pg )

  call particle_group%group(1)%set_common_weight (1.0_f64)


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
  SLL_ALLOCATE(particle_info_check(n_particles,4), i_part)
  particle_info_ref = 0.0_f64
  particle_info_ref = reshape([11.7809724509617_f64,        5.49778714378213_f64,       -1.53412054435254_f64,       0.15731068461017_f64,       0.15731068461017_f64,       -1.53412054435254_f64,        6.86367593760747_f64,        5.7026946767517_f64   ], [n_particles, 4])

  ! Initialize particles from particle_info_ref
  xi = 0.0_f64
  do i_part =1, n_particles
     xi(1) = particle_info_ref(i_part, 1) 
     call particle_group%group(1)%set_x(i_part, xi)
     xi(1:2) = particle_info_ref(i_part, 2:3)
     call particle_group%group(1)%set_v(i_part, xi)
     xi(1) = particle_info_ref(i_part, 4)
     call particle_group%group(1)%set_weights(i_part, xi(1))
  end do

  call particle_group%group(1)%set_common_weight (1.0_f64)

  ! Initialize kernel smoother    
  call sll_s_new_particle_mesh_coupling_spline_1d_ptr(kernel_smoother_1, &
       domain(1:2), num_cells, &
       n_particles, degree_smoother-1, sll_p_galerkin) 
  call sll_s_new_particle_mesh_coupling_spline_1d_ptr(kernel_smoother_0, &
       domain(1:2), num_cells, &
       n_particles, degree_smoother, sll_p_galerkin) 
  
  ! Initialize Maxwell solver
  allocate( sll_t_maxwell_1d_fem :: maxwell_solver )
  select type ( maxwell_solver )
  type is ( sll_t_maxwell_1d_fem )
     call maxwell_solver%init( [eta_min, eta_max], num_cells, &
          degree_smoother, delta_t, solver_tolerance = 1d-14)
  end select
  SLL_ALLOCATE(phi(kernel_smoother_0%n_dofs),ierr)
  SLL_ALLOCATE(efield(kernel_smoother_0%n_dofs,2),ierr)
  SLL_ALLOCATE(bfield(kernel_smoother_0%n_dofs),ierr)
  SLL_ALLOCATE(efield_ref(kernel_smoother_0%n_dofs,2),ierr)
  SLL_ALLOCATE(bfield_ref(kernel_smoother_0%n_dofs),ierr)
  SLL_ALLOCATE(rho(kernel_smoother_0%n_dofs),ierr)
  SLL_ALLOCATE(rho_local(kernel_smoother_0%n_dofs),ierr)

  efield = 1.0_f64

  rho_local = 0.0_f64
  do i_part = 1,n_particles
     xi = particle_group%group(1)%get_x(i_part)
     wi(1) = particle_group%group(1)%get_charge( i_part)
     call kernel_smoother_0%add_charge(xi(1), wi(1), rho_local)
  end do
  ! MPI to sum up contributions from each processor
  rho = 0.0_f64
  call sll_o_collective_allreduce( sll_v_world_collective, &
       rho_local, &
       kernel_smoother_0%n_dofs, MPI_SUM, rho)
  ! Solve Poisson problem
  call maxwell_solver%compute_E_from_rho( rho, efield(:,1) )
  bfield = 1.0_f64

  call filter%init( 0, 0 )

  call propagator%init( maxwell_solver, &
       kernel_smoother_0, kernel_smoother_1, particle_group, &
       phi, efield, bfield, &
       eta_min, eta_max-eta_min, filter, solver_tolerance = 1d-14, adiabatic_electrons = .false. )
   
  call propagator%advect_x( delta_t )

  ! Compare to reference
  ! Particle information after advect_x application 
  particle_info_check(:,1) = [11.62756039652646_f64, 5.513518212243155_f64]
  particle_info_check(:,2) = [-1.53412054435254_f64, 0.157310684610170_f64]
  particle_info_check(:,3) = [ 0.15731068461017_f64,-1.534120544352545_f64]
  particle_info_check(:,4) = [ 6.86367593760747_f64, 5.702694676751700_f64]
  ! Compare computed values to reference values
  do i_part=1,n_particles
     xi = particle_group%group(1)%get_x(i_part)
     if (abs(xi(1)-particle_info_check(i_part,1))> EQV_TOL) then
        passed = .FALSE.
     end if
     xi = particle_group%group(1)%get_v(i_part)
     if (abs(xi(1)-particle_info_check(i_part,2))> EQV_TOL) then
        passed = .FALSE.
     elseif (abs(xi(2)-particle_info_check(i_part,3))> EQV_TOL) then
        passed = .FALSE.
     end if
     xi(1:1) = particle_group%group(1)%get_charge(i_part)
     if (abs(xi(1)-particle_info_check(i_part,4))> EQV_TOL) then
        passed = .FALSE.
     end if
  end do
  
  if (passed .EQV. .FALSE.) then
     print*, 'Error in advect_x.'
  end if

  ! Reset particle info
  ! Initialize particles from particle_info_ref
  xi = 0.0_f64
  do i_part =1, n_particles
     xi(1) = particle_info_ref(i_part, 1) 
     call particle_group%group(1)%set_x(i_part, xi)
     xi(1:2) = particle_info_ref(i_part, 2:3)
     call particle_group%group(1)%set_v(i_part, xi)
     xi(1) = particle_info_ref(i_part, 4)
     call particle_group%group(1)%set_weights(i_part, xi(1))
  end do

  
  call propagator%advect_vb( delta_t )

  ! Compare to reference
  ! Particle information after advect_vb application 
  particle_info_check(:,1) = [11.780972450961723_f64,    5.49778714378213_f64]
  particle_info_check(:,2) = [-1.510751468549690_f64,   3.368290939051138E-003_f64]
  particle_info_check(:,3) = [ 0.309681281920664_f64, -1.5421611947890606_f64]
  particle_info_check(:,4) = [ 6.863675937607472_f64,   5.702694676751700_f64]
  ! Compare computed values to reference values
  do i_part=1,n_particles
     xi = particle_group%group(1)%get_x(i_part)
     if (abs(xi(1)-particle_info_check(i_part,1))> EQV_TOL) then
        passed = .FALSE.
     end if
     xi = particle_group%group(1)%get_v(i_part)
     if (abs(xi(1)-particle_info_check(i_part,2))> EQV_TOL) then
        passed = .FALSE.
     elseif (abs(xi(2)-particle_info_check(i_part,3))> EQV_TOL) then
        passed = .FALSE.
     end if
     xi(1:1) = particle_group%group(1)%get_charge(i_part)
     if (abs(xi(1)-particle_info_check(i_part,4))> EQV_TOL) then
        passed = .FALSE.
     end if
  end do
  
  if (passed .EQV. .FALSE.) then
     print*, 'Error in advect_vb.'
  end if

  

  ! Reset particle info
  ! Initialize particles from particle_info_ref
  xi = 0.0_f64
  do i_part =1, n_particles
     xi(1) = particle_info_ref(i_part, 1) 
     call particle_group%group(1)%set_x(i_part, xi)
     xi(1:2) = particle_info_ref(i_part, 2:3)
     call particle_group%group(1)%set_v(i_part, xi)
     xi(1) = particle_info_ref(i_part, 4)
     call particle_group%group(1)%set_weights(i_part, xi(1))
  end do  

  
  call propagator%advect_e( delta_t )

  efield_ref = reshape([ 0.31627789686666_f64,       0.21472175951316_f64,  &
      -3.278822290154918_f64,       0.91934237457174_f64,   &
       2.422370939448828_f64,      -0.56276020544715_f64, &
      -1.868722023313087E-002_f64, -4.55371557572595_f64, &
       2.762138473069253_f64,       2.51683394103193_f64, &
       1.283272220918888_f64,       0.45132746301539_f64, &
       1.937776139926337_f64,       1.19271986406577_f64, &
       0.845385871848275_f64,       1.04595702707113_f64, &
       1.058484894254336_f64,       0.82797925560061_f64, &
       1.036069666844840_f64,       0.87590353474025_f64  ], [num_cells,2])

  error = maxval(efield-efield_ref)
  if (error> EQV_TOL) then
     passed = .FALSE.
     print*, 'e field wrong in advect_e.', error
  end if
  
  ! Compare to reference
  ! Particle information after advect_vb application 
  particle_info_check(:,1) = [11.780972450961723_f64,  5.497787143782138_f64]
  particle_info_check(:,2) = [-1.449641955904839_f64,  0.182730900323353_f64]
  particle_info_check(:,3) = [ 0.253341837060126_f64, -1.405566770677020_f64]
  particle_info_check(:,4) = [ 6.863675937607472_f64,  5.702694676751700_f64]
  ! Compare computed values to reference values
  do i_part=1,n_particles
     xi = particle_group%group(1)%get_x(i_part)
     if (abs(xi(1)-particle_info_check(i_part,1))> EQV_TOL) then
        passed = .FALSE.
     end if
     xi = particle_group%group(1)%get_v(i_part)
     if (abs(xi(1)-particle_info_check(i_part,2))> EQV_TOL) then
        passed = .FALSE.
     elseif (abs(xi(2)-particle_info_check(i_part,3))> EQV_TOL) then
        passed = .FALSE.
     end if
     xi(1:1) = particle_group%group(1)%get_charge(i_part)
     if (abs(xi(1)-particle_info_check(i_part,4))> EQV_TOL) then
        passed = .FALSE.
     end if
  end do
  
  if (passed .EQV. .FALSE.) then
     print*, 'V wrong in advect_e.'
  end if


  ! Reset particle info
  ! Initialize particles from particle_info_ref
  xi = 0.0_f64
  do i_part =1, n_particles
     xi(1) = particle_info_ref(i_part, 1) 
     call particle_group%group(1)%set_x(i_part, xi)
     xi(1:2) = particle_info_ref(i_part, 2:3)
     call particle_group%group(1)%set_v(i_part, xi)
     xi(1) = particle_info_ref(i_part, 4)
     call particle_group%group(1)%set_weights(i_part, xi(1))
  end do  
  
  call propagator%advect_eb( delta_t*5.0_f64 )
  
  bfield_ref = [ 0.890457153271563_f64,        1.2430359710432_f64, &
       0.5135863209832091_f64,        1.230611521372678_f64, &
       1.1475625443357698_f64,         0.926940887416955_f64, &
       1.0010906711442236_f64,         1.072527573576882_f64, &
       0.9418813699653678_f64,        1.032305986890148_f64 ] 


  efield_ref = reshape( [  0.316277896866717_f64,       0.214721759513144_f64, &
      -3.278822290154864_f64,       0.919342374571724_f64, &
       2.422370939448880_f64,      -0.562760205447174_f64, &
      -1.868722023307678E-002_f64, -4.553715575725981_f64, &
       2.762138473069307_f64,       2.516833941031918_f64, &
       1.113406335317590_f64,       0.723719059320193_f64, &
       1.682252207377827_f64,       1.268128545033283_f64, &
       0.873732288883408_f64,       1.040396287849270_f64, &
       1.022386109540600_f64,       0.888328400462230_f64, &
       0.972374087064816_f64,       0.970152617438486_f64 ], [num_cells,2])

  error = maxval(efield-efield_ref)
  if (error> EQV_TOL) then
     passed = .FALSE.
     print*, 'e field wrong in advect_eb.', error
  end if
  error = maxval(bfield-bfield_ref)

  if (error> EQV_TOL) then
     passed = .FALSE.
     print*, 'b field wrong in advect_eb.', error
  end if
  
  

  if (passed .EQV. .TRUE.) then
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

end program test_time_propagator_pic_1d2v_vm_disgradE
