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
  particle_info_ref = reshape([11.780972450961723_f64,        5.4977871437821380_f64,       -1.5341205443525459_f64,       0.15731068461017067_f64,       0.15731068461017067_f64,       -1.5341205443525459_f64,        6.8636759376074723_f64,        5.7026946767517002_f64   ], [n_particles, 4])

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
  particle_info_check(:,1) = [11.627560396526469_f64, 5.5135182122431550_f64]
  particle_info_check(:,2) = [-1.5341205443525459_f64, 0.15731068461017067_f64]
  particle_info_check(:,3) = [0.15731068461017067_f64, -1.5341205443525459_f64]
  particle_info_check(:,4) = [6.8636759376074723_f64,  5.7026946767517002_f64]
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
  particle_info_check(:,1) = [11.780972450961723_f64,        5.4977871437821380_f64]
  particle_info_check(:,2) = [-1.5107514685496906_f64, 3.3682909390511384E-003_f64]
  particle_info_check(:,3) = [ 0.30968128192066413_f64, -1.5421611947890606_f64]
  particle_info_check(:,4) = [6.8636759376074723_f64,  5.7026946767517002_f64]
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

  efield_ref = reshape([ 0.31627789686666258_f64,       0.21472175951316622_f64,  &
       -3.2788222901549182_f64,       0.91934237457174484_f64,   &
       2.4223709394488284_f64,      -0.56276020544715444_f64, &
       -1.8687220233130877E-002_f64,  -4.5537155757259580_f64, &
       2.7621384730692538_f64,        2.5168339410319378_f64, &
       1.2832722209188887_f64,       0.45132746301539584_f64, &
       1.9377761399263371_f64,        1.1927198640657763_f64, &
       0.84538587184827541_f64,        1.0459570270711387_f64, &
       1.0584848942543363_f64,       0.82797925560061347_f64, &
       1.0360696668448408_f64,       0.87590353474025973_f64  ], [num_cells,2])

  error = maxval(efield-efield_ref)
  if (error> EQV_TOL) then
     passed = .FALSE.
     print*, 'e field wrong in advect_e.', error
  end if
  
  ! Compare to reference
  ! Particle information after advect_vb application 
  particle_info_check(:,1) = [11.780972450961723_f64,        5.4977871437821380_f64]
  particle_info_check(:,2) = [-1.4496419559048399_f64, 0.18273090032335398_f64]
  particle_info_check(:,3) = [ 0.25334183706012608_f64, -1.4055667706770207_f64]
  particle_info_check(:,4) = [6.8636759376074723_f64,  5.7026946767517002_f64]
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
  
  bfield_ref = [ 0.890457153271563946_f64,        1.24303597104320026_f64, &
       0.513586320983209116_f64,        1.23061152137267896_f64, &
       1.14756254433576986_f64,       0.926940887416955916_f64, &
       1.00109067114422366_f64,        1.07252757357688226_f64, &
       0.941881369965367876_f64,        1.03230598689014856_f64 ] 

!!$  [ 0.96816309794804756_f64,         1.0652888181806115_f64, &
!!$       0.88277345057999335_f64,         1.0585999213799744_f64, &
!!$       1.0277955028243342_f64,        0.98405636336018920_f64, &
!!$       0.99909146342661515_f64,        1.0181286601865549_f64, &
!!$       0.98372116282938227_f64,        1.0123815592842973_f64]

  efield_ref = reshape( [  0.31627789686671726_f64,       0.21472175951314446_f64, &
       -3.2788222901548640_f64,       0.91934237457172441_f64, &
       2.4223709394488808_f64,      -0.56276020544717498_f64, &
       -1.8687220233076788E-002_f64,  -4.5537155757259811_f64, &
       2.7621384730693075_f64,        2.5168339410319187_f64, &
       1.1134063353175903_f64,       0.72371905932019331_f64, &
       1.6822522073778270_f64,        1.2681285450332833_f64, &
       0.87373228888340826_f64,        1.0403962878492701_f64, &
       1.0223861095406002_f64,       0.88832840046223072_f64, &
       0.97237408706481698_f64,       0.97015261743848635_f64 ], [num_cells,2])

!!$  [ 0.31627789686666258_f64,       0.21472175951316622_f64, &
!!$       -3.2788222901549182_f64,       0.91934237457174484_f64, &
!!$       2.4223709394488284_f64,      -0.56276020544715444_f64, &
!!$       -1.8687220233130877E-002_f64,  -4.5537155757259580_f64, &
!!$       2.7621384730692538_f64,        2.5168339410319378_f64, &
!!$       1.2741031034377539_f64,       0.46516088887842894_f64, &
!!$       1.9249367436768068_f64,        1.1972163594702570_f64, &
!!$       0.84597317188298016_f64,        1.0461093105715515_f64, &
!!$       1.0564154579858682_f64,       0.83129817134579298_f64, &
!!$       1.0323395622166271_f64,       0.88132316881979644_f64]

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
