! TODO: Use input from file to initialize and compare
! Unit test for antisymmetric splitting with coordinate transformation
! author: Benedikt Perse, IPP
program test_time_propagator_pic_1d2v_vm_trafo
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
    sll_p_pi, sll_p_twopi

  use sll_m_time_propagator_pic_vm_1d2v_trafo, only: &
    sll_t_time_propagator_pic_vm_1d2v_trafo

  use sll_m_particle_mesh_coupling_base_1d, only: &
    sll_p_galerkin, &
    sll_c_particle_mesh_coupling_1d

  use sll_m_particle_mesh_coupling_spline_1d, only: &
    sll_t_particle_mesh_coupling_spline_1d, &
    sll_s_new_particle_mesh_coupling_spline_1d_ptr

  use sll_m_maxwell_1d_base, only: &
    sll_c_maxwell_1d_base

  use sll_m_maxwell_1d_trafo, only: &
    sll_t_maxwell_1d_trafo

  use sll_m_particle_group_1d2v, only: &
    sll_t_particle_group_1d2v

  use sll_m_particle_group_base, only: &
      sll_t_particle_array

    use sll_m_3d_coordinate_transformations

  use sll_m_mapping_3d, only: &
       sll_t_mapping_3d
  

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! Tolerance for comparison of real numbers: set it here!
  sll_real64, parameter :: EQV_TOL = 1.5e-11_f64

  ! Abstract particle group
  class(sll_t_particle_array), pointer :: particle_group
  class(sll_t_particle_group_1d2v), pointer :: pg

  ! Arrays for the fields
  sll_real64, pointer :: efield(:,:), efield_ref(:,:)
  sll_real64, pointer :: bfield(:), bfield_ref(:)
  sll_real64, pointer :: rho(:), rho_local(:)

  ! Abstract kernel smoothers
  class(sll_c_particle_mesh_coupling_1d), pointer :: kernel_smoother_0     
  class(sll_c_particle_mesh_coupling_1d), pointer :: kernel_smoother_1
  
  ! Maxwell solver 
  ! Abstract 
  class(sll_c_maxwell_1d_base), pointer :: maxwell_solver
  
  ! Specific Hamiltonian splitting
  type(sll_t_time_propagator_pic_vm_1d2v_trafo) :: propagator

  type(sll_t_mapping_3d), pointer :: map  
  
  ! Parameters
  sll_int32  :: n_particles
  sll_real64 :: eta_min, eta_max, domain(3)
  sll_int32  :: num_cells
  sll_real64 :: delta_t
  sll_int32  :: degree_smoother
  sll_int64  :: rnd_seed
  sll_real64 :: params(6)
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
  eta_max = 1._f64
  num_cells = 10
  delta_t = 0.1_f64
  degree_smoother = 3
  passed = .TRUE.
  rnd_seed = 10



  !Initialize mapping
  params=0._f64
  params(1)= 2._f64*sll_p_twopi
  params(2)=1._f64
  params(3)=1._f64
  params(4)=0.1_f64

  allocate(map)
  call map%init( params,&
      sll_f_orthogonal_x1,&
      sll_f_orthogonal_x2,&
      sll_f_orthogonal_x3,&
      sll_f_orthogonal_jac11,&
      sll_f_orthogonal_jac12,&
      sll_f_orthogonal_jac13,&
      sll_f_orthogonal_jac21,&
      sll_f_orthogonal_jac22,&
      sll_f_orthogonal_jac23,&
      sll_f_orthogonal_jac31,&
      sll_f_orthogonal_jac32,&
      sll_f_orthogonal_jac33,&
      sll_f_orthogonal_jacobian)

  
  domain = [eta_min, eta_max, eta_max - eta_min]
  
 ! Initialize
  allocate(pg)
  call pg%init(n_particles, &
       n_particles ,1.0_f64, 1.0_f64, 1)
  allocate( particle_group)
  particle_group%n_species = 1
  allocate( particle_group%group(1), source=pg )

  call particle_group%group(1)%set_common_weight (1.0_f64)


  SLL_ALLOCATE(particle_info_ref(n_particles,4), i_part)
  SLL_ALLOCATE(particle_info_check(n_particles,4), i_part)
  particle_info_ref = 0.0_f64
  particle_info_ref = reshape([0.9375_f64,        0.4375_f64,    &
      -1.534120544352545_f64,        0.157310684610170_f64,    &
       0.157310684610170_f64,       -1.534120544352545_f64,    &
       6.863675937607472_f64,        5.702694676751700_f64], [n_particles, 4])

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
  allocate( sll_t_maxwell_1d_trafo :: maxwell_solver )
  select type ( maxwell_solver )
  type is ( sll_t_maxwell_1d_trafo )
    call maxwell_solver%init( [eta_min, eta_min+params(1)], num_cells, degree_smoother, map, solver_tolerance = 1d-14)
  end select
  
  SLL_ALLOCATE(efield(kernel_smoother_0%n_dofs,2),ierr)
  SLL_ALLOCATE(bfield(kernel_smoother_0%n_dofs),ierr)
  SLL_ALLOCATE(efield_ref(kernel_smoother_0%n_dofs,2),ierr)
  SLL_ALLOCATE(bfield_ref(kernel_smoother_0%n_dofs),ierr)
  SLL_ALLOCATE(rho(kernel_smoother_0%n_dofs),ierr)
  SLL_ALLOCATE(rho_local(kernel_smoother_0%n_dofs),ierr)

  efield(:,2) = 1.0_f64

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
  call maxwell_solver%compute_E_from_rho( rho, efield(:,1))
  bfield = params(1)

  call propagator%init( maxwell_solver, &
       kernel_smoother_0, kernel_smoother_1, particle_group, &
       efield, bfield, &
       eta_min, eta_max-eta_min, map, solver_tolerance = 1d-14)

  call propagator%advect_x( delta_t )

  ! Compare to reference
  ! Particle information after advect_x application 
  particle_info_check(:,1) = [11.146655277439535_f64,  5.994411753754585_f64]
  particle_info_check(:,2) = [-1.534120544352545_f64,  0.157310684610170_f64]
  particle_info_check(:,3) = [ 0.157310684610170_f64, -1.534120544352545_f64]
  particle_info_check(:,4) = [ 6.863675937607472_f64,  5.702694676751700_f64]
  ! Compare computed values to reference values
  do i_part=1,n_particles
    xi = particle_group%group(1)%get_x(i_part)
    xi(1)=map%get_x1([xi(1), 0._f64, 0._f64])
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
  particle_info_check(:,1) = [11.300078267054246_f64,  5.97868132768961_f64]
  particle_info_check(:,2) = [-1.521104168249398_f64, -0.20937811683683_f64]
  particle_info_check(:,3) = [ 0.253997254771338_f64, -1.52788523793082_f64]
  particle_info_check(:,4) = [ 6.863675937607472_f64,  5.7026946767517_f64]
  ! Compare computed values to reference values
  do i_part=1,n_particles
    xi = particle_group%group(1)%get_x(i_part)
    xi(1)=map%get_x1([xi(1), 0._f64, 0._f64])
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

  efield_ref = reshape([ 3.52695203249949_f64,        -1.7744293184851_f64, &
     -24.95823230104683_f64,        5.152551580890536_f64, &
      12.05269842145767_f64,       -3.880285811847831_f64, &
       1.417736642890464_f64,      -77.23828211414418_f64, &
      40.68624422581577_f64,        54.41638736872317_f64, &
       1.366286390402017_f64,        0.188037456501155_f64, &
       2.657153387208462_f64,        1.806767239767959_f64, &
       0.530247250362130_f64,        1.166149515026346_f64, &
       1.000163593159168_f64,        0.884746280866144_f64, &
       1.057369533431829_f64,        0.842109126367992_f64 ], [num_cells,2])
  error = maxval(abs(efield-efield_ref))
  if (error> EQV_TOL) then
     passed = .FALSE.
     print*, 'e field wrong in advect_e.', error
  end if
  
  ! Compare to reference
  ! Particle information after advect_vb application 
  particle_info_check(:,1) = [11.300078267054246_f64,   5.978681327689614_f64]
  particle_info_check(:,2) = [-1.470221580781017_f64,   0.163231181152982_f64]
  particle_info_check(:,3) = [ 0.255200247982141_f64,  -1.374148425372217_f64]
  particle_info_check(:,4) = [ 6.863675937607472_f64,   5.702694676751700_f64]
  ! Compare computed values to reference values
  do i_part=1,n_particles
    xi = particle_group%group(1)%get_x(i_part)
    xi(1)=map%get_x1([xi(1), 0._f64, 0._f64])
    if (abs(xi(1)-particle_info_check(i_part,1))> EQV_TOL) then
       passed = .FALSE.
     end if
     xi = particle_group%group(1)%get_v(i_part)
     if (abs(xi(1)-particle_info_check(i_part,2))> EQV_TOL) then
       passed = .FALSE.
       print*,'error in v1', abs(xi(1)-particle_info_check(i_part,2))
     elseif (abs(xi(2)-particle_info_check(i_part,3))> EQV_TOL) then
       passed = .FALSE.
       print*,'error in v2', abs(xi(2)-particle_info_check(i_part,3))
     end if
     xi(1:1) = particle_group%group(1)%get_charge(i_part)
     if (abs(xi(1)-particle_info_check(i_part,4))> EQV_TOL) then
       passed = .FALSE.
     end if
  end do
  
  if (passed .EQV. .FALSE.) then
     print*, 'Error in advect_ev.'
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
  
  bfield_ref = [ 11.29434057882151_f64,        16.22003743893807_f64, &
        5.60930435013833_f64,        13.32541348196329_f64, &
       14.55855274357444_f64,        13.15087343610345_f64, &
       12.95473215619471_f64,        13.45907118715099_f64, &
       11.93626226903481_f64,        13.15511850167207_f64 ]  
  
  efield_ref = reshape( [ 3.52695203249949_f64,       -1.7744293184851_f64, &
     -24.95823230104683_f64,        5.152551580890536_f64, &
      12.05269842145767_f64,       -3.880285811847831_f64, &
       1.41773664289046_f64,      -77.23828211414418_f64, &
      40.68624422581577_f64,       54.41638736872317_f64, &
       0.956974444676979_f64,       0.673756648746287_f64, &
       0.987467223727320_f64,       1.534236224126173_f64, &
       2.013883361845894_f64,       1.144179968483959_f64, &
       1.154821273616908_f64,       0.913158356793186_f64, &
       0.992578442357238_f64,       0.972339694495922_f64 ], [num_cells,2])

  error = maxval(abs(efield-efield_ref))
  if (error> EQV_TOL) then
     passed = .FALSE.
     print*, 'e field wrong in advect_eb.', error
  end if
  error = maxval(abs(bfield-bfield_ref))

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

end program test_time_propagator_pic_1d2v_vm_trafo
