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
  sll_real64, parameter :: EQV_TOL = 1.5e-12_f64

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
  particle_info_ref = reshape([0.937500000000000_f64,        0.437500000000000_f64,       -1.5341205443525459_f64,       0.15731068461017067_f64,       0.15731068461017067_f64,       -1.5341205443525459_f64,        6.8636759376074723_f64,        5.7026946767517002_f64], [n_particles, 4])

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
  particle_info_check(:,1) = [11.146655277439535_f64, 5.9944117537545853_f64]
  particle_info_check(:,2) = [-1.5341205443525459_f64, 0.15731068461017067_f64]
  particle_info_check(:,3) = [0.15731068461017067_f64, -1.5341205443525459_f64]
  particle_info_check(:,4) = [6.8636759376074723_f64,  5.7026946767517002_f64]
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
  particle_info_check(:,1) = [11.300078267054246_f64,  5.9786813276896149_f64]
  particle_info_check(:,2) = [-1.5211041682493989_f64, -0.20937811683683141_f64]
  particle_info_check(:,3) = [0.25399725477133861_f64, -1.5278852379308252_f64]
  particle_info_check(:,4) = [6.8636759376074723_f64,  5.7026946767517002_f64]
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

  efield_ref = reshape([ 3.5269520324994938_f64,        -1.7744293184851003_f64, &
       -24.958232301046834_f64,        5.1525515808905364_f64, &
       12.052698421457679_f64,        -3.8802858118478318_f64, &
       1.4177366428904643_f64,        -77.238282114144184_f64, &
       40.686244225815770_f64,         54.416387368723171_f64, &
       1.3662863904020170_f64,        0.18803745650115597_f64, &
       2.6571533872084623_f64,         1.8067672397679591_f64, &
       0.53024725036213016_f64,        1.1661495150263466_f64, &
       1.0001635931591684_f64,        0.88474628086614449_f64, &
       1.0573695334318294_f64,        0.84210912636799262_f64 ], [num_cells,2])
  error = maxval(abs(efield-efield_ref))
  if (error> EQV_TOL) then
     passed = .FALSE.
     print*, 'e field wrong in advect_e.', error
  end if
  
  ! Compare to reference
  ! Particle information after advect_vb application 
  particle_info_check(:,1) = [ 11.300078267054246_f64,    5.9786813276896149_f64]
  particle_info_check(:,2) = [-1.4702215807810171_f64,   0.16323118115298291_f64]
  particle_info_check(:,3) = [ 0.25520024798214108_f64,  -1.3741484253722172_f64]
  particle_info_check(:,4) = [ 6.8636759376074723_f64,   5.7026946767517002_f64]
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
  
  bfield_ref = [ 11.294340578821515_f64,        16.220037438938071_f64, &
       5.6093043501383377_f64,        13.325413481963295_f64, &
       14.558552743574444_f64,        13.150873436103453_f64, &
       12.954732156194719_f64,        13.459071187150995_f64, &
       11.936262269034817_f64,        13.155118501672078_f64 ]  
  
  efield_ref = reshape( [ 3.5269520324994938_f64,       -1.7744293184851003_f64, &
       -24.958232301046834_f64,       5.1525515808905364_f64, &
       12.052698421457679_f64,       -3.8802858118478318_f64, &
       1.4177366428904643_f64,       -77.238282114144184_f64, &
       40.686244225815770_f64,        54.416387368723171_f64, &
       0.95697444467697945_f64,      0.67375664874628727_f64, &
       0.98746722372732010_f64,       1.5342362241261738_f64, &
       2.0138833618458940_f64,        1.1441799684839593_f64, &
       1.1548212736169088_f64,       0.91315835679318624_f64, &
       0.99257844235723891_f64,       0.97233969449592283_f64 ], [num_cells,2])

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
