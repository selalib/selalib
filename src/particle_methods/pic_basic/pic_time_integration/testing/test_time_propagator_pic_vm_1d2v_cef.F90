! TODO: Use input from file to initialize and compare
! Unit test for symplectic splitting
! author: Katharina Kormann, IPP
program test_time_propagator_pic_1d2v_vm_cef
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

  use mpi, only: &
       MPI_SUM

  use sll_m_collective, only: &
    sll_s_boot_collective, &
    sll_o_collective_allreduce, &
    sll_s_halt_collective, &
    sll_v_world_collective

  use sll_m_constants, only: &
    sll_p_pi

  use sll_m_time_propagator_pic_vm_1d2v_cef, only: &
    sll_t_time_propagator_pic_vm_1d2v_cef

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
  sll_real64, parameter :: EQV_TOL = 1.0e-12_f64

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
  type(sll_t_time_propagator_pic_vm_1d2v_cef) :: propagator
  
  ! Parameters
  sll_int32  :: n_particles
  sll_real64 :: eta_min, eta_max, domain(3)
  sll_int32  :: num_cells
  sll_real64 :: delta_t
  sll_int32  :: degree_smoother
  sll_int64  :: rnd_seed
  sll_real64 :: charge

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

  charge = 1.0_f64 ! Can also be changed to -1.0_f64
  
  ! Initialize
  allocate( pg )
  call pg%init(n_particles, &
       n_particles , charge, 1.0_f64, 1)
  !print*, 'test1',  pg%n_particles
  allocate( particle_group )
  particle_group%n_species = 1

  ! Alternative allocation if we do not need the pg explicitly
  !allocate( sll_t_particle_group_1d2v :: particle_group%group(1) )
  !particle_group%group(1) = cpg
  
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
  particle_info_ref = 0.0_f64
  particle_info_ref = reshape([11.780972450961723_f64,        &
       5.497787143782138_f64,       -1.534120544352545_f64,       &
       0.1573106846101706_f64,       0.1573106846101706_f64,       &
      -1.534120544352545_f64,        6.863675937607472_f64,        &
       5.7026946767517_f64   ], [n_particles, 4])

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
          degree_smoother, delta_t )
  end select
  SLL_ALLOCATE(phi(kernel_smoother_0%n_dofs),ierr)
  SLL_ALLOCATE(efield(kernel_smoother_0%n_dofs,2),ierr)
  SLL_ALLOCATE(bfield(kernel_smoother_0%n_dofs),ierr)
  SLL_ALLOCATE(efield_ref(kernel_smoother_0%n_dofs,2),ierr)
  SLL_ALLOCATE(bfield_ref(kernel_smoother_0%n_dofs),ierr)
  SLL_ALLOCATE(rho(kernel_smoother_0%n_dofs),ierr)
  SLL_ALLOCATE(rho_local(kernel_smoother_0%n_dofs),ierr)

  efield = charge*1.0_f64

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
  call maxwell_solver%compute_E_from_rho( rho,efield(:,1) )
  bfield = charge*1.0_f64

  call propagator%init( maxwell_solver, &
       kernel_smoother_0, kernel_smoother_1, particle_group, &
       efield, bfield, &
       eta_min, eta_max-eta_min)

  !propagator%adiabatic_electrons = .false.

  call propagator%operatorHp(delta_t)


  ! Compare to reference
  ! Particle information after operatorV application 
  particle_info_ref = reshape([  11.62756039652646_f64,        &
       5.513518212243155_f64,       -1.534120544352545_f64,       &
       0.157310684610170_f64,        0.15731068461017059_f64,       &
      -1.5341205443525450_f64,        6.863675937607472_f64,        &
       5.7026946767517_f64    ], [n_particles, 4])
  ! Compare computed values to reference values
  do i_part=1,n_particles
     xi = particle_group%group(1)%get_x(i_part)
     if (abs(xi(1)-particle_info_ref(i_part,1))> EQV_TOL) then
        print*, 'x', xi(1), particle_info_ref(i_part,1)
        passed = .FALSE.
     end if
     xi = particle_group%group(1)%get_v(i_part)
     if (abs(xi(1)-particle_info_ref(i_part,2))> EQV_TOL) then
        print*, 'v1', xi(1), particle_info_ref(i_part,2)
        passed = .FALSE.
     elseif (abs(xi(2)-particle_info_ref(i_part,3))> EQV_TOL) then
        print*, 'v2', xi(2), particle_info_ref(i_part,3)
        passed = .FALSE.
     end if
     xi(1:1) = particle_group%group(1)%get_charge(i_part)
     if (abs(xi(1)-charge*particle_info_ref(i_part,4))> EQV_TOL) then
        print*, 'charge', xi(1), particle_info_ref(i_part,4)        
        passed = .FALSE.
     end if
  end do
  
  if (passed .EQV. .FALSE.) then
     print*, 'Error in operatorHp.'
  end if
  
  call propagator%operatorHE(delta_t)

  ! Compare to reference
  ! Particle information after operatorV application 
  particle_info_ref = reshape([    11.627560396526469_f64,  &
      5.513518212243155_f64,       -1.462218469940975_f64,&
      0.18259649841985_f64,         0.250325316116684_f64,  &
     -1.374927208519881_f64,        6.863675937607472_f64,   &
      5.702694676751700_f64    ], [n_particles, 4])
  ! Compare computed values to reference values
    do i_part=1,n_particles
     xi = particle_group%group(1)%get_x(i_part)
     if (abs(xi(1)-particle_info_ref(i_part,1))> EQV_TOL) then
        passed = .FALSE.
        print*, 'x(1) of particle ', i_part,' wrong'
     end if
     xi = particle_group%group(1)%get_v(i_part)
     if (abs(xi(1)-particle_info_ref(i_part,2))> EQV_TOL) then
        passed = .FALSE.
        print*, 'v(1) of particle ', i_part,' wrong'
        print*, i_part, xi(1), particle_info_ref(i_part,2)
     elseif (abs(xi(2)-particle_info_ref(i_part,3))> EQV_TOL) then
        passed = .FALSE.
        print*, 'v(2) of particle ', i_part,' wrong'
        print*,   i_part, xi(2), particle_info_ref(i_part,3)
     end if
     xi(1:1) = particle_group%group(1)%get_charge(i_part)
     if (abs(xi(1)-charge*particle_info_ref(i_part,4))> EQV_TOL) then
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
  particle_info_ref = reshape([ 11.62756039652646_f64,    &
    5.51351821224315_f64,       -1.43021010930802_f64,   &
    4.12326594633539E-002_f64,   0.39401099091087_f64,    &
   -1.38638601325786_f64,        6.86367593760747_f64,     &
   5.7026946767517_f64      ], [n_particles, 4])
  ! Compare computed values to reference values
    do i_part=1,n_particles
     xi = particle_group%group(1)%get_x(i_part)
     if (abs(xi(1)-particle_info_ref(i_part,1))> EQV_TOL) then
        passed = .FALSE.
        print*, 'x(1) of particle ', i_part,' wrong'
     end if
     xi = particle_group%group(1)%get_v(i_part)
     if (abs(xi(1)-particle_info_ref(i_part,2))> EQV_TOL) then
        passed = .FALSE.
        print*, 'v(1) of particle ', i_part,' wrong', xi(1)
        print*, i_part, xi(1), particle_info_ref(i_part,2)
     elseif (abs(xi(2)-particle_info_ref(i_part,3))> EQV_TOL) then
        passed = .FALSE.
        print*, 'v(2) of particle ', i_part,' wrong'
        print*,   i_part, xi(2), particle_info_ref(i_part,3)
     end if
     xi(1:1) = particle_group%group(1)%get_charge(i_part)
     if (abs(xi(1)-charge*particle_info_ref(i_part,4))> EQV_TOL) then
        passed = .FALSE.
        print*, 'weight of particle ', i_part,' wrong'
     end if
  end do

  if (passed .EQV. .FALSE.) then
     print*, 'Error in operatorHB.'
  end if


  bfield_ref = [ 0.9640913345584765_f64, &
       1.069501357174962_f64,        0.8779765545703390_f64, &
       1.058999200145200_f64,        1.0316973600817902_f64, &
       0.980877334067277_f64,        1.0018819117882489_f64, &
       1.015391028305587_f64,        0.9823914733157797_f64, &
       1.017192445992337_f64 ];

  efield_ref = reshape([ 0.326689673008277_f64, &
       0.211118166653252_f64,       -3.279735301413321_f64, &
       0.931820598697050_f64,        2.409952052992651_f64, &
      -0.540749353016754_f64,       -6.650814373435491E-002_f64, &
      -4.452675001363864_f64,        2.731401214003940_f64, &
       2.495224958698765_f64,        1.283297305601403_f64, &
       0.458638076711317_f64,        1.936607469843452_f64, &
       1.230253464136155_f64,        0.825389850287529_f64, &
       1.063026015198084_f64,        1.036525537754838_f64, &
       0.853098891024324_f64,        1.059464793580823_f64, &
       0.863969591078809_f64], [num_cells,2]);


  error = maxval(bfield-charge*bfield_ref)
  
  if (error> EQV_TOL) then
     passed = .FALSE.
     print*, 'b field wrong.'
  end if

  error = maxval(efield-charge*efield_ref)

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

end program test_time_propagator_pic_1d2v_vm_cef

