! TODO: Use input from file to initialize and compare
! Unit test for symplectic splitting
! author: Katharina Kormann, IPP
program test_time_propagator_pic_vm_1d2v_hs
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

  use sll_m_time_propagator_pic_vm_1d2v_hs, only: &
    sll_t_time_propagator_pic_vm_1d2v_hs

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
  sll_real64, parameter :: EQV_TOL = 1.0e-13_f64

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

  ! Binomial filter object
  type( sll_t_binomial_filter), target :: filter
  
  ! Specific Hamiltonian splitting
  type(sll_t_time_propagator_pic_vm_1d2v_hs) :: propagator
  
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
  particle_info_ref = reshape([11.78097245096172_f64,        &
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

  call filter%init( 0, 0 )

  call propagator%init( maxwell_solver, &
       kernel_smoother_0, kernel_smoother_1, particle_group, &
       phi, efield, bfield, &
       eta_min, eta_max-eta_min, filter)

  propagator%adiabatic_electrons = .false.

  call propagator%operatorHp1(delta_t)


  ! Compare to reference
  ! Particle information after operatorV application 
  particle_info_ref = reshape([  11.6275603965264_f64,        &
5.51351821224315_f64,       -1.53412054435254_f64,       &
0.15731068461017_f64,       0.310722739045425_f64,       &
-1.54985161281356_f64,        6.86367593760747_f64,        &
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
     print*, 'Error in operatorHp1.'
  end if
  
  call propagator%operatorHp2(delta_t)

  ! Compare to reference
  ! Particle information after operatorV application 
  particle_info_ref = reshape([ 11.6275603965264_f64,        &
5.51351821224315_f64,       -1.503048270448_f64,       &
 2.32552332881433D-003,  0.310722739045425_f64,       &
-1.54985161281356_f64,        6.86367593760747_f64,   &
     5.7026946767517_f64   ], [n_particles, 4])
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
     print*, 'Error in operatorHp2.'
  end if

  !print*, 'v(1)', size(pg%particle_array,1), size(pg%particle_array,2),pg%particle_array(2,:)

  call propagator%operatorHE(delta_t)

  ! Compare to reference
  ! Particle information after operatorV application 
  particle_info_ref = reshape([    11.6275603965264_f64,  &
      5.51351821224315_f64,       -1.43114619603643_f64,&
       2.76113371384939E-002_f64,  0.395114887882236_f64,  &
     -1.39031318482687_f64,        6.86367593760747_f64,   &
     5.7026946767517_f64    ], [n_particles, 4])
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
  particle_info_ref = reshape([   11.6275603965264_f64,    &
    5.51351821224315_f64,       -1.431146196036438_f64,   &
    2.76113371384939E-002_f64,  0.395114887882236_f64,    &
   -1.39031318482687_f64,       6.86367593760747_f64,     &
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
     print*, 'Error in operatorHB.'
  end if


  bfield_ref = [ 0.962414253193435_f64,        1.070484285773022_f64, &
       0.877494708070806_f64,        1.0585080035598995_f64, &
       1.03120678997331_f64,         0.9852448637249098_f64, &
       0.991638668895902_f64,        1.0359265946162064_f64, &
       0.965762298883728_f64,        1.0213195333087741_f64 ];

  efield_ref = reshape([ 0.32668967300827_f64,  0.21111816665325_f64, &
      -3.279735301413320_f64,       0.9318205986970499_f64, &
       2.409952052992651_f64,      -0.5407493530167542_f64, &
      -6.650814373435558E-002_f64, -4.4526750013638647_f64, &
       2.731401214003940_f64,        2.495224958698765_f64, &
       1.289653877027072_f64,        0.454284284708069_f64, &
       1.937203186631420_f64,        1.238094676549730_f64, &
       0.837327955152277_f64,        1.024488521395515_f64, &
       1.118798797824109_f64,        0.687893072962867_f64, &
       1.093857657675167_f64,        0.852015078837940_f64 ], [num_cells,2]);


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

end program test_time_propagator_pic_vm_1d2v_hs
