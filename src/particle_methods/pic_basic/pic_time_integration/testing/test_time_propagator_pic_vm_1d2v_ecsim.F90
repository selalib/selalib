! TODO: Use input from file to initialize and compare
! Unit test for symplectic splitting
! author: Benedikt Perse
program test_time_propagator_pic_1d2v_vm_ecsim
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

  use sll_m_time_propagator_pic_vm_1d2v_ecsim, only: &
    sll_t_time_propagator_pic_vm_1d2v_ecsim

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
  sll_real64, parameter :: EQV_TOL = 1.5e-13_f64

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
  sll_real64 :: solver_tolerance = 1d-13
  ! Specific Hamiltonian splitting
  type(sll_t_time_propagator_pic_vm_1d2v_ecsim) :: propagator
  
  ! Parameters
  sll_int32  :: n_particles
  sll_real64 :: eta_min, eta_max, domain(3)
  sll_int32  :: num_cells
  sll_real64 :: delta_t
  sll_int32  :: degree_smoother
  sll_int64  :: rnd_seed

  ! Helper 
  sll_int32  :: i_part, i_sp
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
  allocate( particle_group )
  particle_group%n_species = 1
  i_sp = 1

  
  allocate( particle_group%group(1), source=pg )

  call particle_group%group(i_sp)%set_common_weight (1.0_f64)


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
     call particle_group%group(i_sp)%set_x(i_part, xi)
     xi(1:2) = particle_info_ref(i_part, 2:3)
     call particle_group%group(i_sp)%set_v(i_part, xi)
     xi(1) = particle_info_ref(i_part, 4)
     call particle_group%group(i_sp)%set_weights(i_part, xi(1))
  end do

  call particle_group%group(i_sp)%set_common_weight (1.0_f64)

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
  
  SLL_ALLOCATE(efield(kernel_smoother_0%n_dofs,2),ierr)
  SLL_ALLOCATE(bfield(kernel_smoother_0%n_dofs),ierr)
  SLL_ALLOCATE(efield_ref(kernel_smoother_0%n_dofs,2),ierr)
  SLL_ALLOCATE(bfield_ref(kernel_smoother_0%n_dofs),ierr)
  SLL_ALLOCATE(rho(kernel_smoother_0%n_dofs),ierr)
  SLL_ALLOCATE(rho_local(kernel_smoother_0%n_dofs),ierr)

  efield = 1.0_f64

  rho_local = 0.0_f64
  do i_part = 1,n_particles
     xi = particle_group%group(i_sp)%get_x(i_part)
     wi(1) = particle_group%group(i_sp)%get_charge( i_part)
     call kernel_smoother_0%add_charge(xi(1), wi(1), rho_local)
  end do
  ! MPI to sum up contributions from each processor
  rho = 0.0_f64
  call sll_o_collective_allreduce( sll_v_world_collective, &
       rho_local, &
       kernel_smoother_0%n_dofs, MPI_SUM, rho)
  ! Solve Poisson problem
  call maxwell_solver%compute_E_from_rho(rho,efield(:,1))
  bfield = 1.0_f64

  call propagator%init( kernel_smoother_0, kernel_smoother_1, &
       particle_group, efield, bfield, &
       eta_min, eta_max-eta_min, solver_tolerance)

 

  call propagator%advect_x( delta_t )

  ! Compare to reference
  ! Particle information after advect_x application 
  particle_info_check(:,1) = [11.627560396526469_f64, 5.5135182122431550_f64]
  particle_info_check(:,2) = [-1.5341205443525459_f64, 0.15731068461017067_f64]
  particle_info_check(:,3) = [0.15731068461017067_f64, -1.5341205443525459_f64]
  particle_info_check(:,4) = [6.8636759376074723_f64,  5.7026946767517002_f64]
  ! Compare computed values to reference values
  do i_part=1,n_particles
     xi = particle_group%group(i_sp)%get_x(i_part)
     if (abs(xi(1)-particle_info_check(i_part,1))> EQV_TOL) then
        passed = .FALSE.
     end if
     xi = particle_group%group(i_sp)%get_v(i_part)
     if (abs(xi(1)-particle_info_check(i_part,2))> EQV_TOL) then
        passed = .FALSE.
     elseif (abs(xi(2)-particle_info_check(i_part,3))> EQV_TOL) then
        passed = .FALSE.
     end if
     xi(1:1) = particle_group%group(i_sp)%get_charge(i_part)
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
     call particle_group%group(i_sp)%set_x(i_part, xi)
     xi(1:2) = particle_info_ref(i_part, 2:3)
     call particle_group%group(i_sp)%set_v(i_part, xi)
     xi(1) = particle_info_ref(i_part, 4)
     call particle_group%group(i_sp)%set_weights(i_part, xi(1))
  end do

  call propagator%advect_e( delta_t )

  

  ! Compare to reference
  efield_ref = reshape([0.31199377839865949_f64,        0.22050392415748257_f64, &
       -3.2916675305308023_f64,       0.97695337617560474_f64, &
       2.4003219706846859_f64,       -0.55299036556468861_f64, &
       -2.3834706718177692E-002_f64,  -4.5490281730908615_f64, &
       2.7480332443923041_f64,         2.5229634028470049_f64, &
       1.2739512649332709_f64,        0.45849001425514457_f64, &
       1.9353003807811671_f64,         1.1915325743511542_f64, &
       0.85415086120891093_f64,        1.0287177029790429_f64, &
       1.0900936409342892_f64,        0.77401542604316043_f64, &
       1.0220769816757111_f64,       0.88886534405692808_f64], [num_cells,2])
  
  error = maxval(efield-efield_ref)

  if (error> EQV_TOL) then
     passed = .FALSE.
     print*, 'e field wrong in advect_e.', error
  end if

  bfield_ref = [0.98467791804435667_f64,        1.0324461722363252_f64,  &
      0.94123958253950857_f64,        1.0295935807264864_f64,  &
       1.0134239918388528_f64,       0.99305420605808870_f64,   &
      0.99755792902188478_f64,        1.0125763525758935_f64,  &
      0.99012994430750567_f64,        1.0053003226510977_f64]


  error = maxval(bfield-bfield_ref)

  if (error> EQV_TOL) then
     passed = .FALSE.
     print*, 'b field wrong in advect_e.', error
  end if
  
  ! Particle information after advect_e application 
  particle_info_check(:,1) = [11.780972450961723_f64,        5.4977871437821380_f64]
  particle_info_check(:,2) = [-1.4222539473208027_f64,  3.7169101803288168E-002_f64 ]
  particle_info_check(:,3) = [0.39946152826921183_f64  , -1.4153643512244618_f64]
  particle_info_check(:,4) = [6.8636759376074723_f64,  5.7026946767517002_f64]
  ! Compare computed values to reference values
  do i_part=1,n_particles
     xi = particle_group%group(i_sp)%get_x(i_part)
     if (abs(xi(1)-particle_info_check(i_part,1))> EQV_TOL) then
        passed = .FALSE.
     end if


     xi = particle_group%group(i_sp)%get_v(i_part)
     if (abs(xi(1)-particle_info_check(i_part,2))> EQV_TOL) then
        passed = .FALSE.
     elseif (abs(xi(2)-particle_info_check(i_part,3))> EQV_TOL) then
        passed = .FALSE.
     end if
      
     xi(1:1) = particle_group%group(i_sp)%get_charge(i_part)
     if (abs(xi(1)-particle_info_check(i_part,4))> EQV_TOL) then
        passed = .FALSE.
     end if
  end do
 
  if (passed .EQV. .FALSE.) then
     print*, 'Error in advect_e.'
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

end program test_time_propagator_pic_1d2v_vm_ecsim
