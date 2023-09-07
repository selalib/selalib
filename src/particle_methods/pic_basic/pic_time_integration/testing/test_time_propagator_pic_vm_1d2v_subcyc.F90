! TODO: Use input from file to initialize and compare
! Unit test for symplectic splitting
! author: Katharina Kormann, IPP
program test_time_propagator_pic_1d2v_vm_subcyc
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

  use mpi, only: &
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

  use sll_m_time_propagator_pic_vm_1d2v_subcyc, only: &
    sll_t_time_propagator_pic_vm_1d2v_subcyc

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
  sll_real64, pointer :: efield(:,:), efield_ref(:,:)
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
  type(sll_t_time_propagator_pic_vm_1d2v_subcyc) :: propagator
  
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
  particle_info_ref = reshape([11.7809724509617_f64,        &
       5.497787143782138_f64,      -1.534120544352545_f64,       &
       0.157310684610170_f64,       0.157310684610170_f64,       &
      -1.534120544352545_f64,       6.863675937607472_f64,        &
       5.702694676751700_f64   ], [n_particles, 4])

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
       efield, bfield, &
       eta_min, eta_max-eta_min, 2, 'test_subcyc_', filter)

  call propagator%operatorall(delta_t)


  ! Compare to reference
  ! Particle information after operatorV application
 ! particle_info_ref = reshape([   11.633445350371955_f64,        &
 !5.5055070035890106_f64,        -1.4660295487739201_f64,       &
!4.1128223050050053E-002_f64,       0.40654006008902749_f64,       &
!-1.4447449656976983_f64,        6.8636759376074723_f64,        &
!5.7026946767517002_f64    ], [n_particles, 4]) ! Reference for Picard
  particle_info_ref = reshape([   11.633445350371_f64,        &
       5.505507003588953_f64,        -1.466029548775424_f64,       &
       4.112822304949182E-002_f64,    0.406540060089496_f64,       &
      -1.444744965698378_f64,         6.863675937607472_f64,        &
       5.702694676751700_f64    ], [n_particles, 4]) ! Reference for Newton
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


  bfield_ref = 1.0_f64
  
 ! efield_ref = reshape([  0.31502302980372776d0,       0.22070774991100911d0,       -3.2957522297314017d0,       0.99607238352733296d0,        2.3856942413645363d0,      -0.53021634675511797d0,       -7.1669522418057852d-002,  -4.4472062513548636d0,        2.6802643557198871d0,        2.5178335384449486d0,        1.2792492342544843d0,       0.46113150101229250d0,        1.9186663491519194d0,        1.1853257272543620d0,       0.86966367446528203d0,       0.99709998122863119d0,        1.1513366492981942d0,       0.65647776138716973d0,        1.0599747631911016d0,       0.87392736806525706d0  ], [num_cells,2]); ! Reference for Picard
  efield_ref = reshape([ 0.315023029803914d0,       0.220707749910959d0,    &
      -3.295752229731461d0,       0.996072383527755d0, &
       2.385694241364380d0,      -0.530216346755056d0,  &
      -7.166952241807099d-002,   -4.447206251354928d0,   &
       2.680264355720997d0,       2.517833538444482d0,  &
       1.279249234254694d0,       0.461131501011904d0,  &
       1.918666349152554d0,       1.185325727254327d0,  &
       0.869663674465310d0,       0.997099981228532d0, &
       1.151336649298411d0,       0.656477761386754d0,  &
       1.059974763191194d0,       0.873927368065145d0      ], [num_cells,2]); ! Reference for Newton

 

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

end program test_time_propagator_pic_1d2v_vm_subcyc
