! TODO: Use input from file to initialize and compare
! Unit test for antisymmetric splitting with coordinate transformation
! author: Benedikt Perse, IPP
program test_time_propagator_pic_3d3v_vm_clamped_disgrad_trafo
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

  use sll_m_io_utilities, only : &
       sll_s_concatenate_filename_and_path

  use sll_m_time_propagator_pic_vm_3d3v_trafo_helper, only: &
       sll_t_time_propagator_pic_vm_3d3v_trafo_helper

  use sll_m_time_propagator_pic_vm_3d3v_cl_helper, only: &
       sll_p_boundary_particles_periodic, &
       sll_p_boundary_particles_singular, &
       sll_p_boundary_particles_reflection, &
       sll_p_boundary_particles_absorption

  use sll_m_particle_mesh_coupling_spline_cl_3d_feec, only: &
       sll_t_particle_mesh_coupling_spline_cl_3d_feec

  use sll_m_maxwell_3d_base, only: &
       sll_c_maxwell_3d_base

  use sll_m_maxwell_clamped_3d_trafo, only: &
       sll_t_maxwell_clamped_3d_trafo

  use sll_m_particle_group_3d3v, only: &
       sll_t_particle_group_3d3v

  use sll_m_particle_group_base, only: &
       sll_t_particle_array

  use sll_m_3d_coordinate_transformations

  use sll_m_mapping_3d, only: &
       sll_t_mapping_3d

  use sll_m_splines_pp, only: &
       sll_p_boundary_periodic, &
       sll_p_boundary_clamped, &
       sll_p_boundary_clamped_square, &
       sll_p_boundary_clamped_cubic

  implicit none
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! Tolerance for comparison of real numbers: set it here!
  sll_real64, parameter :: EQV_TOL = 1.5e-12_f64

  ! Abstract particle group
  class(sll_t_particle_array), pointer :: particle_group
  class(sll_t_particle_group_3d3v), pointer :: pg

  ! Arrays for the fields
  sll_real64, allocatable :: phi(:)
  sll_real64, allocatable :: efield(:), efield_ref1(:), efield_ref2(:), efield_ref3(:)
  sll_real64, allocatable :: bfield(:), bfield_ref(:)
  sll_real64, allocatable :: rho(:), rho_local(:)

  ! Abstract kernel smoothers
  type(sll_t_particle_mesh_coupling_spline_cl_3d_feec), pointer :: particle_mesh_coupling   
  ! Maxwell solver 
  ! Abstract 
  class(sll_c_maxwell_3d_base), pointer :: maxwell_solver

  ! Specific Hamiltonian splitting
  type(sll_t_time_propagator_pic_vm_3d3v_trafo_helper) :: propagator

  type(sll_t_mapping_3d), pointer :: map  

  ! Parameters
  sll_int32  :: n_particles
  sll_real64 :: eta_min(3), eta_max(3), domain(3,3), domain_logical(3,3)
  sll_int32  :: num_cells(3)
  sll_real64 :: delta_t
  sll_int32  :: degree_smoother(3), boundary(3)
  sll_int64  :: rnd_seed
  sll_real64 :: params(6)
  ! Helper 
  sll_int32  :: i_part
  sll_real64 :: xi(3), wi(1)
  logical    :: passed
  sll_real64 :: error
  sll_int32  :: ierr, n_total0, n_total1   ! error code for SLL_ALLOCATE

  ! Reference
  sll_real64, allocatable :: particle_info_ref(:,:)
  sll_real64, allocatable :: particle_info_check(:,:)

  character(len=256) :: reference_efield1
  character(len=256) :: reference_efield2
  character(len=256) :: reference_efield3
  character(len=256) :: reference_bfield
  sll_int32  :: file_id
  sll_int32   :: io_stat
  character(len=256) :: reffile_full  

  call sll_s_boot_collective()

  ! Set parameters
  n_particles = 1!10
  eta_min = 0.0_f64
  eta_max = sll_p_twopi
  num_cells = [8, 4, 2]
  delta_t = 0.1_f64
  degree_smoother = [2, 2, 1]
  passed = .TRUE.
  rnd_seed = 10

  if( degree_smoother(1) == 2 ) then
     boundary = [ sll_p_boundary_clamped_square, sll_p_boundary_periodic, sll_p_boundary_periodic]
  else if( degree_smoother(1) == 3 ) then
     boundary = [ sll_p_boundary_clamped_cubic, sll_p_boundary_periodic, sll_p_boundary_periodic]
  else
     boundary = [ sll_p_boundary_clamped, sll_p_boundary_periodic, sll_p_boundary_periodic]
  end if

  !Initialize mapping
  params = 0._f64
  params(1) = sll_p_twopi
  params(2) = sll_p_twopi
  params(3) = sll_p_twopi
  params(4) = 0.05_f64
  params(5) = 0.05_f64

  allocate(map)
  call map%init( params,&
       sll_f_colbound_x1,&
       sll_f_colbound_x2,&
       sll_f_colbound_x3,&
       sll_f_colbound_jac11,&
       sll_f_colbound_jac12,&
       sll_f_colbound_jac13,&
       sll_f_colbound_jac21,&
       sll_f_colbound_jac22,&
       sll_f_colbound_jac23,&
       sll_f_colbound_jac31,&
       sll_f_colbound_jac32,&
       sll_f_colbound_jac33,&
       sll_f_colbound_jacobian, flag2d = .true., Lx=params(1:3))

  domain(:,1) = eta_min
  domain(:,2) = eta_max
  domain(:,3) = eta_max - eta_min

  domain_logical(:,1) = 0._f64
  domain_logical(:,2) = 1._f64
  domain_logical(:,3) = eta_max - eta_min

  ! Initialize
  allocate(pg)
  call pg%init(n_particles, &
       n_particles, -1.0_f64, 1.0_f64, 1)
  allocate( particle_group)
  particle_group%n_species = 1
  allocate( particle_group%group(1), source=pg )

  call particle_group%group(1)%set_common_weight (1.0_f64)


  SLL_ALLOCATE(particle_info_ref(7, n_particles), ierr)
  SLL_ALLOCATE(particle_info_check(7, n_particles), ierr)
  particle_info_ref = 0.0_f64
  particle_info_ref = reshape([0.9375_f64,        0.4375_f64,       0.15731068461017_f64,     -1.53412054435254_f64,    0.15731068461017_f64,       1.62412057435254_f64,   5.7026946767517_f64], [7, n_particles])

  ! Initialize particles from particle_info_ref
  xi = 0.0_f64
  do i_part =1, n_particles
     xi = particle_info_ref(1:3, i_part) 
     call particle_group%group(1)%set_x(i_part, xi)
     xi = particle_info_ref(4:6, i_part)
     call particle_group%group(1)%set_v(i_part, xi)
     xi(1) = particle_info_ref(7, i_part)
     call particle_group%group(1)%set_weights(i_part, xi(1))
  end do

  call particle_group%group(1)%set_common_weight (1.0_f64)

  ! Initialize kernel smoother
  allocate(particle_mesh_coupling)
  call particle_mesh_coupling%init( num_cells, domain_logical(:,1:2), &
       degree_smoother, boundary )


  ! Initialize Maxwell solver
  allocate( sll_t_maxwell_clamped_3d_trafo :: maxwell_solver )
  select type ( maxwell_solver )
  type is ( sll_t_maxwell_clamped_3d_trafo )
     call maxwell_solver%init( domain(:,1:2), num_cells, degree_smoother, boundary, map, mass_tolerance=1d-12, poisson_tolerance=1d-12, solver_tolerance=1d-12 )
  end select

  n_total0 = maxwell_solver%n_total0
  n_total1 = maxwell_solver%n_total1

  SLL_ALLOCATE(phi(n_total0),ierr)
  SLL_ALLOCATE(efield(n_total1+2*n_total0),ierr)
  SLL_ALLOCATE(bfield(n_total0+2*n_total1),ierr)
  SLL_ALLOCATE(efield_ref1(n_total1+2*n_total0),ierr)
  SLL_ALLOCATE(efield_ref2(n_total1+2*n_total0),ierr)
  SLL_ALLOCATE(efield_ref3(n_total1+2*n_total0),ierr)
  SLL_ALLOCATE(bfield_ref(n_total0+2*n_total1),ierr)
  SLL_ALLOCATE(rho(n_total0),ierr)
  SLL_ALLOCATE(rho_local(n_total0),ierr)

  ! Reference fields 
  !------------------------------------------------
  reference_efield1 = "ref_disgrad_trafo_3d3v_efield1"
  reference_efield2 = "ref_disgrad_trafo_3d3v_efield2"
  reference_efield3 = "ref_disgrad_trafo_3d3v_efield3"
  reference_bfield = "ref_disgrad_trafo_3d3v_bfield"

  ! Read reference
  call sll_s_concatenate_filename_and_path( trim(reference_efield1), __FILE__,&
       reffile_full)
  open(newunit=file_id, file=reffile_full, status='old', action='read', IOStat=io_stat)
  if (io_stat /= 0) then
     print*, ' failed to open reference file: ', reference_efield1
     STOP
  end if
  read(unit=file_id,fmt=*) efield_ref1
  close(file_id)

  call sll_s_concatenate_filename_and_path( trim(reference_efield2), __FILE__,&
       reffile_full)
  open(newunit=file_id, file=reffile_full, status='old', action='read', IOStat=io_stat)
  if (io_stat /= 0) then
     print*, ' failed to open reference file: ', reference_efield2
     STOP
  end if
  read(unit=file_id,fmt=*) efield_ref2
  close(file_id)

  call sll_s_concatenate_filename_and_path( trim(reference_efield3), __FILE__,&
       reffile_full)
  open(newunit=file_id, file=reffile_full, status='old', action='read', IOStat=io_stat)
  if (io_stat /= 0) then
     print*, ' failed to open reference file: ', reference_efield3
     STOP
  end if
  read(unit=file_id,fmt=*) efield_ref3
  close(file_id)

  call sll_s_concatenate_filename_and_path( trim(reference_bfield), __FILE__,&
       reffile_full)
  open(newunit=file_id, file=reffile_full, status='old', action='read', IOStat=io_stat)
  if (io_stat /= 0) then
     print*, ' failed to open reference file: ', reference_bfield
     STOP
  end if
  read(unit=file_id,fmt=*) bfield_ref
  close(file_id)

  rho_local = 0.0_f64
  do i_part = 1,n_particles
     xi = particle_group%group(1)%get_x(i_part)
     wi(1) = particle_group%group(1)%get_charge( i_part)
     call particle_mesh_coupling%add_charge(xi, wi(1), degree_smoother, rho_local)
  end do
  ! MPI to sum up contributions from each processor
  rho = 0.0_f64
  call sll_o_collective_allreduce( sll_v_world_collective, &
       rho_local, &
       n_total0, MPI_SUM, rho)

  ! Solve Poisson problem
  call maxwell_solver%compute_E_from_rho( rho, efield)
  bfield = 0._f64
  bfield(n_total0+n_total1+1:n_total0+2*n_total1) = params(1)

  call propagator%init( maxwell_solver, &
       particle_mesh_coupling, particle_group, &
       phi, efield, bfield, &
       domain(:,1), domain(:,3), map=map, boundary_particles=sll_p_boundary_particles_reflection, iter_tolerance=1d-10, max_iter=5 )

  call propagator%advect_x( delta_t )

  ! Compare to reference
  ! Particle information after advect_x application 
  particle_info_check(:,1) = [ 5.856716218387105_f64, 2.884266687693583_f64, 1.150824239640236_f64, -1.53412054435254_f64, 0.15731068461017_f64, 1.62412057435254_f64, -5.7026946767517_f64] 

  ! Compare computed values to reference values
  do i_part=1,n_particles
     xi = particle_group%group(1)%get_x(i_part)
     xi = map%get_x(xi)
     if( maxval(abs(xi-particle_info_check(1:3,i_part))) > EQV_TOL ) then
        passed = .FALSE.
        print*, 'error x', maxval(abs(xi-particle_info_check(1:3,i_part)))
     end if
     xi = particle_group%group(1)%get_v(i_part)
     if( maxval(abs(xi-particle_info_check(4:6,i_part))) > EQV_TOL ) then
        passed = .FALSE.
        print*, 'error v', maxval(abs(xi-particle_info_check(4:6,i_part)))
     end if
     xi(1:1) = particle_group%group(1)%get_charge(i_part)
     if( abs(xi(1)-particle_info_check(7,i_part)) > EQV_TOL) then
        passed = .FALSE.
        print*, 'error wi', abs(xi-particle_info_check(7,i_part))
     end if
  end do

  if (passed .EQV. .FALSE.) then
     print*, 'Error in advect_x.'
  end if

  ! Reset particle info
  ! Initialize particles from particle_info_ref
  xi = 0.0_f64
  do i_part =1, n_particles
     xi = particle_info_ref(1:3, i_part) 
     call particle_group%group(1)%set_x(i_part, xi)
     xi = particle_info_ref(4:6, i_part)
     call particle_group%group(1)%set_v(i_part, xi)
     xi(1) = particle_info_ref(7, i_part)
     call particle_group%group(1)%set_weights(i_part, xi(1))
  end do


  call propagator%advect_vb( delta_t )

  ! Compare to reference
  ! Particle information after advect_vb application 
  particle_info_check(:,1) = [ 6.010130862223416_f64, 2.868538208633622_f64,       0.9884121822049821_f64, -1.5372451992717_f64, 0.1230840908207219_f64,        1.62412057435254_f64, -5.7026946767517_f64]
  ! Compare computed values to reference values
  do i_part=1,n_particles
     xi = particle_group%group(1)%get_x(i_part)
     xi = map%get_x(xi)
     if( maxval(abs(xi-particle_info_check(1:3,i_part))) > EQV_TOL ) then
        passed = .FALSE.
        print*, 'error x', maxval(abs(xi-particle_info_check(1:3,i_part)))
     end if
     xi = particle_group%group(1)%get_v(i_part)
     if( maxval(abs(xi-particle_info_check(4:6,i_part))) > EQV_TOL ) then
        passed = .FALSE.
        print*, 'error v', maxval(abs(xi-particle_info_check(4:6,i_part)))
     end if
     xi(1:1) = particle_group%group(1)%get_charge(i_part)
     if( abs(xi(1)-particle_info_check(7,i_part)) > EQV_TOL) then
        passed = .FALSE.
        print*, 'error wi', abs(xi-particle_info_check(7,i_part))
     end if
  end do

  if (passed .EQV. .FALSE.) then
     print*, 'Error in advect_vb.'
  end if



  ! Reset particle info
  ! Initialize particles from particle_info_ref
  xi = 0.0_f64
  do i_part =1, n_particles
     xi = particle_info_ref(1:3, i_part) 
     call particle_group%group(1)%set_x(i_part, xi)
     xi = particle_info_ref(4:6, i_part)
     call particle_group%group(1)%set_v(i_part, xi)
     xi(1) = particle_info_ref(7, i_part)
     call particle_group%group(1)%set_weights(i_part, xi(1))
  end do

  call propagator%advect_e( delta_t )

  error = maxval(abs(efield-efield_ref1))
  if (error> EQV_TOL) then
     passed = .FALSE.
     print*, 'e field wrong in advect_e.', error
  end if

  ! Compare to reference
  ! Particle information after advect_ev application 
  particle_info_check(:,1) = [ 6.010130862223416_f64, 2.86853820863362_f64,       0.9884121822049821_f64, -1.503258398841246_f64, 0.1735663654187207_f64, 1.611567924439785_f64, -5.7026946767517_f64]
  ! Compare computed values to reference values
  do i_part=1,n_particles
     xi = particle_group%group(1)%get_x(i_part)
     xi = map%get_x(xi)
     if( maxval(abs(xi-particle_info_check(1:3,i_part))) > EQV_TOL ) then
        passed = .FALSE.
        print*, 'error x', maxval(abs(xi-particle_info_check(1:3,i_part)))
     end if
     xi = particle_group%group(1)%get_v(i_part)
     if( maxval(abs(xi-particle_info_check(4:6,i_part))) > EQV_TOL ) then
        passed = .FALSE.
        print*, 'error v', maxval(abs(xi-particle_info_check(4:6,i_part)))
     end if
     xi(1:1) = particle_group%group(1)%get_charge(i_part)
     if( abs(xi(1)-particle_info_check(7,i_part)) > EQV_TOL) then
        passed = .FALSE.
        print*, 'error wi', abs(xi-particle_info_check(7,i_part))
     end if
  end do

  if (passed .EQV. .FALSE.) then
     print*, 'Error in advect_ev.'
  end if

  call propagator%advect_eb( delta_t*5.0_f64 )

  error = maxval(abs(efield-efield_ref2))
  if (error> EQV_TOL) then
     passed = .FALSE.
     print*, 'e field wrong in advect_eb.', error
  end if
  error = maxval(abs(bfield-bfield_ref))

  if (error> EQV_TOL) then
     passed = .FALSE.
     print*, 'b field wrong in advect_eb.', error
  end if

  ! Reset particle info
  ! Initialize particles from particle_info_ref
  xi = 0.0_f64
  do i_part =1, n_particles
     xi = particle_info_ref(1:3, i_part) 
     call particle_group%group(1)%set_x(i_part, xi)
     xi = particle_info_ref(4:6, i_part)
     call particle_group%group(1)%set_v(i_part, xi)
     xi(1) = particle_info_ref(7, i_part)
     call particle_group%group(1)%set_weights(i_part, xi(1))
  end do

  ! Solve Poisson problem
  call maxwell_solver%compute_E_from_rho( rho, efield)
  bfield = 0._f64
  bfield(n_total0+n_total1+1:n_total0+2*n_total1) = params(1)

  propagator%max_iter=10

  ! Reset particle info
  ! Initialize particles from particle_info_ref
  xi = 0.0_f64
  do i_part =1, n_particles
     xi = particle_info_ref(1:3, i_part) 
     call particle_group%group(1)%set_x(i_part, xi)
     xi = particle_info_ref(4:6, i_part)
     call particle_group%group(1)%set_v(i_part, xi)
     xi(1) = particle_info_ref(7, i_part)
     call particle_group%group(1)%set_weights(i_part, xi(1))
  end do


  call propagator%advect_ex( delta_t )

  error = maxval(abs(efield-efield_ref3))
  if (error> EQV_TOL) then
     passed = .FALSE.
     print*, 'e field wrong in advect_ex.', error
  end if

  ! Compare to reference
  ! Particle information after advect_ev application 
  particle_info_check(:,1) = [ 5.857730198896661_f64, 2.884896293145641_f64,        1.150183285518458_f64, -1.5138394086609_f64, 0.1699043191518476_f64, 1.611301491916905_f64, -5.7026946767517_f64]
  ! Compare computed values to reference values
  do i_part=1,n_particles
     xi = particle_group%group(1)%get_x(i_part)
     xi = map%get_x(xi)
     if( maxval(abs(xi-particle_info_check(1:3,i_part))) > EQV_TOL ) then
        passed = .FALSE.
        print*, 'error x', maxval(abs(xi-particle_info_check(1:3,i_part)))
     end if
     xi = particle_group%group(1)%get_v(i_part)
     if( maxval(abs(xi-particle_info_check(4:6,i_part))) > EQV_TOL ) then
        passed = .FALSE.
        print*, 'error v', maxval(abs(xi-particle_info_check(4:6,i_part)))
     end if
     xi(1:1) = particle_group%group(1)%get_charge(i_part)
     if( abs(xi(1)-particle_info_check(7,i_part)) > EQV_TOL) then
        passed = .FALSE.
        print*, 'error wi', abs(xi-particle_info_check(7,i_part))
     end if
  end do

  if (passed .EQV. .FALSE.) then
     print*, 'Error in advect_ex.'
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
  deallocate(efield_ref1)
  deallocate(efield_ref2)
  deallocate(efield_ref3)
  deallocate(bfield)
  deallocate(bfield_ref)
  deallocate(rho)
  deallocate(rho_local)
  call particle_mesh_coupling%free()
  deallocate(particle_mesh_coupling)
  call maxwell_solver%free()
  deallocate(maxwell_solver)

  call sll_s_halt_collective()


end program test_time_propagator_pic_3d3v_vm_clamped_disgrad_trafo
