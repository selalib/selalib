! TODO: Use input from file to initialize and compare
! Unit test for antisymmetric splitting with coordinate transformation
! author: Benedikt Perse, IPP
program test_time_propagator_pic_3d3v_vm_cef
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
       sll_p_pi, sll_p_twopi

  use sll_m_io_utilities, only : &
       sll_s_concatenate_filename_and_path

  use sll_m_time_propagator_pic_vm_3d3v_cef, only: &
       sll_t_time_propagator_pic_vm_3d3v_cef

  use sll_m_particle_mesh_coupling_spline_3d_feec, only: &
       sll_t_particle_mesh_coupling_spline_3d_feec

  use sll_m_maxwell_3d_base, only: &
       sll_c_maxwell_3d_base

  use sll_m_maxwell_3d_fem_fft, only: &
       sll_t_maxwell_3d_fem_fft

  use sll_m_particle_group_3d3v, only: &
       sll_t_particle_group_3d3v

  use sll_m_particle_group_base, only: &
       sll_t_particle_array


  implicit none
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! Tolerance for comparison of real numbers: set it here!
  sll_real64, parameter :: EQV_TOL = 1.5e-12_f64

  ! Abstract particle group
  class(sll_t_particle_array), pointer :: particle_group
  class(sll_t_particle_group_3d3v), pointer :: pg

  ! Arrays for the fields
  sll_real64, allocatable :: phi(:), efield(:), efield_ref(:)
  sll_real64, allocatable :: bfield(:), bfield_ref(:)
  sll_real64, allocatable :: rho(:), rho_local(:)

  ! Abstract kernel smoothers
  type(sll_t_particle_mesh_coupling_spline_3d_feec), pointer :: particle_mesh_coupling     


  ! Maxwell solver 
  ! Abstract 
  class(sll_c_maxwell_3d_base), pointer :: maxwell_solver

  ! Specific Hamiltonian splitting
  type(sll_t_time_propagator_pic_vm_3d3v_cef) :: propagator

  ! Parameters
  sll_int32  :: n_particles
  sll_real64 :: eta_min(3), eta_max(3), domain(3,3), domain_logical(3,3)
  sll_int32  :: num_cells(3)
  sll_real64 :: delta_t
  sll_int32  :: degree_smoother(3)
  sll_int64  :: rnd_seed
  sll_real64 :: params(6)
  ! Helper 
  sll_int32  :: i_part
  sll_real64 :: xi(3), wi(1)
  logical    :: passed
  sll_real64 :: error
  sll_int32  :: ierr,  n_total   ! error code for SLL_ALLOCATE

  ! Reference
  sll_real64, allocatable :: particle_info_ref(:,:)
  sll_real64, allocatable :: particle_info_check(:,:)

  character(len=256) :: reference_efield
  character(len=256) :: reference_bfield
  sll_int32  :: file_id
  sll_int32   :: io_stat
  character(len=256) :: reffile_full  

  call sll_s_boot_collective()

  ! Set parameters
  n_particles = 1!10
  eta_min = 0.0_f64
  eta_max = sll_p_twopi
  num_cells = [4, 4, 2]
  delta_t = 0.1_f64
  degree_smoother = [2, 2, 1]
  passed = .TRUE.
  rnd_seed = 10


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
  particle_info_ref = reshape([5.89048622548_f64, 2.74889357189_f64, &
       0.9884121822_f64,     -1.53412054435254_f64, &
       0.15731068461017_f64,       1.62412057435254_f64,   5.7026946767517_f64], [7, n_particles])

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
       degree_smoother)

  ! Initialize Maxwell solver
  allocate( sll_t_maxwell_3d_fem_fft :: maxwell_solver )
  select type ( maxwell_solver )
  type is ( sll_t_maxwell_3d_fem_fft )
     call maxwell_solver%init( domain(:,1:2), num_cells, degree_smoother )
  end select

  n_total = maxwell_solver%n_total

  SLL_ALLOCATE(phi(n_total),ierr)
  SLL_ALLOCATE(efield(3*n_total),ierr)
  SLL_ALLOCATE(bfield(3*n_total),ierr)
  SLL_ALLOCATE(efield_ref(6*n_total),ierr)
  SLL_ALLOCATE(bfield_ref(3*n_total),ierr)
  SLL_ALLOCATE(rho(n_total),ierr)
  SLL_ALLOCATE(rho_local(n_total),ierr)


  ! Reference fields 
  !------------------------------------------------
  reference_efield = "ref_cef_3d3v_efield"
  reference_bfield = "ref_cef_3d3v_bfield"

  ! Read reference
  call sll_s_concatenate_filename_and_path( trim(reference_efield), __FILE__,&
       reffile_full)
  open(newunit=file_id, file=reffile_full, status='old', action='read', IOStat=io_stat)
  if (io_stat /= 0) then
     print*, ' failed to open reference file: ', reference_efield
     STOP
  end if
  read(unit=file_id,fmt=*) efield_ref
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
       n_total, MPI_SUM, rho)

  ! Solve Poisson problem
  call maxwell_solver%compute_E_from_rho( rho, efield)
  bfield = 0._f64
  bfield(2*n_total+1:3*n_total) = params(1)

  call propagator%init( maxwell_solver, &
       particle_mesh_coupling, particle_group, &
       phi, efield, bfield, &
       domain(:,1), domain(:,3) )

  call propagator%operatorHp( delta_t )

  ! Compare to reference
  ! Particle information after advect_x application 
  particle_info_check(:,1) = [5.73707417104475_f64, 2.76462464035102_f64, 1.15082423963525_f64, -1.534120544353_f64, 0.15731068461_f64, 1.624120574353_f64, -5.7026946767517_f64] 

  ! Compare computed values to reference values
  do i_part=1,n_particles
     xi = particle_group%group(1)%get_x(i_part)
     if( maxval(abs(xi-particle_info_check(1:3,i_part))) > EQV_TOL ) then
        passed = .FALSE.
        print*, 'error x', maxval(abs(xi-particle_info_check(1:3,i_part))), xi, particle_info_check(1:3,i_part)
     end if
     xi = particle_group%group(1)%get_v(i_part)
     if( maxval(abs(xi-particle_info_check(4:6,i_part))) > EQV_TOL ) then
        passed = .FALSE.
        print*, 'error v', maxval(abs(xi-particle_info_check(4:6,i_part))), xi, particle_info_check(4:6,i_part)
     end if
     xi(1:1) = particle_group%group(1)%get_charge(i_part)
     if( abs(xi(1)-particle_info_check(7,i_part)) > EQV_TOL) then
        passed = .FALSE.
        print*, 'error wi', abs(xi-particle_info_check(7,i_part))
     end if
  end do

  error = maxval(abs(efield-efield_ref(1:3*n_total)))
  if (error> EQV_TOL) then
     passed = .FALSE.
     print*, 'e field wrong in Hp.', error
  end if

  if (passed .EQV. .FALSE.) then
     print*, 'Error in Hp.'
  end if

  !Reset field info
  efield = [  -0.1208722522953525_f64, 6.783545860148128E-002_f64,&
       0.356353462528972_f64,     -0.325896148390076_f64,&
      -0.147238956656419_f64,      8.837289195243920E-002_f64,&
       1.216516029272858E-002_f64, 5.043767141305302E-002_f64,&
      -2.391124245395114E-003_f64,-2.450242988238941E-002_f64,&
       1.317175751084487E-003_f64, 2.207490825650721E-002_f64,&
      -0.118229753583657_f64,      1.867939985271847E-002_f64,&
       5.298130008588725E-002_f64, 5.569224368648823E-002_f64,&
       4.213226136544545E-002_f64,-0.325514456051976_f64,&
       2.011891545553868_f64,     -2.832207983740051_f64,&
      -0.543939198589565_f64,      0.157637303950138_f64,&
       0.117989349078764_f64,      0.464995452522253_f64,&
       5.334407967454392E-002_f64,-0.132532237048027_f64,&
       0.140178681909910_f64,     -0.236905549133555_f64,&
      -0.384522909795975_f64,      0.292520752603064_f64,&
      -0.376662023769924_f64,      0.901020553909419_f64,&
       1.156172332682932E-003_f64, 5.217454881824712E-002_f64,&
       0.356270941860991_f64,      1.089637294651349E-002_f64,&
       2.488823502358686E-002_f64, 4.521205715802832E-002_f64,&
      -0.296803584895001_f64,      4.728723467664095E-002_f64,&
      -2.892072397305311E-003_f64,-0.116226985082640_f64,&
      -0.126004474441369_f64,     -0.147083471922241_f64,&
      -2.492625324470412E-002_f64, 1.840148586491818E-002_f64,&
       6.792128663811199E-002_f64, 9.004259906672866E-002_f64,&
       0.123742774067357_f64,     -0.368990922500690_f64,&
       1.924388039281720_f64,      6.192587750524988E-002_f64,&
      -9.722608061579932E-002_f64, 0.388940679545226_f64,&
      -1.381100831957693_f64,      0.305426436043234_f64,&
       2.794412730977026E-002_f64,-0.284813268813339_f64,&
      -0.224879734474267_f64,     -0.534295210410154_f64,&
      -0.153854282206331_f64,      0.276924648487647_f64,&
      -0.290900199602184_f64,      0.237802761922037_f64,&
       0.377718220708566_f64,     -0.632861468982267_f64,&
       0.732194837605924_f64,     -0.412668937110055_f64,&
      -9.042773638068581E-002_f64, 5.316966977446258E-002_f64,&
      -0.258213915785729_f64,     -9.782890334523987E-002_f64,&
       5.994150110561960E-002_f64,-0.13425559310180_f64,&
       0.305452179004565_f64,     -0.10895712803335_f64,&
      -0.162452265692013_f64,      0.25459363323910_f64,&
      -0.638946888089264_f64,      0.11723867843223_f64,&
       5.282467479262681E-002_f64,-0.12857750721543_f64,&
       0.711118053999076_f64,      9.95989116389381E-002_f64,&
       2.241510727790141E-002_f64, 5.78510833263527E-002_f64,&
       9.959891163893774E-002_f64, 0.20235750901088_f64,&
       2.698172991868413E-002_f64,-2.76401247317766E-002_f64,&
       5.282467479262692E-002_f64, 2.24151072779014E-002_f64,&
      -2.764012473177654E-002_f64, 0.10517288277746_f64,&
      -0.128577507215434_f64,      5.78510833263527E-002_f64 ]

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

  call propagator%operatorHB( delta_t )
  ! Compare to reference
  ! Particle information after advect_ev application 
  particle_info_check(:,1) = [5.89048622548_f64, 2.74889357189_f64, 0.9884121822_f64, -1.53412054435254_f64, 0.15731068461017_f64, 1.62412057435254_f64, -5.7026946767517_f64]
  ! Compare computed values to reference values
  do i_part=1,n_particles
     xi = particle_group%group(1)%get_x(i_part)
     if( maxval(abs(xi-particle_info_check(1:3,i_part))) > EQV_TOL ) then
        passed = .FALSE.
        print*, 'error x', maxval(abs(xi-particle_info_check(1:3,i_part))), xi, particle_info_check(1:3,i_part)
     end if
     xi = particle_group%group(1)%get_v(i_part)
     if( maxval(abs(xi-particle_info_check(4:6,i_part))) > EQV_TOL ) then
        passed = .FALSE.
        print*, 'error v', maxval(abs(xi-particle_info_check(4:6,i_part))), xi, particle_info_check(4:6,i_part)
     end if
     xi(1:1) = particle_group%group(1)%get_charge(i_part)
     if( abs(xi(1)-particle_info_check(7,i_part)) > EQV_TOL) then
        passed = .FALSE.
        print*, 'error wi', abs(xi-particle_info_check(7,i_part))
     end if
  end do

  error = maxval(abs(efield-efield_ref(3*n_total+1:6*n_total)))
  if (error> EQV_TOL) then
     passed = .FALSE.
     print*, 'e field wrong in HB.', error
  end if

  if (passed .EQV. .FALSE.) then
     print*, 'Error in HB.'
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


  call propagator%operatorHE( delta_t )

  ! Compare to reference
  ! Particle information after advect_ev application 
  particle_info_check(:,1) = [5.89048622548_f64, 2.74889357189_f64, 0.9884121822_f64, -1.54612349983218_f64, 0.185197469629497_f64, 1.61648252940304_f64, -5.7026946767517_f64]
  ! Compare computed values to reference values
  do i_part=1,n_particles
     xi = particle_group%group(1)%get_x(i_part)
     if( maxval(abs(xi-particle_info_check(1:3,i_part))) > EQV_TOL ) then
        passed = .FALSE.
        print*, 'error x', maxval(abs(xi-particle_info_check(1:3,i_part))), xi, particle_info_check(1:3,i_part)
     end if
     xi = particle_group%group(1)%get_v(i_part)
     if( maxval(abs(xi-particle_info_check(4:6,i_part))) > EQV_TOL ) then
        passed = .FALSE.
        print*, 'error v', maxval(abs(xi-particle_info_check(4:6,i_part))), xi, particle_info_check(4:6,i_part)
     end if
     xi(1:1) = particle_group%group(1)%get_charge(i_part)
     if( abs(xi(1)-particle_info_check(7,i_part)) > EQV_TOL) then
        passed = .FALSE.
        print*, 'error wi', abs(xi-particle_info_check(7,i_part))
     end if
  end do

  error = maxval(abs(bfield-bfield_ref))
  if (error> EQV_TOL) then
     passed = .FALSE.
     print*, 'b field wrong in HE.', error
  end if

  if (passed .EQV. .FALSE.) then
     print*, 'Error in HE.'
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
  call particle_mesh_coupling%free()
  deallocate(particle_mesh_coupling)
  call maxwell_solver%free()
  deallocate(maxwell_solver)

  call sll_s_halt_collective()

end program test_time_propagator_pic_3d3v_vm_cef
