! TODO: Use input from file to initialize and compare
! Unit test for antisymmetric splitting with coordinate transformation
! author: Benedikt Perse, IPP
program test_time_propagator_pic_3d3v_vm_trafo
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

  use sll_m_time_propagator_pic_vm_3d3v_trafo_helper, only: &
       sll_t_time_propagator_pic_vm_3d3v_trafo_helper

  use sll_m_particle_mesh_coupling_spline_3d_feec, only: &
       sll_t_particle_mesh_coupling_spline_3d_feec

  use sll_m_maxwell_3d_base, only: &
       sll_c_maxwell_3d_base

  use sll_m_maxwell_3d_trafo, only: &
       sll_t_maxwell_3d_trafo

  use sll_m_particle_group_3d3v, only: &
       sll_t_particle_group_3d3v

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
  class(sll_t_particle_group_3d3v), pointer :: pg

  ! Arrays for the fields
  sll_real64, allocatable :: phi(:)
  sll_real64, allocatable :: efield(:), efield_ref(:)
  sll_real64, allocatable :: bfield(:), bfield_ref(:)
  sll_real64, allocatable :: rho(:), rho_local(:)

  ! Abstract kernel smoothers
  type(sll_t_particle_mesh_coupling_spline_3d_feec), pointer :: particle_mesh_coupling     


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



  !Initialize mapping
  params = 0._f64
  params(1) = sll_p_twopi
  params(2) = sll_p_twopi
  params(3) = sll_p_twopi
  params(4) = 0.1_f64
  params(5) = 0.1_f64

  allocate(map)
  call map%init( params,&
       sll_f_colella_x1,&
       sll_f_colella_x2,&
       sll_f_colella_x3,&
       sll_f_colella_jac11,&
       sll_f_colella_jac12,&
       sll_f_colella_jac13,&
       sll_f_colella_jac21,&
       sll_f_colella_jac22,&
       sll_f_colella_jac23,&
       sll_f_colella_jac31,&
       sll_f_colella_jac32,&
       sll_f_colella_jac33,&
       sll_f_colella_jacobian, flag2d = .true., Lx=params(1:3))

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
  particle_info_ref = reshape([0.937500000000000_f64,        0.437500000000000_f64,       0.15731068461017067_f64,     -1.5341205443525459_f64,    0.15731068461017067_f64,       1.6241205743525459_f64,   5.7026946767517002_f64], [7, n_particles])

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
  allocate( sll_t_maxwell_3d_trafo :: maxwell_solver )
  select type ( maxwell_solver )
  type is ( sll_t_maxwell_3d_trafo )
     !call maxwell_solver%init( domain(:,1:2), num_cells, degree_smoother, map, solver_tolerance=1d-14)
     call maxwell_solver%init( domain(:,1:2), num_cells, degree_smoother, map, mass_tolerance=1d-12, poisson_tolerance=1d-12, solver_tolerance=1d-14 )
  end select

  n_total = maxwell_solver%n_total
  SLL_ALLOCATE(phi(n_total),ierr)
  SLL_ALLOCATE(efield(3*n_total),ierr)
  SLL_ALLOCATE(bfield(3*n_total),ierr)
  SLL_ALLOCATE(efield_ref(3*n_total),ierr)
  SLL_ALLOCATE(bfield_ref(3*n_total),ierr)
  SLL_ALLOCATE(rho(n_total),ierr)
  SLL_ALLOCATE(rho_local(n_total),ierr)

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
       domain(:,1), domain(:,3), map, solver_tolerance=1d-14, iter_tolerance=1d-10, max_iter=5)

  call propagator%advect_x( delta_t )

  ! Compare to reference
  ! Particle information after advect_x application 
  particle_info_check(:,1) = [5.6450994972205857_f64, 2.6726499665270644_f64, 1.1508242396402408_f64, -1.5341205443525459_f64, 0.15731068461017067_f64, 1.6241205743525460_f64, -5.7026946767517002_f64] 

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
  particle_info_check(:,1) = [5.7984711070298012_f64, 2.6568784534400081_f64, 0.98841218220498617_f64, -1.5357608514931156_f64, 0.14039623612551846_f64, 1.6241205743525460_f64, -5.7026946767517002_f64]
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

  efield_ref = [  -0.14206646549074659_f64,7.0888840061299266E-002_f64,&
       0.35636291469320819_f64,-0.22032526058074112_f64,&
       -0.14405176058429067_f64,8.8449348361050437E-002_f64,&
       1.0610105661923164E-002_f64,3.2973539414827738E-002_f64,&
       -6.2779320039068002E-003_f64,-2.4491788068771760E-002_f64,&
       2.7255166091282369E-003_f64,3.8411280173855762E-002_f64,&
       -0.11033442853248576_f64,1.7787687952167199E-002_f64,&
       5.0900674492875254E-002_f64,1.6861643037516815E-002_f64,&
       8.7443577651757073E-002_f64,-0.35188369662460939_f64,&
       2.0825727656940440_f64,-3.0384036964202714_f64,&
       -0.55710436549198417_f64,0.16168206560689247_f64,&
       0.11014680756036981_f64,0.51139555270288706_f64,&
       6.3811973532424418E-002_f64,-0.13736643591309369_f64,&
       0.15131791913161949_f64,-0.27280875408531885_f64,&
       -0.40239138558748622_f64,0.30244378278183259_f64,&
       -0.41033232490654159_f64,0.97657194986231144_f64,&
       4.1123675937858638E-003_f64,4.8419205876077495E-002_f64,&
       0.36521454684387161_f64,1.3674456579715239E-002_f64,&
       2.9911389920207081E-002_f64,4.7657323790881806E-002_f64,&
       -0.30894759718005510_f64,4.2839864433428840E-002_f64,&
       -4.7115038746653416E-003_f64,-0.11640897792221376_f64,&
       -0.12528868718641192_f64,-0.14472197409144383_f64,&
       -2.3192223558355585E-002_f64,1.8818599800606624E-002_f64,&
       6.7501505422252722E-002_f64,8.3707030478428618E-002_f64,&
       0.12522582294823606_f64,-0.36364618592010223_f64,&
       1.9160389130976065_f64,5.2494370658436845E-002_f64,&
       -0.11289254355198773_f64,0.39718463779747543_f64,&
       -1.3710706115508482_f64,0.32577135500158338_f64,&
       3.4341328262850808E-002_f64,-0.28810533128302956_f64,&
       -0.22821211289806775_f64,-0.54449502570970343_f64,&
       -0.16181633181090088_f64,0.28304829516121460_f64,&
       -0.28815467152736046_f64,0.25090362992211834_f64,&
       0.31073836684068951_f64,-0.61138178038899493_f64,&
       0.67332369492527988_f64,-0.29031936205559944_f64,&
       -9.1808371258457533E-002_f64,6.7819989039837075E-002_f64,&
       -0.29031936205559988_f64,-9.0356355342549086E-002_f64,&
       5.1873358069521966E-002_f64,-0.13712389664449146_f64,&
       0.31073836684068973_f64,-9.1808371258457339E-002_f64,&
       -0.13712389664449132_f64,0.24661611348109788_f64,&
       -0.61138178038899460_f64,6.7819989039836465E-002_f64,&
       5.2824674792626813E-002_f64,-0.12857750721543390_f64,&
       0.71111805399907602_f64,9.9598911638938192E-002_f64,&
       2.2415107277901414E-002_f64,5.7851083326352748E-002_f64,&
       9.9598911638937748E-002_f64,0.20235750901088362_f64,&
       2.6981729918684139E-002_f64,-2.7640124731776600E-002_f64,&
       5.2824674792626924E-002_f64,2.2415107277901414E-002_f64,&
       -2.7640124731776544E-002_f64,0.10517288277746234_f64,&
       -0.12857750721543401_f64,5.7851083326352748E-002_f64 ]

  error = maxval(abs(efield-efield_ref))
  if (error> EQV_TOL) then
     passed = .FALSE.
     print*, 'e field wrong in advect_e.', error
  end if

  ! Compare to reference
  ! Particle information after advect_ev application 
  particle_info_check(:,1) = [5.7984711070298012_f64, 2.6568784534400081_f64, 0.98841218220498617_f64, -1.5278225724228562_f64, 0.15948663055539666_f64, 1.6258095163014388_f64, -5.7026946767517002_f64]
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

  bfield_ref =  [ -0.53271095926650691_f64,1.3308679942987209_f64,&
       -2.9583688867871807_f64,0.29564914839532797_f64,&
       0.49154531342942409_f64,-1.0749935158387103_f64,&
       2.1252266889261620_f64,-0.41532257390093352_f64,&
       -0.16171172064557282_f64,0.32431510100307570_f64,&
       -0.51142373149335374_f64,0.27688518913160776_f64,&
       0.29433835629689986_f64,-0.59798674285987219_f64,&
       1.3283575128269782_f64,-0.25306180692833363_f64,&
       -4.4326446757240123E-002_f64,6.6946107915122732E-002_f64,&
       -0.17804732034477744_f64,-3.3183057558774154E-002_f64,&
       -5.5998389008615243E-002_f64,-3.5731745282422678E-002_f64,&
       0.16293034792121985_f64,7.6912622634296823E-002_f64,&
       1.6784029989544428E-002_f64,5.5767976652454365E-003_f64,&
       1.8121333013699226E-003_f64,-1.2539911636143930E-002_f64,&
       -7.9201840379333199E-003_f64,-1.8993996901159826E-002_f64,&
       2.9513255649581460E-002_f64,6.4660389862952694E-002_f64,&
       0.70937330909728413_f64,-1.3491332475410973_f64,&
       3.0481347467401392_f64,-3.4816169664164502_f64,&
       -0.32619694810917482_f64,0.20040386151880929_f64,&
       -0.26025660711514770_f64,0.50759402354405214_f64,&
       0.18994465527575921_f64,-0.28370562536991395_f64,&
       0.56557419109598728_f64,-0.59100115872734293_f64,&
       -0.39847438186353340_f64,0.60648180515662653_f64,&
       -1.3399708301969131_f64,1.4829450802721924_f64,&
       -0.19976171090530287_f64,-8.7375209354619288E-003_f64,&
       8.8281460391818944E-002_f64,1.1934599295690700_f64,&
       6.1851670613710608E-002_f64,-1.2002444453426797E-002_f64,&
       -2.2094837214089030E-003_f64,-0.16918407227741389_f64,&
       -4.5016964619731245E-002_f64,-2.7125468890530913E-003_f64,&
       1.1463214927761006E-002_f64,0.15545423430653363_f64,&
       6.8582483195213900E-002_f64,1.0498934604405041E-002_f64,&
       -5.7843272016930347E-002_f64,-0.37221981915106106_f64,&
       8.8407916228187791_f64,2.5246958314784589_f64,&
       8.8479061011716844_f64,6.4987105760329502_f64,&
       2.5261865123329779_f64,7.3369643479190394_f64,&
       6.4288400555854066_f64,7.2694867012696198_f64,&
       8.8559597452211030_f64,6.4196175340111052_f64,&
       8.8573546368544118_f64,2.5363390401315007_f64,&
       6.4223729365246589_f64,7.3373726884080419_f64,&
       2.5180014613754818_f64,7.3103651237381602_f64,&
       9.3997667894167432_f64,2.3406236332134669_f64,&
       9.0456434928739853_f64,3.0776225530206811_f64,&
       2.2687817725807751_f64,7.3029609075025839_f64,&
       6.2124977574045772_f64,10.166810155536433_f64,&
       9.0110491324366109_f64,6.4234522035309558_f64,&
       8.8375385682838239_f64,1.9157665168386337_f64,&
       6.3403351886965407_f64,7.3330973511475781_f64,&
       2.5595999301633814_f64,8.2954189622266075_f64]

  efield_ref = [1.0957167250965685_f64,-0.66081218254844432_f64,&
       0.55454806741937990_f64,-0.25227843829315766_f64,&
       -1.3596543394066656_f64,-5.7593575597489122E-002_f64,&
       -8.4151863499605711E-002_f64,0.74418913439631540_f64,&
       8.7455525180043014E-002_f64,0.11865489815668974_f64,&
       1.2250243030980754_f64,-0.67420148689927140_f64,&
       -0.21362256441771568_f64,0.72121350798336648_f64,&
       -1.2002049830939834_f64,-5.9930239847901169E-002_f64,&
       1.2068943626466464_f64,-0.97755285785791801_f64,&
       2.0268727862268814_f64,-2.3803933571994134_f64,&
       -1.7828583010502519_f64,2.2067156159435646E-002_f64,&
       3.6027116369819395E-002_f64,0.97250321656597516_f64,&
       0.13827914193293389_f64,1.9243880882726227E-002_f64,&
       1.3523177213094320_f64,-0.78587100599590165_f64,&
       -0.44509752030479344_f64,0.95365573697438399_f64,&
       -1.5506226186078937_f64,0.50040702409877191_f64,&
       -1.1812555694299687_f64,1.2040025295661885_f64,&
       0.38270168531999310_f64,8.9181404667271036E-002_f64,&
       0.69043264204103105_f64,0.23968678268160920_f64,&
       -0.53381578470354529_f64,-0.62213826772202108_f64,&
       -8.2450523438509102E-002_f64,-4.3877938898123206E-002_f64,&
       -1.3078759817638650_f64,1.0454506916526936_f64,&
       -0.18833074005232917_f64,-0.63969085879213083_f64,&
       0.69975595308201588_f64,0.24909222138229725_f64,&
       -1.0196804390170950_f64,0.82621716957861857_f64,&
       1.8544743803018717_f64,2.7805649319178768E-002_f64,&
       0.56293177374629433_f64,0.54626476832268855_f64,&
       -1.5780640605447160_f64,-0.25583505256544858_f64,&
       -6.5707354597888248E-002_f64,-0.18539936820007327_f64,&
       -1.4421841311784103_f64,0.69893904262409523_f64,&
       -0.30080839640680163_f64,-0.42136299256032500_f64,&
       0.40883476718896788_f64,0.39320154473477187_f64,&
       -5.8350885720122377E-002_f64,-0.11560738003384406_f64,&
       4.8554606487851391E-002_f64,3.6773744251731937E-003_f64,&
       -1.2197060166933662E-002_f64,-4.7868283447802390E-002_f64,&
       -5.9843380351448976E-002_f64,-0.10558041124930717_f64,&
       -2.2065069338408846E-002_f64,-1.0630607435166637E-002_f64,&
       -3.0861590191284344E-002_f64,-7.8712341415777161E-003_f64,&
       -1.8552886177865603E-003_f64,-3.7662655497090747E-002_f64,&
       -0.11381027436724765_f64,-7.0090748909620054E-002_f64,&
       7.1263049029022599E-002_f64,-0.12582889589606971_f64,&
       0.72202428818540820_f64,5.3906682741723284E-002_f64,&
       2.2518604668897491E-002_f64,5.1526908668402767E-002_f64,&
       0.11742743751834474_f64,0.19885281559784454_f64,&
       2.9065952527638845E-002_f64,-2.5167572270096243E-002_f64,&
       4.3773753500185079E-002_f64,1.8192778643541768E-002_f64,&
       -3.3942891087476022E-002_f64,0.10229219791803155_f64,&
       -0.12762600156266607_f64,7.3749374130219897E-002_f64 ]

  error = maxval(abs(efield-efield_ref))
  if (error> EQV_TOL) then
     passed = .FALSE.
     print*, 'e field wrong in advect_eb_avf.', error
  end if
  error = maxval(abs(bfield-bfield_ref))

  if (error> EQV_TOL) then
     passed = .FALSE.
     print*, 'b field wrong in advect_eb_avf.', error
  end if

  !call propagator%free()
  !particle_group => null()
  !call pg%free()
  !deallocate(pg)
  !deallocate(efield)
  !deallocate(efield_ref)
  !deallocate(bfield)
  !deallocate(bfield_ref)
  !deallocate(rho)
  !deallocate(rho_local)
  !call particle_mesh_coupling%free()
  !deallocate(particle_mesh_coupling)
  !call maxwell_solver%free()
  !deallocate(maxwell_solver)


!!$  ! Set parameters
!!$  n_particles = 1!10
!!$  eta_min = 0.0_f64
!!$  eta_max = sll_p_twopi
!!$  num_cells = [4, 4, 2]
!!$  delta_t = 0.1_f64
!!$  degree_smoother = [2, 2, 1]
!!$  passed = .TRUE.
!!$  rnd_seed = 10
!!$
!!$
!!$
!!$  !Initialize mapping
!!$  params = 0._f64
!!$  params(1) = sll_p_twopi
!!$  params(2) = sll_p_twopi
!!$  params(3) = sll_p_twopi
!!$  params(4) = 0.1_f64
!!$  params(5) = 0.1_f64
!!$
!!$  allocate(map)
!!$  call map%init( params,&
!!$       sll_f_colella_x1,&
!!$       sll_f_colella_x2,&
!!$       sll_f_colella_x3,&
!!$       sll_f_colella_jac11,&
!!$       sll_f_colella_jac12,&
!!$       sll_f_colella_jac13,&
!!$       sll_f_colella_jac21,&
!!$       sll_f_colella_jac22,&
!!$       sll_f_colella_jac23,&
!!$       sll_f_colella_jac31,&
!!$       sll_f_colella_jac32,&
!!$       sll_f_colella_jac33,&
!!$       sll_f_colella_jacobian)
!!$
!!$  domain(:,1) = eta_min
!!$  domain(:,2) = eta_max
!!$  domain(:,3) = eta_max - eta_min
!!$
!!$  domain_logical(:,1) = 0._f64
!!$  domain_logical(:,2) = 1._f64
!!$  domain_logical(:,3) = eta_max - eta_min
!!$  ! Initialize
!!$  allocate(pg)
!!$  call pg%init(n_particles, &
!!$       n_particles, -1.0_f64, 1.0_f64, 1)
!!$  allocate( particle_group)
!!$  particle_group%n_species = 1
!!$  allocate( particle_group%group(1), source=pg )
!!$
!!$  call particle_group%group(1)%set_common_weight (1.0_f64)
!!$
!!$
!!$  SLL_ALLOCATE(particle_info_ref(7, n_particles), ierr)
!!$  SLL_ALLOCATE(particle_info_check(7, n_particles), ierr)
!!$  particle_info_ref = 0.0_f64
!!$  particle_info_ref = reshape([0.937500000000000_f64,        0.437500000000000_f64,       0.15731068461017067_f64,     -1.5341205443525459_f64,    0.15731068461017067_f64,       1.6241205743525459_f64,   5.7026946767517002_f64], [7, n_particles])
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

!!$  call particle_group%group(1)%set_common_weight (1.0_f64)
!!$
!!$  ! Initialize kernel smoother
!!$  allocate(particle_mesh_coupling)
!!$  call particle_mesh_coupling%init( num_cells, domain_logical(:,1:2), &
!!$       degree_smoother)
!!$
!!$  ! Initialize Maxwell solver
!!$  allocate( sll_t_maxwell_3d_trafo :: maxwell_solver )
!!$  select type ( maxwell_solver )
!!$  type is ( sll_t_maxwell_3d_trafo )
!!$     call maxwell_solver%init( domain(:,1:2), num_cells, degree_smoother, map)
!!$  end select
!!$
!!$  n_total = maxwell_solver%n_total
!!$
!!$  SLL_ALLOCATE(efield(3*n_total),ierr)
!!$  SLL_ALLOCATE(bfield(3*n_total),ierr)
!!$  SLL_ALLOCATE(efield_ref(3*n_total),ierr)
!!$  SLL_ALLOCATE(bfield_ref(3*n_total),ierr)
!!$  SLL_ALLOCATE(rho(n_total),ierr)
!!$  SLL_ALLOCATE(rho_local(n_total),ierr)
!!$
!!$  rho_local = 0.0_f64
!!$  do i_part = 1,n_particles
!!$     xi = particle_group%group(1)%get_x(i_part)
!!$     wi(1) = particle_group%group(1)%get_charge( i_part)
!!$     call particle_mesh_coupling%add_charge(xi, wi(1), degree_smoother, rho_local)
!!$  end do
!!$  ! MPI to sum up contributions from each processor
!!$  rho = 0.0_f64
!!$  call sll_o_collective_allreduce( sll_v_world_collective, &
!!$       rho_local, &
!!$       n_total, MPI_SUM, rho)

  ! Solve Poisson problem
  call maxwell_solver%compute_E_from_rho( rho, efield)
  bfield = 0._f64
  bfield(2*n_total+1:3*n_total) = params(1)
!!$
!!$  call propagator%init( maxwell_solver, &
!!$       particle_mesh_coupling, particle_group, &
!!$       efield, bfield, &
!!$       domain(:,1), domain(:,3), map, solver_tolerance=1d-14, iter_tolerance=1d-12, max_iter=10)

  propagator%iter_tolerance=1d-12
  propagator%max_iter=10


  call propagator%advect_vb( delta_t )

  ! Compare to reference
  ! Particle information after advect_vb application 
  particle_info_check(:,1) = [5.7984711070298012_f64, 2.6568784534400081_f64, 0.98841218220498617_f64, -1.5357608514931156_f64, 0.14039623612551846_f64, 1.6241205743525460_f64, -5.7026946767517002_f64]
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

  call propagator%advect_eb( delta_t*5.0_f64 )

  bfield_ref =  [ 0._f64,        0._f64,&
       0._f64,        0._f64,&
       0._f64,        0._f64,&
       0._f64,        0._f64,&
       0._f64,        0._f64,&
       0._f64,        0._f64,&
       0._f64,        0._f64,&
       0._f64,        0._f64,&
       0._f64,        0._f64,&
       0._f64,        0._f64,&
       0._f64,        0._f64,&
       0._f64,        0._f64,&
       0._f64,        0._f64,&
       0._f64,        0._f64,&
       0._f64,        0._f64,&
       0._f64,        0._f64,&
       0._f64,        0._f64,&
       0._f64,        0._f64,&
       0._f64,        0._f64,&
       0._f64,        0._f64,&
       0._f64,        0._f64,&
       0._f64,        0._f64,&
       0._f64,        0._f64,&
       0._f64,        0._f64,&
       0._f64,        0._f64,&
       0._f64,        0._f64,&
       0._f64,        0._f64,&
       0._f64,        0._f64,&
       0._f64,        0._f64,&
       0._f64,        0._f64,&
       0._f64,        0._f64,&
       0._f64,        0._f64,&
       8.8561684283915234_f64,       2.5187234515222334_f64,&
       8.8561684283915234_f64,       6.4207758944380418_f64,&
       2.5187234515222343_f64,       7.3370734543665499_f64,&
       6.4207758944380364_f64,       7.3370734543665455_f64,&
       8.8561684283915234_f64,       6.4207758944380382_f64,&
       8.8561684283915252_f64,       2.5187234515222361_f64,&
       6.4207758944380409_f64,       7.3370734543665463_f64,&
       2.5187234515222343_f64,       7.3370734543665472_f64,&
       8.8561684283915199_f64,       2.5187234515222379_f64,&
       8.8561684283915199_f64,       6.4207758944380426_f64,&
       2.5187234515222388_f64,       7.3370734543665446_f64,&
       6.4207758944380435_f64,       7.3370734543665446_f64,&
       8.8561684283915181_f64,       6.4207758944380444_f64,&
       8.8561684283915199_f64,       2.5187234515222370_f64,&
       6.4207758944380453_f64,       7.3370734543665428_f64,&
       2.5187234515222374_f64,       7.3370734543665463_f64]

  efield_ref = [1.0759309337915530_f64,-0.61267588550112562_f64,&
       0.44352462608227139_f64,-0.14723218736020990_f64,&
       -1.3484278882032075_f64,-6.5025588153098113E-002_f64,&
       -7.3379099199841868E-002_f64,0.72728508854365770_f64,&
       8.1597314891300385E-002_f64,0.12891904332475318_f64,&
       1.2037956343442624_f64,-0.65476450554782695_f64,&
       -0.20124379975253059_f64,0.69928628506925106_f64,&
       -1.1504757641293459_f64,-0.10711420819986267_f64,&
       0.98238246009893060_f64,-0.97548024951724632_f64,&
       2.1229157485112897_f64,-1.3702704720804868_f64,&
       -1.7083126916691702_f64,5.8463639438028903E-003_f64,&
       1.0116557425329797E-002_f64,0.93280228328755066_f64,&
       9.0730560172863697E-002_f64,1.9675334023833507E-002_f64,&
       1.3647252333930679_f64,-0.71558364057727841_f64,&
       -0.37222621586878779_f64,0.96491230008772821_f64,&
       -1.6179765441151372_f64,0.26574297288371074_f64,&
       -1.1995899074360632_f64,1.2529097776509128_f64,&
       0.27392704665059181_f64,9.6218480231789333E-002_f64,&
       0.70631180832916407_f64,0.20007395849022780_f64,&
       -0.45442035405033576_f64,-0.63379122533342747_f64,&
       -8.8000264540291695E-002_f64,-3.1646220320892315E-002_f64,&
       -1.3274546079887273_f64,1.0549576535771381_f64,&
       -0.17826912336529799_f64,-0.66179002880775961_f64,&
       0.74840042837598408_f64,0.24216257853698767_f64,&
       -1.0386603083872552_f64,0.78540899766511885_f64,&
       1.9533181690796133_f64,0.17971413685695892_f64,&
       0.64549267329971161_f64,0.57293113957380293_f64,&
       -1.6774586387706150_f64,-0.42827403058953528_f64,&
       -7.8867019258724927E-002_f64,-0.20262863643715223_f64,&
       -1.4210030816813475_f64,0.69507285011117292_f64,&
       -0.28751283266621980_f64,-0.39616401378928207_f64,&
       0.38559606435986188_f64,0.31303453063389086_f64,&
       -5.2824674792626848E-002_f64,0.12857750721543390_f64,&
       -0.71111805399907602_f64,-9.9598911638938151E-002_f64,&
       -2.2415107277901410E-002_f64,-5.7851083326352741E-002_f64,&
       -9.9598911638937707E-002_f64,-0.20235750901088365_f64,&
       -2.6981729918684125E-002_f64,2.7640124731776562E-002_f64,&
       -5.2824674792626966E-002_f64,-2.2415107277901435E-002_f64,&
       2.7640124731776520E-002_f64,-0.10517288277746233_f64,&
       0.12857750721543404_f64,-5.7851083326352741E-002_f64,&
       5.2824674792626848E-002_f64,-0.12857750721543390_f64,&
       0.71111805399907602_f64,9.9598911638938151E-002_f64,&
       2.2415107277901410E-002_f64,5.7851083326352741E-002_f64,&
       9.9598911638937707E-002_f64,0.20235750901088365_f64,&
       2.6981729918684125E-002_f64,-2.7640124731776562E-002_f64,&
       5.2824674792626966E-002_f64,2.2415107277901435E-002_f64,&
       -2.7640124731776520E-002_f64,0.10517288277746233_f64,&
       -0.12857750721543404_f64,5.7851083326352741E-002_f64 ]

  error = maxval(abs(efield-efield_ref))
  if (error> EQV_TOL) then
     passed = .FALSE.
     print*, 'e field wrong in advect_eb_disg.', error
  end if
  error = maxval(abs(bfield-bfield_ref))

  if (error> EQV_TOL) then
     passed = .FALSE.
     print*, 'b field wrong in advect_eb_disg.', error
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


  call propagator%advect_ex( delta_t )

  efield_ref = [ 1.0808255965780844_f64, -0.61270404810453938_f64, &
       0.44115793268610226_f64,      -0.17232711949883883_f64, &
       -1.3489279631422635_f64,       -6.5225865947767178E-002_f64, &
       -7.2631867420990232E-002_f64,  0.73097090368301609_f64, &
       8.2406160263130648E-002_f64,  0.12909308341009085_f64, &
       1.2030118700840577_f64,      -0.65846092765098707_f64, &
       -0.20302991929767436_f64,       0.69921759730699617_f64, &
       -1.1487134919191349_f64,       -9.7894138381176513E-002_f64, &
       1.2436081394008807_f64,       -1.0060438494400965_f64, &
       2.0969536892400344_f64,       -2.6776441253735745_f64, &
       -1.7454650406646872_f64,        3.8477737658581300E-003_f64, &
       3.3251684320049740E-002_f64,   1.1450815456757559_f64, &
       0.13802069000609979_f64,        2.1122004732005362E-002_f64, &
       1.3419733507516234_f64,      -0.91705671419885193_f64, &
       -0.46922577814273875_f64,       0.97304800314807838_f64, &
       -1.5784923429856546_f64,       0.74700059596836055_f64, &
       -1.2005371635414954_f64,        1.2538694527297909_f64, &
       0.27146511591384020_f64,        9.5694627233087873E-002_f64, &
       0.70543875237470643_f64,       0.19877804464625928_f64, &
       -0.45031712040823557_f64,      -0.63326806284216475_f64, &
       -8.7694198973684995E-002_f64,  -3.1422486238521584E-002_f64, &
       -1.3277125955996074_f64,        1.0546139486198074_f64, &
       -0.17852024492006918_f64,      -0.66213801874249079_f64, &
       0.74846286548037910_f64,       0.24363710894384794_f64, &
       -1.0779021489564116_f64,       0.83268199689844535_f64, &
       1.8395374461198795_f64,       0.14674072606091337_f64, &
       0.58390121447047505_f64,       0.54139279401847429_f64, &
       -1.5315965174998536_f64,      -0.37581042558633565_f64, &
       -5.7014204647059280E-002_f64, -0.19977538414156037_f64, &
       -1.4272084530644085_f64,       0.66758354169402989_f64, &
       -0.30741602150451514_f64,      -0.40365324108381018_f64, &
       0.38976448007314285_f64,       0.39130356972906322_f64, &
       0.37785917621188492_f64,      -0.63312505973792632_f64, &
       0.73266455128375430_f64,      -0.41275664215272356_f64, &
       -9.0399941417334692E-002_f64,   5.3105836307543487E-002_f64, &
       -0.25799984505905321_f64,       -9.7809812618065733E-002_f64, &
       5.9945796005926012E-002_f64, -0.13426517842113742_f64, &
       0.30544412899083628_f64,      -0.10896379064552378_f64, &
       -0.16251624157559658_f64,       0.25472258255910196_f64, &
       -0.63920633099139657_f64,       0.11729589573307428_f64, &
       5.2824674792627001E-002_f64, -0.12857750721543440_f64, &
       0.71111805399907657_f64,        9.9598911638937540E-002_f64, &
       2.2415107277900980E-002_f64,   5.7851083326353032E-002_f64, &
       9.9598911638937068E-002_f64,  0.20235750901088417_f64, &
       2.6981729918684326E-002_f64,  -2.7640124731776790E-002_f64, &
       5.2824674792626897E-002_f64,   2.2415107277900966E-002_f64, &
       -2.7640124731776912E-002_f64,  0.10517288277746292_f64, &
       -0.12857750721543451_f64,        5.7851083326353407E-002_f64 ]

  error = maxval(abs(efield-efield_ref))
  if (error> EQV_TOL) then
     passed = .FALSE.
     print*, 'e field wrong in advect_ex.', error
  end if

  ! Compare to reference
  ! Particle information after advect_ev application 
  particle_info_check(:,1) = [5.6451678535987559_f64, 2.6729422370841269_f64, 1.1509024733343229_f64, -1.5327598632956168_f64, 0.16314964924494874_f64, 1.6256852482341670_f64, -5.7026946767517002_f64]
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


end program test_time_propagator_pic_3d3v_vm_trafo
