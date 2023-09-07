! TODO: Use input from file to initialize and compare
! Unit test for antisymmetric splitting with coordinate transformation
! author: Benedikt Perse, IPP
program test_time_propagator_pic_3d3v_vm_disgrad_trafo
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
  particle_info_check(:,1) = [5.64509949722058_f64, 2.67264996652706_f64, 1.15082423964024_f64, -1.53412054435254_f64, 0.15731068461017_f64, 1.62412057435254_f64, -5.7026946767517_f64] 

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
  particle_info_check(:,1) = [5.7984711070298_f64, 2.65687845344_f64, 0.988412182204986_f64, -1.53576085149311_f64, 0.140396236125518_f64, 1.62412057435254_f64, -5.7026946767517_f64]
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

  efield_ref = [  -0.14206646549074_f64,7.08888400612992E-002_f64,&
       0.35636291469320_f64,     -0.22032526058074_f64,&
      -0.14405176058429_f64,      8.84493483610504E-002_f64,&
       1.06101056619231E-002_f64, 3.29735394148277E-002_f64,&
      -6.27793200390680E-003_f64,-2.44917880687717E-002_f64,&
       2.72551660912823E-003_f64, 3.84112801738557E-002_f64,&
      -0.11033442853248_f64,      1.77876879521671E-002_f64,&
       5.09006744928752E-002_f64, 1.68616430375168E-002_f64,&
       8.74435776517570E-002_f64,-0.35188369662460_f64,&
       2.08257276569404_f64,     -3.03840369642027_f64,&
      -0.55710436549198_f64,      0.16168206560689_f64,&
       0.11014680756036_f64,      0.51139555270288_f64,&
       6.38119735324244E-002_f64,-0.13736643591309_f64,&
       0.15131791913161_f64,     -0.27280875408531_f64,&
      -0.40239138558748_f64,      0.30244378278183_f64,&
      -0.41033232490654_f64,      0.97657194986231_f64,&
       4.11236759378586E-003_f64, 4.84192058760774E-002_f64,&
       0.36521454684387_f64,      1.36744565797152E-002_f64,&
       2.99113899202070E-002_f64, 4.76573237908818E-002_f64,&
      -0.30894759718005_f64,      4.28398644334288E-002_f64,&
      -4.71150387466534E-003_f64,-0.11640897792221_f64,&
      -0.12528868718641_f64,     -0.14472197409144_f64,&
      -2.31922235583555E-002_f64, 1.88185998006066E-002_f64,&
       6.75015054222527E-002_f64, 8.37070304784286E-002_f64,&
       0.12522582294823_f64,     -0.36364618592010_f64,&
       1.91603891309760_f64,      5.24943706584368E-002_f64,&
      -0.11289254355198_f64,      0.39718463779747_f64,&
      -1.37107061155084_f64,      0.32577135500158_f64,&
       3.43413282628508E-002_f64,-0.28810533128302_f64,&
      -0.22821211289806_f64,     -0.54449502570970_f64,&
      -0.16181633181090_f64,      0.28304829516121_f64,&
      -0.28815467152736_f64,      0.25090362992211_f64,&
       0.31073836684068_f64,     -0.61138178038899_f64,&
       0.67332369492527_f64,     -0.29031936205559_f64,&
      -9.18083712584575E-002_f64, 6.78199890398370E-002_f64,&
      -0.29031936205559_f64,     -9.03563553425490E-002_f64,&
       5.18733580695219E-002_f64,-0.13712389664449_f64,&
       0.31073836684068_f64,     -9.18083712584573E-002_f64,&
      -0.13712389664449_f64,      0.24661611348109_f64,&
      -0.61138178038899_f64,      6.78199890398364E-002_f64,&
       5.28246747926268E-002_f64,-0.12857750721543_f64,&
       0.71111805399907_f64,      9.95989116389381E-002_f64,&
       2.24151072779014E-002_f64, 5.78510833263527E-002_f64,&
       9.95989116389377E-002_f64, 0.20235750901088_f64,&
       2.69817299186841E-002_f64,-2.76401247317766E-002_f64,&
       5.28246747926269E-002_f64, 2.24151072779014E-002_f64,&
      -2.76401247317765E-002_f64, 0.10517288277746_f64,&
      -0.12857750721543_f64,      5.78510833263527E-002_f64 ]

  error = maxval(abs(efield-efield_ref))
  if (error> EQV_TOL) then
     passed = .FALSE.
     print*, 'e field wrong in advect_e.', error
  end if

  ! Compare to reference
  ! Particle information after advect_ev application 
  particle_info_check(:,1) = [5.7984711070298_f64, 2.65687845344_f64, 0.988412182204986_f64, -1.52782257242285_f64, 0.159486630555396_f64, 1.62580951630143_f64, -5.7026946767517_f64]
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

  bfield_ref =  [ -0.5327109592665069_f64,1.3308679942987209_f64,&
      -2.958368886787180_f64, 0.295649148395327_f64,&
       0.491545313429424_f64,-1.074993515838710_f64,&
       2.125226688926162_f64,-0.415322573900933_f64,&
      -0.161711720645572_f64, 0.324315101003075_f64,&
      -0.511423731493353_f64, 0.276885189131607_f64,&
       0.294338356296899_f64,-0.597986742859872_f64,&
       1.328357512826978_f64,-0.253061806928333_f64,&
      -4.432644675724012E-002_f64, 6.694610791512273E-002_f64,&
      -0.178047320344777_f64,     -3.318305755877415E-002_f64,&
      -5.599838900861524E-002_f64,-3.573174528242267E-002_f64,&
       0.162930347921219_f64,      7.691262263429682E-002_f64,&
       1.678402998954442E-002_f64, 5.576797665245436E-003_f64,&
       1.812133301369922E-003_f64,-1.253991163614393E-002_f64,&
      -7.920184037933319E-003_f64,-1.899399690115982E-002_f64,&
       2.951325564958146E-002_f64, 6.466038986295269E-002_f64,&
       0.709373309097284_f64,-1.349133247541097_f64,&
       3.048134746740139_f64,-3.481616966416450_f64,&
      -0.326196948109174_f64, 0.200403861518809_f64,&
      -0.260256607115147_f64, 0.507594023544052_f64,&
       0.189944655275759_f64,-0.283705625369913_f64,&
       0.565574191095987_f64,-0.591001158727342_f64,&
      -0.398474381863533_f64, 0.606481805156626_f64,&
      -1.339970830196913_f64, 1.482945080272192_f64,&
      -0.199761710905302_f64,-8.737520935461928E-003_f64,&
       8.828146039181894E-002_f64, 1.193459929569070_f64,&
       6.185167061371060E-002_f64,-1.200244445342679E-002_f64,&
      -2.209483721408903E-003_f64,-0.169184072277413_f64,&
      -4.501696461973124E-002_f64,-2.712546889053091E-003_f64,&
       1.146321492776100E-002_f64, 0.155454234306533_f64,&
       6.858248319521390E-002_f64, 1.049893460440504E-002_f64,&
      -5.784327201693034E-002_f64,-0.372219819151061_f64,&
       8.840791622818779_f64,2.524695831478458_f64,&
       8.847906101171684_f64,6.498710576032950_f64,&
       2.526186512332977_f64,7.336964347919039_f64,&
       6.428840055585406_f64,7.269486701269619_f64,&
       8.855959745221103_f64,6.419617534011105_f64,&
       8.857354636854411_f64,2.536339040131500_f64,&
       6.422372936524658_f64,7.337372688408041_f64,&
       2.518001461375481_f64,7.310365123738160_f64,&
       9.399766789416743_f64,2.340623633213466_f64,&
       9.045643492873985_f64,3.077622553020681_f64,&
       2.268781772580775_f64,7.302960907502583_f64,&
       6.212497757404577_f64,10.16681015553643_f64,&
       9.011049132436610_f64,6.423452203530955_f64,&
       8.837538568283823_f64,1.915766516838633_f64,&
       6.340335188696540_f64,7.333097351147578_f64,&
       2.559599930163381_f64,8.295418962226607_f64]

  efield_ref = [1.09571672509656_f64,-0.660812182548444_f64,&
       0.554548067419379_f64,-0.252278438293157_f64,&
      -1.359654339406665_f64,-5.759357559748912E-002_f64,&
      -8.415186349960571E-002_f64,0.744189134396315_f64,&
       8.745552518004301E-002_f64,0.118654898156689_f64,&
       1.225024303098075_f64,-0.674201486899271_f64,&
      -0.213622564417715_f64,0.721213507983366_f64,&
      -1.200204983093983_f64,-5.993023984790116E-002_f64,&
       1.206894362646646_f64,-0.977552857857918_f64,&
       2.026872786226881_f64,-2.380393357199413_f64,&
      -1.782858301050251_f64,2.206715615943564E-002_f64,&
       3.602711636981939E-002_f64,0.972503216565975_f64,&
       0.138279141932933_f64,1.924388088272622E-002_f64,&
       1.352317721309432_f64,-0.785871005995901_f64,&
      -0.445097520304793_f64,0.953655736974383_f64,&
      -1.550622618607893_f64,0.500407024098771_f64,&
      -1.181255569429968_f64,1.20400252956618_f64,&
       0.382701685319993_f64,8.918140466727103E-002_f64,&
       0.690432642041031_f64,0.239686782681609_f64,&
      -0.533815784703545_f64,-0.622138267722021_f64,&
      -8.245052343850910E-002_f64,-4.38779388981232E-002_f64,&
      -1.307875981763865_f64,1.04545069165269_f64,&
      -0.188330740052329_f64,-0.639690858792130_f64,&
       0.699755953082015_f64,0.249092221382297_f64,&
      -1.019680439017095_f64,0.826217169578618_f64,&
       1.854474380301871_f64,2.780564931917876E-002_f64,&
       0.562931773746294_f64,0.546264768322688_f64,&
      -1.578064060544716_f64,-0.255835052565448_f64,&
      -6.570735459788824E-002_f64,-0.185399368200073_f64,&
      -1.442184131178410_f64,0.698939042624095_f64,&
      -0.300808396406801_f64,-0.42136299256032_f64,&
       0.408834767188967_f64,0.393201544734771_f64,&
      -5.835088572012237E-002_f64,-0.115607380033844_f64,&
       4.855460648785139E-002_f64,3.677374425173193E-003_f64,&
      -1.219706016693366E-002_f64,-4.78682834478023E-002_f64,&
      -5.984338035144897E-002_f64,-0.105580411249307_f64,&
      -2.206506933840884E-002_f64,-1.063060743516663E-002_f64,&
      -3.086159019128434E-002_f64,-7.871234141577716E-003_f64,&
      -1.855288617786560E-003_f64,-3.766265549709074E-002_f64,&
      -0.113810274367247_f64,-7.009074890962005E-002_f64,&
       7.126304902902259E-002_f64,-0.125828895896069_f64,&
       0.722024288185408_f64,5.390668274172328E-002_f64,&
       2.251860466889749E-002_f64,5.15269086684027E-002_f64,&
       0.117427437518344_f64,0.198852815597844_f64,&
       2.906595252763884E-002_f64,-2.51675722700962E-002_f64,&
       4.377375350018507E-002_f64,1.819277864354176E-002_f64,&
      -3.394289108747602E-002_f64,0.102292197918031_f64,&
      -0.127626001562666_f64,7.374937413021989E-002_f64 ]

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
  particle_info_check(:,1) = [5.7984711070298_f64, 2.65687845344_f64, 0.988412182204986_f64, -1.53576085149311_f64, 0.140396236125518_f64, 1.62412057435254_f64, -5.7026946767517_f64]
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
       8.85616842839152_f64,       2.51872345152223_f64,&
       8.85616842839152_f64,       6.42077589443804_f64,&
       2.51872345152223_f64,       7.33707345436654_f64,&
       6.42077589443803_f64,       7.33707345436654_f64,&
       8.85616842839152_f64,       6.42077589443803_f64,&
       8.85616842839152_f64,       2.51872345152223_f64,&
       6.42077589443804_f64,       7.33707345436654_f64,&
       2.51872345152223_f64,       7.33707345436654_f64,&
       8.85616842839151_f64,       2.51872345152223_f64,&
       8.85616842839151_f64,       6.42077589443804_f64,&
       2.51872345152223_f64,       7.33707345436654_f64,&
       6.42077589443804_f64,       7.33707345436654_f64,&
       8.85616842839151_f64,       6.42077589443804_f64,&
       8.85616842839151_f64,       2.51872345152223_f64,&
       6.42077589443804_f64,       7.33707345436654_f64,&
       2.51872345152223_f64,       7.33707345436654_f64]

  efield_ref = [1.075930933791_f64,-0.612675885501125_f64,&
       0.443524626082271_f64,    -0.147232187360209_f64,&
      -1.348427888203207_f64,    -6.502558815309811E-002_f64,&
      -7.337909919984186E-002_f64,0.727285088543657_f64,&
       8.159731489130038E-002_f64,0.128919043324753_f64,&
       1.203795634344262_f64,    -0.654764505547826_f64,&
      -0.201243799752530_f64,     0.699286285069251_f64,&
      -1.150475764129345_f64,    -0.107114208199862_f64,&
       0.982382460098930_f64,    -0.975480249517246_f64,&
       2.122915748511289_f64,    -1.370270472080486_f64,&
      -1.708312691669170_f64,     5.846363943802890E-003_f64,&
       1.011655742532979E-002_f64,0.932802283287550_f64,&
       9.073056017286369E-002_f64,1.967533402383350E-002_f64,&
       1.364725233393067_f64,    -0.715583640577278_f64,&
      -0.372226215868787_f64,     0.964912300087728_f64,&
      -1.617976544115137_f64,     0.265742972883710_f64,&
      -1.199589907436063_f64,     1.252909777650912_f64,&
       0.273927046650591_f64,     9.621848023178933E-002_f64,&
       0.706311808329164_f64,     0.200073958490227_f64,&
      -0.454420354050335_f64,    -0.633791225333427_f64,&
      -8.800026454029169E-002_f64,-3.16462203208923E-002_f64,&
      -1.327454607988727_f64, 1.054957653577138_f64,&
      -0.178269123365297_f64,-0.661790028807759_f64,&
       0.748400428375984_f64, 0.242162578536987_f64,&
      -1.038660308387255_f64, 0.785408997665118_f64,&
       1.953318169079613_f64, 0.179714136856958_f64,&
       0.645492673299711_f64, 0.572931139573802_f64,&
      -1.677458638770615_f64,-0.428274030589535_f64,&
      -7.886701925872492E-002_f64,-0.202628636437152_f64,&
      -1.421003081681347_f64, 0.695072850111172_f64,&
      -0.287512832666219_f64,-0.396164013789282_f64,&
       0.385596064359861_f64, 0.313034530633890_f64,&
      -5.282467479262684E-002_f64, 0.128577507215433_f64,&
      -0.711118053999076_f64,     -9.959891163893815E-002_f64,&
      -2.241510727790141E-002_f64,-5.785108332635274E-002_f64,&
      -9.959891163893770E-002_f64,-0.202357509010883_f64,&
      -2.698172991868412E-002_f64, 2.764012473177656E-002_f64,&
      -5.282467479262696E-002_f64,-2.241510727790143E-002_f64,&
       2.764012473177652E-002_f64,-0.105172882777462_f64,&
       0.128577507215434_f64,     -5.785108332635274E-002_f64,&
       5.282467479262684E-002_f64,-0.128577507215433_f64,&
       0.711118053999076_f64,      9.959891163893815E-002_f64,&
       2.241510727790141E-002_f64, 5.785108332635274E-002_f64,&
       9.959891163893770E-002_f64, 0.202357509010883_f64,&
       2.698172991868412E-002_f64,-2.764012473177656E-002_f64,&
       5.282467479262696E-002_f64, 2.241510727790143E-002_f64,&
      -2.764012473177652E-002_f64, 0.105172882777462_f64,&
      -0.128577507215434_f64,      5.785108332635274E-002_f64 ]

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

  efield_ref = [ 1.08082559657808_f64, -0.612704048104539_f64, &
       0.441157932686102_f64,      -0.172327119498838_f64, &
      -1.348927963142263_f64,      -6.522586594776717E-002_f64, &
      -7.263186742099023E-002_f64,  0.730970903683016_f64, &
       8.240616026313064E-002_f64,  0.129093083410090_f64, &
       1.203011870084057_f64,      -0.658460927650987_f64, &
      -0.203029919297674_f64,       0.699217597306996_f64, &
      -1.148713491919134_f64,      -9.789413838117651E-002_f64, &
       1.243608139400880_f64,      -1.006043849440096_f64, &
       2.096953689240034_f64,      -2.677644125373574_f64, &
       -1.74546504066468_f64,       3.847773765858130E-003_f64, &
       3.325168432004974E-002_f64,  1.145081545675755_f64, &
       0.138020690006099_f64,       2.112200473200536E-002_f64, &
       1.341973350751623_f64,      -0.917056714198851_f64, &
      -0.469225778142738_f64,       0.973048003148078_f64, &
      -1.578492342985654_f64,       0.747000595968360_f64, &
      -1.200537163541495_f64,       1.253869452729790_f64, &
       0.271465115913840_f64,       9.569462723308787E-002_f64, &
       0.705438752374706_f64,       0.198778044646259_f64, &
      -0.450317120408235_f64,      -0.633268062842164_f64, &
      -8.769419897368499E-002_f64, -3.142248623852158E-002_f64, &
      -1.327712595599607_f64,       1.054613948619807_f64, &
      -0.178520244920069_f64,      -0.662138018742490_f64, &
       0.748462865480379_f64,       0.243637108943847_f64, &
      -1.077902148956411_f64,       0.832681996898445_f64, &
       1.839537446119879_f64,       0.146740726060913_f64, &
       0.583901214470475_f64,       0.541392794018474_f64, &
      -1.531596517499853_f64,      -0.375810425586335_f64, &
      -5.701420464705928E-002_f64, -0.199775384141560_f64, &
      -1.427208453064408_f64,       0.667583541694029_f64, &
      -0.307416021504515_f64,      -0.403653241083810_f64, &
       0.389764480073142_f64,       0.391303569729063_f64, &
       0.377859176211884_f64,      -0.633125059737926_f64, &
       0.732664551283754_f64,      -0.412756642152723_f64, &
      -9.039994141733469E-002_f64,  5.310583630754348E-002_f64, &
      -0.257999845059053_f64,      -9.780981261806573E-002_f64, &
       5.994579600592601E-002_f64, -0.134265178421137_f64, &
       0.305444128990836_f64,      -0.108963790645523_f64, &
      -0.162516241575596_f64,       0.254722582559101_f64, &
      -0.639206330991396_f64,       0.117295895733074_f64, &
       5.282467479262700E-002_f64, -0.128577507215434_f64, &
       0.711118053999076_f64,       9.959891163893754E-002_f64, &
       2.241510727790098E-002_f64,  5.785108332635303E-002_f64, &
       9.959891163893706E-002_f64,  0.202357509010884_f64, &
       2.698172991868432E-002_f64, -2.764012473177670E-002_f64, &
       5.282467479262689E-002_f64,  2.241510727790096E-002_f64, &
      -2.764012473177691E-002_f64,  0.105172882777462_f64, &
      -0.128577507215434_f64,       5.785108332635340E-002_f64 ]

  error = maxval(abs(efield-efield_ref))
  if (error> EQV_TOL) then
     passed = .FALSE.
     print*, 'e field wrong in advect_ex.', error
  end if

  ! Compare to reference
  ! Particle information after advect_ev application 
  particle_info_check(:,1) = [5.64516785359875_f64, 2.67294223708412_f64, 1.15090247333432_f64, -1.53275986329561_f64, 0.163149649244948_f64, 1.62568524823416_f64, -5.7026946767517_f64]
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


end program test_time_propagator_pic_3d3v_vm_disgrad_trafo
