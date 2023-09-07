! TODO: Use input from file to initialize and compare
! Unit test for antisymmetric splitting with coordinate transformation
! author: Benedikt Perse, IPP
program test_time_propagator_pic_3d3v_vm_hs_trafo
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

  use sll_m_time_propagator_pic_vm_3d3v_hs_trafo, only: &
       sll_t_time_propagator_pic_vm_3d3v_hs_trafo

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
  type(sll_t_time_propagator_pic_vm_3d3v_hs_trafo) :: propagator

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
     call maxwell_solver%init( domain(:,1:2), num_cells, degree_smoother, map, mass_tolerance=1d-12, poisson_tolerance=1d-12)
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
       domain(:,1), domain(:,3), map, iter_tolerance=1d-12, max_iter=10 )

  call propagator%operatorHp( delta_t )

  ! Compare to reference
  ! Particle information after advect_x application 
  particle_info_check(:,1) = [5.64498281855006_f64, 2.67178099953534_f64, 1.15082423964024_f64, -1.53573785045216_f64, 0.140647612086555_f64, 1.62412057435254_f64, -5.7026946767517_f64] 

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
  
  efield_ref =  [ -0.120865084289336_f64, 6.783520304109252E-002_f64,  &
       0.356341879887054_f64,      -0.3259316990565466_f64,  &
       -0.1472516097809145_f64,        8.838483572458752E-002_f64, &
       1.216151274125419E-002_f64,   5.046670844182286E-002_f64,  &
       -2.383567377030787E-003_f64,  -2.450633646487882E-002_f64, &
       1.3117832260736244E-003_f64,   2.205361523101158E-002_f64, &
       -0.1182338850382063_f64,        1.868057011786124E-002_f64, &
       5.298667111751185E-002_f64,   5.570947444446495E-002_f64, &
       4.237009113287387E-002_f64, -0.3255247678116475_f64,   &
       2.011352577271616_f64,       -2.833086605728608_f64,  &
       -0.5443164671662246_f64,       0.1581277277009513_f64, &
       0.1178002183850051_f64,       0.4658802226710472_f64,  &
       5.359458825601224E-002_f64, -0.1326880553172483_f64,  &
       0.1399678307626886_f64,      -0.2376334687908240_f64, &
       -0.3846502340637258_f64,       0.2925475542765447_f64, &
       -0.3764116948186366_f64,       0.9015139792672735_f64, &
       1.155434947588656E-003_f64,   5.217220177063509E-002_f64,&
       0.3562809734686507_f64,        1.089454193651550E-002_f64, &
       2.485225270170377E-002_f64,   4.528096503325593E-002_f64, &
       -0.2969947441034492_f64,        4.733480420007227E-002_f64, &
       -2.881853935422051E-003_f64, -0.1162413822850134_f64, &
       -0.1259666557943634_f64,      -0.1470975370540288_f64, &
       -2.492842487848107E-002_f64,   1.840489120918593E-002_f64, &
       6.791301452220925E-002_f64,   9.004766421048209E-002_f64, &
       0.1236777818981160_f64,      -0.3689973022840232_f64, &
       1.924656573720140_f64,        6.192054095197477E-002_f64, &
       -9.861996799625569E-002_f64,  0.3918211203530287_f64,  &
       -1.388963029501974_f64,       0.3071082445391841_f64,  &
       2.831673277082568E-002_f64, -0.2854129008994981_f64,  &
       -0.2232990935987929_f64,      -0.5347293320528685_f64, &
       -0.1539163315301638_f64,       0.2770397087540344_f64, &
       -0.2912152277687483_f64,       0.2379429524327899_f64, &
       0.3778169175446466_f64,      -0.6331050592758131_f64,  &
       0.7326980579396500_f64,      -0.4126970884239815_f64, &
       -9.058295447004177E-002_f64,   5.348038950779958E-002_f64, &
       -0.2589743810781840_f64,       -9.773617830251420E-002_f64, &
       6.002813709118699E-002_f64, -0.1344315273383278_f64,  &
       0.3059097346991039_f64,      -0.1090284227410261_f64,  &
       -0.1624962707432857_f64,       0.2547009633851950_f64, &
       -0.6391949384798981_f64,       0.1172548258707867_f64, &
       5.282467479262686E-002_f64, -0.1285775072154342_f64,   &
       0.7111180539990762_f64,        9.959891163893769E-002_f64, &
       2.241510727790113E-002_f64,   5.785108332635313E-002_f64, &
       9.959891163893730E-002_f64,  0.2023575090108842_f64,   &
       2.698172991868425E-002_f64,  -2.764012473177668E-002_f64, &
       5.282467479262686E-002_f64,   2.241510727790119E-002_f64, &
       -2.764012473177673E-002_f64,  0.1051728827774628_f64,   &
       -0.1285775072154344_f64,        5.785108332635330E-002_f64 ]
 

  error = maxval(abs(efield-efield_ref))
  if (error> EQV_TOL) then
     passed = .FALSE.
     print*, 'e field wrong in Hp.', error
  end if

  if (passed .EQV. .FALSE.) then
     print*, 'Error in Hp.'
  end if
  
  !Reset field info
  efield =  [ -0.120865137170533_f64,6.78352642859766E-002_f64,&
       0.356341707929176_f64,-0.325931439937101_f64,&
      -0.147251541668634_f64, 8.838475181435975E-002_f64,&
       1.216165423105994E-002_f64,5.04665559516821E-002_f64,&
      -2.383617053938659E-003_f64,-2.45062745338305E-002_f64,&
       1.311675342687885E-003_f64,2.20537285897379E-002_f64,&
      -0.118233822453996_f64,     1.868050512205642E-002_f64,&
       5.298681400137125E-002_f64,5.570927811688920E-002_f64,&
       4.237014401420807E-002_f64,-0.32552482905586_f64,&
       2.011352749228831_f64,-2.833086864848195_f64,&
      -0.544316535278314_f64, 0.158127811611123_f64,&
       0.117800076894635_f64, 0.465880375161616_f64,&
       5.359463793213738E-002_f64,-0.13268811724853_f64,&
       0.139967938646569_f64,-0.237633582149024_f64,&
      -0.384650296648817_f64, 0.292547619273112_f64,&
      -0.376411837702051_f64, 0.901514175594523_f64,&
       1.155327064178219E-003_f64,5.21723446544482E-002_f64,&
       0.356280801510725_f64,1.089468342631703E-002_f64,&
       2.485236606036780E-002_f64,4.52807687056821E-002_f64,&
      -0.296994484984044_f64,4.73346517098916E-002_f64,&
      -2.881903612283355E-003_f64,-0.11624131970079_f64,&
      -0.125966708675535_f64,-0.1470974689417_f64,&
      -2.492836294739257E-002_f64,1.84048262134212E-002_f64,&
       6.791307576715578E-002_f64,9.00475803002520E-002_f64,&
       0.123677889781972_f64,-0.36899744516748_f64,&
       1.924656745677308_f64,6.19203994616005E-002_f64,&
      -9.862008135451856E-002_f64,0.391821316680281_f64,&
      -1.388963288621601_f64,0.307108397029713_f64,&
       2.831678244699745E-002_f64,-0.28541296348458_f64,&
      -0.223299040717434_f64,-0.53472940016491_f64,&
      -0.153916393461410_f64,0.27703977375064_f64,&
      -0.291215289012902_f64,0.23794303634296_f64,&
       0.377816704412492_f64,-0.63310480467487_f64,&
       0.732697774350068_f64,-0.41269683382304_f64,&
      -9.058279425323573E-002_f64,5.34801941698655E-002_f64,&
      -0.258974126477242_f64,-9.77364179643301E-002_f64,&
       6.002801958761451E-002_f64,-0.13443136712152_f64,&
       0.305909521566950_f64,-0.1090282625242_f64,&
      -0.162496110526479_f64,0.25470072372337_f64,&
      -0.639194683878957_f64,0.11725463053285_f64,&
       5.282467479262681E-002_f64,-0.12857750721543_f64,&
       0.7111180539990760_f64,9.95989116389381E-002_f64,&
       2.241510727790141E-002_f64,5.7851083326352E-002_f64,&
       9.959891163893774E-002_f64,0.20235750901088_f64,&
       2.698172991868413E-002_f64,-2.7640124731776E-002_f64,&
       5.282467479262692E-002_f64,2.2415107277901E-002_f64,&
      -2.764012473177654E-002_f64,0.10517288277746_f64,&
      -0.1285775072154340_f64,5.78510833263527E-002_f64 ]


  call propagator%operatorHB( delta_t )

  efield_ref = [    0.28013348235601_f64,-0.1509129419134_f64,&
       0.44587273017667_f64,     -0.338852982219_f64,&
      -0.54825016119518_f64,      0.101306294096_f64,&
      -7.73693680164418E-002_f64, 0.269214762151_f64,&
       8.71474051935627E-002_f64,-3.742781681642E-002_f64,&
       0.40231029486923_f64,-0.196694477609_f64,&
      -0.20776484470149_f64, 0.237428711321_f64,&
      -0.34801180552517_f64, 6.8630820399485E-002_f64,&
       0.44336876354075_f64,-0.5442730352552_f64,&
       2.10088377147633_f64,-2.8460084071307_f64,&
      -0.94531515480486_f64, 0.1710493538937_f64,&
       2.82690546471351E-002_f64,0.684628581361_f64,&
       0.14312566017963_f64,-0.1456096595311_f64,&
       0.54096655817311_f64,-0.4563817883484_f64,&
      -0.47418131889631_f64, 0.5112958254725_f64,&
      -0.77741045722860_f64, 0.9144357178771_f64,&
      -0.39984329246237_f64, 0.4531709641809_f64,&
       0.26674977926322_f64, 0.1004257056738_f64,&
       0.24360057225976_f64, 3.23592264230853E-002_f64,&
      -0.28407294270144_f64,-0.17141355448950_f64,&
      -9.24129258597849E-002_f64,-2.6710297453297E-002_f64,&
      -0.52696532820208_f64,      0.25390115058484_f64,&
      -1.20068206647965E-002_f64,-0.20034337998597_f64,&
       0.28666128196655_f64,7.71260380176559E-002_f64,&
      -0.27732072974457_f64,3.20011743590639E-002_f64,&
       1.83512572342980_f64,0.1514514217091_f64,&
       0.12012812484488_f64,0.37889977439768_f64,&
      -1.37604174633900_f64,8.83601908303128E-002_f64,&
      -6.12142398005033E-002_f64,-0.195881941237_f64,&
      -0.62429766024398_f64,-0.13373078063836_f64,&
      -0.14099485117881_f64, 5.82915675512412E-002_f64,&
      -7.24670828135012E-002_f64,0.2250214940603_f64,&
       0.37781670441249_f64,-0.63310480467487_f64,&
       0.73269777435006_f64,-0.41269683382304_f64,&
      -9.05827942532357E-002_f64,5.34801941698655E-002_f64,&
      -0.25897412647724_f64,-9.77364179643301E-002_f64,&
       6.00280195876145E-002_f64,-0.1344313671215_f64,&
       0.30590952156695_f64,-0.10902826252422_f64,&
      -0.16249611052647_f64, 0.25470072372337_f64,&
      -0.63919468387895_f64, 0.11725463053285_f64,&
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
     print*, 'e field wrong in HB.', error
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
  particle_info_check(:,1) = [5.7984711070298_f64, 2.65687845344_f64, 0.988412182204986_f64, -1.52589193591242_f64, 0.159503586388184_f64, 1.62120606580492_f64, -5.7026946767517_f64]
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

  bfield_ref =  [ -0.240629638519148_f64,0.439356169323687_f64,&
      -0.862432172124927_f64,0.201775442535300_f64,&
       0.212054288949269_f64,-0.343942109132815_f64,&
       0.615062521058436_f64,-0.177938915407448_f64,&
      -6.64840627481964E-002_f64,0.108998953273312_f64,&
      -0.206486992809297_f64,8.204312406859801E-002_f64,&
       0.114807258148441_f64,-0.207379825845404_f64,&
       0.449867355134374_f64,-0.120092248431370_f64,&
      -7.681407266202322E-003_f64, 9.266198032771667E-003_f64,&
      -2.220303565248733E-002_f64,-6.493988117977517E-003_f64,&
      -1.253066247708736E-002_f64,-5.263326621794584E-003_f64,&
       2.621389621654364E-002_f64, 1.085131011518627E-002_f64,&
       4.413088155543229E-003_f64, 3.621544664942466E-004_f64,&
      -7.567716698553939E-004_f64,-5.549425551449105E-003_f64,&
      -3.948864242619282E-003_f64,-1.398213496251421E-003_f64,&
       7.351998472127530E-004_f64, 1.540470078916111E-002_f64,&
       0.348852471531161_f64,-0.48304062230331_f64,&
       0.877323239869907_f64,-0.95958892825146_f64,&
      -7.655154923749796E-002_f64, 7.15738073285931E-002_f64,&
      -0.103854043726127_f64,      0.14757784724715_f64,&
       7.881816384194895E-002_f64,-9.94201232265951E-002_f64,&
       0.203867608136164_f64,-0.21791257578422_f64,&
      -0.165183591262697_f64, 0.22165215653015_f64,&
      -0.443437893381619_f64, 0.47174070526025_f64,&
      -5.135675097547273E-002_f64, 6.11114586514377E-003_f64,&
       4.876016225873103E-003_f64, 0.256823428038163_f64,&
       7.436038028742864E-003_f64, 2.257784600278867E-004_f64,&
      -4.428553207681396E-003_f64,-4.197932489320838E-002_f64,&
      -9.369001940901969E-003_f64,-2.123733172435943E-004_f64,&
       4.454667148985291E-003_f64, 3.977363514186221E-002_f64,&
       1.908681161571242E-002_f64,-1.648219826515429E-003_f64,&
      -7.620425656474161E-003_f64,-9.458954327881197E-002_f64,&
       6.678452237257067_f64,5.786642943228285_f64,&
       6.675307595427438_f64,6.186721415567674_f64,&
       5.785826199059394_f64,6.468569539918412_f64,&
       6.200461335552152_f64,6.481348649643123_f64,&
       6.675869964312938_f64,6.201410611451637_f64,&
       6.675159184633373_f64,5.784475019760513_f64,&
       6.200873550694542_f64,6.468462542163213_f64,&
       5.788254602240805_f64,6.473129523962806_f64,&
       6.821714200735887_f64,5.737229001247018_f64,&
       6.713253179033261_f64,5.452477377864704_f64,&
       5.715004566235510_f64,6.465805603018059_f64,&
       6.156116028742584_f64,7.109679327708581_f64,&
       6.689555016838243_f64,6.210388782384278_f64,&
       6.659630596192739_f64,5.630554407453559_f64,&
       6.182669053644875_f64,6.466232933689021_f64,&
       5.808137961164795_f64,6.712516878920258_f64]

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

end program test_time_propagator_pic_3d3v_vm_hs_trafo
