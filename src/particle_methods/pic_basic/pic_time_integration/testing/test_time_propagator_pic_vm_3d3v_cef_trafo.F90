! TODO: Use input from file to initialize and compare
! Unit test for antisymmetric splitting with coordinate transformation
! author: Benedikt Perse, IPP
program test_time_propagator_pic_3d3v_vm_cef_trafo
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

  use sll_m_time_propagator_pic_vm_3d3v_cef_trafo, only: &
       sll_t_time_propagator_pic_vm_3d3v_cef_trafo

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
  sll_real64, allocatable :: phi(:), efield(:), efield_ref(:)
  sll_real64, allocatable :: bfield(:), bfield_ref(:)
  sll_real64, allocatable :: rho(:), rho_local(:)

  ! Abstract kernel smoothers
  type(sll_t_particle_mesh_coupling_spline_3d_feec), pointer :: particle_mesh_coupling     


  ! Maxwell solver 
  ! Abstract 
  class(sll_c_maxwell_3d_base), pointer :: maxwell_solver

  ! Specific Hamiltonian splitting
  type(sll_t_time_propagator_pic_vm_3d3v_cef_trafo) :: propagator

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
  particle_info_ref = reshape([0.9375_f64, 0.4375_f64, &
       0.15731068461017_f64,     -1.53412054435254_f64, &
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
  allocate( sll_t_maxwell_3d_trafo :: maxwell_solver )
  select type ( maxwell_solver )
  type is ( sll_t_maxwell_3d_trafo )
     call maxwell_solver%init( domain(:,1:2), num_cells, degree_smoother, map, mass_tolerance=1d-12, poisson_tolerance=1d-12 )
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
  particle_info_check(:,1) = [5.64506342246758_f64, 2.67261389177405_f64, 1.15082423964024_f64, -1.53412054435254_f64, 0.15731068461017_f64, 1.62412057435254_f64, -5.7026946767517_f64] 

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

  efield_ref = [  -0.12087219726697_f64,        6.783539486983726E-002_f64, &
       0.356353641468976_f64,       -0.325896418030751_f64, &
      -0.147239027534311_f64,        8.837297926973777E-002_f64, &
       1.216501305789119E-002_f64,   5.043783009487787E-002_f64, &
      -2.391072551433183E-003_f64,  -2.450249432807672E-002_f64, &
       1.317288014952872E-003_f64,   2.207479029499933E-002_f64, &
      -0.118229818709040_f64,        1.867946748761371E-002_f64, &
       5.298115140041303E-002_f64,   5.569244798569784E-002_f64, &
       4.213220633693237E-002_f64,  -0.325514392321000_f64, &
       2.011891366614528_f64,       -2.832207714099234_f64, &
      -0.543939127711864_f64,        0.157637216632894_f64, &
       0.117989496314166_f64,        0.464995293840000_f64, &
       5.334402798136481E-002_f64,  -0.132532172602102_f64, &
       0.140178569645546_f64,       -0.236905431172573_f64, &
      -0.384522844669711_f64,        0.292520684967407_f64, &
      -0.376661875084894_f64,        0.901020349610536_f64, &
       1.156284596575764E-003_f64,   5.217440013281924E-002_f64, &
       0.356271120801042_f64,        1.089622571168057E-002_f64, &
       2.488811706214149E-002_f64,   4.521226145723579E-002_f64, &
      -0.296803854535637_f64,        4.728739335850543E-002_f64, &
      -2.892020703389962E-003_f64,  -0.116227050208028_f64, &
      -0.126004419413018_f64,       -0.147083542800179_f64, &
      -2.492631769043163E-002_f64,   1.840155349977320E-002_f64, &
       6.792122290640563E-002_f64,   9.004268638402947E-002_f64, &
       0.123742661803018_f64,       -0.368990773815614_f64, &
       1.924387860342426_f64,        6.192602474065532E-002_f64, &
      -9.722596265475513E-002_f64,   0.388940475246340_f64, &
      -1.381100562316837_f64,        0.305426277361021_f64, &
       2.794407561654449E-002_f64,  -0.284813203687078_f64, &
      -0.224879789502805_f64,       -0.534295139532499_f64, &
      -0.153854217760446_f64,        0.276924580851949_f64, &
      -0.290900135871270_f64,        0.237802674604796_f64, &
       0.377718442494688_f64,       -0.632861733920966_f64, &
       0.732195132710310_f64,       -0.412669202048753_f64, &
      -9.042790310289820E-002_f64,   5.316987304385317E-002_f64, &
      -0.258214180724428_f64,       -9.782865395225069E-002_f64, &
       5.994162338028012E-002_f64,  -0.134255759824019_f64, &
       0.305452400790688_f64,       -0.108957294755565_f64, &
      -0.162452432414225_f64,        0.254593882632091_f64, &
      -0.638947153027963_f64,        0.117238881701623_f64, &
       5.282467479262686E-002_f64,  -0.128577507215434_f64, &
       0.711118053999076_f64,        9.959891163893769E-002_f64, &
       2.241510727790113E-002_f64,   5.785108332635313E-002_f64, &
       9.959891163893730E-002_f64,   0.202357509010884_f64, &
       2.698172991868425E-002_f64,  -2.764012473177668E-002_f64, &
       5.282467479262686E-002_f64,   2.241510727790119E-002_f64, &
      -2.764012473177673E-002_f64,   0.105172882777462_f64, &
      -0.128577507215434_f64,        5.785108332635330E-002_f64  ]
  
  error = maxval(abs(efield-efield_ref))
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
  particle_info_check(:,1) = [5.79847110702980_f64, 2.65687845344_f64, 0.988412182204986_f64, -1.53576085149311_f64, 0.140396236125518_f64, 1.62412057435254_f64, -5.7026946767517_f64]
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

  efield_ref = [   0.2801263672311972_f64,-0.1509127475979201_f64,&
       0.4458844847764738_f64,-0.3388176906726732_f64,&
       -0.5482375761829695_f64,0.1012944342350359_f64,&
       -7.736586195477324E-002_f64,0.2691858776124547_f64,&
       8.713989800210625E-002_f64,-3.742397216498618E-002_f64,&
       0.4023157952776344_f64,-0.1966732979428942_f64,&
       -0.2077607758311584_f64,0.2374276060521199_f64,&
       -0.3480173194406627_f64,6.861378596908489E-002_f64,&
       0.4431308808919948_f64,-0.5442626622513776_f64,&
       2.101422567801368_f64,-2.845129526022647_f64,&
       -0.9449378181161147_f64,0.1705588462327338_f64,&
       2.845832683126435E-002_f64,0.6837436587216543_f64,&
       0.1428751019220445_f64,-0.1454537793306237_f64,&
       0.5411773014364592_f64,-0.4556537553329561_f64,&
       -0.4740539320434760_f64,0.5112689588024655_f64,&
       -0.7776606432964734_f64,0.9139420961920156_f64,&
       -0.3998424471938668_f64,0.4531731683447969_f64,&
       0.2667399196134904_f64,0.1004273951940146_f64,&
       0.2436364412229886_f64,3.229051487543151E-002_f64,&
       -0.2838820426124053_f64,-0.1714609715227601_f64,&
       -9.242309464480685E-002_f64,-2.669596283513925E-002_f64,&
       -0.5270030939679192_f64,0.2539151476043087_f64,&
       -1.200471096210808E-002_f64,-0.2003467203344826_f64,&
       0.2866694928375130_f64,7.712105678413262E-002_f64,&
       -0.2772558454591920_f64,3.200769702585839E-002_f64,&
       1.834857017034219_f64,0.1514568997527508_f64,&
       0.1215221255836012_f64,0.3760191372626315_f64,&
       -1.368179289675097_f64,8.667822984383404E-002_f64,&
       -6.158689493773050E-002_f64,-0.1952822465658381_f64,&
       -0.6258783540008178_f64,-0.1332965908836054_f64,&
       -0.1409327399237359_f64,5.817644228824629E-002_f64,&
       -7.215199340278304E-002_f64,0.2248812196394413_f64,&
       0.3777182207085665_f64,-0.6328614689822679_f64,&
       0.7321948376059241_f64,-0.4126689371100556_f64,&
       -9.042773638068581E-002_f64,5.316966977446250E-002_f64,&
       -0.2582139157857291_f64,-9.782890334523987E-002_f64,&
       5.994150110561960E-002_f64,-0.1342555931018072_f64,&
       0.3054521790045657_f64,-0.1089571280333536_f64,&
       -0.1624522656920130_f64,0.2545936332391019_f64,&
       -0.6389468880892644_f64,0.1172386784322332_f64,&
       5.282467479262681E-002_f64,-0.128577507215433_f64,&
       0.7111180539990760_f64,9.959891163893819E-002_f64,&
       2.241510727790141E-002_f64,5.785108332635274E-002_f64,&
       9.959891163893774E-002_f64,0.2023575090108836_f64,&
       2.698172991868413E-002_f64,-2.76401247317766E-002_f64,&
       5.282467479262692E-002_f64,2.241510727790141E-002_f64,&
       -2.764012473177654E-002_f64,0.1051728827774634_f64,&
       -0.1285775072154340_f64,5.785108332635274E-002_f64 ]

  error = maxval(abs(efield-efield_ref))
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
  particle_info_check(:,1) = [5.798471107029801_f64, 2.656878453440008_f64, 0.9884121822049861_f64, -1.525889225523391_f64, 0.159470737623141_f64, 1.621207419677591_f64, -5.7026946767517_f64]
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

  bfield_ref =  [ -0.240585514907166_f64,0.4392151351523356_f64,&
       -0.8620801097622212_f64,0.2017571453051683_f64,&
       0.2116812459635784_f64,-0.3431581799801322_f64,&
       0.6130229507691998_f64,-0.1775638537792451_f64,&
       -6.631493493593744E-002_f64,0.1086873618966477_f64,&
       -0.2056913859095383_f64,8.189363757282835E-002_f64,&
       0.1147431125113786_f64,-0.2072443230609094_f64,&
       0.4495239240855913_f64,-0.1200303551572965_f64,&
       -7.668599462826376E-003_f64,9.267061733370786E-003_f64,&
       -2.225480500165826E-002_f64,-6.493230413286946E-003_f64,&
       -1.225903612198732E-002_f64,-5.825711739274652E-003_f64,&
       2.774820753151678E-002_f64,1.052440132454065E-002_f64,&
       4.340590885102181E-003_f64,4.792264771119659E-004_f64,&
       -1.065357268055322E-003_f64,-5.465387004389954E-003_f64,&
       -3.936863932141293E-003_f64,-1.420570479149796E-003_f64,&
       7.965755551651466E-004_f64,1.537764215168121E-002_f64,&
       0.3487557658596084_f64,-0.4829018588070253_f64,&
       0.8771301392402559_f64,-0.9592078769563867_f64,&
       -7.637958160080743E-002_f64,7.129184486159892E-002_f64,&
       -0.1033885964668691_f64,0.1470655611980356_f64,&
       7.870649243957696E-002_f64,-9.928479911609827E-002_f64,&
       0.2036554100743141_f64,-0.2175598142931801_f64,&
       -0.165135008892162_f64,0.2215866301225151_f64,&
       -0.4433448733025087_f64,0.4715398886531853_f64,&
       -5.131059747068407E-002_f64,6.109110127467210E-003_f64,&
       4.770607880824950E-003_f64,0.2566547101259397_f64,&
       7.363087693436166E-003_f64,3.215080198409631E-004_f64,&
       -4.465706432173519E-003_f64,-4.180811727306157E-002_f64,&
       -9.320391727674570E-003_f64,-2.427804270567823E-004_f64,&
       4.413618577996447E-003_f64,3.963226447212218E-002_f64,&
       1.906214801921178E-002_f64,-1.643067546373566E-003_f64,&
       -7.57149122599640E-003_f64,-9.449422582787145E-002_f64,&
       6.678448101359681_f64,5.786642919504104_f64,&
       6.675319328358963_f64,6.186737726290672_f64,&
       5.785800764715619_f64,6.468606550451792_f64,&
       6.200354191482222_f64,6.481418306057779_f64,&
       6.675871593753262_f64,6.201407091895710_f64,&
       6.675180822525661_f64,5.784474340328555_f64,&
       6.200875344744776_f64,6.468462742215378_f64,&
       5.788245576023468_f64,6.473119515165730_f64,&
       6.821544330438551_f64,5.737267241764028_f64,&
       6.713678863615379_f64,5.452916705206308_f64,&
       5.714020269280435_f64,6.467315105901619_f64,&
       6.151678981566636_f64,7.112791573269734_f64,&
       6.689626596816499_f64,6.210258397605485_f64,&
       6.660511339995656_f64,5.630393636310857_f64,&
       6.182739277418648_f64,6.466230729548029_f64,&
       5.807781503562824_f64,6.712210362572685_f64]

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

end program test_time_propagator_pic_3d3v_vm_cef_trafo
