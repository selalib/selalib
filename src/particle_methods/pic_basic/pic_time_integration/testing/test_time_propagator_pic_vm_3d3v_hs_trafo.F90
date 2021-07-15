! TODO: Use input from file to initialize and compare
! Unit test for antisymmetric splitting with coordinate transformation
! author: Benedikt Perse, IPP
program test_time_propagator_pic_3d3v_vm_hs_trafo
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
  particle_info_check(:,1) = [5.6449828185500666_f64, 2.6717809995353443_f64, 1.1508242396402408_f64, -1.5357378504521602_f64, 0.14064761208655593_f64, 1.6241205743525460_f64, -5.7026946767517002_f64] 

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
  
  efield_ref =  [ -0.12086508428933694_f64, 6.7835203041092529E-002_f64,  &
       0.35634187988705468_f64,      -0.32593169905654662_f64,  &
       -0.14725160978091453_f64,        8.8384835724587524E-002_f64, &
       1.2161512741254191E-002_f64,   5.0466708441822869E-002_f64,  &
       -2.3835673770307873E-003_f64,  -2.4506336464878820E-002_f64, &
       1.3117832260736244E-003_f64,   2.2053615231011580E-002_f64, &
       -0.11823388503820638_f64,        1.8680570117861242E-002_f64, &
       5.2986671117511853E-002_f64,   5.5709474444464953E-002_f64, &
       4.2370091132873877E-002_f64, -0.32552476781164752_f64,   &
       2.0113525772716163_f64,       -2.8330866057286084_f64,  &
       -0.54431646716622462_f64,       0.15812772770095138_f64, &
       0.11780021838500510_f64,       0.46588022267104723_f64,  &
       5.3594588256012245E-002_f64, -0.13268805531724839_f64,  &
       0.13996783076268865_f64,      -0.23763346879082406_f64, &
       -0.38465023406372589_f64,       0.29254755427654477_f64, &
       -0.37641169481863668_f64,       0.90151397926727350_f64, &
       1.1554349475886561E-003_f64,   5.2172201770635092E-002_f64,&
       0.35628097346865079_f64,        1.0894541936515508E-002_f64, &
       2.4852252701703777E-002_f64,   4.5280965033255932E-002_f64, &
       -0.29699474410344923_f64,        4.7334804200072279E-002_f64, &
       -2.8818539354220516E-003_f64, -0.11624138228501348_f64, &
       -0.12596665579436345_f64,      -0.14709753705402889_f64, &
       -2.4928424878481070E-002_f64,   1.8404891209185936E-002_f64, &
       6.7913014522209256E-002_f64,   9.0047664210482090E-002_f64, &
       0.12367778189811601_f64,      -0.36899730228402328_f64, &
       1.9246565737201402_f64,        6.1920540951974773E-002_f64, &
       -9.8619967996255692E-002_f64,  0.39182112035302874_f64,  &
       -1.3889630295019748_f64,       0.30710824453918417_f64,  &
       2.8316732770825689E-002_f64, -0.28541290089949811_f64,  &
       -0.22329909359879299_f64,      -0.53472933205286854_f64, &
       -0.15391633153016387_f64,       0.27703970875403444_f64, &
       -0.29121522776874831_f64,       0.23794295243278998_f64, &
       0.37781691754464669_f64,      -0.63310505927581318_f64,  &
       0.73269805793965004_f64,      -0.41269708842398156_f64, &
       -9.0582954470041777E-002_f64,   5.3480389507799583E-002_f64, &
       -0.25897438107818405_f64,       -9.7736178302514201E-002_f64, &
       6.0028137091186998E-002_f64, -0.13443152733832786_f64,  &
       0.30590973469910393_f64,      -0.10902842274102616_f64,  &
       -0.16249627074328574_f64,       0.25470096338519504_f64, &
       -0.63919493847989817_f64,       0.11725482587078673_f64, &
       5.2824674792626869E-002_f64, -0.12857750721543426_f64,   &
       0.71111805399907624_f64,        9.9598911638937693E-002_f64, &
       2.2415107277901136E-002_f64,   5.7851083326353137E-002_f64, &
       9.9598911638937304E-002_f64,  0.20235750901088428_f64,   &
       2.6981729918684250E-002_f64,  -2.7640124731776683E-002_f64, &
       5.2824674792626869E-002_f64,   2.2415107277901192E-002_f64, &
       -2.7640124731776738E-002_f64,  0.10517288277746289_f64,   &
       -0.12857750721543446_f64,        5.7851083326353303E-002_f64 ]
 

  error = maxval(abs(efield-efield_ref))
  if (error> EQV_TOL) then
     passed = .FALSE.
     print*, 'e field wrong in Hp.', error
  end if

  if (passed .EQV. .FALSE.) then
     print*, 'Error in Hp.'
  end if
  
  !Reset field info
  efield =  [ -0.12086513717053321_f64,6.7835264285976665E-002_f64,&
       0.35634170792917602_f64,-0.32593143993710183_f64,&
       -0.14725154166863458_f64,8.8384751814359758E-002_f64,&
       1.2161654231059946E-002_f64,5.0466555951682146E-002_f64,&
       -2.3836170539386593E-003_f64,-2.4506274533830592E-002_f64,&
       1.3116753426878850E-003_f64,2.2053728589737979E-002_f64,&
       -0.11823382245399601_f64,1.8680505122056424E-002_f64,&
       5.2986814001371257E-002_f64,5.5709278116889208E-002_f64,&
       4.2370144014208078E-002_f64,-0.32552482905586400_f64,&
       2.0113527492288310_f64,-2.8330868648481955_f64,&
       -0.54431653527831414_f64,0.15812781161112399_f64,&
       0.11780007689463567_f64,0.46588037516161651_f64,&
       5.3594637932137382E-002_f64,-0.13268811724853505_f64,&
       0.13996793864656948_f64,-0.23763358214902461_f64,&
       -0.38465029664881722_f64,0.29254761927311207_f64,&
       -0.37641183770205133_f64,0.90151417559452318_f64,&
       1.1553270641782199E-003_f64,5.2172344654448290E-002_f64,&
       0.35628080151072533_f64,1.0894683426317034E-002_f64,&
       2.4852366060367809E-002_f64,4.5280768705682171E-002_f64,&
       -0.29699448498404407_f64,4.7334651709891623E-002_f64,&
       -2.8819036122833559E-003_f64,-0.11624131970079891_f64,&
       -0.12596670867553580_f64,-0.14709746894170217_f64,&
       -2.4928362947392573E-002_f64,1.8404826213421276E-002_f64,&
       6.7913075767155787E-002_f64,9.0047580300252034E-002_f64,&
       0.12367788978197283_f64,-0.36899744516748523_f64,&
       1.9246567456773085_f64,6.1920399461600557E-002_f64,&
       -9.8620081354518566E-002_f64,0.39182131668028114_f64,&
       -1.3889632886216019_f64,0.30710839702971365_f64,&
       2.8316782446997455E-002_f64,-0.28541296348458572_f64,&
       -0.22329904071743439_f64,-0.53472940016491177_f64,&
       -0.15391639346141034_f64,0.27703977375064226_f64,&
       -0.29121528901290250_f64,0.23794303634296035_f64,&
       0.37781670441249299_f64,-0.63310480467487196_f64,&
       0.73269777435006866_f64,-0.41269683382304079_f64,&
       -9.0582794253235732E-002_f64,5.3480194169865519E-002_f64,&
       -0.25897412647724283_f64,-9.7736417964330177E-002_f64,&
       6.0028019587614515E-002_f64,-0.13443136712152157_f64,&
       0.30590952156695028_f64,-0.10902826252422022_f64,&
       -0.16249611052647969_f64,0.25470072372337893_f64,&
       -0.63919468387895784_f64,0.11725463053285312_f64,&
       5.2824674792626813E-002_f64,-0.12857750721543390_f64,&
       0.71111805399907602_f64,9.9598911638938192E-002_f64,&
       2.2415107277901414E-002_f64,5.7851083326352748E-002_f64,&
       9.9598911638937748E-002_f64,0.20235750901088362_f64,&
       2.6981729918684139E-002_f64,-2.7640124731776600E-002_f64,&
       5.2824674792626924E-002_f64,2.2415107277901414E-002_f64,&
       -2.7640124731776544E-002_f64,0.10517288277746234_f64,&
       -0.12857750721543401_f64,5.7851083326352748E-002_f64 ]


  call propagator%operatorHB( delta_t )

  efield_ref = [    0.28013348235601659_f64,-0.15091294191342475_f64,&
       0.44587273017667750_f64,-0.33885298221969851_f64,&
       -0.54825016119518477_f64,0.10130629409695655_f64,&
       -7.7369368016441883E-002_f64,0.26921476215108381_f64,&
       8.7147405193562719E-002_f64,-3.7427816816427359E-002_f64,&
       0.40231029486923786_f64,-0.19669447760966349_f64,&
       -0.20776484470149731_f64,0.23742871132145793_f64,&
       -0.34801180552517880_f64,6.8630820399485870E-002_f64,&
       0.44336876354075749_f64,-0.54427303525526505_f64,&
       2.1008837714763318_f64,-2.8460084071307916_f64,&
       -0.94531515480486350_f64,0.17104935389371978_f64,&
       2.8269054647135101E-002_f64,0.68462858136101745_f64,&
       0.14312566017963801_f64,-0.14560965953113086_f64,&
       0.54096655817311845_f64,-0.45638178834842552_f64,&
       -0.47418131889631798_f64,0.51129582547251284_f64,&
       -0.77741045722860069_f64,0.91443571787711919_f64,&
       -0.39984329246237155_f64,0.45317096418099806_f64,&
       0.26674977926322413_f64,0.10042570567381820_f64,&
       0.24360057225976964_f64,3.2359226423085363E-002_f64,&
       -0.28407294270144762_f64,-0.17141355448950943_f64,&
       -9.2412925859784903E-002_f64,-2.6710297453297352E-002_f64,&
       -0.52696532820208553_f64,0.25390115058484769_f64,&
       -1.2006820664796542E-002_f64,-0.20034337998597951_f64,&
       0.28666128196655688_f64,7.7126038017655996E-002_f64,&
       -0.27732072974457644_f64,3.2001174359063932E-002_f64,&
       1.8351257234298075_f64,0.15145142170910150_f64,&
       0.12012812484488201_f64,0.37889977439768574_f64,&
       -1.3760417463390060_f64,8.8360190830312810E-002_f64,&
       -6.1214239800503306E-002_f64,-0.19588194123708483_f64,&
       -0.62429766024398414_f64,-0.13373078063836225_f64,&
       -0.14099485117881433_f64,5.8291567551241243E-002_f64,&
       -7.2467082813501213E-002_f64,0.22502149406036426_f64,&
       0.37781670441249299_f64,-0.63310480467487196_f64,&
       0.73269777435006866_f64,-0.41269683382304079_f64,&
       -9.0582794253235732E-002_f64,5.3480194169865519E-002_f64,&
       -0.25897412647724283_f64,-9.7736417964330177E-002_f64,&
       6.0028019587614515E-002_f64,-0.13443136712152157_f64,&
       0.30590952156695028_f64,-0.10902826252422022_f64,&
       -0.16249611052647969_f64,0.25470072372337893_f64,&
       -0.63919468387895784_f64,0.11725463053285312_f64,&
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
  particle_info_check(:,1) = [5.7984711070298012_f64, 2.6568784534400081_f64, 0.98841218220498617_f64, -1.5258919359124290_f64, 0.15950358638818404_f64, 1.6212060658049259_f64, -5.7026946767517002_f64]
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

  bfield_ref =  [ -0.24062963851914809_f64,0.43935616932368715_f64,&
       -0.86243217212492740_f64,0.20177544253530089_f64,&
       0.21205428894926903_f64,-0.34394210913281509_f64,&
       0.61506252105843628_f64,-0.17793891540744869_f64,&
       -6.6484062748196418E-002_f64,0.10899895327331233_f64,&
       -0.20648699280929753_f64,8.2043124068598011E-002_f64,&
       0.11480725814844124_f64,-0.20737982584540438_f64,&
       0.44986735513437487_f64,-0.12009224843137099_f64,&
       -7.6814072662023227E-003_f64,9.2661980327716671E-003_f64,&
       -2.2203035652487338E-002_f64,-6.4939881179775175E-003_f64,&
       -1.2530662477087366E-002_f64,-5.2633266217945842E-003_f64,&
       2.6213896216543645E-002_f64,1.0851310115186275E-002_f64,&
       4.4130881555432292E-003_f64,3.6215446649424666E-004_f64,&
       -7.5677166985539397E-004_f64,-5.5494255514491053E-003_f64,&
       -3.9488642426192821E-003_f64,-1.3982134962514216E-003_f64,&
       7.3519984721275304E-004_f64,1.5404700789161119E-002_f64,&
       0.34885247153116167_f64,-0.48304062230331407_f64,&
       0.87732323986990723_f64,-0.95958892825146236_f64,&
       -7.6551549237497962E-002_f64,7.1573807328593156E-002_f64,&
       -0.10385404372612794_f64,0.14757784724715181_f64,&
       7.8818163841948952E-002_f64,-9.9420123226595145E-002_f64,&
       0.20386760813616489_f64,-0.21791257578422060_f64,&
       -0.16518359126269730_f64,0.22165215653015446_f64,&
       -0.44343789338161915_f64,0.47174070526025108_f64,&
       -5.1356750975472737E-002_f64,6.1111458651437744E-003_f64,&
       4.8760162258731036E-003_f64,0.25682342803816349_f64,&
       7.4360380287428647E-003_f64,2.2577846002788671E-004_f64,&
       -4.4285532076813964E-003_f64,-4.1979324893208381E-002_f64,&
       -9.3690019409019692E-003_f64,-2.1237331724359432E-004_f64,&
       4.4546671489852918E-003_f64,3.9773635141862211E-002_f64,&
       1.9086811615712421E-002_f64,-1.6482198265154293E-003_f64,&
       -7.6204256564741616E-003_f64,-9.4589543278811972E-002_f64,&
       6.6784522372570674_f64,5.7866429432282853_f64,&
       6.6753075954274381_f64,6.1867214155676749_f64,&
       5.7858261990593940_f64,6.4685695399184127_f64,&
       6.2004613355521521_f64,6.4813486496431238_f64,&
       6.6758699643129384_f64,6.2014106114516379_f64,&
       6.6751591846333733_f64,5.7844750197605137_f64,&
       6.2008735506945429_f64,6.4684625421632136_f64,&
       5.7882546022408050_f64,6.4731295239628066_f64,&
       6.8217142007358875_f64,5.7372290012470186_f64,&
       6.7132531790332619_f64,5.4524773778647040_f64,&
       5.7150045662355105_f64,6.4658056030180591_f64,&
       6.1561160287425842_f64,7.1096793277085819_f64,&
       6.6895550168382432_f64,6.2103887823842783_f64,&
       6.6596305961927396_f64,5.6305544074535598_f64,&
       6.1826690536448750_f64,6.4662329336890219_f64,&
       5.8081379611647952_f64,6.7125168789202583_f64]

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
