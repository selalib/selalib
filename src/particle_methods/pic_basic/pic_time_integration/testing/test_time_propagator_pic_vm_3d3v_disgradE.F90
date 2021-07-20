! TODO: Use input from file to initialize and compare
! Unit test for antisymmetric splitting with coordinate transformation
! author: Benedikt Perse, IPP
program test_time_propagator_pic_3d3v_vm_disgradE
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

  use sll_m_time_propagator_pic_vm_3d3v_helper, only: &
       sll_t_time_propagator_pic_vm_3d3v_helper

  use sll_m_particle_mesh_coupling_spline_3d_feec, only: &
       sll_t_particle_mesh_coupling_spline_3d_feec

  use sll_m_maxwell_3d_base, only: &
       sll_c_maxwell_3d_base

  use sll_m_maxwell_3d_fem_fft, only: &
       sll_t_maxwell_3d_fem_fft

  use sll_m_maxwell_3d_fem, only: &
       sll_t_maxwell_3d_fem

  use sll_m_particle_group_3d3v, only: &
       sll_t_particle_group_3d3v

  use sll_m_particle_group_base, only: &
       sll_t_particle_array

  use sll_m_fft_filter_3d, only: &
       sll_t_fft_filter_3d

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
  type(sll_t_time_propagator_pic_vm_3d3v_helper) :: propagator

  type(sll_t_fft_filter_3d) :: filter

  ! Parameters
  sll_int32  :: n_particles
  sll_real64 :: eta_min(3), eta_max(3), domain(3,3)
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

  domain(:,1) = eta_min
  domain(:,2) = eta_max
  domain(:,3) = eta_max - eta_min

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
  particle_info_ref = reshape([5.937500000000000_f64,        3.437500000000000_f64,       0.65731068461017067_f64,     -1.5341205443525459_f64,    0.15731068461017067_f64,       1.6241205743525459_f64,   5.7026946767517002_f64], [7, n_particles])

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
  call particle_mesh_coupling%init( num_cells, domain(:,1:2), &
       degree_smoother, n_particles)

  ! Initialize Maxwell solver
  allocate( sll_t_maxwell_3d_fem :: maxwell_solver )
  select type ( maxwell_solver )
  type is ( sll_t_maxwell_3d_fem )
     call maxwell_solver%init( domain(:,1:2), num_cells, degree_smoother )
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
  bfield(2*n_total+1:3*n_total) = 1.5_f64

  call propagator%init( maxwell_solver, &
       particle_mesh_coupling, particle_group, &
       phi, efield, bfield, &
       domain(:,1), domain(:,3), filter, solver_tolerance=1d-14)

  call propagator%advect_x( delta_t )

  ! Compare to reference
  ! Particle information after advect_x application 
  particle_info_check(:,1) = [5.7840879455647451_f64, 3.4532310684610170_f64, 0.81972274204542517_f64,     -1.5341205443525459_f64,    0.15731068461017067_f64,       1.6241205743525459_f64,   -5.7026946767517002_f64] 
  
  ! Compare computed values to reference values
  do i_part=1,n_particles
     xi = particle_group%group(1)%get_x(i_part)
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
  particle_info_check(:,1) = [5.937500000000000_f64,        3.437500000000000_f64,       0.65731068461017067_f64, -1.5404228405042519_f64, -7.3280069254089197E-002_f64, 1.6241205743525460_f64,   -5.7026946767517002_f64]
  ! Compare computed values to reference values
  do i_part=1,n_particles
     xi = particle_group%group(1)%get_x(i_part)
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

  efield_ref = [  -1.9540083946034448E-002_f64, 9.5300536002179050E-003_f64, &
       1.5165917054294587E-002_f64,  -6.4619318394243501E-003_f64, &
       -4.7673409216992695E-002_f64,   3.5394079124252957E-002_f64, &
       -3.7154800808068197E-002_f64,  0.10950973660536834_f64, &
       -2.6453206184795068E-003_f64,  -4.4548338140878814E-003_f64, &
       3.4878115163240062E-002_f64,  -5.3537442819442477E-002_f64, &
       -1.2276150464793535E-002_f64,   6.0934805407155748E-003_f64, &
       4.1454390064148295E-003_f64,   1.0878935407110090E-002_f64, &
       -6.8363370312354788E-002_f64,   9.2236344695499854E-003_f64, &
       0.13016914008798994_f64,       -6.6215451102733733E-002_f64, &
       -2.2179614454132593E-002_f64,  -8.7928194914971655E-002_f64, &
       0.46719488872691911_f64,      -0.57851981855925927_f64, &
       -5.6322701345351486E-002_f64,   4.8732095061747721E-002_f64, &
       -7.9865941923783171E-002_f64,  0.18240345073847286_f64, &
       3.0475323671229713E-003_f64,  -5.1515868400788326E-003_f64, &
       2.7137457861086334E-002_f64,  -5.7623051404934747E-002_f64, &
       9.3194848614671210E-003_f64,   1.2232238999430283E-002_f64, &
       2.6486270113071197E-002_f64,   1.9387116320619100E-002_f64, &
       -2.5213163085076639E-003_f64,   2.3251868631464032E-002_f64, &
       -4.0907511427739421E-002_f64,   7.2530364168153115E-003_f64, &
       -2.6208080701286759E-003_f64,  -3.5472422969779302E-002_f64, &
       3.6841669313854138E-002_f64,  -1.9025011277121154E-002_f64, &
       -5.2858351597258359E-003_f64,   2.7432014471160983E-003_f64, &
       -2.8816802454405877E-002_f64,  -7.0641784306400002E-003_f64, &
       7.2970889026639160E-003_f64,   2.3603054934407990E-002_f64, &
       0.11471618046417383_f64,        6.8373838212211088E-002_f64, &
       2.4943054252578659E-002_f64,  -7.1873945726164887E-002_f64, &
       0.30878794054271752_f64,        4.6428654899829773E-002_f64, &
       -3.3337235003209648E-002_f64,   7.7531706662584771E-002_f64, &
       -0.47056551205789204_f64,      -0.10469895783464220_f64, &
       5.1828198429107143E-003_f64,  -3.9415053518195069E-002_f64, &
       7.0637794484083416E-002_f64,  -1.2134330496663976E-002_f64, &
       -4.2396426984670386E-003_f64,  -3.8890913031058285E-004_f64, &
       -6.6048050367910935E-002_f64,  -2.6149793757883781E-002_f64, &
       3.1827061439659785E-002_f64,  -7.6590911989053281E-002_f64, &
       4.6465180839071307E-002_f64,  -7.0489947424012309E-002_f64, &
       -1.7834209101075631E-002_f64,   2.8499594490781370E-002_f64, &
       -7.5021546544068965E-002_f64,   4.6647076594328995E-003_f64, &
       3.3442484373463087E-003_f64,  -1.6064829568481938E-002_f64, &
       2.7664256558511739E-002_f64,  -5.9285884489633433E-003_f64, &
       3.2533315334109837E-003_f64,   2.8401829245096235E-003_f64, &
       6.0356611787884804E-002_f64,   2.6640035888829038E-002_f64, &
       1.3541386594736457E-002_f64,  -3.6163046082459918E-002_f64, &
       0.21533022877103658_f64,        4.7939736668972929E-002_f64, &
       -1.6190733224855196E-003_f64,   1.9847542909126145E-002_f64, &
       -3.7232239092948996E-002_f64,   5.0044707571623576E-003_f64, &
       3.3329112177333080E-003_f64,  -5.2987914218084395E-004_f64, &
       1.0865819356987199E-002_f64,   2.6097321955394943E-003_f64 ]
 
  error = maxval(abs(efield-efield_ref))
  if (error> EQV_TOL) then
     passed = .FALSE.
     print*, 'e field wrong in advect_e.', error
  end if

  ! Compare to reference
  ! Particle information after advect_ev application 
  particle_info_check(:,1) = [5.937500000000000_f64,        3.437500000000000_f64,       0.65731068461017067_f64, -1.5218099966480290_f64, 0.14929160117182377_f64, 1.6307170671191553_f64,   -5.7026946767517002_f64]
  ! Compare computed values to reference values
 do i_part=1,n_particles
     xi = particle_group%group(1)%get_x(i_part)
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

  bfield_ref =  [    9.4290430542918351E-004_f64, -2.3205851270612322E-003_f64, &
       5.2215161469197201E-003_f64,  -5.3731922245685727E-004_f64, &
       -8.9686792293844752E-003_f64,   2.5855940005841453E-002_f64, &
       -6.7010475933814401E-002_f64,   1.6649264389435123E-003_f64, &
       1.1938607141322069E-002_f64,  -3.4569034815472161E-002_f64, &
       8.9893366311919909E-002_f64,  -2.0785754607523053E-003_f64, &
       -4.4260307084426310E-003_f64,   1.2543718673234106E-002_f64, &
       -3.2079006906723564E-002_f64,   1.0183428172061406E-003_f64, &
       -1.4458451801888702E-004_f64,   4.2862287230474604E-004_f64, &
       -1.1427046980600000E-003_f64,   1.2620264051124361E-005_f64, &
       4.6823036347203289E-004_f64,  -1.2586181825458703E-003_f64, &
       3.0823147059512723E-003_f64,  -1.6910418772662049E-004_f64, &
       2.3799215224143584E-004_f64,  -8.6533173105837746E-004_f64, &
       2.6003007056110000E-003_f64,   1.1924561679277010E-004_f64, &
       -4.8439506618727819E-005_f64,   1.8528830475733702E-004_f64, &
       -5.6531033180392404E-004_f64,  -3.0136266057762245E-005_f64, &
       2.8916780442458413E-004_f64,  -1.0962004450399635E-003_f64, &
       3.4931221134973049E-003_f64,  -4.4790767711794296E-003_f64, &
       1.9373698029735244E-002_f64,  -3.4969339675477343E-002_f64, &
       8.7656124556165074E-002_f64, -0.10954357904999798_f64, &
       -7.8133853936671917E-003_f64,   1.3682557584269587E-002_f64, &
       -3.3752504628277766E-002_f64,   4.2097299973314038E-002_f64, &
       2.3975400010201168E-003_f64,  -4.0054538928097844E-003_f64, &
       9.6205189897836905E-003_f64,  -1.1958544717152429E-002_f64, &
       -2.9143762670727665E-004_f64,  -8.0606523606953769E-005_f64, &
       5.1265831708396936E-004_f64,   1.6523731315277630E-003_f64, &
       -9.3796967348886012E-003_f64,   6.9476201762240508E-004_f64, &
       4.8751389255749095E-003_f64,   4.1292891931266291E-002_f64, &
       3.6714575509907920E-003_f64,  -3.4616940203048650E-004_f64, &
       -1.6442654540436800E-003_f64,  -1.5894990230555291E-002_f64, &
       -1.0767910774867788E-003_f64,   1.3836488199608444E-004_f64, &
       3.5603514441379363E-004_f64,   4.5283306702353070E-003_f64, &
       1.5002677537248748_f64,        1.5001438809316225_f64, &
       1.4994142063200953_f64,        1.4991762605511667_f64, &
       1.4966285469273173_f64,        1.4996546265819848_f64, &
       1.5033419963344117_f64,        1.5143130144260302_f64, &
       1.5052638764219801_f64,        1.4991031017918039_f64, &
       1.4989060206809901_f64,        1.4766125457442592_f64, &
       1.4980826667318337_f64,        1.5002988960414709_f64, &
       1.5004728505488953_f64,        1.5083197562412636_f64, &
       1.4990114563874557_f64,        1.4994354089621813_f64, &
       1.5022436151154843_f64,       1.5026175257043597_f64, &
       1.5135303960412827_f64,        1.5015575865915620_f64, &
       1.4859351693404355_f64,        1.4415348146139089_f64, &
       1.4789240747793240_f64,        1.5033916123977094_f64, &
       1.5050135645668885_f64,       1.5959504202455388_f64, &
       1.5076157704699107_f64,        1.4988623718506657_f64, &
       1.4979734466251029_f64,        1.4664027663081900_f64 ]

  efield_ref = [   -2.1201798529421539E-002_f64, 9.2605451743335253E-003_f64, &
       1.7205419315076412E-002_f64,  -2.5469087561400330E-004_f64, &
       -3.8014770906463463E-002_f64,   2.6699898816937441E-002_f64, &
       -1.8821365484890010E-002_f64,   6.1340256435701603E-002_f64, &
       -7.8921079123088177E-003_f64,  -7.9754541567077925E-004_f64, &
       2.8071108540148788E-002_f64,  -2.8725198943500160E-002_f64, &
       -9.9356790406780228E-003_f64,   5.0030965790633045E-003_f64, &
       5.6000010810401575E-003_f64,   4.6428786339568857E-004_f64, &
       -6.1315287000178478E-002_f64,   9.6488485074493063E-003_f64, &
       0.12389061998905636_f64,       -9.5866047582595468E-002_f64, &
       -4.0226932059912957E-002_f64,  -7.8809714984830692E-002_f64, &
       0.45355043485543495_f64,      -0.49282790190695602_f64, &
       -4.3574624183127353E-002_f64,   4.4242484441070894E-002_f64, &
       -7.5923775999949011E-002_f64,  0.12382377134412370_f64, &
       -3.7919183589256653E-003_f64,  -3.8088858910071688E-003_f64, &
       2.8097772872175956E-002_f64,  -2.7520049309398581E-002_f64, &
       9.1205004958755335E-003_f64,   1.2997737457131175E-002_f64, &
       2.5050252868022649E-002_f64,   1.9265268276783577E-002_f64, &
       -2.2627565732623572E-003_f64,   1.7437482838967868E-002_f64, &
       -2.7749022029181662E-002_f64,   1.1725723939289424E-002_f64, &
       -1.3580804255452390E-003_f64,  -3.0217445053070034E-002_f64, &
       2.3267527443205156E-002_f64,  -2.7030717817826368E-002_f64, &
       -5.6826656029908933E-003_f64,   7.5923492356027699E-004_f64, &
       -2.3855552412162851E-002_f64,  -4.2768815769034683E-003_f64, &
       7.3930001339984550E-003_f64,   2.2467419923727389E-002_f64, &
       0.11608177420916821_f64,        6.9039319443337646E-002_f64, &
       3.0679615132270438E-002_f64,  -6.7900382952549798E-002_f64, &
       0.30649523163201986_f64,        2.6935889293439541E-002_f64, &
       -4.6891068160524962E-002_f64,   8.0642350391366524E-002_f64, &
       -0.47929662761751285_f64,       -7.0462509996426559E-002_f64, &
       9.8902905540835043E-003_f64,  -4.0298959845784049E-002_f64, &
       7.3406351198300218E-002_f64,  -2.4094201523623252E-002_f64, &
       -8.6047640379350102E-004_f64,  -1.1639782603172709E-002_f64, &
       -3.4876919679699622E-002_f64,  -2.5235607662985320E-002_f64, &
       -1.0778678116354613E-002_f64,  -7.6636020178766015E-003_f64, &
       -6.8320765888057189E-002_f64,  -1.8549387254143793E-002_f64, &
       1.9948187798670466E-003_f64,  -7.4824194001573158E-003_f64, &
       -6.6650487126523328E-003_f64,  -1.6091098174665447E-002_f64, &
       -4.3948732219979718E-003_f64,  -1.2427933862499928E-004_f64, &
       -6.1507572244712582E-003_f64,   5.4749731336118344E-004_f64, &
       3.2752347488351069E-003_f64,   3.2228597341703615E-003_f64, &
       5.9771769930138233E-002_f64,   2.6188074426042127E-002_f64, &
       2.0248572680676355E-002_f64,  -4.0097620061512101E-002_f64, &
       0.22740125129027910_f64,        3.3143778512186972E-002_f64, &
       -4.2110235303194563E-003_f64,   2.1260793103724428E-002_f64, &
       -4.1360120249711148E-002_f64,   1.0776932748559027E-002_f64, &
       4.1621519055910113E-003_f64,  -1.1342104663398418E-003_f64, &
       1.2452021145616079E-002_f64,   9.1728704791739848E-004_f64 ]
 
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
  !deallocate(maxwell_solver)

  call sll_s_halt_collective()

end program test_time_propagator_pic_3d3v_vm_disgradE
