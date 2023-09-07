! TODO: Use input from file to initialize and compare
! Unit test for antisymmetric splitting with coordinate transformation
! author: Benedikt Perse, IPP
program test_time_propagator_pic_3d3v_vm_disgradE
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
  particle_info_ref = reshape([5.9375_f64,        3.4375_f64,&
       0.65731068461017_f64,     -1.534120544352545_f64,  &
       0.15731068461017_f64,      1.624120574352545_f64, &
       5.70269467675170_f64], [7, n_particles])

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
  particle_info_check(:,1) = [5.784087945564745_f64, 3.453231068461017_f64, &
       0.819722742045425_f64,     -1.534120544352545_f64,    &
       0.157310684610170_f64,      1.624120574352545_f64,   &
      -5.702694676751700_f64] 
  
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
  particle_info_check(:,1) = [5.9375_f64,        3.4375_f64,    &
       0.657310684610170_f64,     -1.540422840504251_f64, &
      -7.328006925408919E-002_f64, 1.624120574352546_f64,  &
      -5.702694676751700_f64]
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

  efield_ref = [  -1.954008394603444E-002_f64, 9.530053600217905E-003_f64, &
       1.516591705429458E-002_f64,  -6.461931839424350E-003_f64, &
      -4.767340921699269E-002_f64,   3.539407912425295E-002_f64, &
      -3.715480080806819E-002_f64,   0.109509736605368_f64, &
      -2.645320618479506E-003_f64,  -4.454833814087881E-003_f64, &
       3.487811516324006E-002_f64,  -5.353744281944247E-002_f64, &
      -1.227615046479353E-002_f64,   6.093480540715574E-003_f64, &
       4.145439006414829E-003_f64,   1.087893540711009E-002_f64, &
      -6.836337031235478E-002_f64,   9.223634469549985E-003_f64, &
       0.130169140087989_f64,       -6.621545110273373E-002_f64, &
      -2.217961445413259E-002_f64,  -8.792819491497165E-002_f64, &
       0.467194888726919_f64,       -0.578519818559259_f64, &
      -5.632270134535148E-002_f64,   4.873209506174772E-002_f64, &
      -7.986594192378317E-002_f64,   0.182403450738472_f64, &
       3.047532367122971E-003_f64,  -5.151586840078832E-003_f64, &
       2.713745786108633E-002_f64,  -5.762305140493474E-002_f64, &
       9.319484861467121E-003_f64,   1.223223899943028E-002_f64, &
       2.648627011307119E-002_f64,   1.938711632061910E-002_f64, &
      -2.521316308507663E-003_f64,   2.325186863146403E-002_f64, &
      -4.090751142773942E-002_f64,   7.253036416815311E-003_f64, &
      -2.620808070128675E-003_f64,  -3.547242296977930E-002_f64, &
       3.684166931385413E-002_f64,  -1.902501127712115E-002_f64, &
      -5.285835159725835E-003_f64,   2.743201447116098E-003_f64, &
      -2.881680245440587E-002_f64,  -7.064178430640000E-003_f64, &
       7.297088902663916E-003_f64,   2.360305493440799E-002_f64, &
       0.114716180464173_f64,        6.837383821221108E-002_f64, &
       2.494305425257865E-002_f64,  -7.187394572616488E-002_f64, &
       0.308787940542717_f64,        4.642865489982977E-002_f64, &
      -3.333723500320964E-002_f64,   7.753170666258477E-002_f64, &
      -0.470565512057892_f64,       -0.104698957834642_f64, &
       5.182819842910714E-003_f64,  -3.941505351819506E-002_f64, &
       7.063779448408341E-002_f64,  -1.213433049666397E-002_f64, &
      -4.239642698467038E-003_f64,  -3.889091303105828E-004_f64, &
      -6.604805036791093E-002_f64,  -2.614979375788378E-002_f64, &
       3.182706143965978E-002_f64,  -7.659091198905328E-002_f64, &
       4.646518083907130E-002_f64,  -7.048994742401230E-002_f64, &
      -1.783420910107563E-002_f64,   2.849959449078137E-002_f64, &
      -7.502154654406896E-002_f64,   4.664707659432899E-003_f64, &
       3.344248437346308E-003_f64,  -1.606482956848193E-002_f64, &
       2.766425655851173E-002_f64,  -5.928588448963343E-003_f64, &
       3.253331533410983E-003_f64,   2.840182924509623E-003_f64, &
       6.035661178788480E-002_f64,   2.664003588882903E-002_f64, &
       1.354138659473645E-002_f64,  -3.616304608245991E-002_f64, &
       0.215330228771036_f64,        4.793973666897292E-002_f64, &
      -1.619073322485519E-003_f64,   1.984754290912614E-002_f64, &
      -3.723223909294899E-002_f64,   5.004470757162357E-003_f64, &
       3.332911217733308E-003_f64,  -5.298791421808439E-004_f64, &
       1.086581935698719E-002_f64,   2.609732195539494E-003_f64 ]
 
  error = maxval(abs(efield-efield_ref))
  if (error> EQV_TOL) then
     passed = .FALSE.
     print*, 'e field wrong in advect_e.', error
  end if

  ! Compare to reference
  ! Particle information after advect_ev application 
  particle_info_check(:,1) = [5.9375_f64,        3.4375_f64, &
       0.657310684610170_f64, -1.521809996648029_f64, &
       0.149291601171823_f64,  1.630717067119155_f64,  &
      -5.702694676751700_f64]
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

  bfield_ref =  [    9.429043054291835E-004_f64, -2.320585127061232E-003_f64, &
       5.221516146919720E-003_f64,  -5.373192224568572E-004_f64, &
      -8.968679229384475E-003_f64,   2.585594000584145E-002_f64, &
      -6.701047593381440E-002_f64,   1.664926438943512E-003_f64, &
       1.193860714132206E-002_f64,  -3.456903481547216E-002_f64, &
       8.989336631191990E-002_f64,  -2.078575460752305E-003_f64, &
      -4.426030708442631E-003_f64,   1.254371867323410E-002_f64, &
      -3.207900690672356E-002_f64,   1.018342817206140E-003_f64, &
      -1.445845180188870E-004_f64,   4.286228723047460E-004_f64, &
      -1.142704698060000E-003_f64,   1.262026405112436E-005_f64, &
       4.682303634720328E-004_f64,  -1.258618182545870E-003_f64, &
       3.082314705951272E-003_f64,  -1.691041877266204E-004_f64, &
       2.379921522414358E-004_f64,  -8.653317310583774E-004_f64, &
       2.600300705611000E-003_f64,   1.192456167927701E-004_f64, &
      -4.843950661872781E-005_f64,   1.852883047573370E-004_f64, &
      -5.653103318039240E-004_f64,  -3.013626605776224E-005_f64, &
       2.891678044245841E-004_f64,  -1.096200445039963E-003_f64, &
       3.493122113497304E-003_f64,  -4.479076771179429E-003_f64, &
       1.937369802973524E-002_f64,  -3.496933967547734E-002_f64, &
       8.765612455616507E-002_f64,  -0.109543579049997_f64, &
      -7.813385393667191E-003_f64,   1.368255758426958E-002_f64, &
      -3.375250462827776E-002_f64,   4.209729997331403E-002_f64, &
       2.397540001020116E-003_f64,  -4.005453892809784E-003_f64, &
       9.620518989783690E-003_f64,  -1.195854471715242E-002_f64, &
      -2.914376267072766E-004_f64,  -8.060652360695376E-005_f64, &
       5.126583170839693E-004_f64,   1.652373131527763E-003_f64, &
      -9.379696734888601E-003_f64,   6.947620176224050E-004_f64, &
       4.875138925574909E-003_f64,   4.129289193126629E-002_f64, &
       3.671457550990792E-003_f64,  -3.461694020304865E-004_f64, &
      -1.644265454043680E-003_f64,  -1.589499023055529E-002_f64, &
      -1.076791077486778E-003_f64,   1.383648819960844E-004_f64, &
       3.560351444137936E-004_f64,   4.528330670235307E-003_f64, &
       1.500267753724874_f64,        1.500143880931622_f64, &
       1.499414206320095_f64,        1.499176260551166_f64, &
       1.496628546927317_f64,        1.499654626581984_f64, &
       1.503341996334411_f64,        1.514313014426030_f64, &
       1.505263876421980_f64,        1.499103101791803_f64, &
       1.498906020680990_f64,        1.476612545744259_f64, &
       1.498082666731833_f64,        1.500298896041470_f64, &
       1.500472850548895_f64,        1.508319756241263_f64, &
       1.499011456387455_f64,        1.499435408962181_f64, &
       1.502243615115484_f64,        1.502617525704359_f64, &
       1.513530396041282_f64,        1.501557586591562_f64, &
       1.485935169340435_f64,        1.441534814613908_f64, &
       1.478924074779324_f64,        1.503391612397709_f64, &
       1.505013564566888_f64,        1.595950420245538_f64, &
       1.507615770469910_f64,        1.498862371850665_f64, &
       1.497973446625102_f64,        1.466402766308190_f64 ]

  efield_ref = [   -2.120179852942153E-002_f64, 9.260545174333525E-003_f64, &
       1.720541931507641E-002_f64,  -2.546908756140033E-004_f64, &
      -3.801477090646346E-002_f64,   2.669989881693744E-002_f64, &
      -1.882136548489001E-002_f64,   6.134025643570160E-002_f64, &
      -7.892107912308817E-003_f64,  -7.975454156707792E-004_f64, &
       2.807110854014878E-002_f64,  -2.872519894350016E-002_f64, &
      -9.935679040678022E-003_f64,   5.003096579063304E-003_f64, &
       5.600001081040157E-003_f64,   4.642878633956885E-004_f64, &
      -6.131528700017847E-002_f64,   9.648848507449306E-003_f64, &
       0.123890619989056_f64,       -9.586604758259546E-002_f64, &
      -4.022693205991295E-002_f64,  -7.880971498483069E-002_f64, &
       0.453550434855434_f64,       -0.492827901906956_f64, &
      -4.357462418312735E-002_f64,   4.424248444107089E-002_f64, &
      -7.592377599994901E-002_f64,   0.123823771344123_f64, &
      -3.791918358925665E-003_f64,  -3.808885891007168E-003_f64, &
       2.809777287217595E-002_f64,  -2.752004930939858E-002_f64, &
       9.120500495875533E-003_f64,   1.299773745713117E-002_f64, &
       2.505025286802264E-002_f64,   1.926526827678357E-002_f64, &
      -2.262756573262357E-003_f64,   1.743748283896786E-002_f64, &
      -2.774902202918166E-002_f64,   1.172572393928942E-002_f64, &
      -1.358080425545239E-003_f64,  -3.021744505307003E-002_f64, &
       2.326752744320515E-002_f64,  -2.703071781782636E-002_f64, &
      -5.682665602990893E-003_f64,   7.592349235602769E-004_f64, &
      -2.385555241216285E-002_f64,  -4.276881576903468E-003_f64, &
       7.393000133998455E-003_f64,   2.246741992372738E-002_f64, &
       0.116081774209168_f64,        6.903931944333764E-002_f64, &
       3.067961513227043E-002_f64,  -6.790038295254979E-002_f64, &
       0.306495231632019_f64,        2.693588929343954E-002_f64, &
      -4.689106816052496E-002_f64,   8.064235039136652E-002_f64, &
      -0.479296627617512_f64,       -7.046250999642655E-002_f64, &
       9.890290554083504E-003_f64,  -4.029895984578404E-002_f64, &
       7.340635119830021E-002_f64,  -2.409420152362325E-002_f64, &
      -8.60476403793501E-004_f64,   -1.163978260317270E-002_f64, &
      -3.487691967969962E-002_f64,  -2.523560766298532E-002_f64, &
      -1.077867811635461E-002_f64,  -7.663602017876601E-003_f64, &
      -6.832076588805718E-002_f64,  -1.854938725414379E-002_f64, &
       1.994818779867046E-003_f64,  -7.482419400157315E-003_f64, &
      -6.665048712652332E-003_f64,  -1.609109817466544E-002_f64, &
      -4.394873221997971E-003_f64,  -1.242793386249992E-004_f64, &
      -6.150757224471258E-003_f64,   5.474973133611834E-004_f64, &
       3.275234748835106E-003_f64,   3.222859734170361E-003_f64, &
       5.977176993013823E-002_f64,   2.618807442604212E-002_f64, &
       2.024857268067635E-002_f64,  -4.009762006151210E-002_f64, &
       0.227401251290279_f64,        3.314377851218697E-002_f64, &
      -4.211023530319456E-003_f64,   2.126079310372442E-002_f64, &
      -4.136012024971114E-002_f64,   1.077693274855902E-002_f64, &
       4.162151905591011E-003_f64,  -1.134210466339841E-003_f64, &
       1.245202114561607E-002_f64,   9.172870479173984E-004_f64 ]
 
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
