! TODO: Use input from file to initialize and compare
! Unit test for antisymmetric splitting with coordinate transformation
! author: Benedikt Perse, IPP
program test_time_propagator_pic_3d3v_vm_cef_trafo
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
  particle_info_check(:,1) = [5.6450634224675804_f64, 2.6726138917740592_f64, 1.1508242396402408_f64, -1.5341205443525459_f64, 0.15731068461017067_f64, 1.6241205743525460_f64, -5.7026946767517002_f64] 

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

  efield_ref = [  -0.12087219726697715_f64,        6.7835394869837265E-002_f64, &
       0.35635364146897641_f64,      -0.32589641803075198_f64, &
       -0.14723902753431117_f64,        8.8372979269737778E-002_f64, &
       1.2165013057891192E-002_f64,   5.0437830094877872E-002_f64, &
       -2.3910725514331836E-003_f64,  -2.4502494328076725E-002_f64, &
       1.3172880149528727E-003_f64,   2.2074790294999334E-002_f64, &
       -0.11822981870904055_f64,        1.8679467487613712E-002_f64, &
       5.2981151400413033E-002_f64,   5.5692447985697842E-002_f64, &
       4.2132206336932376E-002_f64, -0.32551439232100032_f64, &
       2.0118913666145284_f64,       -2.8322077140992343_f64, &
       -0.54393912771186415_f64,       0.15763721663289462_f64, &
       0.11798949631416607_f64,       0.46499529384000043_f64, &
       5.3344027981364818E-002_f64, -0.13253217260210221_f64, &
       0.14017856964554679_f64,      -0.23690543117257348_f64, &
       -0.38452284466971115_f64,       0.29252068496740713_f64, &
       -0.37666187508489452_f64,       0.90102034961053645_f64, &
       1.1562845965757645E-003_f64,   5.2174400132819244E-002_f64, &
       0.35627112080104223_f64,        1.0896225711680573E-002_f64, &
       2.4888117062141495E-002_f64,   4.5212261457235794E-002_f64, &
       -0.29680385453563707_f64,        4.7287393358505438E-002_f64, &
       -2.8920207033899625E-003_f64, -0.11622705020802837_f64, &
       -0.12600441941301821_f64,      -0.14708354280017963_f64, &
       -2.4926317690431639E-002_f64,   1.8401553499773203E-002_f64, &
       6.7921222906405634E-002_f64,   9.0042686384029472E-002_f64, &
       0.12374266180301813_f64,      -0.36899077381561429_f64, &
       1.9243878603424269_f64,        6.1926024740655325E-002_f64, &
       -9.7225962654755138E-002_f64,  0.38894047524634096_f64, &
       -1.3811005623168371_f64,       0.30542627736102157_f64, &
       2.7944075616544491E-002_f64, -0.28481320368707846_f64, &
       -0.22487978950280529_f64,      -0.53429513953249996_f64, &
       -0.15385421776044653_f64,       0.27692458085194915_f64, &
       -0.29090013587127039_f64,       0.23780267460479643_f64, &
       0.37771844249468867_f64,      -0.63286173392096656_f64, &
       0.73219513271031000_f64,      -0.41266920204875362_f64, &
       -9.0427903102898208E-002_f64,   5.3169873043853177E-002_f64, &
       -0.25821418072442848_f64,       -9.7828653952250696E-002_f64, &
       5.9941623380280126E-002_f64, -0.13425575982401985_f64, &
       0.30545240079068858_f64,      -0.10895729475556595_f64, &
       -0.16245243241422536_f64,       0.25459388263209104_f64, &
       -0.63894715302796312_f64,       0.11723888170162328_f64, &
       5.2824674792626869E-002_f64, -0.12857750721543426_f64, &
       0.71111805399907624_f64,        9.9598911638937693E-002_f64, &
       2.2415107277901136E-002_f64,   5.7851083326353137E-002_f64, &
       9.9598911638937304E-002_f64,  0.20235750901088428_f64, &
       2.6981729918684250E-002_f64,  -2.7640124731776683E-002_f64, &
       5.2824674792626869E-002_f64,   2.2415107277901192E-002_f64, &
       -2.7640124731776738E-002_f64,  0.10517288277746289_f64, &
       -0.12857750721543446_f64,        5.7851083326353303E-002_f64  ]
  
  error = maxval(abs(efield-efield_ref))
  if (error> EQV_TOL) then
     passed = .FALSE.
     print*, 'e field wrong in Hp.', error
  end if

  if (passed .EQV. .FALSE.) then
     print*, 'Error in Hp.'
  end if

  !Reset field info
  efield = [  -0.12087225229535253_f64,6.7835458601481283E-002_f64,&
       0.35635346252897238_f64,-0.32589614839007652_f64,&
       -0.14723895665641937_f64,8.8372891952439206E-002_f64,&
       1.2165160292728587E-002_f64,5.0437671413053026E-002_f64,&
       -2.3911242453951145E-003_f64,-2.4502429882389418E-002_f64,&
       1.3171757510844870E-003_f64,2.2074908256507219E-002_f64,&
       -0.11822975358365719_f64,1.8679399852718478E-002_f64,&
       5.2981300085887258E-002_f64,5.5692243686488230E-002_f64,&
       4.2132261365445456E-002_f64,-0.32551445605197660_f64,&
       2.0118915455538682_f64,-2.8322079837400516_f64,&
       -0.54393919858956541_f64,0.15763730395013806_f64,&
       0.11798934907876493_f64,0.46499545252225338_f64,&
       5.3344079674543923E-002_f64,-0.13253223704802794_f64,&
       0.14017868190991031_f64,-0.23690554913355524_f64,&
       -0.38452290979597525_f64,0.29252075260306482_f64,&
       -0.37666202376992414_f64,0.90102055390941960_f64,&
       1.1561723326829329E-003_f64,5.2174548818247124E-002_f64,&
       0.35627094186099162_f64,1.0896372946513497E-002_f64,&
       2.4888235023586863E-002_f64,4.5212057158028320E-002_f64,&
       -0.29680358489500180_f64,4.7287234676640956E-002_f64,&
       -2.8920723973053118E-003_f64,-0.11622698508264082_f64,&
       -0.12600447444136950_f64,-0.14708347192224111_f64,&
       -2.4926253244704120E-002_f64,1.8401485864918184E-002_f64,&
       6.7921286638111991E-002_f64,9.0042599066728665E-002_f64,&
       0.12374277406735726_f64,-0.36899092250069077_f64,&
       1.9243880392817201_f64,6.1925877505249881E-002_f64,&
       -9.7226080615799321E-002_f64,0.38894067954522693_f64,&
       -1.3811008319576938_f64,0.30542643604323488_f64,&
       2.7944127309770264E-002_f64,-0.28481326881333902_f64,&
       -0.22487973447426779_f64,-0.53429521041015493_f64,&
       -0.15385428220633193_f64,0.27692464848764731_f64,&
       -0.29090019960218433_f64,0.23780276192203745_f64,&
       0.37771822070856659_f64,-0.63286146898226792_f64,&
       0.73219483760592419_f64,-0.41266893711005564_f64,&
       -9.0427736380685811E-002_f64,5.3169669774462580E-002_f64,&
       -0.25821391578572916_f64,-9.7828903345239873E-002_f64,&
       5.9941501105619605E-002_f64,-0.13425559310180726_f64,&
       0.30545217900456578_f64,-0.10895712803335365_f64,&
       -0.16245226569201307_f64,0.25459363323910195_f64,&
       -0.63894688808926448_f64,0.11723867843223326_f64,&
       5.2824674792626813E-002_f64,-0.12857750721543390_f64,&
       0.71111805399907602_f64,9.9598911638938192E-002_f64,&
       2.2415107277901414E-002_f64,5.7851083326352748E-002_f64,&
       9.9598911638937748E-002_f64,0.20235750901088362_f64,&
       2.6981729918684139E-002_f64,-2.7640124731776600E-002_f64,&
       5.2824674792626924E-002_f64,2.2415107277901414E-002_f64,&
       -2.7640124731776544E-002_f64,0.10517288277746234_f64,&
       -0.12857750721543401_f64,5.7851083326352748E-002_f64 ]

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

  efield_ref = [   0.28012636723119727_f64,-0.15091274759792014_f64,&
       0.44588448477647386_f64,-0.33881769067267320_f64,&
       -0.54823757618296953_f64,0.10129443423503599_f64,&
       -7.7365861954773241E-002_f64,0.26918587761245472_f64,&
       8.7139898002106253E-002_f64,-3.7423972164986184E-002_f64,&
       0.40231579527763445_f64,-0.19667329794289423_f64,&
       -0.20776077583115848_f64,0.23742760605211999_f64,&
       -0.34801731944066278_f64,6.8613785969084892E-002_f64,&
       0.44313088089199487_f64,-0.54426266225137765_f64,&
       2.1014225678013689_f64,-2.8451295260226477_f64,&
       -0.94493781811611477_f64,0.17055884623273385_f64,&
       2.8458326831264358E-002_f64,0.68374365872165432_f64,&
       0.14287510192204456_f64,-0.14545377933062376_f64,&
       0.54117730143645926_f64,-0.45565375533295616_f64,&
       -0.47405393204347601_f64,0.51126895880246559_f64,&
       -0.77766064329647344_f64,0.91394209619201561_f64,&
       -0.39984244719386686_f64,0.45317316834479693_f64,&
       0.26673991961349042_f64,0.10042739519401465_f64,&
       0.24363644122298869_f64,3.2290514875431511E-002_f64,&
       -0.28388204261240535_f64,-0.17146097152276010_f64,&
       -9.2423094644806858E-002_f64,-2.6695962835139253E-002_f64,&
       -0.52700309396791922_f64,0.25391514760430878_f64,&
       -1.2004710962108089E-002_f64,-0.20034672033448261_f64,&
       0.28666949283751308_f64,7.7121056784132627E-002_f64,&
       -0.27725584545919202_f64,3.2007697025858393E-002_f64,&
       1.8348570170342191_f64,0.15145689975275081_f64,&
       0.12152212558360126_f64,0.37601913726263153_f64,&
       -1.3681792896750979_f64,8.6678229843834040E-002_f64,&
       -6.1586894937730501E-002_f64,-0.19528224656583812_f64,&
       -0.62587835400081748_f64,-0.13329659088360540_f64,&
       -0.14093273992373592_f64,5.8176442288246290E-002_f64,&
       -7.2151993402783041E-002_f64,0.22488121963944135_f64,&
       0.37771822070856659_f64,-0.63286146898226792_f64,&
       0.73219483760592419_f64,-0.41266893711005564_f64,&
       -9.0427736380685811E-002_f64,5.3169669774462580E-002_f64,&
       -0.25821391578572916_f64,-9.7828903345239873E-002_f64,&
       5.9941501105619605E-002_f64,-0.13425559310180726_f64,&
       0.30545217900456578_f64,-0.10895712803335365_f64,&
       -0.16245226569201307_f64,0.25459363323910195_f64,&
       -0.63894688808926448_f64,0.11723867843223326_f64,&
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
  particle_info_check(:,1) = [5.7984711070298012_f64, 2.6568784534400081_f64, 0.98841218220498617_f64, -1.5258892255233918_f64, 0.15947073762314104_f64, 1.6212074196775912_f64, -5.7026946767517002_f64]
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

  bfield_ref =  [ -0.24058551490716684_f64,0.43921513515233568_f64,&
       -0.86208010976222127_f64,0.20175714530516836_f64,&
       0.21168124596357846_f64,-0.34315817998013221_f64,&
       0.61302295076919988_f64,-0.17756385377924513_f64,&
       -6.6314934935937445E-002_f64,0.10868736189664772_f64,&
       -0.20569138590953831_f64,8.1893637572828359E-002_f64,&
       0.11474311251137864_f64,-0.20724432306090948_f64,&
       0.44952392408559139_f64,-0.12003035515729650_f64,&
       -7.6685994628263761E-003_f64,9.2670617333707860E-003_f64,&
       -2.2254805001658263E-002_f64,-6.4932304132869463E-003_f64,&
       -1.2259036121987327E-002_f64,-5.8257117392746527E-003_f64,&
       2.7748207531516789E-002_f64,1.0524401324540656E-002_f64,&
       4.3405908851021817E-003_f64,4.7922647711196590E-004_f64,&
       -1.0653572680553226E-003_f64,-5.4653870043899549E-003_f64,&
       -3.9368639321412935E-003_f64,-1.4205704791497966E-003_f64,&
       7.9657555516514660E-004_f64,1.5377642151681215E-002_f64,&
       0.34875576585960844_f64,-0.48290185880702535_f64,&
       0.87713013924025596_f64,-0.95920787695638676_f64,&
       -7.6379581600807434E-002_f64,7.1291844861598924E-002_f64,&
       -0.10338859646686918_f64,0.14706556119803565_f64,&
       7.8706492439576961E-002_f64,-9.9284799116098277E-002_f64,&
       0.20365541007431417_f64,-0.21755981429318016_f64,&
       -0.16513500889216204_f64,0.22158663012251512_f64,&
       -0.44334487330250871_f64,0.47153988865318530_f64,&
       -5.1310597470684075E-002_f64,6.1091101274672105E-003_f64,&
       4.7706078808249508E-003_f64,0.25665471012593971_f64,&
       7.3630876934361661E-003_f64,3.2150801984096315E-004_f64,&
       -4.4657064321735190E-003_f64,-4.1808117273061579E-002_f64,&
       -9.3203917276745707E-003_f64,-2.4278042705678238E-004_f64,&
       4.4136185779964478E-003_f64,3.9632264472122182E-002_f64,&
       1.9062148019211789E-002_f64,-1.6430675463735668E-003_f64,&
       -7.5714912259964079E-003_f64,-9.4494225827871453E-002_f64,&
       6.6784481013596810_f64,5.7866429195041045_f64,&
       6.6753193283589631_f64,6.1867377262906729_f64,&
       5.7858007647156198_f64,6.4686065504517920_f64,&
       6.2003541914822220_f64,6.4814183060577797_f64,&
       6.6758715937532624_f64,6.2014070918957103_f64,&
       6.6751808225256610_f64,5.7844743403285550_f64,&
       6.2008753447447766_f64,6.4684627422153786_f64,&
       5.7882455760234688_f64,6.4731195151657301_f64,&
       6.8215443304385515_f64,5.7372672417640285_f64,&
       6.7136788636153790_f64,5.4529167052063086_f64,&
       5.7140202692804358_f64,6.4673151059016192_f64,&
       6.1516789815666364_f64,7.1127915732697344_f64,&
       6.6896265968164998_f64,6.2102583976054859_f64,&
       6.6605113399956561_f64,5.6303936363108571_f64,&
       6.1827392774186487_f64,6.4662307295480295_f64,&
       5.8077815035628246_f64,6.7122103625726854_f64]

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
