program test_particle_mesh_coupling_spline_cl_1d

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"

  use sll_m_constants, only: &
       sll_p_twopi
  
  use sll_m_particle_mesh_coupling_base_1d, only: &
    sll_p_collocation, &
    sll_c_particle_mesh_coupling_1d

  use sll_m_particle_mesh_coupling_spline_cl_1d, only: &
    sll_t_particle_mesh_coupling_spline_cl_1d, &
     sll_s_new_particle_mesh_coupling_spline_cl_1d, &
     sll_s_new_particle_mesh_coupling_spline_cl_1d_ptr
    
  use sll_m_particle_group_1d2v, only: &
       sll_t_particle_group_1d2v
  
  use sll_m_splines_pp
  
  use sll_m_maxwell_clamped_1d_fem_sm, only: &
       sll_t_maxwell_clamped_1d_fem_sm

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


  sll_int32, parameter :: n_cells = 10; ! Number of cells
  sll_int32, parameter :: n_particles = 4 ! Number of particles
  sll_int32, parameter :: spline_degree = 3 ! Spline degree
  
  class(sll_c_particle_mesh_coupling_1d), pointer     :: ksp
  type(sll_t_particle_mesh_coupling_spline_cl_1d) :: kernel
  type(sll_t_maxwell_clamped_1d_fem_sm) :: maxwell
  ! Abstract particle group
  type(sll_t_particle_group_1d2v) :: particle_group

  ! Parameters for the test
  sll_int32 :: boundary
  sll_real64 :: domain(2)
  sll_real64 :: dx
  sll_real64 :: x_vec(4)
  sll_real64 :: v_vec(4,2)

  ! helper variables
  sll_int32  :: i_part 
  sll_real64 :: xi(3), wi(1), vi(3), x_new(1), vh
  logical :: passed
  sll_real64 :: error
  sll_real64 :: error2
  sll_real64 :: eval, eval_ref
  
  ! Sparse structure for shape factors (for reference)
  sll_real64 :: values_grid(4,1,4) !< Values of the space factors in each dimesion.

  ! Rho dofs
  sll_real64 :: rho_dofs(n_cells+spline_degree)
  sll_real64 :: rho_dofs_ref(n_cells+spline_degree)
  ! J dofs
  sll_real64 :: j_dofs(n_cells+spline_degree)
  sll_real64 :: j_dofs_ref(n_cells+spline_degree)
  ! 
  passed = .TRUE.

  ! This tests the kernel smoother for a fixed particle and grid and spline degree 3.
  ! Test parameters
  domain = [0.0_f64, 2.0_f64] ! x_min, x_max
  x_vec = [0.1_f64, 0.65_f64, 0.7_f64, 1.5_f64] ! Particle positions
  v_vec(:,1) = [1.5_f64, -3.0_f64, 0.2_f64, 6.0_f64]
  v_vec(:,2) = [0.0_f64, 0.5_f64, 0.0_f64, 0.0_f64]
  dx = (domain(2)-domain(1))/real(n_cells,f64)

  ! We need to initialize the particle group
  call particle_group%init( n_particles, &
       n_particles ,1.0_f64, 1.0_f64, 1)
  
  
  call particle_group%set_common_weight(1.0_f64/real(n_particles,f64))

  do i_part = 1,n_particles
     xi(1) = x_vec(i_part)
     call particle_group%set_x(i_part, xi)
     call particle_group%set_weights(i_part, [1.0_f64])
     xi(1:2) = v_vec(i_part,:)
     call particle_group%set_v(i_part, xi)
  end do
  
  values_grid(:,1,1) = [ 2.0833333333333E-002_f64,  0.4791666666666_f64,    &
       0.47916666666666_f64,        2.0833333333333E-002_f64]
  values_grid(:,1,3) = values_grid(:,1,1)   
  values_grid(:,1,4) = values_grid(:,1,1) 
  values_grid(:,1,2) = [7.03125E-002_f64,  0.61197916666666_f64, &
       0.315104166666666_f64,        2.60416666666666E-003_f64 ]

  if( spline_degree == 2 ) then
     boundary = sll_p_boundary_clamped_square
  else if( spline_degree == 3 ) then
     boundary = sll_p_boundary_clamped_cubic
  else
     boundary = sll_p_boundary_clamped
  end if

  ! Initialize Maxwell solver
  call maxwell%init( domain, n_cells, spline_degree, boundary )

  ! Initialize the kernel
  call kernel%init &
       (domain, n_cells, n_particles, spline_degree, sll_p_collocation, boundary)
  
  ! Check that the constructors for the abstract type are working.
  !call sll_s_new_particle_mesh_coupling_spline_1d(ksa, domain, n_cells, n_particles, spline_degree, sll_p_collocation)
  call sll_s_new_particle_mesh_coupling_spline_cl_1d_ptr(ksp, domain, n_cells, n_particles, spline_degree, sll_p_collocation, boundary)

  ! Accumulate rho
  rho_dofs = 0.0_f64
  do i_part = 1,n_particles
     xi = particle_group%get_x(i_part)
     wi = particle_group%get_charge(i_part)
     call kernel%add_charge(xi(1), wi(1), rho_dofs)
  end do
  
  rho_dofs_ref = [   0.0_f64,       0.7421875_f64, &
       0.3255208333333_f64,       0.1399739583333_f64,&
       1.3639322916666_f64,       0.99283854166666_f64,&
       2.9296874999999E-002_f64,  2.60416666666666E-002_f64, &
       0.5989583333333_f64,       0.59895833333333_f64, &
       2.6041666666666E-002_f64,  0.0_f64,   0.0_f64 ];
  error = maxval(abs(rho_dofs-rho_dofs_ref))

  
  if (error > 1.d-13) then
     passed = .FALSE.
     print*, 'Error in procedure add_charge:', error
  end if

  ! Test j accumulations
  j_dofs = 0.0_f64
  do i_part=1, n_particles 
     xi = particle_group%get_x(i_part)
     wi = particle_group%get_charge(i_part)
     vi = particle_group%get_v(i_part)

     x_new = xi(1) + vi(1)/20.0_f64
     vh = vi(2)/vi(1) 
        
     call kernel%add_current( xi(1), x_new, vh*wi(1), j_dofs )
  end do

  j_dofs_ref = [    0.0_f64,   0.0_f64,  1.08506944444444E-004_f64, &
       7.80571831597222E-003_f64,   1.922607421875E-002_f64, &
       4.10291883680555E-003_f64,   6.78168402777778E-006_f64,   &
       0.0_f64, 0.0_f64, 0.0_f64, 0.0_f64, 0.0_f64, 0.0_f64 ];
  

  error= maxval(abs(j_dofs-j_dofs_ref))  
  if (error > 1.d-14) then
     passed = .FALSE.
     print*, 'Error in procedure add_current:', error2
  end if
   
  ! Test evaluate
  xi = 1.756_f64
  call maxwell%L2projection( test_func, spline_degree, rho_dofs  )
  call kernel%evaluate( xi, rho_dofs, eval )

  eval_ref = test_func( xi(1) )

  error = eval_ref -eval;
  if (error > 1.1d-3) then
     print*, 'Error in procedure evaluate', error
     passed = .false.
  end if

  
  if (passed .EQV. .TRUE.) then
     print*, 'PASSED'
  else
     print*, 'FAILED'
     stop
  end if

  call kernel%free()
  call particle_group%free()

contains
  
  function test_func(x)
    sll_real64             ::test_func
    sll_real64, intent(in) :: x

    test_func = sin(sll_p_twopi*x/domain(2)) 
  end function test_func


end program test_particle_mesh_coupling_spline_cl_1d
