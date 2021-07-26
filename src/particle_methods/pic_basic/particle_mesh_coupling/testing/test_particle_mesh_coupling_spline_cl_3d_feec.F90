program test_particle_mesh_coupling_spline_cl_3d_feec

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"

  use sll_m_constants, only: &
       sll_p_twopi, &
       sll_p_pi
  
  use sll_m_io_utilities, only : &
    sll_s_read_data_real_array

  use sll_m_particle_mesh_coupling_spline_cl_3d_feec, only: &
    sll_t_particle_mesh_coupling_spline_cl_3d_feec

  use sll_m_particle_group_3d3v, only: &
    sll_t_particle_group_3d3v

  use sll_m_maxwell_clamped_3d_fem, only: &
       sll_t_maxwell_clamped_3d_fem

  use sll_m_splines_pp
  
  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  sll_int32, parameter :: n_grid = 8 ! number of grid points per dimension
  sll_int32, parameter :: n_particles = 4 ! number of particles
  sll_int32, parameter :: spline_degree = 3 ! degree of the spline
  
  type(sll_t_particle_mesh_coupling_spline_cl_3d_feec) :: particle_mesh
  type(sll_t_particle_group_3d3v) :: particle_group
  type(sll_t_maxwell_clamped_3d_fem) :: maxwell

  ! Parameters for the test
  sll_int32 :: n_cells(3), boundary(3)
  sll_real64 :: domain(3,2)
  sll_real64 :: x_vec(3,n_particles)
  sll_real64 :: v_vec(3,n_particles)
  sll_real64 :: v_ref(3,n_particles)

  ! helper variables
  sll_int32  :: i_part, i, j, ind 
  sll_real64 :: xi(3), wi(1), vi(3), x_new(3), vh(3)
  logical :: passed
  sll_real64 :: error

  sll_real64 :: eval, eval_ref

  character(len=256) :: reference_add_charge
  character(len=256) :: reference_add_particle_mass
  character(len=256) :: reference_add_current
  logical    :: file_exists
  
  ! Sparse structure for shape factors (for reference)
  sll_real64 :: values_grid(4,1,4) !< Values of the space factors in each dimesion.

  ! Rho dofs
  sll_real64 :: rho_dofs((n_grid+spline_degree)*n_grid**2)
  sll_real64 :: rho_dofs_ref((n_grid+spline_degree)*n_grid**2)
  ! J dofs
  sll_real64 :: j_dofs((n_grid+spline_degree)*n_grid**2*2+(n_grid+spline_degree-1)*n_grid**2)
  sll_real64 :: j_dofs_ref((n_grid+spline_degree)*n_grid**2*2+(n_grid+spline_degree-1)*n_grid**2)
  ! B dofs
  sll_real64 :: b_dofs((n_grid+spline_degree-1)*n_grid**2*2+(n_grid+spline_degree)*n_grid**2)
  sll_real64 :: e_dofs((n_grid+spline_degree-1)*n_grid**2)
  ! Particle mass
  sll_real64 :: particle_mass((2*spline_degree+1)**2*(2*spline_degree-1),(n_grid+spline_degree)*n_grid**2)
  sll_real64 :: particle_mass_ref((2*spline_degree+1)**2*(2*spline_degree-1)*(n_grid+spline_degree)*n_grid**2)


  sll_int32  :: n_dofs0, n_dofs1

  ! Read name of reference file from input argument
  !------------------------------------------------
  call get_command_argument( 1, reference_add_charge )
  call get_command_argument( 2, reference_add_particle_mass )
  call get_command_argument( 3, reference_add_current )

  
  ! Check that file exists    
  !-----------------------
  inquire( file=reference_add_charge, exist=file_exists )
  if (.not. file_exists) then
    write(*,*) &
      "ERROR: reference file '"//trim( reference_add_charge )//"' does not exist"
    stop
 end if
 !-----------------------
  inquire( file=reference_add_particle_mass, exist=file_exists )
  if (.not. file_exists) then
    write(*,*) &
      "ERROR: reference file '"//trim( reference_add_particle_mass )//"' does not exist"
    stop
 end if
 !-----------------------
  inquire( file=reference_add_current, exist=file_exists )
  if (.not. file_exists) then
    write(*,*) &
      "ERROR: reference file '"//trim( reference_add_current )//"' does not exist"
    stop
 end if
  ! 
 passed = .TRUE.

 
  ! This tests the kernel smoother for a fixed particle and grid and spline degree 3.
  ! Test parameters
  n_cells = n_grid; ! Number of cells
  n_dofs0 = (n_grid+spline_degree)*n_grid**2
  n_dofs1 = (n_grid+spline_degree-1)*n_grid**2
  domain(1,:) = [0.0_f64, 2.0_f64] ! x_min, x_max
  domain(2,:) = [0.0_f64, 2.0_f64] ! x_min, x_max
  domain(3,:) = [0.0_f64, 2.0_f64] ! x_min, x_max
  x_vec(1,:) = [0.1_f64, 0.65_f64, 0.7_f64, 1.5_f64] ! Particle positions
  x_vec(2,:) = [0.0_f64, 1.0_f64, 0.5_f64, 1.7_f64 ]
  x_vec(3,:) = [1.0_f64, 0.8_f64, 0.5_f64, 1.8_f64 ]
  v_vec(1,:) = [1.5_f64, -3.0_f64, 0.0_f64, 6.0_f64]
  v_vec(2,:) = [0.0_f64, 0.5_f64, 0.0_f64, 0.0_f64]
  v_vec(3,:) = [0.0_f64, 0.25_f64, 0.0_f64, 0.0_f64]

  
  ! We need to initialize the particle group
  call particle_group%init( n_particles, &
       n_particles ,1.0_f64, 1.0_f64, 1)
  
  
  call particle_group%set_common_weight(1.0_f64/real(n_particles,f64))

  xi = 0.0_f64
  vi = 0.0_f64
  do i_part = 1,n_particles
     xi = x_vec(:,i_part)
     call particle_group%set_x(i_part, xi)
     call particle_group%set_weights(i_part, [1.0_f64])
     xi = v_vec(:,i_part)
     call particle_group%set_v(i_part, xi)
  end do
  
  values_grid(:,1,1) = [ 2.083333333333333E-002_f64,  0.479166666666666_f64,    &
       0.479166666666666_f64,        2.08333333333333E-002_f64]
  values_grid(:,1,3) = values_grid(:,1,1)   
  values_grid(:,1,4) = values_grid(:,1,1) 
  values_grid(:,1,2) = [7.03125E-002_f64,  0.611979166666666_f64, &
       0.315104166666666_f64,        2.60416666666666E-003_f64 ]


  if( spline_degree == 2 ) then
     boundary = [ sll_p_boundary_clamped_square, sll_p_boundary_periodic, sll_p_boundary_periodic]
  else if( spline_degree == 3 ) then
     boundary = [ sll_p_boundary_clamped_cubic, sll_p_boundary_periodic, sll_p_boundary_periodic]
  else
     boundary = [ sll_p_boundary_clamped, sll_p_boundary_periodic, sll_p_boundary_periodic]
  end if

  ! Initialize the particle mesh coupler
  call particle_mesh%init( n_cells, domain, [spline_degree,spline_degree,spline_degree], boundary, n_particles )

  ! Initialize Maxwell solver (for evaluate)
  call maxwell%init( domain, n_cells, [spline_degree,spline_degree,spline_degree], boundary )
  

  ! Test evaluate
  print*, 'Evaluate function'
  xi = [0.01_f64, 0.134_f64, 1.756_f64]
  call maxwell%L2projection( 1, 1, e_dofs, test_func )
  call particle_mesh%evaluate( xi, [spline_degree-1, spline_degree, spline_degree],  &
       e_dofs, eval )

  eval_ref = test_func( xi )

  error = eval_ref -eval;
  if (error > 1.1d-3) then
     print*, 'Error in procedure evaluate', error
     passed = .false.
  end if
  
  ! Accumulate rho
  rho_dofs = 0.0_f64
  do i_part = 1, n_particles
     xi = particle_group%get_x(i_part)
     wi = particle_group%get_charge(i_part)
     call particle_mesh%add_charge(xi, wi(1), &
          [spline_degree, spline_degree, spline_degree], rho_dofs)
  end do

  call sll_s_read_data_real_array( reference_add_charge, rho_dofs_ref)
  error = maxval(abs(rho_dofs-rho_dofs_ref))

  if (error > 1.d-14) then
     passed = .FALSE.
     print*, 'Error in procedure add_charge.'
  end if

  ! add_particle_mass
  particle_mass = 0.0_f64
  do i_part = 1, n_particles
     xi = particle_group%get_x(i_part)
     wi = particle_group%get_charge(i_part)
     call particle_mesh%add_particle_mass(xi, wi(1), &
          [spline_degree, spline_degree-1, spline_degree], particle_mass)
  end do

  
  call sll_s_read_data_real_array( reference_add_particle_mass, particle_mass_ref)

  error = 0.0_f64
  ind=1
  do i = 1, (n_grid+spline_degree)*n_grid**2
     do j = 1, (2*spline_degree+1)**2*(2*spline_degree-1)
        error = max( error,abs( particle_mass(j,i)-particle_mass_ref(ind) ))
        ind = ind+1
     end do
  end do

  if (error > 1.d-14) then
     passed = .FALSE.
     print*, 'Error in procedure add_particle_mass.'
  end if

  ! add_current

  
  j_dofs = 0.0_f64
  do i_part = 1, n_particles
     xi = particle_group%get_x(i_part)
     wi = particle_group%get_charge(i_part)
     vi = particle_group%get_v(i_part)
     x_new = xi + vi/20.0_f64
     vh = (x_new-xi)*wi(1)
     call particle_mesh%add_current( xi, x_new , vh, j_dofs )
  end do

  call sll_s_read_data_real_array( reference_add_current, j_dofs_ref)
  error = maxval(abs(j_dofs-j_dofs_ref))
  print*, error

  if (error > 1.d-14) then
     passed = .false.
     print*, 'Error in procedure add_current.'
  end if


  call particle_mesh%free()
  call particle_group%free()

  if ( passed .eqv. .false.) then
     print*, 'FAILED.'
  else
     print*, 'PASSED.'
  end if

  
contains
  
  function test_func(x)
      sll_real64             ::test_func
      sll_real64, intent(in) :: x(3)

      test_func = cos(sll_p_twopi*x(1)/domain(1,2)) *  cos(sll_p_twopi*x(2)/domain(2,2)) *  &
           cos(sll_p_twopi*x(3)/domain(3,2))
    end function test_func

    
    function beta_cos_k(x)
      sll_real64             :: beta_cos_k
      sll_real64, intent(in) :: x(3)
      sll_int32 :: r

      r=1

      beta_cos_k = 0.1_f64 * cos(2._f64*sll_p_pi*x(r)/domain(r,2)) 
    end function beta_cos_k
    
end program test_particle_mesh_coupling_spline_cl_3d_feec
