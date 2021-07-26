program test_particle_mesh_coupling_spline_2d_feec

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"

  use sll_m_constants, only: &
       sll_p_twopi, &
       sll_p_pi
  
  use sll_m_io_utilities, only : &
    sll_s_read_data_real_array

  use sll_m_particle_mesh_coupling_spline_2d_feec, only: &
    sll_t_particle_mesh_coupling_spline_2d_feec

  use sll_m_particle_group_2d3v, only: &
    sll_t_particle_group_2d3v

  use sll_m_maxwell_2d_fem_fft, only: &
       sll_t_maxwell_2d_fem_fft

  use sll_m_splines_pp
  
  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  sll_int32, parameter :: n_grid = 8 ! number of grid points per dimension
  sll_int32, parameter :: n_particles = 4 ! number of particles
  sll_int32, parameter :: spline_degree = 3 ! degree of the spline
  
  type(sll_t_particle_mesh_coupling_spline_2d_feec) :: particle_mesh
  type(sll_t_particle_group_2d3v) :: particle_group
  type(sll_t_maxwell_2d_fem_fft) :: maxwell

  ! Parameters for the test
  sll_int32 :: n_cells(2)
  sll_real64 :: domain(2,2)
  sll_real64 :: x_vec(3,n_particles)
  sll_real64 :: v_vec(3,n_particles)
  sll_real64 :: v_ref(3,n_particles)

  ! helper variables
  sll_int32  :: i_part 
  sll_real64 :: xi(3), wi(1), vi(3), x_new
  logical :: passed
  sll_real64 :: error

  sll_real64 :: eval, eval_ref

  character(len=256) :: reference_add_charge
  character(len=256) :: reference_add_current
  logical    :: file_exists
  
  ! Sparse structure for shape factors (for reference)
  !sll_int32 :: index_grid(1,4)    !< First dof index with contribution from th
  sll_real64 :: values_grid(4,1,4) !< Values of the space factors in each dimesion.

  ! Rho dofs
  sll_real64 :: rho_dofs(n_grid**2)
  sll_real64 :: rho_dofs_ref(n_grid**2)
  ! J dofs
  sll_real64 :: j_dofs(n_grid**2*2)
  sll_real64 :: j_dofs_ref(n_grid**2*2)
  ! B dofs
  sll_real64 :: b_dofs(n_grid**2*3)
  sll_real64 :: e_dofs(n_grid**2)


  sll_int32  :: n_dofs

  ! Read name of reference file from input argument
  !------------------------------------------------
  call get_command_argument( 1, reference_add_charge )
  call get_command_argument( 2, reference_add_current )

  
  ! Check that file exists    
  !-----------------------
  inquire( file=reference_add_charge, exist=file_exists )
  if (.not. file_exists) then
    write(*,*) &
      "ERROR: reference file '"//trim( reference_add_charge )//"' does not exist"
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
  n_dofs = n_grid**2
  !n_particles = 4 ! Number of particles
  !spline_degree = 3 ! Spline degree
  domain(1,:) = [0.0_f64, 2.0_f64] ! x_min, x_max
  domain(2,:) = [0.0_f64, 2.0_f64] ! x_min, x_max
  x_vec(1,:) = [0.1_f64, 0.65_f64, 0.7_f64, 1.5_f64] ! Particle positions
  x_vec(2,:) = [0.0_f64, 1.0_f64, 0.5_f64, 1.7_f64 ]
  x_vec(3,:) = 0.0_f64![1.0_f64, 0.8_f64, 0.5_f64, 1.8_f64 ]
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
  
  values_grid(:,1,1) = [ 2.0833333333333332E-002_f64,  0.47916666666666663_f64,    &
       0.47916666666666663_f64,        2.0833333333333332E-002_f64]
  values_grid(:,1,3) = values_grid(:,1,1)   
  values_grid(:,1,4) = values_grid(:,1,1) 
  values_grid(:,1,2) = [7.03125E-002_f64,  0.611979166666666_f64, &
       0.315104166666666_f64,        2.60416666666666E-003_f64 ]

  ! Initialize the particle mesh coupler
  call particle_mesh%init( n_cells, domain, n_particles, spline_degree )

  ! Initialize Maxwell solver (for evaluate)
  call maxwell%init( domain, n_cells, spline_degree )
  

  ! Test evaluate
  !print*, 'Evaluate function'
  xi = [0.01_f64, 0.134_f64, 1.756_f64]
  call maxwell%L2projection( test_func, 1, 1, e_dofs )
  call particle_mesh%evaluate( xi, [spline_degree-1, spline_degree],  &
       e_dofs, eval )

  !call sll_s_spline_pp_b_to_pp_2d(particle_mesh%spline_pp_11, n_cells, &
  !     e_dofs(1:n_dofs), e_dofs_pp)
  !call particle_mesh%evaluate_pp( xi, [spline_degree-1, spline_degree],  &
  !     e_dofs_pp, evalpp )
  eval_ref = test_func( xi(1:2) )

  error = eval_ref - eval;
  if (error > 2.1d-3) then
     print*, 'Error in procedure evaluate', error
     passed = .false.
  end if

  !error = eval_ref -evalpp;
  !if (error > 1.1e-3) then
  !   print*, 'Error in procedure evaluate_pp', error
  !   passed = .false.
  !end if
  
  ! Accumulate rho
  rho_dofs = 0.0_f64
  !rho_dofs_pp = 0.0_f64
  do i_part = 1, n_particles
     xi = particle_group%get_x(i_part)
     wi = particle_group%get_charge(i_part)
     call particle_mesh%add_charge(xi, wi(1), &
          [spline_degree, spline_degree, spline_degree], rho_dofs)
    ! call particle_mesh%add_charge_pp(xi, wi(1), &
    !           [spline_degree, spline_degree, spline_degree], rho_dofs_pp)
  end do

  !print*, rho_dofs
  !print*, 'read charge'
  call sll_s_read_data_real_array( reference_add_charge, rho_dofs_ref)
  !write(33,*) rho_dofs

  
  error = maxval(abs(rho_dofs-rho_dofs_ref))

  if (error > 1.d-14) then
     passed = .FALSE.
     print*, 'Error in procedure add_charge.'
  end if

  !error = maxval(abs(rho_dofs_pp-_rho_dofs_ref))

  !if (error > 1.d-14) then
  !   passed = .FALSE.
  !   print*, 'Error in procedure add_charge_pp.'
  !   print*, error
  !end if

  ! add_current_update_v
  ! Set dofs of bfield
  b_dofs = 0.0_f64
  call maxwell%L2projection( beta_cos_k, 1, 2, b_dofs(1:n_dofs) )
  call maxwell%L2projection( beta_cos_k, 2, 2, b_dofs(n_dofs+1:2*n_dofs) )
  call maxwell%L2projection( beta_cos_k, 3, 2, b_dofs(n_dofs*2+1:3*n_dofs) )
  
  j_dofs = 0.0_f64
  do i_part = 1, n_particles
     xi = particle_group%get_x(i_part)
     wi = particle_group%get_charge(i_part)
     vi = particle_group%get_v(i_part)
     x_new = xi(1) + vi(1)/10.0_f64
     call particle_mesh%add_current_update_v_component1( xi, x_new , wi(1), -1.0_f64, b_dofs, vi, &
          j_dofs(1:n_dofs) )
     call particle_group%set_v( i_part, vi )
     
  end do

  !write(32,*) j_dofs(1:n_dofs)


  j_dofs(1+n_dofs:2*n_dofs) = 0.0_f64
  do i_part = 1, n_particles
     xi = particle_group%get_x(i_part)
     wi = particle_group%get_charge(i_part)
     vi = particle_group%get_v(i_part)
     x_new = xi(2) + vi(2)/10.0_f64
     call particle_mesh%add_current_update_v_component2( xi, x_new , wi(1), -1.0_f64, b_dofs, vi, &
          j_dofs(1+n_dofs:2*n_dofs) )
     call particle_group%set_v( i_part, vi )
     !print*, vi
  end do


  !write(32,*) j_dofs
  !print*, 'read current', reference_add_current
  call sll_s_read_data_real_array( reference_add_current, j_dofs_ref)
  error = maxval(abs(j_dofs-j_dofs_ref))
  print*, error

  if (error > 1.d-14) then
     passed = .false.
     print*, 'Error in procedure add_current_update_v.'
  end if

  v_ref(:,1) = [1.49987921014047_f64,        1.2694135490275E-002_f64,  -1.25734555551388E-002_f64]
  v_ref(:,2) = [-2.99773689214078_f64,       0.500_f64,       0.247730406177548_f64]     
  v_ref(:,3) = [0.000_f64,       0.00000_f64,        0.0000_f64]     
  v_ref(:,4) = [6.0000_f64,        4.1677425955085E-002_f64,  -4.1677425955085E-002_f64]

  do i_part = 1, n_particles
     vi = particle_group%get_v( i_part )
     if (maxval(abs(vi-v_ref(:,i_part)))>1E-14_f64) then
        passed = .false.
        print*, 'Error in v in procedure add_current_update_v.'
     end if
  end do

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
      sll_real64, intent(in) :: x(2)

      test_func = cos(sll_p_twopi*x(1)/domain(1,2)) *  cos(sll_p_twopi*x(2)/domain(2,2))
      
    end function test_func

    
    function beta_cos_k(x)
      sll_real64             :: beta_cos_k
      sll_real64, intent(in) :: x(2)
      sll_int32 :: r

      r=1

      beta_cos_k = 0.1_f64 * cos(2*sll_p_pi*x(r)/domain(r,2)) 
    end function beta_cos_k
    
  end program test_particle_mesh_coupling_spline_2d_feec
