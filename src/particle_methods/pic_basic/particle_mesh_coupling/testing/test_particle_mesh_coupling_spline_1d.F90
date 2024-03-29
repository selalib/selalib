program test_particle_mesh_coupling_spline_1d

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"

  use sll_m_particle_mesh_coupling_base_1d, only: &
    sll_p_collocation, &
    sll_c_particle_mesh_coupling_1d

  use sll_m_particle_mesh_coupling_spline_1d, only: &
    sll_t_particle_mesh_coupling_spline_1d, &
     sll_s_new_particle_mesh_coupling_spline_1d, &
     sll_s_new_particle_mesh_coupling_spline_1d_ptr
    
  use sll_m_particle_group_1d2v, only: &
    sll_t_particle_group_1d2v

  use sll_m_splines_pp

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !class(sll_c_particle_mesh_coupling), allocatable :: ksa
  class(sll_c_particle_mesh_coupling_1d), pointer     :: ksp
  type(sll_t_particle_mesh_coupling_spline_1d) :: kernel
  ! Abstract particle group
  type(sll_t_particle_group_1d2v) :: particle_group

  ! Parameters for the test
  sll_int32 :: n_cells
  sll_int32 :: n_particles
  sll_int32 :: spline_degree
  sll_real64 :: domain(2)
  sll_real64 :: dx
  sll_real64 :: x_vec(4)
  sll_real64 :: v_vec(4,2)

  ! helper variables
  sll_int32  :: i_part 
  sll_real64 :: xi(3), wi(1), vi(3), vi1(3), x_new(1)
  logical :: passed
  sll_real64 :: error
  sll_real64 :: error1
  sll_real64 :: error2
  sll_real64 :: error3
  sll_real64 :: error4
  
  ! Sparse structure for shape factors (for reference)
  !sll_int32 :: index_grid(1,4)    !< First dof index with contribution from th
  sll_real64 :: values_grid(4,1,4) !< Values of the space factors in each dimesion.

  ! Rho dofs
  sll_real64 :: rho_dofs(10)
  sll_real64 :: rho_dofs1(10)
  sll_real64 :: rho_dofs_ref(10)
  sll_real64 :: rho_dofs_pp(4,10) 
  ! J dofs
  sll_real64 :: j_dofs(10)
  sll_real64 :: j_dofs1(10)
  sll_real64 :: j_dofs2(10)
  sll_real64 :: j_dofs3(10,2)
  sll_real64 :: j_dofs4(10)
  sll_real64 :: j_dofs_ref(10)
  ! B dofs
  sll_real64 :: b_dofs(10)

  ! For evaluation check
  sll_real64 :: particle_values(4)
  sll_real64 :: particle_values1(4)
  sll_real64 :: particle_values_ref(4)

  ! 
  passed = .TRUE.

  ! This tests the kernel smoother for a fixed particle and grid and spline degree 3.
  ! Test parameters
  n_cells = 10; ! Number of cells
  n_particles = 4 ! Number of particles
  spline_degree = 3 ! Spline degree
  domain = [0.0_f64, 2.0_f64] ! x_min, x_max
  x_vec = [0.1_f64, 0.65_f64, 0.7_f64, 1.5_f64] ! Particle positions
  v_vec(:,1) = [1.5_f64, -3.0_f64, 0.0_f64, 6.0_f64]
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
  
  values_grid(:,1,1) = [ 2.08333333333333E-002_f64,  0.479166666666666_f64,    &
       0.479166666666666_f64,        2.08333333333333E-002_f64]
  values_grid(:,1,3) = values_grid(:,1,1)   
  values_grid(:,1,4) = values_grid(:,1,1) 
  values_grid(:,1,2) = [7.03125E-002_f64,  0.611979166666666_f64, &
       0.315104166666666_f64,        2.60416666666666E-003_f64 ]

  ! Initialize the kernel
  call kernel%init &
       (domain, n_cells, n_particles, spline_degree, sll_p_collocation)
  
  ! Check that the constructors for the abstract type are working.
  !call sll_s_new_particle_mesh_coupling_spline_1d(ksa, domain, n_cells, n_particles, spline_degree, sll_p_collocation)
  call sll_s_new_particle_mesh_coupling_spline_1d_ptr(ksp, domain, n_cells, n_particles, spline_degree, sll_p_collocation)

  ! Accumulate rho
  rho_dofs = 0.0_f64
  rho_dofs1 = 0.0_f64
  rho_dofs_pp = 0.0_f64
  do i_part = 1,n_particles
     xi = particle_group%get_x(i_part)
     wi = particle_group%get_charge(i_part)
     call kernel%add_charge(xi(1), wi(1), rho_dofs)
     !call kernel%add_charge_pp(xi(1), wi(1), rho_dofs_pp)
  end do
  !rho_dofs = rho_dofs
  rho_dofs_ref = 0.0_f64
  rho_dofs_ref(8:10) = values_grid(1:3,1,1)
  rho_dofs_ref(1) = values_grid(4,1,1)
  rho_dofs_ref(1:4) = rho_dofs_ref(1:4) + values_grid(:,1,2) + values_grid(:,1,3)
  rho_dofs_ref(5:8) = rho_dofs_ref(5:8) + values_grid(:,1,4)
  rho_dofs_ref = rho_dofs_ref/real(n_particles, f64) * real(n_cells,f64)/domain(2)
  error = maxval(abs(rho_dofs-rho_dofs_ref))

  
  if (error > 1.d-14) then
     passed = .FALSE.
     print*, 'Error in procedure add_charge .'
  end if
 
!!$  call sll_s_spline_pp_pp_to_b_1d( kernel%spline_pp, n_cells, rho_dofs_pp, rho_dofs1 )
!!$  error = maxval(abs(rho_dofs1-rho_dofs_ref))
!!$  if (error > 1.d-14) then
!!$     passed = .FALSE.
!!$     print*, 'Error in procedure add_charge_pp .'
!!$  end if

  ! Test j accumulations
  j_dofs = 0.0_f64
  j_dofs1 = 0.0_f64
  j_dofs2 = 0.0_f64
  j_dofs3 = 0.0_f64
  j_dofs4 = 0.0_f64
  b_dofs = 0.0_f64
  do i_part=1, n_particles 
     xi = particle_group%get_x(i_part)
     wi = particle_group%get_charge(i_part)
     vi = particle_group%get_v(i_part)
     vi1 = vi
     x_new = xi(1) + vi(1)/10.0_f64
     
     call kernel%add_current_update_v( xi(1), x_new, wi(1), 1.0_f64, b_dofs, vi, j_dofs )
     !call kernel%add_current_update_v_pp( xi(1), x_new, wi(1), 1.0_f64, b_dofs, vi1, j_dofs1 )

     call kernel%add_current( xi(1), x_new, wi(1), j_dofs2 )
     call kernel%add_current_split( xi(1), x_new, 0, 1, wi(1)*vi(1)/10.0_f64, j_dofs3 )
     call kernel%add_charge_int( xi(1), x_new, wi(1)*vi(1)/10.0_f64, j_dofs4 )
  end do
  j_dofs_ref = [ 2.46175130208333D-002,   4.06901041666666D-005, 0.0_f64, &
       0.0_f64, 0.0_f64, 0.0_f64,  0.0_f64,  6.51041666666666D-004, &
       4.68343098958333D-002,    0.115356445312499_f64];
  j_dofs_ref = j_dofs_ref + [ -0.162190755208333_f64, -0.162190755208333_f64, &
       -2.52685546875D-002,  -4.06901041666666D-005 , &
       0.0_f64, 0.0_f64, 0.0_f64, 0.0_f64,  &
       -4.06901041666666D-005,  -2.52685546875D-002];
  j_dofs_ref = j_dofs_ref + [ 6.51041666666666D-004,  0.0_f64, 0.0_f64,  &
       0.0_f64, 6.51041666666666D-004, 5.01302083333333D-002,  &
       0.199869791666666_f64, 0.248697916666666_f64,  &
       0.199869791666666_f64, 5.01302083333333D-002];
 
  error = maxval(abs(j_dofs-j_dofs_ref))
  !error1= maxval(abs(j_dofs1-j_dofs_ref))
  error2= maxval(abs(j_dofs2-j_dofs_ref*dx))  
  error3= maxval(abs(j_dofs3(:,1)+j_dofs3(:,2)-j_dofs_ref*dx))
  error4= maxval(abs(j_dofs4-j_dofs_ref*dx))
  
  if (error > 1.d-14) then
     passed = .FALSE.
     print*, 'Error in procedure add_current_update_v.'
     
!!$  elseif (error1 > 1.d-14) then
!!$     passed = .FALSE.
!!$     print*, 'Error in procedure add_current_update_v_spline_pp.'
  elseif (error2 > 1.d-14) then
     passed = .FALSE.
     print*, 'Error in procedure add_current.'
  elseif (error3 > 1.d-14) then
     passed = .FALSE.
     print*, 'Error in procedure add_current_split.'
  elseif (error4 > 1.d-14) then
     passed = .FALSE.
     print*, 'Error in procedure add_charge_int.'
  end if

  !call sll_s_spline_pp_b_to_pp_1d(kernel%spline_pp,n_cells,rho_dofs,rho_dofs_pp)
   
  ! Test function evaluation
  do i_part = 1, n_particles
     xi = particle_group%get_x(i_part)
     call kernel%evaluate(xi(1), rho_dofs, particle_values(i_part))
     !call kernel%evaluate_pp(xi(1), rho_dofs_pp, particle_values1(i_part))
  end do
  particle_values_ref = [1.1560058593749998_f64,       2.3149278428819446_f64, &
       2.2656250000000000_f64,        1.1512586805555554_f64]/domain(2);
  error = maxval(abs(particle_values-particle_values_ref))
  error1 = maxval(abs(particle_values1-particle_values_ref))
  !print*,'fehler=', maxval(abs(particle_values1-particle_values))
  if (error > 1.d-14) then
     passed = .FALSE.
     print*, 'Error in procedure evaluate_field_single.'
!!$  else  if (error1 > 1.d-14) then
!!$     passed = .FALSE.
!!$     print*, 'Error in procedure evaluate_field_single_spline_pp.'
  end if
  if (passed .EQV. .TRUE.) then
     print*, 'PASSED'
  else
     print*, 'FAILED'
     stop
  end if

  call kernel%free()
  call particle_group%free()


end program test_particle_mesh_coupling_spline_1d
