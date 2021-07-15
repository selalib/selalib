program test_particle_mesh_coupling_spline_strong_1d

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"

  use sll_m_particle_mesh_coupling_spline_strong_1d, only: &
    sll_t_particle_mesh_coupling_spline_strong_1d, &
     sll_s_new_particle_mesh_coupling_spline_strong_1d
    
  use sll_m_particle_group_1d2v, only: &
    sll_t_particle_group_1d2v

  use sll_m_splines_pp

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type(sll_t_particle_mesh_coupling_spline_strong_1d) :: kernel
  type(sll_t_particle_mesh_coupling_spline_strong_1d) :: kernel_int
  type(sll_t_particle_mesh_coupling_spline_strong_1d) :: kernel_mid
  type(sll_t_particle_mesh_coupling_spline_strong_1d) :: kernel_mid_int

  ! Abstract particle group
  type(sll_t_particle_group_1d2v) :: particle_group

  ! Parameters for the test
  sll_int32 :: n_cells
  sll_int32 :: n_particles
  sll_int32 :: spline_degree
  sll_real64 :: domain(2)
  sll_real64 :: x_vec(2)
  sll_real64 :: v_vec(2,2)

  ! helper variables
  sll_int32  :: i_part 
  sll_real64 :: xi(3), wi(1)
  logical :: passed
  sll_real64 :: error
  
  ! Sparse structure for shape factors (for reference)
  !sll_int32 :: index_grid(1,4)    !< First dof index with contribution from th
  sll_real64 :: values_grid(4,1,4) !< Values of the space factors in each dimesion.

  ! Rho dofs
  sll_real64 :: rho_dofs(10)
  sll_real64 :: rho_dofs1(10)
  sll_real64 :: rho_dofs_ref(10)
  sll_real64 :: rho_dofs_pp(10) 
  ! J dofs


  ! 
  passed = .TRUE.

  ! This tests the kernel smoother for a fixed particle and grid and spline degree 3.
  ! Test parameters
  n_cells = 10; ! Number of cells
  n_particles = 2 ! Number of particles
  spline_degree = 3 ! Spline degree
  domain = [0.0_f64, 2.0_f64] ! x_min, x_max
  x_vec = [0.05_f64, 1.75_f64] ! Particle positions
  v_vec(:,1) = [1.5_f64, -3.0_f64]
  v_vec(:,2) = [0.0_f64, 0.5_f64]

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

  

  ! Reference for the second order spline (evaluated at the grid points)
  values_grid(1:3,1,1) = [0.0_f64, 0.5_f64, 0.5_f64]
  
  values_grid(:,1,3) = values_grid(:,1,1)   
  values_grid(:,1,4) = values_grid(:,1,1) 
  values_grid(1:3,1,2) = [  3.1250000000000000d-002,       0.68750000000000000_f64,       0.28125000000000000_f64]
  rho_dofs_ref = 0.0_f64
  rho_dofs_ref(1:2) = values_grid(2:3,1,1)
  rho_dofs_ref(3:5) = rho_dofs_ref(3:5) + values_grid(1:3,1,2) + values_grid(1:3,1,3)
  rho_dofs_ref(7:9) = rho_dofs_ref(7:9) + values_grid(1:3,1,4)
  rho_dofs_ref = rho_dofs_ref/real(n_particles, f64) * real(n_cells,f64)/domain(2)

  
  values_grid(:,1,1) = [ 2.0833333333333332d-002,  0.47916666666666663_f64,    &
       0.47916666666666663_f64,        2.0833333333333332d-002]
  
  values_grid(:,1,3) = values_grid(:,1,1)   
  values_grid(:,1,4) = values_grid(:,1,1) 
  values_grid(:,1,2) = [7.0312500000000000d-002,  0.61197916666666663_f64, &
       0.31510416666666663_f64,        2.6041666666666665d-003 ]

  ! Initialize the kernel
  call kernel%init &
       (domain, n_cells, spline_degree, integ=.false., eval_grid_points=.true.)
  call kernel_int%init &
       (domain, n_cells, spline_degree+1, integ=.true., eval_grid_points=.true.)
  call kernel_mid%init &
       (domain, n_cells, spline_degree, integ=.true., eval_grid_points=.false.)
  call kernel_mid_int%init &
       (domain, n_cells, spline_degree+1, integ=.true., eval_grid_points=.false.)
  

  ! Accumulate rho
  rho_dofs = 0.0_f64
  rho_dofs1 = 0.0_f64
  rho_dofs_pp = 0.0_f64
  do i_part = 1,n_particles
     xi = particle_group%get_x(i_part)
     wi = particle_group%get_charge(i_part)
     call kernel%add_charge(xi(1), wi(1), rho_dofs)
     call kernel_int%add_charge(xi(1), wi(1), rho_dofs_pp)
  end do
  rho_dofs_ref = [ 1.7057291666666665_f64,       0.78776041666666652_f64, &
       6.5104166666666661d-003,   0.0000000000000000_f64, &
       0.0000000000000000_f64,        0.0000000000000000_f64, &
       0.0000000000000000_f64,        6.5104166666666661d-003, &
       0.78776041666666652_f64,        1.7057291666666665_f64   ] 

  
  error = maxval(abs(rho_dofs-rho_dofs_ref))

  if (error > 1.d-14) then
     passed = .FALSE.
     print*, 'Error in procedure add_charge (third order).'
  end if
  rho_dofs_ref = [1.6243489583333330_f64,        1.4351399739583330_f64, &
       0.25227864583333326_f64,         4.0690104166666663d-004, &
       0.0000000000000000_f64,         0.0000000000000000_f64, &
       0.0000000000000000_f64,         4.0690104166666663d-004, &
       0.25227864583333326_f64,       1.4351399739583330_f64    ] 
  
  error = maxval(abs(rho_dofs_pp-rho_dofs_ref))

  if (error > 1.d-14) then
     passed = .FALSE.
     print*, 'Error in procedure add_charge (third order integrated).'
  end if
  
   ! Accumulate rho
  rho_dofs = 0.0_f64
  rho_dofs1 = 0.0_f64
  rho_dofs_pp = 0.0_f64
  do i_part = 1,n_particles
     xi = particle_group%get_x(i_part)
     wi = particle_group%get_charge(i_part)
     call kernel_mid%add_charge(xi(1), wi(1), rho_dofs)
     call kernel_mid_int%add_charge(xi(1), wi(1), rho_dofs_pp)
  end do
  rho_dofs_ref = [  1.7057291666666665_f64,  0.78776041666666652_f64, &
       6.5104166666666661d-003,   0.0000000000000000_f64, &
       0.0000000000000000_f64,        0.0000000000000000_f64,        0.0000000000000000_f64, &
       6.5104166666666661d-003,  0.78776041666666652_f64,        1.7057291666666665_f64 ] 

  
  error = maxval(abs(rho_dofs-rho_dofs_ref))
  print*, rho_dofs
  print*, rho_dofs_ref

  if (error > 1.d-14) then
     passed = .FALSE.
     print*, 'Error in procedure add_charge (third order mid point).'
  end if
  rho_dofs_ref = [1.6544596354166663_f64,        0.81258138020833315_f64,  &
       3.2958984375000000d-002,    0.0000000000000000_f64, &
       0.0000000000000000_f64,         0.0000000000000000_f64, &
       0.0000000000000000_f64,         3.2958984375000000d-002, &
       0.81258138020833315_f64,         1.6544596354166663_f64   ] 
  
  error = maxval(abs(rho_dofs_pp-rho_dofs_ref))

  if (error > 1.d-14) then
     passed = .FALSE.
     print*, 'Error in procedure add_charge (third order mid point integrated).'
  end if
 
 
!!$  ! Test j accumulations
!!$  j_dofs = 0.0_f64
!!$  j_dofs1 = 0.0_f64
!!$  b_dofs = 0.0_f64
!!$  do i_part=1, n_particles 
!!$     xi = particle_group%get_x(i_part)
!!$     wi = particle_group%get_charge(i_part)
!!$     vi = particle_group%get_v(i_part)
!!$     vi1 = vi
!!$     x_new = xi(1) + vi(1)/10.0_f64
!!$
!!$     call kernel%add_current_update_v( xi(1), x_new, wi(1), 1.0_f64, b_dofs, vi, j_dofs )
!!$     call kernel_int%add_current_update_v( xi(1), x_new, wi(1), 1.0_f64, b_dofs, vi, j_dofs1 )
!!$     
!!$ end do
!!$ ! For eval_grid_points=.false. (nur erstes Teilchen)
!!$ ! j_dofs_ref = [ 9.9731445312499986E-002_f64,   7.4096679687500000E-002_f64,   3.2958984375000000E-003_f64,   0.0000000000000000_f64,        0.0000000000000000_f64,        0.0000000000000000_f64,        0.0000000000000000_f64,        0.0000000000000000_f64,        0.0000000000000000_f64,        1.0375976562500000E-002_f64];
!!$ j_dofs_ref = [0.24666341145833334_f64,       0.14021809895833331_f64,      -0.13692220052083331_f64, &
!!$      -0.16215006510416663_f64,       -2.5268554687500000E-002_f64,  -4.0690104166685170E-005_f64, &
!!$      6.5104166666666663E-004_f64,    5.0130208333333329E-002_f64,  0.19986979166666666_f64, &
!!$      0.24934895833333331_f64]    
!!$! ! Nur ein Teilchen
!!$! j_dofs_ref = [ 4.6834309895833329E-002_f64,  0.11535644531249999_f64,        2.4617513020833315E-002_f64, &
!!$!      4.0690104166685170E-005_f64,   0.0000000000000000_f64,        0.0000000000000000_f64, &
!!$!      0.0000000000000000_f64,        0.0000000000000000_f64,        0.0000000000000000_f64, &
!!$!      6.5104166666666663E-004_f64] 
!!$ print*, j_dofs, 'a'
!!$ print*, j_dofs_ref
!!$
!!$ error = maxval(abs(j_dofs-j_dofs_ref))
!!$  if (error > 1.e-14) then
!!$     passed = .FALSE.
!!$     print*, 'Error in procedure add_current_update_v (third order).'
!!$
!!$  end if
!!$ 
!!$ j_dofs_ref = [  0.24934895833333334_f64,       0.15299479166666666_f64,      -0.15234374999999997_f64, &
!!$      -0.16992187500000000_f64,       -1.7578125000000028E-002_f64,   0.0000000000000000_f64, &
!!$      0.0000000000000000_f64,        4.1666666666666664E-002_f64,  0.20833333333333334_f64, &
!!$      0.24999999999999997_f64 ]
!!$! Nur erstes Teilchen
!!$ ! j_dofs_ref = [4.1015625000000000E-002_f64,  0.12890625000000000_f64,        1.7578125000000000E-002_f64, &
!!$!      0.0000000000000000_f64,        0.0000000000000000_f64,        0.0000000000000000_f64, &
!!$!      0.0000000000000000_f64,        0.0000000000000000_f64,        0.0000000000000000_f64, &
!!$!      0.0000000000000000_f64 ]
!!$  error = maxval(abs(j_dofs1-j_dofs_ref))
!!$
!!$  
!!$  if (error > 1.e-14) then
!!$     passed = .FALSE.
!!$     print*, 'Error in procedure add_current_update_v (second order).'
!!$
!!$  end if
!!$
!!$  call sll_s_spline_pp_b_to_pp_1d(kernel%spline_pp,n_cells,rho_dofs,rho_dofs_pp)
!!$  
!!$  print*, '#########################'
!!$  ! Test function evaluation
!!$  do i_part = 1, n_particles
!!$     xi = particle_group%get_x(i_part)
!!$     call kernel%evaluate(xi(1), rho_dofs_ref, particle_values(i_part))
!!$    ! print*, particle_values
!!$     call kernel_int%evaluate(xi(1), rho_dofs_ref, particle_values1(i_part))
!!$    ! print*, particle_values1
!!$    ! stop
!!$  end do
!!$  particle_values_ref = [2.8900146484374996_f64,        5.7873196072048598_f64,        5.6640625000000000_f64,        2.8781467013888884_f64];
!!$  error = maxval(abs(particle_values-particle_values_ref))
!!$  !print*,'fehler=', maxval(abs(particle_values1-particle_values))
!!$
!!$  if (error > 1.e-14) then
!!$     passed = .FALSE.
!!$     print*, 'Error in procedure evaluate_field_single (third order).', error
!!$  end if
!!$
!!$  particle_values_ref = [2.9947916666666661_f64,        6.1065673828124982_f64,        5.8919270833333321_f64,        2.9947916666666661_f64]
!!$  error = maxval(abs(particle_values1-particle_values_ref))
!!$  
!!$  if (error > 1.e-14) then
!!$     passed = .FALSE.
!!$     print*, 'Error in procedure evaluate_field_single (second order).', error
!!$  end if
!!$  
!!$
!!$  if (passed .EQV. .TRUE.) then
!!$     print*, 'PASSED'
!!$  else
!!$     print*, 'FAILED'
!!$     stop
!!$  end if

  call kernel%free()
  call kernel_int%free()
  call kernel_mid%free()
  call kernel_mid_int%free()
  !call particle_group%free()
  
  ! Initialize the kernel
  call kernel%init &
       (domain, n_cells, spline_degree-1, integ=.false., eval_grid_points=.true.)
  call kernel_int%init &
       (domain, n_cells, spline_degree, integ=.true., eval_grid_points=.true.)
  call kernel_mid%init &
       (domain, n_cells, spline_degree-1, integ=.false., eval_grid_points=.false.)
  call kernel_mid_int%init &
       (domain, n_cells, spline_degree, integ=.true., eval_grid_points=.false.)
  
  print*, '#####################'
  ! Accumulate rho
  rho_dofs = 0.0_f64
  rho_dofs1 = 0.0_f64
  rho_dofs_pp = 0.0_f64
  do i_part = 1,n_particles
     xi = particle_group%get_x(i_part)
     wi = particle_group%get_charge(i_part)
     call kernel%add_charge(xi(1), wi(1), rho_dofs)
     call kernel_int%add_charge(xi(1), wi(1), rho_dofs_pp)
  end do
  rho_dofs_ref = [ 1.7968750000000000_f64,       0.70312500000000000_f64, &
       0.0000000000000000_f64,        0.0000000000000000_f64, &
       0.0000000000000000_f64,        0.0000000000000000_f64, &
       0.0000000000000000_f64,        0.0000000000000000_f64, &
       0.70312500000000000_f64,        1.7968750000000000_f64   ] 

  
  error = maxval(abs(rho_dofs-rho_dofs_ref))

  if (error > 1.d-14) then
     passed = .FALSE.
     print*, 'Error in procedure add_charge (second order).'
  end if
  rho_dofs_ref = [1.5755208333333330_f64,        1.5364583333333333_f64, &
       0.17578125000000000_f64,        0.0000000000000000_f64, &
       0.0000000000000000_f64,        0.0000000000000000_f64, &
       0.0000000000000000_f64,        0.0000000000000000_f64, &
       0.17578125000000000_f64,        1.5364583333333333_f64  ] 
  
  error = maxval(abs(rho_dofs_pp-rho_dofs_ref))

  if (error > 1.d-14) then
     passed = .FALSE.
     print*, 'Error in procedure add_charge (second order integrated).'
  end if
  
   ! Accumulate rho
  rho_dofs = 0.0_f64
  rho_dofs1 = 0.0_f64
  rho_dofs_pp = 0.0_f64
  do i_part = 1,n_particles
     xi = particle_group%get_x(i_part)
     wi = particle_group%get_charge(i_part)
     call kernel_mid%add_charge(xi(1), wi(1), rho_dofs)
     call kernel_mid_int%add_charge(xi(1), wi(1), rho_dofs_pp)
  end do
  rho_dofs_ref = [ 1.7187500000000000_f64,        7.8125000000000000d-002, &
       0.0000000000000000_f64,        0.0000000000000000_f64, &
       0.0000000000000000_f64,        0.0000000000000000_f64, &
       0.0000000000000000_f64,        7.8125000000000000d-002, &
       1.7187500000000000_f64,        1.4062500000000000_f64    ] 

  
  error = maxval(abs(rho_dofs-rho_dofs_ref))

  if (error > 1.d-14) then
     passed = .FALSE.
     print*, 'Error in procedure add_charge (second order mid point).'
  end if
  rho_dofs_ref = [1.7057291666666665_f64,       0.78776041666666652_f64, &
       6.5104166666666661d-003,  0.0000000000000000_f64, &
       0.0000000000000000_f64,        0.0000000000000000_f64, &
       0.0000000000000000_f64,        6.5104166666666661d-003, &
       0.78776041666666652_f64,        1.7057291666666665_f64   ] 
  
  error = maxval(abs(rho_dofs_pp-rho_dofs_ref))

  if (error > 1.d-14) then
     passed = .FALSE.
     print*, 'Error in procedure add_charge (second order mid point integrated).'
  end if
  
  call kernel%free()
  call kernel_int%free()
  call kernel_mid%free()
  call kernel_mid_int%free()
  call particle_group%free()

  if (passed .eqv. .true. ) then
     print*, 'PASSED.'
  else
     print*, 'FAILED.'
  end if

end program test_particle_mesh_coupling_spline_strong_1d
