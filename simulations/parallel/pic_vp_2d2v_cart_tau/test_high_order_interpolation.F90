program test_high_order_interpolation

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"

  use sll_m_particle_mesh_coupling_base, only: &
    sll_c_particle_mesh_coupling, &
    sll_p_collocation

  use sll_m_particle_mesh_coupling_spline_2d, only: &
    sll_t_particle_mesh_coupling_spline_2d, &
    sll_s_new_particle_mesh_coupling_spline_2d, &
    sll_s_new_particle_mesh_coupling_spline_2d_ptr

  use sll_m_particle_group_2d2v, only: &
    sll_t_particle_group_2d2v

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  class(sll_c_particle_mesh_coupling), pointer     :: ksp
  type(sll_t_particle_mesh_coupling_spline_2d) :: kernel
  ! Abstract particle group
  type(sll_t_particle_group_2d2v) :: particle_group
  ! Parameters for the test
  sll_int32 :: n_cells
  sll_int32 :: n_particles
  sll_int32 :: spline_degree
  sll_real64 :: domain(2,2)
  sll_real64 :: x_vec(4,2)
  sll_real64 :: v_vec(4,2)

  ! helper variables
  sll_int32  :: i_part, i, j, i1, i2
  sll_real64 :: xi(3), wi(1)
  logical :: passed
  sll_real64 :: error

  ! Sparse structure for shape factors (for reference)
  sll_int32 :: index_grid(2,4)    !< First dof index with contribution from th
  sll_real64 :: values_grid(4,2,4) !< Values of the space factors in each dimesion.

  ! Rho dofs
  sll_real64 :: rho_dofs(100)
  sll_real64 :: rho_dofs_ref(100)

  ! For evaluation check
  sll_real64 :: particle_values(4)
  sll_real64 :: particle_values_ref(4)

  sll_real64 :: volume

  !
  passed = .TRUE.

  ! This tests the kernel smoother for a fixed particle and grid and spline degree 3.
  ! Test parameters
  n_cells = 10; ! Number of cells
  n_particles = 4 ! Number of particles
  spline_degree = 3 ! Spline degree
  domain(1,:) = [0.0_f64, 2.0_f64] ! x1_min, x1_max
  domain(2,:) = [0.0_f64, 1.0_f64] ! x2_min, x2_max
  volume = (domain(1,2)-domain(1,1))*(domain(2,2)-domain(2,1))
  x_vec(:,1) = [0.1_f64, 0.65_f64, 0.7_f64, 1.5_f64] ! Particle positions
  x_vec(:,2) = [0.0_f64, 0.0_f64, 0.0_f64, 0.0_f64] ! Particle positions
  v_vec(:,1) = [1.5_f64, 0.0_f64, 0.0_f64, 0.0_f64]
  v_vec(:,2) = [0.0_f64, 0.5_f64, 0.0_f64, 0.0_f64]

  ! We need to initialize the particle group
  call particle_group%init( n_particles, &
       n_particles ,1.0_f64, 1.0_f64, 1)

  call particle_group%set_common_weight(1.0_f64/real(n_particles,f64))

  do i_part = 1,n_particles
     xi(1:2) = x_vec(i_part,:)
     call particle_group%set_x(i_part, xi)
     call particle_group%set_weights(i_part, [1.0_f64])
     xi(1:2) = v_vec(i_part,:)
     call particle_group%set_v(i_part, xi)
  end do

  ! Initialize the kernel
  call kernel%init &
       (domain, [n_cells, n_cells], &
       n_particles, spline_degree, sll_p_collocation)

  ! Check that the constructors for the abstract type are working.
  call sll_s_new_particle_mesh_coupling_spline_2d_ptr(ksp, &
          domain, [n_cells, n_cells], n_particles, &
          spline_degree, sll_p_collocation)

  ! Reference values of the shape factors
  index_grid(1,:) = [-2, 1, 1, 5]
  index_grid(2,:) = [-3, -3, -3, -3]
  values_grid(:,1,1) = [ 2.0833333333333332E-002_f64, &
                         0.47916666666666663_f64,     &
                         0.47916666666666663_f64,     &
                         2.0833333333333332E-002_f64]
  values_grid(:,1,3) = values_grid(:,1,1)
  values_grid(:,1,4) = values_grid(:,1,1)
  values_grid(:,1,2) = [ 7.0312500000000000E-002_f64, &
                         0.6119791666666666_f64,      &
                         0.3151041666666666_f64,      &
                         2.6041666666666665E-003_f64 ]
  values_grid(1,2,:) = 0.0_f64
  values_grid(2,2,:) = 1.0_f64/6.0_f64
  values_grid(3,2,:) = 2.0_f64/3.0_f64
  values_grid(4,2,:) = 1.0_f64/6.0_f64


  ! Accumulate rho
  rho_dofs = 0.0_f64

  do i_part = 1, n_particles
     xi = particle_group%get_x(i_part)
     wi = particle_group%get_charge(i_part)
     call kernel%add_charge(xi(1:2), wi(1), rho_dofs)
  end do
  rho_dofs = rho_dofs

  rho_dofs_ref = 0.0_f64
  rho_dofs_ref(8:10) = values_grid(1:3,1,1)
  rho_dofs_ref(1) = values_grid(4,1,1)
  rho_dofs_ref(1:4) = rho_dofs_ref(1:4) + values_grid(:,1,2) + values_grid(:,1,3)
  rho_dofs_ref(5:8) = rho_dofs_ref(5:8) + values_grid(:,1,4)

  rho_dofs_ref(71:80) = rho_dofs_ref(1:10)/6.0_f64
  rho_dofs_ref(81:90) = rho_dofs_ref(1:10)*2.0_f64/3.0_f64
  rho_dofs_ref(91:100)= rho_dofs_ref(1:10)/6.0_f64
  rho_dofs_ref(1:10) = 0.0_f64

  rho_dofs_ref = rho_dofs_ref  *&
       real(n_cells**2,f64)/volume/real(n_particles, f64)
  error = maxval(abs(rho_dofs-rho_dofs_ref))
  if (error > 1.d-14) then
     passed = .FALSE.
     print*, 'Error in procedure add_charge.'
  end if

  ! Test function evaluation
   do i_part = 1, n_particles
     xi = particle_group%get_x(i_part)
     call  kernel%evaluate(xi(1), rho_dofs, particle_values(i_part))
  end do
  particle_values_ref = [1.1560058593749998_f64, &
                         2.3149278428819446_f64, &
                         2.2656250000000000_f64, &
                         1.1512586805555554_f64]/volume;

  particle_values_ref = 0.0_f64
  do i_part = 1, n_particles
     do i=1,4
        i1 = modulo(index_grid(1,i_part)+i-2, n_cells)
        do j=1,4
           i2 = modulo(index_grid(2,i_part)+j-2, n_cells)
           particle_values_ref(i_part ) = particle_values_ref( i_part) +&
                values_grid(i,1,i_part)*values_grid(j,2,i_part)*&
                rho_dofs_ref(i1+i2*n_cells+1)
        end do
     end do
  end do

  error = maxval(abs(particle_values-particle_values_ref))
  if (error > 1.d-14) then
     passed = .FALSE.
     print*, 'Error in procedure evaluate_field_single.'
  end if

  if (passed .EQV. .TRUE.) then
     print*, 'PASSED'
  else
     print*, 'FAILED'
     stop
  end if

  call kernel%free()
  call particle_group%free()

end program test_high_order_interpolation
