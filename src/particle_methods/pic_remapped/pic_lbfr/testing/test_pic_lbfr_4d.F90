!**************************************************************
!  Copyright INRIA
!  Authors : MCP, ALH
!     CALVI project team
!  
!  This code SeLaLib (for Semi-Lagrangian-Library) 
!  is a parallel library for simulating the plasma turbulence 
!  in a tokamak.
!  
!  This software is governed by the CeCILL-B license 
!  under French law and abiding by the rules of distribution 
!  of free software.  You can  use, modify and redistribute 
!  the software under the terms of the CeCILL-B license as 
!  circulated by CEA, CNRS and INRIA at the following URL
!  "http://www.cecill.info". 
!**************************************************************

!> @file
!> @brief Unit test for initialization and remapping of pic_lbfr_4d particles

! Goal of this test: Start from a particle distribution initialized from a pw affine
! hat function and move it along an affine flow.
! The initialization and remapping algorithms should give exact values on grid points.

! todo: finish the unit test (some functions not written), and add the test in CMakeLists

program test_pic_lbfr_4d

  ! [[file:../working_precision/sll_m_working_precision.h]]

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_pic_lbfr_4d_group, only: &
    SLL_PIC_LBFR_STRUCTURED, &
    SLL_PIC_LBFR_PUSHED, &
    SLL_PIC_LBFR_HAT_F0, &
    SLL_PIC_LBFR_REMAP_WITH_SPARSE_GRIDS, &
    SLL_PIC_LBFR_REMAP_WITH_SPLINES,  &
    SLL_PIC_LBFR_STRUCTURED,  &
    sll_t_pic_lbfr_4d_group, &
    sll_f_pic_lbfr_4d_group_new

  use sll_m_pic_lbfr_4d_utilities, only: &
    sll_f_eval_hat_function

  use sll_m_cartesian_meshes, only: &
    sll_f_new_cartesian_mesh_2d, &
    sll_t_cartesian_mesh_2d

  use sll_m_constants, only: &
    sll_p_pi

  use sll_m_remapped_pic_utilities, only: &
    sll_s_apply_periodic_bc_on_cartesian_mesh_2d, &
    sll_f_x_is_in_domain_2d

  use sll_m_timer, only: &
    sll_s_set_time_mark, &
    sll_f_time_elapsed_since, &
    sll_t_time_mark

  implicit none

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
#define SPECIES_CHARGE  1._f64
#define SPECIES_MASS    1._f64
#define particle_group_id 1_i32

#define NUM_MARKERS_X  10_i32
#define NUM_MARKERS_Y  10_i32
#define NUM_MARKERS_VX 10_i32
#define NUM_MARKERS_VY 10_i32

#define N_UNSTRUCT_MARKERS_PER_CELL 7

#define REMAP_NC_X   20_i32
#define REMAP_NC_Y   20_i32
#define REMAP_NC_VX  20_i32
#define REMAP_NC_VY  20_i32
#define REMAP_DEGREE 1_i32

#define REMAP_SPARSE_GRID_LEVEL_X  6_i32
#define REMAP_SPARSE_GRID_LEVEL_Y  6_i32
#define REMAP_SPARSE_GRID_LEVEL_VX 6_i32
#define REMAP_SPARSE_GRID_LEVEL_VY 6_i32
!#define NUM_DEPOSITION_PARTS    1000_i32

#define DOMAIN_IS_X_PERIODIC .true.
#define DOMAIN_IS_Y_PERIODIC .true. 

#define REMAP_GRID_VX_MIN -6._f64
#define REMAP_GRID_VX_MAX  6._f64
#define REMAP_GRID_VY_MIN -6._f64
#define REMAP_GRID_VY_MAX  6._f64

#define NC_X 256_i32
#define KX   0.5_f64
#define OMEGA   1.323_f64
#define GAMMA  -0.151_f64
#define X_MIN 0._f64
#define X_MAX (2._f64*sll_p_pi/KX)
#define NC_Y 64_i32
#define Y_MIN 0._f64
#define Y_MAX 1._f64

! MCP parameters for 2d plotting of the lt_pic density
#define PLOT_DIM1 1
#define X1_MIN_PLOT X_MIN
#define X1_MAX_PLOT X_MAX
! note that for the plotting grid, n_points = n_cells + 1
#define NC_X1_PLOT 20_i32

#define PLOT_DIM2 3
#define X2_MIN_PLOT REMAP_GRID_VX_MIN
#define X2_MAX_PLOT REMAP_GRID_VX_MAX
! note that for the plotting grid, n_points = n_cells + 1
#define NC_X2_PLOT 20_i32   

#define PLOT_CST_DIM3 2
#define X3_PLOT_CST 0.5_f64

#define PLOT_CST_DIM4 4
#define X4_PLOT_CST 0._f64

  type(sll_t_pic_lbfr_4d_group),        pointer     :: p_group
  type(sll_t_cartesian_mesh_2d),        pointer     :: mesh_2d

  character(5)      :: ncx_name, ncy_name

  sll_int32 :: k
  sll_int32 :: j
  sll_int32 :: j_x, j_y, j_vx, j_vy
  sll_int32 :: n_nodes_x
  sll_int32 :: n_nodes_y
  sll_int32 :: n_nodes_vx
  sll_int32 :: n_nodes_vy
  sll_real64 :: error, tolerance, f_target
  sll_real64 :: f_j
  sll_real64 :: h_x
  sll_real64 :: h_y
  sll_real64 :: h_vx
  sll_real64 :: h_vy
  sll_real64 :: d_vol
  sll_real64 :: x_min
  sll_real64 :: y_min
  sll_real64 :: vx_min
  sll_real64 :: vy_min
  sll_real64 :: x_j
  sll_real64 :: y_j
  sll_real64 :: vx_j
  sll_real64 :: vy_j
  sll_real64 :: new_x
  sll_real64 :: new_y
  sll_real64 :: new_vx
  sll_real64 :: new_vy
  sll_real64 :: x0  
  sll_real64 :: y0  
  sll_real64 :: vx0 
  sll_real64 :: vy0 
  sll_real64 :: r_x 
  sll_real64 :: r_y 
  sll_real64 :: r_vx
  sll_real64 :: r_vy
  sll_real64 :: shift
  sll_real64 :: basis_height

  sll_int32 :: plot_np_x
  sll_int32 :: plot_np_y
  sll_int32 :: plot_np_vx
  sll_int32 :: plot_np_vy

  sll_int32 :: deposition_particles_pos_type
  sll_int32 :: deposition_particles_move_type
  sll_int32 :: n_deposition_particles
  sll_int32 :: n_deposition_particles_per_cell_x
  sll_int32 :: n_deposition_particles_per_cell_y
  sll_int32 :: n_deposition_particles_vx
  sll_int32 :: n_deposition_particles_vy

  sll_int32 :: flow_markers_type

  ! Benchmarking remap performance
  type(sll_t_time_mark) :: remapstart
  sll_real64 :: remaptime

  sll_real64, dimension(3) :: coords
  sll_real64, dimension(4) :: eta
  sll_int32,  dimension(4) :: remapping_sparse_grid_max_levels
  sll_real64 :: x_k,  new_x_k
  sll_real64 :: y_k,  new_y_k
  sll_real64 :: vx_k, new_vx_k
  sll_real64 :: vy_k, new_vy_k

  character(len=10) :: remap_type
  logical :: evaluate_error_on_sparse_grid
  logical :: testing

  sll_real64    :: target_total_charge
  logical       :: enforce_total_charge

  logical :: DEBUG_MODE

  ! --- end of declarations

  DEBUG_MODE = .false.

  plot_np_x  = 21
  plot_np_y  = 10
  plot_np_vx = 21
  plot_np_vy = 10

  open(82,file='error_test_lt_pic_init4D.dat')


  ! A -- build and initialize the particle group

  mesh_2d =>  sll_f_new_cartesian_mesh_2d( NC_X, NC_Y, X_MIN, X_MAX, Y_MIN, Y_MAX )

  if( DEBUG_MODE )then
    print*, "[test_pic_lbfr_4d - DEBUG] -- AA"
  end if

  remapping_sparse_grid_max_levels(1) = REMAP_SPARSE_GRID_LEVEL_X
  remapping_sparse_grid_max_levels(2) = REMAP_SPARSE_GRID_LEVEL_Y
  remapping_sparse_grid_max_levels(3) = REMAP_SPARSE_GRID_LEVEL_VX
  remapping_sparse_grid_max_levels(4) = REMAP_SPARSE_GRID_LEVEL_VY

  deposition_particles_pos_type = SLL_PIC_LBFR_STRUCTURED
  deposition_particles_move_type = SLL_PIC_LBFR_PUSHED
  n_deposition_particles = 1000
  n_deposition_particles_per_cell_x = 2
  n_deposition_particles_per_cell_y = 2
  n_deposition_particles_vx = 50
  n_deposition_particles_vy = 50

  flow_markers_type = SLL_PIC_LBFR_STRUCTURED

  p_group => sll_f_pic_lbfr_4d_group_new(             &
        SPECIES_CHARGE,                               &
        SPECIES_MASS,                                 &
        particle_group_id,                            &
        DOMAIN_IS_X_PERIODIC,                         &
        DOMAIN_IS_Y_PERIODIC,                         &
        SLL_PIC_LBFR_REMAP_WITH_SPARSE_GRIDS,       &       !! other option is SLL_PIC_LBFR_REMAP_WITH_SPLINES
        REMAP_DEGREE,                                 &
        REMAP_GRID_VX_MIN,                            &
        REMAP_GRID_VX_MAX,                            &
        REMAP_GRID_VY_MIN,                            &
        REMAP_GRID_VY_MAX,                            &
        REMAP_NC_X,                                   &
        REMAP_NC_Y,                                   &
        REMAP_NC_VX,                                  &
        REMAP_NC_VY,                                  &
        remapping_sparse_grid_max_levels,             &
        deposition_particles_pos_type,                &
        deposition_particles_move_type,               &
        n_deposition_particles,                  &
        n_deposition_particles_per_cell_x,           &
        n_deposition_particles_per_cell_y,           &
        n_deposition_particles_vx,                   &
        n_deposition_particles_vy,                   &
        flow_markers_type,                            &
        NUM_MARKERS_X,                                &
        NUM_MARKERS_Y,                                &
        NUM_MARKERS_VX,                               &
        NUM_MARKERS_VY,                               &
        N_UNSTRUCT_MARKERS_PER_CELL,                 &
        int(NUM_MARKERS_X/2),                         &
        int(NUM_MARKERS_Y/2),                         &
        int(NUM_MARKERS_VX/2),                        &
        int(NUM_MARKERS_VY/2),                        &
        mesh_2d )

  ! MCP: parameters for f_target (a hat function, with given centers and radius in every dimension)
  x0    =  0.5 * (X_MAX +X_MIN )
  y0    =  0.5 * (Y_MAX +Y_MIN )
  vx0   =  0.5 * ((REMAP_GRID_VX_MAX)+(REMAP_GRID_VX_MIN))
  vy0   =  0.5 * ((REMAP_GRID_VY_MAX)+(REMAP_GRID_VY_MIN))
  r_x   =  0.5 * (X_MAX -X_MIN )
  r_y   =  0.5 * (Y_MAX -Y_MIN )                     
  r_vx  =  0.5 * ((REMAP_GRID_VX_MAX)-(REMAP_GRID_VX_MIN))
  r_vy  =  0.5 * ((REMAP_GRID_VY_MAX)-(REMAP_GRID_VY_MIN))

  !MCP: for a constant (1) function, take shift = 0 and basis_height = 1
  !MCP: for the std tensor-product hat function, take shift = 1 and basis_height = 0
  shift = 1.0_f64
  basis_height = 0.0_f64

  target_total_charge = 0._f64
  enforce_total_charge = .false.

  if( DEBUG_MODE )then
    print*, "[test_pic_lbfr_4d - DEBUG] -- AA B"
  end if

  ! This initializes the particles [[?? file:../pic_particle_initializers/lt_pic_4d_init.F90::sll_lt_pic_4d_init_hat_f]]

  call p_group%set_hat_f0_parameters( x0, y0, vx0, vy0, r_x, r_y, r_vx, r_vy, basis_height, shift )

  if( DEBUG_MODE )then
    print*, "[test_pic_lbfr_4d - DEBUG] -- AA C"
  end if

  call p_group%initializer( SLL_PIC_LBFR_HAT_F0, target_total_charge, enforce_total_charge )


  print *, "plotting initial f slice in gnuplot format "

  call p_group%pic_lbfr_4d_visualize_f_slice_x_vx( "f_initial", plot_np_x, plot_np_y, plot_np_vx, plot_np_vy, 1 )

  if( DEBUG_MODE )then
    print*, "[test_pic_lbfr_4d - DEBUG] -- AA D"
  end if



  ! B --- measure of the error after initialization
  !       (error should be almost 0 on the sparse grid nodes if interpolation is exact)

  open(70,file='f_values_after_init_test4D.dat')

  write(70,*) "# x_j, y_j, vx_j, vy_j, f_j, f_target, abs( f_j - f_target )"

  error = 0.0_f64
  if( p_group%remapped_f_interpolation_type == SLL_PIC_LBFR_REMAP_WITH_SPLINES )then

    n_nodes_x  = p_group%remapping_cart_grid_n_nodes_x()
    n_nodes_y  = p_group%remapping_cart_grid_n_nodes_y()
    n_nodes_vx = p_group%remapping_cart_grid_n_nodes_vx()
    n_nodes_vy = p_group%remapping_cart_grid_n_nodes_vy()

    h_x    = p_group%remapping_cart_grid%delta_eta1
    h_y    = p_group%remapping_cart_grid%delta_eta2
    h_vx   = p_group%remapping_cart_grid%delta_eta3
    h_vy   = p_group%remapping_cart_grid%delta_eta4

    x_min    = p_group%remapping_cart_grid%eta1_min
    y_min    = p_group%remapping_cart_grid%eta2_min
    vx_min   = p_group%remapping_cart_grid%eta3_min
    vy_min   = p_group%remapping_cart_grid%eta4_min

    d_vol = h_x * h_y * h_vx * h_vy


    do j_x = 1, n_nodes_x
      x_j = x_min + (j_x-1) * h_x

      do j_y = 1, n_nodes_y
        y_j = y_min + (j_y-1) * h_y

        do j_vx = 1, n_nodes_vx
          vx_j = vx_min + (j_vx-1) * h_vx

          do j_vy = 1, n_nodes_vy
            vy_j = vy_min + (j_vy-1) * h_vy

            f_j = real( p_group%remapped_f_splines_coefficients(j_x,j_y,j_vx,j_vy)/d_vol ,f64)

            SLL_ASSERT( sll_f_x_is_in_domain_2d( x_j, y_j, p_group%space_mesh_2d,  DOMAIN_IS_X_PERIODIC, DOMAIN_IS_Y_PERIODIC ) )

            f_target = sll_f_eval_hat_function(x0,y0,vx0,vy0,r_x,r_y,r_vx,r_vy, basis_height, shift, x_j, y_j, vx_j, vy_j)
            error = max( error, abs( f_j - f_target ) )
            write(70,*) x_j, y_j, vx_j, vy_j,  f_j, f_target, abs( f_j - f_target )

          end do
        end do
      end do
    end do

  else

    if( DEBUG_MODE )then
      print*, "[test_pic_lbfr_4d - DEBUG] -- AA  12"
    end if

    SLL_ASSERT( p_group%remapped_f_interpolation_type == SLL_PIC_LBFR_REMAP_WITH_SPARSE_GRIDS )

    if( DEBUG_MODE )then
      print*, "[test_pic_lbfr_4d - DEBUG] -- AA  24"
    end if

    evaluate_error_on_sparse_grid = .false.

    if( evaluate_error_on_sparse_grid )then

      do j=1, p_group%sparse_grid_interpolator%size_basis

          f_j = p_group%pic_lbfr_4d_interpolate_value_of_remapped_f( p_group%sparse_grid_interpolator%hierarchy(j)%coordinate )

          x_j  = p_group%sparse_grid_interpolator%hierarchy(j)%coordinate(1)
          y_j  = p_group%sparse_grid_interpolator%hierarchy(j)%coordinate(2)
          vx_j = p_group%sparse_grid_interpolator%hierarchy(j)%coordinate(3)
          vy_j = p_group%sparse_grid_interpolator%hierarchy(j)%coordinate(4)

          SLL_ASSERT( sll_f_x_is_in_domain_2d( x_j, y_j, p_group%space_mesh_2d, DOMAIN_IS_X_PERIODIC, DOMAIN_IS_Y_PERIODIC ) )

          f_target = sll_f_eval_hat_function(x0,y0,vx0,vy0,r_x,r_y,r_vx,r_vy, basis_height, shift, x_j, y_j, vx_j, vy_j)
          error = max( error, abs( f_j - f_target ) )
          write(70,*) x_j, y_j, vx_j, vy_j,  f_j, f_target, abs( f_j - f_target )

      end do

    else

      n_nodes_x  = plot_np_x - 1
      n_nodes_y  = plot_np_y - 1
      n_nodes_vx = plot_np_vx - 1
      n_nodes_vy = plot_np_vy - 1

      x_min    = p_group%remapping_grid_eta_min(1)
      y_min    = p_group%remapping_grid_eta_min(2)
      vx_min   = p_group%remapping_grid_eta_min(3)
      vy_min   = p_group%remapping_grid_eta_min(4)

      h_x    = (p_group%remapping_grid_eta_max(1) - x_min)  / n_nodes_x
      h_y    = (p_group%remapping_grid_eta_max(2) - y_min)  / n_nodes_y
      h_vx   = (p_group%remapping_grid_eta_max(3) - vx_min) / n_nodes_vx
      h_vy   = (p_group%remapping_grid_eta_max(4) - vy_min) / n_nodes_vy

      do j_x = 1, n_nodes_x
        x_j = x_min + (j_x-1) * h_x

        do j_y = 1, n_nodes_y
          y_j = y_min + (j_y-1) * h_y

          do j_vx = 1, n_nodes_vx
            vx_j = vx_min + (j_vx-1) * h_vx

            do j_vy = 1, n_nodes_vy
              vy_j = vy_min + (j_vy-1) * h_vy

              eta(1) = x_j
              eta(2) = y_j
              eta(3) = vx_j
              eta(4) = vy_j

              f_j = p_group%pic_lbfr_4d_interpolate_value_of_remapped_f( eta )

              SLL_ASSERT( sll_f_x_is_in_domain_2d( x_j, y_j, p_group%space_mesh_2d,  DOMAIN_IS_X_PERIODIC, DOMAIN_IS_Y_PERIODIC ) )

              f_target = sll_f_eval_hat_function(x0,y0,vx0,vy0,r_x,r_y,r_vx,r_vy, basis_height, shift, x_j, y_j, vx_j, vy_j)
              error = max( error, abs( f_j - f_target ) )
              write(70,*) x_j, y_j, vx_j, vy_j,  f_j, f_target, abs( f_j - f_target )

            end do
          end do
        end do
      end do

    end if


  end if

  write(82,*) "(maximum) error after initialization: ", error

  close(70)


  if( DEBUG_MODE )then
    print*, "[test_pic_lbfr_4d - DEBUG] -- AA 56"
  end if




  ! C -- push the markers with a measure-preserving affine flow

  ! loop over all particles (taken from [[file:../simulation/simulation_4d_vp_lt_pic_cartesian.F90]])

  do k = 1, p_group%n_particles

     ! -- --  push the k-th particle [begin]  -- --

     ! get particle position
     coords = p_group%get_x(k) ! [[selalib:src/particle_methods/pic_remapped/sll_m_remapped_pic_base.F90::get_v]]
     x_k = coords(1)
     y_k = coords(2)

     ! get particle speed
     coords = p_group%get_v(k) ! [[selalib:src/particle_methods/pic_remapped/sll_m_remapped_pic_base.F90::get_v]]
     vx_k = coords(1)
     vy_k = coords(2)

     ! Flow is cte + A*(x,y,vx,vy), where det A = 1, cf [[test_forward_push]]
     call test_forward_push(x_k, y_k, vx_k, vy_k, new_x_k, new_y_k, new_vx_k, new_vy_k)

     if( .not. sll_f_x_is_in_domain_2d( new_x_k, new_y_k, p_group%space_mesh_2d, DOMAIN_IS_X_PERIODIC, DOMAIN_IS_Y_PERIODIC ) )then
        !! -- -- put outside particles back in domain
        ! [[selalib:src/particle_methods/pic_remapped/remapped_pic_utilities/sll_m_remapped_pic_utilities.F90::apply_periodic_bc_on_cartesian_mesh_2d]]
        call sll_s_apply_periodic_bc_on_cartesian_mesh_2d( p_group%space_mesh_2d, new_x_k, new_y_k)
     end if

     ! set particle position
     coords(1) = new_x_k
     coords(2) = new_y_k
     coords(3) = 0.0_f64
     call p_group%set_x(k, coords)

     ! set particle speed
     coords(1) = new_vx_k
     coords(2) = new_vy_k
     coords(3) = 0.0_f64
     call p_group%set_v(k, coords)

  end do

  testing = .true.
  if( testing )then
    remap_type = 'bsl_ltp'
  else
    ! called as an executable, with an argument
    call get_command_argument(1,remap_type)
  end if

  print *, "plotting transported f slice in gnuplot format "
  call p_group%pic_lbfr_4d_visualize_f_slice_x_vx( "f_transported", plot_np_x, plot_np_y, plot_np_vx, plot_np_vy, 1 )

  ! D --- remap the particle group to get the values of transported f on the remapping grid

  call sll_s_set_time_mark(remapstart)
  if(remap_type == 'ltp') then
    print*, "[test_pic_lbfr_4d]  Error (875454367554242): ltp remapping not implemented yet, stop."
    stop
    print*, "[lt_pic_4d_init_tester]  OLD VERSION: calling sll_lt_pic_4d_write_f_on_remap_grid..."
    ! OLD: call sll_lt_pic_4d_write_f_on_remap_grid( p_group )

  else if (remap_type == 'bsl_ltp') then
    ! remap with [[??? file:lt_pic_4d_utilities.F90::sll_lt_pic_4d_write_bsl_f_on_remap_grid]]
    call p_group%resample( target_total_charge, enforce_total_charge )

  else
    print*, 'ERROR (code=656536756757657): option is ltp (WARNING: not implemented yet) or bsl_ltp'
    stop
  end if

  print *, "plotting remapped f slice in gnuplot format "
  call p_group%pic_lbfr_4d_visualize_f_slice_x_vx( "f_remapped", plot_np_x, plot_np_y, plot_np_vx, plot_np_vy, 1 )

  remaptime = sll_f_time_elapsed_since(remapstart)

  ! formats at [[http://www.cs.mtu.edu/~shene/COURSES/cs201/NOTES/chap05/format.html]]
  write(*,'(A,A,ES8.2,A,A,ES8.2,A)') trim(remap_type),' remap time = ',remaptime,' sec',&
       ' ie ',p_group%n_particles/remaptime,' remapped-ptc/sec'



  ! E --- measure of the error after transport (and remapping):
  !       (result should still be exact on the new grid because all the transformations have been linear)

  open(80,file='f_values_after_remap_test4D.dat')

  write(80,*) "# x_j, y_j, vx_j, vy_j, f_j, f_target, abs( f_j - f_target )"

  error = 0.0_f64
  if( p_group%remapped_f_interpolation_type == SLL_PIC_LBFR_REMAP_WITH_SPLINES )then

    n_nodes_x  = p_group%remapping_cart_grid_n_nodes_x()
    n_nodes_y  = p_group%remapping_cart_grid_n_nodes_y()
    n_nodes_vx = p_group%remapping_cart_grid_n_nodes_vx()
    n_nodes_vy = p_group%remapping_cart_grid_n_nodes_vy()

    h_x    = p_group%remapping_cart_grid%delta_eta1
    h_y    = p_group%remapping_cart_grid%delta_eta2
    h_vx   = p_group%remapping_cart_grid%delta_eta3
    h_vy   = p_group%remapping_cart_grid%delta_eta4

    x_min    = p_group%remapping_cart_grid%eta1_min
    y_min    = p_group%remapping_cart_grid%eta2_min
    vx_min   = p_group%remapping_cart_grid%eta3_min
    vy_min   = p_group%remapping_cart_grid%eta4_min

    d_vol = h_x * h_y * h_vx * h_vy


    do j_x = 1, n_nodes_x
      x_j = x_min + (j_x-1) * h_x

      do j_y = 1, n_nodes_y
        y_j = y_min + (j_y-1) * h_y

        do j_vx = 1, n_nodes_vx
          vx_j = vx_min + (j_vx-1) * h_vx

          do j_vy = 1, n_nodes_vy
            vy_j = vy_min + (j_vy-1) * h_vy

            if( DEBUG_MODE )then
              print*, "[test_pic_lbfr_4d - DEBUG] -- AA -- ddd"
            end if

            f_j = real( p_group%remapped_f_splines_coefficients(j_x,j_y,j_vx,j_vy)/d_vol ,f64)
            call sll_s_test_backward_push(x_j, y_j, vx_j, vy_j, new_x, new_y, new_vx, new_vy)

            if( .not. sll_f_x_is_in_domain_2d( new_x, new_y, p_group%space_mesh_2d,   &
                                         DOMAIN_IS_X_PERIODIC, DOMAIN_IS_Y_PERIODIC ) &
                 ) then
              !! -- -- put outside particles back in domain
              call sll_s_apply_periodic_bc_on_cartesian_mesh_2d( p_group%space_mesh_2d, new_x, new_y)
            end if

            f_target = sll_f_eval_hat_function(x0,y0,vx0,vy0,r_x,r_y,r_vx,r_vy, basis_height, shift, new_x, new_y, new_vx, new_vy)
            error = max( error, abs( f_j - f_target ) )
            write(80,*) x_j, y_j, vx_j, vy_j,  f_j, f_target, abs( f_j - f_target )

          end do
        end do
      end do
    end do

  else

    if( DEBUG_MODE )then
      print*, "[test_pic_lbfr_4d - DEBUG] -- SG -- 1"
    end if

    SLL_ASSERT( p_group%remapped_f_interpolation_type == SLL_PIC_LBFR_REMAP_WITH_SPARSE_GRIDS )

    if( DEBUG_MODE )then
      print*, "[test_pic_lbfr_4d - DEBUG] -- SG -- 2"
    end if

    do j=1, p_group%sparse_grid_interpolator%size_basis

        if( DEBUG_MODE )then
          print*, "[test_pic_lbfr_4d - DEBUG] -- SG -- 3a - j =", j
        end if

        f_j = p_group%pic_lbfr_4d_interpolate_value_of_remapped_f( p_group%sparse_grid_interpolator%hierarchy(j)%coordinate )

        if( DEBUG_MODE )then
          print*, "[test_pic_lbfr_4d - DEBUG] -- SG -- 3b - f_j =", f_j
        end if

        x_j  = p_group%sparse_grid_interpolator%hierarchy(j)%coordinate(1)
        y_j  = p_group%sparse_grid_interpolator%hierarchy(j)%coordinate(2)
        vx_j = p_group%sparse_grid_interpolator%hierarchy(j)%coordinate(3)
        vy_j = p_group%sparse_grid_interpolator%hierarchy(j)%coordinate(4)
        call sll_s_test_backward_push(x_j, y_j, vx_j, vy_j, new_x, new_y, new_vx, new_vy)

        if( .not. sll_f_x_is_in_domain_2d( new_x, new_y, p_group%space_mesh_2d,   &
                                     DOMAIN_IS_X_PERIODIC, DOMAIN_IS_Y_PERIODIC ) &
             ) then
          !! -- -- put outside particles back in domain
          call sll_s_apply_periodic_bc_on_cartesian_mesh_2d( p_group%space_mesh_2d, new_x, new_y)
        end if


        f_target = sll_f_eval_hat_function(x0,y0,vx0,vy0,r_x,r_y,r_vx,r_vy, basis_height, shift, new_x, new_y, new_vx, new_vy)
        error = max( error, abs( f_j - f_target ) )

        write(80,*) x_j, y_j, vx_j, vy_j,  f_j, f_target, abs( f_j - f_target )

        if( DEBUG_MODE )then
          print*, "[test_pic_lbfr_4d - DEBUG] -- SG -- 3c - j, f_target =", j, f_target
        end if

    end do


  end if

  close(80)

  ! uses [[selalib:src/data_structures/fields/sll_m_array_plotting.F90::write_projection_2d]] developed by PN (cf
  ! example in [[selalib:src/data_structures/fields/testing/test_plot_array_4d.F90::write_projection_2d]])

    !  plotting_mesh_2d =>  new_cartesian_mesh_2d( NC_X1_PLOT, NC_X2_PLOT, &
    !         X1_MIN_PLOT, X1_MAX_PLOT, X2_MIN_PLOT, X2_MAX_PLOT )
    !
    !  file_name = "lt_pic_density_test_init4D.dat"
    !
    !  call sll_lt_pic_4d_plot_f_on_2d_grid( &
    !    file_name,        &
    !    PLOT_DIM1,        &
    !    PLOT_DIM2,        &
    !    PLOT_CST_DIM3,    &
    !    PLOT_CST_DIM4,    &
    !    X3_PLOT_CST,      &
    !    X4_PLOT_CST,      &
    !    plotting_mesh_2d, &
    !    p_group    &
    !    )


  write(ncx_name,'(i3)') NC_X
  write(ncy_name,'(i3)') NC_Y
  open(83,file='positions_and_speeds_'//trim(adjustl(ncx_name))//'x'//trim(adjustl(ncy_name))//'.dat')
  do k = 1, p_group%n_particles
     coords = p_group%get_x(k)
     x_k = coords(1)
     y_k = coords(2)

     coords = p_group%get_v(k)
     vx_k = coords(1)
     vy_k = coords(2)

     write(83,*) vx_k, vy_k, x_k, y_k

  enddo
  close(83)

  write(82,*) "(maximum) error after remapping: ", error
  close(82)
  
!  call sll_o_delete( particle_group )
!  call sll_o_delete( mesh_2d )

  tolerance = real( 1e-6, f64)
  if( error < tolerance )then
    ! print*, "TEST PASSED: error = ", error
    print *, 'PASSED'
   else
    print*, "FAILED: error = ", error
  endif

!  print *, 'PASSED'

contains

  ! Only translations on (x,y) for the moment (rotations need non-periodic boundaries)

!#define TRANSLATION_X 1.0
!#define TRANSLATION_Y 1.0
#define TRANSLATION_X 0.0
#define TRANSLATION_Y 0.0

#define COEFF_X_VX 3.0
#define COEFF_Y_VY 0.05555555555555555

!#define COEFF_X_VX 0.0
!#define COEFF_Y_VY 0.0

  ! <<test_forward_push>>
  subroutine test_forward_push(x,y,vx,vy,new_x,new_y,new_vx,new_vy)

    ! [[file:../working_precision/sll_m_working_precision.h]]
    use sll_m_working_precision

    sll_real64,intent(in) :: x
    sll_real64,intent(in) :: y
    sll_real64,intent(in) :: vx
    sll_real64,intent(in) :: vy
    sll_real64,intent(out) :: new_x
    sll_real64,intent(out) :: new_y
    sll_real64,intent(out) :: new_vx
    sll_real64,intent(out) :: new_vy

    new_x = x + TRANSLATION_X + COEFF_X_VX * vx
    new_y = y + TRANSLATION_Y + COEFF_Y_VY * vy
    new_vx = vx
    new_vy = vy

  end subroutine test_forward_push

  subroutine sll_s_test_backward_push(x, y, vx, vy, new_x, new_y, new_vx, new_vy)

    ! [[file:../working_precision/sll_m_working_precision.h]]
    use sll_m_working_precision

    sll_real64,intent(in) :: x
    sll_real64,intent(in) :: y
    sll_real64,intent(in) :: vx
    sll_real64,intent(in) :: vy
    sll_real64,intent(out) :: new_x
    sll_real64,intent(out) :: new_y
    sll_real64,intent(out) :: new_vx
    sll_real64,intent(out) :: new_vy

    new_x = x - TRANSLATION_X - COEFF_X_VX * vx
    new_y = y - TRANSLATION_Y - COEFF_Y_VY * vy
    new_vx = vx
    new_vy = vy

  end subroutine sll_s_test_backward_push

end program test_pic_lbfr_4d
