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
!> @brief Unit test for initialization and remapping of bsl_lt_pic_4d particles

! Goal of this test: Start from a particle distribution initialized from a pw affine
! hat function and move it along an affine flow.
! The initialization and remapping algorithms should give exact values on grid points.

! todo: finish the unit test (some functions not written), and add the test in CMakeLists

program unit_test_bsl_lt_pic_4d

  ! [[file:../working_precision/sll_working_precision.h]]
#include "sll_working_precision.h"

  ! [[file:../memory/sll_memory.h]]
#include "sll_memory.h"

  ! [[file:../assert/sll_assert.h]]
#include "sll_assert.h"

  use sll_constants, only: sll_pi
  use sll_bsl_lt_pic_4d_group_module
  use sll_cartesian_meshes
!  use sll_representation_conversion_module
  use sll_timer

#define SPECIES_CHARGE  1._f64
#define SPECIES_MASS    1._f64
#define particle_group_id 1_i32

#define NUM_PARTS_X 20_i32
#define NUM_PARTS_Y 20_i32
#define NUM_PARTS_VX 21_i32
#define NUM_PARTS_VY 21_i32
#define SPLINE_DEGREE 1_i32

#define DOMAIN_IS_X_PERIODIC .true.
#define DOMAIN_IS_Y_PERIODIC .true. 

#define REMAP_GRID_VX_MIN -6._f64
#define REMAP_GRID_VX_MAX  6._f64
#define REMAP_GRID_VY_MIN -6._f64
#define REMAP_GRID_VY_MAX  6._f64

#define N_virtual_particles_per_deposition_cell_x  2_i32
#define N_virtual_particles_per_deposition_cell_y  2_i32
#define N_virtual_particles_per_deposition_cell_vx 2_i32
#define N_virtual_particles_per_deposition_cell_vy 2_i32

#define NC_X 256_i32
#define KX   0.5_f64
#define OMEGA   1.323_f64
#define GAMMA  -0.151_f64
#define X_MIN 0._f64
#define X_MAX (2._f64*sll_pi/KX)
#define NC_Y 64_i32
#define Y_MIN 0._f64
#define Y_MAX 1._f64
!#define SIMDT 0.01_f64
!#define NUMBER_ITERATIONS 1600_i32

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


  implicit none
  
  type(sll_bsl_lt_pic_4d_group),      pointer     :: part_group
  !  type(sll_lt_pic_4d_particle), pointer :: particle


  !  type(sll_lt_pic_4d_group),  pointer :: part_group
  type(sll_cartesian_mesh_2d),       pointer :: mesh_2d
  type(sll_cartesian_mesh_2d),       pointer :: plotting_mesh_2d

  character(len=30) :: file_name
  character(5)      :: ncx_name, ncy_name

  sll_int64 :: i
  sll_int64 :: j, j_x, j_y, j_vx, j_vy
!  sll_int32 :: part_array_size, part_guard_size
  sll_real64 :: x, y, error, tolerance, f_target
  sll_real64 :: f_j
  sll_int64 :: number_nodes_x
  sll_int64 :: number_nodes_y
  sll_int64 :: number_nodes_vx
  sll_int64 :: number_nodes_vy
  sll_real64 :: h_nodes_x    
  sll_real64 :: h_nodes_y    
  sll_real64 :: h_nodes_vx   
  sll_real64 :: h_nodes_vy   
  sll_real64 :: nodes_x_min  
  sll_real64 :: nodes_y_min  
  sll_real64 :: nodes_vx_min 
  sll_real64 :: nodes_vy_min 
  sll_real64 :: x_j
  sll_real64 :: y_j
  sll_real64 :: vx_j
  sll_real64 :: vy_j
  sll_real64 :: x0  
  sll_real64 :: y0  
  sll_real64 :: vx0 
  sll_real64 :: vy0 
  sll_real64 :: r_x 
  sll_real64 :: r_y 
  sll_real64 :: r_vx
  sll_real64 :: r_vy
  sll_real64 :: hat_shift
  sll_real64 :: basis_height


  sll_real64 :: x_j_bsl_before_per
  sll_real64 :: y_j_bsl_before_per
  sll_real64 :: vx_j_bsl_before_per
  sll_real64 :: vy_j_bsl_before_per

  sll_real64 :: x_j_bsl
  sll_real64 :: y_j_bsl
  sll_real64 :: vx_j_bsl
  sll_real64 :: vy_j_bsl

  sll_real64 :: x_j_bsl_end
  sll_real64 :: y_j_bsl_end
  sll_real64 :: vx_j_bsl_end
  sll_real64 :: vy_j_bsl_end

  ! Benchmarking remap performance
  type(sll_time_mark) :: remapstart
  sll_real64 :: remaptime

!  sll_real64 :: max_f
!  sll_real64 :: inv_r_x  
!  sll_real64 :: inv_r_y  
!  sll_real64 :: inv_r_vx 
!  sll_real64 :: inv_r_vy 
!  sll_real32 :: d_vol


  sll_real64 :: new_x,new_y,new_vx,new_vy

  character(len=10) :: remap_type

  ! --- end of declarations


  mesh_2d =>  new_cartesian_mesh_2d( NC_X, NC_Y, X_MIN, X_MAX, Y_MIN, Y_MAX )

!  part_group => new_particle_2d_group( &
!       NUM_PARTICLES, &
!       PARTICLE_ARRAY_SIZE, &
!       GUARD_SIZE, QoverM, mesh_2d )

!  part_array_size = NUM_PARTS_X * NUM_PARTS_Y * NUM_PARTS_VX * NUM_PARTS_VY
!  part_guard_size = part_array_size / 10

  part_group => sll_bsl_lt_pic_4d_group_new( &
        SPECIES_CHARGE,                                             &
        SPECIES_MASS,                                               &
        particle_group_id,                                          &
        DOMAIN_IS_X_PERIODIC,                                       &
        DOMAIN_IS_Y_PERIODIC,                                       &
        SPLINE_DEGREE,          &
        NUM_PARTS_X,            &
        NUM_PARTS_Y,            &
        NUM_PARTS_VX,           &
        NUM_PARTS_VY,           &
        REMAP_GRID_VX_MIN,      &
        REMAP_GRID_VX_MAX,      &
        REMAP_GRID_VY_MIN,      &
        REMAP_GRID_VY_MAX,      &
        N_virtual_particles_per_deposition_cell_x,   &
        N_virtual_particles_per_deposition_cell_y,   &
        N_virtual_particles_per_deposition_cell_vx,  &
        N_virtual_particles_per_deposition_cell_vy,  &
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

  !MCP: for a constant (1) function, take hat_shift = 0 and basis_height = 1
  !MCP: for the std tensor-product hat function, take hat_shift = 1 and basis_height = 0
  hat_shift = 1.
  basis_height = 0.

  ! This initializes the particles [[?? file:../pic_particle_initializers/lt_pic_4d_init.F90::sll_lt_pic_4d_init_hat_f]]

  call part_group%set_hat_f0_parameters( x0,y0,vx0,vy0,r_x,r_y,r_vx,r_vy, basis_height, hat_shift )   ! todo

  initial_density_identifier = 1  ! 1: initialize with a hat function
  call sim%particle_group%initializer( initial_density_identifier ) !todo with initial_density_identifier = 1


    !  call sll_lt_pic_4d_init_hat_f (                                           &
    !             x0,y0,vx0,vy0,r_x,r_y,r_vx,r_vy, basis_height, hat_shift,      &
    !             part_group )


  ! push particles with a measure-preserving affine flow

  ! loop over all particles (taken from [[file:../simulation/simulation_4d_vp_lt_pic_cartesian.F90]])

  do k = 1, part_group%number_particles

     ! -- --  push the k-th particle [begin]  -- --

     ! get particle position
     coords = particle_group%get_x(k) ! [[selalib:src/particle_methods/sll_pic_base.F90::get_v]]
     x_k = coords(1)
     y_k = coords(2)

     ! get particle speed
     coords = particle_group%get_v(k) ! [[selalib:src/particle_methods/sll_pic_base.F90::get_v]]
     vx_k = coords(1)
     vy_k = coords(2)

     ! Flow is cte + A*(x,y,vx,vy), where det A = 1, cf [[test_forward_push]]
     call test_forward_push(x_k, y_k, vx_k, vy_k, new_x_k, new_y_k, new_vx_k, new_vy_k)

     ! set particle position
     coords(1) = new_x_k
     coords(2) = new_y_k
     coords(3) = 0
     call particle_group%set_x(k, coords)


     ! set particle speed
     coords(1) = new_vx_k
     coords(2) = new_vy_k
     coords(3) = 0
     call particle_group%set_v(k, coords)

    ! SHOULD WE PUT THE PARTICLES BACK INTO THE DOMAIN ???
    ! (see the LD simulation if that is needed)

  end do

  ! remap au choix
  call getarg(1,remap_type)
  call sll_set_time_mark(remapstart)
  if(remap_type == 'ltp') then
     ! remap with [[file:lt_pic_4d_utilities.F90::sll_lt_pic_4d_write_f_on_remap_grid]]
     print*, "[unit_test_bsl_lt_pic_4d]  Error (875454367554242): ltp remapping not implemented yet, stop."
     stop
     print*, "[lt_pic_4d_init_tester]  OLD VERSION: calling sll_lt_pic_4d_write_f_on_remap_grid..."
     ! OLD: call sll_lt_pic_4d_write_f_on_remap_grid( part_group )

  else if (remap_type == 'bsl_ltp') then
     ! remap with [[??? file:lt_pic_4d_utilities.F90::sll_lt_pic_4d_write_bsl_f_on_remap_grid]]
     print*, "[unit_test_bsl_lt_pic_4d]  Error (8767848745765): please write down the f*$%ing remapping routine!"
     stop

    ! OLD    call sll_lt_pic_4d_write_f_on_remapping_grid( part_group, 2, 2, 2, 2)

  else
     print*, 'ERROR (code=656536756757657): option is ltp (WARNING: not implemented yet) or bsl_ltp'
     stop
  end if
  remaptime = sll_time_elapsed_since(remapstart)

  ! formats at [[http://www.cs.mtu.edu/~shene/COURSES/cs201/NOTES/chap05/format.html]]
  write(*,'(A,A,ES8.2,A,A,ES8.2,A)') trim(remap_type),' remap time = ',remaptime,' sec',&
       ' ie ',part_group%number_particles/remaptime,' remapped-ptc/sec'
  
  ! result should still be exact on the new grid because all the transformations have been linear

  number_nodes_x  = part_group%number_parts_x
  number_nodes_y  = part_group%number_parts_y
  number_nodes_vx = part_group%number_parts_vx
  number_nodes_vy = part_group%number_parts_vy

  h_nodes_x    = part_group%remapping_grid%delta_eta1
  h_nodes_y    = part_group%remapping_grid%delta_eta2
  h_nodes_vx   = part_group%remapping_grid%delta_eta3
  h_nodes_vy   = part_group%remapping_grid%delta_eta4
    
  nodes_x_min    = part_group%remapping_grid%eta1_min
  nodes_y_min    = part_group%remapping_grid%eta2_min
  nodes_vx_min   = part_group%remapping_grid%eta3_min
  nodes_vy_min   = part_group%remapping_grid%eta4_min

  open(80,file='target_values_test_init4D.dat')
  error = 0
  x_j = nodes_x_min
  do j_x = 1, number_nodes_x
    y_j = nodes_y_min
    do j_y = 1, number_nodes_y
      vx_j = nodes_vx_min
      do j_vx = 1, number_nodes_vx
        vy_j = nodes_vy_min
        do j_vy = 1, number_nodes_vy
          f_j = real( part_group%target_values(j_x,j_y,j_vx,j_vy) ,f32)
          call test_backward_push(x_j, y_j, vx_j, vy_j, new_x, new_y, new_vx, new_vy)
          if(.not.in_bounds_periodic( new_x, new_y, part_group%mesh,DOMAIN_IS_X_PERIODIC,DOMAIN_IS_Y_PERIODIC )) then
             call apply_periodic_bc( part_group%mesh, new_x, new_y)
          end if

          f_target = eval_hat_function(x0,y0,vx0,vy0,r_x,r_y,r_vx,r_vy, basis_height, hat_shift, new_x, new_y, new_vx, new_vy)

          !! MCP: old implementation
          !          f_target = max_f * max(0._f64, 1.-inv_r_x*abs(new_x-x0) )       &
          !                           * max(0._f64, 1.-inv_r_y*abs(new_y-y0) )       &
          !                           * max(0._f64, 1.-inv_r_vx*abs(new_vx-vx0) )    &
          !                           * max(0._f64, 1.-inv_r_vy*abs(new_vy-vy0) )

          error = max( error, abs( f_j - f_target ) )

          ! MCP: [DEBUG] print the (computed) absolute initial position of the virtual particles (used in the bsl reconstruction)
          !x_j_bsl_before_per = part_group%debug_bsl_remap(j_x,j_y,j_vx,j_vy,1,1)
          !y_j_bsl_before_per = part_group%debug_bsl_remap(j_x,j_y,j_vx,j_vy,2,1)
          !vx_j_bsl_before_per = part_group%debug_bsl_remap(j_x,j_y,j_vx,j_vy,3,1)
          !vy_j_bsl_before_per = part_group%debug_bsl_remap(j_x,j_y,j_vx,j_vy,4,1)

          !x_j_bsl = part_group%debug_bsl_remap(j_x,j_y,j_vx,j_vy,1,2)
          !y_j_bsl = part_group%debug_bsl_remap(j_x,j_y,j_vx,j_vy,2,2)
          !vx_j_bsl = part_group%debug_bsl_remap(j_x,j_y,j_vx,j_vy,3,2)
          !vy_j_bsl = part_group%debug_bsl_remap(j_x,j_y,j_vx,j_vy,4,2)

          !x_j_bsl_end = part_group%debug_bsl_remap(j_x,j_y,j_vx,j_vy,1,3)
          !y_j_bsl_end = part_group%debug_bsl_remap(j_x,j_y,j_vx,j_vy,2,3)
          !vx_j_bsl_end = part_group%debug_bsl_remap(j_x,j_y,j_vx,j_vy,3,3)
          !vy_j_bsl_end = part_group%debug_bsl_remap(j_x,j_y,j_vx,j_vy,4,3)

          ! write(80,*) x_j, y_j, vx_j, vy_j, f_j, f_target, abs( f_j - f_target ) &
          !              , x_j_bsl_before_per, y_j_bsl_before_per, vx_j_bsl_before_per, vy_j_bsl_before_per &      ! MCP : this line to DEBUG only
          !              , x_j_bsl, y_j_bsl, vx_j_bsl, vy_j_bsl &      ! MCP : this line to DEBUG only
          !              , x_j_bsl_end, y_j_bsl_end, vx_j_bsl_end, vy_j_bsl_end       ! MCP : this line to DEBUG only


          vy_j = vy_j + h_nodes_vy
        end do
        vx_j = vx_j + h_nodes_vx
      end do
      y_j = y_j + h_nodes_y
    end do
    x_j = x_j + h_nodes_x
  end do

  close(80)

  ! uses [[file:~/mcp/selalib/src/fields/sll_array_plotting_module.F90::write_projection_2d]] developed by PN (cf
  ! example in [[file:~/mcp/selalib/src/fields/unit_test_4d.F90::write_projection_2d]])

  plotting_mesh_2d =>  new_cartesian_mesh_2d( NC_X1_PLOT, NC_X2_PLOT, &
         X1_MIN_PLOT, X1_MAX_PLOT, X2_MIN_PLOT, X2_MAX_PLOT )

  file_name = "lt_pic_density_test_init4D.dat"

  ! [[file:lt_pic_4d_utilities.F90::sll_lt_pic_4d_plot_f_on_2d_grid]]
  call sll_lt_pic_4d_plot_f_on_2d_grid( &
    file_name,        &
    PLOT_DIM1,        &
    PLOT_DIM2,        &
    PLOT_CST_DIM3,    &
    PLOT_CST_DIM4,    &
    X3_PLOT_CST,      &
    X4_PLOT_CST,      &      
    plotting_mesh_2d,     &
    part_group        &
    )


  write(ncx_name,'(i3)') NC_X
  write(ncy_name,'(i3)') NC_Y
  open(83,file='vit_pos_'//trim(adjustl(ncx_name))//'x'//trim(adjustl(ncy_name))//'.dat')
  do j = 1, part_group%number_particles
    ! call cell_offset_to_global ( part_group%p_list(j)%dx, &
    !                            part_group%p_list(j)%dy, &
    !                            part_group%p_list(j)%ic, &
    !                            mesh_2d, x, y )
     call cell_offset_to_global_extended (  part_group%p_list(j)%dx, &
                                            part_group%p_list(j)%dy, &
                                            part_group%p_list(j)%ic_x, &
                                            part_group%p_list(j)%ic_y, &
                                            mesh_2d, x, y )

     write(83,*) part_group%p_list(j)%vx, part_group%p_list(j)%vy, &
     part_group%p_list(j)%ic_x, part_group%p_list(j)%ic_y, part_group%p_list(j)%dx, part_group%p_list(j)%dy, &
     x, y
  enddo
  close(83)

  open(82,file='error_test_lt_pic_init4D.dat')
  write(82,*) "(maximum) error : ", error
  close(82)
  
  call sll_delete( part_group )
  call sll_delete( mesh_2d )
  tolerance = 1e-6
  if( error < tolerance )then
    print*, "TEST PASSED: error = ", error
  else
    print*, "TEST FAILED: error = ", error
  endif

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

  subroutine test_forward_push(x,y,vx,vy,new_x,new_y,new_vx,new_vy)

    ! [[file:../working_precision/sll_working_precision.h]]
    use sll_working_precision

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

  subroutine test_backward_push(x, y, vx, vy, new_x, new_y, new_vx, new_vy)

    ! [[file:../working_precision/sll_working_precision.h]]
    use sll_working_precision

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

  end subroutine test_backward_push
end program
