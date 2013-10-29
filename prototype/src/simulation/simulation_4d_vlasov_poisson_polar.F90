!**************************************************************
!  Copyright INRIA
!  Authors : 
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
!> @author Pierre Navaro
!> @brief 
!> Simulation class to solve vlasov poisson system in polar coordinates
!> (2d space 2d phase)
!> @details
!> Example of use in test program 
!> 
!> \code
!>
!>  use sll_simulation_4d_vlasov_poisson_polar
!>  type(sll_simulation_4d_vp_polar)    :: simulation
!>  class(sll_coordinate_transformation_2d_base), pointer :: transformation
!>  call simulation%init_from_file(trim(filename))
!>
!>  call initialize_vp4d_polar( &
!>       simulation, &
!>       mx, &
!>       mv, &
!>       transformation, &
!>       sll_periodic_gaussian_initializer_4d, &
!>       params )
!>
!>  call simulation%run()
!>  call delete(simulation)
!> \endcode

#define MPI_MASTER 0

module sll_simulation_4d_vlasov_poisson_polar

#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_field_2d.h"
#include "sll_utilities.h"

use sll_collective
use sll_remapper
use sll_constants
use sll_cubic_spline_interpolator_2d
use sll_poisson_2d_polar
use sll_cubic_spline_interpolator_1d
use sll_simulation_base
use sll_logical_meshes
use sll_parallel_array_initializer_module
use sll_coordinate_transformation_2d_base_module
use sll_gnuplot_parallel
use sll_poisson_polar_parallel

implicit none

!> vp4d polar simulation class extended from sll_simulation_base_class
type, extends(sll_simulation_base_class) :: sll_simulation_4d_vp_polar
  
 sll_int32  :: world_size !< Parallel environment parameters
 sll_int32  :: my_rank    !< Processor id
 sll_real64 :: dt              !< time step
 sll_int32  :: num_iterations  !< steps number
 sll_int32  :: nc_x1, nc_x2, nc_x3, nc_x4 !< Mesh parameters
 type(sll_logical_mesh_2d), pointer :: mesh2d_x !< the logical mesh for space
 type(sll_logical_mesh_2d), pointer :: mesh2d_v !< the logical mesh for velocity
 class(sll_coordinate_transformation_2d_base), pointer :: transfx !< coordinate transformation
 
 sll_real64, dimension(:,:,:,:), pointer :: f_x1x2 !< sequential in x1 and x2
 sll_real64, dimension(:,:,:,:), pointer :: f_x3x4 !< sequential in x3 and x4
 
 sll_real64, dimension(:,:), pointer :: proj_f_x1x2 !< f projection to x1x2
 sll_real64, dimension(:,:), pointer :: proj_f_x3x4 !< f projection to x3x4
      
 type(layout_4D), pointer :: sequential_x1x2 !< layout 4d sequential in x1x2
 type(layout_4D), pointer :: sequential_x3x4 !< layout 4d sequential in x3x4
 type(remap_plan_4D_real64), pointer :: seqx1x2_to_seqx3x4 !< transpose x to v
 type(remap_plan_4D_real64), pointer :: seqx3x4_to_seqx1x2 !< transpose v to x 
 
 type(cubic_spline_2d_interpolator) :: interp_x1x2 !< interpolator 2d in xy
 type(cubic_spline_1d_interpolator) :: interp_x3   !< interpolator 1d in vx
 type(cubic_spline_1d_interpolator) :: interp_x4   !< interpolator 1d in vx
 procedure(sll_scalar_initializer_4d), nopass, pointer :: init_func !< for distribution function initializer:
 sll_real64, dimension(:),     pointer :: params !< function initializer parameters
 sll_real64, dimension(:,:),   pointer :: x1     !< x1 mesh mapped coordinates
 sll_real64, dimension(:,:),   pointer :: x2     !< x2 mesh mapped coordinates
 sll_real64, dimension(:,:),   pointer :: phi_x1    !< potential
 sll_real64, dimension(:,:),   pointer :: phi_x2    !< potential
 sll_real64, dimension(:,:),   pointer :: rho    !< charge density
 sll_real64, dimension(:,:,:), pointer :: partial_reduction
 sll_real64, dimension(:,:,:), pointer :: efields_x1
 sll_real64, dimension(:,:,:), pointer :: efields_x2


 type(sll_poisson_polar)  :: poisson   ! Poisson solver in polar coordinates
 type(layout_2D), pointer :: layout_x1 ! sequential in r direction
 type(layout_2D), pointer :: layout_x2 ! sequential in theta direction
 type(remap_plan_2D_real64), pointer :: rmp_x1x2   !< remap r->theta 
 type(remap_plan_2D_real64), pointer :: rmp_x2x1   !< remap theta->r

contains
 
 procedure, pass(sim) :: run => run_vp4d_polar !< run the simulation
 procedure, pass(sim) :: init_from_file => init_vp4d_par_polar !< init the simulation

end type sll_simulation_4d_vp_polar

interface delete
   module procedure delete_vp4d_par_polar
end interface delete

interface initialize
   module procedure initialize_vp4d_polar
end interface initialize

!> Local variables
sll_int32,  private :: i, j, k, l
sll_int32,  private :: loc_sz_x1, loc_sz_x2, loc_sz_x3, loc_sz_x4 
sll_int32,  private :: global_indices(4)
sll_int32,  private :: gi, gj, gk, gl
sll_real64, private :: delta_eta1, delta_eta2, delta_eta3, delta_eta4
sll_real64, private :: alpha1, alpha2, alpha3, alpha4
sll_real64, private :: eta1, eta2, eta3, eta4
sll_int32,  private :: nc_eta1, nc_eta2, nc_eta3, nc_eta4
sll_real64, private :: jac_m(2,2), inv_j(2,2)
sll_int32,  private :: itime, error
sll_real64, private :: eta1_min, eta2_min, eta3_min, eta4_min
sll_real64, private :: eta1_max, eta2_max, eta3_max, eta4_max

character(len=4), private :: ctime

contains

  !> Function to initialize the simulation object 'manually'.
  subroutine initialize_vp4d_polar( &
   sim, &
   mesh2d_x, &
   mesh2d_v, &
   transformation_x, &
   init_func, &
   params )

   type(sll_simulation_4d_vp_polar), intent(inout)       :: sim
   type(sll_logical_mesh_2d), pointer                    :: mesh2d_x
   type(sll_logical_mesh_2d), pointer                    :: mesh2d_v
   class(sll_coordinate_transformation_2d_base), pointer :: transformation_x
   procedure(sll_scalar_initializer_4d)                  :: init_func
   sll_real64, dimension(:), target                      :: params
   sim%mesh2d_x  => mesh2d_x
   sim%mesh2d_v  => mesh2d_v
   sim%transfx   => transformation_x
   sim%init_func => init_func
   sim%params    => params

   sim%nc_x1 = mesh2d_x%num_cells1
   sim%nc_x2 = mesh2d_x%num_cells2
   sim%nc_x3 = mesh2d_v%num_cells1
   sim%nc_x4 = mesh2d_v%num_cells2

  end subroutine initialize_vp4d_polar

  !> Function to initialize the simulation object from a file.
  subroutine init_vp4d_par_polar( sim, filename )
    class(sll_simulation_4d_vp_polar), intent(inout) :: sim !< simulation class
    character(len=*), intent(in)                     :: filename !< input file name
    sll_int32  :: IO_stat
    sll_real64 :: dt
    sll_int32  :: number_iterations
    sll_int32  :: num_cells_x1
    sll_int32  :: num_cells_x2
    sll_int32  :: num_cells_x3
    sll_int32  :: num_cells_x4
    sll_int32  :: input_file 

    namelist /sim_params/ dt, number_iterations
    namelist /grid_dims/ num_cells_x1, num_cells_x2, num_cells_x3, num_cells_x4
    
    call sll_new_file_id(input_file, error)
    open(unit = input_file, file=filename,IOStat=IO_stat)
    if( IO_stat /= 0 ) then
       print *, 'init_vp4d_par_cart() failed to open file ', filename
       stop
    end if
    read(input_file, sim_params)
    read(input_file,grid_dims)
    close(input_file)

    sim%dt = dt
    sim%num_iterations = number_iterations
    ! In this particular simulation, since the system is periodic, the number
    ! of points is the same as the number of cells in all directions.
    sim%nc_x1 = num_cells_x1
    sim%nc_x2 = num_cells_x2
    sim%nc_x3 = num_cells_x3
    sim%nc_x4 = num_cells_x4

  end subroutine init_vp4d_par_polar

  !> run simulation
  subroutine run_vp4d_polar(sim)
    class(sll_simulation_4d_vp_polar), intent(inout) :: sim
    sll_int32  :: error
    sll_int32  :: prank
    sll_int32  :: psize

    psize  = sll_get_collective_size(sll_world_collective)
    prank  = sll_get_collective_rank(sll_world_collective)

    sim%my_rank = prank

    ! allocate the layouts...
    sim%sequential_x1x2  => new_layout_4D( sll_world_collective )
    sim%sequential_x3x4  => new_layout_4D( sll_world_collective )

    nc_eta1    = sim%nc_x1
    nc_eta2    = sim%nc_x2 
    nc_eta3    = sim%nc_x3
    nc_eta4    = sim%nc_x4 

    delta_eta1 = sim%mesh2d_x%delta_eta1
    delta_eta2 = sim%mesh2d_x%delta_eta2
    delta_eta3 = sim%mesh2d_v%delta_eta1
    delta_eta4 = sim%mesh2d_v%delta_eta2

    eta1_min   = sim%mesh2d_x%eta1_min
    eta2_min   = sim%mesh2d_x%eta2_min
    eta3_min   = sim%mesh2d_v%eta1_min
    eta4_min   = sim%mesh2d_v%eta2_min

    eta1_max   = sim%mesh2d_x%eta1_max
    eta2_max   = sim%mesh2d_x%eta2_max
    eta3_max   = sim%mesh2d_v%eta1_max
    eta4_max   = sim%mesh2d_v%eta2_max

    call initialize_layout_with_distributed_4D_array( &
         sim%nc_x1+1, sim%nc_x2+1, sim%nc_x3+1, sim%nc_x4+1, &
         psize, 1, 1, 1, &
         sim%sequential_x3x4 )

    call initialize_layout_with_distributed_4D_array( &
         sim%nc_x1+1, sim%nc_x2+1, sim%nc_x3+1, sim%nc_x4+1, &
         1, 1, psize, 1, &
         sim%sequential_x1x2 )
    
    call compute_local_sizes_4d( sim%sequential_x1x2, &
         loc_sz_x1, loc_sz_x2, loc_sz_x3, loc_sz_x4 )

    SLL_ALLOCATE(sim%f_x1x2(loc_sz_x1,loc_sz_x2,loc_sz_x3,loc_sz_x4),error)
    SLL_ALLOCATE(sim%proj_f_x3x4(loc_sz_x3,loc_sz_x4),error)
    
    call compute_local_sizes_4d( sim%sequential_x3x4, &
         loc_sz_x1, loc_sz_x2, loc_sz_x3, loc_sz_x4 )

    SLL_ALLOCATE(sim%f_x3x4(loc_sz_x1,loc_sz_x2,loc_sz_x3,loc_sz_x4),error)
    SLL_ALLOCATE(sim%proj_f_x1x2(loc_sz_x1,loc_sz_x2),error)

    
    if ( prank == MPI_MASTER ) call sll_view_lims_4D( sim%sequential_x1x2 )
    call flush(6)
    if ( prank == MPI_MASTER ) call sll_view_lims_4D( sim%sequential_x3x4 )
    call flush(6)
    
    call sll_4d_parallel_array_initializer( &
         sim%sequential_x3x4, sim%mesh2d_x, sim%mesh2d_v, &
         sim%f_x3x4, sim%init_func, sim%params, &
         transf_x1_x2=sim%transfx )

    sim%seqx3x4_to_seqx1x2 => &
         NEW_REMAP_PLAN(sim%sequential_x3x4,sim%sequential_x1x2,sim%f_x3x4)

    sim%seqx1x2_to_seqx3x4 => &
         NEW_REMAP_PLAN(sim%sequential_x1x2,sim%sequential_x3x4,sim%f_x1x2)

    call sim%interp_x1x2%initialize( &
         sim%nc_x1+1, sim%nc_x2+1,   &
         eta1_min, eta1_max,         &
         eta2_min, eta2_max,         &
         SLL_PERIODIC, SLL_PERIODIC )

    call sim%interp_x3%initialize( &
         sim%nc_x3+1, eta3_min, eta3_max, SLL_PERIODIC)

    call sim%interp_x4%initialize( &
         sim%nc_x4+1, eta4_min, eta4_max, SLL_PERIODIC)

    sim%layout_x1 => new_layout_2D( sll_world_collective )

    call initialize_layout_with_distributed_2D_array( sim%nc_x1+1, &
                                                      sim%nc_x2+1, &
                                                      1,           &
                                                      psize,       &
                                                      sim%layout_x1 )

    call compute_local_sizes_2d(sim%layout_x1, loc_sz_x1, loc_sz_x2)
    SLL_CLEAR_ALLOCATE(sim%phi_x1(1:loc_sz_x1,1:loc_sz_x2),error)
    SLL_CLEAR_ALLOCATE(sim%efields_x1(1:loc_sz_x1,1:loc_sz_x2,2),error)

    sim%layout_x2 => new_layout_2D( sll_world_collective )

    call initialize_layout_with_distributed_2D_array( sim%nc_x1+1, &
                                                      sim%nc_x2+1, &
                                                      psize,       &
                                                      1,           &
                                                      sim%layout_x2 )

    call compute_local_sizes_2d(sim%layout_x2, loc_sz_x1, loc_sz_x2)
    SLL_CLEAR_ALLOCATE(sim%rho(1:loc_sz_x1,1:loc_sz_x2),error)
    SLL_CLEAR_ALLOCATE(sim%phi_x2(1:loc_sz_x1,1:loc_sz_x2),error)
    SLL_CLEAR_ALLOCATE(sim%efields_x2(1:loc_sz_x1,1:loc_sz_x2,2),error)

    call initialize( sim%poisson,   &
                     sim%layout_x1, &
                     sim%layout_x2, &
                     eta1_min,      &
                     eta1_max,      &
                     sim%nc_x1,     &
                     sim%nc_x2,     &
                     SLL_DIRICHLET, &
                     SLL_DIRICHLET)

    sim%rmp_x1x2 => new_remap_plan(sim%layout_x1, sim%layout_x2, sim%phi_x1)
    sim%rmp_x2x1 => new_remap_plan(sim%layout_x2, sim%layout_x1, sim%phi_x2)

    SLL_CLEAR_ALLOCATE(sim%x1(loc_sz_x1,loc_sz_x2),error)
    SLL_CLEAR_ALLOCATE(sim%x2(loc_sz_x1,loc_sz_x2),error)

    do j = 1, loc_sz_x2
    do i = 1, loc_sz_x1
       global_indices(1:2) = local_to_global_2D(sim%layout_x2,(/i,j/))
       sim%x1(i,j) = sim%transfx%x1_at_node(global_indices(1),global_indices(2))
       sim%x2(i,j) = sim%transfx%x2_at_node(global_indices(1),global_indices(2))
    end do
    end do
    
    do itime=1, sim%num_iterations !Loop over time

       if(sim%my_rank == 0) then
          print *, 'Starting iteration ', itime, ' of ', sim%num_iterations
       end if

       call plot_f_x1x2(sim)

       !transpose to x space
       call apply_remap_4D( sim%seqx3x4_to_seqx1x2, sim%f_x3x4, sim%f_x1x2 )

       call advection_x1x2(sim,sim%dt)


       call compute_charge_density( sim )

       call plot_rho(sim)

       call solve(sim%poisson, sim%rho, sim%phi_x2)

       call plot_phi(sim)

       !transpose to v space
       call apply_remap_4D( sim%seqx1x2_to_seqx3x4, sim%f_x1x2, sim%f_x3x4 )

       !call apply_remap_2D( sim%rmp_x2x1, sim%phi_x2, sim%phi_x1 )
       !call compute_electric_fields_eta1( sim )
       !call apply_remap_2D( sim%rmp_x1x2, sim%phi_x1, sim%phi_x2 )
       !call compute_electric_fields_eta2( sim )
!
!       call advection_x3(sim,sim%dt)

!       call advection_x4(sim,sim%dt)

       call plot_f_x3x4(sim)

    end do ! next time step 

  end subroutine run_vp4d_polar

  subroutine advection_x1x2(sim,deltat)
    class(sll_simulation_4d_vp_polar) :: sim
    sll_real64, intent(in) :: deltat

    call compute_local_sizes_4d( sim%sequential_x1x2, &
         loc_sz_x1, loc_sz_x2, loc_sz_x3, loc_sz_x4 )

    do l=1,loc_sz_x4
       do k=1,loc_sz_x3
          call sim%interp_x1x2%compute_interpolants(sim%f_x1x2(:,:,k,l))
          do j=1,loc_sz_x2
             do i=1,loc_sz_x1
                global_indices = local_to_global_4D(sim%sequential_x1x2,(/i,j,k,l/))
                gi = global_indices(1)
                gj = global_indices(2)
                gk = global_indices(3)
                gl = global_indices(4)
                eta1 = eta1_min + (gi-1)*delta_eta1
                eta2 = eta2_min + (gj-1)*delta_eta2
                eta3 =  eta1*sin(eta2) !sim%mesh2d_v%eta1_min + (gk-1)*delta_eta3
                eta4 = -eta1*cos(eta2) !sim%mesh2d_v%eta2_min + (gl-1)*delta_eta4
                inv_j  = sim%transfx%inverse_jacobian_matrix(eta1,eta2)
                alpha1 = -deltat*(inv_j(1,1)*eta3 + inv_j(1,2)*eta4)
                alpha2 = -deltat*(inv_j(2,1)*eta3 + inv_j(2,2)*eta4)

                eta1 = eta1+alpha1
                ! This is hardwiring the periodic BC, please improve this...
                if( eta1 < eta1_min ) then
                   eta1 = eta1+eta1_max-eta1_min
                else if( eta1 > eta1_max ) then
                   eta1 = eta1+eta1_min-eta1_max
                end if

                eta2 = eta2+alpha2
                if( eta2 < eta2_min ) then
                   eta2 = eta2+eta2_max-eta2_min
                else if( eta2 > eta2_max ) then
                   eta2 = eta2+eta2_min-eta2_max
                end if
                
                sim%f_x1x2(i,j,k,l) = sim%interp_x1x2%interpolate_value(eta1,eta2)

             end do
          end do
       end do
    end do

  end subroutine advection_x1x2

  subroutine advection_x3(sim,deltat)

    class(sll_simulation_4d_vp_polar) :: sim
    sll_real64, intent(in) :: deltat
    sll_real64 :: ex, ey

    call compute_local_sizes_4d( sim%sequential_x3x4, &
            loc_sz_x1, loc_sz_x2, loc_sz_x3, loc_sz_x4 ) 

    do l=1,loc_sz_x4
       do j=1,loc_sz_x2
          do i=1,loc_sz_x1
             global_indices = local_to_global_4D( sim%sequential_x3x4, (/i,j,1,1/))
             eta1   =  eta1_min + (global_indices(1)-1)*delta_eta1
             eta2   =  eta2_min + (global_indices(2)-1)*delta_eta2
             inv_j  =  sim%transfx%inverse_jacobian_matrix(eta1,eta2)
             jac_m  =  sim%transfx%jacobian_matrix(eta1,eta2)
             ex     =  sim%efields_x2(i,j,1)
             ey     =  sim%efields_x2(i,j,2)
             alpha3 = -deltat*(inv_j(1,1)*ex + inv_j(2,1)*ey)
             sim%f_x3x4(i,j,:,l) = sim%interp_x3%interpolate_array_disp( &
                                   sim%nc_x3+1, sim%f_x3x4(i,j,:,l), alpha3 )
          end do
       end do
    end do

  end subroutine advection_x3

  subroutine advection_x4(sim,deltat)

    class(sll_simulation_4d_vp_polar) :: sim
    sll_real64, intent(in) :: deltat
    sll_real64 :: ex, ey

    call compute_local_sizes_4d( sim%sequential_x3x4, &
         loc_sz_x1, loc_sz_x2, loc_sz_x3, loc_sz_x4 ) 

    do j=1,loc_sz_x2
       do i=1,loc_sz_x1
          do k=1,loc_sz_x3
             global_indices = local_to_global_4D( sim%sequential_x3x4, (/i,j,1,1/))
             eta1   =  eta1_min+(global_indices(1)-1)*delta_eta1
             eta2   =  eta2_min+(global_indices(2)-1)*delta_eta2
             inv_j  =  sim%transfx%inverse_jacobian_matrix(eta1,eta2)
             ex     =  sim%efields_x2(i,j,1)
             ey     =  sim%efields_x2(i,j,2)
             alpha4 = -deltat*(inv_j(1,2)*ex + inv_j(2,2)*ey)
             sim%f_x3x4(i,j,k,:) = sim%interp_x4%interpolate_array_disp( &
                                   sim%nc_x4+1, sim%f_x3x4(i,j,k,:), alpha4 )
          end do
       end do
    end do

  end subroutine advection_x4

  subroutine plot_f_x1x2(sim)
    class(sll_simulation_4d_vp_polar) :: sim


    call compute_local_sizes_4d( sim%sequential_x3x4, &
         loc_sz_x1, loc_sz_x2, loc_sz_x3, loc_sz_x4 )

    do i = 1, loc_sz_x1
       do j = 1, loc_sz_x2
          sim%proj_f_x1x2(i,j) = sum(sim%f_x3x4(i,j,:,:))
       end do
    end do

    call sll_gnuplot_2d_parallel( sim%x1, sim%x2, sim%proj_f_x1x2, &
                                  "fxy", itime, error )

  end subroutine plot_f_x1x2

  subroutine plot_f_x3x4(sim)
    class(sll_simulation_4d_vp_polar) :: sim

    call compute_local_sizes_4d( sim%sequential_x1x2, &
         loc_sz_x1, loc_sz_x2, loc_sz_x3, loc_sz_x4 )

    do l = 1, loc_sz_x4
       do k = 1, loc_sz_x3
          sim%proj_f_x3x4(k,l) = sum(sim%f_x1x2(:,:,k,l))
       end do
    end do

    call sll_gnuplot_2d_parallel( &
        eta3_min+(global_indices(1)-1)*delta_eta3, delta_eta3, &
        eta4_min+(global_indices(2)-1)*delta_eta4, delta_eta4, &
        sim%proj_f_x3x4, "fvxvy", itime, error )

  end subroutine plot_f_x3x4

  subroutine delete_vp4d_par_polar( sim )
    type(sll_simulation_4d_vp_polar) :: sim
    sll_int32 :: error
    SLL_DEALLOCATE( sim%f_x1x2, error )
    SLL_DEALLOCATE( sim%f_x3x4, error )
    call delete( sim%sequential_x1x2 )
    call delete( sim%sequential_x3x4 )
    call delete( sim%seqx1x2_to_seqx3x4 )
    call delete( sim%seqx3x4_to_seqx1x2 )
    call delete( sim%interp_x1x2 )
    call delete( sim%interp_x3 )
    call delete( sim%interp_x4 )
  end subroutine delete_vp4d_par_polar

  subroutine compute_charge_density( sim )

    class(sll_simulation_4d_vp_polar) :: sim

    call compute_local_sizes_4d(sim%sequential_x3x4, loc_sz_x1, &
                                                     loc_sz_x2, &
                                                     loc_sz_x3, &
                                                     loc_sz_x4)
    sim%rho(:,:) = 0.0
    
    do j=1,loc_sz_x2
       do i=1,loc_sz_x1
          sim%rho(i,j) = sum(sim%f_x3x4(i,j,:,:))
       end do
    end do
    sim%rho = sim%rho *delta_eta3*delta_eta4

  end subroutine compute_charge_density

  
  subroutine plot_rho(sim)

    class(sll_simulation_4d_vp_polar) :: sim
    call compute_local_sizes_2d(sim%layout_x2, loc_sz_x1, loc_sz_x2)

    call sll_gnuplot_2d_parallel(sim%x1, sim%x2, sim%rho, 'rho', itime, error)

  end subroutine plot_rho
 
  subroutine plot_phi(sim)

    class(sll_simulation_4d_vp_polar) :: sim
    call compute_local_sizes_2d(sim%layout_x2, loc_sz_x1, loc_sz_x2)

    call sll_gnuplot_2d_parallel(sim%x1, sim%x2, sim%phi_x2, 'phi', itime, error)

  end subroutine plot_phi

  subroutine compute_electric_fields_eta1( sim )

    class(sll_simulation_4d_vp_polar) :: sim
    sll_real64                        :: r_delta
    
    r_delta = 1.0_f64/delta_eta1
    
    i = 1
    sim%efields_x1(i,:,1) = -r_delta*(- 1.5_f64*sim%phi_x1(i  ,:) &
                                   + 2.0_f64*sim%phi_x1(i+1,:) &
                                   - 0.5_f64*sim%phi_x1(i+2,:) )
    i = sim%nc_x1+1
    sim%efields_x1(i,:,1) = -r_delta*(  0.5_f64*sim%phi_x1(i-2,:) &
                                   - 2.0_f64*sim%phi_x1(i-1,:) &
                                   + 1.5_f64*sim%phi_x1(i  ,:) )
    do i=2, sim%nc_x1
       sim%efields_x1(i,:,1) = -r_delta*0.5_f64*(sim%phi_x1(i+1,:) &
                                            - sim%phi_x1(i-1,:))
    end do

    call apply_remap_2D( sim%rmp_x1x2, sim%efields_x1(:,:,1), sim%efields_x2(:,:,1) )

  end subroutine compute_electric_fields_eta1

  subroutine compute_electric_fields_eta2( sim )

    class(sll_simulation_4d_vp_polar) :: sim
    sll_real64                        :: r_delta
  
    r_delta = 1.0_f64/delta_eta2
    
    j=1 
    sim%efields_x2(:,j,2) = -r_delta*(- 1.5_f64*sim%phi_x2(:,j)   &
                                     + 2.0_f64*sim%phi_x2(:,j+1) &
                                     - 0.5_f64*sim%phi_x2(:,j+2))
    j=sim%nc_x2+1
    sim%efields_x2(:,j,2) = -r_delta*(  0.5_f64*sim%phi_x2(:,j-2) &
                                     - 2.0_f64*sim%phi_x2(:,j-1) &
                                     + 1.5_f64*sim%phi_x2(:,j))
    
    do j=2,sim%nc_x2
       sim%efields_x2(:,j,2) = -r_delta*0.5_f64*(sim%phi_x2(:,j+1) &
                                              - sim%phi_x2(:,j-1))
    end do

    call apply_remap_2D( sim%rmp_x2x1, sim%efields_x2(:,:,2), sim%efields_x1(:,:,2) )

  end subroutine compute_electric_fields_eta2


end module sll_simulation_4d_vlasov_poisson_polar
