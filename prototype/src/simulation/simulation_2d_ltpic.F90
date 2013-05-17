!> \file simulation_2d_ltpic.F90
!> \authors                    
!> Martin CAMPOS PINTO (campos@ann.jussieu.fr) 

module sll_simulation_2d_ltpic
#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"

  use sll_constants
  use sll_tsi_2d_initializer
  use sll_simulation_base
  use sll_poisson_1d_periodic
  use sll_lin_trans_pic_2d
  
  implicit none

  type, extends(sll_simulation_base_class) :: sll_simulation_2d_ltpic

     ! Testcase parameters
     sll_real64 :: final_time
     sll_real64 :: kx
     sll_real64 :: epsilon

     ! Numerical parameters
     sll_real64 :: dt
     sll_real64 :: time_before_remap 

     ! Collection of particles with linear transformations 
     type(lin_trans_pic_2d), pointer    :: particles

   contains
     procedure, pass(sim) :: run => run_2d_ltpic
  end type sll_simulation_2d_ltpic


contains


  subroutine run_2d_ltpic(sim)
    implicit none
    class(sll_simulation_2d_ltpic),     intent(inout) :: sim
    type(init_tsi_2d),                  target        :: init_tsi
    type(sll_mapped_mesh_2d_cartesian), target        :: mesh
    class(sll_mapped_mesh_2d_base),     pointer       :: mesh_base
    type(leap_frog_1st_flow_2d),        target        :: lf_1st_flow
    type(leap_frog_2nd_flow_2d),        target        :: lf_2nd_flow
    class(flow_base_class),             pointer       :: flow
    type(poisson_1d_periodic)                         :: poisson
    sll_real64, dimension(:),           allocatable   :: rho
                                                      
    sll_int32                                         :: n_time
    sll_int32                                         :: num_iterations
    sll_int32                                         :: num_particles_x
    sll_int32                                         :: num_particles_v
    sll_real64                                        :: xmin
    sll_real64                                        :: xmax
    sll_real64                                        :: vmin
    sll_real64                                        :: vmax
    sll_int32                                         :: bspline_degree
    sll_int32                                         :: ierr
    sll_int32                                         :: nc_poisson_mesh
    sll_real64                                        :: xmin_poisson_mesh
    sll_real64                                        :: xmax_poisson_mesh
    sll_real64                                        :: elementary_charge
                                           
    ! testcase parameters
    sim%kx = 0.5_f64
    sim%epsilon = 0.01_f64
    sim%final_time = 45
    elementary_charge = 1.
        
    ! numerical parameters (particles)
    num_particles_x = 500
    num_particles_v = 500
    xmin = 0.0
    xmax = 2*sll_pi/kx
    vmin = -6.0
    vmax =  6.0
    bspline_degree = 3

    ! numerical parameters (poisson mesh)
    nc_poisson_mesh = 128 ! number of cells
    xmin_poisson_mesh = xmin
    xmax_poisson_mesh = xmax
    
    ! numerical parameters (time resolution)
    num_iterations = 100
    sim%time_before_remap = 10 
    sim%dt = sim%final_time/num_iterations

    print *, 'initializing the particles...'
    sim%particles => new_ltpic_2d( &
                                    num_particles_x,                        &
                                    num_particles_v,                        &
                                    xmin,                                   &
                                    xmax,                                   &
                                    vmin,                                   &
                                    vmax,                                   &
                                    bspline_degree,                         &
                                    nc_poisson_mesh,                        &
                                    xmin_poisson_mesh,                      &
                                    xmax_poisson_mesh,                      &
                                    elementary_charge,                      &
                                    PERIODIC_LTPIC                          &
                                  )

    print *, 'building a mesh for the initializer...'
    call mesh%initialize( &
         "map_cartesian", &
         sim%particles%qi_grid_xmin, &
         sim%particles%qi_grid_xmax, &
         sim%particles%qi_grid_npx,  &
         sim%particles%qi_grid_vmin, &
         sim%particles%qi_grid_vmax, &
         sim%particles%qi_grid_npv    )         
    mesh_base => mesh
    call init_tsi%initialize(mesh_base,NODE_CENTERED_FIELD,sim%epsilon,sim%kx,0_f64)

    print *, 'initializing the leap-frog flows...'
    call init_lf_1st_flow( lf_1st_flow, sim%dt )
    ! MCP -- here, the array describing the electric field coefficients is a member of the lf_2nd_flow object. 
    ! but I think it should rather have a specific type -- and be a member of some poisson solver class ? 
    ! (then passed as an argument to the flow initializer?)
    call init_lf_2nd_flow( lf_2nd_flow, sim%dt, sim%nc_poisson_mesh, sim%xmin_poisson_mesh, sim%xmax_poisson_mesh)

    print *, 'initializing the poisson solver...'
    SLL_ALLOCATE(rho(sim%nc_poisson_mesh+1),  ierr)
    call new(poisson, sim%xmin_poisson_mesh, sim%xmax_poisson_mesh, sim%nc_poisson_mesh,  ierr) 

    ! ------------------------------------------------------------------------
    !
    !                           COMPUTE INITIAL DATA
    !
    ! ------------------------------------------------------------------------

    print *, 'loading the particles...'
    call load_ltpic_2d( init_tsi,sim%particles )
    
    n_time_last_remap = 0
    n_time = 0
    
    print *, 'plotting the initial particle density...'
    call plot_density( 'initial_density',n_time,sim%particles )

    ! ------------------------------------------------------------------------
    !
    !                                MAIN LOOP
    !
    ! ------------------------------------------------------------------------

    do n_time=1, num_iterations
       print *, 'time step ', n_time, '/', num_iterations

       print *, 'transporting the particles (1st substep in leap-frog scheme)...'
       flow => lf_1st_flow
       call transport_ltpic_2d( flow,sim%particles )
       
       print *, 'computing the E field for the intermediate solution...'
       call deposit_charges( rho,sim%particles )
       call solve( poisson, lf_2nd_flow%electric_field, rho )

       print *, 'transporting the particles (2nd substep in leap-frog scheme)...'
       flow => lf_2nd_flow
       call transport_ltpic_2d( flow,sim%particles )

      ! conditional remapping
      if( n_time < num_iterations .and. (n_time - n_time_last_remap)*sim%dt > sim%time_before_remap ) then
         call remap_ltpic_2d( sim%particles )
         n_time_last_remap = n_time
      end if
    end do ! main loop

    print *, 'plotting the final particle density...'
    call plot_density( 'final_density',n_time,sim%particles )

    print *, 'done.'
    
  end subroutine run_2d_ltpic



end module sll_simulation_2d_ltpic
