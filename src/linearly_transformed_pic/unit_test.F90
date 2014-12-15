
!*****************************************************************************
!
! Selalib      
! Module: unit_test.F90
!
!> @brief 
!> sll_lin_trans_pic_2d unit test
!   
!> @authors                    
!> Martin CAMPOS PINTO (campos@ann.jussieu.fr) 
!                                  
!*****************************************************************************

program ltpic_tester
#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"
  
  use sll_lin_trans_pic_2d
  use sll_poly_2d_initializer
  use sll_module_mapped_meshes_2d_cartesian
  use sll_leap_frog_1st_flow_2d
  use sll_leap_frog_2nd_flow_2d
  use sll_poisson_1d_periodic
  implicit none

  logical                                :: test_passed
  logical                                :: test_flag
  test_passed = .true.
  

  print *, '***************************************************'
  print *, 'Test of the 2D ltpic partices (using f0(x,v) = v) '
  call test_poly_ltpic_2d(test_flag)        
  print *, '********************************'
  call test_error_flag( test_flag, test_passed, 'test_poly_ltpic_2d' )


  if ( test_passed ) then
     print *, ' '
     print *, 'ltpic_2d unit test: PASSED'
     print *, ' '
  else
     print *, ' '
     print *, 'ltpic_2d unit test: FAILED'
     print *, ' '
  endif
  
contains

  subroutine test_error_flag( individual_flag, general_flag, message )
    logical, intent(in) :: individual_flag 
    logical, intent(inout) :: general_flag
    character(len=*) :: message
    ! print *, individual_flag
    general_flag = general_flag .and. individual_flag
    if( individual_flag .eqv. .false. ) then
       print *, 'FAILURE IN FUNCTION: ', message
    end if
  end subroutine test_error_flag


  ! testing the ltpic particles with a polynomial function (and stationary solution)
  subroutine test_poly_ltpic_2d( test_passed )
    logical, intent(out)                              :: test_passed
    
    type(lin_trans_pic_2d),             pointer       :: particles
    type(init_poly_2d),                 target        :: init_poly
    class(scalar_field_2d_initializer_base), pointer  :: p_init_f    
    type(init_poly_2d),                 target        :: init_poly_test
    class(scalar_field_2d_initializer_base), pointer  :: p_init_f_test    
    type(sll_mapped_mesh_2d_cartesian), target        :: mesh
    class(sll_mapped_mesh_2d_base),     pointer       :: mesh_base
    type(sll_mapped_mesh_2d_cartesian), target        :: mesh_test
    class(sll_mapped_mesh_2d_base),     pointer       :: mesh_test_base
    type(leap_frog_1st_flow_2d),        target        :: lf_1st_flow
    type(leap_frog_2nd_flow_2d),        target        :: lf_2nd_flow
    class(flow_base_class),             pointer       :: flow
    type(poisson_1d_periodic)                         :: poisson
    sll_real64, dimension(:),           allocatable   :: rho
    sll_real64, dimension(:,:),         allocatable   :: f0_coefs_xv
    sll_real64, dimension(:,:),         allocatable   :: f_h_data
    sll_real64, dimension(:,:),         allocatable   :: f_ref_data
    sll_int32                                         :: ierr
    sll_int32                                         :: f0_degree      
    sll_real64                                        :: final_time
    sll_real64                                        :: elementary_charge
    sll_int32                                         :: num_particles_x
    sll_int32                                         :: num_particles_v
    sll_real64                                        :: xmin 
    sll_real64                                        :: xmax 
    sll_real64                                        :: vmin 
    sll_real64                                        :: vmax 
    sll_int32                                         :: bspline_degree
    sll_int32                                         :: nc_poisson_mesh 
    sll_real64                                        :: xmin_poisson_mesh 
    sll_real64                                        :: xmax_poisson_mesh 
    sll_int32                                         :: num_iterations 
    sll_real64                                        :: time_before_remap 
    sll_real64                                        :: dt
    sll_int32                                         :: n_time_last_remap
    sll_int32                                         :: n_time
    sll_real64                                        :: xmin_test
    sll_real64                                        :: xmax_test
    sll_int32                                         :: npx_test
    sll_real64                                        :: vmin_test
    sll_real64                                        :: vmax_test
    sll_int32                                         :: npv_test
    

    test_passed = .true.

    ! testcase parameters: defining f0(x,v) = v
    f0_degree = 2
    SLL_ALLOCATE( f0_coefs_xv(f0_degree+1,f0_degree+1), ierr )
    f0_coefs_xv = 0
    f0_coefs_xv(1,2) = 1.

    elementary_charge = 1.
    
    xmin_test = 0.4
    xmax_test = 0.6
    npx_test  = 50
    vmin_test = 0.4
    vmax_test = 0.6
    npv_test  = 50
        
    ! numerical parameters (particles)
    num_particles_x = 100
    num_particles_v = 100
    xmin = 0.0
    xmax = 1.0
    vmin = 0.0
    vmax = 1.0
    bspline_degree = 3

    ! numerical parameters (poisson mesh)
    nc_poisson_mesh = 50 ! number of cells
    xmin_poisson_mesh = xmin
    xmax_poisson_mesh = xmax
    
    ! numerical parameters (time resolution)
    final_time = 10
    num_iterations = 10
    time_before_remap = 5 
    dt = final_time/num_iterations

    print *, 'initializing the particles...'
    particles =>      new_ltpic_2d( &
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
    call mesh%initialize(         &
          "map_cartesian",        &
          particles%qi_grid_xmin, &
          particles%qi_grid_xmax, &
          particles%qi_grid_npx,  &
          particles%qi_grid_vmin, &
          particles%qi_grid_vmax, &
          particles%qi_grid_npv    )         
    mesh_base => mesh
    call init_poly%initialize(mesh_base,NODE_CENTERED_FIELD,f0_degree,f0_coefs_xv)
    p_init_f => init_poly
        
    print *, 'initializing the leap-frog flows...'
    call init_lf_1st_flow( lf_1st_flow, dt )
    ! MCP -- here, the array describing the electric field coefficients is a member of the lf_2nd_flow object. 
    ! but I think it should rather have a specific type -- and be a member of some poisson solver class ? 
    ! (then passed as an argument to the flow initializer?)
    call init_lf_2nd_flow( lf_2nd_flow, dt, nc_poisson_mesh, xmin_poisson_mesh, xmax_poisson_mesh, PERIODIC_E_FLOW )

    print *, 'initializing the poisson solver...'
    SLL_ALLOCATE(rho(nc_poisson_mesh+1),  ierr)
    call new(poisson, xmin_poisson_mesh, xmax_poisson_mesh, nc_poisson_mesh,  ierr) 


    print *, 'preparing the mesh and ref data for the tests...'
    call mesh_test%initialize(  &
              "map_cartesian",  &
              xmin_test,        &
              xmax_test,        &
              npx_test,         &
              vmin_test,        &
              vmax_test,        &
              npv_test          )
    mesh_test_base => mesh_test
    call init_poly_test%initialize(mesh_test_base,NODE_CENTERED_FIELD,f0_degree,f0_coefs_xv)
    p_init_f_test => init_poly_test
    
    SLL_ALLOCATE( f_h_data  (npx_test,npv_test),   ierr )
    SLL_ALLOCATE( f_ref_data(npx_test,npv_test),   ierr )

    call p_init_f_test%f_of_x1x2(f_ref_data)


    ! ------------------------------------------------------------------------
    !
    !                           COMPUTE INITIAL DATA
    !
    ! ------------------------------------------------------------------------

    print *, 'loading the particles...'
    call load_ltpic_2d( p_init_f,particles )
    
    n_time_last_remap = 0
    n_time = 0
    
    call measure_error( xmin_test,xmax_test,npx_test,vmin_test,vmax_test,npv_test,f_h_data,f_ref_data,particles,test_passed )      

    ! ------------------------------------------------------------------------
    !
    !                                MAIN LOOP
    !
    ! ------------------------------------------------------------------------

    print *, 'starting time loop : '
    print *, 'final_time        = ', final_time
    print *, 'time_before_remap = ', time_before_remap

    do n_time=1, num_iterations
       print *, 'time step ', n_time, '/', num_iterations

       print *, 'transporting the particles (1st substep in leap-frog scheme)...'
       flow => lf_1st_flow
       call transport_ltpic_2d( flow,particles )
       
!       print *, 'computing the E field for the intermediate solution...'
!       call deposit_charges( rho,particles )
!       call solve( poisson, lf_2nd_flow%electric_field, rho )
       print *, 'setting E = 0 (by hand)...'
       lf_2nd_flow%electric_field = 0

       print *, 'transporting the particles (2nd substep in leap-frog scheme)...'
       flow => lf_2nd_flow
       call transport_ltpic_2d( flow,particles )

      
      call measure_error( xmin_test,xmax_test,npx_test,vmin_test,vmax_test,npv_test,f_h_data,f_ref_data,particles,test_passed )      
      
      ! conditional remapping
      if( n_time < num_iterations .and. (n_time - n_time_last_remap)*dt > time_before_remap ) then
         print *, 'remapping the particles...'
         call remap_ltpic_2d( particles )
         n_time_last_remap = n_time
      end if
    end do ! main loop

    call measure_error( xmin_test,xmax_test,npx_test,vmin_test,vmax_test,npv_test,f_h_data,f_ref_data,particles,test_passed )      

    print *, 'done.'
  end subroutine test_poly_ltpic_2d
  

  ! measure_error
  subroutine measure_error( xmin_test,xmax_test,npx_test,vmin_test,vmax_test,npv_test,f_h_data,f_ref_data,ltpic_object,test_passed )
    sll_real64                                      :: xmin_test
    sll_real64                                      :: xmax_test
    sll_int32                                       :: npx_test
    sll_real64                                      :: vmin_test
    sll_real64                                      :: vmax_test
    sll_int32                                       :: npv_test
    sll_real64, dimension(:,:),     intent(inout)   :: f_h_data
    sll_real64, dimension(:,:),     intent(inout)   :: f_ref_data
    type(lin_trans_pic_2d),         pointer         :: ltpic_object
    logical,                        intent(inout)   :: test_passed

    sll_real64                                        :: error
    sll_int32                                         :: i
    sll_int32                                         :: j
    
    print *, 'measuring the error...'
    call write_f_on_grid (                                                                                &
                          xmin_test,xmax_test,npx_test,                                                   &
                          vmin_test,vmax_test,npv_test,                                                   &
                          f_h_data, ltpic_object                                                          &
                         )
    error = 0
    do j = 1, npv_test
      do i = 1, npx_test
        error = max( error, abs(f_h_data(i,j)-f_ref_data(i,j)) )
      end do
    end do    

    print *, 'measured error (L_inf) = ', error
    if( error .ge. 1.0e-15 ) then
       test_passed = .false.
       print *, 'test_poly_ltpic_2d(): TEST FAILED'
    end if

  end subroutine measure_error

end program ltpic_tester
