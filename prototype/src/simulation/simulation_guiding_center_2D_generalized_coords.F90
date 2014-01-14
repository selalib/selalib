module sll_simulation_2d_guiding_center_generalized_coords_module

#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_field_2d.h"
#include "sll_utilities.h"

  use sll_constants
  use sll_cubic_spline_interpolator_2d
  use sll_cubic_spline_interpolator_1d
  use sll_module_interpolators_2d_base
  use sll_arbitrary_degree_spline_interpolator_2d_module
  use sll_simulation_base
  use sll_logical_meshes
  use sll_coordinate_transformation_2d_base_module
  use sll_general_coordinate_elliptic_solver_module
  use sll_module_scalar_field_2d_base
  use sll_module_scalar_field_2d_alternative
  implicit none



  abstract interface
     function sll_scalar_initializer_2d( x1, x2, params )
       use sll_working_precision
       sll_real64                                     :: sll_scalar_initializer_2d
       sll_real64, intent(in)                         :: x1
       sll_real64, intent(in)                         :: x2
       sll_real64, dimension(:), intent(in), optional :: params
     end function sll_scalar_initializer_2d
   end interface



  type, extends(sll_simulation_base_class) :: sll_simulation_2d_guiding_center_generalized
     
     ! Physics/numerical parameters
     sll_real64 :: dt
     sll_int32  :: num_iterations
     sll_int32  :: carac_case
     sll_int32  :: time_scheme
     sll_int32  :: visu_step
     ! Mesh parameters
     sll_int32  :: nc_x1
     sll_int32  :: nc_x2
     ! for QNS spline_degre in each direction
     sll_int32  :: spline_degree_eta1
     sll_int32  :: spline_degree_eta2
     ! for QNS boundary conditions
     sll_int32  :: bc_left
     sll_int32  :: bc_right
     sll_int32  :: bc_bottom
     sll_int32  :: bc_top
     ! the logical meshes are split in two one for space, one for velocity
     type(sll_logical_mesh_2d), pointer    :: mesh2d
     ! This simulation only applies a coordinate transformation to the spatial
     ! coordinates.
     class(sll_coordinate_transformation_2d_base), pointer :: transf
     type(general_coordinate_elliptic_solver), pointer           :: qns
     
     sll_real64, dimension(:,:), pointer     :: rho_n
     sll_real64, dimension(:,:), pointer     :: rho_np1
     sll_real64, dimension(:,:), pointer     :: rho_nm1
     
     type(arb_deg_2d_interpolator)     :: interp_rho
     type(arb_deg_2d_interpolator)     :: interp_phi
      
     ! for distribution function initializer:
     procedure(sll_scalar_initializer_2d), nopass, pointer :: init_func
     sll_real64, dimension(:), pointer :: params
     
     
     ! for general coordinate QNS, analytical fields
     procedure(two_var_parametrizable_function),nopass,pointer :: a11_f
     procedure(two_var_parametrizable_function),nopass,pointer :: a12_f
     procedure(two_var_parametrizable_function),nopass,pointer :: a21_f
     procedure(two_var_parametrizable_function),nopass,pointer :: a22_f
     procedure(two_var_parametrizable_function),nopass,pointer :: c_f
     sll_real64, dimension(:), pointer :: params_field
   contains
     procedure, pass(sim) :: run => run_2d_gc_general
     procedure, pass(sim) :: init_from_file => init_2d_gc_general
  end type sll_simulation_2d_guiding_center_generalized

  interface delete
     module procedure delete_2d_gc_general
  end interface delete

  interface initialize
     module procedure initialize_2d_gc_general
  end interface initialize

contains

  ! Tentative function to initialize the simulation object 'manually'.
  subroutine initialize_2d_gc_general( &
   sim, &
   mesh2d, &
   transformation, &
   init_func, &
   params,&
   params_field,&
   a11_f,&
   a12_f,&
   a21_f,&
   a22_f,&
   c_f,&
   spline_degre1,&
   spline_degre2,&
   bc_left,&
   bc_right,&
   bc_bottom,&
   bc_top)

   type(sll_simulation_2d_guiding_center_generalized), intent(inout)     :: sim
   type(sll_logical_mesh_2d), pointer                    :: mesh2d
   class(sll_coordinate_transformation_2d_base), pointer :: transformation
   procedure(sll_scalar_initializer_2d)                  :: init_func !! see it 
   sll_real64, dimension(:),target                      :: params
   sll_real64, dimension(:),target                      :: params_field
   procedure(two_var_parametrizable_function) :: a11_f
   procedure(two_var_parametrizable_function) :: a12_f
   procedure(two_var_parametrizable_function) :: a21_f
   procedure(two_var_parametrizable_function) :: a22_f
   procedure(two_var_parametrizable_function) :: c_f
   sll_int32  :: spline_degre1
   sll_int32  :: spline_degre2
   sll_int32  :: bc_left
   sll_int32  :: bc_right
   sll_int32  :: bc_bottom
   sll_int32  :: bc_top
   sll_int32  :: err


   !SLL_ALLOCATE(sim%params(size(params)),err)
   !SLL_ALLOCATE(sim%params_field(size(params_field)),err)
   
   sim%params    => params
   sim%params_field  => params_field
   
   sim%mesh2d  => mesh2d
   sim%transf  => transformation
   sim%init_func => init_func
   
   
   sim%a11_f     => a11_f
   sim%a12_f     => a12_f
   sim%a21_f     => a21_f
   sim%a22_f     => a22_f
   sim%c_f       => c_f
   sim%spline_degree_eta1 = spline_degre1
   sim%spline_degree_eta2 = spline_degre2

   sim%bc_left   = bc_left
   sim%bc_right  = bc_right
   sim%bc_bottom = bc_bottom
   sim%bc_top    = bc_top

   call sim%interp_phi%initialize( &
        sim%mesh2d%num_cells1 +1 , &
        sim%mesh2d%num_cells2 +1, &
        sim%mesh2d%eta1_min, &
        sim%mesh2d%eta1_max, &
        sim%mesh2d%eta2_min, &
        sim%mesh2d%eta2_max, &
        sim%bc_left, &
        sim%bc_right, &
        sim%bc_bottom, &
        sim%bc_top, &
        sim%spline_degree_eta1, &
        sim%spline_degree_eta2)

   call sim%interp_rho%initialize( &
        sim%mesh2d%num_cells1 +1, &
        sim%mesh2d%num_cells2 +1, &
        sim%mesh2d%eta1_min, &
        sim%mesh2d%eta1_max, &
        sim%mesh2d%eta2_min, &
        sim%mesh2d%eta2_max, &
        sim%bc_left, &
        sim%bc_right, &
        sim%bc_bottom, &
        sim%bc_top, &
        sim%spline_degree_eta1, &
        sim%spline_degree_eta2)
        
  end subroutine initialize_2d_gc_general


  subroutine init_2d_gc_general( sim, filename )
    intrinsic :: trim
    class(sll_simulation_2d_guiding_center_generalized), intent(inout) :: sim
    character(len=*), intent(in)                                   :: filename
    sll_int32             :: IO_stat
    sll_real64            :: dt
    sll_int32             :: number_iterations
    sll_int32             :: time_scheme
    sll_int32             :: carac_case
    sll_int32             :: visu_step
    sll_int32             :: num_cells_x1
    sll_int32             :: num_cells_x2
    sll_int32, parameter  :: input_file = 99

    namelist /sim_params/ dt, number_iterations,visu_step,time_scheme,carac_case
    namelist /grid_dims/ num_cells_x1, num_cells_x2
    ! Try to add here other parameters to initialize the mesh values like
    ! xmin, xmax and also for the distribution function initializer.
    open(unit = input_file, file=trim(filename),IOStat=IO_stat)
    if( IO_stat /= 0 ) then
       print *, 'init_par_cart() failed to open file ', filename
       STOP
    end if
    read(input_file, sim_params)
    read(input_file,grid_dims)
    close(input_file)

    sim%dt = dt
    sim%num_iterations = number_iterations
    sim%carac_case = carac_case
    sim%time_scheme = time_scheme
    sim%visu_step = visu_step
    ! In this particular simulation, since the system is periodic, the number
    ! of points is the same as the number of cells in all directions.
    sim%nc_x1 = num_cells_x1
    sim%nc_x2 = num_cells_x2

  end subroutine init_2d_gc_general

  ! Note that the following function has no local variables, which is silly...
  ! This just happened since the guts of the unit test were transplanted here
  ! directly, but this should be cleaned up.
  subroutine run_2d_gc_general(sim)
    class(sll_simulation_2d_guiding_center_generalized), intent(inout) :: sim
    sll_int32  :: i
    sll_int32  :: j
    sll_real64 :: delta1
    sll_real64 :: delta2
    sll_int32  :: itemp
    sll_int32  :: ierr
    sll_int32  :: itime
    sll_int32  :: nc_x1
    sll_int32  :: nc_x2  
    sll_real64 :: eta1
    sll_real64 :: eta2 
    sll_real64 :: eta1_min
    sll_real64 :: eta2_min  
    sll_real64 :: eta1_max
    sll_real64 :: eta2_max
    sll_real64 :: x
    sll_real64 :: y
    sll_real64 :: val_intg_L1
    sll_real64 :: val_intg_L2
    sll_real64 :: val_intg_Linf
    sll_real64 :: val_mass
    sll_real64 :: val_energy
    sll_real64 :: node_val,ref_val
    ! The following could probably be abstracted for convenience
    sll_int32 :: iplot
    character(len=4) :: cplot
    class(sll_scalar_field_2d_base), pointer              :: a11_field_mat
    class(sll_scalar_field_2d_base), pointer              :: a12_field_mat
    class(sll_scalar_field_2d_base), pointer              :: a21_field_mat
    class(sll_scalar_field_2d_base), pointer              :: a22_field_mat
    class(sll_scalar_field_2d_base), pointer              :: b1_field
    class(sll_scalar_field_2d_base), pointer              :: b2_field
    class(sll_scalar_field_2d_base), pointer              :: c_field
    class(sll_scalar_field_2d_discrete_alt), pointer      :: rho_n_ptr
    class(sll_scalar_field_2d_discrete_alt), pointer      :: rho_np1_ptr
    class(sll_scalar_field_2d_discrete_alt), pointer      :: rho_nm1_ptr
    class(sll_scalar_field_2d_discrete_alt), pointer      :: rho_tmp_ptr
    type(sll_scalar_field_2d_discrete_alt), pointer       :: phi
    sll_real64, dimension(:,:), allocatable :: phi_values

    
    ! Start with the fields  
    a11_field_mat => new_scalar_field_2d_analytic_alt( &
         sim%a11_f, &
         "a11", &
         sim%transf, &
         sim%bc_left, &
         sim%bc_right, &
         sim%bc_bottom, &
         sim%bc_top, &
         sim%params_field) 
       

    a12_field_mat => new_scalar_field_2d_analytic_alt( &
         sim%a12_f, &
         "a12", &
         sim%transf, &
         sim%bc_left, &
         sim%bc_right, &
         sim%bc_bottom, &
         sim%bc_top, &
         sim%params_field) 
    
    
    a21_field_mat => new_scalar_field_2d_analytic_alt( &
         sim%a21_f, &
         "a21", &
         sim%transf, &
         sim%bc_left, &
         sim%bc_right, &
         sim%bc_bottom, &
         sim%bc_top, &
         sim%params_field)
       
    
    a22_field_mat => new_scalar_field_2d_analytic_alt( &
         sim%a22_f, &
         "a22", &
         sim%transf, &
         sim%bc_left, &
         sim%bc_right, &
         sim%bc_bottom, &
         sim%bc_top, &
         sim%params_field) 
     

    b1_field => new_scalar_field_2d_analytic_alt( &
          sim%c_f, &
         "c_field", &
         sim%transf, &
         sim%bc_left, &
         sim%bc_right, &
         sim%bc_bottom, &
         sim%bc_top, &
         sim%params_field)

    b2_field => new_scalar_field_2d_analytic_alt( &
          sim%c_f, &
         "c_field", &
         sim%transf, &
         sim%bc_left, &
         sim%bc_right, &
         sim%bc_bottom, &
         sim%bc_top, &
         sim%params_field)
       
    !print*,'pass 1'
       
    !print*,'pass 1'
    c_field => new_scalar_field_2d_analytic_alt( &
          sim%c_f, &
         "c_field", &
         sim%transf, &
         sim%bc_left, &
         sim%bc_right, &
         sim%bc_bottom, &
         sim%bc_top, &
         sim%params_field)
       
    !print*,'pass 1'

    SLL_ALLOCATE(phi_values(sim%mesh2d%num_cells1+1,sim%mesh2d%num_cells2+1),ierr)
    
    phi_values(:,:) = 0.0_f64
 
    print*,'phi => new_scalar_field_2d_discrete_alt ' 
    phi => new_scalar_field_2d_discrete_alt( &
         "phi_check", &
         sim%interp_phi, &
         sim%transf, &
         sim%bc_left, &
         sim%bc_right, &
         sim%bc_bottom, &
         sim%bc_top)
    call phi%set_field_data(phi_values)
    call phi%update_interpolation_coefficients( )
    !print*,'pass 3'
    
    nc_x1 = sim%mesh2d%num_cells1
    nc_x2 = sim%mesh2d%num_cells2
  
    delta1 = sim%mesh2d%delta_eta1
    delta2 = sim%mesh2d%delta_eta2
  
    eta1_min = sim%mesh2d%eta1_min
    eta2_min = sim%mesh2d%eta2_min
   
    eta1_max = sim%mesh2d%eta1_max
    eta2_max = sim%mesh2d%eta2_max
   
    
    SLL_ALLOCATE(sim%rho_n(nc_x1+1,nc_x2+1),ierr) 
    SLL_ALLOCATE(sim%rho_np1(nc_x1+1,nc_x2+1),ierr)
    SLL_ALLOCATE(sim%rho_nm1(nc_x1+1,nc_x2+1),ierr)
    
    sim%rho_n = 0._f64
    sim%rho_np1 = 0._f64
    sim%rho_nm1 = 0._f64
  
    ! this only works because there is no transformation applied in the
    ! velocity space...
    
    
     do j=1,nc_x2+1
        eta2=eta2_min+real(j-1,f64)*delta2
        do i=1,nc_x1+1
          eta1=eta1_min+real(i-1,f64)*delta1
          x = sim%transf%x1(eta1,eta2)
          y = sim%transf%x2(eta1,eta2)
          sim%rho_n(i,j) =  sim%init_func(x,y,sim%params) !0.001*cos(2*sll_pi*eta1)
        end do
     end do
     
   rho_n_ptr => new_scalar_field_2d_discrete_alt( &
         "rho_n", &
         sim%interp_rho, &     
         sim%transf, &
         sim%bc_left, &
         sim%bc_right, &
         sim%bc_bottom, &
         sim%bc_top)     
   call rho_n_ptr%set_field_data(sim%rho_n)
   call rho_n_ptr%update_interpolation_coefficients( )
    
    rho_np1_ptr => new_scalar_field_2d_discrete_alt( &
         "rho_np1", &
         sim%interp_rho, &     
         sim%transf, &
         sim%bc_left, &
         sim%bc_right, &
         sim%bc_bottom, &
         sim%bc_top)                
    call rho_np1_ptr%set_field_data(sim%rho_n)
    call rho_np1_ptr%update_interpolation_coefficients( )
   rho_nm1_ptr => new_scalar_field_2d_discrete_alt( &
         "rho_nm1", &
         sim%interp_rho, &     
         sim%transf, &
         sim%bc_left, &
         sim%bc_right, &
         sim%bc_bottom, &
         sim%bc_top)
    call rho_nm1_ptr%set_field_data(sim%rho_n)     
     call rho_nm1_ptr%update_interpolation_coefficients( )    
    ! Initialize the poisson plan before going into the main loop.
    print *,'Initialize the poisson plan before going into the main loop.'
    sim%qns => new_general_elliptic_solver( &
         sim%spline_degree_eta1, & 
         sim%spline_degree_eta2, & 
         sim%mesh2d%num_cells1, &
         sim%mesh2d%num_cells2, &
         ES_GAUSS_LEGENDRE, &  ! put in arguments
         ES_GAUSS_LEGENDRE, &  ! put in arguments
         sim%bc_left, &
         sim%bc_right, &
         sim%bc_bottom, &
         sim%bc_top, &
         sim%mesh2d%eta1_min, &  
         sim%mesh2d%eta1_max, & 
         sim%mesh2d%eta2_min, & 
         sim%mesh2d%eta2_max ) 
 
    
     ! compute matrix the field
    print *,'Compute matrix the field'
    call factorize_mat_es(&
       sim%qns, &
       a11_field_mat, &
       a12_field_mat,&
       a21_field_mat,&
       a22_field_mat,&
       b1_field, &
       b2_field, &
       c_field)
       
  
    print *, 'started solve_quasi_neutral_eq_general_coords before loop ...'
    call solve_general_coordinates_elliptic_eq( &
            sim%qns, & 
            rho_n_ptr, &
            phi )
          
    do i=1,sim%mesh2d%num_cells1+1
     do j=1,sim%mesh2d%num_cells2+1
     
        eta1       = real(i-1,f64)*delta1 + eta1_min
        eta2       = real(j-1,f64)*delta2 + eta2_min
        node_val   = phi%value_at_point(eta1,eta2)
        ref_val    = -0.001/((2*sll_pi)**2)*cos(2*sll_pi*eta1)
      write(32,*) eta1,eta2,node_val,ref_val,abs(ref_val-node_val)
     enddo
    enddo
    

             
    call calcul_integral(rho_n_ptr,phi,&
       sim%spline_degree_eta1, &
       val_intg_L1,&
       val_intg_L2,&
       val_intg_Linf,&
       val_mass,&
       val_energy)
                        
     write(30,*) 0.*sim%dt, val_intg_L1,val_intg_L2,val_intg_Linf,val_mass,val_energy  
    
    print*, ' ... finished initialization, entering main loop.'
    
   
    ! ------------------------------------------------------------------------
    !
    !                                MAIN LOOP
    !
    ! ------------------------------------------------------------------------
    
    do itime=1,sim%num_iterations

          print *, 'Starting iteration ', itime, ' of ', sim%num_iterations
   
         
       
      select case (sim%time_scheme)

      case(1) 
            
            !Classical semi-Lagrangian scheme (order 1)
             ! compute matrix the field
            
            call solve_general_coordinates_elliptic_eq( &
            sim%qns, & 
            rho_n_ptr, &
            phi )
            !print *, 'advection started ...'
            call advect_CG_curvilinear(rho_n_ptr,sim%rho_np1,phi,&
            sim%dt,&
            sim%carac_case,&
            sim%bc_left, &
            sim%bc_bottom)
            
            call rho_n_ptr%set_field_data(sim%rho_np1)
            call rho_n_ptr%update_interpolation_coefficients( )
            

   case(2) 
            !!'Semi-Lagrangian predictor-corrector scheme'  
            call solve_general_coordinates_elliptic_eq( &
            sim%qns, & 
            rho_n_ptr, &
            phi )
       
            call advect_CG_curvilinear(rho_n_ptr,sim%rho_np1,phi,&
            sim%dt/2_f64,&
            sim%carac_case,&
            sim%bc_left, &
            sim%bc_bottom)
            
            call rho_np1_ptr%set_field_data(sim%rho_np1)
            call rho_np1_ptr%update_interpolation_coefficients( )
            !!we just obtained f^(n+1/2)    
            
            call solve_general_coordinates_elliptic_eq( &
            sim%qns, & 
            rho_np1_ptr, &
            phi )
       
            call advect_CG_curvilinear(rho_n_ptr,sim%rho_np1,phi,&
            sim%dt,&
            sim%carac_case,&
            sim%bc_left, &
            sim%bc_bottom)      
            
            call rho_n_ptr%set_field_data(sim%rho_np1)
            call rho_n_ptr%update_interpolation_coefficients( )    

      case(3)
           !Leap-frog scheme
             if (itime==1) then
                call solve_general_coordinates_elliptic_eq( &
                sim%qns, & 
                rho_n_ptr, &
                phi )
       
                call advect_CG_curvilinear(rho_n_ptr,sim%rho_np1,phi,&
                sim%dt/2_f64,&
                sim%carac_case,&
                sim%bc_left, &
                sim%bc_bottom)
            
                call rho_np1_ptr%set_field_data(sim%rho_np1)
                call rho_np1_ptr%update_interpolation_coefficients( )
             !!we just obtained f^(n+1/2)    
            
                call solve_general_coordinates_elliptic_eq( &
                sim%qns, & 
                rho_np1_ptr, &
                phi )
       
                call advect_CG_curvilinear(rho_n_ptr,sim%rho_np1,phi,&
                sim%dt,&
                sim%carac_case,&
                sim%bc_left, &
                sim%bc_bottom)       
                
             else 
             
               call solve_general_coordinates_elliptic_eq( &
                sim%qns, & 
                rho_n_ptr, &
                phi )
       
               call advect_CG_curvilinear(rho_nm1_ptr,sim%rho_np1,phi,&
               sim%dt*2_f64,&
               sim%carac_case,&
               sim%bc_left, &
               sim%bc_bottom)
               
             end if
                
               call rho_nm1_ptr%set_field_data(sim%rho_n)
               call rho_nm1_ptr%update_interpolation_coefficients( )
               sim%rho_n=sim%rho_np1
               call rho_n_ptr%set_field_data(sim%rho_np1)
               call rho_n_ptr%update_interpolation_coefficients( )
       case default
           print*,'#no scheme defined'
    end select
    
    
  
    call calcul_integral(rho_n_ptr,phi,&
       sim%spline_degree_eta1, &
       val_intg_L1,&
       val_intg_L2,&
       val_intg_Linf,&
       val_mass,&
       val_energy)
                        
     write(30,*) itime*sim%dt, val_intg_L1,val_intg_L2,val_intg_Linf,val_mass,val_energy
     
    
#ifndef NOHDF5
     if (itime==1 .or. ((itime/sim%visu_step)*sim%visu_step==itime)) then
     call plot_f1(rho_n_ptr,sim,itime)
     endif
#endif
     
    end do ! main loop
    
  call rho_n_ptr%delete()
  call rho_np1_ptr%delete()
  call rho_nm1_ptr%delete()
  call c_field%delete()
  call phi%delete()
  call a11_field_mat%delete()
  call a12_field_mat%delete()
  call a21_field_mat%delete()
  call a22_field_mat%delete()
  SLL_DEALLOCATE_ARRAY(phi_values, ierr)
  
  end subroutine run_2d_gc_general
  

  subroutine delete_2d_gc_general( sim )
    type(sll_simulation_2d_guiding_center_generalized) :: sim
    sll_int32 :: ierr
    SLL_DEALLOCATE( sim%rho_n, ierr )
    SLL_DEALLOCATE( sim%rho_np1, ierr )
    SLL_DEALLOCATE( sim%rho_nm1, ierr )
  end subroutine delete_2d_gc_general


 subroutine advect_CG_curvilinear(rho_n,rho_np1,phi,&
    dt,&
    carac_case,&
    bc1_type,bc2_type)

    implicit none
    
    sll_real64, dimension(:,:), pointer ::  rho_np1
    class(sll_scalar_field_2d_discrete_alt), pointer     :: rho_n
    type(sll_scalar_field_2d_discrete_alt) , pointer      :: phi
    class(sll_coordinate_transformation_2d_base), pointer :: T
    type(sll_logical_mesh_2D), pointer :: M
    sll_real64 :: eta1_loc,eta2_loc,eta1,eta1n,eta2,eta20,eta2n,tolr
    sll_real64 :: a_eta1,a_eta2,eta10,eta2_min,eta2_max 
    sll_real64 :: dt, delta_eta1, delta_eta2, eta1_min, eta1_max
    sll_int32  :: N_eta1, N_eta2,carac_case
    sll_int32  :: i,j,maxiter,iter,k_eta1,k_eta2,ii,jj
    sll_int32,intent(in),optional :: bc1_type,bc2_type
    sll_real64 :: phi_loc(2,2)
  
   
    T => phi%get_transformation() 
    M => phi%get_logical_mesh()
    
    N_eta1 = M%num_cells1
    N_eta2 = M%num_cells2
    eta1_min = M%eta1_min
    eta2_min = M%eta2_min
    eta1_max = M%eta1_max
    eta2_max = M%eta2_max
    delta_eta1 = M%delta_eta1
    delta_eta2 = M%delta_eta2
    !print*,N_eta1 ,N_eta2,eta1_min,eta1_max,eta2_min,eta2_max
    
    rho_np1=0._f64
    
    if (carac_case==1) then
       !explicit Euler with linear interpolation  
       do j=1,N_eta2+1
          do i=1,N_eta1+1
             eta1=eta1_min+real(i-1,f64)*delta_eta1
             eta2=eta2_min+real(j-1,f64)*delta_eta2
             
             eta2=eta2+dt*phi%first_deriv_eta1_value_at_point(eta1,eta2)/T%jacobian(eta1,eta2)
             eta1=eta1-dt*phi%first_deriv_eta2_value_at_point(eta1,eta2)/T%jacobian(eta1,eta2)
            
             call correction_BC(bc1_type,bc2_type,eta1_min,eta1_max,eta2_min,eta2_max,&
                &eta1,eta2)
                eta1_loc=(eta1-eta1_min)/(eta1_max-eta1_min)
                eta1_loc=eta1_loc*real(N_eta1,f64)
                k_eta1=floor(eta1_loc)+1
                eta1_loc=eta1_loc-real(k_eta1-1,f64)
                if(((k_eta1-1).gt.(N_eta1)).or.((k_eta1-1).lt.0))then
                  print *,"#bad value of k_eta1=",k_eta1,N_eta1,eta1_loc,eta1,i,j,iter
                endif
                if((k_eta1-1)==N_eta1)then
                  k_eta1=N_eta1
                  if (abs(eta1_loc)>1.e-13) print *,'#eta1_loc=',eta1_loc
                  eta1_loc=1._f64
                endif
                
                eta2_loc=(eta2-eta2_min)/(eta2_max-eta2_min)
                eta2_loc=eta2_loc*real(N_eta2,f64)
                k_eta2=floor(eta2_loc)+1
                eta2_loc=eta2_loc-real(k_eta2-1,f64)
                if(((k_eta2-1).gt.(N_eta2)).or.((k_eta2-1).lt.0))then
                  print *,"#bad value of k_eta2=",k_eta2,N_eta2
                endif
                if((k_eta2-1)==N_eta2)then
                  k_eta2=N_eta2
                  if (abs(eta2_loc)>1.e-13) print *,'#eta2_loc=',eta2_loc
                  eta2_loc=1._f64
                endif

             
             rho_np1(i,j)=(1.0_f64-eta2_loc)*((1.0_f64-eta1_loc)*rho_n%values(k_eta1,k_eta2) &
             &+eta1_loc*rho_n%values(k_eta1+1,k_eta2)) +eta2_loc*((1.0_f64-eta1_loc)* &
             & rho_n%values(k_eta1,k_eta2+1)+eta1_loc*rho_n%values(k_eta1+1,k_eta2+1))
            
          end do

       end do
  end if

  if (carac_case==2) then
      !explicit Euler with "Spline interpolation"
       do j=1,N_eta2+1
          do i=1,N_eta1+1
            
             eta10=eta1_min+real(i-1,f64)*delta_eta1
             eta20=eta2_min+real(j-1,f64)*delta_eta2
             eta2=eta20+(dt*phi%first_deriv_eta1_value_at_point(eta10,eta20)) !/T%jacobian(eta10,eta20))
             eta1=eta10-(dt*phi%first_deriv_eta2_value_at_point(eta10,eta20)) !/T%jacobian(eta10,eta20))

             call correction_BC(bc1_type,bc2_type,eta1_min,eta1_max,eta2_min,eta2_max,&
                &eta1,eta2)
                
             rho_np1(i,j)=rho_n%value_at_point(eta1,eta2)

          end do
       end do
       
  end if

    
    if (carac_case==3) then
       !using fixed point method
    
       !initialization
       maxiter=40
       tolr=1e-10
       eta1n=0.0_f64
       eta2n=0.0_f64

       do j=1,N_eta2+1 !N_eta2+1
          eta20=eta2_min+real(j-1,f64)*delta_eta2
          eta2=eta20
          k_eta2=j
          do i=1,N_eta1+1 !N_eta1+1
             eta10=eta1_min+real(i-1,f64)*delta_eta1
             eta1=eta10
             eta1_loc=0.0_f64
             eta2_loc=0.0_f64
             k_eta1=i
             a_eta1=0.0_f64
             a_eta2=0.0_f64
             iter=0
 
           do while (((iter<maxiter) .and. (abs((eta1n-eta1))+abs((eta2n-eta2))>tolr)).or.(iter==0))    
      
                eta1_loc=(eta1-eta1_min)/(eta1_max-eta1_min)
                eta1_loc=eta1_loc*real(N_eta1,f64)
                k_eta1=floor(eta1_loc)+1
                eta1_loc=eta1_loc-real(k_eta1-1,f64)
                if(((k_eta1-1).gt.(N_eta1)).or.((k_eta1-1).lt.0))then
                  print *,"#bad value of k_eta1=",k_eta1,N_eta1,eta1_loc,eta1,i,j,iter
                endif
                if((k_eta1-1)==N_eta1)then
                  k_eta1= N_eta1
                  if (abs(eta1_loc)>1.e-13) print *,'#eta1_loc=',eta1_loc
                  eta1_loc=1._f64 
                endif  

             
                eta2_loc=(eta2-eta2_min)/(eta2_max-eta2_min)
                eta2_loc=eta2_loc*real(N_eta2,f64)
                k_eta2=floor(eta2_loc)+1
                eta2_loc=eta2_loc-real(k_eta2-1,f64)
                if(((k_eta2-1).gt.(N_eta2)).or.((k_eta2-1).lt.0))then
                  print *,"#bad value of k_eta2=",k_eta2,N_eta2
                endif
                if((k_eta2-1)==N_eta2)then
                  k_eta2=N_eta2
                  if (abs(eta2_loc)>1.e-13) print *,'#eta2_loc=',eta2_loc
                  eta2_loc=1._f64 
                endif
             
               do jj=0,1
                 do ii=0,1               
                 phi_loc(1+ii,1+jj) = phi%first_deriv_eta2_value_at_indices(k_eta1+ii,k_eta2+jj)/T%jacobian_at_node(k_eta1+ii,k_eta2+jj)
                 enddo
               enddo
               phi_loc(1,1) = (1.0_f64-eta1_loc)*phi_loc(1,1) + eta1_loc*phi_loc(2,1)
               phi_loc(1,2) = (1.0_f64-eta1_loc)*phi_loc(1,2) + eta1_loc*phi_loc(2,2)
               a_eta1 = (1.0_f64-eta2_loc)*phi_loc(1,1) + eta2_loc*phi_loc(1,2)
               a_eta1 = 0.5_f64*dt*a_eta1

               do jj=0,1
                 do ii=0,1               
                   phi_loc(1+ii,1+jj) = -phi%first_deriv_eta1_value_at_indices(k_eta1+ii,k_eta2+jj)/T%jacobian_at_node(k_eta1+ii,k_eta2+jj)
                 enddo
               enddo
               phi_loc(1,1) = (1.0_f64-eta1_loc)*phi_loc(1,1) + eta1_loc*phi_loc(2,1)
               phi_loc(1,2) = (1.0_f64-eta1_loc)*phi_loc(1,2) + eta1_loc*phi_loc(2,2)
               a_eta2 = (1.0_f64-eta2_loc)*phi_loc(1,1) + eta2_loc*phi_loc(1,2)
               a_eta2 = 0.5_f64*dt*a_eta2

                eta1n=eta1
                eta2n=eta2
                eta1=eta10-a_eta1
                eta2=eta20-a_eta2              
                call correction_BC(bc1_type,bc2_type,eta1_min,eta1_max,& 
                                 & eta2_min,eta2_max,eta1,eta2)                           
                iter=iter+1
                
             end do

             if (iter==maxiter .and. abs((eta1n-eta1))+abs((eta2n-eta2))>tolr) then
                print*,'#no convergence in fixed point method',i,j !,step
             end if
             eta1=eta10-2.0_f64*a_eta1
             eta2=eta20-2.0_f64*a_eta2
            call correction_BC(bc1_type,bc2_type,eta1_min,eta1_max, & 
                             & eta2_min,eta2_max,eta1,eta2)   
            rho_np1(i,j)=rho_n%value_at_point(eta1,eta2)

          end do

       end do
  
  end if




    if(bc2_type==sll_PERIODIC) rho_np1(:,N_eta2+1)=rho_np1(:,1) 
    if(bc1_type==sll_PERIODIC) rho_np1(N_eta1+1,:)=rho_np1(1,:)
 

  end subroutine advect_CG_curvilinear

!****************************************************************************

subroutine correction_BC(bc1_type,bc2_type,eta1_min,eta1_max,eta2_min,eta2_max,eta1,eta2)

   implicit none

    sll_real64 :: eta1_min,eta1_max,eta2_min,eta2_max,eta1,eta2  
    sll_int32  :: bc1_type,bc2_type
     
  
! --- Corrections on the BC ---
        if (bc1_type.eq.SLL_DIRICHLET) then
          eta1 = min(max(eta1,eta1_min),eta1_max)
        endif
        if (bc2_type.eq.SLL_DIRICHLET) then
          eta2 = min(max(eta2,eta2_min),eta2_max)
        endif
        if (bc1_type==SLL_PERIODIC) then
          do while (eta1>eta1_max)
            eta1 = eta1-(eta1_max-eta1_min)
          enddo
          do while (eta1<eta1_min)
            eta1 = eta1+(eta1_max-eta1_min)
          enddo
        endif
        if (bc2_type==SLL_PERIODIC) then
          do while (eta2>eta2_max)
            eta2 = eta2-(eta2_max-eta2_min)
          enddo
          do while (eta2<eta2_min)
            eta2 = eta2+(eta2_max-eta2_min)
          enddo
        endif
  end subroutine correction_BC 

!****************************************************************************

subroutine calcul_integral(rho_n,phi,&
       spline_degree, &
       val_intg_L1,&
       val_intg_L2,&
       val_intg_Linf,&
       val_mass,&
       val_energy)
    implicit none
    
    ! IMPUT VARIABLES
    class(sll_scalar_field_2d_discrete_alt), pointer      :: rho_n
    type(sll_scalar_field_2d_discrete_alt), pointer       :: phi
    integer :: spline_degree
    
    !OUTPUT VARIABLES
     sll_real64, intent(out)  :: val_intg_L1,val_mass,val_energy
     sll_real64, intent(out)  :: val_intg_L2
     sll_real64, intent(out)  :: val_intg_Linf
    !LOCAL VARIABLES
    integer :: ig,jg, mx,my
    sll_real64 :: ptgaussx, ptgaussy, wgaussptx,wgausspty
    sll_real64 :: val1,val2
    sll_real64 :: valjac,et1i,et2i,valf,yg,valEx,valEy
    sll_real64 :: valdx_1,valdx_2,valdy_1,valdy_2
    sll_real64 :: delta_eta1, delta_eta2, eta1_min, eta1_max
    sll_real64:: eta2_min, eta2_max
    sll_real64, dimension(spline_degree+2) :: gauss_x,gauss_w
    integer :: N_eta1,N_eta2,size_ptgauss
    
    
    class(sll_coordinate_transformation_2d_base), pointer :: T
    type(sll_logical_mesh_2D), pointer :: M
    sll_real64, dimension(1:2,1:2) :: jac_m
    
  
   
    T => phi%get_transformation() 
    M => phi%get_logical_mesh()
    
    N_eta1 = M%num_cells1
    N_eta2 = M%num_cells2
    eta1_min = M%eta1_min
    eta2_min = M%eta2_min
    eta1_max = M%eta1_max
    eta2_max = M%eta2_max
    delta_eta1 = M%delta_eta1
    delta_eta2 = M%delta_eta2
    
    val_intg_L1 = 0.0_f64
    val_intg_L2 = 0.0_f64
    val_intg_Linf = 0.0_f64 
    val_mass = 0.0_f64
    val_energy=0.0_f64 
   
    !call pointgauss(spline_degree)    
    gauss_x = 0.0
    gauss_w = 0.0
   
      ! set Gauss points and weights            
      select case(spline_degree+1)
      
      case(1) 
         gauss_x(1) = -1.0_8/sqrt(3.0_8)
         gauss_x(2) =  1.0_8/sqrt(3.0_8)
         gauss_w(1) =  1.0_8 
         gauss_w(2) =  1.0_8
      case(2)
         gauss_x(1) = -sqrt(3.0_8/5.0_8)
         gauss_x(2) = 0.0_8 
         gauss_x(3) = sqrt(3.0_8/5.0_8)
         gauss_w(1) = 5.0_8/9.0_8
         gauss_w(2) = 8.0_8/9.0_8
         gauss_w(3) = gauss_w(1)
      case(3)
         gauss_x(4) = sqrt((3.0_8+2.0_8*sqrt(6.0_8/5.0_8))/7.0_8)
         gauss_x(3) = sqrt((3.0_8-2.0_8*sqrt(6.0_8/5.0_8))/7.0_8)
         gauss_x(2) = -gauss_x(3) 
         gauss_x(1) = -gauss_x(4) 
         gauss_w(1) = (18.0_8-sqrt(30.0_8))/36.0_8
         gauss_w(2) = (18.0_8+sqrt(30.0_8))/36.0_8  
         gauss_w(3) = gauss_w(2)  
         gauss_w(4) = gauss_w(1)
      case(4)
         gauss_x(1) = -0.90617984593866374_8
         gauss_x(2) = -0.53846931010568311_8
         gauss_x(3) =  0._8
         gauss_x(4) =  0.53846931010568311_8
         gauss_x(5) =  0.90617984593866374_8
         gauss_w(1) =  0.23692688505618875_8  
         gauss_w(2) =  0.47862867049936653_8
         gauss_w(3) =  0.568888888888888888_8
         gauss_w(4) =  0.47862867049936653_8
         gauss_w(5) =  0.23692688505618875_8
      case default
         print*, 'spline degree ', spline_degree, ' not implemented'
         stop
      end select
      
    size_ptgauss = size(gauss_x)

    do my=1,N_eta2

       et2i  = eta2_min + (my-1)*delta_eta2
       
       do ig = 1, size_ptgauss
          
          ptgaussy = et2i  + 0.5_8 * delta_eta2*(gauss_x(ig)+1.0_8)
          wgausspty= 0.5_8 * delta_eta2*gauss_w(ig)
          
          do mx = 1, N_eta1

             et1i  = eta1_min + (mx-1)*delta_eta1
             do jg=1,size_ptgauss 
                
                ptgaussx=  0.5_8*delta_eta1*(gauss_x(jg)+1)+ et1i
                wgaussptx= 0.5_8*delta_eta1*gauss_w(jg)
                
                jac_m =  T%jacobian_matrix(ptgaussx,ptgaussy)
                valdx_1 = jac_m(1,1)    !dx/deta1
                valdx_2 = jac_m(1,2)    !dx/deta2
                valdy_1 = jac_m(2,1)    !dy/deta1
                valdy_2 = jac_m(2,2)    !dy/deta2
                
                valjac = 1. !T%jacobian(ptgaussx,ptgaussy)
        
                valf = rho_n%value_at_point(ptgaussx,ptgaussy)
               
                valEx = phi%first_deriv_eta1_value_at_point(ptgaussx,ptgaussy)
                valEy = phi%first_deriv_eta2_value_at_point(ptgaussx,ptgaussy)
                        
                val_mass = val_mass + &
                               valf*wgausspty*wgaussptx*abs(valjac)
                val_intg_L1 = val_intg_L1 + &
                               abs(valf)*wgausspty*wgaussptx*abs(valjac)
                val_intg_L2 = val_intg_L2 + &
                               (valf)**2*wgausspty*wgaussptx*abs(valjac)
                val_intg_Linf= max(val_intg_Linf,abs(valf))
                
                val_energy = val_energy + &
                               ((valdy_2*valEx -valdy_1*valEy )**2 + (-valdx_2*valEx +valdx_1*valEy )**2) &
                               *wgausspty*wgaussptx/abs(valjac)      
               
             end do
          end do
          
       end do
    end do
   val_intg_L2=sqrt(val_intg_L2)  
   val_energy =sqrt(val_energy)
  end subroutine calcul_integral

!****************************************************************************

#ifndef NOHDF5
subroutine plot_f1(rho,sim,itime)!

  use sll_xdmf
  use sll_hdf5_io

  sll_int32 :: file_id, hfile_id
  sll_int32 :: error
  class(sll_simulation_2d_guiding_center_generalized) :: sim
  !sll_real64, dimension (:,:), intent(in):: rho
  class(sll_scalar_field_2d_discrete_alt), pointer      :: rho
  sll_real64, dimension(:,:),allocatable :: f
  sll_real64, dimension(:,:),allocatable :: x1
  sll_real64, dimension(:,:),allocatable :: x2
  sll_int32 :: i, j
  sll_int32, intent(in) :: itime
  character(len=4)      :: cplot
  sll_int32             :: nnodes_x1, nnodes_x2
  

  nnodes_x1 = sim%nc_x1+1
  nnodes_x2 = sim%nc_x2+1
  call int2string(itime,cplot)
   

  if (itime == 1) then
     SLL_ALLOCATE(x1(nnodes_x1,nnodes_x2),error)
     SLL_ALLOCATE(x2(nnodes_x1,nnodes_x2),error)
     do j=1,nnodes_x2
      do i=1,nnodes_x1
        x1(i,j) = sim%transf%x1_at_node(i,j)
        x2(i,j) = sim%transf%x2_at_node(i,j)
      enddo
     enddo   
     call sll_hdf5_file_create("curvilinear_mesh-x1.h5",hfile_id,error)
     call sll_hdf5_write_array(hfile_id,x1,"/x1",error)
     call sll_hdf5_file_close(hfile_id, error)
     call sll_hdf5_file_create("curvilinear_mesh-x2.h5",hfile_id,error)
     call sll_hdf5_write_array(hfile_id,x2,"/x2",error)
     call sll_hdf5_file_close(hfile_id, error)
     deallocate(x1)
     deallocate(x2)
  end if
  
  SLL_ALLOCATE(f(nnodes_x1,nnodes_x2),error)
  do j=1,nnodes_x2
    do i=1,nnodes_x1
      f(i,j) = rho%value_at_indices(i,j)
    enddo
  enddo  
   
  call int2string(itime,cplot)
  call sll_xdmf_open("rho"//cplot//".xmf","curvilinear_mesh",nnodes_x1,nnodes_x2,file_id,error)
  call sll_xdmf_write_array("rho"//cplot,f,"values",error,file_id,"Node")
  call sll_xdmf_close(file_id,error)!

  deallocate(f)
 end subroutine plot_f1

#endif
end module sll_simulation_2d_guiding_center_generalized_coords_module
