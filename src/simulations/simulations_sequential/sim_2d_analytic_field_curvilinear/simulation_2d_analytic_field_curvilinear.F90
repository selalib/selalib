module sll_simulation_2d_analytic_field_curvilinear_module

!2d guiding center cartesian simulation
!contact: Adnane Hamiaz (hamiaz@math.unistra.fr
!         Michel Mehrenberger (mehrenbe@math.unistra.fr)

#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_fields.h"
#include "sll_utilities.h"
#include "sll_poisson_solvers.h"
  use sll_constants
  use sll_cartesian_meshes  
  use sll_module_advection_1d_periodic
  use sll_module_advection_2d_BSL
  use sll_module_advection_2d_tensor_product
  use sll_module_characteristics_2d_explicit_euler
  use sll_module_characteristics_2d_verlet
  use sll_module_advection_1d_BSL
  use sll_module_advection_1d_CSL
  use sll_module_advection_1d_PSM
  use sll_module_characteristics_1d_explicit_euler
  use sll_module_characteristics_1d_trapezoid
  use sll_module_characteristics_1d_explicit_euler_conservative
  use sll_reduction_module
  use sll_simulation_base
  use sll_module_cubic_spline_interpolator_2d
  use sll_module_cubic_spline_interpolator_1d
  use sll_coordinate_transformation_2d_base_module
  use sll_module_coordinate_transformations_2d
  use sll_common_coordinate_transformations
  use sll_common_array_initializers_module
  !use sll_mudpack_curvilinear
#ifdef MUDPACK
  use sll_module_poisson_2d_mudpack_curvilinear_solver_old
#endif
!  use sll_module_poisson_2d_elliptic_solver
!  use sll_timer
!  use sll_fft
!  use sll_module_poisson_2d_periodic_solver
  use sll_parallel_array_initializer_module
  
  implicit none
  
  
  sll_int32, parameter :: SLL_EULER = 0 
  sll_int32, parameter :: SLL_PREDICTOR_CORRECTOR = 1 
  sll_int32, parameter :: SLL_LEAP_FROG = 2 
  sll_int32, parameter :: SLL_COMPUTE_FIELD_FROM_PHI = 0 
  sll_int32, parameter :: SLL_COMPUTE_FIELD_FROM_ANALYTIC = 1 

  type, extends(sll_simulation_base_class) :: &
    sll_simulation_2d_analytic_field_curvilinear

   !geometry
   type(sll_cartesian_mesh_2d), pointer :: mesh_2d
  
   !transformation 
    class(sll_coordinate_transformation_2d_base), pointer :: transformation
  
   !initial function
   procedure(sll_scalar_initializer_2d), nopass, pointer :: init_func
   sll_real64, dimension(:), pointer :: params
   !set boundary functions for characteristics
   procedure(signature_process_outside_point), nopass, pointer :: process_outside_point1_func
   procedure(signature_process_outside_point), nopass, pointer :: process_outside_point2_func
   
      
    !advector
   class(sll_advection_2d_base), pointer    :: advect_2d
   procedure(sll_scalar_initializer_2d), nopass, pointer :: phi_func
   procedure(sll_scalar_initializer_2d), nopass, pointer :: A1_func
   procedure(sll_scalar_initializer_2d), nopass, pointer :: A2_func
   procedure(sll_scalar_initializer_4d), nopass, pointer :: A1_exact_charac_func
   procedure(sll_scalar_initializer_4d), nopass, pointer :: A2_exact_charac_func
   sll_real64, dimension(:), pointer :: A_func_params
   procedure(sll_scalar_initializer_1d), nopass, pointer :: A_time_func
   sll_real64, dimension(:), pointer :: A_time_func_params
   
   
   !interpolator for derivatives
   class(sll_interpolator_2d_base), pointer   :: phi_interp2d

    
   !time_iterations
   sll_real64 :: dt
   sll_int32  :: num_iterations
   sll_int32  :: freq_diag
   sll_int32  :: freq_diag_time

   !time_loop
   sll_int32 :: time_loop_case
   
   sll_int32 :: compute_field_case
     
   !boundaries conditions 
   sll_int32  :: bc_eta1_left
   sll_int32  :: bc_eta1_right
   sll_int32  :: bc_eta2_left
   sll_int32  :: bc_eta2_right
   sll_int32  :: bc_interp2d_eta1
   sll_int32  :: bc_interp2d_eta2   
   sll_int32  :: bc_charac2d_eta1
   sll_int32  :: bc_charac2d_eta2 
  contains
    procedure, pass(sim) :: run => run_af2d_curvilinear
    procedure, pass(sim) :: init_from_file => init_fake
     
  end type sll_simulation_2d_analytic_field_curvilinear





contains

  function new_analytic_field_2d_curvilinear(filename) result(sim)
    type(sll_simulation_2d_analytic_field_curvilinear), pointer :: sim 
    character(len=*), intent(in), optional :: filename   
    sll_int32 :: ierr
    
    SLL_ALLOCATE(sim,ierr)
    
    call initialize_analytic_field_2d_curvilinear(sim,filename)
    
  
  
  end function new_analytic_field_2d_curvilinear
  
  subroutine initialize_analytic_field_2d_curvilinear(sim,filename)
    class(sll_simulation_2d_analytic_field_curvilinear), intent(inout) :: sim
    character(len=*), intent(in), optional :: filename
    sll_int32             :: IO_stat
    sll_int32, parameter  :: input_file = 99
        
    !geometry
    character(len=256) :: mesh_case_eta1
    character(len=256) :: mesh_case_eta2
    character(len=256) :: transf_case
    sll_int32 :: num_cells_eta1
    sll_int32 :: num_cells_eta2
    sll_real64 :: eta1_min
    sll_real64 :: eta1_max
    sll_real64 :: eta2_min
    sll_real64 :: eta2_max
    sll_real64 :: alpha1
    sll_real64 :: alpha2
    
     !initial_function
    character(len=256) :: initial_function_case
    sll_real64 :: kmode_eta1
    sll_real64 :: kmode_eta2
    sll_real64 :: eps
    sll_real64 :: xc_1
    sll_real64 :: xc_2
    sll_real64 :: sigma_1
    sll_real64 :: sigma_2

    
     !time_iterations
    sll_real64 :: dt
    sll_int32 :: number_iterations
    sll_int32 :: freq_diag
    sll_int32 :: freq_diag_time
    character(len=256) :: time_loop_case


    !advector
    character(len=256) :: advect2d_case 
    character(len=256) :: f_interp2d_case
    character(len=256) :: phi_interp2d_case
    character(len=256) ::  charac2d_case
    character(len=256) ::  A_interp_case
    character(len=256) ::  advection_field_case
    character(len=256) :: advect1d_x1_case 
    character(len=256) :: advect1d_x2_case 
    character(len=256) :: charac1d_x1_case
    character(len=256) :: charac1d_x2_case
    character(len=256) :: f_interp1d_x1_case
    character(len=256) :: f_interp1d_x2_case
    character(len=256) ::  compute_field_case
    sll_real64 :: time_period

 
    !boundaries conditions
    !character(len=256) :: boundaries_conditions
    character(len=256) :: bc_interp2d_eta1
    character(len=256) :: bc_interp2d_eta2   
    character(len=256) :: bc_charac2d_eta1
    character(len=256) :: bc_charac2d_eta2  
    character(len=256) :: bc_eta1_left
    character(len=256) :: bc_eta1_right
    character(len=256) :: bc_eta2_left
    character(len=256) :: bc_eta2_right  
    
    !local variables
    sll_int32 :: Nc_eta1
    sll_int32 :: Nc_eta2
    type(sll_cartesian_mesh_1d), pointer :: mesh_x1
    type(sll_cartesian_mesh_1d), pointer :: mesh_x2
    class(sll_interpolator_2d_base), pointer :: f_interp2d
    class(sll_interpolator_2d_base), pointer :: phi_interp2d
    class(sll_characteristics_2d_base), pointer :: charac2d
    class(sll_characteristics_1d_base), pointer :: charac1d_x1
    class(sll_characteristics_1d_base), pointer :: charac1d_x2
    class(sll_interpolator_2d_base), pointer   :: A1_interp2d
    class(sll_interpolator_2d_base), pointer   :: A2_interp2d
    class(sll_interpolator_1d_base), pointer   :: A1_interp1d_x1
    class(sll_interpolator_1d_base), pointer   :: A2_interp1d_x1
    class(sll_interpolator_1d_base), pointer   :: A2_interp1d_x2
    class(sll_interpolator_1d_base), pointer :: f_interp1d_x1
    class(sll_interpolator_1d_base), pointer :: f_interp1d_x2
    class(sll_advection_1d_base), pointer    :: advect_1d_x1
    class(sll_advection_1d_base), pointer    :: advect_1d_x2
    sll_int32 :: ierr
    sll_real64, dimension(4) :: params_mesh
    sll_real64 :: eta1_min_bis
    sll_real64 :: eta1_max_bis
    sll_real64 :: eta2_min_bis
    sll_real64 :: eta2_max_bis
    sll_int32  :: Nc_eta1_bis
    sll_int32  :: Nc_eta2_bis
    sll_real64 :: A1 !parameter for translation
    sll_real64 :: A2 !parameter for translation

    !here we do all the initialization
    !in future, we will use namelist file
    
    namelist /geometry/ &
      mesh_case_eta1, &
      mesh_case_eta2, &
      num_cells_eta1, &
      num_cells_eta2, &
      eta1_min, &
      eta1_max, &
      eta2_min, &
      eta2_max, &
      transf_case,&
      alpha1  , &
      alpha2

     namelist /initial_function/ &
      initial_function_case, &
      kmode_eta1, &
      kmode_eta2, &
      eps, &
      xc_1, &
      xc_2, &
      sigma_1, &
      sigma_2
      
    namelist /time_iterations/ &
      dt, &
      number_iterations, &
      freq_diag, &
      freq_diag_time, &
      time_loop_case

     namelist /advector/ &
      advect2d_case, &   
      f_interp2d_case, &
      phi_interp2d_case, &
      charac2d_case, &
      A_interp_case, &
      advection_field_case, &
      charac1d_x1_case, &
      charac1d_x2_case, &
      advect1d_x1_case, &   
      advect1d_x2_case, & 
      time_period, &
      A1, &
      A2, &
      compute_field_case
   
      
     namelist /boundaries/ &
      bc_interp2d_eta1, &
      bc_interp2d_eta2, &    
      bc_charac2d_eta1, &
      bc_charac2d_eta2, &
      bc_eta1_left,  &
      bc_eta1_right, &
      bc_eta2_left,  &
      bc_eta2_right   
    
        !! set default parameters
    
    !geometry
    transf_case = "SLL_CARTESIAN"
    mesh_case_eta1="SLL_CARTESIAN_MESH"
    mesh_case_eta2="SLL_CARTESIAN_MESH"
    num_cells_eta1 = 128
    num_cells_eta2 = 128
    eta1_min = 0.0_f64
    eta1_max = 2._f64*sll_pi
    eta2_min = 0.0_f64
    eta2_max = 2._f64*sll_pi
    alpha1 = 0._f64
    alpha2 = 0._f64
    params_mesh = (/ alpha1, alpha2, eta1_max-eta1_min, eta2_max-eta2_min/)
    
    !initial function
    initial_function_case="SLL_KHP1"
    kmode_eta1 = 0.5_f64
    kmode_eta2 = 1._f64
    eps = 0.015_f64
    
    !time_iterations
    dt = 0.1_f64
    number_iterations  = 600
    freq_diag = 100
    freq_diag_time = 1
    !time_loop_case = "SLL_EULER"
    time_loop_case = "SLL_PREDICTOR_CORRECTOR" 

    !advector
    advect2d_case = "SLL_BSL"    
    f_interp2d_case = "SLL_CUBIC_SPLINES"
    phi_interp2d_case = "SLL_CUBIC_SPLINES"
    !charac2d_case = "SLL_EULER"
    charac2d_case = "SLL_VERLET"
    A_interp_case = "SLL_CUBIC_SPLINES"
    advection_field_case = "SLL_SWIRLING_DEFORMATION_FLOW"
    !advection_field_case = "SLL_TRANSLATION_FLOW"
    advect1d_x1_case = "SLL_BSL"
    advect1d_x2_case = "SLL_BSL"
    charac1d_x1_case = "SLL_EULER"
    !charac1d_x1_case = "SLL_TRAPEZOID"
    charac1d_x2_case = "SLL_EULER"
    !charac1d_x2_case = "SLL_TRAPEZOID"
    f_interp1d_x1_case = "SLL_CUBIC_SPLINES"
    f_interp1d_x2_case = "SLL_CUBIC_SPLINES"
    A1 = 1._f64
    A2 = 1._f64
    compute_field_case = "SLL_COMPUTE_FIELD_FROM_ANALYTIC"
    
    !boundaries conditions
    sim%bc_eta1_left = SLL_PERIODIC
    sim%bc_eta1_right= SLL_PERIODIC
    sim%bc_eta2_left = SLL_PERIODIC
    sim%bc_eta2_right= SLL_PERIODIC 
    sim%bc_interp2d_eta1 = SLL_PERIODIC
    sim%bc_interp2d_eta2 = SLL_PERIODIC
    sim%bc_charac2d_eta1 = SLL_PERIODIC
    sim%bc_charac2d_eta2 = SLL_PERIODIC

    if(present(filename))then
      open(unit = input_file, file=trim(filename)//'.nml',IOStat=IO_stat)
        if( IO_stat /= 0 ) then
          print *, '#in initialize_analytic_field_2d_curvilinear() failed to open file ', &
          trim(filename)//'.nml'
          STOP
        end if
      print *,'#initialization with filename:'
      print *,'#',trim(filename)//'.nml'
      read(input_file, geometry) 
      read(input_file, initial_function)
      read(input_file, time_iterations)
      read(input_file, advector)
      read(input_file, boundaries)
      close(input_file)
    else
      print *,'#initialization with default parameters'    
    endif

 

    Nc_eta1 =  num_cells_eta1
    Nc_eta2 =  num_cells_eta2
    sim%dt = dt
    sim%num_iterations = number_iterations
    sim%freq_diag = freq_diag 
    sim%freq_diag_time = freq_diag_time
   
    select case(bc_eta1_left)
      case ("SLL_PERIODIC")
        print*,"#bc_eta1_left = SLL_PERIODIC" 
        sim%bc_eta1_left = SLL_PERIODIC
      case ("SLL_DIRICHLET")
        print*,"#bc_eta1_left = SLL_DIRICHLET"  
        sim%bc_eta1_left = SLL_DIRICHLET
      case default
        print *,'#bad bc_eta1_left',bc_eta1_left
        print *,'#not implemented'
        print *,'#in initialize_analytic_field_2d_curvilinear'
        stop
    end select
    
    select case(bc_eta1_right)
      case ("SLL_PERIODIC")
        print*,"#bc_eta1_right = SLL_PERIODIC" 
        sim%bc_eta1_right = SLL_PERIODIC
      case ("SLL_DIRICHLET")
        print*,"#bc_eta1_right = SLL_DIRICHLET"  
        sim%bc_eta1_right = SLL_DIRICHLET
      case default
        print *,'#bad bc_eta1_right',bc_eta1_right
        print *,'#not implemented'
        print *,'#in initialize_analytic_field_2d_curvilinear'
        stop
    end select
    
       select case(bc_eta2_left)
      case ("SLL_PERIODIC")
        print*,"#bc_eta2_left = SLL_PERIODIC" 
        sim%bc_eta2_left = SLL_PERIODIC
      case ("SLL_DIRICHLET")
        print*,"#bc_eta2_left = SLL_DIRICHLET"  
        sim%bc_eta2_left = SLL_DIRICHLET
      case default
        print *,'#bad bc_eta2_left',bc_eta2_left
        print *,'#not implemented'
        print *,'#in initialize_analytic_field_2d_curvilinear'
        stop
    end select
    
    select case(bc_eta2_right)
      case ("SLL_PERIODIC")
        print*,"#bc_eta2_right = SLL_PERIODIC" 
        sim%bc_eta2_right = SLL_PERIODIC
      case ("SLL_DIRICHLET")
        print*,"#bc_eta2_right = SLL_DIRICHLET"  
        sim%bc_eta2_right = SLL_DIRICHLET
      case default
        print *,'#bad bc_eta2_right',bc_eta2_right
        print *,'#not implemented'
        print *,'#in initialize_analytic_field_2d_curvilinear'
        stop
    end select
    
    select case(bc_interp2d_eta1)
      case ("SLL_PERIODIC")
        print*,"#bc_interp2d_eta1= SLL_PERIODIC" 
        sim%bc_interp2d_eta1 = SLL_PERIODIC
      case ("SLL_DIRICHLET")
        print*,"#bc_interp2d_eta1 = SLL_DIRICHLET"  
        sim%bc_interp2d_eta1= SLL_DIRICHLET
      case ("SLL_HERMITE")
        print*,"#bc_interp2d_eta1 = SLL_HERMITE"  
        sim%bc_interp2d_eta1= SLL_HERMITE 
      case default
        print *,'#bad bc_interp2d_eta1',bc_interp2d_eta1
        print *,'#not implemented'
        print *,'#in initialize_analytic_field_2d_curvilinear'
        stop
    end select

    select case(bc_interp2d_eta2)
      case ("SLL_PERIODIC")
        print*,"#bc_interp2d_eta2= SLL_PERIODIC" 
        sim%bc_interp2d_eta2 = SLL_PERIODIC
      case ("SLL_DIRICHLET")
        print*,"#bc_interp2d_eta2 = SLL_DIRICHLET"  
        sim%bc_interp2d_eta2= SLL_DIRICHLET
      case ("SLL_HERMITE")
        print*,"#bc_interp2d_eta2 = SLL_HERMITE"  
        sim%bc_interp2d_eta2= SLL_HERMITE 
      case default
        print *,'#bad bc_interp2d_eta2',bc_interp2d_eta2
        print *,'#not implemented'
        print *,'#in initialize_analytic_field_2d_curvilinear'
        stop
    end select
    
    select case(bc_charac2d_eta1)
      case ("SLL_PERIODIC")
        print*,"#bc_charac2d_eta1= SLL_PERIODIC" 
        sim%bc_charac2d_eta1= SLL_PERIODIC
        sim%process_outside_point1_func => process_outside_point_periodic1            
      !case ("SLL_DIRICHLET")
      !  print*,"#bc_charac2d_eta1 = SLL_DIRICHLET"  
      !  sim%bc_charac2d_eta1= SLL_DIRICHLET
      !case ("SLL_HERMITE")
      !  print*,"#bc_charac2d_eta1 = SLL_HERMITE"  
      !  sim%bc_charac2d_eta1= SLL_HERMITE 
      case ("SLL_SET_TO_LIMIT")
        print*,"#bc_charac2d_eta1 = SLL_SET_TO_LIMIT"  
        sim%bc_charac2d_eta1= SLL_SET_TO_LIMIT
        sim%process_outside_point1_func => process_outside_point_set_to_limit1    
      case default
        print *,'#bad bc_charac2d_eta1',bc_charac2d_eta1
        print *,'#not implemented'
        print *,'#in initialize_analytic_field_2d_curvilinear'
        stop
    end select

    select case(bc_charac2d_eta2)
      case ("SLL_PERIODIC")
        print*,"#bc_charac2d_eta2= SLL_PERIODIC" 
        sim%bc_charac2d_eta2= SLL_PERIODIC
        sim%process_outside_point2_func => process_outside_point_periodic1
      !case ("SLL_DIRICHLET")
      !  print*,"#bc_charac2d_eta2= SLL_DIRICHLET"  
      !  sim%bc_charac2d_eta2= SLL_DIRICHLET
      !case ("SLL_HERMITE")
      !  print*,"#bc_charac2d_eta2 = SLL_HERMITE"  
      !  sim%bc_charac2d_eta2= SLL_HERMITE 
      case ("SLL_SET_TO_LIMIT")
        print*,"#bc_charac2d_eta2 = SLL_SET_TO_LIMIT"  
        sim%bc_charac2d_eta2= SLL_SET_TO_LIMIT
        sim%process_outside_point2_func => process_outside_point_set_to_limit1   
      case default
        print *,'#bad bc_charac2d_eta2',bc_charac2d_eta2
        print *,'#not implemented'
        print *,'#in initialize_analytic_field_2d_curvilinear'
        stop
    end select
    
    
    select case (mesh_case_eta1)
      case ("SLL_CARTESIAN_MESH")
        mesh_x1 => new_cartesian_mesh_1d(Nc_eta1,eta_min=eta1_min, eta_max=eta1_max)  
      case default
        print*,'#mesh_case_eta1', mesh_case_eta1, ' not implemented'
        stop 
    end select
    select case (mesh_case_eta2)
      case ("SLL_CARTESIAN_MESH")
        mesh_x2 => new_cartesian_mesh_1d(Nc_eta2,eta_min=eta2_min, eta_max=eta2_max)
      case default
        print*,'#mesh_case_eta2', mesh_case_eta2, ' not implemented'
        stop 
    end select
    sim%mesh_2d =>  mesh_x1 * mesh_x2 ! tensor product
    !  In collela  mesh params_mesh =( alpha1, alpha2, L1, L2 ) such that :
    !  x1= eta1 + alpha1*sin(2*pi*eta1/L1)*sin(2*pi*eta2/L2)
    params_mesh = (/ alpha1, alpha2, eta1_max-eta1_min, eta2_max-eta2_min/)
          
!   sim%mesh_2d => new_cartesian_mesh_2d( &
!      Nc_eta1, &
!      Nc_eta2, &
!      eta1_min , &
!      eta1_max , &
!      eta2_min , &
!      eta2_max ) 
      
    select case (transf_case)
      case ("SLL_CARTESIAN")
        print*,'#transf_case =SLL_CARTESIAN'       
        sim%transformation => new_coordinate_transformation_2d_analytic( &
         "analytic_identity_transformation", &
         sim%mesh_2d, &
         identity_x1, &
         identity_x2, &
         identity_jac11, &
         identity_jac12, &
         identity_jac21, &
         identity_jac22, &
         params_mesh   )  
      case ("SLL_POLAR") 
        print*,'#transf_case =SLL_POLAR'  
        sim%transformation => new_coordinate_transformation_2d_analytic( &
         "analytic_polar_transformation", &
         sim%mesh_2d, &
         polar_x1, &
         polar_x2, &
         polar_jac11, &
         polar_jac12, &
         polar_jac21, &
         polar_jac22, &
         params_mesh  )     
      case ("SLL_COLLELA")
        print*,'#transf_case =SLL_COLLELA'  
        sim%transformation => new_coordinate_transformation_2d_analytic( &
         "analytic_collela_transformation", &
         sim%mesh_2d, &
         sinprod_x1, &
         sinprod_x2, &
         sinprod_jac11, &
         sinprod_jac12, &
         sinprod_jac21, &
         sinprod_jac22, &
         params_mesh  )  
        case default
        print *,'#bad transf_case',transf_case
        print *,'#not implemented'
        print *,'#in initialize_analytic_field_2d_curvilinear'
        stop
    end select  
     
    select case(advect1d_x1_case)
      case ("SLL_BSL")
        eta1_min_bis = eta1_min
        eta1_max_bis = eta1_max
        Nc_eta1_bis = Nc_eta1
      case ("SLL_PSM")
        eta1_min_bis = eta1_min
        eta1_max_bis = eta1_max
        Nc_eta1_bis = Nc_eta1
      case ("SLL_CSL")
        eta1_min_bis = eta1_min-0.5_f64*mesh_x1%delta_eta
        eta1_max_bis = eta1_max-0.5_f64*mesh_x1%delta_eta
        Nc_eta1_bis = Nc_eta1
      case default
        print *,'#bad value of advect1d_x1_case'
        stop  
    end select

    select case(advect1d_x2_case)
      case ("SLL_BSL")
        eta2_min_bis = eta2_min
        eta2_max_bis = eta2_max
        Nc_eta2_bis = Nc_eta2
      case ("SLL_PSM")
        eta2_min_bis = eta2_min
        eta2_max_bis = eta2_max
        Nc_eta2_bis = Nc_eta2
      case ("SLL_CSL")
        eta2_min_bis = eta2_min-0.5_f64*mesh_x2%delta_eta
        eta2_max_bis = eta2_max-0.5_f64*mesh_x2%delta_eta        
        Nc_eta2_bis = Nc_eta2
      case default
        print *,'#bad value of advect1d_x2_case'
        stop  
    end select
  
    select case (f_interp2d_case)
      case ("SLL_CUBIC_SPLINES")
        print*,"#f interpolation SLL_CUBIC_SPLINES"
        f_interp2d => new_cubic_spline_interpolator_2d( &
          Nc_eta1+1, &
          Nc_eta2+1, &
          eta1_min, &
          eta1_max, &
          eta2_min, &
          eta2_max, &
          sim%bc_interp2d_eta1, &
          sim%bc_interp2d_eta2)
      case default
        print *,'#bad f_interp2d_case',f_interp2d_case
        print *,'#not implemented'
        print *,'#in initialize_analytic_field_2d_curvilinear'
        stop
    end select




    select case (A_interp_case)
      case ("SLL_CUBIC_SPLINES")
       print*,"#A1_2d interpolation SLL_CUBIC_SPLINES"
        A1_interp2d => new_cubic_spline_interpolator_2d( &
          Nc_eta1+1, &
          Nc_eta2+1, &
          eta1_min, &
          eta1_max, &
          eta2_min, &
          eta2_max, &
          sim%bc_interp2d_eta1, &
          sim%bc_interp2d_eta2)
       print*,"#A2_2d interpolation SLL_CUBIC_SPLINES"   
        A2_interp2d => new_cubic_spline_interpolator_2d( &
          Nc_eta1+1, &
          Nc_eta2+1, &
          eta1_min, &
          eta1_max, &
          eta2_min, &
          eta2_max, &
          sim%bc_interp2d_eta1, &
          sim%bc_interp2d_eta2)  
       print*,"#A1_1d interpolation SLL_CUBIC_SPLINES"   
        A1_interp1d_x1 => new_cubic_spline_interpolator_1d( &
          Nc_eta1+1, &
          eta1_min, &
          eta1_max, &
          sim%bc_interp2d_eta1)
       print*,"#A2_1d interpolation SLL_CUBIC_SPLINES"     
        A2_interp1d_x1 => new_cubic_spline_interpolator_1d( &
          Nc_eta1+1, &
          eta1_min, &
          eta1_max, &
          sim%bc_interp2d_eta1)
      case default
        print *,'#bad A_interp_case',A_interp_case
        print *,'#not implemented'
        print *,'#in initialize_analytic_field_2d_curvilinear'
        stop
    end select

    select case (phi_interp2d_case)
      case ("SLL_CUBIC_SPLINES")
      print*,"#phi interpolation SLL_CUBIC_SPLINES"  
        phi_interp2d => new_cubic_spline_interpolator_2d( &
          Nc_eta1+1, &
          Nc_eta2+1, &
          eta1_min, &
          eta1_max, &
          eta2_min, &
          eta2_max, &
          sim%bc_interp2d_eta1, &
          sim%bc_interp2d_eta2)         
      case default
        print *,'#bad phi_interp2d_case',phi_interp2d_case
        print *,'#not implemented'
        print *,'#in initialize_analytic_field_2d_curvilinear'
        stop
    end select

    select case (f_interp1d_x1_case)
      case ("SLL_CUBIC_SPLINES")
        f_interp1d_x1 => new_cubic_spline_interpolator_1d( &
          Nc_eta1_bis+1, &
          eta1_min_bis, &
          eta1_max_bis, &
          sim%bc_interp2d_eta1)
      case default
        print *,'#bad f_interp1d_x1_case',f_interp1d_x1_case
        print *,'#not implemented'
        print *,'#in initialize_guiding_center_2d_curvilinear'
        stop
    end select


    select case (f_interp1d_x2_case)
      case ("SLL_CUBIC_SPLINES")
        f_interp1d_x2 => new_cubic_spline_interpolator_1d( &
          Nc_eta2_bis+1, &
          eta2_min_bis, &
          eta2_max_bis, &
          sim%bc_interp2d_eta2)
      case default
        print *,'#bad f_interp1d_x2_case',f_interp1d_x2_case
        print *,'#not implemented'
        print *,'#in initialize_guiding_center_2d_curvilinear'
        stop
    end select


    select case(charac1d_x1_case)
      case ("SLL_EULER")
        charac1d_x1 => new_explicit_euler_1d_charac(&
          Nc_eta1_bis+1, &
          eta_min=eta1_min_bis, &
          eta_max=eta1_max_bis, &
          bc_type= sim%bc_charac2d_eta1)    
      case ("SLL_TRAPEZOID")
        charac1d_x1 => &
          new_trapezoid_1d_charac(&
          Nc_eta1_bis+1, &
          A1_interp1d_x1, &
          bc_type= sim%bc_charac2d_eta1, &
          eta_min=eta1_min_bis, &
          eta_max=eta1_max_bis)
      case ("SLL_EULER_CONSERVATIVE")
        charac1d_x1 => new_explicit_euler_conservative_1d_charac(&
          Nc_eta1_bis+1, &
          eta_min=eta1_min_bis, &
          eta_max=eta1_max_bis, &
          bc_type=sim%bc_charac2d_eta1)    
      case default
        print *,'#bad charac1d_x1_case',charac1d_x1_case
        print *,'#not implemented'
        print *,'#in initialize_guiding_center_2d_curvilinear'
        stop
    end select

    select case(charac1d_x2_case)
      case ("SLL_EULER")
        charac1d_x2 => new_explicit_euler_1d_charac(&
          Nc_eta2_bis+1, &
          eta_min=eta2_min_bis, &
          eta_max=eta2_max_bis, &
          bc_type= sim%bc_charac2d_eta2)    
      case ("SLL_TRAPEZOID")
        charac1d_x2 => &
          new_trapezoid_1d_charac(&
          Nc_eta2_bis+1, &
          A2_interp1d_x2, &
          bc_type= sim%bc_charac2d_eta2, &
          eta_min=eta2_min_bis, &
          eta_max=eta2_max_bis)
      case ("SLL_EULER_CONSERVATIVE")
        charac1d_x2 => new_explicit_euler_conservative_1d_charac(&
          Nc_eta2_bis+1, &
          eta_min=eta2_min_bis, &
          eta_max=eta2_max_bis, &
          bc_type= sim%bc_charac2d_eta2)    
      case default
        print *,'#bad charac1d_x2_case',charac1d_x2_case
        print *,'#not implemented'
        print *,'#in initialize_guiding_center_2d_curvilinear'
        stop
    end select

    select case(advect1d_x1_case)
      case ("SLL_BSL")
        advect_1d_x1 => new_BSL_1d_advector(&
          f_interp1d_x1, &
          charac1d_x1, &
          Nc_eta1_bis+1, &
          eta_min = eta1_min_bis, &
          eta_max = eta1_max_bis)
      case ("SLL_CSL")
        advect_1d_x1 => new_CSL_1d_advector(&
          f_interp1d_x1, &
          charac1d_x1, &
          Nc_eta1_bis+1, &
          eta_min = eta1_min_bis, &
          eta_max = eta1_max_bis)
      case ("SLL_PSM")
        advect_1d_x1 => new_PSM_1d_advector(&
          Nc_eta1+1, &
          eta_min = eta1_min, &
          eta_max = eta1_max)
      case default
        print *,'#bad advect_case',advect1d_x1_case
        print *,'#not implemented'
        print *,'#in initialize_guiding_center_2d_curvilinear'
        stop
    end select

    select case(advect1d_x2_case)
      case ("SLL_BSL")
        advect_1d_x2 => new_BSL_1d_advector(&
          f_interp1d_x2, &
          charac1d_x2, &
          Nc_eta2_bis+1, &
          eta_min = eta2_min_bis, &
          eta_max = eta2_max_bis)
      case ("SLL_CSL")
        advect_1d_x2 => new_CSL_1d_advector(&
          f_interp1d_x2, &
          charac1d_x2, &
          Nc_eta2_bis+1, &
          eta_min = eta2_min_bis, &
          eta_max = eta2_max_bis)
      case ("SLL_PSM")
        advect_1d_x2 => new_PSM_1d_advector(&
          Nc_eta2+1, &
          eta_min = eta2_min, &
          eta_max = eta2_max)
      case default
        print *,'#bad advect_case',advect1d_x2_case
        print *,'#not implemented'
        print *,'#in initialize_guiding_center_2d_curvilinear'
        stop
    end select


    select case(charac2d_case)
      case ("SLL_EULER")
         print*,"#charac = SLL_EULER"  
        charac2d => new_explicit_euler_2d_charac(&
          Nc_eta1+1, &
          Nc_eta2+1, &
          eta1_min=eta1_min, &
          eta1_max=eta1_max, &
          eta2_min=eta2_min, &
          eta2_max=eta2_max, &
          bc_type_1=sim%bc_charac2d_eta1, & !SLL_SET_TO_LIMIT, &
          bc_type_2=sim%bc_charac2d_eta2)    
      case ("SLL_VERLET")    
          print*,"#charac =SLL_VERLET"   
        charac2d => new_verlet_2d_charac(&
          Nc_eta1+1, &
          Nc_eta2+1, &
          A1_interp2d, &
          A2_interp2d, &
          A1_interp1d_x1, &
          A2_interp1d_x1, &
          bc_type_1=sim%bc_charac2d_eta1, & !SLL_SET_TO_LIMIT, &
          bc_type_2=sim%bc_charac2d_eta2, &
          eta1_min=eta1_min, &
          eta1_max=eta1_max, &
          eta2_min=eta2_min, &
          eta2_max=eta2_max)
      case default
        print *,'#bad charac2d_case',charac2d_case
        print *,'#not implemented'
        print *,'#in initialize_analytic_field_2d_curvilinear'
        stop
    end select

  
    sim%phi_interp2d => phi_interp2d

    select case(advect2d_case)
      case ("SLL_BSL")
       print*,"#advect2d = SLL_BSL "  
        sim%advect_2d => new_BSL_2d_advector(&
          f_interp2d, &
          charac2d, &
          Nc_eta1+1, &
          Nc_eta2+1, &
          eta1_min = eta1_min, &
          eta1_max = eta1_max, &
          eta2_min = eta2_min, &
          eta2_max = eta2_max)
      case ("SLL_TENSOR_PRODUCT")
       print*,"#advect2d = SLL_SPLITING " 
        sim%advect_2d => new_tensor_product_2d_advector(&
          advect_1d_x1, &
          advect_1d_x2, &
          Nc_eta1+1, &
          Nc_eta2+1)    
      case default
        print *,'#bad advect_case',advect2d_case
        print *,'#not implemented'
        print *,'#in initialize_analytic_field_2d__curvilinear'
        stop
    end select



    
    select case(advection_field_case)
      case ("SLL_SWIRLING_DEFORMATION_FLOW")
        sim%phi_func => sll_SDF_phi_initializer_2d 
        sim%A1_func => sll_SDF_A1_initializer_2d 
        sim%A2_func => sll_SDF_A2_initializer_2d 
        sim%A1_exact_charac_func => sll_SDF_A1_exact_charac_2d 
        sim%A2_exact_charac_func => sll_SDF_A2_exact_charac_2d 
        SLL_ALLOCATE(sim%A_func_params(2),ierr)
        sim%A_time_func => sll_SDF_time_initializer_1d 
        SLL_ALLOCATE(sim%A_time_func_params(1),ierr)
        sim%A_time_func_params(1) = time_period
      case ("SLL_ROTATION_FLOW")
        sim%phi_func => sll_rotation_phi_initializer_2d 
        sim%A1_func => sll_rotation_A1_initializer_2d 
        sim%A2_func => sll_rotation_A2_initializer_2d 
        sim%A1_exact_charac_func => sll_rotation_A1_exact_charac_2d 
        sim%A2_exact_charac_func => sll_rotation_A2_exact_charac_2d 
        SLL_ALLOCATE(sim%A_func_params(2),ierr)
        sim%A_time_func => sll_constant_time_initializer_1d 
        SLL_ALLOCATE(sim%A_time_func_params(1),ierr)
        sim%A_time_func_params(1) = 1._f64
      case ("SLL_TRANSLATION_FLOW")
        sim%phi_func => sll_translation_phi_initializer_2d 
        sim%A1_func => sll_translation_A1_initializer_2d 
        sim%A2_func => sll_translation_A2_initializer_2d 
        sim%A1_exact_charac_func => sll_translation_A1_exact_charac_2d 
        sim%A2_exact_charac_func => sll_translation_A2_exact_charac_2d 
        SLL_ALLOCATE(sim%A_func_params(2),ierr)
        sim%A_func_params(1) = A1
        sim%A_func_params(2) = A2
        sim%A_time_func => sll_constant_time_initializer_1d 
        SLL_ALLOCATE(sim%A_time_func_params(1),ierr)
        sim%A_time_func_params(1) = 1._f64
      case default
        print *,'#bad advect_case',advection_field_case
        print *,'#not implemented'
        print *,'#in initialize_analytic_field_2d_curvilinear'
        stop
    end select
    
   
    select case(initial_function_case)
      case ("SLL_KHP1")
        print*,"#f0 = SLL_KHP1" 
        sim%init_func => sll_KHP1_2d
        SLL_ALLOCATE(sim%params(3),ierr)
        sim%params(1) = eps
        sim%params(2) = kmode_eta1
        sim%params(3) = kmode_eta2
      case ("SLL_GAUSSIAN")
        print*,"#f0 = SLL_GAUSSIAN " 
        sim%init_func => sll_gaussian_initializer_2d
        SLL_ALLOCATE(sim%params(4),ierr)
        sim%params(1) = xc_1
        sim%params(2) = xc_2
        sim%params(3) = sigma_1
        sim%params(4) = sigma_2
      case ("SLL_COS_BELL")
        print*,"#f0 = SLL__COS_BELL " 
        sim%init_func => sll_cos_bell_initializer_2d
        SLL_ALLOCATE(sim%params(2),ierr)
        sim%params(1) = xc_1
        sim%params(2) = xc_2  
      case default
        print *,'#bad initial_function_case',initial_function_case
        print *,'#not implemented'
        print *,'#in initialize_analytic_field_2d_curvilinear'
        stop
    end select
    
    !time_loop
    select case(time_loop_case)
      case ("SLL_EULER")
        print*,"#time_loop = SLL_EULER " 
        sim%time_loop_case = SLL_EULER
      case ("SLL_PREDICTOR_CORRECTOR")
       print*,"#time_loop = SLL_PREDICTOR_CORRECTOR " 
        sim%time_loop_case = SLL_PREDICTOR_CORRECTOR
      case default
        print *,'#bad time_loop_case',time_loop_case
        print *,'#not implemented'
        print *,'#in initialize_analytic_field_2d_curvilinear'
        stop
    end select


    select case(compute_field_case)
      case ("SLL_COMPUTE_FIELD_FROM_ANALYTIC")
        print*,"#compute_field_case = SLL_COMPUTE_FIELD_FROM_ANALYTIC " 
        sim%compute_field_case = SLL_COMPUTE_FIELD_FROM_ANALYTIC
      case ("SLL_COMPUTE_FIELD_FROM_PHI")
        print*,"#compute_field_case = SLL_COMPUTE_FIELD_FROM_PHI " 
        sim%compute_field_case = SLL_COMPUTE_FIELD_FROM_PHI
      case default
        SLL_ERROR("bad compute_field_case")
        stop
    end select

    
   
  end subroutine initialize_analytic_field_2d_curvilinear
  


  subroutine init_fake(sim, filename)
    class(sll_simulation_2d_analytic_field_curvilinear), intent(inout) :: sim
    character(len=*), intent(in)                                :: filename
  
    print *,'# Do not use the routine init_vp4d_fake'
    print *,'#use instead initialize_vlasov_par_poisson_seq_curv'
    stop
  
  end subroutine init_fake
  
  subroutine run_af2d_curvilinear(sim)
    class(sll_simulation_2d_analytic_field_curvilinear), intent(inout) :: sim
    sll_int32 :: Nc_eta1
    sll_int32 :: Nc_eta2
    sll_real64 :: delta_eta1
    sll_real64 :: delta_eta2
    sll_real64 :: eta1_min,eta1_max
    sll_real64 :: eta2_min,eta2_max 
    sll_real64 :: eta1
    sll_real64 :: eta2  
    sll_real64 :: x1
    sll_real64 :: x2
    sll_int32 :: i1 
    sll_int32 :: i2
    sll_real64,dimension(:,:),  pointer :: f
    sll_real64,dimension(:,:),  pointer :: f_exact
    sll_real64,dimension(:,:),  pointer :: f_old
    sll_real64,dimension(:,:),  pointer :: f_init
    sll_real64,dimension(:,:),  pointer :: phi
    sll_real64,dimension(:,:),  pointer :: A1 !advection fields
    sll_real64,dimension(:,:),  pointer :: A2
    sll_real64,dimension(:,:),  pointer :: A1_init !advection fields
    sll_real64,dimension(:,:),  pointer :: A2_init
    sll_real64,dimension(:,:),  pointer :: div
    sll_int32  :: ierr
    sll_int32  :: nb_step
    sll_int32  :: step
    sll_real64 :: dt
    !sll_int32  :: diag_id = 77
    sll_int32 :: thdiag_id
    sll_int32  :: iplot
    sll_real64 :: time_factor
    sll_real64 :: err
    sll_real64 :: feet1
    sll_real64 :: feet2
    sll_real64 :: jac_m(2,2)

    Nc_eta1 = sim%mesh_2d%num_cells1
    Nc_eta2 = sim%mesh_2d%num_cells2
    delta_eta1 = sim%mesh_2d%delta_eta1
    delta_eta2 = sim%mesh_2d%delta_eta2
    eta1_min = sim%mesh_2d%eta1_min
    eta2_min = sim%mesh_2d%eta2_min
    eta1_max = sim%mesh_2d%eta1_max
    eta2_max = sim%mesh_2d%eta2_max
    nb_step = sim%num_iterations
    dt = sim%dt
    
    
    !allocation
    SLL_ALLOCATE(f(Nc_eta1+1,Nc_eta2+1),ierr)
    SLL_ALLOCATE(f_old(Nc_eta1+1,Nc_eta2+1),ierr)
    SLL_ALLOCATE(f_init(Nc_eta1+1,Nc_eta2+1),ierr)
    SLL_ALLOCATE(f_exact(Nc_eta1+1,Nc_eta2+1),ierr)
    SLL_ALLOCATE(phi(Nc_eta1+1,Nc_eta2+1),ierr)
    SLL_ALLOCATE(A1(Nc_eta1+1,Nc_eta2+1),ierr)
    SLL_ALLOCATE(A2(Nc_eta1+1,Nc_eta2+1),ierr)
    SLL_ALLOCATE(A1_init(Nc_eta1+1,Nc_eta2+1),ierr)
    SLL_ALLOCATE(A2_init(Nc_eta1+1,Nc_eta2+1),ierr)
    SLL_ALLOCATE(div(Nc_eta1+1,Nc_eta2+1),ierr)
    

    
    !initialisation of distribution function    
     do i2=1,Nc_eta2+1
        eta2=eta2_min+real(i2-1,f64)*delta_eta2
        do i1=1,Nc_eta1+1
          eta1=eta1_min+real(i1-1,f64)*delta_eta1
          x1 = sim%transformation%x1(eta1,eta2)
          x2 = sim%transformation%x2(eta1,eta2)
          f(i1,i2) =  sim%init_func(x1,x2,sim%params) 
          f_init(i1,i2)  =  sim%init_func(x1,x2,sim%params)
          phi(i1,i2) = sim%phi_func(x1,x2,sim%A_func_params)
          !warning specific change for colella mesh
          !phi(i1,i2) = phi(i1,i2)-sim%A_func_params(1)*eta2+sim%A_func_params(2)*eta1
          phi(i1,i2) = exp(1._f64/(cos(x1)-1.001_f64))*exp(1._f64/(cos(x2)-1.001_f64))-sim%A_func_params(1)*eta2+sim%A_func_params(2)*eta1
          jac_m  =  sim%transformation%jacobian_matrix(eta1,eta2)          
          call compute_curvilinear_field_2d( &
            sim%A1_func(x1,x2,sim%A_func_params), &
            sim%A2_func(x1,x2,sim%A_func_params), &
            A1_init(i1,i2), &
            A2_init(i1,i2), &
            jac_m, &
            sim%transformation%jacobian(eta1,eta2))
          !A1_init(i1,i2) =  &
          !  sim%A1_func(x1,x2,sim%A_func_params) !/sim%transformation%jacobian(eta1,eta2)
          !A2_init(i1,i2) =  &
          !  sim%A2_func(x1,x2,sim%A_func_params) !/sim%transformation%jacobian(eta1,eta2)
        end do
     end do
    
    if(sim%compute_field_case==SLL_COMPUTE_FIELD_FROM_PHI) then   
      call compute_field_from_phi_2d_curvilinear( &
        phi, &
        sim%mesh_2d, &
        sim%transformation, &
        A1_init, &
        A2_init, &
        sim%phi_interp2d)
      do i2=1,Nc_eta2+1
        eta2=eta2_min+real(i2-1,f64)*delta_eta2
        do i1=1,Nc_eta1+1
          eta1=eta1_min+real(i1-1,f64)*delta_eta1
          A1_init(i1,i2) =  A1_init(i1,i2)+sim%A_func_params(1)/sim%transformation%jacobian(eta1,eta2)
          A2_init(i1,i2) =  A2_init(i1,i2)+sim%A_func_params(2)/sim%transformation%jacobian(eta1,eta2)
        enddo
      enddo    
    endif           

    call sll_ascii_file_create('thdiag.dat', thdiag_id, ierr)
    
    iplot = 0

    err = 0._f64
    f_exact = f

 


    do step=1,nb_step
      print*,"step= ", step
      f_old = f
    




      select case (sim%time_loop_case)
        case (SLL_EULER)
          time_factor = sim%A_time_func( &
            real(step-1,f64)*sim%dt, &
            sim%A_time_func_params )
          A1 = time_factor*A1_init
          A2 = time_factor*A2_init          
          call sim%advect_2d%advect_2d(A1, A2, sim%dt, f_old, f)
        case (SLL_PREDICTOR_CORRECTOR)
          time_factor = sim%A_time_func( &
            (real(step-1,f64)+0.5_f64)*sim%dt, &
            sim%A_time_func_params )
          A1 = time_factor*A1_init
          A2 = time_factor*A2_init          
          call sim%advect_2d%advect_2d(A1, A2, sim%dt, f_old, f)
        case default  
          print *,'#bad time_loop_case',sim%time_loop_case
          print *,'#not implemented'
          print *,'#in run_af2d_curvilinear'
          print *,'#available options are:'
          print *,'#SLL_EULER=',SLL_EULER
          print *,'#SLL_PREDICTOR_CORRECTOR=',SLL_PREDICTOR_CORRECTOR
          
      end select
      
      do i2=1,Nc_eta2+1
        eta2=eta2_min+real(i2-1,f64)*delta_eta2
        do i1=1,Nc_eta1+1
          eta1=eta1_min+real(i1-1,f64)*delta_eta1
          x1 = sim%transformation%x1(eta1,eta2)
          x2 = sim%transformation%x2(eta1,eta2)
          feet1 = sim%A1_exact_charac_func( &
            0._f64, &
            real(step,f64)*sim%dt, &
            x1, &
            x2, &
            sim%A_func_params) 
          feet2 = sim%A2_exact_charac_func( &
            0._f64, &
            real(step,f64)*sim%dt, &
            x1, &
            x2, &
            sim%A_func_params)
          feet1 = sim%process_outside_point1_func( feet1, eta1_min, eta1_max )   
          feet2 = sim%process_outside_point2_func( feet2, eta2_min, eta2_max )   
          f_exact(i1,i2) =  sim%init_func(feet1,feet2,sim%params) 
        end do
      end do
      
      err = max(maxval(abs(f_exact-f)),err)
      
      
      
      call compute_divergence_2d_curvilinear( &
        div, &
        sim%mesh_2d, &
        sim%transformation, &
        A1, &
        A2, &
        sim%phi_interp2d)


#ifndef NOHDF5
      if(modulo(step,sim%freq_diag)==0)then
        call plot_f_curvilinear(iplot,div,sim%mesh_2d,sim%transformation)
        iplot = iplot+1  
      endif            
#endif  


      if(modulo(step,sim%freq_diag_time)==0)then
        call time_history_diagnostic_curvilinear( &
          thdiag_id , &    
          step, &
          dt, &
          sim%mesh_2d, &
          sim%transformation, &
          f, &
          phi, &
          A1, &
          A2)
      endif            
      
      
         
    enddo
    
    close(thdiag_id)   
    print *,err    
    print *,'#run_af2d_curvilinear PASSED'
    
  end subroutine run_af2d_curvilinear
  
  

#ifndef NOHDF5
!*********************
!*********************

  !---------------------------------------------------
  ! Save the mesh structure
  !---------------------------------------------------
  subroutine plot_f_curvilinear(iplot,f,mesh_2d,transf)
    use sll_xdmf
    use sll_hdf5_io_serial
    sll_int32 :: file_id
    sll_int32 :: error
    sll_real64, dimension(:,:), allocatable :: x1
    sll_real64, dimension(:,:), allocatable :: x2
    sll_int32 :: i, j
    sll_int32, intent(in) :: iplot
    character(len=4)      :: cplot
    sll_int32             :: nnodes_x1, nnodes_x2
    type(sll_cartesian_mesh_2d), pointer :: mesh_2d
    class(sll_coordinate_transformation_2d_base), pointer :: transf
    sll_real64, dimension(:,:), intent(in) :: f
    sll_real64 :: eta1
    sll_real64 :: eta2
    sll_real64 ::  eta1_min, eta2_min
    sll_real64 ::  eta1_max, eta2_max  
    sll_real64 :: deta1
    sll_real64 :: deta2
    
    
    nnodes_x1 = mesh_2d%num_cells1+1
    nnodes_x2 = mesh_2d%num_cells2+1
    eta1_min = mesh_2d%eta1_min
    eta1_max = mesh_2d%eta1_max
    eta2_min = mesh_2d%eta2_min
    eta2_max = mesh_2d%eta2_max
    deta1 = mesh_2d%delta_eta1
    deta2 = mesh_2d%delta_eta2
    
    !print *,'#maxf=',iplot,maxval(f),minval(f)
    

    
    if (iplot == 1) then

      SLL_ALLOCATE(x1(nnodes_x1,nnodes_x2), error)
      SLL_ALLOCATE(x2(nnodes_x1,nnodes_x2), error)
      do j = 1,nnodes_x2
        do i = 1,nnodes_x1
          eta1 = eta1_min+real(i-1,f32)*deta1
          eta2 = eta2_min+real(j-1,f32)*deta2
          x1(i,j) = transf%x1(eta1,eta2)
          x2(i,j) = transf%x2(eta1,eta2)
        end do
      end do
      call sll_hdf5_file_create("curvilinear_mesh-x1.h5",file_id,error)
      call sll_hdf5_write_array(file_id,x1,"/x1",error)
      call sll_hdf5_file_close(file_id, error)
      call sll_hdf5_file_create("curvilinear_mesh-x2.h5",file_id,error)
      call sll_hdf5_write_array(file_id,x2,"/x2",error)
      call sll_hdf5_file_close(file_id, error)
      deallocate(x1)
      deallocate(x2)

    end if

    call int2string(iplot,cplot)
    call sll_xdmf_open("f"//cplot//".xmf","curvilinear_mesh", &
      nnodes_x1,nnodes_x2,file_id,error)
    call sll_xdmf_write_array("f"//cplot,f,"values", &
      error,file_id,"Node")
    call sll_xdmf_close(file_id,error)
  end subroutine plot_f_curvilinear

#endif


  subroutine compute_curvilinear_field_2d( &
    input1, &
    input2, &
    output1, &
    output2, &
    jac_m, &
    jacobian)
    sll_real64, intent(in) :: input1
    sll_real64, intent(in) :: input2
    sll_real64, intent(out) :: output1
    sll_real64, intent(out) :: output2
    sll_real64, intent(in) :: jac_m(2,2)
    sll_real64, intent(in), optional :: jacobian
    sll_real64 :: inv_jac
    sll_real64 :: inv_j11
    sll_real64 :: inv_j12
    sll_real64 :: inv_j21
    sll_real64 :: inv_j22
    
    if(present(jacobian))then
      inv_jac = jacobian
    else
      inv_jac = jac_m(1,1)*jac_m(2,2)-jac_m(1,2)*jac_m(2,1)  
    endif
    inv_jac = 1._f64/inv_jac
    
    inv_j11 =  jac_m(2,2)*inv_jac
    inv_j12 = -jac_m(1,2)*inv_jac
    inv_j21 = -jac_m(2,1)*inv_jac
    inv_j22 =  jac_m(1,1)*inv_jac
 
    output1 = inv_j11*input1+inv_j12*input2 
    output2 = inv_j21*input1+inv_j22*input2 
    
  end subroutine compute_curvilinear_field_2d

  subroutine compute_field_from_phi_2d_curvilinear(phi,mesh_2d,transformation,A1,A2,interp2d)
    sll_real64, dimension(:,:), intent(in) :: phi
    sll_real64, dimension(:,:), intent(out) :: A1
    sll_real64, dimension(:,:), intent(out) :: A2
    type(sll_cartesian_mesh_2d), pointer :: mesh_2d
    class(sll_coordinate_transformation_2d_base), pointer :: transformation
    class(sll_interpolator_2d_base), pointer   :: interp2d
    sll_int32 :: Nc_eta1
    sll_int32 :: Nc_eta2
    sll_real64 :: eta1_min
    sll_real64 :: eta2_min
    sll_real64 :: delta_eta1
    sll_real64 :: delta_eta2
    sll_real64 :: eta1
    sll_real64 :: eta2
    sll_int32 :: i1
    sll_int32 :: i2
    
    Nc_eta1 = mesh_2d%num_cells1
    Nc_eta2 = mesh_2d%num_cells2
    eta1_min = mesh_2d%eta1_min
    eta2_min = mesh_2d%eta2_min
    delta_eta1 = mesh_2d%delta_eta1
    delta_eta2 = mesh_2d%delta_eta2

    call interp2d%compute_interpolants(phi)
    A1 = 0._f64
    A2 = 0._f64
    do i2=1,Nc_eta2+1
      eta2=eta2_min+real(i2-1,f64)*delta_eta2
      do i1=1,Nc_eta1+1
        eta1=eta1_min+real(i1-1,f64)*delta_eta1
        A1(i1,i2) = phi(i1,modulo(i2+1-1+Nc_eta2,Nc_eta2)+1)-phi(i1,modulo(i2-1-1+Nc_eta2,Nc_eta2)+1)
        A1(i1,i2) = A1(i1,i2)/(2._f64*delta_eta2)
        A1(i1,i2) = A1(i1,i2)/transformation%jacobian(eta1,eta2)
        A2(i1,i2) = phi(modulo(i1+1-1+Nc_eta1,Nc_eta1)+1,i2)-phi(modulo(i1-1-1+Nc_eta1,Nc_eta1)+1,i2)
        A2(i1,i2) = A2(i1,i2)/(2._f64*delta_eta1)
        A2(i1,i2) = -A2(i1,i2)/transformation%jacobian(eta1,eta2)
        !A1(i1,i2)=interp2d%interpolate_derivative_eta2(eta1,eta2)/transformation%jacobian(eta1,eta2)
        !A2(i1,i2)=-interp2d%interpolate_derivative_eta1(eta1,eta2)/transformation%jacobian(eta1,eta2)
      end do
    end do
   
    
  end subroutine compute_field_from_phi_2d_curvilinear


  subroutine compute_divergence_2d_curvilinear(div,mesh_2d,transformation,A1,A2,interp2d)
    sll_real64, dimension(:,:), intent(out) :: div
    sll_real64, dimension(:,:), intent(in) :: A1
    sll_real64, dimension(:,:), intent(in) :: A2
    sll_real64, dimension(:,:), allocatable :: tmp
    type(sll_cartesian_mesh_2d), pointer :: mesh_2d
    class(sll_coordinate_transformation_2d_base), pointer :: transformation
    class(sll_interpolator_2d_base), pointer   :: interp2d
    sll_int32 :: Nc_eta1
    sll_int32 :: Nc_eta2
    sll_real64 :: eta1_min
    sll_real64 :: eta2_min
    sll_real64 :: delta_eta1
    sll_real64 :: delta_eta2
    sll_real64 :: eta1
    sll_real64 :: eta2
    sll_int32 :: i1
    sll_int32 :: i2
    sll_int32 :: ierr
    
    Nc_eta1 = mesh_2d%num_cells1
    Nc_eta2 = mesh_2d%num_cells2
    eta1_min = mesh_2d%eta1_min
    eta2_min = mesh_2d%eta2_min
    delta_eta1 = mesh_2d%delta_eta1
    delta_eta2 = mesh_2d%delta_eta2
    
    SLL_ALLOCATE(tmp(Nc_eta1+1,Nc_eta2+1),ierr)

    do i2=1,Nc_eta2+1
      eta2=eta2_min+real(i2-1,f64)*delta_eta2
      do i1=1,Nc_eta1+1
        eta1=eta1_min+real(i1-1,f64)*delta_eta1
        tmp(i1,i2)=A1(i1,i2)*transformation%jacobian(eta1,eta2)
      end do
    end do
    call interp2d%compute_interpolants(tmp)
    do i2=1,Nc_eta2+1
      eta2=eta2_min+real(i2-1,f64)*delta_eta2
      do i1=1,Nc_eta1+1
        eta1=eta1_min+real(i1-1,f64)*delta_eta1
        div(i1,i2)=interp2d%interpolate_derivative_eta1(eta1,eta2)
      end do
    end do

    do i2=1,Nc_eta2+1
      eta2=eta2_min+real(i2-1,f64)*delta_eta2
      do i1=1,Nc_eta1+1
        eta1=eta1_min+real(i1-1,f64)*delta_eta1
        tmp(i1,i2)=A2(i1,i2)*transformation%jacobian(eta1,eta2)
      end do
    end do
    call interp2d%compute_interpolants(tmp)
    do i2=1,Nc_eta2+1
      eta2=eta2_min+real(i2-1,f64)*delta_eta2
      do i1=1,Nc_eta1+1
        eta1=eta1_min+real(i1-1,f64)*delta_eta1
        div(i1,i2)=div(i1,i2)+interp2d%interpolate_derivative_eta2(eta1,eta2)
      end do
    end do

   
    
  end subroutine compute_divergence_2d_curvilinear


  subroutine time_history_diagnostic_curvilinear( &
    file_id, &    
    step, &
    dt, &
    mesh_2d, &
    transformation,&
    f, &
    phi, &
    A1, &
    A2)
    sll_int32, intent(in) :: file_id
    sll_int32, intent(in) :: step
    sll_real64, intent(in) :: dt
    type(sll_cartesian_mesh_2d), pointer :: mesh_2d
    class(sll_coordinate_transformation_2d_base), pointer :: transformation
    sll_real64, dimension(:,:), intent(in) :: f
    sll_real64, dimension(:,:), intent(in) :: phi
    sll_real64, dimension(:,:), intent(in) :: A1
    sll_real64, dimension(:,:), intent(in) :: A2 
    sll_real64 :: mass
    sll_real64 :: linf
    sll_real64 :: l1
    sll_real64 :: l2
    sll_real64 :: e
    
    sll_real64, dimension(:), allocatable :: mass_array
    sll_real64, dimension(:), allocatable  :: l1_array
    sll_real64, dimension(:), allocatable  :: l2_array
    sll_real64, dimension(:), allocatable  :: e_array
    
    sll_real64 :: eta1
    sll_real64 :: eta2
    sll_real64, dimension(:),allocatable :: data
    sll_real64, dimension(1:2,1:2) :: jac_m
    sll_int32 :: i1
    sll_int32 :: i2
    sll_int32 :: Nc_eta1
    sll_int32 :: Nc_eta2
    sll_real64 :: eta1_min
    sll_real64 :: eta1_max
    sll_real64 :: eta2_min
    sll_real64 :: eta2_max
    sll_real64 :: delta_eta1
    sll_real64 :: delta_eta2
    sll_real64 :: dphi_eta1
    sll_real64 :: dphi_eta2
    sll_int32 :: ierr 

    
    Nc_eta1 = mesh_2d%num_cells1
    Nc_eta2 = mesh_2d%num_cells2
    
    
    eta1_min = mesh_2d%eta1_min
    eta1_max = mesh_2d%eta1_max
    eta2_min = mesh_2d%eta2_min
    eta2_max = mesh_2d%eta2_max
    delta_eta1 = mesh_2d%delta_eta1
    delta_eta2 = mesh_2d%delta_eta2

    SLL_ALLOCATE(data(Nc_eta1+1),ierr)
    SLL_ALLOCATE(mass_array(Nc_eta2+1),ierr)
    SLL_ALLOCATE(l1_array(Nc_eta2+1),ierr)
    SLL_ALLOCATE(l2_array(Nc_eta2+1),ierr)
    SLL_ALLOCATE(e_array(Nc_eta2+1),ierr)
 
    linf  = 0.0_f64
    !l1    = 0.0_f64
    !l2    = 0.0_f64
    !mass  = 0.0_f64
     !e     = 0.0_f64
    
    do i2 = 1, Nc_eta2+1
      eta2 = eta2_min + (i2-1)* delta_eta2 
      do i1=1,Nc_eta1+1
        eta1 = eta1_min + (i1-1)* delta_eta1
        data(i1) = f(i1,i2)*abs(transformation%jacobian(eta1,eta2))
      enddo
      mass_array(i2) = compute_integral_trapezoid_1d(data, Nc_eta1+1, delta_eta1)

      do i1=1,Nc_eta1+1
        eta1 = eta1_min + (i1-1)* delta_eta1
        data(i1) = abs(f(i1,i2))*abs(transformation%jacobian(eta1,eta2))
      enddo
      l1_array(i2) = compute_integral_trapezoid_1d(data, Nc_eta1+1, delta_eta1)

      do i1=1,Nc_eta1+1
        eta1 = eta1_min + (i1-1)* delta_eta1
        data(i1) = (f(i1,i2))**2 *abs(transformation%jacobian(eta1,eta2))
      enddo
      l2_array(i2) = compute_integral_trapezoid_1d(data, Nc_eta1+1, delta_eta1)

      do i1=1,Nc_eta1+1
        eta1 = eta1_min + (i1-1)* delta_eta1
        jac_m  =  transformation%jacobian_matrix(eta1,eta2)
        dphi_eta1 = -A2(i1,i2)* transformation%jacobian(eta1,eta2)
        dphi_eta2 = A1(i1,i2)* transformation%jacobian(eta1,eta2)
        data(i1) = (( jac_m(2,2)*dphi_eta1 - jac_m(2,1)*dphi_eta2 )**2 + &
        ( -jac_m(1,2)*dphi_eta1 + jac_m(1,1)*dphi_eta2 )**2) &
        /abs(transformation%jacobian(eta1,eta2)) 
      enddo
      e_array(i2) = compute_integral_trapezoid_1d(data, Nc_eta1+1, delta_eta1)

      do i1=1,Nc_eta1+1
       linf = max(linf,abs(f(i1,i2)))
      enddo
         
    enddo     

    mass = compute_integral_trapezoid_1d(mass_array, Nc_eta2+1, delta_eta2)
    l1 = compute_integral_trapezoid_1d(l1_array, Nc_eta2+1, delta_eta2)
    l2 = compute_integral_trapezoid_1d(l2_array, Nc_eta2+1, delta_eta2)
    l2 = sqrt(l2)
    e = compute_integral_trapezoid_1d(e_array, Nc_eta2+1, delta_eta2)

    !mass = mass*delta_eta2
    !l1 = l1*delta_eta2
    !l2 = sqrt(l2*delta_eta2)
    !e  = e*delta_eta2
    
    write(file_id,*) &
      dt*real(step,f64), &
      linf, &
      l1, &
      l2, &
      mass, &
      e, &
      maxval(abs(phi(1:Nc_eta1+1,1:Nc_eta2+1)))
   
    
  end subroutine time_history_diagnostic_curvilinear

  ! periodic case
  ! called when bc_type = SLL_PERIODIC
  function process_outside_point_periodic1( eta, eta_min, eta_max ) result(eta_out)
      use sll_working_precision
      sll_real64, intent(in)  :: eta
      sll_real64, intent(in) :: eta_min
      sll_real64, intent(in) :: eta_max
      sll_real64 :: eta_out

      eta_out = (eta-eta_min)/(eta_max-eta_min)      
      eta_out = eta_out-floor(eta_out)
      if(eta_out==1._f64)then
        eta_out = 0._f64
      endif      
      if(.not.((eta_out>=0).and.(eta_out<1)))then
        print *,'#eta=',eta
        print *,'#eta_min=',eta_min
        print *,'#eta_max=',eta_max
        print *,'#(eta-eta_min)/(eta_max-eta_min)=',(eta-eta_min)/(eta_max-eta_min)
        print *,'#floor(-1e-19)',floor(-1e-19)
        print *,'#eta_out=',eta_out
      endif      
      SLL_ASSERT((eta_out>=0).and.(eta_out<1))
      eta_out = eta_min+eta_out*(eta_max-eta_min) 
      SLL_ASSERT((eta_out>=eta_min).and.(eta_out<eta_max))      
      
  end function process_outside_point_periodic1

  ! set to limit case
  ! called when bc_type = SLL_SET_TO_LIMIT
  
  function process_outside_point_set_to_limit1( eta, eta_min, eta_max ) result(eta_out)
      use sll_working_precision
      sll_real64, intent(in)  :: eta
      sll_real64, intent(in) :: eta_min
      sll_real64, intent(in) :: eta_max
      sll_real64 :: eta_out
      
      eta_out = (eta-eta_min)/(eta_max-eta_min)      
      if(eta_out>1)then
        eta_out = 1._f64
      endif
      if(eta_out<0)then
        eta_out = 0._f64
      endif
      SLL_ASSERT((eta_out>=0).and.(eta_out<=1))      
      eta_out = eta_min+eta_out*(eta_max-eta_min) 
      SLL_ASSERT((eta_out>=eta_min).and.(eta_out<=eta_max))      
      
  end function process_outside_point_set_to_limit1





end module sll_simulation_2d_analytic_field_curvilinear_module
