module sll_simulation_2d_guiding_center_curvilinear_module

!2d guiding center cartesian simulation
!contact: Adnane Hamiaz (hamiaz@math.unistra.fr
!         Michel Mehrenberger (mehrenbe@math.unistra.fr)

#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"
!#include "sll_field_2d.h"
#include "sll_utilities.h"
#include "sll_poisson_solvers.h"
!  use sll_constants
!  use sll_logical_meshes  
!  use sll_module_advection_1d_periodic
!  use sll_module_interpolators_1d_base
  use sll_module_advection_2d_base
  use sll_module_characteristics_1d_explicit_euler
  use sll_module_advection_2d_BSL
  use sll_module_advection_2d_tensor_product
  use sll_module_characteristics_2d_explicit_euler
  use sll_module_characteristics_2d_verlet
  use sll_module_advection_1d_BSL
  use sll_module_advection_1d_CSL
  use sll_module_advection_1d_PSM
!  use sll_module_characteristics_1d_explicit_euler
  use sll_module_characteristics_1d_trapezoid
  use sll_module_characteristics_1d_explicit_euler_conservative
  use sll_module_characteristics_1d_trapezoid_conservative
  use sll_reduction_module
  use sll_simulation_base
 
  use sll_cubic_spline_interpolator_1d
  use sll_cubic_spline_interpolator_2d
!  use sll_coordinate_transformation_2d_base_module
  use sll_module_coordinate_transformations_2d
  use sll_common_coordinate_transformations
  use sll_common_array_initializers_module
  
#ifdef MUDPACK
 !use sll_mudpack_curvilinear
  use sll_module_poisson_2d_mudpack_curvilinear_solver_old
#endif
  use sll_module_poisson_2d_elliptic_solver
!  use sll_module_scalar_field_2d_base
!  use sll_module_scalar_field_2d_alternative
!  use sll_timer
!  use sll_fft
  implicit none
  
  
  sll_int32, parameter :: SLL_EULER = 0 
  sll_int32, parameter :: SLL_PREDICTOR_CORRECTOR = 1 
  sll_int32, parameter :: SLL_ADVECTIVE = 0 
  sll_int32, parameter :: SLL_CONSERVATIVE = 1

  type, extends(sll_simulation_base_class) :: &
    sll_simulation_2d_guiding_center_curvilinear

   !geometry
   type(sll_logical_mesh_2d), pointer :: mesh_2d
  
   !transformation 
   class(sll_coordinate_transformation_2d_base), pointer :: transformation
  
   ! Simulation should carry the pointers to the heap-allocated objects itself
    class(sll_interpolator_2d_base), pointer :: f_interp2d => null()
    class(sll_interpolator_2d_base), pointer :: phi_interp2d => null()
    class(sll_characteristics_2d_base), pointer :: charac2d => null()
    class(sll_characteristics_1d_base), pointer :: charac1d_x1 => null()
    class(sll_characteristics_1d_base), pointer :: charac1d_x2 => null()
    class(sll_interpolator_2d_base), pointer   :: A1_interp2d => null()
    class(sll_interpolator_2d_base), pointer   :: A2_interp2d => null()
    class(sll_interpolator_1d_base), pointer   :: A1_interp1d_x1 => null()
    class(sll_interpolator_1d_base), pointer   :: A2_interp1d_x1 => null()
    class(sll_interpolator_1d_base), pointer   :: A1_interp1d_x2 => null()
    class(sll_interpolator_1d_base), pointer   :: A2_interp1d_x2 => null()
    class(sll_interpolator_1d_base), pointer :: f_interp1d_x1 => null()
    class(sll_interpolator_1d_base), pointer :: f_interp1d_x2 => null()
    class(sll_advection_1d_base), pointer    :: advect_1d_x1 => null()
    class(sll_advection_1d_base), pointer    :: advect_1d_x2 => null()
   
   !initial function
   procedure(sll_scalar_initializer_2d), nopass, pointer :: init_func
   sll_real64, dimension(:), pointer :: params
      
   !advector
   class(sll_advection_2d_base), pointer    :: advect_2d
   
   !interpolator for derivatives
!   class(sll_interpolator_2d_base), pointer   :: phi_interp2d
   !coef
   sll_real64, dimension(:,:), pointer :: b11
   sll_real64, dimension(:,:), pointer :: b12
   sll_real64, dimension(:,:), pointer :: b21
   sll_real64, dimension(:,:), pointer :: b22
   sll_real64, dimension(:,:), pointer :: b1
   sll_real64, dimension(:,:), pointer :: b2
   sll_real64, dimension(:,:), pointer :: c
   
   !poisson solver
   class(sll_poisson_2d_base), pointer   :: poisson
   !type(poisson_2d_periodic), pointer   :: poisson
   !type(sll_plan_poisson_polar), pointer :: poisson 
#ifdef MUDPACK
    type(mudpack_2d) :: poisson2
#endif    
   !time_iterations
   sll_real64 :: dt
   sll_int32  :: num_iterations
   sll_int32  :: freq_diag
   sll_int32  :: freq_diag_time

   !time_loop
   sll_int32 :: time_loop_case
   ! quadrature 
   sll_int32  :: quadrature_type1
   sll_int32  :: quadrature_type2
   !boundaries conditions 
   sll_int32  :: bc_eta1_left
   sll_int32  :: bc_eta1_right
   sll_int32  :: bc_eta2_left
   sll_int32  :: bc_eta2_right
   sll_int32  :: bc_interp2d_eta1
   sll_int32  :: bc_interp2d_eta2   
   sll_int32  :: bc_charac2d_eta1
   sll_int32  :: bc_charac2d_eta2
   ! for QNS spline_degre in each direction
   sll_int32  :: spline_degree_eta1
   sll_int32  :: spline_degree_eta2   
   sll_int32  :: advection_form
  contains
    procedure, pass(sim) :: run => run_gc2d_curvilinear
    procedure, pass(sim) :: init_from_file => init_fake
     
  end type sll_simulation_2d_guiding_center_curvilinear


  abstract interface
    function sll_scalar_initializer_2d( x1, x2, params )
      use sll_working_precision
      sll_real64                                     :: sll_scalar_initializer_2d
      sll_real64, intent(in)                         :: x1
      sll_real64, intent(in)                         :: x2
      sll_real64, dimension(:), intent(in), optional :: params
    end function sll_scalar_initializer_2d
  end interface



contains

  function new_guiding_center_2d_curvilinear(filename) result(sim)
    type(sll_simulation_2d_guiding_center_curvilinear), pointer :: sim 
    character(len=*), intent(in), optional :: filename   
    sll_int32 :: ierr
    
    SLL_ALLOCATE(sim,ierr)
    
    call initialize_guiding_center_2d_curvilinear(sim,filename)
    
  
  
  end function new_guiding_center_2d_curvilinear
  
  subroutine initialize_guiding_center_2d_curvilinear(sim,filename)
    class(sll_simulation_2d_guiding_center_curvilinear), intent(inout) :: sim
    character(len=*), intent(in), optional :: filename
    sll_int32             :: IO_stat
    sll_int32, parameter  :: input_file = 99
        
    !geometry
    character(len=256) :: mesh_case
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
    character(len=256) :: advect1d_x1_case 
    character(len=256) :: advect1d_x2_case 
    character(len=256) :: charac1d_x1_case
    character(len=256) :: charac1d_x2_case
    character(len=256) :: f_interp1d_x1_case
    character(len=256) :: f_interp1d_x2_case
    character(len=256) :: advection_form
    
    !poisson
    character(len=256) :: poisson_solver
    sll_int32 :: mudpack_method    
    sll_int32 :: spline_degree_eta1
    sll_int32 :: spline_degree_eta2
    
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
    sll_real64 :: r_minus
    sll_real64 :: r_plus
    sll_int32 :: visu_step
!!$    class(sll_interpolator_2d_base), pointer :: f_interp2d
!!$    class(sll_interpolator_2d_base), pointer :: phi_interp2d
!!$    class(sll_characteristics_2d_base), pointer :: charac2d
!!$    class(sll_characteristics_1d_base), pointer :: charac1d_x1
!!$    class(sll_characteristics_1d_base), pointer :: charac1d_x2
!!$    class(sll_interpolator_2d_base), pointer   :: A1_interp2d
!!$    class(sll_interpolator_2d_base), pointer   :: A2_interp2d
!!$    class(sll_interpolator_1d_base), pointer   :: A1_interp1d_x1
!!$    class(sll_interpolator_1d_base), pointer   :: A2_interp1d_x1
!!$    class(sll_interpolator_1d_base), pointer   :: A1_interp1d_x2
!!$    class(sll_interpolator_1d_base), pointer   :: A2_interp1d_x2
!!$    class(sll_interpolator_1d_base), pointer :: f_interp1d_x1
!!$    class(sll_interpolator_1d_base), pointer :: f_interp1d_x2
!!$    class(sll_advection_1d_base), pointer    :: advect_1d_x1
!!$    class(sll_advection_1d_base), pointer    :: advect_1d_x2
    sll_int32 :: ierr
    sll_real64, dimension(4) :: params_mesh
    sll_real64 :: eta1_min_bis
    sll_real64 :: eta1_max_bis
    sll_real64 :: eta2_min_bis
    sll_real64 :: eta2_max_bis
    sll_int32  :: Nc_eta1_bis
    sll_int32  :: Nc_eta2_bis
    !here we do all the initialization
    !in future, we will use namelist file
    
    namelist /geometry/ &
      mesh_case, &
      num_cells_eta1, &
      num_cells_eta2, &
      eta1_min, &
      eta1_max, &
      eta2_min, &
      eta2_max, &
      alpha1  , &
      alpha2

    namelist /initial_function/ &
      initial_function_case, &
      kmode_eta1, &
      kmode_eta2, &
      r_minus  ,&
      r_plus   ,&
      eps

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
      charac1d_x1_case, &
      charac1d_x2_case, &
      advect1d_x1_case, &   
      advect1d_x2_case, &
      advection_form  

    namelist /poisson/ &
      poisson_solver, &
      mudpack_method, &    
      spline_degree_eta1, &
      spline_degree_eta2    
      
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
    mesh_case="SLL_LANDAU_MESH"
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
    
    advection_form = "SLL_CONSERVATIVE"
    advect1d_x1_case = "SLL_BSL"
    advect1d_x2_case = "SLL_BSL"
    charac1d_x1_case = "SLL_EULER"
    !charac1d_x1_case = "SLL_TRAPEZOID"
    charac1d_x2_case = "SLL_EULER"
    !charac1d_x2_case = "SLL_TRAPEZOID"
    f_interp1d_x1_case = "SLL_CUBIC_SPLINES"
    f_interp1d_x2_case = "SLL_CUBIC_SPLINES"
    
    !poisson 
    !poisson_solver = "SLL_ELLIPTIC_FINITE_ELEMENT_SOLVER" !use with "SLL_PHI_FROM_RHO"
    poisson_solver = "SLL_MUDPACK_CURVILINEAR"   !use with "SLL_PHI_FROM_RHO"    
    !mudpack_method = SLL_NON_SEPARABLE_WITH_CROSS_TERMS  
#ifdef MUDPACK

    mudpack_method = SLL_NON_SEPARABLE_WITHOUT_CROSS_TERMS  
#else
    mudpack_method = 0
#endif
    spline_degree_eta1 = 3
    spline_degree_eta2 = 3    
     
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
          print *, '#initialize_guiding_center_2d_curvilinear() failed to open file ', &
          trim(filename)//'.nml'
          STOP
        end if
      print *,'#initialization with filename:'
      print *,'#',trim(filename)//'.nml'
      read(input_file, geometry) 
      read(input_file, initial_function)
      read(input_file, time_iterations)
      read(input_file, advector)
      read(input_file, poisson)
      read(input_file, boundaries)
      close(input_file)
    else
      print *,'#initialization with default parameters'    
    endif

 
!   initial_function_case = "SLL_DIOCOTRON"
!    eta1_min = 1._f64
!    eta1_max = 10._f64
!    eta2_min = 0._f64
!    eta2_max = 2._f64*sll_pi
!    r_minus = 4._f64
!    r_plus = 5._f64
!    k_mode = 3
!    eps = 1.e-6_f64

    Nc_eta1 =  num_cells_eta1
    Nc_eta2 =  num_cells_eta2
    sim%dt = dt
    sim%num_iterations = number_iterations
    sim%freq_diag = freq_diag 
    sim%freq_diag_time = freq_diag_time
    sim%spline_degree_eta1 = spline_degree_eta1
    sim%spline_degree_eta2 = spline_degree_eta2
    sim%quadrature_type1 = ES_GAUSS_LEGENDRE
    sim%quadrature_type2 = ES_GAUSS_LEGENDRE
      
    SLL_ALLOCATE(sim%b11(Nc_eta1+1,Nc_eta2+1),ierr)
    SLL_ALLOCATE(sim%b12(Nc_eta1+1,Nc_eta2+1),ierr)
    SLL_ALLOCATE(sim%b21(Nc_eta1+1,Nc_eta2+1),ierr)
    SLL_ALLOCATE(sim%b22(Nc_eta1+1,Nc_eta2+1),ierr)
    SLL_ALLOCATE(sim%b1(Nc_eta1+1,Nc_eta2+1),ierr)
    SLL_ALLOCATE(sim%b2(Nc_eta1+1,Nc_eta2+1),ierr)
    SLL_ALLOCATE(sim%c(Nc_eta1+1,Nc_eta2+1),ierr)
    
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
      case ("SLL_DIRICHLET")
        print*,"#bc_charac2d_eta1 = SLL_DIRICHLET"  
        sim%bc_charac2d_eta1= SLL_DIRICHLET
      case ("SLL_HERMITE")
        print*,"#bc_charac2d_eta1 = SLL_HERMITE"  
        sim%bc_charac2d_eta1= SLL_HERMITE 
      case ("SLL_SET_TO_LIMIT")
        print*,"#bc_charac2d_eta1 = SLL_SET_TO_LIMIT"  
        sim%bc_charac2d_eta1= SLL_SET_TO_LIMIT    
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
      case ("SLL_DIRICHLET")
        print*,"#bc_charac2d_eta2= SLL_DIRICHLET"  
        sim%bc_charac2d_eta2= SLL_DIRICHLET
      case ("SLL_HERMITE")
        print*,"#bc_charac2d_eta2 = SLL_HERMITE"  
        sim%bc_charac2d_eta2= SLL_HERMITE 
      case ("SLL_SET_TO_LIMIT")
        print*,"#bc_charac2d_eta2 = SLL_SET_TO_LIMIT"  
        sim%bc_charac2d_eta2= SLL_SET_TO_LIMIT   
      case default
        print *,'#bad bc_charac2d_eta2',bc_charac2d_eta2
        print *,'#not implemented'
        print *,'#in initialize_analytic_field_2d_curvilinear'
        stop
    end select
    
    select case (mesh_case)
      case ("SLL_LANDAU_MESH")
        eta1_max = 2._f64*sll_pi/kmode_eta1
        eta2_max = 2._f64*sll_pi/kmode_eta2
      case ("SLL_POLAR_MESH")
        eta2_max = 2._f64*sll_pi
      case ("SLL_COLLELA_MESH")  
        eta1_max = 2._f64*sll_pi/kmode_eta1
        eta2_max = 2._f64*sll_pi/kmode_eta2
         !  In collela  mesh params_mesh =( alpha1, alpha2, L1, L2 ) such that :
         !  x1= eta1 + alpha1*sin(2*pi*eta1/L1)*sin(2*pi*eta2/L2)
         params_mesh = (/ alpha1, alpha2, eta1_max-eta1_min, eta2_max-eta2_min/)
      case default
        print *,'#bad mesh_case',mesh_case
        print *,'#not implemented'
        print *,'#in initialize_guiding_center_2d_curvilinear'
        stop
    end select

        
    sim%mesh_2d => new_logical_mesh_2d( &
      Nc_eta1, &
      Nc_eta2, &
      eta1_min , &
      eta1_max , &
      eta2_min , &
      eta2_max ) 
      
    select case (mesh_case)
      case ("SLL_LANDAU_MESH")       
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
      case ("SLL_POLAR_MESH") 
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
      case ("SLL_COLLELA_MESH")
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
        print *,'#bad mesh_case',mesh_case
        print *,'#not implemented'
        print *,'#in initialize_guiding_center_2d_curvilinear'
        stop
    end select  
     
    select case(advection_form)
      case ("SLL_CONSERVATIVE")
        sim%advection_form = SLL_CONSERVATIVE
      case ("SLL_ADVECTIVE")  
        sim%advection_form = SLL_ADVECTIVE
      case default
        print *,'#bad advection_form',advection_form
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
        eta1_min_bis = eta1_min-0.5_f64*sim%mesh_2d%delta_eta1
        eta1_max_bis = eta1_max-0.5_f64*sim%mesh_2d%delta_eta1
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
        eta2_min_bis = eta2_min-0.5_f64*sim%mesh_2d%delta_eta2
        eta2_max_bis = eta2_max-0.5_f64*sim%mesh_2d%delta_eta2        
        Nc_eta2_bis = Nc_eta2
      case default
        print *,'#bad value of advect1d_x2_case'
        stop  
    end select
      
    select case (f_interp2d_case)
      case ("SLL_CUBIC_SPLINES")
        print*,"#f interpolation SLL_CUBIC_SPLINES"
        sim%f_interp2d => new_cubic_spline_2d_interpolator( &
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
        print *,'#in initialize_guiding_center_2d_curvilinear'
        stop
    end select




    select case (A_interp_case)
      case ("SLL_CUBIC_SPLINES")
       print*,"#A1_2d interpolation SLL_CUBIC_SPLINES"
        sim%A1_interp2d => new_cubic_spline_2d_interpolator( &
          Nc_eta1+1, &
          Nc_eta2+1, &
          eta1_min, &
          eta1_max, &
          eta2_min, &
          eta2_max, &
          sim%bc_interp2d_eta1, &
          sim%bc_interp2d_eta2)
       print*,"#A2_2d interpolation SLL_CUBIC_SPLINES"   
        sim%A2_interp2d => new_cubic_spline_2d_interpolator( &
          Nc_eta1+1, &
          Nc_eta2+1, &
          eta1_min, &
          eta1_max, &
          eta2_min, &
          eta2_max, &
          sim%bc_interp2d_eta1, &
          sim%bc_interp2d_eta2)  
       print*,"#A1_1d interpolation SLL_CUBIC_SPLINES"   
        sim%A1_interp1d_x1 => new_cubic_spline_1d_interpolator( &
          Nc_eta1+1, &
          eta1_min, &
          eta1_max, &
          sim%bc_interp2d_eta1)
        sim%A1_interp1d_x2 => new_cubic_spline_1d_interpolator( &
          Nc_eta2+1, &
          eta2_min, &
          eta2_max, &
          sim%bc_interp2d_eta2)
       print*,"#A2_1d interpolation SLL_CUBIC_SPLINES"     
        sim%A2_interp1d_x1 => new_cubic_spline_1d_interpolator( &
          Nc_eta1+1, &
          eta1_min, &
          eta1_max, &
          sim%bc_interp2d_eta1)
        sim%A2_interp1d_x2 => new_cubic_spline_1d_interpolator( &
          Nc_eta2+1, &
          eta2_min, &
          eta2_max, &
          sim%bc_interp2d_eta2)
      case default
        print *,'#bad A_interp_case',A_interp_case
        print *,'#not implemented'
        print *,'#in initialize_guiding_center_2d_curvilinear'
        stop
    end select

    select case (phi_interp2d_case)
      case ("SLL_CUBIC_SPLINES")
      print*,"#phi interpolation SLL_CUBIC_SPLINES"  
        sim%phi_interp2d => new_cubic_spline_2d_interpolator( &
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
        print *,'#in initialize_guiding_center_2d_curvilinear'
        stop
    end select

     select case (f_interp1d_x1_case)
      case ("SLL_CUBIC_SPLINES")
        sim%f_interp1d_x1 => new_cubic_spline_1d_interpolator( &
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
        sim%f_interp1d_x2 => new_cubic_spline_1d_interpolator( &
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
        sim%charac1d_x1 => new_explicit_euler_1d_charac(&
          Nc_eta1_bis+1, &
          eta_min=eta1_min_bis, &
          eta_max=eta1_max_bis, &
          bc_type= sim%bc_charac2d_eta1)    
      case ("SLL_TRAPEZOID")
        sim%charac1d_x1 => &
          new_trapezoid_1d_charac(&
          Nc_eta1_bis+1, &
          sim%A1_interp1d_x1, &
          bc_type= sim%bc_charac2d_eta1, &
          eta_min=eta1_min_bis, &
          eta_max=eta1_max_bis)
      case ("SLL_TRAPEZOID_CONSERVATIVE")
        sim%charac1d_x1 => &
          new_trapezoid_conservative_1d_charac(&
          Nc_eta1_bis+1, &
          sim%A1_interp1d_x1, &
          bc_type=SLL_PERIODIC, &
          eta_min=eta1_min_bis, &
          eta_max=eta1_max_bis)
      case ("SLL_EULER_CONSERVATIVE")
        sim%charac1d_x1 => new_explicit_euler_conservative_1d_charac(&
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
        sim%charac1d_x2 => new_explicit_euler_1d_charac(&
          Nc_eta2_bis+1, &
          eta_min=eta2_min_bis, &
          eta_max=eta2_max_bis, &
          bc_type= sim%bc_charac2d_eta2)    
      case ("SLL_TRAPEZOID")
        sim%charac1d_x2 => &
          new_trapezoid_1d_charac(&
          Nc_eta2_bis+1, &
          sim%A2_interp1d_x2, &
          bc_type= sim%bc_charac2d_eta2, &
          eta_min=eta2_min_bis, &
          eta_max=eta2_max_bis)
      case ("SLL_EULER_CONSERVATIVE")
        sim%charac1d_x2 => new_explicit_euler_conservative_1d_charac(&
          Nc_eta2_bis+1, &
          eta_min=eta2_min_bis, &
          eta_max=eta2_max_bis, &
          bc_type= sim%bc_charac2d_eta2)    
      case ("SLL_TRAPEZOID_CONSERVATIVE")
        sim%charac1d_x2 => &
          new_trapezoid_conservative_1d_charac(&
          Nc_eta2_bis+1, &
          sim%A2_interp1d_x2, &
          bc_type=SLL_PERIODIC, &
          eta_min=eta2_min_bis, &
          eta_max=eta2_max_bis)
      case default
        print *,'#bad charac1d_x2_case',charac1d_x2_case
        print *,'#not implemented'
        print *,'#in initialize_guiding_center_2d_curvilinear'
        stop
    end select

    select case(advect1d_x1_case)
      case ("SLL_BSL")
        sim%advect_1d_x1 => new_BSL_1d_advector(&
          sim%f_interp1d_x1, &
          sim%charac1d_x1, &
          Nc_eta1_bis+1, &
          eta_min = eta1_min_bis, &
          eta_max = eta1_max_bis)
      case ("SLL_CSL")
        sim%advect_1d_x1 => new_CSL_1d_advector(&
          sim%f_interp1d_x1, &
          sim%charac1d_x1, &
          Nc_eta1_bis+1, &
          eta_min = eta1_min_bis, &
          eta_max = eta1_max_bis, &
          bc_type = sim%bc_charac2d_eta1)
      case ("SLL_PSM")
        sim%advect_1d_x1 => new_PSM_1d_advector(&
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
        sim%advect_1d_x2 => new_BSL_1d_advector(&
          sim%f_interp1d_x2, &
          sim%charac1d_x2, &
          Nc_eta2_bis+1, &
          eta_min = eta2_min_bis, &
          eta_max = eta2_max_bis)
      case ("SLL_CSL")
        sim%advect_1d_x2 => new_CSL_1d_advector(&
          sim%f_interp1d_x2, &
          sim%charac1d_x2, &
          Nc_eta2_bis+1, &
          eta_min = eta2_min_bis, &
          eta_max = eta2_max_bis, &
          bc_type = sim%bc_charac2d_eta2)
      case ("SLL_PSM")
        sim%advect_1d_x2 => new_PSM_1d_advector(&
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
        sim%charac2d => new_explicit_euler_2d_charac(&
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
        sim%charac2d => new_verlet_2d_charac(&
          Nc_eta1+1, &
          Nc_eta2+1, &
          sim%A1_interp2d, &
          sim%A2_interp2d, &
          sim%A1_interp1d_x1, &
          sim%A2_interp1d_x1, &
          bc_type_1=sim%bc_charac2d_eta1, & !SLL_SET_TO_LIMIT, &
          bc_type_2=sim%bc_charac2d_eta2, &
          eta1_min=eta1_min, &
          eta1_max=eta1_max, &
          eta2_min=eta2_min, &
          eta2_max=eta2_max)
      case default
        print *,'#bad charac2d_case',charac2d_case
        print *,'#not implemented'
        print *,'#in initialize_guiding_center_2d_curvilinear'
        stop
    end select

  
    select case(advect2d_case)
      case ("SLL_BSL")
       print*,"#advect2d = SLL_BSL "  
        sim%advect_2d => new_BSL_2d_advector(&
          sim%f_interp2d, &
          sim%charac2d, &
          Nc_eta1+1, &
          Nc_eta2+1, &
          eta1_min = eta1_min, &
          eta1_max = eta1_max, &
          eta2_min = eta2_min, &
          eta2_max = eta2_max)
      case ("SLL_TENSOR_PRODUCT")
       print*,"#advect2d = SLL_SPLITING " 
        sim%advect_2d => new_tensor_product_2d_advector(&
          sim%advect_1d_x1, &
          sim%advect_1d_x2, &
          Nc_eta1+1, &
          Nc_eta2+1)        
      case default
        print *,'#bad advect_case',advect2d_case
        print *,'#not implemented'
        print *,'#in initialize_guiding_center_2d_curvilinear'
        stop
    end select
    
    
    select case(initial_function_case)
      case ("SLL_KHP1")
        print*,"#f0 = SLL_KHP1"  
        sim%init_func => sll_KHP1_2d
        SLL_ALLOCATE(sim%params(2),ierr)
        sim%params(1) = eps
        sim%params(2) = kmode_eta1
      case("SLL_DIOCOTRON")
        print*,"#f0 = SLL_DIOCOTRON " 
        sim%init_func => sll_diocotron_initializer_2d2
        SLL_ALLOCATE(sim%params(4),ierr)
        sim%params(1) = r_minus
        sim%params(2) = r_plus
        sim%params(3) = eps
        sim%params(4) = kmode_eta2  
      case default
        print *,'#bad initial_function_case',initial_function_case
        print *,'#not implemented'
        print *,'#in initialize_guiding_center_2d_curvilinear'
        stop
    end select
    
    
    !time_loop
    select case(trim(time_loop_case))
      case ("SLL_EULER")
        print*,"#time_loop = SLL_EULER " 
        sim%time_loop_case = SLL_EULER
      case ("SLL_PREDICTOR_CORRECTOR")
       print*,"#time_loop = SLL_PREDICTOR_CORRECTOR " 
        sim%time_loop_case = SLL_PREDICTOR_CORRECTOR
      case default
        print *,'#bad time_loop_case',time_loop_case
        print *,'#not implemented'
        print *,'#in initialize_guiding_center_2d_curvilinear'
        stop
    end select
    
    
    
     !poisson solver
    select case(poisson_solver)    
#ifdef MUDPACK
      case ("SLL_MUDPACK_CURVILINEAR")    
        print *,'#poisson = MUDPACK_CURVILINEAR', mudpack_method
        sim%b11 = 1._f64
        sim%b22 = 1._f64
        sim%b12 = 0._f64
        sim%b21 = 0._f64
        sim%c   = 0._f64 
         
        sim%poisson => new_poisson_2d_mudpack_curvilinear_solver( &
         sim%transformation, &
         eta1_min,&
         eta1_max,&
         Nc_eta1,&
         eta2_min,&
         eta2_max,&
         Nc_eta2,&
         sim%bc_eta1_left, &
         sim%bc_eta1_right, &
         sim%bc_eta2_left, &
         sim%bc_eta2_right, &
         sim%bc_interp2d_eta1, &
         sim%bc_interp2d_eta2, &
         sim%b11,&
         sim%b12,&
         sim%b21,&
         sim%b22,&
         sim%c, &
         mudpack_curvilinear_case = mudpack_method)
#endif
       case("SLL_ELLIPTIC_FINITE_ELEMENT_SOLVER")
        print *,'#poisson = ELLIPTIC_FINITE_ELEMENT_SOLVER '
        sim%b11 = 1._f64
        sim%b22 = 1._f64
        sim%b12 = 0._f64
        sim%b21 = 0._f64
        sim%b1  = 0._f64
        sim%b2  = 0._f64
        sim%c   = 0._f64 
        
        sim%poisson => new_poisson_2d_elliptic_solver( &
         sim%transformation,&
         sim%spline_degree_eta1, &
         sim%spline_degree_eta2, &
         Nc_eta1, &
         Nc_eta2, &
         sim%quadrature_type1, &
         sim%quadrature_type2, &
         sim%bc_eta1_left, &
         sim%bc_eta1_right, &
         sim%bc_eta2_left, &
         sim%bc_eta2_right, &
         eta1_min, &
         eta1_max, &
         eta2_min, &
         eta2_max, &
         sim%b11, & 
         sim%b12, & 
         sim%b21, & 
         sim%b22, & 
         sim%b1, & 
         sim%b2, & 
         sim%c ) 
      case default
        print *,'#bad poisson_case',poisson_solver
        print *,'#not implemented'
        print *,'#in  initialize_guiding_center_2d_curvilinear'
        stop
    end select

   
  end subroutine initialize_guiding_center_2d_curvilinear
  


  subroutine init_fake(sim, filename)
    class(sll_simulation_2d_guiding_center_curvilinear), intent(inout) :: sim
    character(len=*), intent(in)                                :: filename
  
    print *,'# Do not use the routine init_vp4d_fake'
    print *,'#use instead initialize_vlasov_par_poisson_seq_curv'
    stop
  
  end subroutine init_fake
  
  subroutine run_gc2d_curvilinear(sim)
    class(sll_simulation_2d_guiding_center_curvilinear), intent(inout) :: sim
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
    sll_real64,dimension(:,:), pointer :: f
    sll_real64,dimension(:,:), pointer :: f_conserv
    sll_real64,dimension(:,:), pointer :: f_old
    sll_real64,dimension(:,:), pointer :: rho
    sll_real64,dimension(:,:), pointer :: phi
    sll_real64,dimension(:,:), pointer :: A1 !advection fields
    sll_real64,dimension(:,:), pointer :: A2
    sll_int32 :: ierr
    sll_int32 :: nb_step
    sll_int32 :: step
    sll_real64 :: dt
    sll_int32 :: thdiag_id
    sll_int32 :: thdiagp_id
    sll_int32 :: iplot
    
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
    SLL_ALLOCATE(rho(Nc_eta1+1,Nc_eta2+1),ierr)
    SLL_ALLOCATE(phi(Nc_eta1+1,Nc_eta2+1),ierr)
    SLL_ALLOCATE(A1(Nc_eta1+1,Nc_eta2+1),ierr)
    SLL_ALLOCATE(A2(Nc_eta1+1,Nc_eta2+1),ierr)
    SLL_ALLOCATE(f_conserv(Nc_eta1+1,Nc_eta2+1),ierr)
    

    
    !initialisation of distribution function    
     do i2=1,Nc_eta2+1
        eta2=eta2_min+real(i2-1,f64)*delta_eta2
        do i1=1,Nc_eta1+1
          eta1=eta1_min+real(i1-1,f64)*delta_eta1
          x1 = sim%transformation%x1(eta1,eta2)
          x2 = sim%transformation%x2(eta1,eta2)
          f(i1,i2) =  sim%init_func(x1,x2,sim%params) 
        end do
     end do
        
    !solve poisson
    call sim%poisson%compute_phi_from_rho(phi, f)
    call compute_field_from_phi_2d_curvilinear(phi,sim%mesh_2d,sim%transformation,A1,A2,sim%phi_interp2d)
    
    !print *,A1
    !print *,A2
    thdiagp_id=thdiag_id+1
    call sll_ascii_file_create('thdiag.dat', thdiag_id, ierr)
    call sll_ascii_file_create('thdiagp.dat', thdiagp_id, ierr)
!    open(unit = diag_id, file='diag_curvilinear.dat',IOStat=IO_stat)
!    if( IO_stat /= 0 ) then
!       print *, '#run_gc2d_curvilinear (sim) failed to open file diag_curvilinear.dat'
!       STOP
!    end if
!    open(unit = diag_id+1, file='diag2_curvilinear.dat',IOStat=IO_stat)
!    if( IO_stat /= 0 ) then
!       print *, '#run_gc2d_curvilinear (sim) failed to open file diag_curvilinear.dat'
!       STOP
!    end if
    
    iplot = 0

    do step=1,nb_step+1
      print*,"step= ", step
      f_old = f
      
      call sim%poisson%compute_phi_from_rho(phi, f_old) 
      call compute_field_from_phi_2d_curvilinear(phi,sim%mesh_2d,sim%transformation,A1,A2,sim%phi_interp2d)      
      if(modulo(step-1,sim%freq_diag_time)==0)then
        call time_history_diagnostic_gc( &
          thdiag_id , &    
          step-1, &
          dt, &
          sim%mesh_2d, &
          f, &
          phi, &
          A1, &
          A2)
        call time_history_diagnostic_gc2( &
          thdiagp_id, &    
          step-1, &
          dt, &
          sim%mesh_2d, &
          f, &
          phi, &
          A1, &
          A2)
      endif            
#ifndef NOHDF5
      if(modulo(step-1,sim%freq_diag)==0)then
        call plot_f_curvilinear(iplot,f,sim%mesh_2d,sim%transformation)
        iplot = iplot+1  
      endif            
#endif  
      select case (sim%time_loop_case)
        case (SLL_EULER)
          select case (sim%advection_form)
            case(SLL_CONSERVATIVE)
              call compute_f_to_f_conserv(f_old,f_conserv,sim%transformation,sim%mesh_2d)
              call sim%advect_2d%advect_2d(A1, A2, dt, f_conserv, f)
              f_conserv =f
              call compute_f_conserv_to_f(f_conserv,f,sim%transformation,sim%mesh_2d)
            case(SLL_ADVECTIVE)  
              call sim%advect_2d%advect_2d(A1, A2, dt, f_old, f)
           case default 
              print *,'#bad advection_form',sim%advection_form
          end select    
        case (SLL_PREDICTOR_CORRECTOR)
          select case (sim%advection_form)
            case(SLL_CONSERVATIVE)
              call compute_f_to_f_conserv(f_old,f_conserv,sim%transformation,sim%mesh_2d)
              call sim%advect_2d%advect_2d(A1, A2, 0.5_f64*dt, f_conserv, f)
              f_conserv =f
              call compute_f_conserv_to_f(f_conserv,f,sim%transformation,sim%mesh_2d)
            case(SLL_ADVECTIVE)  
              call sim%advect_2d%advect_2d(A1, A2, 0.5_f64*dt, f_old, f)
           case default 
              print *,'#bad advection_form',sim%advection_form
          end select    
          !call sim%advect_2d%advect_2d(A1, A2, 0.5_f64*dt, f_old, f)
          call sim%poisson%compute_phi_from_rho(phi, f)
          call compute_field_from_phi_2d_curvilinear(phi,sim%mesh_2d,sim%transformation,A1,A2,sim%phi_interp2d)      
          !call sim%advect_2d%advect_2d(A1, A2, dt, f_old, f)
          select case (sim%advection_form)
            case(SLL_CONSERVATIVE)
              call compute_f_to_f_conserv(f_old,f_conserv,sim%transformation,sim%mesh_2d)
              call sim%advect_2d%advect_2d(A1, A2, dt, f_conserv, f)
              f_conserv =f
              call compute_f_conserv_to_f(f_conserv,f,sim%transformation,sim%mesh_2d)
            case(SLL_ADVECTIVE)  
              call sim%advect_2d%advect_2d(A1, A2, dt, f_old, f)
           case default 
              print *,'#bad advection_form',sim%advection_form
          end select    
        case default  
          print *,'#bad time_loop_case',sim%time_loop_case
          print *,'#not implemented'
          print *,'#in run_gc2d_curvilinear'
          print *,'#available options are:'
          print *,'#SLL_EULER=',SLL_EULER
          print *,'#SLL_PREDICTOR_CORRECTOR=',SLL_PREDICTOR_CORRECTOR
          
      end select
    enddo
    
    close(thdiag_id)
    close(thdiagp_id)
    !print *,'#not implemented for the moment!'
  end subroutine run_gc2d_curvilinear
  
  
  subroutine compute_field_from_phi_2d_curvilinear(phi,mesh_2d,transformation,A1,A2,interp2d)
    sll_real64, dimension(:,:), intent(in) :: phi
    sll_real64, dimension(:,:), intent(out) :: A1
    sll_real64, dimension(:,:), intent(out) :: A2
    type(sll_logical_mesh_2d), pointer :: mesh_2d
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
        A1(i1,i2)=interp2d%interpolate_derivative_eta2(eta1,eta2)/transformation%jacobian(eta1,eta2)
        A2(i1,i2)=-interp2d%interpolate_derivative_eta1(eta1,eta2)/transformation%jacobian(eta1,eta2)
      end do
    end do
   
    
  end subroutine compute_field_from_phi_2d_curvilinear


  subroutine compute_f_to_f_conserv(f,f_conserv,transformation,mesh_2d)
    sll_real64, dimension(:,:), intent(in) :: f
    sll_real64, dimension(:,:), intent(out) :: f_conserv
    type(sll_logical_mesh_2d), pointer :: mesh_2d
    class(sll_coordinate_transformation_2d_base), pointer :: transformation
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

    f_conserv(1:Nc_eta1+1,1:Nc_eta2+1) =  0._f64
    do i2=1,Nc_eta2+1
      eta2=eta2_min+real(i2-1,f64)*delta_eta2
      do i1=1,Nc_eta1+1
        eta1=eta1_min+real(i1-1,f64)*delta_eta1
        f_conserv(i1,i2)=f(i1,i2)*transformation%jacobian(eta1,eta2)
      end do
    end do
    
      
  end subroutine compute_f_to_f_conserv
  
  
   subroutine compute_f_conserv_to_f(f_conserv,f,transformation,mesh_2d)
    sll_real64, dimension(:,:), intent(out) :: f
    sll_real64, dimension(:,:), intent(in) :: f_conserv
    type(sll_logical_mesh_2d), pointer :: mesh_2d
    class(sll_coordinate_transformation_2d_base), pointer :: transformation
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

    f(1:Nc_eta1+1,1:Nc_eta2+1) =  0._f64
    do i2=1,Nc_eta2+1
      eta2=eta2_min+real(i2-1,f64)*delta_eta2
      do i1=1,Nc_eta1+1
        eta1=eta1_min+real(i1-1,f64)*delta_eta1
        f(i1,i2)=f_conserv(i1,i2)/transformation%jacobian(eta1,eta2)
      end do
    end do
    
      
  end subroutine compute_f_conserv_to_f
  subroutine time_history_diagnostic_gc( &
    file_id, &    
    step, &
    dt, &
    mesh_2d, &
    f, &
    phi, &
    A1, &
    A2)
    sll_int32, intent(in) :: file_id
    sll_int32, intent(in) :: step
    sll_real64, intent(in) :: dt
    type(sll_logical_mesh_2d), pointer :: mesh_2d
    sll_real64, dimension(:,:), intent(in) :: f
    sll_real64, dimension(:,:), intent(in) :: phi
    sll_real64, dimension(:,:), intent(in) :: A1
    sll_real64, dimension(:,:), intent(in) :: A2 
    sll_real64 :: mass
    sll_real64 :: linf
    sll_real64 :: l1
    sll_real64 :: l2
    sll_real64 :: e
    sll_real64, dimension(:),allocatable :: data
    sll_int32 :: i1
    sll_int32 :: i2
    sll_int32 :: Nc_eta1
    sll_int32 :: Nc_eta2
    sll_real64 ::eta1_min
    sll_real64 ::eta1_max
    sll_real64 :: delta_eta1
    sll_real64 :: delta_eta2
    sll_int32 :: ierr 

    
    Nc_eta1 = mesh_2d%num_cells1
    Nc_eta2 = mesh_2d%num_cells2
    
    
    eta1_min = mesh_2d%eta1_min
    eta1_max = mesh_2d%eta1_max

    delta_eta1 = mesh_2d%delta_eta1
    delta_eta2 = mesh_2d%delta_eta2

    SLL_ALLOCATE(data(Nc_eta1+1),ierr)
 
    linf  = 0.0_f64
    l1    = 0.0_f64
    l2    = 0.0_f64
    mass  = 0.0_f64
     e     = 0.0_f64
    
    do i2 = 1, Nc_eta2+1
      do i1=1,Nc_eta1+1
        data(i1) = f(i1,i2)
      enddo
      mass = mass + compute_integral_trapezoid_1d(data, Nc_eta1+1, delta_eta1)

      do i1=1,Nc_eta1+1
        data(i1) = abs(f(i1,i2))
      enddo
      l1 = l1 + compute_integral_trapezoid_1d(data, Nc_eta1+1, delta_eta1)

      do i1=1,Nc_eta1+1
        data(i1) = (f(i1,i2))**2
      enddo
      l2 = l2 + compute_integral_trapezoid_1d(data, Nc_eta1+1, delta_eta1)

      do i1=1,Nc_eta1+1
        data(i1) = A2(i1,i2)**2+A1(i1,i2)**2
      enddo
      e = e + compute_integral_trapezoid_1d(data, Nc_eta1+1, delta_eta1)

      do i1=1,Nc_eta1+1
       linf = max(linf,abs(f(i1,i2)))
      enddo
         
    enddo     

    mass = mass*delta_eta2
    l1 = l1*delta_eta2
    l2 = sqrt(l2*delta_eta2)
    e  = e*delta_eta2
    
    write(file_id,*) &
      dt*real(step,f64), &
      linf, &
      l1, &
      l2, &
      mass, &
      e, &
      maxval(abs(phi(1:Nc_eta1+1,1:Nc_eta2+1)))
   
    
  end subroutine time_history_diagnostic_gc

  subroutine time_history_diagnostic_gc2( &
    file_id, &    
    step, &
    dt, &
    mesh_2d, &
    f, &
    phi, &
    A1, &
    A2)
    sll_int32, intent(in) :: file_id
    sll_int32, intent(in) :: step
    sll_real64, intent(in) :: dt
    type(sll_logical_mesh_2d), pointer :: mesh_2d
    sll_real64, dimension(:,:), intent(in) :: f
    sll_real64, dimension(:,:), intent(in) :: phi
    sll_real64, dimension(:,:), intent(in) :: A1
    sll_real64, dimension(:,:), intent(in) :: A2
    sll_real64 :: time_mode(8) 
    !sll_real64 :: mode_slope(8) 
    sll_real64 :: w
    sll_real64 :: l1
    sll_real64 :: l2
    sll_real64 :: e
    sll_real64, dimension(:),allocatable :: int_r
    sll_real64, dimension(:),allocatable :: data
    sll_int32 :: i1
    sll_int32 :: i2
    sll_int32 :: Nc_eta1
    sll_int32 :: Nc_eta2
    sll_real64 ::eta1_min
    sll_real64 :: delta_eta1
    sll_real64 :: delta_eta2
    sll_real64 :: eta1
    sll_int32 :: ierr 
    type(sll_fft_plan), pointer         :: pfwd
    
    Nc_eta1 = mesh_2d%num_cells1
    Nc_eta2 = mesh_2d%num_cells2
    
    
    eta1_min = mesh_2d%eta1_min

    delta_eta1 = mesh_2d%delta_eta1
    delta_eta2 = mesh_2d%delta_eta2

    
    SLL_ALLOCATE(int_r(Nc_eta2),ierr)
    SLL_ALLOCATE(data(Nc_eta1+1),ierr)
    pfwd => fft_new_plan(Nc_eta2,int_r,int_r,FFT_FORWARD,FFT_NORMALIZE)
 
    w     = 0.0_f64
    l1    = 0.0_f64
    l2    = 0.0_f64
    e     = 0.0_f64
    int_r = 0.0_f64
    
    
    do i2 = 1, Nc_eta2
      do i1=1,Nc_eta1+1
        eta1 = eta1_min+real(i1-1,f64)*delta_eta1
        data(i1) = eta1*f(i1,i2)
      enddo
      w = w + compute_integral_trapezoid_1d(data, Nc_eta1+1, delta_eta1)

      do i1=1,Nc_eta1+1
        eta1 = eta1_min+real(i1-1,f64)*delta_eta1
        data(i1) = eta1*abs(f(i1,i2))
      enddo
      l1 = l1 + compute_integral_trapezoid_1d(data, Nc_eta1+1, delta_eta1)

      do i1=1,Nc_eta1+1
        eta1 = eta1_min+real(i1-1,f64)*delta_eta1
        data(i1) = eta1*(f(i1,i2))**2
      enddo
      l2 = l2 + compute_integral_trapezoid_1d(data, Nc_eta1+1, delta_eta1)

      do i1=1,Nc_eta1+1
        eta1 = eta1_min+real(i1-1,f64)*delta_eta1
        data(i1) = eta1*((eta1*A2(i1,i2))**2+A1(i1,i2)**2)
      enddo
      e = e + compute_integral_trapezoid_1d(data, Nc_eta1+1, delta_eta1)

      do i1=1,Nc_eta1+1
        eta1 = eta1_min+real(i1-1,f64)*delta_eta1
        data(i1) = eta1*phi(i1,i2)
      enddo
      int_r(i2) = compute_integral_trapezoid_1d(data, Nc_eta1+1, delta_eta1)     
    enddo     

    w = w*delta_eta2
    l1 = l1*delta_eta2
    l2 = sqrt(l2*delta_eta2)
    e  = 0.5_f64*e*delta_eta2
    call fft_apply_plan(pfwd,int_r,int_r)
    do i1=1,8
      !mode_slope(i1) = time_mode(i1)
      time_mode(i1) = abs(fft_get_mode(pfwd,int_r,i1-1))**2
      !mode_slope(i1) = &
      !  (log(0*time_mode(i1)+1.e-40_f64)-log(0*mode_slope(i1)+1.e-40_f64))/(dt+1.e-40_f64)
    enddo
    
    write(file_id,*) &
      dt*real(step,f64), &
      w, &
      l1, &
      l2, &
      e, &
      maxval(abs(phi(1:Nc_eta1+1,1:Nc_eta2+1))), &
      time_mode(1:8)!,mode_slope

    call fft_delete_plan(pfwd)

    
  end subroutine time_history_diagnostic_gc2

#ifndef NOHDF5
!*********************
!*********************

  !---------------------------------------------------
  ! Save the mesh structure
  !---------------------------------------------------
  subroutine plot_f_curvilinear(iplot,f,mesh_2d,transf)
    use sll_xdmf
    use sll_hdf5_io
    sll_int32 :: file_id
    sll_int32 :: error
    sll_real64, dimension(:,:), allocatable :: x1
    sll_real64, dimension(:,:), allocatable :: x2
    sll_int32 :: i, j
    sll_int32, intent(in) :: iplot
    character(len=4)      :: cplot
    sll_int32             :: nnodes_x1, nnodes_x2
    type(sll_logical_mesh_2d), pointer :: mesh_2d
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


end module sll_simulation_2d_guiding_center_curvilinear_module
