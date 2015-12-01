module sll_m_sim_bsl_gc_2d0v_curv

!2d guiding center cartesian simulation
!contact: Adnane Hamiaz (hamiaz@math.unistra.fr
!         Michel Mehrenberger (mehrenbe@math.unistra.fr)

#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_errors.h"

  use sll_m_fft
  use sll_m_advection_2d_base
  use sll_m_characteristics_1d_explicit_euler
  use sll_m_advection_2d_bsl
  use sll_m_advection_2d_tensor_product
  use sll_m_characteristics_2d_explicit_euler
  use sll_m_characteristics_2d_verlet
  use sll_m_advection_1d_BSL
  use sll_m_advection_1d_CSL_periodic
  use sll_m_characteristics_1d_trapezoid
  use sll_m_reduction
  use sll_m_sim_base
 
  use sll_m_cubic_spline_interpolator_1d
  use sll_m_cubic_spline_interpolator_2d
  use sll_m_hermite_interpolation_2d
  use sll_m_hermite_interpolator_2d
  use sll_m_hermite_interpolation_1d
  use sll_m_hermite_interpolator_1d
  use sll_m_arbitrary_degree_spline_interpolator_2d

  use sll_m_coordinate_transformations_2d
  use sll_m_common_coordinate_transformations
  use sll_m_common_array_initializers
  use sll_m_parallel_array_initializer

  use sll_m_xdmf

#ifdef MUDPACK
 !use sll_m_mudpack_curvilinear
  use sll_m_poisson_2d_mudpack_curvilinear_solver_old
#endif
  use sll_m_poisson_2d_base
  use sll_m_poisson_2d_elliptic_solver, only: new_poisson_2d_elliptic_solver
  use sll_m_general_coordinate_elliptic_solver, only: es_gauss_legendre

  implicit none

  sll_int32, parameter :: SLL_COMPUTE_FIELD_FROM_PHI = 0 
  sll_int32, parameter :: SLL_COMPUTE_FIELD_FROM_PHI_FD = 1 

  
  sll_int32, parameter :: SLL_EULER = 0 
  sll_int32, parameter :: SLL_PREDICTOR_CORRECTOR = 1 
  sll_int32, parameter :: SLL_ADVECTIVE = 0 
  sll_int32, parameter :: SLL_CONSERVATIVE = 1
  sll_int32, parameter :: SLL_X_Y = 0
  sll_int32, parameter :: SLL_ETA1_ETA2 = 1
  type, extends(sll_simulation_base_class) :: &
    sll_simulation_2d_guiding_center_curvilinear

   !geometry
   type(sll_cartesian_mesh_2d), pointer :: mesh_2d
  
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
    class(sll_c_interpolator_1d), pointer   :: A1_interp1d_x1 => null()
    class(sll_c_interpolator_1d), pointer   :: A2_interp1d_x1 => null()
    class(sll_c_interpolator_1d), pointer   :: A1_interp1d_x2 => null()
    class(sll_c_interpolator_1d), pointer   :: A2_interp1d_x2 => null()
    class(sll_c_interpolator_1d), pointer :: f_interp1d_x1 => null()
    class(sll_c_interpolator_1d), pointer :: f_interp1d_x2 => null()
    class(sll_advection_1d_base), pointer    :: advect_1d_x1 => null()
    class(sll_advection_1d_base), pointer    :: advect_1d_x2 => null()
   
   !initial function
   procedure(sll_scalar_initializer_2d), nopass, pointer :: init_func
   sll_int32  :: type_var
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
#ifdef MUDPACK
    type(mudpack_2d) :: poisson2
#endif    
   !time_iterations
   sll_real64 :: dt
   sll_int32  :: num_iterations
   sll_int32  :: freq_diag
   sll_int32  :: freq_diag_time
   character(len=256)      :: thdiag_filename
   character(len=256)      :: thdiagp_filename
   character(len=256)      :: mesh_name
   character(len=256)      :: f_name
   character(len=256)      :: phi_name

   !time_loop
   sll_int32 :: time_loop_case
   sll_int32 :: compute_field_case
   sll_int32 :: fd_degree1
   sll_int32 :: fd_degree2
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

contains

  function new_guiding_center_2d_curvilinear( &
    filename, &
    num_run &
    ) result(sim)
    type(sll_simulation_2d_guiding_center_curvilinear), pointer :: sim 
    character(len=*), intent(in), optional :: filename   
    sll_int32, intent(in), optional :: num_run   
    sll_int32 :: ierr

    
    SLL_ALLOCATE(sim,ierr)
    
    call initialize_guiding_center_2d_curvilinear(sim,filename,num_run)
  
  end function new_guiding_center_2d_curvilinear
  
  subroutine initialize_guiding_center_2d_curvilinear( &
    sim, &
    filename, &
    num_run)
    class(sll_simulation_2d_guiding_center_curvilinear), intent(inout) :: sim
    character(len=*), intent(in), optional :: filename
    sll_int32, intent(in), optional :: num_run
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
    sll_real64 :: alpha3
    sll_real64 :: alpha4
    sll_real64 :: alpha5
    sll_real64 :: num_sides
    
     !initial_function
    character(len=256) :: initial_function_case
    character(len=256) :: type_var
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
    character(len=256) ::  compute_field_case
    sll_int32 :: hermite_degree1
    sll_int32 :: hermite_degree2
    sll_int32 :: fd_degree1
    sll_int32 :: fd_degree2
    
    !poisson
    character(len=256) :: poisson_solver
    sll_int32 :: mudpack_method    
    sll_int32 :: spline_degree_eta1
    sll_int32 :: spline_degree_eta2
    character(len=256) ::  es_control_case
    character(len=256) ::  interp_rho_case
    sll_int32 :: rho_degree1
    sll_int32 :: rho_degree2
    class(sll_interpolator_2d_base), pointer   :: interp_rho
    logical :: precompute_rhs
    logical :: with_constraint
    logical :: with_constraint_loc    
    logical :: zero_mean
    logical :: zero_mean_loc    
    sll_real64 :: eps_penalization

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
    sll_int32 :: ierr
    sll_real64, dimension(4) :: params_mesh
    sll_real64, dimension(9) :: params_mesh_DSG

    character(len=256)      :: str_num_run
    character(len=256)      :: filename_loc
    logical :: feet_inside1 
    logical :: feet_inside2

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
      alpha2  , &
      alpha3  , &
      alpha4  , & 
      alpha5, &
      num_sides

    namelist /initial_function/ &
      type_var, &
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
      f_interp1d_x1_case, &
      f_interp1d_x2_case, &
      charac1d_x1_case, &
      charac1d_x2_case, &
      advect1d_x1_case, &   
      advect1d_x2_case, &
      advection_form, &
      hermite_degree1, &
      hermite_degree2, &
      fd_degree1, &
      fd_degree2, &
      compute_field_case  

    namelist /poisson/ &
      poisson_solver, &
      mudpack_method, &    
      spline_degree_eta1, &
      spline_degree_eta2, &
      es_control_case, &
      precompute_rhs, &
      interp_rho_case, &
      rho_degree1, &
      rho_degree2, &
      with_constraint, &
      eps_penalization, &
      zero_mean    
      
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
    alpha3 = 0._f64
    alpha4 = 0._f64
    alpha5 = 0._f64
    num_sides = 6._f64
    params_mesh_DSG = (/ alpha1, alpha2, alpha3,alpha4,alpha5,eta1_min,eta2_min,eta1_max,eta2_max/)
    !initial function
    type_var = "SLL_X_Y"
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
    compute_field_case = "SLL_COMPUTE_FIELD_FROM_PHI"
    hermite_degree1 = 4
    hermite_degree2 = 4
    fd_degree1 = 4
    fd_degree2 = 4
    
    advection_form = "SLL_ADVECTIVE"
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
    es_control_case = "SLL_SOLVE_ELLIPTIC_SOLVER"
    interp_rho_case = "SLL_CUBIC_SPLINES"
    rho_degree1 = 3
    rho_degree2 = 3
    precompute_rhs = .false.
    with_constraint = .true.
    eps_penalization = 0._f64
    zero_mean = .true.
     
#ifdef MUDPACK
    print *,'#MUDPACK IS ON'
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

    if(present(num_run))then
      !call int2string(num_run, str_num_run)
      write(str_num_run, *) num_run
      str_num_run = adjustl(str_num_run) 
      sim%thdiag_filename = "thdiag_"//trim(str_num_run)//".dat"
      sim%thdiagp_filename = "thdiagp_"//trim(str_num_run)//".dat"
      sim%mesh_name = "curvilinear_"//trim(str_num_run)
      sim%f_name = "f_"//trim(str_num_run)//"_"
      sim%phi_name = "phi_"//trim(str_num_run)//"_"
    else      
      sim%thdiag_filename = "thdiag.dat"
      sim%thdiagp_filename = "thdiagp.dat"
      sim%mesh_name = "curvilinear"
      sim%f_name = "f_"
      sim%phi_name = "phi_"
    endif



    if(present(filename))then
      filename_loc = filename
      filename_loc = adjustl(filename_loc)
      if(present(num_run)) then
        filename_loc = trim(filename)//"_"//trim(str_num_run)
        !filename_loc = adjustl(filename_loc)
        !print *,'filename_loc=',filename_loc
      endif
      open(unit = input_file, file=trim(filename_loc)//'.nml',IOStat=IO_stat)
        if( IO_stat /= 0 ) then
          print *, '#initialize_guiding_center_2d_curvilinear() failed to open file ', &
          trim(filename)//'.nml'
          STOP
        end if
      print *,'#initialization with filename:'
      print *,'#',trim(filename_loc)//'.nml'
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

    sim%fd_degree1 = fd_degree1
    sim%fd_degree2 = fd_degree2

      
    SLL_ALLOCATE(sim%b11(Nc_eta1+1,Nc_eta2+1),ierr)
    SLL_ALLOCATE(sim%b12(Nc_eta1+1,Nc_eta2+1),ierr)
    SLL_ALLOCATE(sim%b21(Nc_eta1+1,Nc_eta2+1),ierr)
    SLL_ALLOCATE(sim%b22(Nc_eta1+1,Nc_eta2+1),ierr)
    SLL_ALLOCATE(sim%b1(Nc_eta1+1,Nc_eta2+1),ierr)
    SLL_ALLOCATE(sim%b2(Nc_eta1+1,Nc_eta2+1),ierr)
    SLL_ALLOCATE(sim%c(Nc_eta1+1,Nc_eta2+1),ierr)



    select case(compute_field_case)
      case ("SLL_COMPUTE_FIELD_FROM_PHI")
        print*,"#compute_field_case = SLL_COMPUTE_FIELD_FROM_PHI " 
        sim%compute_field_case = SLL_COMPUTE_FIELD_FROM_PHI
      case ("SLL_COMPUTE_FIELD_FROM_PHI_FD")
        print*,"#compute_field_case = SLL_COMPUTE_FIELD_FROM_PHI_FD" 
        sim%compute_field_case = SLL_COMPUTE_FIELD_FROM_PHI_FD
      case default
        SLL_ERROR("initialize_analytic_field_2d_curvilinear","bad compute_field_case")
        stop
    end select


    
    select case(bc_eta1_left)
      case ("SLL_PERIODIC")
        print*,"#bc_eta1_left = SLL_PERIODIC" 
        sim%bc_eta1_left = SLL_PERIODIC
      case ("SLL_DIRICHLET")
        print*,"#bc_eta1_left = SLL_DIRICHLET"  
        sim%bc_eta1_left = SLL_DIRICHLET
      case ("SLL_NEUMANN")
        print*,"#bc_eta1_left = SLL_NEUMANN"  
        sim%bc_eta1_left = SLL_NEUMANN
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
      case ("SLL_NEUMANN")
        print*,"#bc_eta1_right = SLL_NEUMANN"  
        sim%bc_eta1_right = SLL_NEUMANN
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
      case ("SLL_NEUMANN")
        print*,"#bc_eta2_left = SLL_NEUMANN"  
        sim%bc_eta2_left = SLL_NEUMANN
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
      case ("SLL_NEUMANN")
        print*,"#bc_eta2_right = SLL_NEUMANN"  
        sim%bc_eta2_right = SLL_NEUMANN
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
!      case ("SLL_DIRICHLET")
!        print*,"#bc_interp2d_eta1 = SLL_DIRICHLET"  
!        sim%bc_interp2d_eta1= SLL_DIRICHLET
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
!      case ("SLL_DIRICHLET")
!        print*,"#bc_interp2d_eta2 = SLL_DIRICHLET"  
!        sim%bc_interp2d_eta2= SLL_DIRICHLET
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
      case ("SLL_CARTESIAN_MESH")
        eta1_max = eta1_max
        eta2_max = eta2_max
      case ("SLL_LANDAU_MESH")
        eta1_max = 2._f64 *sll_pi/kmode_eta1
        eta2_max = 2._f64 *sll_pi/kmode_eta2
      case ("SLL_POLAR_MESH")
        eta2_max = 2._f64*sll_pi
      case ("SLL_POLAR_SHEAR_MESH")
        eta2_max = 2._f64*sll_pi
        params_mesh = (/alpha1,alpha2,0._f64,0._f64/)
      case ("SLL_POLYGONAL_MESH")
        eta2_max = 2._f64*sll_pi
        params_mesh = (/alpha1,num_sides,alpha2,0._f64/)        
      case ("SLL_COLLELA_MESH")  
        eta1_max = 2._f64*sll_pi/kmode_eta1
        eta2_max = 2._f64*sll_pi/kmode_eta2
         !  In collela  mesh params_mesh =( alpha1, alpha2, L1, L2 ) such that :
         !  x1= eta1 + alpha1*sin(2*pi*eta1/L1)*sin(2*pi*eta2/L2)
         params_mesh = (/ alpha1, alpha2, eta1_max-eta1_min, eta2_max-eta2_min/)
      case ("SLL_D_SHAPED_MESH")
        eta1_max = eta1_max
        eta2_max = eta2_max
         !        x1 = alpha1+(alpha2*(2._f64*eta1-1)+alpha3)* &
         !             cos(2._f64*sll_pi*eta2+alpha4*sin(2._f64*sll_p*eta2))
         !        x2 = alpha5*(alpha2*(2._f64*eta1-1)+alpha3)*sin(2._f64*sll_pi*eta2)
         ! Domain: [0,1] X [0,1]
         ! By default the values of the alpha parameters are:
         !             alpha1  = 1.7_f64
         !             alpha2  = 0.074_f64
         !             alpha3  = 0.536_f64
         !             alpha4  = 0.4290421957_f64
         !             alpha5  = 1.66_f64
        params_mesh_DSG = (/ alpha1, alpha2, alpha3,alpha4,alpha5,eta1_min,eta2_min,eta1_max,eta2_max/)
      case default
        print *,'#bad mesh_case',mesh_case
        print *,'#not implemented'
        print *,'#in initialize_guiding_center_2d_curvilinear'
        stop
    end select

        
    sim%mesh_2d => new_cartesian_mesh_2d( &
      Nc_eta1, &
      Nc_eta2, &
      eta1_min , &
      eta1_max , &
      eta2_min , &
      eta2_max ) 
      
    select case (mesh_case)
      case ("SLL_CARTESIAN_MESH")       
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
      case ("SLL_POLAR_SHEAR_MESH") 
        sim%transformation => new_coordinate_transformation_2d_analytic( &
         "analytic_polar_shear_transformation", &
         sim%mesh_2d, &
         polar_shear_x1, &
         polar_shear_x2, &
         polar_shear_jac11, &
         polar_shear_jac12, &
         polar_shear_jac21, &
         polar_shear_jac22, &
         params_mesh  )     
      case ("SLL_POLYGONAL_MESH") 
        sim%transformation => new_coordinate_transformation_2d_analytic( &
         "analytic_polygonal_transformation", &
         sim%mesh_2d, &
         polygonal_x1, &
         polygonal_x2, &
         polygonal_jac11, &
         polygonal_jac12, &
         polygonal_jac21, &
         polygonal_jac22, &
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
       case ("SLL_D_SHAPED_MESH")
        sim%transformation => new_coordinate_transformation_2d_analytic( &
         "analytic_D_SHAPED_transformation", &
         sim%mesh_2d, &
         D_sharped_Geo_x1, &
         D_sharped_Geo_x2, &
         D_sharped_Geo_jac11, &
         D_sharped_Geo_jac12, &
         D_sharped_Geo_jac21, &
         D_sharped_Geo_jac22, &
         params_mesh_DSG  )    
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
        feet_inside1 = .true.
      case ("SLL_CSL_PERIODIC")
        feet_inside1 = .false.
      case default
        print *,'#bad value of advect1d_x1_case'
        stop  
    end select

    select case(advect1d_x2_case)
      case ("SLL_BSL")
        feet_inside2 = .true.
      case ("SLL_CSL_PERIODIC")
        feet_inside2 = .false.
      case default
        print *,'#bad value of advect1d_x2_case'
        stop  
    end select


      
    select case (f_interp2d_case)
      case ("SLL_CUBIC_SPLINES")
        print*,"#f interpolation SLL_CUBIC_SPLINES"
        sim%f_interp2d => new_cubic_spline_interpolator_2d( &
          Nc_eta1+1, &
          Nc_eta2+1, &
          eta1_min, &
          eta1_max, &
          eta2_min, &
          eta2_max, &
          sim%bc_interp2d_eta1, &
          sim%bc_interp2d_eta2, &
          const_eta1_min_slope = 0._f64, & 
          const_eta1_max_slope = 0._f64, &
          const_eta2_min_slope = 0._f64, &
          const_eta2_max_slope = 0._f64 )
      case ("SLL_HERMITE")
        print*,"#f interpolation SLL_HERMITE"
        sim%f_interp2d => new_hermite_interpolator_2d( &
          Nc_eta1+1, &
          Nc_eta2+1, &
          eta1_min, &
          eta1_max, &
          eta2_min, &
          eta2_max, &
          hermite_degree1, &          
          hermite_degree2, &          
          SLL_HERMITE_C0, &
          SLL_HERMITE_C0, &
          SLL_HERMITE_PERIODIC, &
          SLL_HERMITE_PERIODIC)
      case default
        print *,'#bad f_interp2d_case',f_interp2d_case
        print *,'#not implemented'
        print *,'#in initialize_guiding_center_2d_curvilinear'
        stop
    end select


    select case (interp_rho_case)
      case ("SLL_CUBIC_SPLINES")
        print*,"#rhs interpolation SLL_CUBIC_SPLINES"
        interp_rho => new_cubic_spline_interpolator_2d( &
          Nc_eta1+1, &
          Nc_eta2+1, &
          eta1_min, &
          eta1_max, &
          eta2_min, &
          eta2_max, &
          sim%bc_interp2d_eta1, &
          sim%bc_interp2d_eta2, &
          const_eta1_min_slope = 0._f64, & 
          const_eta1_max_slope = 0._f64, &
          const_eta2_min_slope = 0._f64, &
          const_eta2_max_slope = 0._f64 )
      case ("SLL_HERMITE")
        print*,"#rho interpolation SLL_HERMITE"
        interp_rho => new_hermite_interpolator_2d( &
          Nc_eta1+1, &
          Nc_eta2+1, &
          eta1_min, &
          eta1_max, &
          eta2_min, &
          eta2_max, &
          rho_degree1, &          
          rho_degree2, &          
          SLL_HERMITE_C0, &
          SLL_HERMITE_C0, &
          SLL_HERMITE_PERIODIC, &
          SLL_HERMITE_PERIODIC)
      case ("SLL_ARBITRARY_DEGREE_SPLINES")
        print*,"#rho interpolation SLL_ARBITRARY_DEGREE_SPLINES"
        interp_rho => new_arbitrary_degree_spline_interp2d( &
          Nc_eta1+1, &
          Nc_eta2+1, &
          eta1_min, &
          eta1_max, &
          eta2_min, &
          eta2_max, &
          sim%bc_eta1_left, &
          sim%bc_eta1_right, &
          sim%bc_eta2_left, &
          sim%bc_eta2_right, &
          rho_degree1, &
          rho_degree2)          
      case default
        print *,'#bad interp_rho_case',interp_rho_case
        print *,'#not implemented'
        print *,'#in initialize_guiding_center_2d_curvilinear'
        stop
    end select





    select case (A_interp_case)
      case ("SLL_CUBIC_SPLINES")
       print*,"#A1_2d interpolation SLL_CUBIC_SPLINES"
        sim%A1_interp2d => new_cubic_spline_interpolator_2d( &
          Nc_eta1+1, &
          Nc_eta2+1, &
          eta1_min, &
          eta1_max, &
          eta2_min, &
          eta2_max, &
          sim%bc_interp2d_eta1, &
          sim%bc_interp2d_eta2, &
          const_eta1_min_slope = 0._f64, & 
          const_eta1_max_slope = 0._f64, &
          const_eta2_min_slope = 0._f64, &
          const_eta2_max_slope = 0._f64 )
       print*,"#A2_2d interpolation SLL_CUBIC_SPLINES"   
        sim%A2_interp2d => new_cubic_spline_interpolator_2d( &
          Nc_eta1+1, &
          Nc_eta2+1, &
          eta1_min, &
          eta1_max, &
          eta2_min, &
          eta2_max, &
          sim%bc_interp2d_eta1, &
          sim%bc_interp2d_eta2, &  
          const_eta1_min_slope = 0._f64, & 
          const_eta1_max_slope = 0._f64, &
          const_eta2_min_slope = 0._f64, &
          const_eta2_max_slope = 0._f64 )
       print*,"#A1_1d interpolation SLL_CUBIC_SPLINES"   
        sim%A1_interp1d_x1 => new_cubic_spline_interpolator_1d( &
          Nc_eta1+1, &
          eta1_min, &
          eta1_max, &
          sim%bc_interp2d_eta1, &
          slope_left = 0._f64, &
          slope_right = 0._f64)
        sim%A1_interp1d_x2 => new_cubic_spline_interpolator_1d( &
          Nc_eta2+1, &
          eta2_min, &
          eta2_max, &
          sim%bc_interp2d_eta2, &
          slope_left = 0._f64, &
          slope_right = 0._f64)
       print*,"#A2_1d interpolation SLL_CUBIC_SPLINES"     
        sim%A2_interp1d_x1 => new_cubic_spline_interpolator_1d( &
          Nc_eta1+1, &
          eta1_min, &
          eta1_max, &
          sim%bc_interp2d_eta1, &
          slope_left = 0._f64, &
          slope_right = 0._f64)
        sim%A2_interp1d_x2 => new_cubic_spline_interpolator_1d( &
          Nc_eta2+1, &
          eta2_min, &
          eta2_max, &
          sim%bc_interp2d_eta2, &
          slope_left = 0._f64, &
          slope_right = 0._f64)
      case default
        print *,'#bad A_interp_case',A_interp_case
        print *,'#not implemented'
        print *,'#in initialize_guiding_center_2d_curvilinear'
        stop
    end select

    select case (phi_interp2d_case)
      case ("SLL_CUBIC_SPLINES")
      print*,"#phi interpolation SLL_CUBIC_SPLINES"  
        sim%phi_interp2d => new_cubic_spline_interpolator_2d( &
          Nc_eta1+1, &
          Nc_eta2+1, &
          eta1_min, &
          eta1_max, &
          eta2_min, &
          eta2_max, &
          sim%bc_interp2d_eta1, &
          sim%bc_interp2d_eta2, &
          const_eta1_min_slope = 0._f64, & 
          const_eta1_max_slope = 0._f64, &
          const_eta2_min_slope = 0._f64, &
          const_eta2_max_slope = 0._f64 )
      case default
        print *,'#bad phi_interp2d_case',phi_interp2d_case
        print *,'#not implemented'
        print *,'#in initialize_guiding_center_2d_curvilinear'
        stop
    end select

    select case (f_interp1d_x1_case)
      case ("SLL_CUBIC_SPLINES")
        sim%f_interp1d_x1 => new_cubic_spline_interpolator_1d( &
          Nc_eta1+1, &
          eta1_min, &
          eta1_max, &
          sim%bc_interp2d_eta1, &
          slope_left = 0._f64, &
          slope_right = 0._f64)
      case ("SLL_HERMITE")
        sim%f_interp1d_x1 => new_hermite_interpolator_1d( &
          Nc_eta1+1, &
          eta1_min, &
          eta1_max, &
          hermite_degree1, &          
          SLL_HERMITE_1d_C0, &
          sim%bc_interp2d_eta1) 
      case default
        print *,'#bad f_interp1d_x1_case',f_interp1d_x1_case
        print *,'#not implemented'
        print *,'#in initialize_guiding_center_2d_curvilinear'
        stop
    end select


    select case (f_interp1d_x2_case)
      case ("SLL_CUBIC_SPLINES")
        sim%f_interp1d_x2 => new_cubic_spline_interpolator_1d( &
          Nc_eta2+1, &
          eta2_min, &
          eta2_max, &
          sim%bc_interp2d_eta2, &
          slope_left = 0._f64, &
          slope_right = 0._f64)
      case ("SLL_HERMITE")
        sim%f_interp1d_x2 => new_hermite_interpolator_1d( &
          Nc_eta2+1, &
          eta2_min, &
          eta2_max, &
          hermite_degree1, &          
          SLL_HERMITE_1d_C0, &
          sim%bc_interp2d_eta2 )
      case default
        print *,'#bad f_interp1d_x2_case',f_interp1d_x2_case
        print *,'#not implemented'
        print *,'#in initialize_guiding_center_2d_curvilinear'
        stop
    end select


    select case(charac1d_x1_case)
      case ("SLL_EULER")
        sim%charac1d_x1 => new_explicit_euler_1d_charac(&
          Nc_eta1+1, &
          eta_min=eta1_min, &
          eta_max=eta1_max, &
          bc_type= sim%bc_charac2d_eta1, &
          feet_inside = feet_inside1)    
      case ("SLL_TRAPEZOID")
        sim%charac1d_x1 => &
          new_trapezoid_1d_charac(&
          Nc_eta1+1, &
          sim%A1_interp1d_x1, &
          bc_type= sim%bc_charac2d_eta1, &
          eta_min=eta1_min, &
          eta_max=eta1_max, &
          feet_inside = feet_inside1)
      case default
        print *,'#bad charac1d_x1_case',charac1d_x1_case
        print *,'#not implemented'
        print *,'#in initialize_guiding_center_2d_curvilinear'
        stop
    end select

    select case(charac1d_x2_case)
      case ("SLL_EULER")
        sim%charac1d_x2 => new_explicit_euler_1d_charac(&
          Nc_eta2+1, &
          eta_min=eta2_min, &
          eta_max=eta2_max, &
          bc_type= sim%bc_charac2d_eta2, &
          feet_inside = feet_inside2)    
      case ("SLL_TRAPEZOID")
        sim%charac1d_x2 => &
          new_trapezoid_1d_charac(&
          Nc_eta2+1, &
          sim%A2_interp1d_x2, &
          bc_type= sim%bc_charac2d_eta2, &
          eta_min=eta2_min, &
          eta_max=eta2_max, &
          feet_inside = feet_inside2)
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
          Nc_eta1+1, &
          eta_min = eta1_min, &
          eta_max = eta1_max)
      case ("SLL_CSL_PERIODIC")
        sim%advect_1d_x1 => new_CSL_periodic_1d_advector(&
          sim%f_interp1d_x1, &
          sim%charac1d_x1, &
          Nc_eta1+1, &
          eta_min = eta1_min, &
          eta_max = eta1_max, &
          csl_degree = hermite_degree1+1)
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
          Nc_eta2+1, &
          eta_min = eta2_min, &
          eta_max = eta2_max)
      case ("SLL_CSL_PERIODIC")
        sim%advect_1d_x2 => new_CSL_periodic_1d_advector(&
          sim%f_interp1d_x2, &
          sim%charac1d_x2, &
          Nc_eta2+1, &
          eta_min = eta2_min, &
          eta_max = eta2_max, &
          csl_degree = hermite_degree2+1)
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
    
    select case(type_var)
      case ("SLL_X_Y")
        print*,"# type_var = SLL_X_Y" 
        sim%type_var = SLL_X_Y
      case ("SLL_ETA1_ETA2")
        print*,"## tupe_var =  SLL_ETA1_ETA2"  
        sim%type_var = SLL_ETA1_ETA2
      case default
        print *,'#bad type_var',type_var
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
      case("SLL_DIOCOTRON")
        print*,"#f0 = SLL_DIOCOTRON (X,Y)" 
        sim%init_func => sll_diocotron_initializer_2d2
        SLL_ALLOCATE(sim%params(4),ierr)
        sim%params(1) = r_minus
        sim%params(2) = r_plus
        sim%params(3) = eps
        sim%params(4) = kmode_eta2
       case ("SLL_DSG_2D")
        print*,"#f0 = SLL_DSG_2D"  
        sim%init_func => SLL_DSG_2D
        SLL_ALLOCATE(sim%params(4),ierr)
        sim%params(1) = eta1_min
        sim%params(2) = eta2_min
        sim%params(3) = eta1_max
        sim%params(4) = eta2_max 
        case("SLL_DIOCOTRON_ETA1_ETA2")
        print*,"#f0 = SLL_DIOCOTRON (ETA1,ETA2)" 
        sim%init_func => sll_diocotron_initializer_2d
        SLL_ALLOCATE(sim%params(4),ierr)
        sim%params(1) = r_minus
        sim%params(2) = r_plus
        sim%params(3) = eps
        sim%params(4) = kmode_eta2 *2._f64*sll_pi 
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
        if ((sim%bc_eta1_left == SLL_PERIODIC .and. sim%bc_eta1_right == SLL_PERIODIC .and. &
             sim%bc_eta2_left== SLL_PERIODIC .and. sim%bc_eta2_right== SLL_PERIODIC).or. &
             (sim%bc_eta1_left == SLL_NEUMANN .and. sim%bc_eta1_right == SLL_NEUMANN .and. &
             sim%bc_eta2_left== SLL_NEUMANN .and. sim%bc_eta2_right== SLL_NEUMANN)) then
            sim%c   = eps_penalization  ! entre 1e-8 et 1e-7
            with_constraint_loc = with_constraint
            zero_mean_loc = zero_mean
        else
            sim%c   = 0._f64
            with_constraint_loc = .false.
            zero_mean_loc = .false.
        endif    
        
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
         sim%bc_interp2d_eta1, &
         sim%bc_interp2d_eta2, &
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
         sim%c, &
         interp_rho=interp_rho, &
         precompute_rhs=precompute_rhs, &
         with_constraint = with_constraint_loc, &
         zero_mean = zero_mean_loc ) 
      case default
        print *,'#bad poisson_case',poisson_solver
        print *,'#not implemented'
        print *,'#in  initialize_guiding_center_2d_curvilinear'
        stop
    end select

    select case(compute_field_case)
      case ("SLL_COMPUTE_FIELD_FROM_PHI")
        print*,"#compute_field_case = SLL_COMPUTE_FIELD_FROM_PHI " 
        sim%compute_field_case = SLL_COMPUTE_FIELD_FROM_PHI
      case ("SLL_COMPUTE_FIELD_FROM_PHI_FD")
        print*,"#compute_field_case = SLL_COMPUTE_FIELD_FROM_PHI_FD" 
        sim%compute_field_case = SLL_COMPUTE_FIELD_FROM_PHI_FD
      case default
        SLL_ERROR("initialize_guiding_center_2d_curvilinear","bad compute_field_case")
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
    sll_real64,dimension(:,:), pointer :: fone_conserv
    sll_real64,dimension(:,:), pointer :: fone_old
    sll_real64,dimension(:,:), pointer :: fone
    
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
    sll_real64 :: time
    
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
    SLL_ALLOCATE(phi(Nc_eta1+1,Nc_eta2+1),ierr); phi = 0._f64
    SLL_ALLOCATE(A1(Nc_eta1+1,Nc_eta2+1),ierr)
    SLL_ALLOCATE(A2(Nc_eta1+1,Nc_eta2+1),ierr)
    SLL_ALLOCATE(f_conserv(Nc_eta1+1,Nc_eta2+1),ierr)

    SLL_ALLOCATE(fone_conserv(Nc_eta1+1,Nc_eta2+1),ierr)
    SLL_ALLOCATE(fone_old(Nc_eta1+1,Nc_eta2+1),ierr)
    SLL_ALLOCATE(fone(Nc_eta1+1,Nc_eta2+1),ierr)
  

    
    !initialisation of distribution function
    select case(sim%type_var)    
       case(SLL_X_Y)
         do i2=1,Nc_eta2+1
           eta2=eta2_min+real(i2-1,f64)*delta_eta2
           do i1=1,Nc_eta1+1
             eta1=eta1_min+real(i1-1,f64)*delta_eta1
             x1 = sim%transformation%x1(eta1,eta2)
             x2 = sim%transformation%x2(eta1,eta2)
             f(i1,i2) =  sim%init_func(x1,x2,sim%params) 
           end do
         end do 
       case(SLL_ETA1_ETA2)
          do i2=1,Nc_eta2+1
           eta2=eta2_min+real(i2-1,f64)*delta_eta2
           do i1=1,Nc_eta1+1
             eta1=eta1_min+real(i1-1,f64)*delta_eta1
             f(i1,i2) =  sim%init_func(eta1,eta2,sim%params) 
             write(201,*) eta1,eta2,f(i1,i2)
           end do
         end do 
       case default
        print *,'#bad type_var',sim%type_var
        print *,'#not implemented'
        print *,'#in  run_gc2d_curvilinear'
        stop
    end select    

    fone = 1._f64    


    !solve poisson
    call sim%poisson%compute_phi_from_rho(phi, f)

    select case (sim%compute_field_case)
      case (SLL_COMPUTE_FIELD_FROM_PHI)
        call compute_field_from_phi_2d_curvilinear( &
          phi, &
          sim%mesh_2d, &
          sim%transformation, &
          A1, &
          A2, &
          sim%phi_interp2d)
      case (SLL_COMPUTE_FIELD_FROM_PHI_FD)
            call compute_field_from_phi_2d_fd_curvilinear( &
              phi, &
              sim%mesh_2d, &
              sim%transformation, &
              A1, &
              A2, &
              sim%phi_interp2d, &
              sim%fd_degree1, &
              sim%fd_degree2)      
 
        
      case default
        SLL_ERROR("run_af2d_curvilinear","bad value of sim%compute_field_case")
    end select    







    !call compute_field_from_phi_2d_curvilinear(phi,sim%mesh_2d,sim%transformation,A1,A2,sim%phi_interp2d)  
!    do i2=1,Nc_eta2+1
!        eta2=eta2_min+real(i2-1,f64)*delta_eta2
!        do i1=1,Nc_eta1+1
!          eta1=eta1_min+real(i1-1,f64)*delta_eta1
!         write(201,*) eta1,eta2,phi(i1,i2),sin(2.0*sll_pi*eta1)*sin(2.0*sll_pi*eta2)
!        end do
!     end do   
!     stop  
    print *,maxval(phi),minval(phi)
    print *,maxval(A1),minval(A1)
    print *,maxval(A2),maxval(A2)

    thdiagp_id=thdiag_id+1
    call sll_ascii_file_create(sim%thdiag_filename, thdiag_id, ierr)
    !call sll_ascii_file_create('thdiag.dat', thdiag_id, ierr)
    call sll_ascii_file_create(sim%thdiagp_filename, thdiagp_id, ierr)
    iplot = 0


    !initialize mesh
    call sll_plot_curvilinear_init( &
      sim%mesh_2d, &
      sim%transformation, &
      sim%mesh_name )    




    do step=1,nb_step+1
      if(modulo(step-1,sim%freq_diag)==0)then 
         print*,"step= ", step
      endif   
      f_old = f
      fone_old = fone
      
      call sim%poisson%compute_phi_from_rho(phi, f_old) 

    select case (sim%compute_field_case)
      case (SLL_COMPUTE_FIELD_FROM_PHI)
        call compute_field_from_phi_2d_curvilinear( &
          phi, &
          sim%mesh_2d, &
          sim%transformation, &
          A1, &
          A2, &
          sim%phi_interp2d)
      case (SLL_COMPUTE_FIELD_FROM_PHI_FD)
            call compute_field_from_phi_2d_fd_curvilinear( &
              phi, &
              sim%mesh_2d, &
              sim%transformation, &
              A1, &
              A2, &
              sim%phi_interp2d, &
              sim%fd_degree1, &
              sim%fd_degree2)      
 
        
      case default
        SLL_ERROR("run_af2d_curvilinear","bad value of sim%compute_field_case")
    end select    


      !call compute_field_from_phi_2d_curvilinear(phi,sim%mesh_2d,sim%transformation,A1,A2,sim%phi_interp2d)      
      
      
      
      time = real(step-1,f64)*dt
      if(modulo(step-1,sim%freq_diag_time)==0)then
        call time_history_diagnostic_collela( &
          thdiag_id , &    
          step-1, &
          dt, &
          sim%mesh_2d, &
          sim%transformation, &
          f, &
          phi, &
          A1, &
          A2, &
          fone)
        call time_history_diagnostic_polar( &
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
        call sll_plot_f( &
          iplot, &
          f, &  
          Nc_eta1+1, &
          Nc_eta2+1,  &
          sim%f_name, &
          sim%mesh_name, &
          time )    
        call sll_plot_f( &
          iplot, &
          phi, &  
          Nc_eta1+1, &
          Nc_eta2+1,  &
          sim%phi_name, &
          sim%mesh_name, &
          time )    
        !call plot_f_curvilinear(iplot,f,sim%mesh_2d,sim%transformation)
        !call plot_phi_curvilinear(iplot,phi,sim%mesh_2d,sim%transformation)
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


              call compute_f_to_f_conserv(fone_old,fone_conserv,sim%transformation,sim%mesh_2d)
              call sim%advect_2d%advect_2d(A1, A2, dt, fone_conserv, f)
              fone_conserv =fone
              call compute_f_conserv_to_f(fone_conserv,fone,sim%transformation,sim%mesh_2d)
            
            
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
          call sim%poisson%compute_phi_from_rho(phi, f)


    select case (sim%compute_field_case)
      case (SLL_COMPUTE_FIELD_FROM_PHI)
        call compute_field_from_phi_2d_curvilinear( &
          phi, &
          sim%mesh_2d, &
          sim%transformation, &
          A1, &
          A2, &
          sim%phi_interp2d)
      case (SLL_COMPUTE_FIELD_FROM_PHI_FD)
            call compute_field_from_phi_2d_fd_curvilinear( &
              phi, &
              sim%mesh_2d, &
              sim%transformation, &
              A1, &
              A2, &
              sim%phi_interp2d, &
              sim%fd_degree1, &
              sim%fd_degree2)      
 
        
      case default
        SLL_ERROR("run_af2d_curvilinear","bad value of sim%compute_field_case")
    end select    


          
          !call compute_field_from_phi_2d_curvilinear(phi,sim%mesh_2d,sim%transformation,A1,A2,sim%phi_interp2d)      
          
          
          
          
          select case (sim%advection_form)
            case(SLL_CONSERVATIVE)
              
              call compute_f_to_f_conserv(f_old,f_conserv,sim%transformation,sim%mesh_2d)
              call sim%advect_2d%advect_2d(A1, A2, dt, f_conserv, f)
              f_conserv =f
              call compute_f_conserv_to_f(f_conserv,f,sim%transformation,sim%mesh_2d)

              call compute_f_to_f_conserv(fone_old,fone_conserv,sim%transformation,sim%mesh_2d)
              call sim%advect_2d%advect_2d(A1, A2, dt, fone_conserv, fone)
              fone_conserv =fone
              call compute_f_conserv_to_f(fone_conserv,fone,sim%transformation,sim%mesh_2d)

            
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
        A1(i1,i2)=interp2d%interpolate_derivative_eta2(eta1,eta2)/transformation%jacobian(eta1,eta2)
        A2(i1,i2)=-interp2d%interpolate_from_interpolant_derivative_eta1(eta1,eta2)/transformation%jacobian(eta1,eta2)
      end do
    end do
   
    
  end subroutine compute_field_from_phi_2d_curvilinear



subroutine compute_field_from_phi_2d_fd_curvilinear(phi,mesh_2d,transformation,A1,A2,interp2d,d1,d2)
    sll_real64, dimension(:,:), intent(in) :: phi
    sll_real64, dimension(:,:), intent(out) :: A1
    sll_real64, dimension(:,:), intent(out) :: A2
    type(sll_cartesian_mesh_2d), pointer :: mesh_2d
    class(sll_coordinate_transformation_2d_base), pointer :: transformation
    class(sll_interpolator_2d_base), pointer   :: interp2d
    sll_int32, intent(in) :: d1
    sll_int32, intent(in) :: d2
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
    sll_real64, dimension(:), allocatable :: w1
    sll_real64, dimension(:), allocatable :: w2
    sll_int32 :: ierr
    sll_int32 :: r1
    sll_int32 :: s1
    sll_int32 :: r2
    sll_int32 :: s2
    sll_int32 :: ii
    
    Nc_eta1 = mesh_2d%num_cells1
    Nc_eta2 = mesh_2d%num_cells2
    eta1_min = mesh_2d%eta1_min
    eta2_min = mesh_2d%eta2_min
    delta_eta1 = mesh_2d%delta_eta1
    delta_eta2 = mesh_2d%delta_eta2
    
    r1 = -d1/2
    r2 = -d2/2
    s1 = (d1+1)/2
    s2 = (d2+1)/2
    SLL_ALLOCATE(w1(r1:s1),ierr)
    SLL_ALLOCATE(w2(r2:s2),ierr)
     
    call compute_w_hermite(w1,r1,s1)
    call compute_w_hermite(w2,r2,s2)

    A1 = 0._f64
    A2 = 0._f64
    do i2=1,Nc_eta2+1
      eta2=eta2_min+real(i2-1,f64)*delta_eta2
      do i1=1,Nc_eta1+1
        eta1=eta1_min+real(i1-1,f64)*delta_eta1
        A1(i1,i2) = 0._f64
        do ii=r1,s1
          A1(i1,i2) = A1(i1,i2)+w1(ii)*phi(i1,modulo(i2+ii-1+Nc_eta2,Nc_eta2)+1)
        enddo
        A1(i1,i2) = A1(i1,i2)/(delta_eta2)
        A1(i1,i2) = A1(i1,i2)/transformation%jacobian(eta1,eta2)
        A2(i1,i2) = 0._f64
        do ii=r2,s2
          A2(i1,i2) = A2(i1,i2)+w2(ii)*phi(modulo(i1+ii-1+Nc_eta1,Nc_eta1)+1,i2)
        enddo
        A2(i1,i2) = A2(i1,i2)/(delta_eta1)
        A2(i1,i2) = -A2(i1,i2)/transformation%jacobian(eta1,eta2)
        
      end do
    end do
   
    
  end subroutine compute_field_from_phi_2d_fd_curvilinear



  subroutine compute_f_to_f_conserv(f,f_conserv,transformation,mesh_2d)
    sll_real64, dimension(:,:), intent(in) :: f
    sll_real64, dimension(:,:), intent(out) :: f_conserv
    type(sll_cartesian_mesh_2d), pointer :: mesh_2d
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
    type(sll_cartesian_mesh_2d), pointer :: mesh_2d
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
  subroutine time_history_diagnostic_collela( &
    file_id, &    
    step, &
    dt, &
    mesh_2d, &
    transformation,&
    f, &
    phi, &
    A1, &
    A2, &
    fone)
    sll_int32, intent(in) :: file_id
    sll_int32, intent(in) :: step
    sll_real64, intent(in) :: dt
    type(sll_cartesian_mesh_2d), pointer :: mesh_2d
    class(sll_coordinate_transformation_2d_base), pointer :: transformation
    sll_real64, dimension(:,:), intent(in) :: f
    sll_real64, dimension(:,:), intent(in) :: phi
    sll_real64, dimension(:,:), intent(in) :: A1
    sll_real64, dimension(:,:), intent(in) :: A2 
    sll_real64, dimension(:,:), intent(in) :: fone 
    sll_real64 :: mass
    sll_real64 :: linf
    sll_real64 :: l1
    sll_real64 :: l2
    sll_real64 :: e
    
    sll_real64, dimension(:), allocatable :: mass_array
    sll_real64, dimension(:), allocatable  :: l1_array
    sll_real64, dimension(:), allocatable  :: l2_array
    sll_real64, dimension(:), allocatable  :: e_array
    sll_real64, dimension(:), allocatable  :: array
    
    sll_real64 :: eta1
    sll_real64 :: eta2
    sll_real64, dimension(:),allocatable :: data
    sll_real64, dimension(:,:),allocatable :: data_2d
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
    sll_real64 :: val
    sll_real64 :: int_phi 
    sll_real64 :: int_jac 

    
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
    SLL_ALLOCATE(array(Nc_eta2+1),ierr)
    SLL_ALLOCATE(data_2d(Nc_eta1+1,Nc_eta2+1),ierr)
 
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
        eta1 = eta1_min + (i1-1)* delta_eta1
        jac_m  =  transformation%jacobian_matrix(eta1,eta2)
        dphi_eta1 = -A2(i1,i2)* transformation%jacobian(eta1,eta2)
        dphi_eta2 = A1(i1,i2)* transformation%jacobian(eta1,eta2)
        data(i1) = abs( jac_m(2,2)*dphi_eta1 - jac_m(2,1)*dphi_eta2 ) + &
          abs( -jac_m(1,2)*dphi_eta1 + jac_m(1,1)*dphi_eta2 )
        !/abs(transformation%jacobian(eta1,eta2)) 
      enddo
      array(i2) = compute_integral_trapezoid_1d(data, Nc_eta1+1, delta_eta1)


      do i1=1,Nc_eta1+1
       linf = max(linf,abs(f(i1,i2)))
      enddo
         
    enddo     

    mass = compute_integral_trapezoid_1d(mass_array, Nc_eta2+1, delta_eta2)
    l1 = compute_integral_trapezoid_1d(l1_array, Nc_eta2+1, delta_eta2)
    l2 = compute_integral_trapezoid_1d(l2_array, Nc_eta2+1, delta_eta2)
    l2 = sqrt(l2)
    e = compute_integral_trapezoid_1d(e_array, Nc_eta2+1, delta_eta2)
    val = compute_integral_trapezoid_1d(array, Nc_eta2+1, delta_eta2)


    !new diag
    do i2 = 1, Nc_eta2+1
      eta2 = eta2_min + (i2-1)* delta_eta2 
      do i1=1,Nc_eta1+1
        eta1 = eta1_min + (i1-1)* delta_eta1
        data_2d(i1,i2) = phi(i1,i2)*abs(transformation%jacobian(eta1,eta2))
      enddo
    enddo 
    
    int_phi = compute_integral_trapezoid_2d( &
      data_2d, &
      Nc_eta1+1, &
      Nc_eta2+1, &
      delta_eta1, &
      delta_eta2)

    do i2 = 1, Nc_eta2+1
      eta2 = eta2_min + (i2-1)* delta_eta2 
      do i1=1,Nc_eta1+1
        eta1 = eta1_min + (i1-1)* delta_eta1
        data_2d(i1,i2) = abs(transformation%jacobian(eta1,eta2))
      enddo
    enddo 


    int_jac = compute_integral_trapezoid_2d( &
      data_2d, &
      Nc_eta1+1, &
      Nc_eta2+1, &
      delta_eta1, &
      delta_eta2)


    do i2 = 1, Nc_eta2+1
      do i1=1,Nc_eta1+1
        data_2d(i1,i2) = phi(i1,i2) - int_phi/int_jac
      enddo
    enddo 


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
      maxval(abs(phi(1:Nc_eta1+1,1:Nc_eta2+1))), &
      val, &
      maxval(abs(fone(1:Nc_eta1+1,1:Nc_eta2+1)-1._f64)), &
      maxval(abs(data_2d)), &
      int_phi, &
      int_jac
   
    
  end subroutine time_history_diagnostic_collela

  subroutine time_history_diagnostic_polar( &
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
    type(sll_cartesian_mesh_2d), pointer :: mesh_2d
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

    
  end subroutine time_history_diagnostic_polar
  
  subroutine time_history_diagnostic_gc3( &
    file_id, &    
    step, &
    dt, &
    mesh_2d, &
    f, &
    phi, &
    A1, &
    A2, &
    params)
    sll_int32, intent(in) :: file_id
    sll_int32, intent(in) :: step
    sll_real64, intent(in) :: dt
    type(sll_cartesian_mesh_2d), pointer :: mesh_2d
    sll_real64, dimension(:,:), intent(in) :: f
    sll_real64, dimension(:,:), intent(in) :: phi
    sll_real64, dimension(:,:), intent(in) :: A1
    sll_real64, dimension(:,:), intent(in) :: A2
    sll_real64 :: params(9)
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
    sll_real64 :: alpha1,alpha2,alpha3,alpha4,alpha5
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
    alpha1 = params(5)
    alpha2 = params(6)
    alpha3 = params(7) 
    alpha4 = params(8) 
    alpha5 = params(9) 
   
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

    
  end subroutine time_history_diagnostic_gc3
#ifndef NOHDF5
!*********************
!*********************

  !---------------------------------------------------
  ! Save the mesh structure
  !---------------------------------------------------
  subroutine plot_f_curvilinear(iplot,f,mesh_2d,transf)
    use sll_m_xdmf
    use sll_m_hdf5_io_serial
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
          eta1 = eta1_min+real(i-1,f64)*deta1
          eta2 = eta2_min+real(j-1,f64)*deta2
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


  subroutine plot_phi_curvilinear(iplot,f,mesh_2d,transf)
    use sll_m_xdmf
    use sll_m_hdf5_io_serial
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
          eta1 = eta1_min+real(i-1,f64)*deta1
          eta2 = eta2_min+real(j-1,f64)*deta2
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
    call sll_xdmf_open("phi"//cplot//".xmf","curvilinear_mesh", &
      nnodes_x1,nnodes_x2,file_id,error)
    call sll_xdmf_write_array("phi"//cplot,f,"values", &
      error,file_id,"Node")
    call sll_xdmf_close(file_id,error)
  end subroutine plot_phi_curvilinear


#endif

subroutine sll_DSG( eta1_min,eta1_max, eta2_min,eta2_max,n_eta1,n_eta2, f ) 
    sll_real64, intent(in)   :: eta1_min,eta1_max
    sll_real64, intent(in)   :: eta2_min,eta2_max  
    sll_int32, intent(in)    :: n_eta1,n_eta2  
    !type(sll_arbitrary_degree_spline_interpolator_2d)   :: a11_interp
    !type(sll_arbitrary_degree_spline_interpolator_2d)   :: a22_interp
    !type(sll_arbitrary_degree_spline_interpolator_2d)   :: a12_interp
    type(sll_cubic_spline_interpolator_2d)   :: a11_interp
    type(sll_cubic_spline_interpolator_2d)   :: a22_interp
    type(sll_cubic_spline_interpolator_2d)   :: a12_interp
    sll_real64,dimension(:,:),allocatable :: cxx_array
    sll_real64,dimension(:,:),allocatable :: cyy_array
    sll_real64,dimension(:,:),allocatable :: cxy_array
    sll_real64,dimension(:,:),allocatable :: cx_array
    sll_real64,dimension(:,:),allocatable :: cy_array
    sll_real64,dimension(:,:),intent(out)  :: f
    sll_real64  :: pi2
    sll_real64  :: eta1
    sll_real64  :: eta2
    sll_real64  :: delta1
    sll_real64  :: delta2
    sll_real64  :: jac
    sll_real64  :: jac_11
    sll_real64  :: jac_12
    sll_real64  :: jac_21 
    sll_real64  :: jac_22
    sll_real64  :: dphi_eta1
    sll_real64  :: dphi_eta2
    sll_real64  :: ddphi_eta11
    sll_real64  :: ddphi_eta22
    sll_real64  :: ddphi_eta12
    sll_real64  :: rho
    sll_int32   :: i,j,spline_degree_eta1, spline_degree_eta2
    sll_real64, dimension(:), allocatable      :: eta1_pos
    sll_real64, dimension(:), allocatable      :: eta2_pos
    
    allocate(cxx_array(n_eta1+1,n_eta2+1)) 
    allocate(cyy_array(n_eta1+1,n_eta2+1)) 
    allocate(cxy_array(n_eta1+1,n_eta2+1)) 
    allocate(cx_array(n_eta1+1,n_eta2+1)) 
    allocate(cy_array(n_eta1+1,n_eta2+1))
    allocate(eta1_pos(n_eta1+1))
    allocate(eta2_pos(n_eta2+1))
    pi2 = 2.0_f64*sll_pi
    
    delta1 = (eta1_max - eta1_min)/n_eta1
    delta2 = (eta2_max - eta2_min)/n_eta2
    do j= 1, n_eta2+1
      eta2= eta2_min + (j-1)*delta2
      do i= 1, n_eta1+1
         eta1= eta1_min + (i-1)*delta1
         jac_11 = 0.148_f64*cos(pi2*eta2 + 0.4290421957_f64*sin(pi2*eta2))
         jac_12 = -(0.148*eta1 + 0.462_f64)*sin(pi2*eta2 + 0.4290421957_f64* &
              sin(pi2*eta2))*(pi2 + 0.8580843914_f64*cos(pi2*eta2)*sll_pi)
         jac_21 = 0.24568_f64*sin(pi2*eta2) 
         jac_22 = 3.32_f64*sll_pi*(0.148_f64*eta1 + 0.462_f64)*cos(pi2*eta2)
         jac = jac_11*jac_22 - jac_12*jac_21   
 
         cxx_array(i,j) = (jac_12*jac_12 + jac_22*jac_22)/jac
         cyy_array(i,j) = (jac_11*jac_11 + jac_21*jac_21)/jac
         cxy_array(i,j) =-2._f64*(jac_11*jac_12 + jac_21*jac_22)/jac
         
         eta1_pos(i) = eta1
         eta2_pos(j) = eta2
       enddo
     enddo  
     
   spline_degree_eta1 = 3
   spline_degree_eta2 = 3  
   call a11_interp%initialize( &
         n_eta1+1, &
         n_eta2+1, &
         eta1_min, &
         eta1_max, &
         eta2_min, &
         eta2_max, &
         SLL_DIRICHLET, &
         SLL_PERIODIC )            
   call a22_interp%initialize( &
         n_eta1+1, &
         n_eta2+1, &
         eta1_min, &
         eta1_max, &
         eta2_min, &
         eta2_max, &
         SLL_DIRICHLET, &
         SLL_PERIODIC )    
           
          
    call a12_interp%initialize( &
         n_eta1+1, &
         n_eta2+1, &
         eta1_min, &
         eta1_max, &
         eta2_min, &
         eta2_max, &
         SLL_DIRICHLET, &
         SLL_PERIODIC )   
    call a11_interp%compute_interpolants(cxx_array, &
       eta1_pos,&
       n_eta1+1,&
       eta2_pos,&
       n_eta2+1)   
    call a22_interp%compute_interpolants(cyy_array,&
       eta1_pos,&
       n_eta1+1,&
       eta2_pos,&
       n_eta2+1)
    call a12_interp%compute_interpolants(cxy_array/2., &
       eta1_pos,&
       n_eta1+1,&
       eta2_pos,&
       n_eta2+1)
    do j= 1, n_eta2+1
      eta2= eta2_min + (j-1)*delta2
      do i= 1, n_eta1+1 
         eta1= eta1_min + (i-1)*delta1 
         jac_11 = 0.148_f64*cos(pi2*eta2 + 0.4290421957_f64*sin(pi2*eta2))
         jac_12 = -(0.148*eta1 + 0.462_f64)*sin(pi2*eta2 + 0.4290421957_f64* &
                  sin(pi2*eta2))*(pi2 + 0.8580843914_f64*cos(pi2*eta2)*sll_pi)
         jac_21 = 0.24568_f64*sin(pi2*eta2) 
         jac_22 = 3.32_f64*sll_pi*(0.148_f64*eta1 + 0.462_f64)*cos(pi2*eta2)
         jac    = jac_11*jac_22 - jac_12*jac_21
    
         dphi_eta1   = 4._f64*(1._f64-eta1)*(1._f64+0.1_f64*sin(4._f64*pi2*eta2))
         dphi_eta2   = 4._f64*eta1*(1._f64-eta1)*(pi2*0.4*cos(4._f64*pi2*eta2))
         ddphi_eta11 = -4._f64*(1._f64+0.1*sin(4._f64*pi2*eta2))
         ddphi_eta22 = -4._f64*eta1*(1._f64-eta1)*(pi2*pi2*1.6_f64*sin(4._f64*pi2*eta2))
         ddphi_eta12 = 4._f64*(1._f64-eta1)*(pi2*0.4_f64*cos(4._f64*pi2*eta2))   
                  
         cx_array(i,j)  = (a11_interp%interpolate_from_interpolant_derivative_eta1(eta1,eta2)+ &
                           a12_interp%interpolate_derivative_eta2(eta1,eta2))
         cy_array(i,j)  = (a22_interp%interpolate_derivative_eta2(eta1,eta2)+ &
                           a12_interp%interpolate_from_interpolant_derivative_eta1(eta1,eta2))  
         write(202,*) eta1,eta2,cx_array(i,j),cy_array(i,j)                                     
         rho = (cxx_array(i,j)*ddphi_eta11 + cyy_array(i,j)*ddphi_eta22 + &
                 cxy_array(i,j)*ddphi_eta12 + cx_array(i,j) *dphi_eta1   + &
                 cy_array(i,j) *dphi_eta2)
         f(i,j) = rho/jac
     enddo
    enddo 
  end subroutine sll_DSG 


  subroutine delete_guiding_center_2d_curvilinear( sim )

    class(sll_simulation_2d_guiding_center_curvilinear) :: sim
            
  end subroutine delete_guiding_center_2d_curvilinear


  function compute_integral_trapezoid_2d( &
    input, &
    Npts1, &
    Npts2, &
    delta1, &
    delta2) &
    result (res)
    sll_real64, dimension(:,:), intent(in) :: input
    sll_int32, intent(in) :: Npts1
    sll_int32, intent(in) :: Npts2
    sll_real64, intent(in) :: delta1
    sll_real64, intent(in) :: delta2
    sll_real64 :: res
    
    sll_int32 :: i
    sll_real64, dimension(:), allocatable :: array
    sll_int32 :: ierr
    
    
    if(size(input,1)<Npts1)then
      print *,'size(input,1)=',size(input,1)
      print *,'Npts1=',Npts1
      SLL_ERROR('compute_integral_trapezoid_2d&
      &','bad size1')
    endif
    if(size(input,2)<Npts2)then
      print *,'size(input,2)=',size(input,2)
      print *,'Npts2=',Npts2
      SLL_ERROR('compute_integral_trapezoid_2d&
      &','bad size2')
    endif

    SLL_ALLOCATE(array(Npts2),ierr)    

    do i=1,Npts2
      array(i) = compute_integral_trapezoid_1d( &
        input(1:Npts1,i),&
        Npts1, &
        delta1)
    enddo    
    res = compute_integral_trapezoid_1d( &
      array,&
      Npts2, &
      delta2)
  end function compute_integral_trapezoid_2d

  subroutine sll_plot_curvilinear_init( &
    mesh_2d, &
    transf, &
    mesh_name )
    
    type(sll_cartesian_mesh_2d), pointer :: mesh_2d
    class(sll_coordinate_transformation_2d_base), pointer :: transf
    character(len=*), intent(in) :: mesh_name 
    
    sll_real64, allocatable :: x1(:,:)
    sll_real64, allocatable :: x2(:,:)
    sll_real64, allocatable :: f(:,:)
    sll_int32 :: ierr
    sll_int32 :: i
    sll_int32 :: j
    sll_int32 :: num_pts1
    sll_int32 :: num_pts2
    sll_real64 :: eta1_min
    sll_real64 :: eta1_max
    sll_real64 :: eta2_min
    sll_real64 :: eta2_max
    sll_real64 :: delta1
    sll_real64 :: delta2
    sll_real64 :: eta1
    sll_real64 :: eta2

    num_pts1 = mesh_2d%num_cells1+1
    num_pts2 = mesh_2d%num_cells2+1
    eta1_min = mesh_2d%eta1_min
    eta1_max = mesh_2d%eta1_max
    eta2_min = mesh_2d%eta2_min
    eta2_max = mesh_2d%eta2_max
    delta1 = mesh_2d%delta_eta1
    delta2 = mesh_2d%delta_eta2
    
    SLL_ALLOCATE(x1(num_pts1,num_pts2),ierr)
    SLL_ALLOCATE(x2(num_pts1,num_pts2),ierr)
    SLL_ALLOCATE(f(num_pts1,num_pts2),ierr)
    
    f = 0._f64
    
    do j=1,num_pts2
      eta2 = eta2_min+real(j-1,f64)*delta2
      do i=1,num_pts1
        eta1 = eta1_min+real(i-1,f64)*delta1
        x1(i,j) = transf%x1(eta1,eta2)
        x2(i,j) = transf%x2(eta1,eta2) 
      enddo
    enddo
    call sll_plot_f( &
      0, &
      f, &  
      num_pts1, &
      num_pts2,  &
      "f", & !dummy (for sll_plt_f, we should be able to
      !initialize only the mesh TODO)
      mesh_name, &
      0._f64, &
      x1, &
      x2)    
        
  end subroutine sll_plot_curvilinear_init


end module sll_m_sim_bsl_gc_2d0v_curv
