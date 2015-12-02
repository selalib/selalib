module sll_m_sim_bsl_gc_2d0v_polar

!2d guiding center polar simulation
!related to simulation_2d_guiding_center_curvilinear.F90
!but here geometry and test are specifically polar

!see ../selalib/src/simulation/gcsim2d_polar_input.nml
!for example of use

!contact: Michel Mehrenberger (mehrenbe@math.unistra.fr)
!         Adnane Hamiaz (hamiaz@math.unistra.fr)
#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"

  use sll_m_cartesian_meshes  
  use sll_m_advection_1d_periodic
  use sll_m_advection_2d_bsl
  use sll_m_characteristics_2d_explicit_euler
  use sll_m_characteristics_2d_verlet
  !use sll_poisson_2d_periodic  
!  use sll_m_fft
  use sll_m_reduction
  use sll_m_sim_base
  use sll_m_hermite_interpolator_2d
  use sll_m_cubic_spline_interpolator_1d
  use sll_m_cubic_spline_interpolator_2d
  use sll_m_coordinate_transformation_2d_base
  use sll_m_coordinate_transformations_2d
  use sll_m_common_coordinate_transformations
  use sll_m_common_array_initializers
  use sll_m_parallel_array_initializer
#ifdef MUDPACK
  use sll_m_poisson_2d_mudpack
  use sll_m_poisson_2d_mudpack_curvilinear_solver_old
#endif
!  use sll_m_poisson_2d_base
  use sll_m_poisson_2d_polar_wrapper
  use sll_m_general_coordinate_elliptic_solver, only: es_gauss_legendre
  use sll_m_poisson_2d_elliptic_solver, only: new_poisson_2d_elliptic_solver

  use sll_m_boundary_condition_descriptors
  use sll_m_hermite_interpolation_2d
  use sll_m_xdmf
  
  !use sll_m_parallel_array_initializer

  implicit none

!#define OLD_POISSON  
!#define NEW_POISSON  
  
  sll_int32, parameter :: SLL_EULER = 0 
  sll_int32, parameter :: SLL_PREDICTOR_CORRECTOR = 1 
  sll_int32, parameter :: SLL_PHI_FROM_RHO = 0
  sll_int32, parameter :: SLL_E_FROM_RHO = 1


  type, extends(sll_simulation_base_class) :: &
    sll_simulation_2d_guiding_center_polar


   !geometry
   type(sll_cartesian_mesh_2d), pointer :: mesh_2d
   class(sll_coordinate_transformation_2d_base), pointer :: transformation


   !initial function
   procedure(sll_scalar_initializer_2d), nopass, pointer :: init_func
   sll_real64, dimension(:), pointer :: params
      
   !advector
   class(sll_advection_2d_base), pointer    :: advect_2d
   
   !interpolator for derivatives
   class(sll_c_interpolator_2d), pointer   :: phi_interp2d

   
   !poisson solver
   class(sll_poisson_2d_base), pointer   :: poisson
   sll_int32 :: poisson_case

   !time_iterations
   sll_real64 :: dt
   sll_int32  :: num_iterations
   sll_int32  :: freq_diag
   sll_int32  :: freq_diag_time
   character(len=256)      :: thdiag_filename
   character(len=256)      :: mesh_name
   character(len=256)      :: f_name
   character(len=256)      :: phi_name

   !time_loop
   sll_int32 :: time_loop_case

       
  contains
    procedure, pass(sim) :: run => run_gc2d_polar
    procedure, pass(sim) :: init_from_file => init_fake
     
  end type sll_simulation_2d_guiding_center_polar

contains

  function new_guiding_center_2d_polar(filename,num_run) result(sim)
    type(sll_simulation_2d_guiding_center_polar), pointer :: sim    
    character(len=*), intent(in), optional :: filename
    sll_int32, intent(in), optional :: num_run
    sll_int32 :: ierr
    
    SLL_ALLOCATE(sim,ierr)
    
    call initialize_guiding_center_2d_polar(sim,filename,num_run)
    
  
  
  end function new_guiding_center_2d_polar
  
  subroutine initialize_guiding_center_2d_polar(sim,filename,num_run)
    class(sll_simulation_2d_guiding_center_polar), intent(inout) :: sim
    
    character(len=*), intent(in), optional :: filename
    sll_int32, intent(in), optional :: num_run
    sll_int32             :: IO_stat
    sll_int32, parameter  :: input_file = 99
    
    !geometry
    character(len=256) :: mesh_case
    sll_int32 :: num_cells_x1
    sll_real64 :: r_min
    sll_real64 :: r_max
    sll_int32 :: num_cells_x2
    
    !initial_function
    character(len=256) :: initial_function_case
    sll_real64 :: r_minus
    sll_real64 :: r_plus
    sll_real64 :: kmode_x2
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
    sll_int32 :: hermite_degree_eta1 
    sll_int32 :: hermite_degree_eta2 
 
    !poisson
    character(len=256) :: poisson_case
    character(len=256) :: poisson_solver
    !character(len=256) :: mudpack_method    
    sll_int32 :: spline_degree_eta1
    sll_int32 :: spline_degree_eta2
    character(len=256) :: bc_min_case
    character(len=256) :: bc_max_case
    sll_int32 :: bc_min
    sll_int32 :: bc_max
    
    !local variables
    sll_int32 :: Nc_x1
    sll_int32 :: Nc_x2
    sll_real64 :: x1_min
    sll_real64 :: x1_max     
    sll_real64 :: x2_min
    sll_real64 :: x2_max     
    class(sll_c_interpolator_2d), pointer :: f_interp2d
    class(sll_c_interpolator_2d), pointer :: phi_interp2d
    class(sll_characteristics_2d_base), pointer :: charac2d
    class(sll_c_interpolator_2d), pointer   :: A1_interp2d
    class(sll_c_interpolator_2d), pointer   :: A2_interp2d
    class(sll_c_interpolator_1d), pointer   :: A1_interp1d_x1
    class(sll_c_interpolator_1d), pointer   :: A2_interp1d_x1
    sll_real64, dimension(:,:), pointer :: b11
    sll_real64, dimension(:,:), pointer :: b12
    sll_real64, dimension(:,:), pointer :: b21
    sll_real64, dimension(:,:), pointer :: b22
    sll_real64, dimension(:,:), pointer :: b1
    sll_real64, dimension(:,:), pointer :: b2
    sll_real64, dimension(:,:), pointer :: c
    class(sll_coordinate_transformation_2d_base), pointer :: transformation
!    sll_real64, dimension(:,:), allocatable :: cxx_2d
!    sll_real64, dimension(:,:), allocatable :: cxy_2d
!    sll_real64, dimension(:,:), allocatable :: cyy_2d
!    sll_real64, dimension(:,:), allocatable :: cx_2d
!    sll_real64, dimension(:,:), allocatable :: cy_2d
!    sll_real64, dimension(:,:), allocatable :: ce_2d
    sll_int32 :: ierr
    character(len=256)      :: str_num_run
    character(len=256)      :: filename_loc



    namelist /geometry/ &
      mesh_case, &
      num_cells_x1, &
      r_min, &
      r_max, &
      num_cells_x2

    namelist /initial_function/ &
      initial_function_case, &
      r_minus, &
      r_plus, &
      kmode_x2, &
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
      hermite_degree_eta1, &
      hermite_degree_eta2
      

    namelist /poisson/ &
      poisson_case, &
      poisson_solver, &
      spline_degree_eta1, &
      spline_degree_eta2, &
      bc_min_case, &
      bc_max_case    


    !set default parameters
    
    !geometry
    mesh_case="SLL_POLAR_MESH"
    num_cells_x1 = 32
    r_min = 1.0_f64
    r_max = 10._f64
    num_cells_x2 = 32
    
    !initial function
    initial_function_case="SLL_DIOCOTRON"
    r_minus = 4._f64
    r_plus = 5._f64
    kmode_x2 = 3._f64
    eps = 1.e-6_f64
    bc_min_case = "SLL_NEUMANN_MODE0"
    bc_max_case = "SLL_DIRICHLET"
    
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
    hermite_degree_eta1 = 9
    hermite_degree_eta2 = 9
    
    !poisson
    poisson_case = "SLL_PHI_FROM_RHO"
    !poisson_case = "SLL_E_FROM_RHO"    
    poisson_solver = "SLL_POLAR_FFT"  
    !poisson_solver = "SLL_ELLIPTIC_FINITE_ELEMENT_SOLVER" !use with "SLL_PHI_FROM_RHO"
    !poisson_solver = "SLL_MUDPACK_CURVILINEAR"   !use with "SLL_PHI_FROM_RHO"    
    spline_degree_eta1 = 3
    spline_degree_eta2 = 3
    bc_min_case = "SLL_NEUMANN_MODE_0"    
    bc_max_case = "SLL_DIRICHLET"    

    if(present(num_run))then
      write(str_num_run, *) num_run
      str_num_run = adjustl(str_num_run) 
      sim%thdiag_filename = "thdiag_"//trim(str_num_run)//".dat"
      sim%mesh_name = "polar_"//trim(str_num_run)
      sim%f_name = "f_"//trim(str_num_run)//"_"
      sim%phi_name = "phi_"//trim(str_num_run)//"_"
    else      
      sim%thdiag_filename = "thdiag.dat"
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
          print *, '#initialize_guiding_center_2d_polar() failed to open file ', &
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
      close(input_file)
    else
      print *,'#initialization with default parameters'    
    endif




    select case (mesh_case)
      case ("SLL_POLAR_MESH")
        x1_min = r_min
        x1_max = r_max
        x2_min = 0._f64
        x2_max = 2._f64*sll_pi
      case default
        print *,'#bad mesh_case',mesh_case
        print *,'#not implemented'
        print *,'#in initialize_guiding_center_2d_polar'
        stop
    end select

    Nc_x1 = num_cells_x1
    Nc_x2 = num_cells_x2
        
    sim%dt = dt
    sim%num_iterations = number_iterations
    sim%freq_diag = freq_diag
    sim%freq_diag_time = freq_diag_time

    sim%mesh_2d => new_cartesian_mesh_2d( &
      Nc_x1, &
      Nc_x2, &
      eta1_min = x1_min, &
      eta1_max = x1_max, &
      eta2_min = x2_min, &
      eta2_max = x2_max)      
      
      
      
    select case (f_interp2d_case)
      case ("SLL_CUBIC_SPLINES")
        f_interp2d => new_cubic_spline_interpolator_2d( &
          Nc_x1+1, &
          Nc_x2+1, &
          x1_min, &
          x1_max, &
          x2_min, &
          x2_max, &
          SLL_HERMITE, &
          SLL_PERIODIC, &
          const_eta1_min_slope = 0._f64, & !to prevent problem on the boundary
          const_eta1_max_slope = 0._f64)
      case ("SLL_HERMITE")
        f_interp2d => new_hermite_interpolator_2d( &
          Nc_x1+1, &
          Nc_x2+1, &
          x1_min, &
          x1_max, &
          x2_min, &
          x2_max, &
          hermite_degree_eta1, &          
          hermite_degree_eta2, &          
          SLL_HERMITE_C0, &
          SLL_HERMITE_C0, &
          SLL_HERMITE_DIRICHLET, &
          SLL_HERMITE_PERIODIC)
      case default
        print *,'#bad f_interp2d_case',f_interp2d_case
        print *,'#not implemented'
        print *,'#in initialize_guiding_center_2d_polar'
        stop
    end select




    select case (A_interp_case)
      case ("SLL_CUBIC_SPLINES")
        A1_interp2d => new_cubic_spline_interpolator_2d( &
          Nc_x1+1, &
          Nc_x2+1, &
          x1_min, &
          x1_max, &
          x2_min, &
          x2_max, &
          SLL_HERMITE, &
          SLL_PERIODIC, &
          const_eta1_min_slope = 0._f64, & !to prevent problem on the boundary
          const_eta1_max_slope = 0._f64)
        A2_interp2d => new_cubic_spline_interpolator_2d( &
          Nc_x1+1, &
          Nc_x2+1, &
          x1_min, &
          x1_max, &
          x2_min, &
          x2_max, &
          SLL_HERMITE, &
          SLL_PERIODIC)  
        A1_interp1d_x1 => new_cubic_spline_interpolator_1d( &
          Nc_x1+1, &
          x1_min, &
          x1_max, &
          SLL_HERMITE)
        A2_interp1d_x1 => new_cubic_spline_interpolator_1d( &
          Nc_x1+1, &
          x1_min, &
          x1_max, &
          SLL_HERMITE, &
          slope_left = 0._f64, &
          slope_right = 0._f64)
      case default
        print *,'#bad A_interp_case',A_interp_case
        print *,'#not implemented'
        print *,'#in initialize_guiding_center_2d_polar'
        stop
    end select

    select case (phi_interp2d_case)
      case ("SLL_CUBIC_SPLINES")
        phi_interp2d => new_cubic_spline_interpolator_2d( &
          Nc_x1+1, &
          Nc_x2+1, &
          x1_min, &
          x1_max, &
          x2_min, &
          x2_max, &
          SLL_HERMITE, &
          SLL_PERIODIC, &
          const_eta1_min_slope = 0._f64, & !to prevent problem on the boundary
          const_eta1_max_slope = 0._f64)         
      case default
        print *,'#bad phi_interp2d_case',phi_interp2d_case
        print *,'#not implemented'
        print *,'#in initialize_guiding_center_2d_polar'
        stop
    end select


    select case(charac2d_case)
      case ("SLL_EULER")
        charac2d => new_explicit_euler_2d_charac(&
          Nc_x1+1, &
          Nc_x2+1, &
          eta1_min=x1_min, &
          eta1_max=x1_max, &
          eta2_min=x2_min, &
          eta2_max=x2_max, &
          bc_type_1=SLL_SET_TO_LIMIT, &
          bc_type_2=SLL_PERIODIC)    
      case ("SLL_VERLET")      
        charac2d => new_verlet_2d_charac(&
          Nc_x1+1, &
          Nc_x2+1, &
          A1_interp2d, &
          A2_interp2d, &
          A1_interp1d_x1, &
          A2_interp1d_x1, &
          bc_type_1=SLL_SET_TO_LIMIT, &
          bc_type_2=SLL_PERIODIC, &
          eta1_min=x1_min, &
          eta1_max=x1_max, &
          eta2_min=x2_min, &
          eta2_max=x2_max )
      case default
        print *,'#bad charac2d_case',charac2d_case
        print *,'#not implemented'
        print *,'#in initialize_guiding_center_2d_polar'
        stop
    end select

  
    sim%phi_interp2d => phi_interp2d

    select case(advect2d_case)
      case ("SLL_BSL")
        sim%advect_2d => new_BSL_2d_advector(&
          f_interp2d, &
          charac2d, &
          Nc_x1+1, &
          Nc_x2+1, &
          eta1_min = x1_min, &
          eta1_max = x1_max, &
          eta2_min = x2_min, &
          eta2_max = x2_max)
      case default
        print *,'#bad advect_case',advect2d_case
        print *,'#not implemented'
        print *,'#in initialize_guiding_center_2d_polar'
        stop
    end select
    
    
    select case(initial_function_case)
      case ("SLL_DIOCOTRON")
        sim%init_func => sll_diocotron_initializer_2d
        SLL_ALLOCATE(sim%params(4),ierr)
        sim%params(1) = r_minus
        sim%params(2) = r_plus
        sim%params(3) = eps
        sim%params(4) = kmode_x2
      case default
        print *,'#bad initial_function_case',initial_function_case
        print *,'#not implemented'
        print *,'#in initialize_guiding_center_2d_polar'
        stop
    end select
    
    
    !time_loop
    select case(time_loop_case)
      case ("SLL_EULER")
        sim%time_loop_case = SLL_EULER
      case ("SLL_PREDICTOR_CORRECTOR")
        sim%time_loop_case = SLL_PREDICTOR_CORRECTOR
      case default
        print *,'#bad time_loop_case',time_loop_case
        print *,'#not implemented'
        print *,'#in initialize_guiding_center_2d_polar'
        stop
    end select
    
    
    select case(poisson_case)    
      case ("SLL_PHI_FROM_RHO")     
        sim%poisson_case = SLL_PHI_FROM_RHO
      !case ("SLL_E_FROM_RHO")
      !  sim%poisson_case = SLL_E_FROM_RHO     
      case default
        print *,'#bad poisson_case',poisson_case
        print *,'#not implemented'
        print *,'#in initialize_guiding_center_2d_cartesian'
        stop
    end select
    select case(bc_min_case)
      case("SLL_DIRICHLET")
        bc_min = SLL_DIRICHLET
      case("SLL_NEUMANN")
        bc_min = SLL_NEUMANN
      case("SLL_NEUMANN_MODE_0")
        bc_min = SLL_NEUMANN_MODE_0
      case default
        print *,'#bad bc_min',bc_min
        print *,'#not implemented'
        print *,'#in initialize_guiding_center_2d_cartesian'
        stop
    end select
    select case(bc_max_case)
      case("SLL_DIRICHLET")
        bc_max = SLL_DIRICHLET
      case("SLL_NEUMANN")
        bc_max = SLL_NEUMANN
      case("SLL_NEUMANN_MODE_0")
        bc_max = SLL_NEUMANN_MODE_0
      case default
        print *,'#bad bc_max',bc_max
        print *,'#not implemented'
        print *,'#in initialize_guiding_center_2d_cartesian'
        stop
    end select

        sim%transformation => new_coordinate_transformation_2d_analytic( &
          "analytic_polar_transformation", &
          sim%mesh_2d, &
          polar_x1, &
          polar_x2, &
          polar_jac11, &
          polar_jac12, &
          polar_jac21, &
          polar_jac22, &
          params=(/0._f64,0._f64,0._f64,0._f64/))  


    select case(poisson_solver)    
      case ("SLL_POLAR_FFT")     
        sim%poisson =>new_poisson_2d_polar( &
          x1_min, &
          x1_max, &
          Nc_x1, &
          Nc_x2, &
          (/bc_min,bc_max/))
          !(/SLL_NEUMANN_MODE_0, SLL_DIRICHLET/))
          !(/SLL_DIRICHLET, SLL_DIRICHLET/))
      case ("SLL_ELLIPTIC_FINITE_ELEMENT_SOLVER")
        transformation => new_coordinate_transformation_2d_analytic( &
          "analytic_polar_transformation", &
          sim%mesh_2d, &
          polar_x1, &
          polar_x2, &
          polar_jac11, &
          polar_jac12, &
          polar_jac21, &
          polar_jac22, &
          params=(/0._f64,0._f64,0._f64,0._f64/))  

        SLL_ALLOCATE(b11(Nc_x1+1,Nc_x2+1),ierr)
        SLL_ALLOCATE(b12(Nc_x1+1,Nc_x2+1),ierr)
        SLL_ALLOCATE(b21(Nc_x1+1,Nc_x2+1),ierr)
        SLL_ALLOCATE(b22(Nc_x1+1,Nc_x2+1),ierr)
        SLL_ALLOCATE(b1(Nc_x1+1,Nc_x2+1),ierr)
        SLL_ALLOCATE(b2(Nc_x1+1,Nc_x2+1),ierr)
        SLL_ALLOCATE(c(Nc_x1+1,Nc_x2+1),ierr)
        
        b11 = 1._f64
        b22 = 1._f64
        b12 = 0._f64
        b21 = 0._f64
        c = 0._f64
        
        sim%poisson => new_poisson_2d_elliptic_solver( &
         transformation,&
         spline_degree_eta1, &
         spline_degree_eta2, &
         Nc_x1, &
         Nc_x2, &
         ES_GAUSS_LEGENDRE, &
         ES_GAUSS_LEGENDRE, &
         SLL_DIRICHLET, &
         SLL_DIRICHLET, &
         SLL_PERIODIC, &
         SLL_PERIODIC, &
         SLL_HERMITE, &
         SLL_PERIODIC, &
         x1_min, &
         x1_max, &
         x2_min, &
         x2_max, &
         b11, & 
         b12, & 
         b21, & 
         b22, & 
         b1, & 
         b2, & 
         c ) 
#ifdef MUDPACK
      case ("SLL_MUDPACK_CURVILINEAR")     
        transformation => new_coordinate_transformation_2d_analytic( &
          "analytic_polar_transformation", &
          sim%mesh_2d, &
          polar_x1, &
          polar_x2, &
          polar_jac11, &
          polar_jac12, &
          polar_jac21, &
          polar_jac22, &
          params=(/0._f64,0._f64,0._f64,0._f64/))  

        SLL_ALLOCATE(b11(Nc_x1+1,Nc_x2+1),ierr)
        SLL_ALLOCATE(b12(Nc_x1+1,Nc_x2+1),ierr)
        SLL_ALLOCATE(b21(Nc_x1+1,Nc_x2+1),ierr)
        SLL_ALLOCATE(b22(Nc_x1+1,Nc_x2+1),ierr)
        SLL_ALLOCATE(c(Nc_x1+1,Nc_x2+1),ierr)
        
        b11 = 1._f64
        b22 = 1._f64
        b12 = 0._f64
        b21 = 0._f64
        c = 0._f64


        sim%poisson => new_poisson_2d_mudpack_curvilinear_solver( &
         transformation, &
         x1_min,&
         x1_max,&
         Nc_x1,&
         x2_min,&
         x2_max,&
         Nc_x2,&
         SLL_DIRICHLET, &
         SLL_DIRICHLET, &
         SLL_PERIODIC, &
         SLL_PERIODIC, &
         SLL_HERMITE, &
         SLL_PERIODIC, &
         b11,&
         b12,&
         b21,&
         b22,&
         c)

#endif      
          
      case default
        print *,'#bad poisson_solver',poisson_solver
        print *,'#not implemented'
        print *,'#in initialize_guiding_center_2d_polar'
        stop
    end select


    
  
  end subroutine initialize_guiding_center_2d_polar
  


  subroutine init_fake(sim, filename)
    class(sll_simulation_2d_guiding_center_polar), intent(inout) :: sim
    character(len=*), intent(in)                                :: filename
  
    print *,'# Do not use the routine init_vp4d_fake'
    print *,'#use instead initialize_vlasov_par_poisson_seq_cart'
    print *,sim%dt
    print *,filename
    stop
  
  end subroutine init_fake
  
  subroutine run_gc2d_polar(sim)
    class(sll_simulation_2d_guiding_center_polar), intent(inout) :: sim
    sll_int32 :: Nc_x1
    sll_int32 :: Nc_x2
    sll_real64 :: delta_x1
    sll_real64 :: delta_x2
    sll_real64 :: x1_min
    sll_real64 :: x2_min    
    sll_real64 :: x1
    sll_real64 :: x2
    sll_int32 :: i1 
    sll_int32 :: i2
    sll_real64,dimension(:,:), pointer :: f
    sll_real64,dimension(:,:), pointer :: f_old
    sll_real64,dimension(:,:), pointer :: phi
    sll_real64,dimension(:,:), pointer :: A1 !advection fields
    sll_real64,dimension(:,:), pointer :: A2
    sll_int32 :: ierr
    sll_int32 :: nb_step
    sll_int32 :: step
    sll_real64 :: dt
    sll_int32 :: thdiag_id = 99 
    sll_int32 :: iplot
    sll_real64 :: time
    
    Nc_x1 = sim%mesh_2d%num_cells1
    Nc_x2 = sim%mesh_2d%num_cells2
    delta_x1 = sim%mesh_2d%delta_eta1
    delta_x2 = sim%mesh_2d%delta_eta2
    x1_min = sim%mesh_2d%eta1_min
    x2_min = sim%mesh_2d%eta2_min
    nb_step = sim%num_iterations
    dt = sim%dt
    
    
    !allocation
    SLL_ALLOCATE(f(Nc_x1+1,Nc_x2+1),ierr)
    SLL_ALLOCATE(f_old(Nc_x1+1,Nc_x2+1),ierr)
    SLL_ALLOCATE(phi(Nc_x1+1,Nc_x2+1),ierr)
    SLL_ALLOCATE(A1(Nc_x1+1,Nc_x2+1),ierr)
    SLL_ALLOCATE(A2(Nc_x1+1,Nc_x2+1),ierr)

    

    
    !initialisation of distribution function
    do i2=1,Nc_x2+1
      x2=x2_min+real(i2-1,f64)*delta_x2
      do i1=1,Nc_x1+1
        x1=x1_min+real(i1-1,f64)*delta_x1
        f(i1,i2) =  sim%init_func(x1,x2,sim%params)
      end do
    end do

        
    call sim%poisson%compute_phi_from_rho( phi, f )



    call compute_field_from_phi_2d_polar(phi,sim%mesh_2d,A1,A2,sim%phi_interp2d)
    
    
    call sll_ascii_file_create(sim%thdiag_filename, thdiag_id, ierr)
    
    
!    open(unit = thdiag_id, file='thdiag.dat',IOStat=IO_stat)
!    if( IO_stat /= 0 ) then
!       print *, '#run_gc2d_polar(sim) failed to open file thdiag.dat'
!       STOP
!    end if

    call sll_plot_polar_init( &
      sim%mesh_2d, &
      sim%transformation, &
      sim%mesh_name )    

    
    iplot = 0

    do step=1,nb_step+1
      f_old = f
      call sim%poisson%compute_phi_from_rho( phi, f_old )
      call compute_field_from_phi_2d_polar(phi,sim%mesh_2d,A1,A2,sim%phi_interp2d)      

      time = real(step-1,f64)*dt
      
      if(modulo(step-1,sim%freq_diag_time)==0)then
        call time_history_diagnostic_gc_polar( &
          thdiag_id, &    
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
        print*,"#step= ", step
        call sll_plot_f( &
          iplot, &
          f, &  
          Nc_x1+1, &
          Nc_x2+1,  &
          sim%f_name, &
          sim%mesh_name, &
          time )    
        call sll_plot_f( &
          iplot, &
          phi, &  
          Nc_x1+1, &
          Nc_x2+1,  &
          sim%phi_name, &
          sim%mesh_name, &
          time )    
!        call plot_f_polar(iplot,f,sim%mesh_2d)
        iplot = iplot+1  
      endif            
#endif
      
      select case (sim%time_loop_case)
        case (SLL_EULER)
          call sim%advect_2d%advect_2d(A1, A2, sim%dt, f_old, f)
        case (SLL_PREDICTOR_CORRECTOR)
          call sim%advect_2d%advect_2d(A1, A2, 0.5_f64*sim%dt, f_old, f)
          call sim%poisson%compute_phi_from_rho( phi, f )
          call compute_field_from_phi_2d_polar(phi,sim%mesh_2d,A1,A2,sim%phi_interp2d)      
          call sim%advect_2d%advect_2d(A1, A2, sim%dt, f_old, f)
        case default  
          print *,'#bad time_loop_case',sim%time_loop_case
          print *,'#not implemented'
          print *,'#in run_gc2d_polar'
          print *,'#available options are:'
          print *,'#SLL_EULER=',SLL_EULER
          print *,'#SLL_PREDICTOR_CORRECTOR=',SLL_PREDICTOR_CORRECTOR
          
      end select
         
    enddo
    
    close(thdiag_id)

    !print *,'#not implemented for the moment!'
  end subroutine run_gc2d_polar    
  
  
  subroutine compute_field_from_phi_2d_polar(phi,mesh_2d,A1,A2,interp2d)
    sll_real64, dimension(:,:), intent(in) :: phi
    sll_real64, dimension(:,:), intent(out) :: A1
    sll_real64, dimension(:,:), intent(out) :: A2
    type(sll_cartesian_mesh_2d), pointer :: mesh_2d
    class(sll_c_interpolator_2d), pointer   :: interp2d
    sll_int32 :: Nc_x1
    sll_int32 :: Nc_x2
    sll_real64 :: x1_min
    sll_real64 :: x2_min
    sll_real64 :: delta_x1
    sll_real64 :: delta_x2
    sll_real64 :: x1
    sll_real64 :: x2
    sll_int32 :: i1
    sll_int32 :: i2
    
    Nc_x1 = mesh_2d%num_cells1
    Nc_x2 = mesh_2d%num_cells2
    x1_min = mesh_2d%eta1_min
    x2_min = mesh_2d%eta2_min
    delta_x1 = mesh_2d%delta_eta1
    delta_x2 = mesh_2d%delta_eta2

    call interp2d%compute_interpolants(phi)

    do i2=1,Nc_x2+1
      x2=x2_min+real(i2-1,f64)*delta_x2
      do i1=1,Nc_x1+1
        x1=x1_min+real(i1-1,f64)*delta_x1
        A1(i1,i2)=interp2d%interpolate_from_interpolant_derivative_eta2(x1,x2)/x1
        A2(i1,i2)=-interp2d%interpolate_from_interpolant_derivative_eta1(x1,x2)/x1
      end do
    end do
    
    
    
  end subroutine compute_field_from_phi_2d_polar
  
  subroutine time_history_diagnostic_gc_polar( &
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
    sll_int32 :: Nc_x1
    sll_int32 :: Nc_x2
    sll_real64 ::x1_min
    sll_real64 ::x1_max
    sll_real64 :: delta_x1
    sll_real64 :: delta_x2
    sll_real64 :: x1
    sll_int32 :: ierr 
    type(sll_fft_plan), pointer         :: pfwd
    
    Nc_x1 = mesh_2d%num_cells1
    Nc_x2 = mesh_2d%num_cells2
    
    
    x1_min = mesh_2d%eta1_min
    x1_max = mesh_2d%eta1_max

    delta_x1 = mesh_2d%delta_eta1
    delta_x2 = mesh_2d%delta_eta2

    
    SLL_ALLOCATE(int_r(Nc_x2),ierr)
    SLL_ALLOCATE(data(Nc_x1+1),ierr)
    pfwd => fft_new_plan(Nc_x2,int_r,int_r,FFT_FORWARD,FFT_NORMALIZE)
 
    w     = 0.0_f64
    l1    = 0.0_f64
    l2    = 0.0_f64
    e     = 0.0_f64
    int_r = 0.0_f64
    
    
    do i2 = 1, Nc_x2
      do i1=1,Nc_x1+1
        x1 = x1_min+real(i1-1,f64)*delta_x1
        data(i1) = x1*f(i1,i2)
      enddo
      w = w + compute_integral_trapezoid_1d(data, Nc_x1+1, delta_x1)

      do i1=1,Nc_x1+1
        x1 = x1_min+real(i1-1,f64)*delta_x1
        data(i1) = x1*abs(f(i1,i2))
      enddo
      l1 = l1 + compute_integral_trapezoid_1d(data, Nc_x1+1, delta_x1)

      do i1=1,Nc_x1+1
        x1 = x1_min+real(i1-1,f64)*delta_x1
        data(i1) = x1*(f(i1,i2))**2
      enddo
      l2 = l2 + compute_integral_trapezoid_1d(data, Nc_x1+1, delta_x1)

      do i1=1,Nc_x1+1
        x1 = x1_min+real(i1-1,f64)*delta_x1
        data(i1) = x1*((x1*A2(i1,i2))**2+A1(i1,i2)**2)
      enddo
      e = e + compute_integral_trapezoid_1d(data, Nc_x1+1, delta_x1)

      do i1=1,Nc_x1+1
        x1 = x1_min+real(i1-1,f64)*delta_x1
        data(i1) = x1*phi(i1,i2)
      enddo
      int_r(i2) = compute_integral_trapezoid_1d(data, Nc_x1+1, delta_x1)      
    enddo     

    w = w*delta_x2
    l1 = l1*delta_x2
    l2 = sqrt(l2*delta_x2)
    e  = 0.5_f64*e*delta_x2
    call fft_apply_plan(pfwd,int_r,int_r)
    do i1=1,8
      !mode_slope(i1) = time_mode(i1)
      time_mode(i1) = abs(fft_get_mode(pfwd,int_r,i1-1))
      !mode_slope(i1) = &
      !  (log(0*time_mode(i1)+1.e-40_f64)-log(0*mode_slope(i1)+1.e-40_f64))/(dt+1.e-40_f64)
    enddo
    
    write(file_id,*) &
      dt*real(step,f64), &
      w, &
      l1, &
      l2, &
      e, &
      maxval(abs(phi(1:Nc_x1+1,1:Nc_x2+1))), &
      time_mode(1:8)!,mode_slope



    call fft_delete_plan(pfwd)
    
!    call fft_apply_plan(plan_sl%poisson%pfwd,int_r,int_r)
!
!
    
    
    
  end subroutine time_history_diagnostic_gc_polar


#ifndef NOHDF5
!*********************
!*********************

  !---------------------------------------------------
  ! Save the mesh structure
  !---------------------------------------------------
  subroutine plot_f_polar(iplot,f,mesh_2d)
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
    sll_real64, dimension(:,:), intent(in) :: f
    sll_real64 :: r
    sll_real64 :: theta
    sll_real64 :: rmin
    sll_real64 :: rmax
    sll_real64 :: dr
    sll_real64 :: dtheta
    
    
    nnodes_x1 = mesh_2d%num_cells1+1
    nnodes_x2 = mesh_2d%num_cells2+1
    rmin = mesh_2d%eta1_min
    rmax = mesh_2d%eta1_max
    dr = mesh_2d%delta_eta1
    dtheta = mesh_2d%delta_eta2
    
    !print *,'#maxf=',iplot,maxval(f),minval(f)
    

    
    if (iplot == 1) then

      SLL_ALLOCATE(x1(nnodes_x1,nnodes_x2), error)
      SLL_ALLOCATE(x2(nnodes_x1,nnodes_x2), error)
      do j = 1,nnodes_x2
        
        do i = 1,nnodes_x1
          r       = rmin+real(i-1,f64)*dr
          theta   = real(j-1,f64)*dtheta
          x1(i,j) = r*cos(theta)
          x2(i,j) = r*sin(theta)
        end do
      end do
      call sll_hdf5_file_create("polar_mesh-x1.h5",file_id,error)
      call sll_hdf5_write_array(file_id,x1,"/x1",error)
      call sll_hdf5_file_close(file_id, error)
      call sll_hdf5_file_create("polar_mesh-x2.h5",file_id,error)
      call sll_hdf5_write_array(file_id,x2,"/x2",error)
      call sll_hdf5_file_close(file_id, error)
      deallocate(x1)
      deallocate(x2)

    end if

    call int2string(iplot,cplot)
    call sll_xdmf_open("f"//cplot//".xmf","polar_mesh", &
      nnodes_x1,nnodes_x2,file_id,error)
    call sll_xdmf_write_array("f"//cplot,f,"values", &
      error,file_id,"Node")
    call sll_xdmf_close(file_id,error)
  end subroutine plot_f_polar

#endif

  subroutine sll_plot_polar_init( &
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
        
  end subroutine sll_plot_polar_init



end module sll_m_sim_bsl_gc_2d0v_polar
