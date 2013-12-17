module sll_simulation_2d_guiding_center_cartesian_module

!2d guiding center cartesian simulation
!related to simulation_2d_guiding_center_curvilinear.F90
!but here geometry and test are specifically cartesian

!see ../selalib/prototype/src/simulation/gcsim2d_cartesian_input.nml
!for example of use

!contact: Michel Mehrenberger (mehrenbe@math.unistra.fr)
!         Adnane Hamiaz (hamiaz@math.unistra.fr)


#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_field_2d.h"
#include "sll_utilities.h"
#include "sll_poisson_solvers.h"
  use sll_constants
  use sll_logical_meshes  
  use sll_module_advection_1d_periodic
  use sll_module_advection_2d_BSL
  use sll_module_characteristics_2d_explicit_euler
  use sll_module_characteristics_2d_verlet
  use sll_reduction_module
  use sll_simulation_base
  use sll_cubic_spline_interpolator_2d
  use sll_cubic_spline_interpolator_1d
  use sll_coordinate_transformation_2d_base_module
  use sll_module_coordinate_transformations_2d
  use sll_common_coordinate_transformations
  use sll_common_array_initializers_module
  !use sll_mudpack_curvilinear
  use sll_module_poisson_2d_mudpack_solver
  use sll_module_poisson_2d_mudpack_curvilinear_solver_old
  use sll_module_poisson_2d_elliptic_solver
  use sll_module_scalar_field_2d_base
  use sll_module_scalar_field_2d_alternative
  use sll_timer
  use sll_fft
  use sll_module_poisson_2d_periodic_solver

!#include "sll_working_precision.h"
!#include "sll_assert.h"
!#include "sll_memory.h"
!#include "sll_field_2d.h"
!#include "sll_utilities.h"
!#include "sll_poisson_solvers.h"
!  use sll_constants
!  use sll_logical_meshes  
!  use sll_module_advection_1d_periodic
!  use sll_module_advection_2d_BSL
!  use sll_module_characteristics_2d_explicit_euler
!  use sll_module_characteristics_2d_verlet
!  use sll_reduction_module
!  use sll_simulation_base
!  use sll_cubic_spline_interpolator_2d
!  use sll_cubic_spline_interpolator_1d
!  use sll_coordinate_transformation_2d_base_module
!  use sll_module_coordinate_transformations_2d
!  use sll_common_coordinate_transformations
!  use sll_common_array_initializers_module
!  use sll_mudpack_cartesian
!  use sll_module_poisson_2d_mudpack_solver
!  use sll_module_poisson_2d_periodic_solver
  !use sll_parallel_array_initializer_module

  implicit none

  
  sll_int32, parameter :: SLL_EULER = 0 
  sll_int32, parameter :: SLL_PREDICTOR_CORRECTOR = 1 
  sll_int32, parameter :: SLL_PHI_FROM_RHO = 0
  sll_int32, parameter :: SLL_E_FROM_RHO = 1


  type, extends(sll_simulation_base_class) :: &
    sll_simulation_2d_guiding_center_cartesian

   !geometry
   type(sll_logical_mesh_2d), pointer :: mesh_2d


   !initial function
   procedure(sll_scalar_initializer_2d), nopass, pointer :: init_func
   sll_real64, dimension(:), pointer :: params
      
   !advector
   class(sll_advection_2d_base), pointer    :: advect_2d
   
   !interpolator for derivatives
   class(sll_interpolator_2d_base), pointer   :: phi_interp2d

   
   !poisson solver
   class(sll_poisson_2d_base), pointer   :: poisson
   sll_int32 :: poisson_case

   !time_iterations
   sll_real64 :: dt
   sll_int32  :: num_iterations
   sll_int32  :: freq_diag
   sll_int32  :: freq_diag_time

   !time_loop
   sll_int32 :: time_loop_case
   
       
  contains
    procedure, pass(sim) :: run => run_gc2d_cartesian
    procedure, pass(sim) :: init_from_file => init_fake
     
  end type sll_simulation_2d_guiding_center_cartesian


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

  function new_guiding_center_2d_cartesian(filename) result(sim)
    type(sll_simulation_2d_guiding_center_cartesian), pointer :: sim    
    character(len=*), intent(in), optional :: filename
    sll_int32 :: ierr
    
    SLL_ALLOCATE(sim,ierr)
    
    call initialize_guiding_center_2d_cartesian(sim,filename)
    
  
  
  end function new_guiding_center_2d_cartesian
  
  subroutine initialize_guiding_center_2d_cartesian(sim, filename)
    class(sll_simulation_2d_guiding_center_cartesian), intent(inout) :: sim
    character(len=*), intent(in), optional :: filename
    sll_int32             :: IO_stat
    sll_int32, parameter  :: input_file = 99
    
    !geometry
    character(len=256) :: mesh_case_x1
    sll_int32 :: num_cells_x1
    sll_real64 :: x1_min
    sll_int32 :: nbox_x1
    character(len=256) :: mesh_case_x2
    sll_int32 :: num_cells_x2
    sll_real64 :: x2_min
    sll_int32 :: nbox_x2
    
    !initial_function
    character(len=256) :: initial_function_case
    sll_real64 :: kmode_x1
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

 
    !poisson
    character(len=256) :: poisson_case
    character(len=256) :: poisson_solver
    character(len=256) :: mudpack_method    
    sll_int32 :: spline_degree_eta1
    sll_int32 :: spline_degree_eta2

    !local variables
    sll_int32 :: Nc_x1
    sll_int32 :: Nc_x2
    sll_real64 :: x1_max
    sll_real64 :: x2_max     
    class(sll_interpolator_2d_base), pointer :: f_interp2d
    class(sll_interpolator_2d_base), pointer :: phi_interp2d
    class(sll_characteristics_2d_base), pointer :: charac2d
    class(sll_interpolator_2d_base), pointer   :: A1_interp2d
    class(sll_interpolator_2d_base), pointer   :: A2_interp2d
    class(sll_interpolator_1d_base), pointer   :: A1_interp1d_x1
    class(sll_interpolator_1d_base), pointer   :: A2_interp1d_x1
    sll_real64, dimension(:,:), pointer :: b11
    sll_real64, dimension(:,:), pointer :: b12
    sll_real64, dimension(:,:), pointer :: b21
    sll_real64, dimension(:,:), pointer :: b22
    sll_real64, dimension(:,:), pointer :: c
    class(sll_coordinate_transformation_2d_base), pointer :: transformation
    sll_real64, dimension(:,:), allocatable :: cxx_2d
    sll_real64, dimension(:,:), allocatable :: cxy_2d
    sll_real64, dimension(:,:), allocatable :: cyy_2d
    sll_real64, dimension(:,:), allocatable :: cx_2d
    sll_real64, dimension(:,:), allocatable :: cy_2d
    sll_real64, dimension(:,:), allocatable :: ce_2d
    sll_int32 :: ierr

    namelist /geometry/ &
      mesh_case_x1, &
      num_cells_x1, &
      x1_min, &
      nbox_x1, &
      mesh_case_x2, &
      num_cells_x2, &
      x2_min, &
      nbox_x2

    namelist /initial_function/ &
      initial_function_case, &
      kmode_x1, &
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
      A_interp_case

    namelist /poisson/ &
      poisson_case, &
      poisson_solver, &
      mudpack_method, &    
      spline_degree_eta1, &
      spline_degree_eta2    


    !set default parameters
    
    !geometry
    mesh_case_x1="SLL_LANDAU_MESH"
    num_cells_x1 = 32
    x1_min = 0.0_f64
    nbox_x1 = 1
    mesh_case_x2="SLL_LANDAU_MESH"
    num_cells_x2 = 32
    x2_min = 0.0_f64
    nbox_x2 = 1
    
    !initial function
    initial_function_case="SLL_KHP1"
    kmode_x1 = 0.5_f64
    kmode_x2 = 1._f64
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
    
    !poisson
    poisson_case = "SLL_PHI_FROM_RHO"
    !poisson_case = "SLL_E_FROM_RHO"    
    poisson_solver = "SLL_POISSON_FFT"  
    !poisson_solver = "SLL_ELLIPTIC_FINITE_ELEMENT_SOLVER" !use with "SLL_PHI_FROM_RHO"
    !poisson_solver = "SLL_MUDPACK"   !use with "SLL_PHI_FROM_RHO"
    !poisson_solver = "SLL_MUDPACK_CURVILINEAR"   !use with "SLL_PHI_FROM_RHO"    
    !mudpack_method = "SLL_SEPARABLE"
    !mudpack_method = "SLL_NON_SEPARABLE_WITHOUT_CROSS_TERMS"
    mudpack_method = "SLL_NON_SEPARABLE_WITH_CROSS_TERMS"    
    spline_degree_eta1 = 3
    spline_degree_eta2 = 3    


    if(present(filename))then
      open(unit = input_file, file=trim(filename)//'.nml',IOStat=IO_stat)
        if( IO_stat /= 0 ) then
          print *, '#initialize_guiding_center_2d_cartesian() failed to open file ', &
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




    x1_max = 2._f64*sll_pi/kmode_x1
    x2_max = 2._f64*sll_pi/kmode_x2
    Nc_x1 = num_cells_x1
    Nc_x2 = num_cells_x2
     
    
!    Nc_x1 = 32
!    Nc_x2 = 32
!    k_mode = 0.5_f64
!    eps = 0.015_f64
!    x1_min = 0._f64
!    x1_max = 2._f64*sll_pi/k_mode
!    x2_min = 0._f64
!    x2_max = 2._f64*sll_pi
!    nb_step = 600
!    
!    dt = 0.1_f64
!    visu_step = 100
!    f_interp2d_case = "SLL_CUBIC_SPLINES"
!    phi_interp2d_case = "SLL_CUBIC_SPLINES"
!    A_interp_case = "SLL_CUBIC_SPLINES"
!    charac2d_case = "SLL_VERLET"
!    !charac2d_case = "SLL_EULER"
!    advect2d_case = "SLL_BSL"
!    initial_function_case = "SLL_KHP1" 
!    !time_loop_case = "SLL_EULER"
!    time_loop_case = "SLL_PREDICTOR_CORRECTOR" 
!    !poisson_solver = "SLL_MUDPACK"   !use with "SLL_PHI_FROM_RHO"
!    poisson_solver = "SLL_POISSON_FFT"  
!    !poisson_solver = "SLL_ELLIPTIC_FINITE_ELEMENT_SOLVER" !use with "SLL_PHI_FROM_RHO"
!    !poisson_solver = "SLL_MUDPACK_CURVILINEAR"   !use with "SLL_PHI_FROM_RHO"
!    
!    poisson_case = "SLL_PHI_FROM_RHO"
!    !poisson_case = "SLL_E_FROM_RHO"
!    !mudpack_method = "SLL_SEPARABLE"
!    !mudpack_method = "SLL_NON_SEPARABLE_WITHOUT_CROSS_TERMS"
!    mudpack_method = "SLL_NON_SEPARABLE_WITH_CROSS_TERMS"
!    
!    spline_degree_eta1 = 3
!    spline_degree_eta2 = 3
    
    
    
    sim%dt = dt
    sim%num_iterations = number_iterations
    sim%freq_diag = freq_diag
    sim%freq_diag_time = freq_diag_time

    sim%mesh_2d => new_logical_mesh_2d( &
      num_cells_x1, &
      num_cells_x2, &
      eta1_min = x1_min, &
      eta1_max = x1_max, &
      eta2_min = x2_min, &
      eta2_max = x2_max)      
      
      
      
    select case (f_interp2d_case)
      case ("SLL_CUBIC_SPLINES")
        f_interp2d => new_cubic_spline_2d_interpolator( &
          Nc_x1+1, &
          Nc_x2+1, &
          x1_min, &
          x1_max, &
          x2_min, &
          x2_max, &
          SLL_PERIODIC, &
          SLL_PERIODIC)
      case default
        print *,'#bad f_interp2d_case',f_interp2d_case
        print *,'#not implemented'
        print *,'#in initialize_guiding_center_2d_polar'
        stop
    end select




    select case (A_interp_case)
      case ("SLL_CUBIC_SPLINES")
        A1_interp2d => new_cubic_spline_2d_interpolator( &
          Nc_x1+1, &
          Nc_x2+1, &
          x1_min, &
          x1_max, &
          x2_min, &
          x2_max, &
          SLL_PERIODIC, &
          SLL_PERIODIC)
        A2_interp2d => new_cubic_spline_2d_interpolator( &
          Nc_x1+1, &
          Nc_x2+1, &
          x1_min, &
          x1_max, &
          x2_min, &
          x2_max, &
          SLL_PERIODIC, &
          SLL_PERIODIC)  
        A1_interp1d_x1 => new_cubic_spline_1d_interpolator( &
          Nc_x1+1, &
          x1_min, &
          x1_max, &
          SLL_PERIODIC)
        A2_interp1d_x1 => new_cubic_spline_1d_interpolator( &
          Nc_x1+1, &
          x1_min, &
          x1_max, &
          SLL_PERIODIC)
      case default
        print *,'#bad A_interp_case',A_interp_case
        print *,'#not implemented'
        print *,'#in initialize_guiding_center_2d_polar'
        stop
    end select

    select case (phi_interp2d_case)
      case ("SLL_CUBIC_SPLINES")
        phi_interp2d => new_cubic_spline_2d_interpolator( &
          Nc_x1+1, &
          Nc_x2+1, &
          x1_min, &
          x1_max, &
          x2_min, &
          x2_max, &
          SLL_PERIODIC, &
          SLL_PERIODIC)         
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
          bc_type_1=SLL_PERIODIC, &!&SLL_SET_TO_LIMIT, &
          bc_type_2=SLL_PERIODIC)    
      case ("SLL_VERLET")      
        charac2d => new_verlet_2d_charac(&
          Nc_x1+1, &
          Nc_x2+1, &
          A1_interp2d, &
          A2_interp2d, &
          A1_interp1d_x1, &
          A2_interp1d_x1, &
          bc_type_1=SLL_PERIODIC, &!&SLL_SET_TO_LIMIT, &
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
      case ("SLL_KHP1")
        sim%init_func => sll_KHP1_2d
        SLL_ALLOCATE(sim%params(2),ierr)
        sim%params(1) = eps
        sim%params(2) = kmode_x1
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
    

    
    !poisson solver
    
    select case(poisson_case)    
      case ("SLL_PHI_FROM_RHO")     
        sim%poisson_case = SLL_PHI_FROM_RHO
      case ("SLL_E_FROM_RHO")
        sim%poisson_case = SLL_E_FROM_RHO     
      case default
        print *,'#bad poisson_case',poisson_case
        print *,'#not implemented'
        print *,'#in initialize_guiding_center_2d_cartesian'
        stop
    end select
    
    select case(poisson_solver)    
      case ("SLL_MUDPACK")     
        !stop  
        
        select case(mudpack_method)
          case ("SLL_SEPARABLE")
            sim%poisson => new_poisson_2d_mudpack_solver( &
              x1_min,&
              x1_max,&
              Nc_x1,&
              x2_min,&
              x2_max,&
              Nc_x2,&
              SLL_PERIODIC,& 
              SLL_PERIODIC,& 
              SLL_PERIODIC,& 
              SLL_PERIODIC,&
              SLL_SEPARABLE, &
              cxx = 1._f64, &
              cyy = 1._f64)
          case ("SLL_NON_SEPARABLE_WITHOUT_CROSS_TERMS")   
            SLL_ALLOCATE(cxx_2d(Nc_x1+1,Nc_x2+1),ierr)
            SLL_ALLOCATE(cyy_2d(Nc_x1+1,Nc_x2+1),ierr)
            SLL_ALLOCATE(cx_2d(Nc_x1+1,Nc_x2+1),ierr)
            SLL_ALLOCATE(cy_2d(Nc_x1+1,Nc_x2+1),ierr)
            SLL_ALLOCATE(ce_2d(Nc_x1+1,Nc_x2+1),ierr)
            
            cxx_2d = 1._f64
            cyy_2d = 1._f64
            cx_2d = 0._f64
            cy_2d = 0._f64
            ce_2d = 0._f64
             
            sim%poisson => new_poisson_2d_mudpack_solver( &
              x1_min,&
              x1_max,&
              Nc_x1,&
              x2_min,&
              x2_max,&
              Nc_x2,&
              SLL_PERIODIC,& 
              SLL_PERIODIC,& 
              SLL_PERIODIC,& 
              SLL_PERIODIC,&
              SLL_NON_SEPARABLE_WITHOUT_CROSS_TERMS, &
              cxx_2d = cxx_2d, &
              cyy_2d = cyy_2d, &
              cx_2d = cx_2d, &
              cy_2d = cy_2d, &
              ce_2d = ce_2d)
              
          case ("SLL_NON_SEPARABLE_WITH_CROSS_TERMS")   
            SLL_ALLOCATE(cxx_2d(Nc_x1+1,Nc_x2+1),ierr)
            SLL_ALLOCATE(cxy_2d(Nc_x1+1,Nc_x2+1),ierr)
            SLL_ALLOCATE(cyy_2d(Nc_x1+1,Nc_x2+1),ierr)
            SLL_ALLOCATE(cx_2d(Nc_x1+1,Nc_x2+1),ierr)
            SLL_ALLOCATE(cy_2d(Nc_x1+1,Nc_x2+1),ierr)
            SLL_ALLOCATE(ce_2d(Nc_x1+1,Nc_x2+1),ierr)
            
            cxx_2d = 1._f64
            cxy_2d = 0._f64
            cyy_2d = 1._f64
            cx_2d = 0._f64
            cy_2d = 0._f64
            ce_2d = 0._f64
             
            sim%poisson => new_poisson_2d_mudpack_solver( &
              x1_min,&
              x1_max,&
              Nc_x1,&
              x2_min,&
              x2_max,&
              Nc_x2,&
              SLL_PERIODIC,& 
              SLL_PERIODIC,& 
              SLL_PERIODIC,& 
              SLL_PERIODIC,&
              SLL_NON_SEPARABLE_WITH_CROSS_TERMS, &
              cxx_2d = cxx_2d, &
              cxy_2d = cxy_2d, &
              cyy_2d = cyy_2d, &
              cx_2d = cx_2d, &
              cy_2d = cy_2d, &
              ce_2d = ce_2d)
              
          case default
            print *,'#bad mudpack_method',mudpack_method
            print *,'#in initialize_guiding_center_2d_cartesian'
            stop
        end select    
      case ("SLL_MUDPACK_CURVILINEAR")     
!        transformation => new_coordinate_transformation_2d_analytic( &
!          "analytic_identity_transformation", &
!          sim%mesh_2d, &
!          identity_x1, &
!          identity_x2, &
!          identity_jac11, &
!          identity_jac12, &
!          identity_jac21, &
!          identity_jac22, &
!          params=(/0._f64,0._f64,0._f64,0._f64/))  
        transformation => new_coordinate_transformation_2d_analytic( &
          "analytic_collela_transformation", &
          sim%mesh_2d, &
          sinprod_x1, &
          sinprod_x2, &
          sinprod_jac11, &
          sinprod_jac12, &
          sinprod_jac21, &
          sinprod_jac22, &
          params=(/ 0._f64, 0._f64, x1_max-x1_min, x2_max-x2_min/)  )  
        !  In collela  mesh params_mesh =( alpha1, alpha2, L1, L2 ) such that :
        !  x1= eta1 + alpha1*sin(2*pi*eta1/L1)*sin(2*pi*eta2/L2)



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
         SLL_PERIODIC, &
         SLL_PERIODIC, &
         SLL_PERIODIC, &
         SLL_PERIODIC, &
         b11,&
         b12,&
         b21,&
         b22,&
         c)

      case ("SLL_POISSON_FFT")     
        !stop  
        sim%poisson => new_poisson_2d_periodic_solver( &
          x1_min,&
          x1_max,&
          Nc_x1,&
          x2_min,&
          x2_max,&
          Nc_x2)
      case ("SLL_ELLIPTIC_FINITE_ELEMENT_SOLVER")
        transformation => new_coordinate_transformation_2d_analytic( &
          "analytic_identity_transformation", &
          sim%mesh_2d, &
          identity_x1, &
          identity_x2, &
          identity_jac11, &
          identity_jac12, &
          identity_jac21, &
          identity_jac22, &
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
        
        sim%poisson => new_poisson_2d_elliptic_solver( &
         transformation,&
         spline_degree_eta1, &
         spline_degree_eta2, &
         Nc_x1, &
         Nc_x2, &
         ES_GAUSS_LEGENDRE, &
         ES_GAUSS_LEGENDRE, &
         SLL_PERIODIC, &
         SLL_PERIODIC, &
         SLL_PERIODIC, &
         SLL_PERIODIC, &
         x1_min, &
         x1_max, &
         x2_min, &
         x2_max, &
         b11, & 
         b12, & 
         b21, & 
         b22, & 
         c ) 

         
      case default
        print *,'#bad poisson_solver',poisson_solver
        print *,'#not implemented'
        print *,'#in initialize_guiding_center_2d_cartesian'
        stop
    end select





   
  end subroutine initialize_guiding_center_2d_cartesian
  


  subroutine init_fake(sim, filename)
    class(sll_simulation_2d_guiding_center_cartesian), intent(inout) :: sim
    character(len=*), intent(in)                                :: filename
  
    print *,'# Do not use the routine init_vp4d_fake'
    print *,'#use instead initialize_vlasov_par_poisson_seq_cart'
    stop
  
  end subroutine init_fake
  
  subroutine run_gc2d_cartesian(sim)
    class(sll_simulation_2d_guiding_center_cartesian), intent(inout) :: sim
    sll_int32 :: Nc_x1
    sll_int32 :: Nc_x2
    sll_real64 :: delta_x1
    sll_real64 :: delta_x2
    sll_real64 :: x1_min,x1_max
    sll_real64 :: x2_min,x2_max   
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
    sll_int32 :: diag_id = 88 
    sll_int32             :: IO_stat
    sll_int32 :: iplot
    
    Nc_x1 = sim%mesh_2d%num_cells1
    Nc_x2 = sim%mesh_2d%num_cells2
    delta_x1 = sim%mesh_2d%delta_eta1
    delta_x2 = sim%mesh_2d%delta_eta2
    x1_min = sim%mesh_2d%eta1_min
    x2_min = sim%mesh_2d%eta2_min
    x1_max = sim%mesh_2d%eta1_max
    x2_max = sim%mesh_2d%eta2_max
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
        

    select case(sim%poisson_case) 
      case (SLL_PHI_FROM_RHO)        
        call sim%poisson%compute_phi_from_rho( phi, f )            
        call compute_field_from_phi_2d_cartesian(phi,sim%mesh_2d,A1,A2,sim%phi_interp2d)      
      case (SLL_E_FROM_RHO)
        call sim%poisson%compute_E_from_rho( A2, A1, f )        
        A1 = -A1
     end select

     call sll_gnuplot_write(A1(1,:),'A1_init',ierr)
     call sll_gnuplot_write(A2(1,:),'A2_init',ierr)
    
    
    open(unit = diag_id, file='thdiag.dat',IOStat=IO_stat)
    if( IO_stat /= 0 ) then
       print *, '#run_gc2d_cartesian (sim) failed to open file thdiag.dat'
       STOP
    end if
    
    iplot = 0

    do step=1,nb_step+1
      f_old = f

      select case(sim%poisson_case) 
        case (SLL_PHI_FROM_RHO)
          call sim%poisson%compute_phi_from_rho( phi, f_old )    
          call compute_field_from_phi_2d_cartesian(phi,sim%mesh_2d,A1,A2,sim%phi_interp2d)      
        case (SLL_E_FROM_RHO)
          call sim%poisson%compute_E_from_rho( A2, A1, f_old )              
          A1 = -A1
      end select
      !call poisson_solve_cartesian(sim%poisson,f_old,phi)
      !call solve_mudpack_cartesian(sim%poisson, phi, -f_old)
      
      
      if(modulo(step-1,sim%freq_diag_time)==0)then
        if(sim%poisson_case==SLL_E_FROM_RHO) then
          call sim%poisson%compute_phi_from_rho( phi, f_old )
        endif
        call time_history_diagnostic_gc_cartesian( &
          diag_id, &    
          step-1, &
          dt, &
          sim%mesh_2d, &
          f, &
          phi, &
          A1, &
          A2)
      endif            
     
      if(modulo(step-1,sim%freq_diag)==0)then
        print*,"#step= ", step
        call plot_f_cartesian(iplot,f,sim%mesh_2d)
        iplot = iplot+1  
      endif            
      
      select case (sim%time_loop_case)
        case (SLL_EULER)
          call sim%advect_2d%advect_2d(A1, A2, sim%dt, f_old, f)
        case (SLL_PREDICTOR_CORRECTOR)
          call sim%advect_2d%advect_2d(A1, A2, 0.5_f64*sim%dt, f_old, f)

          select case(sim%poisson_case) 
            case (SLL_PHI_FROM_RHO)
              call sim%poisson%compute_phi_from_rho( phi, f )    
              call compute_field_from_phi_2d_cartesian( &
                phi, &
                sim%mesh_2d, &
                A1, &
                A2, &
                sim%phi_interp2d)      
            case (SLL_E_FROM_RHO)
              call sim%poisson%compute_E_from_rho( A2, A1, f )
              A1 = -A1
          end select
          call sim%advect_2d%advect_2d(A1, A2, sim%dt, f_old, f)
        case default  
          print *,'#bad time_loop_case',sim%time_loop_case
          print *,'#not implemented'
          print *,'#in run_gc2d_cartesian'
          print *,'#available options are:'
          print *,'#SLL_EULER=',SLL_EULER
          print *,'#SLL_PREDICTOR_CORRECTOR=',SLL_PREDICTOR_CORRECTOR
          
      end select
         
    enddo
    
    close(diag_id)

    print *,'#run_gc2d_cartesian PASSED'
  end subroutine run_gc2d_cartesian  
  
  
  subroutine compute_field_from_phi_2d_cartesian(phi,mesh_2d,A1,A2,interp2d)
    sll_real64, dimension(:,:), intent(in) :: phi
    sll_real64, dimension(:,:), intent(out) :: A1
    sll_real64, dimension(:,:), intent(out) :: A2
    type(sll_logical_mesh_2d), pointer :: mesh_2d
    class(sll_interpolator_2d_base), pointer   :: interp2d
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
        A1(i1,i2)=interp2d%interpolate_derivative_eta2(x1,x2)
        A2(i1,i2)=-interp2d%interpolate_derivative_eta1(x1,x2)
      end do
    end do
    
    
    
  end subroutine compute_field_from_phi_2d_cartesian
  
  subroutine time_history_diagnostic_gc_cartesian( &
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

    SLL_ALLOCATE(data(Nc_x1+1),ierr)
 
    linf  = 0.0_f64
    l1    = 0.0_f64
    l2    = 0.0_f64
    mass  = 0.0_f64
     e     = 0.0_f64
    
    do i2 = 1, Nc_x2+1
      do i1=1,Nc_x1+1
        data(i1) = f(i1,i2)
      enddo
      mass = mass + compute_integral_trapezoid_1d(data, Nc_x1+1, delta_x1)

      do i1=1,Nc_x1+1
        data(i1) = abs(f(i1,i2))
      enddo
      l1 = l1 + compute_integral_trapezoid_1d(data, Nc_x1+1, delta_x1)

      do i1=1,Nc_x1+1
        data(i1) = (f(i1,i2))**2
      enddo
      l2 = l2 + compute_integral_trapezoid_1d(data, Nc_x1+1, delta_x1)

      do i1=1,Nc_x1+1
        data(i1) = A2(i1,i2)**2+A1(i1,i2)**2
      enddo
      e = e + compute_integral_trapezoid_1d(data, Nc_x1+1, delta_x1)

      do i1=1,Nc_x1+1
       linf = max(linf,abs(f(i1,i2)))
      enddo
         
    enddo     

    mass = mass*delta_x2
    l1 = l1*delta_x2
    l2 = sqrt(l2*delta_x2)
    e  = e*delta_x2


    
    write(file_id,*) &
      dt*real(step,f64), &
      linf, &
      l1, &
      l2, &
      mass, &
      e, &
      maxval(abs(A1(1:Nc_x1+1,1:Nc_x2+1)**2+A2(1:Nc_x1+1,1:Nc_x2+1)**2)), &
      maxval(abs(phi(1:Nc_x1+1,1:Nc_x2+1)))


    
    
  end subroutine time_history_diagnostic_gc_cartesian


#ifndef NOHDF5
!*********************
!*********************

  !---------------------------------------------------
  ! Save the mesh structure
  !---------------------------------------------------
  subroutine plot_f_cartesian(iplot,f,mesh_2d)
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
    sll_real64, dimension(:,:), intent(in) :: f
    sll_real64 :: r
    sll_real64 :: theta
    sll_real64 ::  x1_min, x2_min
    sll_real64 ::  x1_max, x2_max  
    sll_real64 :: dx1
    sll_real64 :: dx2
    
    
    nnodes_x1 = mesh_2d%num_cells1+1
    nnodes_x2 = mesh_2d%num_cells2+1
    x1_min = mesh_2d%eta1_min
    x1_max = mesh_2d%eta1_max
    x2_min = mesh_2d%eta2_min
    x2_max = mesh_2d%eta2_max
    dx1 = mesh_2d%delta_eta1
    dx2 = mesh_2d%delta_eta2
    
    !print *,'#maxf=',iplot,maxval(f),minval(f)
    

    
    if (iplot == 1) then

      SLL_ALLOCATE(x1(nnodes_x1,nnodes_x2), error)
      SLL_ALLOCATE(x2(nnodes_x1,nnodes_x2), error)
      do j = 1,nnodes_x2
        do i = 1,nnodes_x1
          x1(i,j) = x1_min+real(i-1,f32)*dx1
          x2(i,j) = x2_min+real(j-1,f32)*dx2
        end do
      end do
      call sll_hdf5_file_create("cartesian_mesh-x1.h5",file_id,error)
      call sll_hdf5_write_array(file_id,x1,"/x1",error)
      call sll_hdf5_file_close(file_id, error)
      call sll_hdf5_file_create("cartesian_mesh-x2.h5",file_id,error)
      call sll_hdf5_write_array(file_id,x2,"/x2",error)
      call sll_hdf5_file_close(file_id, error)
      deallocate(x1)
      deallocate(x2)

    end if

    call int2string(iplot,cplot)
    call sll_xdmf_open("f"//cplot//".xmf","cartesian_mesh", &
      nnodes_x1,nnodes_x2,file_id,error)
    call sll_xdmf_write_array("f"//cplot,f,"values", &
      error,file_id,"Node")
    call sll_xdmf_close(file_id,error)
  end subroutine plot_f_cartesian

#endif


end module sll_simulation_2d_guiding_center_cartesian_module
