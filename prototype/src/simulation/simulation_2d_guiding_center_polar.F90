module sll_simulation_2d_guiding_center_polar_module

!the aim is to translate CG_polar in simulation class
!related to
!simulation_guiding_center_2D_generalized_coords.F90
!but here geometry and test is specifically polar

#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_field_2d.h"
#include "sll_utilities.h"
#include "sll_poisson_solvers.h"
  use sll_logical_meshes  
  use sll_module_advection_1d_periodic
  use sll_module_advection_2d_BSL
  use sll_module_characteristics_2d_explicit_euler
  use sll_module_characteristics_2d_verlet
  !use sll_poisson_2d_periodic  
  use sll_fft
  use sll_simulation_base
  use sll_cubic_spline_interpolator_2d
  use sll_cubic_spline_interpolator_1d
  use sll_coordinate_transformation_2d_base_module
  use sll_module_coordinate_transformations_2d
  use sll_common_coordinate_transformations
  use sll_common_array_initializers_module
  !use sll_parallel_array_initializer_module

  implicit none


  type, extends(sll_simulation_base_class) :: &
    sll_simulation_2d_guiding_center_polar

   !geometry
   type(sll_logical_mesh_2d), pointer :: mesh_2d


   !initial function
   procedure(sll_scalar_initializer_2d), nopass, pointer :: init_func
   sll_real64, dimension(:), pointer :: params
      
   !advector (should replace interpolator)
   class(sll_advection_2d_base), pointer    :: advect_2d
   
   !poisson solver
   !class(sll_poisson_2d_base), pointer   :: poisson
   !type(poisson_2d_periodic), pointer   :: poisson
   type(sll_plan_poisson_polar), pointer :: poisson 
   
   !time_iterations
   sll_real64 :: dt
   sll_int32  :: num_iterations
   sll_int32  :: freq_diag
   sll_int32  :: freq_diag_time


       
  contains
    procedure, pass(sim) :: run => run_gc2d_polar
    procedure, pass(sim) :: init_from_file => init_fake
     
  end type sll_simulation_2d_guiding_center_polar


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

  function new_guiding_center_2d_polar() result(sim)
    type(sll_simulation_2d_guiding_center_polar), pointer :: sim    
    sll_int32 :: ierr
    
    SLL_ALLOCATE(sim,ierr)
    
    call initialize_guiding_center_2d_polar(sim)
    
  
  
  end function new_guiding_center_2d_polar
  
  subroutine initialize_guiding_center_2d_polar(sim)
    class(sll_simulation_2d_guiding_center_polar), intent(inout) :: sim
    sll_int32 :: Nc_x1
    sll_int32 :: Nc_x2
    sll_real64 :: x1_min
    sll_real64 :: x1_max
    sll_real64 :: x2_min
    sll_real64 :: x2_max
    sll_real64 :: r_minus
    sll_real64 :: r_plus
    sll_real64 :: eps
    sll_real64 :: k_mode
    sll_int32  :: nb_step
    sll_real64 :: dt
    sll_int32 :: visu_step
    class(sll_interpolator_2d_base), pointer :: interp2d
    class(sll_characteristics_2d_base), pointer :: charac2d
    class(sll_interpolator_2d_base), pointer   :: A1_interp2d
    class(sll_interpolator_2d_base), pointer   :: A2_interp2d
    class(sll_interpolator_1d_base), pointer   :: A1_interp1d_x1
    class(sll_interpolator_1d_base), pointer   :: A2_interp1d_x1
    character(len=256)      :: advect2d_case 
    character(len=256)      :: charac2d_case
    character(len=256)      :: interp2d_case 
    character(len=256)      :: interp1d_x1_case 
    character(len=256)      :: initial_function_case 
    sll_int32 :: ierr
    !character(len=256)      :: interp1d_x2_case 
    
    !here we do all the initialization
    !in future, we will use namelist file

    x1_min = 1._f64
    x1_max = 10._f64
    Nc_x1 = 128
    Nc_x2 = 128
    r_minus = 5._f64
    r_plus = 6._f64
    k_mode = 3
    nb_step = 100
    dt = 0.1_f64
    visu_step = 20
    interp1d_x1_case = "SLL_CUBIC_SPLINES"
    interp2d_case = "SLL_CUBIC_SPLINES"
    charac2d_case = "SLL_VERLET"
    !charac2d_case = "SLL_EULER"
    advect2d_case = "SLL_BSL"
    initial_function_case = "SLL_DIOCOTRON" 
    
    
    
    sim%dt = dt
    sim%num_iterations = nb_step
    sim%freq_diag = visu_step
    sim%freq_diag_time = 1

    sim%mesh_2d => new_logical_mesh_2d( &
      Nc_x1, &
      Nc_x2, &
      eta1_min = x1_min, &
      eta1_max = x1_max, &
      eta2_min = x2_min, &
      eta2_max = x2_max)
      
    select case (interp2d_case)
      case ("SLL_CUBIC_SPLINES")
        interp2d => new_cubic_spline_2d_interpolator( &
          Nc_x1+1, &
          Nc_x2+1, &
          x1_min, &
          x1_max, &
          x2_min, &
          x2_max, &
          SLL_HERMITE, &
          SLL_PERIODIC)
        A1_interp2d => new_cubic_spline_2d_interpolator( &
          Nc_x1+1, &
          Nc_x2+1, &
          x1_min, &
          x1_max, &
          x2_min, &
          x2_max, &
          SLL_HERMITE, &
          SLL_PERIODIC)
        A2_interp2d => new_cubic_spline_2d_interpolator( &
          Nc_x1+1, &
          Nc_x2+1, &
          x1_min, &
          x1_max, &
          x2_min, &
          x2_max, &
          SLL_HERMITE, &
          SLL_PERIODIC)       
      case default
        print *,'#bad interp2d_case',interp2d_case
        print *,'#not implemented'
        print *,'#in initialize_guiding_center_2d_polar'
        stop
    end select

    select case (interp1d_x1_case)
      case ("SLL_CUBIC_SPLINES")
        A1_interp1d_x1 => new_cubic_spline_1d_interpolator( &
          Nc_x1+1, &
          x1_min, &
          x1_max, &
          SLL_HERMITE)
        A2_interp1d_x1 => new_cubic_spline_1d_interpolator( &
          Nc_x1+1, &
          x1_min, &
          x1_max, &
          SLL_HERMITE)
      case default
        print *,'#bad interp1d_x1_case',interp1d_x1_case
        print *,'#not implemented'
        print *,'#in initialize_guiding_center_2d_polar'
        stop
    end select


    select case(charac2d_case)
      case ("SLL_EULER")
        charac2d => new_explicit_euler_2d_charac(&
          Nc_x1+1, &
          Nc_x2+1, &
          SLL_SET_TO_LIMIT, &
          SLL_PERIODIC)    
      case ("SLL_VERLET")      
        charac2d => new_verlet_2d_charac(&
          Nc_x1+1, &
          Nc_x2+1, &
          A1_interp2d, &
          A2_interp2d, &
          A1_interp1d_x1, &
          A2_interp1d_x1, &
          bc_type_1=SLL_SET_TO_LIMIT, &
          bc_type_2=SLL_PERIODIC)
      case default
        print *,'#bad charac2d_case',charac2d_case
        print *,'#not implemented'
        print *,'#in initialize_guiding_center_2d_polar'
        stop
    end select

  




    select case(advect2d_case)
      case ("SLL_BSL")
        sim%advect_2d => new_BSL_2d_advector(&
          interp2d, &
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
        sim%params(4) = k_mode
      case default
        print *,'#bad initial_function_case',initial_function_case
        print *,'#not implemented'
        print *,'#in initialize_guiding_center_2d_polar'
        stop
    end select
    !poisson solver
    ! for the moment no choice
    
    sim%poisson => new_plan_poisson_polar( &
      sim%mesh_2d%delta_eta1,& 
      x1_min, &
      Nc_x1, &
      Nc_x2, &
      (/ SLL_NEUMANN_MODE_0,SLL_DIRICHLET/))

    
  
  end subroutine initialize_guiding_center_2d_polar
  


  subroutine init_fake(sim, filename)
    class(sll_simulation_2d_guiding_center_polar), intent(inout) :: sim
    character(len=*), intent(in)                                :: filename
  
    print *,'# Do not use the routine init_vp4d_fake'
    print *,'#use instead initialize_vlasov_par_poisson_seq_cart'
    stop
  
  end subroutine init_fake
  
  subroutine run_gc2d_polar(sim)
    class(sll_simulation_2d_guiding_center_polar), intent(inout) :: sim
    
    
    
    
    
    
    
    
    
    print *,'#not implemented for the moment!'
  end subroutine run_gc2d_polar    



end module sll_simulation_2d_guiding_center_polar_module
