module sll_m_sim_bsl_ad_2d0v_curv

!2d analytic field curvilinear simulation
!contact: Adnane Hamiaz (hamiaz@math.unistra.fr
!         Michel Mehrenberger (mehrenbe@math.unistra.fr)

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_errors.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_advection_1d_base, only: &
    sll_c_advection_1d_base

  use sll_m_advection_1d_bsl, only: &
    sll_f_new_bsl_1d_advector

  use sll_m_advection_1d_csl_periodic, only: &
    sll_f_new_csl_periodic_1d_advector

  use sll_m_advection_2d_base, only: &
    sll_c_advection_2d_base

  use sll_m_advection_2d_bsl, only: &
    sll_f_new_bsl_2d_advector

  use sll_m_advection_2d_tensor_product, only: &
    sll_f_new_tensor_product_2d_advector

  use sll_m_ascii_io, only: &
    sll_s_ascii_file_create

  use sll_m_boundary_condition_descriptors, only: &
    sll_p_dirichlet, &
    sll_p_hermite, &
    sll_p_periodic, &
    sll_p_set_to_limit

  use sll_m_cartesian_meshes, only: &
    sll_o_get_node_positions, &
    sll_f_new_cartesian_mesh_1d, &
    sll_t_cartesian_mesh_1d, &
    sll_t_cartesian_mesh_2d, &
    operator(*)

  use sll_m_characteristics_1d_base, only: &
    sll_c_characteristics_1d_base

  use sll_m_characteristics_1d_explicit_euler, only: &
    sll_f_new_explicit_euler_1d_charac

  use sll_m_characteristics_1d_trapezoid, only: &
    sll_f_new_trapezoid_1d_charac

  use sll_m_characteristics_2d_base, only: &
    sll_i_signature_process_outside_point, &
    sll_c_characteristics_2d_base

  use sll_m_characteristics_2d_explicit_euler, only: &
    sll_f_new_explicit_euler_2d_charac

  use sll_m_characteristics_2d_verlet, only: &
    sll_f_new_verlet_2d_charac

  use sll_m_common_array_initializers, only: &
    sll_f_constant_time_initializer_1d, &
    sll_f_cos_bell0_initializer_2d, &
    sll_f_cos_bell_initializer_2d, &
    sll_f_cos_sin_initializer_2d, &
    sll_f_gaussian_initializer_2d, &
    sll_f_khp1_2d, &
    sll_f_one_initializer_2d, &
    sll_f_rotation_a1_exact_charac_2d, &
    sll_f_rotation_a1_initializer_2d, &
    sll_f_rotation_a2_exact_charac_2d, &
    sll_f_rotation_a2_initializer_2d, &
    sll_f_rotation_phi_initializer_2d, &
    sll_i_scalar_initializer_1d, &
    sll_i_scalar_initializer_2d, &
    sll_i_scalar_initializer_4d, &
    sll_f_sdf_a1_exact_charac_2d, &
    sll_f_sdf_a1_initializer_2d, &
    sll_f_sdf_a2_exact_charac_2d, &
    sll_f_sdf_a2_initializer_2d, &
    sll_f_sdf_phi_initializer_2d, &
    sll_f_sdf_time_initializer_1d, &
    sll_f_translation_a1_exact_charac_2d, &
    sll_f_translation_a1_initializer_2d, &
    sll_f_translation_a2_exact_charac_2d, &
    sll_f_translation_a2_initializer_2d, &
    sll_f_translation_phi_initializer_2d

  use sll_m_common_coordinate_transformations, only: &
    sll_f_identity_jac11, &
    sll_f_identity_jac12, &
    sll_f_identity_jac21, &
    sll_f_identity_jac22, &
    sll_f_identity_x1, &
    sll_f_identity_x2, &
    sll_f_polar_jac11, &
    sll_f_polar_jac12, &
    sll_f_polar_jac21, &
    sll_f_polar_jac22, &
    sll_f_polar_x1, &
    sll_f_polar_x2, &
    sll_f_sinprod_jac11, &
    sll_f_sinprod_jac12, &
    sll_f_sinprod_jac21, &
    sll_f_sinprod_jac22, &
    sll_f_sinprod_x1, &
    sll_f_sinprod_x2

  use sll_m_constants, only: &
    sll_p_pi

  use sll_m_coordinate_transformation_2d_base, only: &
    sll_c_coordinate_transformation_2d_base

  use sll_m_coordinate_transformations_2d, only: &
    sll_f_new_coordinate_transformation_2d_analytic

  use sll_m_cubic_spline_interpolator_1d, only: &
    sll_f_new_cubic_spline_interpolator_1d

  use sll_m_cubic_spline_interpolator_2d, only: &
    sll_f_new_cubic_spline_interpolator_2d

  use sll_m_hdf5_io_serial, only: &
    sll_o_hdf5_file_close, &
    sll_o_hdf5_file_create, &
    sll_o_hdf5_write_array

  use sll_m_hermite_interpolation_1d, only: &
    sll_p_hermite_1d_c0

  use sll_m_hermite_interpolation_2d, only: &
    sll_s_compute_w_hermite, &
    sll_p_hermite_c0, &
    sll_p_hermite_periodic

  use sll_m_hermite_interpolator_1d, only: &
    sll_f_new_hermite_interpolator_1d

  use sll_m_hermite_interpolator_2d, only: &
    sll_f_new_hermite_interpolator_2d

  use sll_m_interpolators_1d_base, only: &
    sll_c_interpolator_1d

  use sll_m_interpolators_2d_base, only: &
    sll_c_interpolator_2d

  use sll_m_operator_splitting, only: &
    sll_s_do_split_steps, &
    sll_p_strang_tvt

  use sll_m_reduction, only: &
    sll_f_compute_integral_trapezoid_1d

  use sll_m_sim_base, only: &
    sll_c_simulation_base_class

  use sll_m_split_advection_2d, only: &
    sll_f_new_split_advection_2d, &
    sll_p_advective, &
    sll_p_conservative, &
    sll_t_split_advection_2d

  use sll_m_utilities, only: &
    sll_s_int2string

  use sll_m_xdmf, only: &
    sll_s_xdmf_close, &
    sll_o_xdmf_open, &
    sll_o_xdmf_write_array

#ifdef MUDPACK
  use sll_m_poisson_2d_mudpack_curvilinear_solver_old
#endif
  implicit none

  public :: &
    sll_s_delete_analytic_field_2d_curvilinear, &
    sll_f_new_analytic_field_2d_curvilinear, &
    sll_t_simulation_2d_analytic_field_curvilinear

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  
  sll_int32, parameter :: SLL_EULER = 0 
  sll_int32, parameter :: SLL_PREDICTOR_CORRECTOR = 1 
  sll_int32, parameter :: SLL_LEAP_FROG = 2 
  sll_int32, parameter :: SLL_SPLITTING = 3 
  sll_int32, parameter :: SLL_COMPUTE_FIELD_FROM_PHI = 0 
  sll_int32, parameter :: SLL_COMPUTE_FIELD_FROM_ANALYTIC = 1 
  sll_int32, parameter :: SLL_COMPUTE_FIELD_FROM_PHI_FD2 = 2 
  sll_int32, parameter :: SLL_COMPUTE_FIELD_FROM_PHI_FD3 = 3 
  sll_int32, parameter :: SLL_COMPUTE_FIELD_FROM_PHI_FD4 = 4 
  sll_int32, parameter :: SLL_COMPUTE_FIELD_FROM_PHI_FD5 = 5 
  sll_int32, parameter :: SLL_COMPUTE_FIELD_FROM_PHI_FD6 = 6 
  sll_int32, parameter :: SLL_COMPUTE_FIELD_FROM_PHI_FD = 7 

  type, extends(sll_c_simulation_base_class) :: &
    sll_t_simulation_2d_analytic_field_curvilinear

   
   !geometry
   type(sll_t_cartesian_mesh_2d), pointer :: mesh_2d
   sll_real64, dimension(:), pointer :: eta1_array
   sll_real64, dimension(:), pointer :: eta2_array
   sll_int32 :: num_dof1  
   sll_int32 :: num_dof2
   sll_real64, dimension(:), pointer :: integration_weight1
   sll_real64, dimension(:), pointer :: integration_weight2
   
     
   !transformation 
    class(sll_c_coordinate_transformation_2d_base), pointer :: transformation
  
   !initial function
   procedure(sll_i_scalar_initializer_2d), nopass, pointer :: init_func
   sll_real64, dimension(:), pointer :: params
   !set boundary functions for characteristics
   procedure(sll_i_signature_process_outside_point), nopass, pointer :: process_outside_point1_func
   procedure(sll_i_signature_process_outside_point), nopass, pointer :: process_outside_point2_func
   
      
    !advector
   sll_int32 :: advection_form
   class(sll_c_advection_2d_base), pointer    :: advect_2d
   class(sll_c_advection_1d_base), pointer    :: advect1_1d
   class(sll_c_advection_1d_base), pointer    :: advect2_1d
   class(sll_c_characteristics_1d_base), pointer :: charac1
   class(sll_c_characteristics_1d_base), pointer :: charac2
   class(sll_c_interpolator_1d), pointer   :: interp1
   class(sll_c_interpolator_1d), pointer   :: interp2


   procedure(sll_i_scalar_initializer_2d), nopass, pointer :: phi_func
   procedure(sll_i_scalar_initializer_2d), nopass, pointer :: A1_func
   procedure(sll_i_scalar_initializer_2d), nopass, pointer :: A2_func
   procedure(sll_i_scalar_initializer_4d), nopass, pointer :: A1_exact_charac_func
   procedure(sll_i_scalar_initializer_4d), nopass, pointer :: A2_exact_charac_func
   sll_real64, dimension(:), pointer :: A_func_params
   procedure(sll_i_scalar_initializer_1d), nopass, pointer :: A_time_func
   sll_real64, dimension(:), pointer :: A_time_func_params
   class(sll_t_split_advection_2d), pointer :: split
   logical :: csl_2012
   
   !interpolator for derivatives
   class(sll_c_interpolator_2d), pointer   :: phi_interp2d

    
   !time_iterations
   sll_real64 :: dt
   sll_int32  :: num_iterations
   sll_int32  :: freq_diag
   sll_int32  :: freq_diag_time
   character(len=256)      :: thdiag_filename

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
   
   sll_int32 :: fd_degree1 
   sll_int32 :: fd_degree2 
  contains
    procedure, pass(sim) :: run => run_af2d_curvilinear
    procedure, pass(sim) :: init_from_file => init_fake
     
  end type sll_t_simulation_2d_analytic_field_curvilinear





contains

  function sll_f_new_analytic_field_2d_curvilinear( &
    filename, &
    num_run &
    ) result(sim)
    type(sll_t_simulation_2d_analytic_field_curvilinear), pointer :: sim 
    character(len=*), intent(in), optional :: filename
    sll_int32, intent(in), optional :: num_run   
    sll_int32 :: ierr
    
    SLL_ALLOCATE(sim,ierr)
    
    call initialize_analytic_field_2d_curvilinear(sim,filename,num_run)
    
  
  
  end function sll_f_new_analytic_field_2d_curvilinear
  
  subroutine initialize_analytic_field_2d_curvilinear( &
    sim, &
    filename, &
    num_run)
    class(sll_t_simulation_2d_analytic_field_curvilinear), intent(inout) :: sim
    character(len=*), intent(in), optional :: filename
    sll_int32, intent(in), optional :: num_run
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
    sll_real64 :: cutx_value
    sll_real64 :: cutf_value

    
     !time_iterations
    sll_real64 :: dt
    sll_int32 :: number_iterations
    sll_int32 :: freq_diag
    sll_int32 :: freq_diag_time
    character(len=256) :: time_loop_case


    !advector
    character(len=256) :: advection_form    
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
    type(sll_t_cartesian_mesh_1d), pointer :: mesh_x1
    type(sll_t_cartesian_mesh_1d), pointer :: mesh_x2
    class(sll_c_interpolator_2d), pointer :: f_interp2d
    class(sll_c_interpolator_2d), pointer :: phi_interp2d
    class(sll_c_characteristics_2d_base), pointer :: charac2d
    !class(sll_c_characteristics_1d_base), pointer :: charac1d_x1
    !class(sll_c_characteristics_1d_base), pointer :: charac1d_x2
    class(sll_c_interpolator_2d), pointer   :: A1_interp2d
    class(sll_c_interpolator_2d), pointer   :: A2_interp2d
    class(sll_c_interpolator_1d), pointer   :: A1_interp1d_x1
    class(sll_c_interpolator_1d), pointer   :: A2_interp1d_x1
    class(sll_c_interpolator_1d), pointer   :: A2_interp1d_x2
    !class(sll_c_interpolator_1d), pointer :: f_interp1d_x1
    !class(sll_c_interpolator_1d), pointer :: f_interp1d_x2
    !class(sll_c_advection_1d_base), pointer    :: advect_1d_x1
    !class(sll_c_advection_1d_base), pointer    :: advect_1d_x2
    sll_int32 :: ierr
    sll_real64, dimension(4) :: params_mesh
    !sll_real64 :: eta1_min_bis
    !sll_real64 :: eta1_max_bis
    !sll_real64 :: eta2_min_bis
    !sll_real64 :: eta2_max_bis
    !sll_int32  :: Nc_eta1_bis
    !sll_int32  :: Nc_eta2_bis
    sll_real64 :: A1 !parameter for translation
    sll_real64 :: A2 !parameter for translation
    sll_int32 :: hermite_degree1
    sll_int32 :: hermite_degree2
    sll_int32 :: fd_degree1
    sll_int32 :: fd_degree2
    character(len=256)      :: str_num_run
    character(len=256)      :: filename_loc
    logical :: csl_2012
    logical :: feet_inside1
    logical :: feet_inside2

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
      sigma_2, &
      cutx_value, &
      cutf_value
      
    namelist /time_iterations/ &
      dt, &
      number_iterations, &
      freq_diag, &
      freq_diag_time, &
      time_loop_case

     namelist /advector/ &
      advection_form, &
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
      f_interp1d_x1_case, & 
      f_interp1d_x2_case, & 
      time_period, &
      A1, &
      A2, &
      compute_field_case, &
      hermite_degree1, &
      hermite_degree2, &
      fd_degree1, &
      fd_degree2, &
      csl_2012
   
      
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
    eta1_max = 2._f64*sll_p_pi
    eta2_min = 0.0_f64
    eta2_max = 2._f64*sll_p_pi
    alpha1 = 0._f64
    alpha2 = 0._f64
    params_mesh = (/ alpha1, alpha2, eta1_max-eta1_min, eta2_max-eta2_min/)
    
    !initial function
    initial_function_case="SLL_KHP1"
    kmode_eta1 = 0.5_f64
    kmode_eta2 = 1._f64
    eps = 1._f64 !0.015_f64
    cutx_value = 0.25_f64
    cutf_value = 1._f64
    
    !time_iterations
    dt = 0.1_f64
    number_iterations  = 600
    freq_diag = 100
    freq_diag_time = 1
    !time_loop_case = "SLL_EULER"
    time_loop_case = "SLL_PREDICTOR_CORRECTOR" 

    !advector
    advection_form = "sll_p_advective"
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
    hermite_degree1 = 4
    hermite_degree2 = 4
    fd_degree1 = 4
    fd_degree2 = 4
    csl_2012 = .false.
    
    !boundaries conditions
    sim%bc_eta1_left = sll_p_periodic
    sim%bc_eta1_right= sll_p_periodic
    sim%bc_eta2_left = sll_p_periodic
    sim%bc_eta2_right= sll_p_periodic 
    sim%bc_interp2d_eta1 = sll_p_periodic
    sim%bc_interp2d_eta2 = sll_p_periodic
    sim%bc_charac2d_eta1 = sll_p_periodic
    sim%bc_charac2d_eta2 = sll_p_periodic



    if(present(num_run))then
      !call sll_s_int2string(num_run, str_num_run)
      write(str_num_run, *) num_run
      str_num_run = adjustl(str_num_run) 
      sim%thdiag_filename = "thdiag_"//trim(str_num_run)//".dat"
    else      
      sim%thdiag_filename = "thdiag.dat"
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
          print *, '#in initialize_analytic_field_2d_curvilinear() failed to open file ', &
          trim(filename_loc)//'.nml'
          STOP
        end if
      print *,'#initialization with filename:'
      print *,'#',trim(filename_loc)//'.nml'
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
   
    sim%fd_degree1 = fd_degree1
    sim%fd_degree2 = fd_degree2
    
    sim%csl_2012 = csl_2012
   
        sim%num_dof1 =  num_cells_eta1+1
        sim%num_dof2 =  num_cells_eta2+1
    select case(advection_form)
      case ("sll_p_advective")
        sim%advection_form = sll_p_advective
      case ("sll_p_conservative")
        sim%advection_form = sll_p_conservative
        !sim%num_dof1 =  num_cells_eta1
        !sim%num_dof2 =  num_cells_eta2
      case default
        print *,'#bad advection_form',advection_form
        print *,'#not implemented'
        print *,'#in initialize_analytic_field_2d_curvilinear'
        stop
    end select
    
      
   
    select case(bc_eta1_left)
      case ("SLL_PERIODIC")
        print*,"#bc_eta1_left = sll_p_periodic" 
        sim%bc_eta1_left = sll_p_periodic
      case ("SLL_DIRICHLET")
        print*,"#bc_eta1_left = sll_p_dirichlet"  
        sim%bc_eta1_left = sll_p_dirichlet
      case default
        print *,'#bad bc_eta1_left',bc_eta1_left
        print *,'#not implemented'
        print *,'#in initialize_analytic_field_2d_curvilinear'
        stop
    end select
    
    select case(bc_eta1_right)
      case ("SLL_PERIODIC")
        print*,"#bc_eta1_right = sll_p_periodic" 
        sim%bc_eta1_right = sll_p_periodic
      case ("SLL_DIRICHLET")
        print*,"#bc_eta1_right = sll_p_dirichlet"  
        sim%bc_eta1_right = sll_p_dirichlet
      case default
        print *,'#bad bc_eta1_right',bc_eta1_right
        print *,'#not implemented'
        print *,'#in initialize_analytic_field_2d_curvilinear'
        stop
    end select
    
       select case(bc_eta2_left)
      case ("SLL_PERIODIC")
        print*,"#bc_eta2_left = sll_p_periodic" 
        sim%bc_eta2_left = sll_p_periodic
      case ("SLL_DIRICHLET")
        print*,"#bc_eta2_left = sll_p_dirichlet"  
        sim%bc_eta2_left = sll_p_dirichlet
      case default
        print *,'#bad bc_eta2_left',bc_eta2_left
        print *,'#not implemented'
        print *,'#in initialize_analytic_field_2d_curvilinear'
        stop
    end select
    
    select case(bc_eta2_right)
      case ("SLL_PERIODIC")
        print*,"#bc_eta2_right = sll_p_periodic" 
        sim%bc_eta2_right = sll_p_periodic
      case ("SLL_DIRICHLET")
        print*,"#bc_eta2_right = sll_p_dirichlet"  
        sim%bc_eta2_right = sll_p_dirichlet
      case default
        print *,'#bad bc_eta2_right',bc_eta2_right
        print *,'#not implemented'
        print *,'#in initialize_analytic_field_2d_curvilinear'
        stop
    end select
    
    select case(bc_interp2d_eta1)
      case ("SLL_PERIODIC")
        print*,"#bc_interp2d_eta1= sll_p_periodic" 
        sim%bc_interp2d_eta1 = sll_p_periodic
      case ("SLL_DIRICHLET")
        print*,"#bc_interp2d_eta1 = sll_p_dirichlet"  
        sim%bc_interp2d_eta1= sll_p_dirichlet
      case ("SLL_HERMITE")
        print*,"#bc_interp2d_eta1 = sll_p_hermite"  
        sim%bc_interp2d_eta1= sll_p_hermite 
      case default
        print *,'#bad bc_interp2d_eta1',bc_interp2d_eta1
        print *,'#not implemented'
        print *,'#in initialize_analytic_field_2d_curvilinear'
        stop
    end select

    select case(bc_interp2d_eta2)
      case ("SLL_PERIODIC")
        print*,"#bc_interp2d_eta2= sll_p_periodic" 
        sim%bc_interp2d_eta2 = sll_p_periodic
      case ("SLL_DIRICHLET")
        print*,"#bc_interp2d_eta2 = sll_p_dirichlet"  
        sim%bc_interp2d_eta2= sll_p_dirichlet
      case ("SLL_HERMITE")
        print*,"#bc_interp2d_eta2 = sll_p_hermite"  
        sim%bc_interp2d_eta2= sll_p_hermite 
      case default
        print *,'#bad bc_interp2d_eta2',bc_interp2d_eta2
        print *,'#not implemented'
        print *,'#in initialize_analytic_field_2d_curvilinear'
        stop
    end select
    
    select case(bc_charac2d_eta1)
      case ("SLL_PERIODIC")
        print*,"#bc_charac2d_eta1= sll_p_periodic" 
        sim%bc_charac2d_eta1= sll_p_periodic
        sim%process_outside_point1_func => process_outside_point_periodic1            
      !case ("SLL_DIRICHLET")
      !  print*,"#bc_charac2d_eta1 = sll_p_dirichlet"  
      !  sim%bc_charac2d_eta1= sll_p_dirichlet
      !case ("SLL_HERMITE")
      !  print*,"#bc_charac2d_eta1 = sll_p_hermite"  
      !  sim%bc_charac2d_eta1= sll_p_hermite 
      case ("sll_p_set_to_limit")
        print*,"#bc_charac2d_eta1 = sll_p_set_to_limit"  
        sim%bc_charac2d_eta1= sll_p_set_to_limit
        sim%process_outside_point1_func => process_outside_point_set_to_limit1    
      case default
        print *,'#bad bc_charac2d_eta1',bc_charac2d_eta1
        print *,'#not implemented'
        print *,'#in initialize_analytic_field_2d_curvilinear'
        stop
    end select

    select case(bc_charac2d_eta2)
      case ("SLL_PERIODIC")
        print*,"#bc_charac2d_eta2= sll_p_periodic" 
        sim%bc_charac2d_eta2= sll_p_periodic
        sim%process_outside_point2_func => process_outside_point_periodic1
      !case ("SLL_DIRICHLET")
      !  print*,"#bc_charac2d_eta2= sll_p_dirichlet"  
      !  sim%bc_charac2d_eta2= sll_p_dirichlet
      !case ("SLL_HERMITE")
      !  print*,"#bc_charac2d_eta2 = sll_p_hermite"  
      !  sim%bc_charac2d_eta2= sll_p_hermite 
      case ("sll_p_set_to_limit")
        print*,"#bc_charac2d_eta2 = sll_p_set_to_limit"  
        sim%bc_charac2d_eta2= sll_p_set_to_limit
        sim%process_outside_point2_func => process_outside_point_set_to_limit1   
      case default
        print *,'#bad bc_charac2d_eta2',bc_charac2d_eta2
        print *,'#not implemented'
        print *,'#in initialize_analytic_field_2d_curvilinear'
        stop
    end select
    
    
    select case (mesh_case_eta1)
      case ("SLL_CARTESIAN_MESH")
        mesh_x1 => sll_f_new_cartesian_mesh_1d(Nc_eta1,eta_min=eta1_min, eta_max=eta1_max)  
        call sll_o_get_node_positions( mesh_x1, sim%eta1_array )
      case default
        print*,'#mesh_case_eta1', mesh_case_eta1, ' not implemented'
        stop 
    end select
    select case (mesh_case_eta2)
      case ("SLL_CARTESIAN_MESH")
        mesh_x2 => sll_f_new_cartesian_mesh_1d(Nc_eta2,eta_min=eta2_min, eta_max=eta2_max)
        call sll_o_get_node_positions( mesh_x2, sim%eta2_array )
      case default
        print*,'#mesh_case_eta2', mesh_case_eta2, ' not implemented'
        stop 
    end select
    sim%mesh_2d =>  mesh_x1 * mesh_x2 ! tensor product
    !  In collela  mesh params_mesh =( alpha1, alpha2, L1, L2 ) such that :
    !  x1= eta1 + alpha1*sin(2*pi*eta1/L1)*sin(2*pi*eta2/L2)
    params_mesh = (/ alpha1, alpha2, eta1_max-eta1_min, eta2_max-eta2_min/)
          
!   sim%mesh_2d => sll_f_new_cartesian_mesh_2d( &
!      Nc_eta1, &
!      Nc_eta2, &
!      eta1_min , &
!      eta1_max , &
!      eta2_min , &
!      eta2_max ) 
      
    select case (transf_case)
      case ("SLL_CARTESIAN")
        print*,'#transf_case =SLL_CARTESIAN'       
        sim%transformation => sll_f_new_coordinate_transformation_2d_analytic( &
         "analytic_identity_transformation", &
         sim%mesh_2d, &
         sll_f_identity_x1, &
         sll_f_identity_x2, &
         sll_f_identity_jac11, &
         sll_f_identity_jac12, &
         sll_f_identity_jac21, &
         sll_f_identity_jac22, &
         params_mesh   )  
      case ("SLL_POLAR") 
        print*,'#transf_case =SLL_POLAR'  
        sim%transformation => sll_f_new_coordinate_transformation_2d_analytic( &
         "analytic_polar_transformation", &
         sim%mesh_2d, &
         sll_f_polar_x1, &
         sll_f_polar_x2, &
         sll_f_polar_jac11, &
         sll_f_polar_jac12, &
         sll_f_polar_jac21, &
         sll_f_polar_jac22, &
         params_mesh  )     
      case ("SLL_COLLELA")
        print*,'#transf_case =SLL_COLLELA'  
        sim%transformation => sll_f_new_coordinate_transformation_2d_analytic( &
         "analytic_collela_transformation", &
         sim%mesh_2d, &
         sll_f_sinprod_x1, &
         sll_f_sinprod_x2, &
         sll_f_sinprod_jac11, &
         sll_f_sinprod_jac12, &
         sll_f_sinprod_jac21, &
         sll_f_sinprod_jac22, &
         params_mesh  )  
        case default
        print *,'#bad transf_case',transf_case
        print *,'#not implemented'
        print *,'#in initialize_analytic_field_2d_curvilinear'
        stop
    end select  


    SLL_ALLOCATE(sim%integration_weight1(sim%num_dof1),ierr)
    SLL_ALLOCATE(sim%integration_weight2(sim%num_dof2),ierr)
    
    call compute_integration_weight(sim%eta1_array,sim%integration_weight1)
    call compute_integration_weight(sim%eta2_array,sim%integration_weight2)


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



     
!    select case(advect1d_x1_case)
!      case ("SLL_BSL")
!        eta1_min_bis = eta1_min
!        eta1_max_bis = eta1_max
!        Nc_eta1_bis = Nc_eta1
!      case ("SLL_PSM")
!        eta1_min_bis = eta1_min
!        eta1_max_bis = eta1_max
!        Nc_eta1_bis = Nc_eta1
!      case ("SLL_CSL")
!        eta1_min_bis = eta1_min-0.5_f64*mesh_x1%delta_eta
!        eta1_max_bis = eta1_max-0.5_f64*mesh_x1%delta_eta
!        Nc_eta1_bis = Nc_eta1
!      case default
!        print *,'#bad value of advect1d_x1_case'
!        stop  
!    end select

!    select case(advect1d_x2_case)
!      case ("SLL_BSL")
!        eta2_min_bis = eta2_min
!        eta2_max_bis = eta2_max
!        Nc_eta2_bis = Nc_eta2
!      case ("SLL_PSM")
!        eta2_min_bis = eta2_min
!        eta2_max_bis = eta2_max
!        Nc_eta2_bis = Nc_eta2
!      case ("SLL_CSL")
!        eta2_min_bis = eta2_min-0.5_f64*mesh_x2%delta_eta
!        eta2_max_bis = eta2_max-0.5_f64*mesh_x2%delta_eta        
!        Nc_eta2_bis = Nc_eta2
!      case default
!        print *,'#bad value of advect1d_x2_case'
!        stop  
!    end select
  
    select case (f_interp2d_case)
      case ("SLL_CUBIC_SPLINES")
        print*,"#f interpolation SLL_CUBIC_SPLINES"
        f_interp2d => sll_f_new_cubic_spline_interpolator_2d( &
          Nc_eta1+1, &
          Nc_eta2+1, &
          eta1_min, &
          eta1_max, &
          eta2_min, &
          eta2_max, &
          sim%bc_interp2d_eta1, &
          sim%bc_interp2d_eta2)
      case ("SLL_HERMITE")
        print*,"#f interpolation sll_p_hermite"
        f_interp2d => sll_f_new_hermite_interpolator_2d( &
          Nc_eta1+1, &
          Nc_eta2+1, &
          eta1_min, &
          eta1_max, &
          eta2_min, &
          eta2_max, &
          hermite_degree1, &          
          hermite_degree2, &          
          sll_p_hermite_c0, &
          sll_p_hermite_c0, &
          sll_p_hermite_periodic, &
          sll_p_hermite_periodic)

      case default
        print *,'#bad f_interp2d_case',f_interp2d_case
        print *,'#not implemented'
        print *,'#in initialize_analytic_field_2d_curvilinear'
        stop
    end select




    select case (A_interp_case)
      case ("SLL_CUBIC_SPLINES")
       print*,"#A1_2d interpolation SLL_CUBIC_SPLINES"
        A1_interp2d => sll_f_new_cubic_spline_interpolator_2d( &
          Nc_eta1+1, &
          Nc_eta2+1, &
          eta1_min, &
          eta1_max, &
          eta2_min, &
          eta2_max, &
          sim%bc_interp2d_eta1, &
          sim%bc_interp2d_eta2)
       print*,"#A2_2d interpolation SLL_CUBIC_SPLINES"   
        A2_interp2d => sll_f_new_cubic_spline_interpolator_2d( &
          Nc_eta1+1, &
          Nc_eta2+1, &
          eta1_min, &
          eta1_max, &
          eta2_min, &
          eta2_max, &
          sim%bc_interp2d_eta1, &
          sim%bc_interp2d_eta2)  
       print*,"#A1_1d interpolation SLL_CUBIC_SPLINES"   
        A1_interp1d_x1 => sll_f_new_cubic_spline_interpolator_1d( &
          Nc_eta1+1, &
          eta1_min, &
          eta1_max, &
          sim%bc_interp2d_eta1)
       print*,"#A2_1d interpolation SLL_CUBIC_SPLINES"     
        A2_interp1d_x1 => sll_f_new_cubic_spline_interpolator_1d( &
          Nc_eta1+1, &
          eta1_min, &
          eta1_max, &
          sim%bc_interp2d_eta1)
        A2_interp1d_x2 => sll_f_new_cubic_spline_interpolator_1d( &
          Nc_eta2+1, &
          eta2_min, &
          eta2_max, &
          sim%bc_interp2d_eta2)
      case default
        print *,'#bad A_interp_case',A_interp_case
        print *,'#not implemented'
        print *,'#in initialize_analytic_field_2d_curvilinear'
        stop
    end select

    select case (phi_interp2d_case)
      case ("SLL_CUBIC_SPLINES")
      print*,"#phi interpolation SLL_CUBIC_SPLINES"  
        phi_interp2d => sll_f_new_cubic_spline_interpolator_2d( &
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
!        f_interp1d_x1 => sll_f_new_cubic_spline_interpolator_1d( &
!          Nc_eta1+1, &
!          eta1_min, &
!          eta1_max, &
!          !Nc_eta1_bis+1, &
!          !eta1_min_bis, &
!          !eta1_max_bis, &
!          sim%bc_interp2d_eta1)
        sim%interp1 => sll_f_new_cubic_spline_interpolator_1d( &
          Nc_eta1+1, &
          eta1_min, &
          eta1_max, &
          !Nc_eta1_bis+1, &
          !eta1_min_bis, &
          !eta1_max_bis, &
          sim%bc_interp2d_eta1)
      case ("SLL_HERMITE")
        sim%interp1 => sll_f_new_hermite_interpolator_1d( &
          Nc_eta1+1, &
          eta1_min, &
          eta1_max, &
          hermite_degree1, &          
          sll_p_hermite_1d_c0, &
          sim%bc_interp2d_eta1) 
    
      case default
        print *,'#bad f_interp1d_x1_case',f_interp1d_x1_case
        print *,'#not implemented'
        print *,'#in initialize_guiding_center_2d_curvilinear'
        stop
    end select


    select case (f_interp1d_x2_case)
      case ("SLL_CUBIC_SPLINES")
!        f_interp1d_x2 => sll_f_new_cubic_spline_interpolator_1d( &
!          Nc_eta2+1, &
!          eta2_min, &
!          eta2_max, &
!          !Nc_eta2_bis+1, &
!          !eta2_min_bis, &
!          !eta2_max_bis, &
!          sim%bc_interp2d_eta2)
        sim%interp2 => sll_f_new_cubic_spline_interpolator_1d( &
          Nc_eta2+1, &
          eta2_min, &
          eta2_max, &
          !Nc_eta2_bis+1, &
          !eta2_min_bis, &
          !eta2_max_bis, &
          sim%bc_interp2d_eta2)
      case ("SLL_HERMITE")
        sim%interp2 => sll_f_new_hermite_interpolator_1d( &
          Nc_eta2+1, &
          eta2_min, &
          eta2_max, &
          hermite_degree2, &          
          sll_p_hermite_1d_c0, &
          sim%bc_interp2d_eta2) 
      case default
        print *,'#bad f_interp1d_x2_case',f_interp1d_x2_case
        print *,'#not implemented'
        print *,'#in initialize_guiding_center_2d_curvilinear'
        stop
    end select


    select case(charac1d_x1_case)
      case ("SLL_EULER")
        sim%charac1 => sll_f_new_explicit_euler_1d_charac(&
          Nc_eta1+1, &
          eta_min=eta1_min, &
          eta_max=eta1_max, &
          !Nc_eta1_bis+1, &
          !eta_min=eta1_min_bis, &
          !eta_max=eta1_max_bis, &
          bc_type= sim%bc_charac2d_eta1, &
          feet_inside = feet_inside1)    
      case ("SLL_TRAPEZOID")
        sim%charac1 => &
          sll_f_new_trapezoid_1d_charac(&
          !Nc_eta1_bis+1, &
          Nc_eta1+1, &
          A1_interp1d_x1, &
          bc_type= sim%bc_charac2d_eta1, &
          !eta_min=eta1_min_bis, &
          !eta_max=eta1_max_bis)
          eta_min=eta1_min, &
          eta_max=eta1_max, &
          feet_inside = feet_inside1)
!      case ("SLL_EULER_CONSERVATIVE")
!        charac1d_x1 => new_explicit_euler_conservative_1d_charac(&
!          Nc_eta1_bis+1, &
!          eta_min=eta1_min_bis, &
!          eta_max=eta1_max_bis, &
!          bc_type=sim%bc_charac2d_eta1)    
!      case ("SLL_TRAPEZOID_CONSERVATIVE")
!        charac1d_x1 => &
!          new_trapezoid_conservative_1d_charac(&
!          Nc_eta1_bis+1, &
!          A1_interp1d_x1, &
!          eta_min=eta1_min_bis, &
!          eta_max=eta1_max_bis, &
!          bc_type=sim%bc_charac2d_eta1)
      case default
        print *,'#bad charac1d_x1_case',charac1d_x1_case
        print *,'#not implemented'
        print *,'#in initialize_guiding_center_2d_curvilinear'
        stop
    end select

    select case(charac1d_x2_case)
      case ("SLL_EULER")
        sim%charac2 => sll_f_new_explicit_euler_1d_charac(&
          Nc_eta2+1, &
          eta_min=eta2_min, &
          eta_max=eta2_max, &
          !Nc_eta2_bis+1, &
          !eta_min=eta2_min_bis, &
          !eta_max=eta2_max_bis, &
          bc_type= sim%bc_charac2d_eta2, &
          feet_inside = feet_inside2)    
      case ("SLL_TRAPEZOID")
        sim%charac2 => &
          sll_f_new_trapezoid_1d_charac(&
          !Nc_eta2_bis+1, &
          Nc_eta2+1, &
          A2_interp1d_x2, &
          bc_type= sim%bc_charac2d_eta2, &
          eta_min=eta2_min, &
          eta_max=eta2_max, &
          feet_inside = feet_inside2)
          !eta_min=eta2_min_bis, &
          !eta_max=eta2_max_bis)
!      case ("SLL_EULER_CONSERVATIVE")
!        charac1d_x2 => new_explicit_euler_conservative_1d_charac(&
!          Nc_eta2_bis+1, &
!          eta_min=eta2_min_bis, &
!          eta_max=eta2_max_bis, &
!          bc_type= sim%bc_charac2d_eta2)    
!      case ("SLL_TRAPEZOID_CONSERVATIVE")
!        charac1d_x2 => &
!          new_trapezoid_conservative_1d_charac(&
!          Nc_eta2_bis+1, &
!          A2_interp1d_x2, &
!          eta_min=eta2_min_bis, &
!          eta_max=eta2_max_bis, &
!          bc_type=sim%bc_charac2d_eta2)
      case default
        print *,'#bad charac1d_x2_case',charac1d_x2_case
        print *,'#not implemented'
        print *,'#in initialize_guiding_center_2d_curvilinear'
        stop
    end select

    select case(advect1d_x1_case)
      case ("SLL_BSL")
        sim%advect1_1d => sll_f_new_bsl_1d_advector(&
          sim%interp1, &
          sim%charac1, &
          Nc_eta1+1, &
          eta_min = eta1_min, &
          eta_max = eta1_max)
      case ("SLL_CSL_PERIODIC")
        select case (f_interp1d_x1_case)
          case ("SLL_CUBIC_SPLINES")
            sim%advect1_1d => sll_f_new_csl_periodic_1d_advector(&
              sim%interp1, &
              sim%charac1, &
              Nc_eta1+1, &
              eta_min = eta1_min, &
              eta_max = eta1_max)
          case ("SLL_HERMITE")   
            sim%advect1_1d => sll_f_new_csl_periodic_1d_advector(&
              sim%interp1, &
              sim%charac1, &
              Nc_eta1+1, &
              eta_min = eta1_min, &
              eta_max = eta1_max, &
              csl_degree = hermite_degree1+1)
          case default
            print *,'#bad f_interp1d_x1_case',f_interp1d_x1_case
        end select       
      case default
        print *,'#bad advect_case',advect1d_x1_case
        print *,'#not implemented'
        print *,'#in initialize_analytic_field_2d_curvilinear'
        stop
    end select

    select case(advect1d_x2_case)
      case ("SLL_BSL")
        sim%advect2_1d => sll_f_new_bsl_1d_advector(&
          sim%interp2, &
          sim%charac2, &
          Nc_eta2+1, &
          eta_min = eta2_min, &
          eta_max = eta2_max)
      case ("SLL_CSL_PERIODIC")
        select case (f_interp1d_x2_case)
          case ("SLL_CUBIC_SPLINES")
            sim%advect2_1d => sll_f_new_csl_periodic_1d_advector(&
              sim%interp2, &
              sim%charac2, &
              Nc_eta2+1, &
              eta_min = eta2_min, &
              eta_max = eta2_max)
          case ("SLL_HERMITE")
            sim%advect2_1d => sll_f_new_csl_periodic_1d_advector(&
              sim%interp2, &
              sim%charac2, &
              Nc_eta2+1, &
              eta_min = eta2_min, &
              eta_max = eta2_max, &
              csl_degree = hermite_degree2+1)
          case default
            print *,'#bad f_interp1d_x2_case',f_interp1d_x2_case
        end select       
      case default
        print *,'#bad advect_case',advect1d_x2_case
        print *,'#not implemented'
        print *,'#in initialize_analytic_field_2d_curvilinear'
        stop
    end select


    select case(charac2d_case)
      case ("SLL_EULER")
         print*,"#charac = SLL_EULER"  
        charac2d => sll_f_new_explicit_euler_2d_charac(&
          Nc_eta1+1, &
          Nc_eta2+1, &
          eta1_min=eta1_min, &
          eta1_max=eta1_max, &
          eta2_min=eta2_min, &
          eta2_max=eta2_max, &
          bc_type_1=sim%bc_charac2d_eta1, & !sll_p_set_to_limit, &
          bc_type_2=sim%bc_charac2d_eta2)    
      case ("SLL_VERLET")    
          print*,"#charac =SLL_VERLET"   
        charac2d => sll_f_new_verlet_2d_charac(&
          Nc_eta1+1, &
          Nc_eta2+1, &
          A1_interp2d, &
          A2_interp2d, &
          A1_interp1d_x1, &
          A2_interp1d_x1, &
          bc_type_1=sim%bc_charac2d_eta1, & !sll_p_set_to_limit, &
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
        sim%advect_2d => sll_f_new_bsl_2d_advector(&
          f_interp2d, &
          charac2d, &
          Nc_eta1+1, &
          Nc_eta2+1, &
          eta1_min = eta1_min, &
          eta1_max = eta1_max, &
          eta2_min = eta2_min, &
          eta2_max = eta2_max)
      case ("SLL_TENSOR_PRODUCT")
       print*,"#advect2d = SLL_TENSOR_PRODUCT" 
        sim%advect_2d => sll_f_new_tensor_product_2d_advector(&
          sim%advect1_1d, &
          sim%advect2_1d, &
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
        sim%phi_func => sll_f_sdf_phi_initializer_2d 
        sim%A1_func => sll_f_sdf_a1_initializer_2d 
        sim%A2_func => sll_f_sdf_a2_initializer_2d 
        sim%A1_exact_charac_func => sll_f_sdf_a1_exact_charac_2d 
        sim%A2_exact_charac_func => sll_f_sdf_a2_exact_charac_2d 
        SLL_ALLOCATE(sim%A_func_params(2),ierr)
        sim%A_time_func => sll_f_sdf_time_initializer_1d 
        SLL_ALLOCATE(sim%A_time_func_params(1),ierr)
        sim%A_time_func_params(1) = time_period
      case ("SLL_ROTATION_FLOW")
        sim%phi_func => sll_f_rotation_phi_initializer_2d 
        sim%A1_func => sll_f_rotation_a1_initializer_2d 
        sim%A2_func => sll_f_rotation_a2_initializer_2d 
        sim%A1_exact_charac_func => sll_f_rotation_a1_exact_charac_2d 
        sim%A2_exact_charac_func => sll_f_rotation_a2_exact_charac_2d 
        SLL_ALLOCATE(sim%A_func_params(2),ierr)
        sim%A_time_func => sll_f_constant_time_initializer_1d 
        SLL_ALLOCATE(sim%A_time_func_params(1),ierr)
        sim%A_time_func_params(1) = 1._f64
      case ("SLL_TRANSLATION_FLOW")
        sim%phi_func => sll_f_translation_phi_initializer_2d 
        sim%A1_func => sll_f_translation_a1_initializer_2d 
        sim%A2_func => sll_f_translation_a2_initializer_2d 
        sim%A1_exact_charac_func => sll_f_translation_a1_exact_charac_2d 
        sim%A2_exact_charac_func => sll_f_translation_a2_exact_charac_2d 
        SLL_ALLOCATE(sim%A_func_params(2),ierr)
        sim%A_func_params(1) = A1
        sim%A_func_params(2) = A2
        sim%A_time_func => sll_f_constant_time_initializer_1d 
        SLL_ALLOCATE(sim%A_time_func_params(1),ierr)
        sim%A_time_func_params(1) = 1._f64
      case ("SLL_COSSIN_FLOW")
        sim%phi_func => sll_f_cos_sin_initializer_2d
        SLL_ALLOCATE(sim%A_func_params(2),ierr)
        sim%A_func_params(1) = 0._f64
        sim%A_func_params(2) = 0._f64
        sim%A_time_func => sll_f_constant_time_initializer_1d 
        SLL_ALLOCATE(sim%A_time_func_params(1),ierr)
        sim%A_time_func_params(1) = 1._f64  
        
        sim%A1_func => sll_f_rotation_a1_initializer_2d 
        sim%A2_func => sll_f_rotation_a2_initializer_2d 
        sim%A1_exact_charac_func => sll_f_rotation_a1_exact_charac_2d 
        sim%A2_exact_charac_func => sll_f_rotation_a2_exact_charac_2d 
      case default
        print *,'#bad advect_case',advection_field_case
        print *,'#not implemented'
        print *,'#in initialize_analytic_field_2d_curvilinear'
        stop
    end select
    
   
    select case(initial_function_case)
      case ("SLL_KHP1")
        print*,"#f0 = SLL_KHP1" 
        sim%init_func => sll_f_khp1_2d
        SLL_ALLOCATE(sim%params(3),ierr)
        sim%params(1) = eps
        sim%params(2) = kmode_eta1
        sim%params(3) = kmode_eta2
      case ("SLL_GAUSSIAN")
        print*,"#f0 = SLL_GAUSSIAN " 
        sim%init_func => sll_f_gaussian_initializer_2d
        SLL_ALLOCATE(sim%params(4),ierr)
        sim%params(1) = xc_1
        sim%params(2) = xc_2
        sim%params(3) = sigma_1
        sim%params(4) = sigma_2
      case ("SLL_COS_BELL")
        print*,"#f0 = SLL__COS_BELL " 
        sim%init_func => sll_f_cos_bell_initializer_2d
        SLL_ALLOCATE(sim%params(2),ierr)
        sim%params(1) = xc_1
        sim%params(2) = xc_2  
      case ("SLL_COS_BELL0")
        print*,"#f0 = SLL__COS_BELL0 " 
        sim%init_func => sll_f_cos_bell0_initializer_2d
        SLL_ALLOCATE(sim%params(4),ierr)
        sim%params(1) = xc_1
        sim%params(2) = xc_2  
        sim%params(3) = cutx_value  
        sim%params(4) = cutf_value  
      case ("SLL_ONE")
        print*,"#f0 = SLL__ONE " 
        sim%init_func => sll_f_one_initializer_2d
        SLL_ALLOCATE(sim%params(1),ierr)
        sim%params(1) = eps
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
      case ("SLL_SPLITTING")
       print*,"#time_loop = SLL_SPLITTING " 
        sim%time_loop_case = SLL_SPLITTING
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
      case ("SLL_COMPUTE_FIELD_FROM_PHI_FD2")
        print*,"#compute_field_case = SLL_COMPUTE_FIELD_FROM_PHI_FD2" 
        sim%compute_field_case = SLL_COMPUTE_FIELD_FROM_PHI_FD2
      case ("SLL_COMPUTE_FIELD_FROM_PHI_FD3")
        print*,"#compute_field_case = SLL_COMPUTE_FIELD_FROM_PHI_FD3" 
        sim%compute_field_case = SLL_COMPUTE_FIELD_FROM_PHI_FD3
      case ("SLL_COMPUTE_FIELD_FROM_PHI_FD4")
        print*,"#compute_field_case = SLL_COMPUTE_FIELD_FROM_PHI_FD4" 
        sim%compute_field_case = SLL_COMPUTE_FIELD_FROM_PHI_FD4
      case ("SLL_COMPUTE_FIELD_FROM_PHI_FD5")
        print*,"#compute_field_case = SLL_COMPUTE_FIELD_FROM_PHI_FD5" 
        sim%compute_field_case = SLL_COMPUTE_FIELD_FROM_PHI_FD5
      case ("SLL_COMPUTE_FIELD_FROM_PHI_FD6")
        print*,"#compute_field_case = SLL_COMPUTE_FIELD_FROM_PHI_FD6" 
        sim%compute_field_case = SLL_COMPUTE_FIELD_FROM_PHI_FD6
      case ("SLL_COMPUTE_FIELD_FROM_PHI_FD")
        print*,"#compute_field_case = SLL_COMPUTE_FIELD_FROM_PHI_FD" 
        sim%compute_field_case = SLL_COMPUTE_FIELD_FROM_PHI_FD
      case default
        SLL_ERROR("initialize_analytic_field_2d_curvilinear","bad compute_field_case")
        stop
    end select

    
   
  end subroutine initialize_analytic_field_2d_curvilinear
  


  subroutine init_fake(sim, filename)
    class(sll_t_simulation_2d_analytic_field_curvilinear), intent(inout) :: sim
    character(len=*), intent(in)                                :: filename
  
    print *,'# Do not use the routine init_vp4d_fake'
    print *,'#use instead initialize_vlasov_par_poisson_seq_curv'
    stop
  
  end subroutine init_fake
  
  subroutine run_af2d_curvilinear(sim)
    class(sll_t_simulation_2d_analytic_field_curvilinear), intent(inout) :: sim
    sll_int32 :: Nc_eta1
    sll_int32 :: Nc_eta2
    sll_int32 :: num_dof1
    sll_int32 :: num_dof2
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
    sll_real64, dimension(:), allocatable :: dof_positions1
    sll_real64, dimension(:), allocatable :: dof_positions2
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
    !sll_real64 :: jac_m(2,2)

    Nc_eta1 = sim%mesh_2d%num_cells1
    Nc_eta2 = sim%mesh_2d%num_cells2
    num_dof1 = sim%num_dof1
    num_dof2 = sim%num_dof2
    delta_eta1 = sim%mesh_2d%delta_eta1
    delta_eta2 = sim%mesh_2d%delta_eta2
    eta1_min = sim%mesh_2d%eta1_min
    eta2_min = sim%mesh_2d%eta2_min
    eta1_max = sim%mesh_2d%eta1_max
    eta2_max = sim%mesh_2d%eta2_max
    nb_step = sim%num_iterations
    dt = sim%dt
    
    
    !allocation
    SLL_ALLOCATE(f(num_dof1,num_dof2),ierr)
    SLL_ALLOCATE(f_old(num_dof1,num_dof2),ierr)
    SLL_ALLOCATE(f_init(num_dof1,num_dof2),ierr)
    SLL_ALLOCATE(f_exact(num_dof1,num_dof2),ierr)
    SLL_ALLOCATE(phi(num_dof1,num_dof2),ierr)
    !SLL_ALLOCATE(phi(Nc_eta1+1,Nc_eta1+1),ierr)
    !SLL_ALLOCATE(A1(Nc_eta1+1,Nc_eta2+1),ierr)
    !SLL_ALLOCATE(A2(Nc_eta1+1,Nc_eta2+1),ierr)
    !SLL_ALLOCATE(A1_init(Nc_eta1+1,Nc_eta2+1),ierr)
    !SLL_ALLOCATE(A2_init(Nc_eta1+1,Nc_eta2+1),ierr)

    !SLL_ALLOCATE(A1(Nc_eta1+1,num_dof2),ierr)
    !SLL_ALLOCATE(A2(num_dof1,Nc_eta2+1),ierr)
    !SLL_ALLOCATE(A1_init(Nc_eta1+1,num_dof2),ierr)
    !SLL_ALLOCATE(A2_init(num_dof1,Nc_eta2+1),ierr)
    
    SLL_ALLOCATE(A1(num_dof1,num_dof2),ierr)
    SLL_ALLOCATE(A2(num_dof1,num_dof2),ierr)
    SLL_ALLOCATE(A1_init(num_dof1,num_dof2),ierr)
    SLL_ALLOCATE(A2_init(num_dof1,num_dof2),ierr)


    SLL_ALLOCATE(div(num_dof1,num_dof2),ierr)
    SLL_ALLOCATE(dof_positions1(num_dof1),ierr)
    SLL_ALLOCATE(dof_positions2(num_dof2),ierr)

    call compute_dof_positions(sim%eta1_array,dof_positions1)
    call compute_dof_positions(sim%eta2_array,dof_positions2)



    do i2=1,num_dof2
      eta2=dof_positions2(i2)
      do i1=1,num_dof1
        eta1=dof_positions1(i1)
        x1 = sim%transformation%x1(eta1,eta2)
        x2 = sim%transformation%x2(eta1,eta2)
        f(i1,i2) =  sim%init_func(x1,x2,sim%params) 
        f_init(i1,i2)  =  sim%init_func(x1,x2,sim%params)
        phi(i1,i2) = sim%phi_func(x1,x2,sim%A_func_params)
        !phi(i1,i2) = sim%phi_func(x1,x2,sim%A_func_params)
        !phi(i1,i2)= sin(x1)*cos(x2) 
        !-0.5_f64*(x1**2+x2**2) &
        !  +sim%A_func_params(1)*x2-sim%A_func_params(2)*x1
        !warning specific change for colella mesh
        !phi(i1,i2) = phi(i1,i2)-sim%A_func_params(1)*eta2+sim%A_func_params(2)*eta1
      end do
    end do

!    do i2=1,num_dof2
!      eta2=dof_positions2(i2)
!      do i1=1,num_dof1
!        eta1=dof_positions1(i1)
!    do i2=1,Nc_eta2+1
!      eta2=sim%eta2_array(i2)
!      do i1=1,Nc_eta1+1
!        eta1=sim%eta1_array(i1)
!        x1 = sim%transformation%x1(eta1,eta2)
!        x2 = sim%transformation%x2(eta1,eta2)
!        phi(i1,i2) = sim%phi_func(x1,x2,sim%A_func_params)
!      end do
!    end do
    
    select case (sim%compute_field_case)
      case (SLL_COMPUTE_FIELD_FROM_ANALYTIC)
        call compute_curvilinear_field_2d( &
          sim%A1_func, &
          sim%A2_func, &
          sim%A_func_params, &
          dof_positions1, &
          num_dof1, &
          dof_positions2, &
          num_dof2, &
          sim%transformation, &
          A1_init, &
          A2_init)  
      case (SLL_COMPUTE_FIELD_FROM_PHI)
        call compute_field_from_phi_2d_curvilinear( &
          phi, &
          sim%mesh_2d, &
          sim%transformation, &
          A1_init, &
          A2_init, &
          sim%phi_interp2d)
      case (SLL_COMPUTE_FIELD_FROM_PHI_FD)
            call compute_field_from_phi_2d_fd_curvilinear( &
              phi, &
              sim%mesh_2d, &
              sim%transformation, &
              A1_init, &
              A2_init, &
              sim%phi_interp2d, &
              sim%fd_degree1, &
              sim%fd_degree2)      
!        select case (sim%advection_form)
!          case (sll_p_advective)
!          case (sll_p_conservative)
!            call compute_field_from_phi_2d_fd_conservative_curvilinear2( &
!              phi, &
!              sim%mesh_2d, &
!              sim%transformation, &
!              A1_init, &
!              A2_init, &
!              sim%phi_interp2d, &
!              sim%fd_degree1, &
!              sim%fd_degree2)                
!        case default
!          print *,'#bad choice of advection_form:',sim%advection_form
!          stop
!        end select

        
      case default
        SLL_ERROR("run_af2d_curvilinear","bad value of sim%compute_field_case")
    end select    

    if(sim%time_loop_case==SLL_SPLITTING)then
      sim%split => sll_f_new_split_advection_2d( &
        f, &
        A1_init, &
        A2_init, &
        sim%interp1, &
        sim%charac1, &
        sim%process_outside_point1_func, &
        sim%interp2, &
        sim%charac2, &
        sim%process_outside_point2_func, &
        sim%mesh_2d, &
        sim%advection_form, &
        sll_p_strang_tvt, &
        transformation=sim%transformation, &
        csl_2012=sim%csl_2012) 

      A1 = A1_init
      A2 = A2_init 

!      sim%split => new_advection_2d( &
!        f, &
!        Nc_eta1+1, &
!        Nc_eta2+1, &
!        num_dof1, &
!        num_dof2, &
!        A1_init, &
!        A2_init, &
!        sim%advect1_1d, &
!        sim%advect2_1d, &
!        sim%advection_form, &
!        sll_p_strang_tvt)      
    endif



    

    call sll_s_ascii_file_create(sim%thdiag_filename, thdiag_id, ierr)
    
    iplot = 0

    err = 0._f64
    f_exact = f

 


    do step=1,nb_step
      !print*,"step= ", step
      f_old = f
    

      if(sim%advection_form==sll_p_conservative) then 
        call advective_to_conservative_2d_curvilinear( &
          f_old, &
          dof_positions1, &
          num_dof1, &
          dof_positions2, &
          num_dof2, &
          sim%transformation)
      endif



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
        case (SLL_SPLITTING)
          f = f_old
          !print *,'f_old=',minval(f),maxval(f)
           call sll_s_do_split_steps(sim%split, dt, 1)
          !print *,'f=',minval(f),maxval(f)
          !stop
        case default  
          print *,'#bad time_loop_case',sim%time_loop_case
          print *,'#not implemented'
          print *,'#in run_af2d_curvilinear'
          print *,'#available options are:'
          print *,'#SLL_EULER=',SLL_EULER
          print *,'#SLL_PREDICTOR_CORRECTOR=',SLL_PREDICTOR_CORRECTOR
          
      end select




      if(sim%advection_form==sll_p_conservative) then 
        call conservative_to_advective_2d_curvilinear( &
          f, &
          dof_positions1, &
          num_dof1, &
          dof_positions2, &
          num_dof2, &
          sim%transformation)
      endif

      
      do i2=1,num_dof2
        eta2=dof_positions2(i2)
        do i1=1,num_dof1
          eta1=dof_positions1(i1)
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
          !feet1 = sim%process_outside_point1_func( feet1, eta1_min, eta1_max )   
          !feet2 = sim%process_outside_point2_func( feet2, eta2_min, eta2_max )   
          f_exact(i1,i2) =  sim%init_func(feet1,feet2,sim%params) 
        end do
      end do
      
      err = max(maxval(abs(f_exact-f)),err)
      print *,"#",step,maxval(abs(f_exact-f))
      
     
          call compute_divergence_2d_curvilinear( &
            div, &
            sim%mesh_2d, &
            sim%transformation, &
            A1, &
            A2, &
            sim%phi_interp2d)


!      select case (sim%advection_form)
!        case (sll_p_advective)
!          call compute_divergence_2d_curvilinear( &
!            div, &
!            sim%mesh_2d, &
!            sim%transformation, &
!            A1, &
!            A2, &
!            sim%phi_interp2d)
!        case (sll_p_conservative)
!          call compute_divergence_2d_conservative_curvilinear( &
!            div, &
!            sim%mesh_2d, &
!            sim%transformation, &
!            A1, &
!            A2, &
!            sim%phi_interp2d)
!        case default
!          print *,'#bad choice of advection_form:',sim%advection_form
!          stop
!        end select
        
#ifndef NOHDF5
      if(modulo(step,sim%freq_diag)==0)then
        if(iplot==0)then
          call plot_f_curvilinear( &
            "divf", &
            iplot, &
            div, &
            dof_positions1, &
            dof_positions2, &
            num_dof1, &
            num_dof2, &
            sim%transformation)
        endif
        call plot_f_curvilinear( &
          "f", &
          iplot, &
          f, &
          dof_positions1, &
          dof_positions2, &
          num_dof1, &
          num_dof2, &
          sim%transformation)
        call plot_f_curvilinear( &
          "errf", &
          iplot, &
          f-f_exact, &
          dof_positions1, &
          dof_positions2, &
          num_dof1, &
          num_dof2, &
          sim%transformation)
        iplot = iplot+1  
      endif            
! stuff of Adnane: commented for the moment to fix conflict
!      
!      if(step==nb_step)then
!          call plot_visit_curvilinear(iplot-1,abs(f_exact-f),sim%mesh_2d,sim%transformation,"Error")
!          do i2=1,Nc_eta2+1
!            eta2=eta2_min+real(i2-1,f64)*delta_eta2
!            do i1=1,Nc_eta1+1
!              eta1=eta1_min+real(i1-1,f64)*delta_eta1
!              x1 = sim%transformation%x1(eta1,eta2)
!              x2 = sim%transformation%x2(eta1,eta2)
!              write(12,*) x1,x2,f(i1,i2) !,f_exact(i1,i2)
!            enddo
!          enddo    
!      endif        
#endif  

      if(modulo(step,sim%freq_diag_time)==0)then
        call time_history_diagnostic_curvilinear2( &
          thdiag_id , &    
          step, &
          dt, &
          sim%mesh_2d, &
          sim%transformation, &
          f, &
          dof_positions1, &
          dof_positions2, &
          num_dof1, &
          num_dof2, &
          sim%integration_weight1, &
          sim%integration_weight2, &    
          f_exact)

!        call time_history_diagnostic_curvilinear( &
!          thdiag_id , &    
!          step, &
!          dt, &
!          sim%mesh_2d, &
!          sim%transformation, &
!          f, &
!          phi, &
!          A1, &
!          A2, &
!          f_exact)
      endif            
      
      
         
    enddo
    
    close(thdiag_id)   
    !print *,err    
    print *,'#run_af2d_curvilinear PASSED'
    
  end subroutine run_af2d_curvilinear
  

  subroutine compute_dof_positions(input, output)
    sll_real64, dimension(:), intent(in) :: input
    sll_real64, dimension(:), intent(out) :: output
    sll_int32 :: Npts
    sll_int32 :: num_dof
    sll_int32 :: i
    
    Npts = size(input)
    num_dof = size(output)
    
    
    if(num_dof==Npts)then
      output(1:num_dof) = input(1:num_dof)
    else
      if(num_dof/=Npts-1)then
        SLL_ERROR('compute_node_positions','#bad value for Npts')
      endif
      do i=1,num_dof
        output(i) = 0.5_f64*(input(i)+input(i+1))
      enddo
    endif
  
  end subroutine compute_dof_positions


  subroutine compute_integration_weight(input,output)
    sll_real64, dimension(:), intent(in) :: input
    sll_real64, dimension(:), intent(out) :: output
    sll_int32 :: Npts
    sll_int32 :: num_dof
    sll_int32 :: i
    
    Npts = size(input)
    num_dof = size(output)
        
    if(num_dof==Npts)then
    !trapezoidal rule
      output(1)=0.5_f64*(input(2)-input(1))
      do i=2,Npts-1
        output(i) = 0.5_f64*(input(i+1)-input(i-1))
      enddo  
      output(Npts) = &
          0.5_f64*(input(Npts)-input(Npts-1))
    else
      if(num_dof/=Npts-1)then
        SLL_ERROR('compute_integration_weight','#bad value for Npts')
      endif
      !conservative rule
      do i=1,num_dof
        output(i) = input(i+1)-input(i)
      enddo
    endif  
  end subroutine compute_integration_weight

  

#ifndef NOHDF5
!*********************
!*********************

  !---------------------------------------------------
  ! Save the mesh structure
  !---------------------------------------------------
  subroutine plot_f_curvilinear( &
    filename, &
    iplot, &
    f, &
    node_positions1, &
    node_positions2, &
    n1, &
    n2, &
    transf)
    use sll_m_xdmf
    use sll_m_hdf5_io_serial
    character(len=*), intent(in) :: filename  !< file name
    sll_int32, intent(in) :: iplot
    sll_real64, dimension(:,:), intent(in) :: f
    sll_real64, dimension(:), intent(in) :: node_positions1
    sll_real64, dimension(:), intent(in) :: node_positions2
    sll_int32, intent(in) :: n1
    sll_int32, intent(in) :: n2
    class(sll_c_coordinate_transformation_2d_base), pointer :: transf    
    sll_int32 :: file_id
    sll_int32 :: error
    sll_real64, dimension(:,:), allocatable :: x1
    sll_real64, dimension(:,:), allocatable :: x2
    sll_int32 :: i, j
    character(len=4)      :: cplot
    !sll_int32             :: nnodes_x1, nnodes_x2
    !type(sll_t_cartesian_mesh_2d), pointer :: mesh_2d
    sll_real64 :: eta1
    sll_real64 :: eta2
    !sll_real64 ::  eta1_min, eta2_min
    !sll_real64 ::  eta1_max, eta2_max  
    !sll_real64 :: deta1
    !sll_real64 :: deta2
    
    
    
    !print *,'#maxf=',iplot,maxval(f),minval(f)
    

    
    if (iplot == 1) then

      SLL_ALLOCATE(x1(n1,n2), error)
      SLL_ALLOCATE(x2(n1,n2), error)
      do j = 1,n2
        do i = 1,n1
          eta1 = node_positions1(i) 
          eta2 = node_positions2(j)
          x1(i,j) = transf%x1(eta1,eta2)
          x2(i,j) = transf%x2(eta1,eta2)
        end do
      end do
      call sll_o_hdf5_file_create("curvilinear_mesh-x1.h5",file_id,error)
      call sll_o_hdf5_write_array(file_id,x1,"/x1",error)
      call sll_o_hdf5_file_close(file_id, error)
      call sll_o_hdf5_file_create("curvilinear_mesh-x2.h5",file_id,error)
      call sll_o_hdf5_write_array(file_id,x2,"/x2",error)
      call sll_o_hdf5_file_close(file_id, error)
      deallocate(x1)
      deallocate(x2)

    end if

    call sll_s_int2string(iplot,cplot)
    call sll_o_xdmf_open( &
      trim(filename)//cplot//".xmf", &
      "curvilinear_mesh", &
      n1, &
      n2, &
      file_id, &
      error)
    call sll_o_xdmf_write_array( &
      trim(filename)//cplot, &
      f, &
      "values", &
      error, &
      file_id, &
      "Node")
    call sll_s_xdmf_close(file_id,error)
  end subroutine plot_f_curvilinear

#endif






#ifndef NOHDF5
!*********************
!*********************

  subroutine plot_visit_curvilinear(iplot,f,mesh_2d,transf,file_name)
    use sll_m_xdmf
    use sll_m_hdf5_io_serial
    sll_int32 :: file_id
    sll_int32 :: error
    sll_real64, dimension(:,:), allocatable :: x1
    sll_real64, dimension(:,:), allocatable :: x2
    sll_int32 :: i, j
    sll_int32, intent(in)               :: iplot
     character(len=*),  intent(in)      :: file_name
    character(len=4)                    :: cplot
    sll_int32             :: nnodes_x1, nnodes_x2
    type(sll_t_cartesian_mesh_2d), pointer :: mesh_2d
    class(sll_c_coordinate_transformation_2d_base), pointer :: transf
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
      call sll_o_hdf5_file_create("curvilinear_mesh-x1.h5",file_id,error)
      call sll_o_hdf5_write_array(file_id,x1,"/x1",error)
      call sll_o_hdf5_file_close(file_id, error)
      call sll_o_hdf5_file_create("curvilinear_mesh-x2.h5",file_id,error)
      call sll_o_hdf5_write_array(file_id,x2,"/x2",error)
      call sll_o_hdf5_file_close(file_id, error)
      deallocate(x1)
      deallocate(x2)

    end if

    call sll_s_int2string(iplot,cplot)
    call sll_o_xdmf_open(file_name//cplot//".xmf","curvilinear_mesh", &
      nnodes_x1,nnodes_x2,file_id,error)
    call sll_o_xdmf_write_array(file_name//cplot,f,"values", &
      error,file_id,"Node")
    call sll_s_xdmf_close(file_id,error)
  end subroutine plot_visit_curvilinear

#endif



  subroutine compute_curvilinear_field_scalar_2d( &
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
    
  end subroutine compute_curvilinear_field_scalar_2d


  subroutine compute_curvilinear_field_2d( &
    A1_func, &
    A2_func, &
    A_func_params, &
    dof_positions1, &
    num_dof1, &
    dof_positions2, &
    num_dof2, &
    transformation, &
    A1, &
    A2)
    procedure(sll_i_scalar_initializer_2d), pointer :: A1_func
    procedure(sll_i_scalar_initializer_2d), pointer :: A2_func
    sll_real64, dimension(:), intent(in) :: A_func_params
    sll_real64, dimension(:), intent(in) :: dof_positions1
    sll_int32, intent(in) :: num_dof1
    sll_real64, dimension(:), intent(in) :: dof_positions2
    sll_int32, intent(in) :: num_dof2
    class(sll_c_coordinate_transformation_2d_base), pointer :: transformation
    sll_real64, dimension(:,:), intent(out) :: A1
    sll_real64, dimension(:,:), intent(out) :: A2
    sll_real64 :: input1
    sll_real64 :: input2
    sll_real64 :: output1
    sll_real64 :: output2
    sll_real64 :: jac_m(2,2)
    sll_real64 :: inv_jac
    sll_real64 :: inv_j11
    sll_real64 :: inv_j12
    sll_real64 :: inv_j21
    sll_real64 :: inv_j22
    sll_real64 :: eta1
    sll_real64 :: eta2
    sll_int32 :: i1
    sll_int32 :: i2
    sll_real64 :: x1
    sll_real64 :: x2


        do i2=1,num_dof2
          eta2=dof_positions2(i2)
          do i1=1,num_dof1
            eta1=dof_positions1(i1)
            x1 = transformation%x1(eta1,eta2)
            x2 = transformation%x2(eta1,eta2)
            jac_m  =  transformation%jacobian_matrix(eta1,eta2)          
            inv_jac = jac_m(1,1)*jac_m(2,2)-jac_m(1,2)*jac_m(2,1)  
            inv_jac = 1._f64/inv_jac
    
            inv_j11 =  jac_m(2,2)*inv_jac
            inv_j12 = -jac_m(1,2)*inv_jac
            inv_j21 = -jac_m(2,1)*inv_jac
            inv_j22 =  jac_m(1,1)*inv_jac
            input1 = A1_func(x1,x2,A_func_params)
            input2 = A2_func(x1,x2,A_func_params) 
            output1 = inv_j11*input1+inv_j12*input2 
            output2 = inv_j21*input1+inv_j22*input2 
            A1(i1,i2) = output1
            A2(i1,i2) = output2
          end do
        end do    
    
  end subroutine compute_curvilinear_field_2d




  subroutine compute_curvilinear_field_2d_old( &
    A1_func, &
    A2_func, &
    A_func_params, &
    eta1_array, &
    n1, &
    eta2_array, &
    n2, &
    dof_positions1, &
    num_dof1, &
    dof_positions2, &
    num_dof2, &
    transformation, &
    A1, &
    A2)
    procedure(sll_i_scalar_initializer_2d), pointer :: A1_func
    procedure(sll_i_scalar_initializer_2d), pointer :: A2_func
    sll_real64, dimension(:), intent(in) :: A_func_params
    sll_real64, dimension(:), intent(in) :: eta1_array
    sll_int32, intent(in) :: n1
    sll_real64, dimension(:), intent(in) :: eta2_array
    sll_int32, intent(in) :: n2
    sll_real64, dimension(:), intent(in) :: dof_positions1
    sll_int32, intent(in) :: num_dof1
    sll_real64, dimension(:), intent(in) :: dof_positions2
    sll_int32, intent(in) :: num_dof2
    class(sll_c_coordinate_transformation_2d_base), pointer :: transformation
    sll_real64, dimension(:,:), intent(out) :: A1
    sll_real64, dimension(:,:), intent(out) :: A2
    sll_real64 :: input1
    sll_real64 :: input2
    sll_real64 :: output1
    sll_real64 :: output2
    sll_real64 :: jac_m(2,2)
    sll_real64 :: inv_jac
    sll_real64 :: inv_j11
    sll_real64 :: inv_j12
    sll_real64 :: inv_j21
    sll_real64 :: inv_j22
    sll_real64 :: eta1
    sll_real64 :: eta2
    sll_int32 :: i1
    sll_int32 :: i2
    sll_real64 :: x1
    sll_real64 :: x2


        do i2=1,n2
          eta2=eta2_array(i2)
          do i1=1,num_dof1
            eta1=dof_positions1(i1)
            x1 = transformation%x1(eta1,eta2)
            x2 = transformation%x2(eta1,eta2)
            jac_m  =  transformation%jacobian_matrix(eta1,eta2)          
            inv_jac = jac_m(1,1)*jac_m(2,2)-jac_m(1,2)*jac_m(2,1)  
            inv_jac = 1._f64/inv_jac
    
            inv_j11 =  jac_m(2,2)*inv_jac
            inv_j12 = -jac_m(1,2)*inv_jac
            inv_j21 = -jac_m(2,1)*inv_jac
            inv_j22 =  jac_m(1,1)*inv_jac
            input1 = A1_func(x1,x2,A_func_params)
            input2 = A2_func(x1,x2,A_func_params) 
            output1 = inv_j11*input1+inv_j12*input2 
            output2 = inv_j21*input1+inv_j22*input2 
            !A1_init(i1,i2) = output1
            A2(i1,i2) = output2
          end do
        end do


        do i2=1,num_dof2
          eta2=dof_positions2(i2)
          do i1=1,n1
            eta1=eta1_array(i1)
            x1 = transformation%x1(eta1,eta2)
            x2 = transformation%x2(eta1,eta2)
            jac_m  =  transformation%jacobian_matrix(eta1,eta2)          
            inv_jac = jac_m(1,1)*jac_m(2,2)-jac_m(1,2)*jac_m(2,1)  
            inv_jac = 1._f64/inv_jac
    
            inv_j11 =  jac_m(2,2)*inv_jac
            inv_j12 = -jac_m(1,2)*inv_jac
            inv_j21 = -jac_m(2,1)*inv_jac
            inv_j22 =  jac_m(1,1)*inv_jac
            input1 = A1_func(x1,x2,A_func_params)
            input2 = A2_func(x1,x2,A_func_params) 
            output1 = inv_j11*input1+inv_j12*input2 
            output2 = inv_j21*input1+inv_j22*input2 
            A1(i1,i2) = output1
            !A2_init(i1,i2) = output2
          end do
        end do



    
    
  end subroutine compute_curvilinear_field_2d_old



  subroutine compute_field_from_phi_2d_curvilinear(phi,mesh_2d,transformation,A1,A2,interp2d)
    sll_real64, dimension(:,:), intent(in) :: phi
    sll_real64, dimension(:,:), intent(out) :: A1
    sll_real64, dimension(:,:), intent(out) :: A2
    type(sll_t_cartesian_mesh_2d), pointer :: mesh_2d
    class(sll_c_coordinate_transformation_2d_base), pointer :: transformation
    class(sll_c_interpolator_2d), pointer   :: interp2d
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
        A1(i1,i2)=interp2d%interpolate_from_interpolant_derivative_eta2(eta1,eta2)/transformation%jacobian(eta1,eta2)
        A2(i1,i2)=-interp2d%interpolate_from_interpolant_derivative_eta1(eta1,eta2)/transformation%jacobian(eta1,eta2)
      end do
    end do
   
    
  end subroutine compute_field_from_phi_2d_curvilinear


  subroutine compute_field_from_phi_2d_fd2_curvilinear(phi,mesh_2d,transformation,A1,A2,interp2d)
    sll_real64, dimension(:,:), intent(in) :: phi
    sll_real64, dimension(:,:), intent(out) :: A1
    sll_real64, dimension(:,:), intent(out) :: A2
    type(sll_t_cartesian_mesh_2d), pointer :: mesh_2d
    class(sll_c_coordinate_transformation_2d_base), pointer :: transformation
    class(sll_c_interpolator_2d), pointer   :: interp2d
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
        
      end do
    end do
   
    
  end subroutine compute_field_from_phi_2d_fd2_curvilinear


  subroutine compute_field_from_phi_2d_fd4_curvilinear(phi,mesh_2d,transformation,A1,A2,interp2d)
    sll_real64, dimension(:,:), intent(in) :: phi
    sll_real64, dimension(:,:), intent(out) :: A1
    sll_real64, dimension(:,:), intent(out) :: A2
    type(sll_t_cartesian_mesh_2d), pointer :: mesh_2d
    class(sll_c_coordinate_transformation_2d_base), pointer :: transformation
    class(sll_c_interpolator_2d), pointer   :: interp2d
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
        A1(i1,i2) = (2._f64/3._f64)*(phi(i1,modulo(i2+1-1+Nc_eta2,Nc_eta2)+1) &
          -phi(i1,modulo(i2-1-1+Nc_eta2,Nc_eta2)+1))
        A1(i1,i2) = A1(i1,i2)-(1._f64/12._f64)*(phi(i1,modulo(i2+2-1+Nc_eta2,Nc_eta2)+1) &
          -phi(i1,modulo(i2-2-1+Nc_eta2,Nc_eta2)+1))
        A1(i1,i2) = A1(i1,i2)/(delta_eta2)
        A1(i1,i2) = A1(i1,i2)/transformation%jacobian(eta1,eta2)
        A2(i1,i2) = (2._f64/3._f64)*(phi(modulo(i1+1-1+Nc_eta1,Nc_eta1)+1,i2) &
          -phi(modulo(i1-1-1+Nc_eta1,Nc_eta1)+1,i2))
        A2(i1,i2) = A2(i1,i2)-(1._f64/12._f64)*(phi(modulo(i1+2-1+Nc_eta1,Nc_eta1)+1,i2) &
          -phi(modulo(i1-2-1+Nc_eta1,Nc_eta1)+1,i2))
        A2(i1,i2) = A2(i1,i2)/(delta_eta1)
        A2(i1,i2) = -A2(i1,i2)/transformation%jacobian(eta1,eta2)
        
      end do
    end do
   
    
  end subroutine compute_field_from_phi_2d_fd4_curvilinear

  subroutine compute_field_from_phi_2d_fd6_curvilinear(phi,mesh_2d,transformation,A1,A2,interp2d)
    sll_real64, dimension(:,:), intent(in) :: phi
    sll_real64, dimension(:,:), intent(out) :: A1
    sll_real64, dimension(:,:), intent(out) :: A2
    type(sll_t_cartesian_mesh_2d), pointer :: mesh_2d
    class(sll_c_coordinate_transformation_2d_base), pointer :: transformation
    class(sll_c_interpolator_2d), pointer   :: interp2d
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
        A1(i1,i2) = (3._f64/4._f64)*(phi(i1,modulo(i2+1-1+Nc_eta2,Nc_eta2)+1) &
          -phi(i1,modulo(i2-1-1+Nc_eta2,Nc_eta2)+1))
        A1(i1,i2) = A1(i1,i2)-(3._f64/20._f64)*(phi(i1,modulo(i2+2-1+Nc_eta2,Nc_eta2)+1) &
          -phi(i1,modulo(i2-2-1+Nc_eta2,Nc_eta2)+1))
        A1(i1,i2) = A1(i1,i2)+(1._f64/60._f64)*(phi(i1,modulo(i2+3-1+Nc_eta2,Nc_eta2)+1) &
          -phi(i1,modulo(i2-3-1+Nc_eta2,Nc_eta2)+1))
        A1(i1,i2) = A1(i1,i2)/(delta_eta2)
        A1(i1,i2) = A1(i1,i2)/transformation%jacobian(eta1,eta2)
        A2(i1,i2) = (3._f64/4._f64)*(phi(modulo(i1+1-1+Nc_eta1,Nc_eta1)+1,i2) &
          -phi(modulo(i1-1-1+Nc_eta1,Nc_eta1)+1,i2))
        A2(i1,i2) = A2(i1,i2)-(3._f64/20._f64)*(phi(modulo(i1+2-1+Nc_eta1,Nc_eta1)+1,i2) &
          -phi(modulo(i1-2-1+Nc_eta1,Nc_eta1)+1,i2))
        A2(i1,i2) = A2(i1,i2)+(1._f64/60._f64)*(phi(modulo(i1+3-1+Nc_eta1,Nc_eta1)+1,i2) &
          -phi(modulo(i1-3-1+Nc_eta1,Nc_eta1)+1,i2))
        A2(i1,i2) = A2(i1,i2)/(delta_eta1)
        A2(i1,i2) = -A2(i1,i2)/transformation%jacobian(eta1,eta2)
        
      end do
    end do
   
    
  end subroutine compute_field_from_phi_2d_fd6_curvilinear





  subroutine compute_field_from_phi_2d_fd3_curvilinear(phi,mesh_2d,transformation,A1,A2,interp2d)
    sll_real64, dimension(:,:), intent(in) :: phi
    sll_real64, dimension(:,:), intent(out) :: A1
    sll_real64, dimension(:,:), intent(out) :: A2
    type(sll_t_cartesian_mesh_2d), pointer :: mesh_2d
    class(sll_c_coordinate_transformation_2d_base), pointer :: transformation
    class(sll_c_interpolator_2d), pointer   :: interp2d
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
        A1(i1,i2) = (-1._f64/3._f64)*phi(i1,modulo(i2-1-1+Nc_eta2,Nc_eta2)+1)
        A1(i1,i2) = A1(i1,i2)+(-1._f64/2._f64)*phi(i1,modulo(i2+0-1+Nc_eta2,Nc_eta2)+1)
        A1(i1,i2) = A1(i1,i2)+(1._f64)*phi(i1,modulo(i2+1-1+Nc_eta2,Nc_eta2)+1)
        A1(i1,i2) = A1(i1,i2)+(-1._f64/6._f64)*phi(i1,modulo(i2+2-1+Nc_eta2,Nc_eta2)+1)
        A1(i1,i2) = A1(i1,i2)/(delta_eta2)
        A1(i1,i2) = A1(i1,i2)/transformation%jacobian(eta1,eta2)
        A2(i1,i2) = (-1._f64/3._f64)*phi(modulo(i1-1-1+Nc_eta1,Nc_eta1)+1,i2)
        A2(i1,i2) = A2(i1,i2)+(-1._f64/2._f64)*phi(modulo(i1+0-1+Nc_eta1,Nc_eta1)+1,i2)
        A2(i1,i2) = A2(i1,i2)+(1._f64)*phi(modulo(i1+1-1+Nc_eta1,Nc_eta1)+1,i2)
        A2(i1,i2) = A2(i1,i2)+(-1._f64/6._f64)*phi(modulo(i1+2-1+Nc_eta1,Nc_eta1)+1,i2)

        A2(i1,i2) = A2(i1,i2)/(delta_eta1)
        A2(i1,i2) = -A2(i1,i2)/transformation%jacobian(eta1,eta2)
        
      end do
    end do
   
    
  end subroutine compute_field_from_phi_2d_fd3_curvilinear


subroutine compute_field_from_phi_2d_fd5_curvilinear(phi,mesh_2d,transformation,A1,A2,interp2d)
    sll_real64, dimension(:,:), intent(in) :: phi
    sll_real64, dimension(:,:), intent(out) :: A1
    sll_real64, dimension(:,:), intent(out) :: A2
    type(sll_t_cartesian_mesh_2d), pointer :: mesh_2d
    class(sll_c_coordinate_transformation_2d_base), pointer :: transformation
    class(sll_c_interpolator_2d), pointer   :: interp2d
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
        A1(i1,i2) = (1._f64/20._f64)*phi(i1,modulo(i2-2-1+Nc_eta2,Nc_eta2)+1)
        A1(i1,i2) = A1(i1,i2)+(-1._f64/2._f64)*phi(i1,modulo(i2-1-1+Nc_eta2,Nc_eta2)+1)
        A1(i1,i2) = A1(i1,i2)+(-1._f64/3._f64)*phi(i1,modulo(i2+0-1+Nc_eta2,Nc_eta2)+1)
        A1(i1,i2) = A1(i1,i2)+(1._f64)*phi(i1,modulo(i2+1-1+Nc_eta2,Nc_eta2)+1)
        A1(i1,i2) = A1(i1,i2)+(-1._f64/4._f64)*phi(i1,modulo(i2+2-1+Nc_eta2,Nc_eta2)+1)
        A1(i1,i2) = A1(i1,i2)+(1._f64/30._f64)*phi(i1,modulo(i2+3-1+Nc_eta2,Nc_eta2)+1)
        A1(i1,i2) = A1(i1,i2)/(delta_eta2)
        A1(i1,i2) = A1(i1,i2)/transformation%jacobian(eta1,eta2)
        A2(i1,i2) = (1._f64/20._f64)*phi(modulo(i1-2-1+Nc_eta1,Nc_eta1)+1,i2)
        A2(i1,i2) = A2(i1,i2)+(-1._f64/2._f64)*phi(modulo(i1-1-1+Nc_eta1,Nc_eta1)+1,i2)
        A2(i1,i2) = A2(i1,i2)+(-1._f64/3._f64)*phi(modulo(i1+0-1+Nc_eta1,Nc_eta1)+1,i2)
        A2(i1,i2) = A2(i1,i2)+(1._f64)*phi(modulo(i1+1-1+Nc_eta1,Nc_eta1)+1,i2)
        A2(i1,i2) = A2(i1,i2)+(-1._f64/4._f64)*phi(modulo(i1+2-1+Nc_eta1,Nc_eta1)+1,i2)
        A2(i1,i2) = A2(i1,i2)+(1._f64/30._f64)*phi(modulo(i1+3-1+Nc_eta1,Nc_eta1)+1,i2)

        A2(i1,i2) = A2(i1,i2)/(delta_eta1)
        A2(i1,i2) = -A2(i1,i2)/transformation%jacobian(eta1,eta2)
        
      end do
    end do
   
    
  end subroutine compute_field_from_phi_2d_fd5_curvilinear

subroutine compute_field_from_phi_2d_fd_curvilinear(phi,mesh_2d,transformation,A1,A2,interp2d,d1,d2)
    sll_real64, dimension(:,:), intent(in) :: phi
    sll_real64, dimension(:,:), intent(out) :: A1
    sll_real64, dimension(:,:), intent(out) :: A2
    type(sll_t_cartesian_mesh_2d), pointer :: mesh_2d
    class(sll_c_coordinate_transformation_2d_base), pointer :: transformation
    class(sll_c_interpolator_2d), pointer   :: interp2d
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
     
    call sll_s_compute_w_hermite(w1,r1,s1)
    call sll_s_compute_w_hermite(w2,r2,s2)

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




subroutine compute_field_from_phi_2d_fd_conservative_curvilinear(phi,mesh_2d,transformation,A1,A2,interp2d,d1,d2)
    sll_real64, dimension(:,:), intent(in) :: phi
    sll_real64, dimension(:,:), intent(out) :: A1
    sll_real64, dimension(:,:), intent(out) :: A2
    type(sll_t_cartesian_mesh_2d), pointer :: mesh_2d
    class(sll_c_coordinate_transformation_2d_base), pointer :: transformation
    class(sll_c_interpolator_2d), pointer   :: interp2d
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
    sll_real64, dimension(:), allocatable :: tmp1
    sll_real64, dimension(:), allocatable :: tmp2
    sll_real64, dimension(:), allocatable :: jac1
    sll_real64, dimension(:), allocatable :: jac2
    sll_real64, dimension(:,:), allocatable :: dphi
    sll_int32 :: ierr
    
    Nc_eta1 = mesh_2d%num_cells1
    Nc_eta2 = mesh_2d%num_cells2
    eta1_min = mesh_2d%eta1_min
    eta2_min = mesh_2d%eta2_min
    delta_eta1 = mesh_2d%delta_eta1
    delta_eta2 = mesh_2d%delta_eta2
    
    SLL_ALLOCATE(tmp1(Nc_eta1+1),ierr)
    SLL_ALLOCATE(tmp2(Nc_eta2+1),ierr)
    SLL_ALLOCATE(jac1(Nc_eta1),ierr)
    SLL_ALLOCATE(jac2(Nc_eta2),ierr)
    SLL_ALLOCATE(dphi(Nc_eta1+1,Nc_eta2+1),ierr)
    
    do i2=1,Nc_eta2
      call compute_deriv_csl_periodic( &
        phi(1:Nc_eta1,i2), &
        Nc_eta1, &
        tmp1(1:Nc_eta1+1), &
        d1)
      dphi(1:Nc_eta1+1,i2) = tmp1(1:Nc_eta1+1)   
    enddo

    do i1=1,Nc_eta1+1
      call compute_deriv_csl_periodic( &
        dphi(i1,1:Nc_eta2), &
        Nc_eta2, &
        tmp2(1:Nc_eta2+1), &
        d2)
      dphi(i1,1:Nc_eta2+1) = tmp2(1:Nc_eta2+1)   
    enddo

    do i2=1,Nc_eta2
      eta2=eta2_min+(real(i2-1,f64)+0.5_f64)*delta_eta2
      do i1=1,Nc_eta1
        eta1=eta1_min+(real(i1-1,f64)+0.5_f64)*delta_eta1
        jac1(i1) = transformation%jacobian(eta1,eta2)
      enddo
      call compute_deriv_csl_periodic( &
        jac1(1:Nc_eta1), &
        Nc_eta1, &
        tmp1(1:Nc_eta1+1), &
        d1)
      do i1=1,Nc_eta1+1
        A1(i1,i2) = (dphi(i1,i2+1)-dphi(i1,i2))/delta_eta2
        !A1(i1,i2) = A1(i1,i2)/tmp1(i1)
        eta1=eta1_min+(real(i1-1,f64))*delta_eta1
        eta2=eta2_min+(real(i2-1,f64)+0.5_f64)*delta_eta2
        A1(i1,i2) = A1(i1,i2)/transformation%jacobian(eta1,eta2)
      enddo
    enddo


    do i1=1,Nc_eta1
      eta1=eta1_min+(real(i1-1,f64)+0.5_f64)*delta_eta1
      do i2=1,Nc_eta2
        eta2=eta2_min+(real(i2-1,f64)+0.5_f64)*delta_eta2
        jac2(i2) = transformation%jacobian(eta1,eta2)
      enddo
      call compute_deriv_csl_periodic( &
        jac2(1:Nc_eta2), &
        Nc_eta2, &
        tmp2(1:Nc_eta2+1), &
        d2)
      do i2=1,Nc_eta2+1
        A2(i1,i2) = (dphi(i1+1,i2)-dphi(i1,i2))/delta_eta1
        !A2(i1,i2) = -A2(i1,i2)/tmp2(i2)
        eta1=eta1_min+(real(i1-1,f64)+0.5_f64)*delta_eta1
        eta2=eta2_min+(real(i2-1,f64))*delta_eta2
        A2(i1,i2) = -A2(i1,i2)/transformation%jacobian(eta1,eta2)
      enddo
    enddo
    
    
    !print *,'#dphi',maxval(dphi),minval(dphi)
    !print *,'#A1',maxval(A1),minval(A1)
    !print *,'#A2',maxval(A2),minval(A2)
    !stop   
  end subroutine compute_field_from_phi_2d_fd_conservative_curvilinear

subroutine compute_field_from_phi_2d_fd_conservative_curvilinear2(phi,mesh_2d,transformation,A1,A2,interp2d,d1,d2)
    sll_real64, dimension(:,:), intent(in) :: phi
    sll_real64, dimension(:,:), intent(out) :: A1
    sll_real64, dimension(:,:), intent(out) :: A2
    type(sll_t_cartesian_mesh_2d), pointer :: mesh_2d
    class(sll_c_coordinate_transformation_2d_base), pointer :: transformation
    class(sll_c_interpolator_2d), pointer   :: interp2d
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
    sll_real64, dimension(:), allocatable :: tmp1
    sll_real64, dimension(:), allocatable :: tmp2
    sll_real64, dimension(:), allocatable :: jac1
    sll_real64, dimension(:), allocatable :: jac2
    !sll_real64, dimension(:,:), allocatable :: dphi
    sll_int32 :: ierr
    
    Nc_eta1 = mesh_2d%num_cells1
    Nc_eta2 = mesh_2d%num_cells2
    eta1_min = mesh_2d%eta1_min
    eta2_min = mesh_2d%eta2_min
    delta_eta1 = mesh_2d%delta_eta1
    delta_eta2 = mesh_2d%delta_eta2
    
    SLL_ALLOCATE(tmp1(Nc_eta1+1),ierr)
    SLL_ALLOCATE(tmp2(Nc_eta2+1),ierr)
    SLL_ALLOCATE(jac1(Nc_eta1),ierr)
    SLL_ALLOCATE(jac2(Nc_eta2),ierr)
    !SLL_ALLOCATE(dphi(Nc_eta1+1,Nc_eta2+1),ierr)
    

    do i2=1,Nc_eta2
      eta2=eta2_min+(real(i2-1,f64)+0.5_f64)*delta_eta2
      do i1=1,Nc_eta1
        eta1=eta1_min+(real(i1-1,f64)+0.5_f64)*delta_eta1
        jac1(i1) = transformation%jacobian(eta1,eta2)
      enddo
      call compute_deriv_csl_periodic( &
        jac1(1:Nc_eta1), &
        Nc_eta1, &
        tmp1(1:Nc_eta1+1), &
        d1)
      do i1=1,Nc_eta1+1
        A1(i1,i2) = (phi(i1,i2+1)-phi(i1,i2))/delta_eta2
        A1(i1,i2) = A1(i1,i2)/tmp1(i1)
        eta1=eta1_min+(real(i1-1,f64))*delta_eta1
        eta2=eta2_min+(real(i2-1,f64)+0.5_f64)*delta_eta2
        !A1(i1,i2) = A1(i1,i2)/transformation%jacobian(eta1,eta2)
      enddo
    enddo


    do i1=1,Nc_eta1
      eta1=eta1_min+(real(i1-1,f64)+0.5_f64)*delta_eta1
      do i2=1,Nc_eta2
        eta2=eta2_min+(real(i2-1,f64)+0.5_f64)*delta_eta2
        jac2(i2) = transformation%jacobian(eta1,eta2)
      enddo
      call compute_deriv_csl_periodic( &
        jac2(1:Nc_eta2), &
        Nc_eta2, &
        tmp2(1:Nc_eta2+1), &
        d2)
      do i2=1,Nc_eta2+1
        A2(i1,i2) = (phi(i1+1,i2)-phi(i1,i2))/delta_eta1
        A2(i1,i2) = -A2(i1,i2)/tmp2(i2)
        eta1=eta1_min+(real(i1-1,f64)+0.5_f64)*delta_eta1
        eta2=eta2_min+(real(i2-1,f64))*delta_eta2
        !A2(i1,i2) = -A2(i1,i2)/transformation%jacobian(eta1,eta2)
      enddo
    enddo
    
    
    !print *,'#dphi',maxval(dphi),minval(dphi)
    !print *,'#A1',maxval(A1),minval(A1)
    !print *,'#A2',maxval(A2),minval(A2)
    !stop   
  end subroutine compute_field_from_phi_2d_fd_conservative_curvilinear2






  subroutine advective_to_conservative_2d_curvilinear( &
    f, &
    dof_positions1, &
    num_dof1, &
    dof_positions2, &
    num_dof2, &
    transformation)
    sll_real64, dimension(:,:), intent(inout) :: f
    sll_real64, dimension(:), intent(in) :: dof_positions1
    sll_int32, intent(in) :: num_dof1
    sll_real64, dimension(:), intent(in) :: dof_positions2
    sll_int32, intent(in) :: num_dof2
    class(sll_c_coordinate_transformation_2d_base), pointer :: transformation
    sll_real64 :: eta1
    sll_real64 :: eta2
    sll_int32 :: i1
    sll_int32 :: i2
    !sll_int32 :: ierr
    

    do i2=1,num_dof2
      eta2=dof_positions2(i2)
      do i1=1,num_dof1
        eta1=dof_positions1(i1)
        f(i1,i2)=f(i1,i2)*transformation%jacobian(eta1,eta2)
      end do
    end do
    
  end subroutine advective_to_conservative_2d_curvilinear






  subroutine conservative_to_advective_2d_curvilinear( &
    f, &
    dof_positions1, &
    num_dof1, &
    dof_positions2, &
    num_dof2, &
    transformation)
    sll_real64, dimension(:,:), intent(inout) :: f
    sll_real64, dimension(:), intent(in) :: dof_positions1
    sll_int32, intent(in) :: num_dof1
    sll_real64, dimension(:), intent(in) :: dof_positions2
    sll_int32, intent(in) :: num_dof2
    class(sll_c_coordinate_transformation_2d_base), pointer :: transformation
    sll_real64 :: eta1
    sll_real64 :: eta2
    sll_int32 :: i1
    sll_int32 :: i2
    !sll_int32 :: ierr
    

    do i2=1,num_dof2
      eta2=dof_positions2(i2)
      do i1=1,num_dof1
        eta1=dof_positions1(i1)
        f(i1,i2)=f(i1,i2)/transformation%jacobian(eta1,eta2)
      end do
    end do
    
  end subroutine conservative_to_advective_2d_curvilinear







  subroutine compute_divergence_2d_curvilinear(div,mesh_2d,transformation,A1,A2,interp2d)
    sll_real64, dimension(:,:), intent(out) :: div
    sll_real64, dimension(:,:), intent(in) :: A1
    sll_real64, dimension(:,:), intent(in) :: A2
    sll_real64, dimension(:,:), allocatable :: tmp
    type(sll_t_cartesian_mesh_2d), pointer :: mesh_2d
    class(sll_c_coordinate_transformation_2d_base), pointer :: transformation
    class(sll_c_interpolator_2d), pointer   :: interp2d
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
        div(i1,i2)=interp2d%interpolate_from_interpolant_derivative_eta1(eta1,eta2)
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
        div(i1,i2)=div(i1,i2)+interp2d%interpolate_from_interpolant_derivative_eta2(eta1,eta2)
      end do
    end do

   
    
  end subroutine compute_divergence_2d_curvilinear


  subroutine compute_divergence_2d_conservative_curvilinear(div,mesh_2d,transformation,A1,A2,interp2d)
    sll_real64, dimension(:,:), intent(out) :: div
    sll_real64, dimension(:,:), intent(in) :: A1
    sll_real64, dimension(:,:), intent(in) :: A2
    sll_real64, dimension(:,:), allocatable :: tmp
    type(sll_t_cartesian_mesh_2d), pointer :: mesh_2d
    class(sll_c_coordinate_transformation_2d_base), pointer :: transformation
    class(sll_c_interpolator_2d), pointer   :: interp2d
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
    sll_int32 :: ip
    sll_int32 :: im
    sll_real64 :: tmp2
    
    Nc_eta1 = mesh_2d%num_cells1
    Nc_eta2 = mesh_2d%num_cells2
    eta1_min = mesh_2d%eta1_min
    eta2_min = mesh_2d%eta2_min
    delta_eta1 = mesh_2d%delta_eta1
    delta_eta2 = mesh_2d%delta_eta2
    
    SLL_ALLOCATE(tmp(Nc_eta1+1,Nc_eta2+1),ierr)

    do i2=1,Nc_eta2
      eta2=eta2_min+(real(i2-1,f64)+0.5_f64)*delta_eta2
      do i1=1,Nc_eta1+1
        eta1=eta1_min+real(i1-1,f64)*delta_eta1
        tmp(i1,i2)=A1(i1,i2)*transformation%jacobian(eta1,eta2)
      end do
    end do
    do i1=1,Nc_eta1
      do i2=1,Nc_eta2
        div(i1,i2) = tmp(i1+1,i2)-tmp(i1,i2)
        div(i1,i2) = div(i1,i2)/delta_eta1
      enddo
    enddo     

    do i2=1,Nc_eta2+1
      eta2=eta2_min+real(i2-1,f64)*delta_eta2
      do i1=1,Nc_eta1
        eta1=eta1_min+(real(i1-1,f64)+0.5_f64)*delta_eta1
        tmp(i1,i2)=A2(i1,i2)*transformation%jacobian(eta1,eta2)
      end do
    end do


    do i1=1,Nc_eta1
      do i2=1,Nc_eta2
        ip = modulo(i2+1+Nc_eta2-1,Nc_eta2)+1
        im = modulo(i2-1+Nc_eta2-1,Nc_eta2)+1
        tmp2 = tmp(i1,i2+1)-tmp(i1,i2)
        tmp2 = tmp2/delta_eta2
        div(i1,i2) = div(i1,i2)+tmp2
      enddo
    enddo     


   
    
  end subroutine compute_divergence_2d_conservative_curvilinear



  subroutine time_history_diagnostic_curvilinear2( &
    file_id, &    
    step, &
    dt, &
    mesh_2d, &
    transformation,&
    f, &
    dof_positions1, &
    dof_positions2, &
    num_dof1, &
    num_dof2, &
    integration_weight1, &
    integration_weight2, &    
    !phi, &
    !A1, &
    !A2, &
    f_exact)
    sll_int32, intent(in) :: file_id
    sll_int32, intent(in) :: step
    sll_real64, intent(in) :: dt
    type(sll_t_cartesian_mesh_2d), pointer :: mesh_2d
    class(sll_c_coordinate_transformation_2d_base), pointer :: transformation
    sll_real64, dimension(:,:), intent(in) :: f
    sll_real64, dimension(:), intent(in) :: dof_positions1
    sll_real64, dimension(:), intent(in) :: dof_positions2
    sll_real64, dimension(:), intent(in) :: integration_weight1
    sll_real64, dimension(:), intent(in) :: integration_weight2
    sll_int32, intent(in) :: num_dof1
    sll_int32, intent(in) :: num_dof2
    !sll_real64, dimension(:,:), intent(in) :: phi
    !sll_real64, dimension(:,:), intent(in) :: A1
    !sll_real64, dimension(:,:), intent(in) :: A2 
    sll_real64, dimension(:,:), intent(in) :: f_exact
    sll_real64 :: mass
    sll_real64 :: linf
    sll_real64 :: l1
    sll_real64 :: l2
    !sll_real64 :: e
    
    sll_real64, dimension(:), allocatable :: mass_array
    sll_real64, dimension(:), allocatable  :: l1_array
    sll_real64, dimension(:), allocatable  :: l2_array
    !sll_real64, dimension(:), allocatable  :: e_array
    
    sll_real64 :: eta1
    sll_real64 :: eta2
    sll_real64, dimension(:),allocatable :: data
    !sll_real64, dimension(1:2,1:2) :: jac_m
    sll_int32 :: i1
    sll_int32 :: i2
!    sll_int32 :: Nc_eta1
!    sll_int32 :: Nc_eta2
!    sll_real64 :: eta1_min
!    sll_real64 :: eta1_max
!    sll_real64 :: eta2_min
!    sll_real64 :: eta2_max
!    sll_real64 :: delta_eta1
!    sll_real64 :: delta_eta2
!    sll_real64 :: dphi_eta1
!    sll_real64 :: dphi_eta2
    sll_int32 :: ierr
    sll_real64 :: err_linf 
    sll_real64 :: err_l1 
    sll_real64 :: err_l2 

    
!    Nc_eta1 = mesh_2d%num_cells1
!    Nc_eta2 = mesh_2d%num_cells2
    
    
!    eta1_min = mesh_2d%eta1_min
!    eta1_max = mesh_2d%eta1_max
!    eta2_min = mesh_2d%eta2_min
!    eta2_max = mesh_2d%eta2_max
!    delta_eta1 = mesh_2d%delta_eta1
!    delta_eta2 = mesh_2d%delta_eta2


!    SLL_ALLOCATE(data(Nc_eta1+1),ierr)
!    SLL_ALLOCATE(mass_array(Nc_eta2+1),ierr)
!    SLL_ALLOCATE(l1_array(Nc_eta2+1),ierr)
!    SLL_ALLOCATE(l2_array(Nc_eta2+1),ierr)
!    SLL_ALLOCATE(e_array(Nc_eta2+1),ierr)


    SLL_ALLOCATE(data(num_dof2),ierr)
    SLL_ALLOCATE(mass_array(num_dof2),ierr)
    SLL_ALLOCATE(l1_array(num_dof2),ierr)
    SLL_ALLOCATE(l2_array(num_dof2),ierr)
    linf  = 0.0_f64
    !l1    = 0.0_f64
    !l2    = 0.0_f64
    !mass  = 0.0_f64
     !e     = 0.0_f64
    
    do i2 = 1, num_dof2
      eta2 = dof_positions2(i2) 
      do i1=1, num_dof1
        eta1 = dof_positions1(i1)
        data(i1) = f(i1,i2)*abs(transformation%jacobian(eta1,eta2))
      enddo
      mass_array(i2) = sum(data(1:num_dof1)*integration_weight1(1:num_dof1))
      !sll_f_compute_integral_trapezoid_1d(data, Nc_eta1+1, delta_eta1)

      do i1=1,num_dof1
        eta1 = dof_positions1(i1)
        data(i1) = abs(f(i1,i2))*abs(transformation%jacobian(eta1,eta2))
      enddo
      l1_array(i2) = sum(data(1:num_dof1)*integration_weight1(1:num_dof1))
      !sll_f_compute_integral_trapezoid_1d(data, Nc_eta1+1, delta_eta1)

      do i1=1,num_dof1
        eta1 = dof_positions1(i1)
        data(i1) = (f(i1,i2))**2 *abs(transformation%jacobian(eta1,eta2))
      enddo
      l2_array(i2) = sum(data(1:num_dof1)*integration_weight1(1:num_dof1))
      !sll_f_compute_integral_trapezoid_1d(data, Nc_eta1+1, delta_eta1)


      do i1=1,num_dof1
       linf = max(linf,abs(f(i1,i2)))
      enddo
         
    enddo     

!    do i2 = 1,Nc_eta2+1
!      do i1 = 1,Nc_eta1+1
!        eta1 = eta1_min + (i1-1)* delta_eta1
!        jac_m  =  transformation%jacobian_matrix(eta1,eta2)
!        dphi_eta1 = -A2(i1,i2)* transformation%jacobian(eta1,eta2)
!        dphi_eta2 = A1(i1,i2)* transformation%jacobian(eta1,eta2)
!        data(i1) = (( jac_m(2,2)*dphi_eta1 - jac_m(2,1)*dphi_eta2 )**2 + &
!        ( -jac_m(1,2)*dphi_eta1 + jac_m(1,1)*dphi_eta2 )**2) &
!        /abs(transformation%jacobian(eta1,eta2)) 
!      enddo
!      e_array(i2) = sll_f_compute_integral_trapezoid_1d(data, Nc_eta1+1, delta_eta1)
!    enddo
!    e = sll_f_compute_integral_trapezoid_1d(e_array, Nc_eta2+1, delta_eta2)

    mass = sum(mass_array(1:num_dof2)*integration_weight2(1:num_dof2))
    l1 = sum(l1_array(1:num_dof2)*integration_weight2(1:num_dof2))
    l2 = sum(l2_array(1:num_dof2)*integration_weight2(1:num_dof2))
    !sll_f_compute_integral_trapezoid_1d(mass_array, Nc_eta2+1, delta_eta2)
    !l1 = sll_f_compute_integral_trapezoid_1d(l1_array, Nc_eta2+1, delta_eta2)
    !l2 = sll_f_compute_integral_trapezoid_1d(l2_array, Nc_eta2+1, delta_eta2)
    l2 = sqrt(l2)


    !mass = mass*delta_eta2
    !l1 = l1*delta_eta2
    !l2 = sqrt(l2*delta_eta2)
    !e  = e*delta_eta2
    
    
    !now, compute errors


    err_linf  = 0.0_f64
    
    do i2 = 1, num_dof2
      eta2 = dof_positions2(i2) 

      do i1=1,num_dof1
        eta1 = dof_positions1(i1) 
        data(i1) = abs(f(i1,i2)-f_exact(i1,i2))*abs(transformation%jacobian(eta1,eta2))
      enddo
      l1_array(i2) = sum(data(1:num_dof1)*integration_weight1(1:num_dof1))
      !sll_f_compute_integral_trapezoid_1d(data, Nc_eta1+1, delta_eta1)

      do i1=1,num_dof1
        eta1 = dof_positions1(i1)
        data(i1) = (f(i1,i2)-f_exact(i1,i2))**2 *abs(transformation%jacobian(eta1,eta2))
      enddo
      l2_array(i2) = sum(data(1:num_dof1)*integration_weight1(1:num_dof1))
      !sll_f_compute_integral_trapezoid_1d(data, Nc_eta1+1, delta_eta1)


      do i1=1,num_dof1
       err_linf = max(err_linf,abs(f(i1,i2)-f_exact(i1,i2)))
      enddo
         
    enddo     

    err_l1 = sum(l1_array(1:num_dof2)*integration_weight2(1:num_dof2))
    !sll_f_compute_integral_trapezoid_1d(l1_array, Nc_eta2+1, delta_eta2)
    err_l2 = sum(l2_array(1:num_dof2)*integration_weight2(1:num_dof2))
    !sll_f_compute_integral_trapezoid_1d(l2_array, Nc_eta2+1, delta_eta2)
    err_l2 = sqrt(err_l2)
    
    write(file_id,*) &
      dt*real(step,f64), &
      linf, &
      l1, &
      l2, &
      mass, &
      0._f64, &
      0._f64, &
!      e, &
!      maxval(abs(phi(1:Nc_eta1+1,1:Nc_eta2+1))), &
      err_linf, &
      err_l1, &
      err_l2
      
   
    
  end subroutine time_history_diagnostic_curvilinear2




  subroutine time_history_diagnostic_curvilinear( &
    file_id, &    
    step, &
    dt, &
    mesh_2d, &
    transformation,&
    f, &
    phi, &
    A1, &
    A2, &
    f_exact)
    sll_int32, intent(in) :: file_id
    sll_int32, intent(in) :: step
    sll_real64, intent(in) :: dt
    type(sll_t_cartesian_mesh_2d), pointer :: mesh_2d
    class(sll_c_coordinate_transformation_2d_base), pointer :: transformation
    sll_real64, dimension(:,:), intent(in) :: f
    sll_real64, dimension(:,:), intent(in) :: phi
    sll_real64, dimension(:,:), intent(in) :: A1
    sll_real64, dimension(:,:), intent(in) :: A2 
    sll_real64, dimension(:,:), intent(in) :: f_exact
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
    sll_real64 :: err_linf 
    sll_real64 :: err_l1 
    sll_real64 :: err_l2 

    
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
      mass_array(i2) = sll_f_compute_integral_trapezoid_1d(data, Nc_eta1+1, delta_eta1)

      do i1=1,Nc_eta1+1
        eta1 = eta1_min + (i1-1)* delta_eta1
        data(i1) = abs(f(i1,i2))*abs(transformation%jacobian(eta1,eta2))
      enddo
      l1_array(i2) = sll_f_compute_integral_trapezoid_1d(data, Nc_eta1+1, delta_eta1)

      do i1=1,Nc_eta1+1
        eta1 = eta1_min + (i1-1)* delta_eta1
        data(i1) = (f(i1,i2))**2 *abs(transformation%jacobian(eta1,eta2))
      enddo
      l2_array(i2) = sll_f_compute_integral_trapezoid_1d(data, Nc_eta1+1, delta_eta1)

      do i1=1,Nc_eta1+1
        eta1 = eta1_min + (i1-1)* delta_eta1
        jac_m  =  transformation%jacobian_matrix(eta1,eta2)
        dphi_eta1 = -A2(i1,i2)* transformation%jacobian(eta1,eta2)
        dphi_eta2 = A1(i1,i2)* transformation%jacobian(eta1,eta2)
        data(i1) = (( jac_m(2,2)*dphi_eta1 - jac_m(2,1)*dphi_eta2 )**2 + &
        ( -jac_m(1,2)*dphi_eta1 + jac_m(1,1)*dphi_eta2 )**2) &
        /abs(transformation%jacobian(eta1,eta2)) 
      enddo
      e_array(i2) = sll_f_compute_integral_trapezoid_1d(data, Nc_eta1+1, delta_eta1)

      do i1=1,Nc_eta1+1
       linf = max(linf,abs(f(i1,i2)))
      enddo
         
    enddo     

    mass = sll_f_compute_integral_trapezoid_1d(mass_array, Nc_eta2+1, delta_eta2)
    l1 = sll_f_compute_integral_trapezoid_1d(l1_array, Nc_eta2+1, delta_eta2)
    l2 = sll_f_compute_integral_trapezoid_1d(l2_array, Nc_eta2+1, delta_eta2)
    l2 = sqrt(l2)
    e = sll_f_compute_integral_trapezoid_1d(e_array, Nc_eta2+1, delta_eta2)

    !mass = mass*delta_eta2
    !l1 = l1*delta_eta2
    !l2 = sqrt(l2*delta_eta2)
    !e  = e*delta_eta2
    
    
    !now, compute errors


    err_linf  = 0.0_f64
    
    do i2 = 1, Nc_eta2+1
      eta2 = eta2_min + (i2-1)* delta_eta2 

      do i1=1,Nc_eta1+1
        eta1 = eta1_min + (i1-1)* delta_eta1
        data(i1) = abs(f(i1,i2)-f_exact(i1,i2))*abs(transformation%jacobian(eta1,eta2))
      enddo
      l1_array(i2) = sll_f_compute_integral_trapezoid_1d(data, Nc_eta1+1, delta_eta1)

      do i1=1,Nc_eta1+1
        eta1 = eta1_min + (i1-1)* delta_eta1
        data(i1) = (f(i1,i2)-f_exact(i1,i2))**2 *abs(transformation%jacobian(eta1,eta2))
      enddo
      l2_array(i2) = sll_f_compute_integral_trapezoid_1d(data, Nc_eta1+1, delta_eta1)


      do i1=1,Nc_eta1+1
       err_linf = max(err_linf,abs(f(i1,i2)-f_exact(i1,i2)))
      enddo
         
    enddo     

    err_l1 = sll_f_compute_integral_trapezoid_1d(l1_array, Nc_eta2+1, delta_eta2)
    err_l2 = sll_f_compute_integral_trapezoid_1d(l2_array, Nc_eta2+1, delta_eta2)
    err_l2 = sqrt(err_l2)
    
    write(file_id,*) &
      dt*real(step,f64), &
      linf, &
      l1, &
      l2, &
      mass, &
      e, &
      maxval(abs(phi(1:Nc_eta1+1,1:Nc_eta2+1))), &
      err_linf, &
      err_l1, &
      err_l2
      
   
    
  end subroutine time_history_diagnostic_curvilinear

  ! periodic case
  ! called when bc_type = sll_p_periodic
  function process_outside_point_periodic1( eta, eta_min, eta_max ) result(eta_out)
      use sll_m_working_precision
      sll_real64, intent(in)  :: eta
      sll_real64, intent(in) :: eta_min
      sll_real64, intent(in) :: eta_max
      sll_real64 :: eta_out

      eta_out = (eta-eta_min)/(eta_max-eta_min)      
      eta_out = eta_out-real(floor(eta_out),f64)
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
  ! called when bc_type = sll_p_set_to_limit
  
  function process_outside_point_set_to_limit1( eta, eta_min, eta_max ) result(eta_out)
      use sll_m_working_precision
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

  subroutine compute_deriv_csl_periodic(input,num_cell,output,p)
    sll_real64, dimension(:), intent(in) :: input
    sll_int32, intent(in) :: num_cell
    sll_real64, dimension(:), intent(out) :: output
    sll_int32, intent(in) :: p
    sll_int32 :: r
    sll_int32 :: s
    sll_int32 :: ierr
    sll_real64, dimension(:), allocatable :: ww
    sll_int32 :: i
    sll_int32 :: ii
    sll_int32 :: ind
    sll_real64 :: tmp
    !sll_real64, dimension(:,:), allocatable :: buf

    r=-p/2
    s=(p+1)/2


    SLL_ALLOCATE(ww(r:s-1),ierr)
    !SLL_ALLOCATE(buf(2,num_cell),ierr)
    
    call compute_ww_test(ww,r,s)
    
    !print *,'ww=',ww
    !print *,'r,s=',r,s
    
    !stop
    
    do i=1,num_cell
      tmp = 0._f64
      do ii=r,s-1
        ind = modulo(i+ii-1,num_cell)+1
        tmp=tmp+ww(r+s-1-ii)*input(ind)
      enddo
      output(i) = tmp
      !buf(1,i) = tmp
      !do ii=r,s-1
      !  ind = modulo(i+ii-1,num_cell)+1
      !  tmp=tmp+ww(ii)*input(ind)
      !enddo
      !buf(2,i) = tmp      
    enddo
    output(num_cell+1) = output(1)
    !output(1) = 0.5_f64*(buf(1,1)+buf(2,num_cell))
    !output(num_cell+1) = output(1)
    !do i=2,num_cell
    !  output(i) = 0.5_f64*(buf(2,i-1)+buf(1,i))
    !enddo
    
  end subroutine compute_deriv_csl_periodic
!        tmp=0._f64
!        do ii=r,s-1
!          i3=i+ii-1;if(i3<=0)i3=0;if(i3>=N0)i3=N0          
!          tmp=tmp+ww(r+s-1-ii)*f(i3,j)
!        enddo
!        aretesvg(i,j)=tmp
!        tmp=0._f64
!        do ii=r,s-1
!          i3=i+ii;if(i3<=0)i3=0;if(i3>=N0)i3=N0          
!          tmp=tmp+ww(ii)*f(i3,j)
!        enddo
!        aretesvd(i,j)=tmp



  subroutine compute_ww_test(ww,r,s)
    sll_int32, intent(in) :: r
    sll_int32, intent(in) :: s
    sll_real64, dimension(r:s-1), intent(out) :: ww
    sll_real64, dimension(:), allocatable :: w
    sll_real64 :: tmp
    sll_int32 :: i
    sll_int32 :: j
    !sll_int32 :: ii
    sll_int32 :: ierr

    SLL_ALLOCATE(w(r:s),ierr)

     !maple code for generation of w
    !for k from r to -1 do
    !  C[k]:=product((k-j),j=r..k-1)*product((k-j),j=k+1..s):
    !  C[k]:=1/C[k]*product((-j),j=r..k-1)*product((-j),j=k+1..-1)*product((-j),j=1..s):
    !od:
    !for k from 1 to s do
    !  C[k]:=product((k-j),j=r..k-1)*product((k-j),j=k+1..s):
    !  C[k]:=1/C[k]*product((-j),j=r..-1)*product((-j),j=1..k-1)*product((-j),j=k+1..s):
    !od:
    !C[0]:=-add(C[k],k=r..-1)-add(C[k],k=1..s):
    
    do i=r,-1
      tmp=1._f64
      do j=r,i-1
        tmp=tmp*real(i-j,f64)
      enddo
      do j=i+1,s
        tmp=tmp*real(i-j,f64)
      enddo
      tmp=1._f64/tmp
      do j=r,i-1
        tmp=tmp*real(-j,f64)
      enddo
      do j=i+1,-1
        tmp=tmp*real(-j,f64)
      enddo
      do j=1,s
        tmp=tmp*real(-j,f64)
      enddo
      w(i)=tmp      
     enddo

!    do i=r,-1
!      tmp=1._f64
!      !do j=r,i-1
!      !  tmp=tmp*real(i-j,f64)
!      !enddo
!      !do j=i+1,s
!      !  tmp=tmp*real(i-j,f64)
!      !enddo
!      !tmp=1._f64/tmp
!      do j=r,i-1 !-j/(i-j)=j/(j-i)=1/(1-i/j)
!        tmp=tmp*(1._f64-real(i,f64)/real(j,f64))
!      enddo
!      do j=i+1,-1
!        tmp=tmp*(1._f64-real(i,f64)/real(j,f64))
!      enddo
!      do j=1,s
!        tmp=tmp*(1._f64-real(i,f64)/real(j,f64))
!      enddo
!      tmp=tmp*real(i,f64)
!      w(i)=1._f64/tmp      
!    enddo


    do i=1,s
      tmp=1._f64
      do j=r,i-1
        tmp=tmp*real(i-j,f64)
      enddo
      do j=i+1,s
        tmp=tmp*real(i-j,f64)
      enddo
      tmp=1._f64/tmp
      do j=r,-1
        tmp=tmp*real(-j,f64)
      enddo
      do j=1,i-1
        tmp=tmp*real(-j,f64)
      enddo
      do j=i+1,s
        tmp=tmp*real(-j,f64)
      enddo
      w(i)=tmp      
    enddo

    tmp=0._f64
    do i=r,-1
      tmp=tmp+w(i)
    enddo
    do i=1,s
      tmp=tmp+w(i)
    enddo
    w(0)=-tmp
    
    
    
    !print *,'w',w
    !do ii=r,s
    !  print *,ii,w(r+s-ii)
    !enddo
    
    !compute now ww
    !maple code
    !#for conservative formulation
    !tmp:=0:
    !for k from r to -1 do
    !tmp:=tmp+C[k]:
    !CC[k]:=-tmp:
    !od:
    !tmp:=0:
    !for k from s to 1 by -1 do
    !  tmp:=tmp+C[k]:
    !  CC[k-1]:=tmp:
    !od:
    !seq(CC[k],k=r..s-1);
    !evalf(%);

    tmp=0._f64
    do i=r,-1
      tmp=tmp+w(i)
      ww(i)=-tmp
    enddo
    tmp=0._f64
    do i=s,1,-1
      tmp=tmp+w(i)
      ww(i-1)=tmp
    enddo

    !print *,'ww',ww(r:s-1)
    !do ii=r,s-1
    !  print *,ii,ww(r+s-1-ii)
    !enddo
    !stop

    SLL_DEALLOCATE_ARRAY(w,ierr)
    
  end subroutine compute_ww_test



  subroutine sll_s_delete_analytic_field_2d_curvilinear( sim )

    class(sll_t_simulation_2d_analytic_field_curvilinear) :: sim
    !sll_int32 :: ierr
    
            
  end subroutine sll_s_delete_analytic_field_2d_curvilinear




end module sll_m_sim_bsl_ad_2d0v_curv
