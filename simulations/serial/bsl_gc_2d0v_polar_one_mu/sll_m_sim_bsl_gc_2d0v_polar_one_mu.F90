module sll_m_sim_bsl_gc_2d0v_polar_one_mu

!2d guiding center polar simulation
!related to simulation_2d_guiding_center_curvilinear.F90
!but here geometry and test are specifically polar

!see ../selalib/prototype/src/simulation/gcsim2d_polar_input.nml
!for example of use

!contact: Michel Mehrenberger (mehrenbe@math.unistra.fr)
!         Adnane Hamiaz (hamiaz@math.unistra.fr)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

   use sll_m_advection_2d_base, only: &
      sll_c_advector_2d

   use sll_m_advection_2d_bsl, only: &
      sll_f_new_advector_2d_bsl

   use sll_m_boundary_condition_descriptors, only: &
      sll_p_dirichlet, &
      sll_p_hermite, &
      sll_p_periodic, &
      sll_p_set_to_limit

   use sll_m_cartesian_meshes, only: &
      sll_f_new_cartesian_mesh_2d, &
      sll_t_cartesian_mesh_2d

   use sll_m_characteristics_2d_base, only: &
      sll_c_characteristics_2d_base

   use sll_m_characteristics_2d_explicit_euler, only: &
      sll_f_new_explicit_euler_2d_charac

   use sll_m_characteristics_2d_verlet, only: &
      sll_f_new_verlet_2d_charac

   use sll_m_common_array_initializers, only: &
      sll_f_diocotron_initializer_2d, &
      sll_i_scalar_initializer_2d

   use sll_m_common_coordinate_transformations, only: &
      sll_f_polar_jac11, &
      sll_f_polar_jac12, &
      sll_f_polar_jac21, &
      sll_f_polar_jac22, &
      sll_f_polar_x1, &
      sll_f_polar_x2

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

   use sll_m_fft, only: &
      sll_s_fft_exec_r2r_1d, &
      sll_s_fft_free, &
      sll_p_fft_forward, &
      sll_f_fft_get_mode_r2c_1d, &
      sll_s_fft_init_r2r_1d, &
      sll_t_fft

   use sll_m_gyroaverage_2d_base, only: &
      sll_c_gyroaverage_2d_base

   use sll_m_gyroaverage_2d_polar_hermite_solver, only: &
      sll_f_new_gyroaverage_2d_polar_hermite_solver

   use sll_m_gyroaverage_2d_polar_pade_solver, only: &
      sll_f_new_gyroaverage_2d_polar_pade_solver

   use sll_m_gyroaverage_2d_polar_splines_solver, only: &
      sll_f_new_gyroaverage_2d_polar_splines_solver

   use sll_m_hdf5_io_serial, only: &
      sll_t_hdf5_ser_handle, &
      sll_s_hdf5_ser_file_create, &
      sll_s_hdf5_ser_file_close, &
      sll_o_hdf5_ser_write_array

   use sll_m_interpolators_1d_base, only: &
      sll_c_interpolator_1d

   use sll_m_interpolators_2d_base, only: &
      sll_c_interpolator_2d

   use sll_m_poisson_2d_base, only: &
      sll_c_poisson_2d_base

   use sll_m_poisson_2d_polar, only: &
      sll_f_new_poisson_2d_polar

   use sll_m_reduction, only: &
      sll_f_compute_integral_trapezoid_1d

   use sll_m_sim_base, only: &
      sll_c_simulation_base_class

   use sll_m_utilities, only: &
      sll_s_int2string

   use sll_m_xdmf, only: &
      sll_s_xdmf_close, &
      sll_o_xdmf_open, &
      sll_o_xdmf_write_array

#ifdef MUDPACK
   use sll_m_poisson_2d_mudpack_curvilinear
   use sll_m_poisson_2d_mudpack_curvilinear_old, only: &
      sll_f_new_poisson_2d_mudpack_curvilinear_old

#endif
   implicit none

   public :: &
      sll_f_new_guiding_center_2d_polar_one_mu

   private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!#define OLD_POISSON
!#define NEW_POISSON

   sll_int32, parameter :: SLL_EULER = 0
   sll_int32, parameter :: SLL_PREDICTOR_CORRECTOR = 1
   sll_int32, parameter :: SLL_PHI_FROM_RHO = 0
   sll_int32, parameter :: SLL_E_FROM_RHO = 1

   type, extends(sll_c_simulation_base_class) :: &
      sll_simulation_2d_guiding_center_polar_one_mu

      !geometry
      type(sll_t_cartesian_mesh_2d), pointer :: mesh_2d

      !initial function
      procedure(sll_i_scalar_initializer_2d), nopass, pointer :: init_func
      sll_real64, dimension(:), pointer :: params

      !advector
      class(sll_c_advector_2d), pointer    :: advect_2d

      !interpolator for derivatives
      class(sll_c_interpolator_2d), pointer   :: phi_interp2d

      !poisson solver
      class(sll_c_poisson_2d_base), pointer   :: poisson
      sll_int32 :: poisson_case

      !gyroaverage
      class(sll_c_gyroaverage_2d_base), pointer :: gyroaverage
      sll_real64  :: mu

      !time_iterations
      sll_real64 :: dt
      sll_int32  :: num_iterations
      sll_int32  :: freq_diag
      sll_int32  :: freq_diag_time

      !time_loop
      sll_int32 :: time_loop_case

   contains
      procedure, pass(sim) :: run => run_gc2d_polar_one_mu
      procedure, pass(sim) :: init_from_file => init_fake

   end type sll_simulation_2d_guiding_center_polar_one_mu

!  abstract interface
!    function sll_i_scalar_initializer_2d( x1, x2, params )
!      use sll_m_working_precision
!      sll_real64                                     :: sll_i_scalar_initializer_2d
!      sll_real64, intent(in)                         :: x1
!      sll_real64, intent(in)                         :: x2
!      sll_real64, dimension(:), intent(in), optional :: params
!    end function sll_i_scalar_initializer_2d
!  end interface

contains

   function sll_f_new_guiding_center_2d_polar_one_mu(filename) result(sim)
      type(sll_simulation_2d_guiding_center_polar_one_mu), pointer :: sim
      character(len=*), intent(in), optional :: filename
      sll_int32 :: ierr

      SLL_ALLOCATE(sim, ierr)

      call initialize_guiding_center_2d_polar_one_mu(sim, filename)

   end function sll_f_new_guiding_center_2d_polar_one_mu

   subroutine initialize_guiding_center_2d_polar_one_mu(sim, filename)
      class(sll_simulation_2d_guiding_center_polar_one_mu), intent(inout) :: sim

      character(len=*), intent(in), optional :: filename
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

      !poisson
      character(len=256) :: poisson_case
      character(len=256) :: poisson_solver
      !character(len=256) :: mudpack_method
      sll_int32 :: spline_degree_eta1
      sll_int32 :: spline_degree_eta2

      !gyroaverage
      character(len=256)      :: gyroaverage_case
      sll_real64              :: mu
      sll_int32               :: gyroaverage_N_points
      sll_int32               :: gyroaverage_interp_degree_x1
      sll_int32               :: gyroaverage_interp_degree_x2
      sll_real64              :: eta_min_gyro(2)
      sll_real64              :: eta_max_gyro(2)
      sll_int32               :: Nc_gyro(2)
      sll_int32               :: interp_degree_gyro(2)

      !local variables
      sll_int32 :: Nc_x1
      sll_int32 :: Nc_x2
      sll_real64 :: x1_min
      sll_real64 :: x1_max
      sll_real64 :: x2_min
      sll_real64 :: x2_max
      class(sll_c_interpolator_2d), pointer :: f_interp2d
      class(sll_c_interpolator_2d), pointer :: phi_interp2d
      class(sll_c_characteristics_2d_base), pointer :: charac2d
      class(sll_c_interpolator_2d), pointer   :: A1_interp2d
      class(sll_c_interpolator_2d), pointer   :: A2_interp2d
      class(sll_c_interpolator_1d), pointer   :: A1_interp1d_x1
      class(sll_c_interpolator_1d), pointer   :: A2_interp1d_x1
#ifdef MUDPACK
      sll_real64, dimension(:, :), pointer :: b11
      sll_real64, dimension(:, :), pointer :: b12
      sll_real64, dimension(:, :), pointer :: b21
      sll_real64, dimension(:, :), pointer :: b22
      !sll_real64, dimension(:,:), pointer :: b1
      !sll_real64, dimension(:,:), pointer :: b2
      sll_real64, dimension(:, :), pointer :: c
      class(sll_c_coordinate_transformation_2d_base), pointer :: transformation
#endif /* MUDPACK */
!    sll_real64, dimension(:,:), allocatable :: cxx_2d
!    sll_real64, dimension(:,:), allocatable :: cxy_2d
!    sll_real64, dimension(:,:), allocatable :: cyy_2d
!    sll_real64, dimension(:,:), allocatable :: cx_2d
!    sll_real64, dimension(:,:), allocatable :: cy_2d
!    sll_real64, dimension(:,:), allocatable :: ce_2d
      sll_int32 :: ierr

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
         A_interp_case

      namelist /poisson/ &
         poisson_case, &
         poisson_solver, &
         spline_degree_eta1, &
         spline_degree_eta2

      namelist /gyroaverage/ &
         gyroaverage_case, &
         mu, &
         gyroaverage_N_points, &
         gyroaverage_interp_degree_x1, &
         gyroaverage_interp_degree_x2

      !set default parameters

      !geometry
      mesh_case = "SLL_POLAR_MESH"
      num_cells_x1 = 32
      r_min = 1.0_f64
      r_max = 10._f64
      num_cells_x2 = 32

      !initial function
      initial_function_case = "SLL_DIOCOTRON"
      r_minus = 4._f64
      r_plus = 5._f64
      kmode_x2 = 3._f64
      eps = 1.e-6_f64

      !time_iterations
      dt = 0.1_f64
      number_iterations = 600
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
      poisson_solver = "SLL_POLAR_FFT"
      !poisson_solver = "SLL_ELLIPTIC_FINITE_ELEMENT_SOLVER" !use with "SLL_PHI_FROM_RHO"
      !poisson_solver = "SLL_MUDPACK_CURVILINEAR"   !use with "SLL_PHI_FROM_RHO"
      spline_degree_eta1 = 3
      spline_degree_eta2 = 3

      !gyroaverage
      gyroaverage_case = "HERMITE_C1_PRECOMPUTE"
      !gyroaverage_case = "SPLINES"
      !gyroaverage_case = "PADE"
      mu = 1.0_f64
      gyroaverage_N_points = 16
      gyroaverage_interp_degree_x1 = 3
      gyroaverage_interp_degree_x2 = 3

      if (present(filename)) then
         open (unit=input_file, file=trim(filename)//'.nml', IOStat=IO_stat)
         if (IO_stat /= 0) then
            print *, '#initialize_guiding_center_2d_cartesian() failed to open file ', &
               trim(filename)//'.nml'
            STOP
         end if
         print *, '#initialization with filename:'
         print *, '#', trim(filename)//'.nml'
         read (input_file, geometry)
         read (input_file, initial_function)
         read (input_file, time_iterations)
         read (input_file, advector)
         read (input_file, poisson)
         read (input_file, gyroaverage)
         close (input_file)
      else
         print *, '#initialization with default parameters'
      end if

      select case (mesh_case)
      case ("SLL_POLAR_MESH")
         x1_min = r_min
         x1_max = r_max
         x2_min = 0._f64
         x2_max = 2._f64*sll_p_pi
      case default
         print *, '#bad mesh_case', mesh_case
         print *, '#not implemented'
         print *, '#in initialize_guiding_center_2d_polar'
         stop
      end select

      Nc_x1 = num_cells_x1
      Nc_x2 = num_cells_x2

      sim%dt = dt
      sim%num_iterations = number_iterations
      sim%freq_diag = freq_diag
      sim%freq_diag_time = freq_diag_time

      sim%mesh_2d => sll_f_new_cartesian_mesh_2d( &
                     Nc_x1, &
                     Nc_x2, &
                     eta1_min=x1_min, &
                     eta1_max=x1_max, &
                     eta2_min=x2_min, &
                     eta2_max=x2_max)

      select case (f_interp2d_case)
      case ("SLL_CUBIC_SPLINES")
         f_interp2d => sll_f_new_cubic_spline_interpolator_2d( &
                       Nc_x1 + 1, &
                       Nc_x2 + 1, &
                       x1_min, &
                       x1_max, &
                       x2_min, &
                       x2_max, &
                       sll_p_hermite, &
                       sll_p_periodic)
      case default
         print *, '#bad f_interp2d_case', f_interp2d_case
         print *, '#not implemented'
         print *, '#in initialize_guiding_center_2d_polar'
         stop
      end select

      select case (A_interp_case)
      case ("SLL_CUBIC_SPLINES")
         A1_interp2d => sll_f_new_cubic_spline_interpolator_2d( &
                        Nc_x1 + 1, &
                        Nc_x2 + 1, &
                        x1_min, &
                        x1_max, &
                        x2_min, &
                        x2_max, &
                        sll_p_hermite, &
                        sll_p_periodic)
         A2_interp2d => sll_f_new_cubic_spline_interpolator_2d( &
                        Nc_x1 + 1, &
                        Nc_x2 + 1, &
                        x1_min, &
                        x1_max, &
                        x2_min, &
                        x2_max, &
                        sll_p_hermite, &
                        sll_p_periodic)
         A1_interp1d_x1 => sll_f_new_cubic_spline_interpolator_1d( &
                           Nc_x1 + 1, &
                           x1_min, &
                           x1_max, &
                           sll_p_hermite)
         A2_interp1d_x1 => sll_f_new_cubic_spline_interpolator_1d( &
                           Nc_x1 + 1, &
                           x1_min, &
                           x1_max, &
                           sll_p_hermite)
      case default
         print *, '#bad A_interp_case', A_interp_case
         print *, '#not implemented'
         print *, '#in initialize_guiding_center_2d_polar'
         stop
      end select

      select case (phi_interp2d_case)
      case ("SLL_CUBIC_SPLINES")
         phi_interp2d => sll_f_new_cubic_spline_interpolator_2d( &
                         Nc_x1 + 1, &
                         Nc_x2 + 1, &
                         x1_min, &
                         x1_max, &
                         x2_min, &
                         x2_max, &
                         sll_p_hermite, &
                         sll_p_periodic)
      case default
         print *, '#bad phi_interp2d_case', phi_interp2d_case
         print *, '#not implemented'
         print *, '#in initialize_guiding_center_2d_polar'
         stop
      end select

      select case (charac2d_case)
      case ("SLL_EULER")
         charac2d => sll_f_new_explicit_euler_2d_charac( &
                     Nc_x1 + 1, &
                     Nc_x2 + 1, &
                     eta1_min=x1_min, &
                     eta1_max=x1_max, &
                     eta2_min=x2_min, &
                     eta2_max=x2_max, &
                     bc_type_1=sll_p_set_to_limit, &
                     bc_type_2=sll_p_periodic)
      case ("SLL_VERLET")
         charac2d => sll_f_new_verlet_2d_charac( &
                     Nc_x1 + 1, &
                     Nc_x2 + 1, &
                     A1_interp2d, &
                     A2_interp2d, &
                     A1_interp1d_x1, &
                     A2_interp1d_x1, &
                     bc_type_1=sll_p_set_to_limit, &
                     bc_type_2=sll_p_periodic, &
                     eta1_min=x1_min, &
                     eta1_max=x1_max, &
                     eta2_min=x2_min, &
                     eta2_max=x2_max)
      case default
         print *, '#bad charac2d_case', charac2d_case
         print *, '#not implemented'
         print *, '#in initialize_guiding_center_2d_polar'
         stop
      end select

      sim%phi_interp2d => phi_interp2d

      select case (advect2d_case)
      case ("SLL_BSL")
         sim%advect_2d => sll_f_new_advector_2d_bsl( &
                          f_interp2d, &
                          charac2d, &
                          Nc_x1 + 1, &
                          Nc_x2 + 1, &
                          eta1_min=x1_min, &
                          eta1_max=x1_max, &
                          eta2_min=x2_min, &
                          eta2_max=x2_max)
      case default
         print *, '#bad advect_case', advect2d_case
         print *, '#not implemented'
         print *, '#in initialize_guiding_center_2d_polar'
         stop
      end select

      select case (initial_function_case)
      case ("SLL_DIOCOTRON")
         sim%init_func => sll_f_diocotron_initializer_2d
         SLL_ALLOCATE(sim%params(4), ierr)
         sim%params(1) = r_minus
         sim%params(2) = r_plus
         sim%params(3) = eps
         sim%params(4) = kmode_x2
      case default
         print *, '#bad initial_function_case', initial_function_case
         print *, '#not implemented'
         print *, '#in initialize_guiding_center_2d_polar'
         stop
      end select

      !time_loop
      select case (time_loop_case)
      case ("SLL_EULER")
         sim%time_loop_case = SLL_EULER
      case ("SLL_PREDICTOR_CORRECTOR")
         sim%time_loop_case = SLL_PREDICTOR_CORRECTOR
      case default
         print *, '#bad time_loop_case', time_loop_case
         print *, '#not implemented'
         print *, '#in initialize_guiding_center_2d_polar'
         stop
      end select

      select case (poisson_case)
      case ("SLL_PHI_FROM_RHO")
         sim%poisson_case = SLL_PHI_FROM_RHO
         !case ("SLL_E_FROM_RHO")
         !  sim%poisson_case = SLL_E_FROM_RHO
      case default
         print *, '#bad poisson_case', poisson_case
         print *, '#not implemented'
         print *, '#in initialize_guiding_center_2d_cartesian'
         stop
      end select

      select case (poisson_solver)
      case ("SLL_POLAR_FFT")
         sim%poisson => sll_f_new_poisson_2d_polar( &
                        x1_min, &
                        x1_max, &
                        Nc_x1, &
                        Nc_x2, &
                        !(/sll_p_neumann_mode_0, sll_p_dirichlet/))
                        (/sll_p_dirichlet, sll_p_dirichlet/))
!      case ("SLL_ELLIPTIC_FINITE_ELEMENT_SOLVER")
!        transformation => sll_f_new_coordinate_transformation_2d_analytic( &
!          "analytic_polar_transformation", &
!          sim%mesh_2d, &
!          sll_f_polar_x1, &
!          sll_f_polar_x2, &
!          sll_f_polar_jac11, &
!          sll_f_polar_jac12, &
!          sll_f_polar_jac21, &
!          sll_f_polar_jac22, &
!          params=(/0._f64,0._f64,0._f64,0._f64/))
!
!        SLL_ALLOCATE(b11(Nc_x1+1,Nc_x2+1),ierr)
!        SLL_ALLOCATE(b12(Nc_x1+1,Nc_x2+1),ierr)
!        SLL_ALLOCATE(b21(Nc_x1+1,Nc_x2+1),ierr)
!        SLL_ALLOCATE(b22(Nc_x1+1,Nc_x2+1),ierr)
!        SLL_ALLOCATE(b1(Nc_x1+1,Nc_x2+1),ierr)
!        SLL_ALLOCATE(b2(Nc_x1+1,Nc_x2+1),ierr)
!        SLL_ALLOCATE(c(Nc_x1+1,Nc_x2+1),ierr)
!
!        b11 = 1._f64
!        b22 = 1._f64
!        b12 = 0._f64
!        b21 = 0._f64
!        c = 0._f64
!
!        sim%poisson => new_poisson_2d_elliptic_solver( &
!         transformation,&
!         spline_degree_eta1, &
!         spline_degree_eta2, &
!         Nc_x1, &
!         Nc_x2, &
!         ES_GAUSS_LEGENDRE, &
!         ES_GAUSS_LEGENDRE, &
!         sll_p_dirichlet, &
!         sll_p_dirichlet, &
!         sll_p_periodic, &
!         sll_p_periodic, &
!         x1_min, &
!         x1_max, &
!         x2_min, &
!         x2_max, &
!         b11, &
!         b12, &
!         b21, &
!         b22, &
!         b1, &
!         b2, &
!         c )
#ifdef MUDPACK
      case ("SLL_MUDPACK_CURVILINEAR")
         transformation => sll_f_new_coordinate_transformation_2d_analytic( &
                           "analytic_polar_transformation", &
                           sim%mesh_2d, &
                           sll_f_polar_x1, &
                           sll_f_polar_x2, &
                           sll_f_polar_jac11, &
                           sll_f_polar_jac12, &
                           sll_f_polar_jac21, &
                           sll_f_polar_jac22, &
                           params=(/0._f64, 0._f64, 0._f64, 0._f64/))

         SLL_ALLOCATE(b11(Nc_x1 + 1, Nc_x2 + 1), ierr)
         SLL_ALLOCATE(b12(Nc_x1 + 1, Nc_x2 + 1), ierr)
         SLL_ALLOCATE(b21(Nc_x1 + 1, Nc_x2 + 1), ierr)
         SLL_ALLOCATE(b22(Nc_x1 + 1, Nc_x2 + 1), ierr)
         SLL_ALLOCATE(c(Nc_x1 + 1, Nc_x2 + 1), ierr)

         b11 = 1._f64
         b22 = 1._f64
         b12 = 0._f64
         b21 = 0._f64
         c = 0._f64

         sim%poisson => sll_f_new_poisson_2d_mudpack_curvilinear_old( &
                        transformation, &
                        x1_min, &
                        x1_max, &
                        Nc_x1, &
                        x2_min, &
                        x2_max, &
                        Nc_x2, &
                        sll_p_dirichlet, &
                        sll_p_dirichlet, &
                        sll_p_periodic, &
                        sll_p_periodic, &
                        sll_p_hermite, &
                        sll_p_periodic, &
                        b11, &
                        b12, &
                        b21, &
                        b22, &
                        c)

#endif

      case default
         print *, '#bad poisson_solver', poisson_solver
         print *, '#not implemented'
         print *, '#in initialize_guiding_center_2d_polar'
         stop
      end select

      !gyroaverage

      sim%mu = mu

      eta_min_gyro(1) = sim%mesh_2d%eta1_min
      eta_max_gyro(1) = sim%mesh_2d%eta1_min + real(sim%mesh_2d%num_cells1, f64)*sim%mesh_2d%delta_eta1
      eta_min_gyro(2) = 0._f64
      eta_max_gyro(2) = 2._f64*sll_p_pi

      Nc_gyro(1) = sim%mesh_2d%num_cells1
      Nc_gyro(2) = sim%mesh_2d%num_cells2
      interp_degree_gyro(1) = gyroaverage_interp_degree_x1
      interp_degree_gyro(2) = gyroaverage_interp_degree_x2

      select case (gyroaverage_case)
      case ("HERMITE")

         sim%gyroaverage => sll_f_new_gyroaverage_2d_polar_hermite_solver( &
                            eta_min_gyro, &
                            eta_max_gyro, &
                            Nc_gyro, &
                            gyroaverage_N_points, &
                            interp_degree_gyro, &
                            sqrt(2*sim%mu), &
                            1)

      case ("HERMITE_C1")

         sim%gyroaverage => sll_f_new_gyroaverage_2d_polar_hermite_solver( &
                            eta_min_gyro, &
                            eta_max_gyro, &
                            Nc_gyro, &
                            gyroaverage_N_points, &
                            interp_degree_gyro, &
                            sqrt(2*sim%mu), &
                            2)

      case ("HERMITE_C1_PRECOMPUTE")

         sim%gyroaverage => sll_f_new_gyroaverage_2d_polar_hermite_solver( &
                            eta_min_gyro, &
                            eta_max_gyro, &
                            Nc_gyro, &
                            gyroaverage_N_points, &
                            interp_degree_gyro, &
                            sqrt(2*sim%mu), &
                            3)

      case ("HERMITE_C1_INVARIANCE")

         sim%gyroaverage => sll_f_new_gyroaverage_2d_polar_hermite_solver( &
                            eta_min_gyro, &
                            eta_max_gyro, &
                            Nc_gyro, &
                            gyroaverage_N_points, &
                            interp_degree_gyro, &
                            sqrt(2*sim%mu), &
                            4)

      case ("SPLINES")

         sim%gyroaverage => sll_f_new_gyroaverage_2d_polar_splines_solver( &
                            eta_min_gyro, &
                            eta_max_gyro, &
                            Nc_gyro, &
                            gyroaverage_N_points, &
                            1)

      case ("SPLINES_INVARIANCE")

         sim%gyroaverage => sll_f_new_gyroaverage_2d_polar_splines_solver( &
                            eta_min_gyro, &
                            eta_max_gyro, &
                            Nc_gyro, &
                            gyroaverage_N_points, &
                            2)

      case ("SPLINES_PRECOMPUTE")

         sim%gyroaverage => sll_f_new_gyroaverage_2d_polar_splines_solver( &
                            eta_min_gyro, &
                            eta_max_gyro, &
                            Nc_gyro, &
                            gyroaverage_N_points, &
                            3)

      case ("SPLINES_PRECOMPUTE_FFT")

         sim%gyroaverage => sll_f_new_gyroaverage_2d_polar_splines_solver( &
                            eta_min_gyro, &
                            eta_max_gyro, &
                            Nc_gyro, &
                            gyroaverage_N_points, &
                            4)

      case ("PADE")

         sim%gyroaverage => sll_f_new_gyroaverage_2d_polar_pade_solver( &
                            eta_min_gyro, &
                            eta_max_gyro, &
                            Nc_gyro)

      case default
         print *, '#bad gyroaverage_case', gyroaverage_case
         print *, '#not implemented'
         print *, '#in initialize_guiding_center_2d_polar'
         stop
      end select

      print *, '##Initial condition'
      print *, '#ell=', kmode_x2
      print *, '#r-=', r_minus
      print *, '#r+=', r_plus
      print *, '#eps=', eps
      print *, '##Gyroaverage'
      print *, '#gyroaverage_case=', gyroaverage_case
      print *, '#mu=', mu
      print *, '#gyroaverage_N_points=', gyroaverage_N_points
      print *, '#gyroaverage_interp_degree=', gyroaverage_interp_degree_x1, gyroaverage_interp_degree_x2

   end subroutine initialize_guiding_center_2d_polar_one_mu

   subroutine init_fake(sim, filename)
      class(sll_simulation_2d_guiding_center_polar_one_mu), intent(inout) :: sim
      character(len=*), intent(in)                                :: filename

      print *, '# Do not use the routine init_vp4d_fake'
      print *, '#use instead initialize_vlasov_par_poisson_seq_cart'
      print *, sim%dt
      print *, filename
      stop

   end subroutine init_fake

   subroutine run_gc2d_polar_one_mu(sim)
      class(sll_simulation_2d_guiding_center_polar_one_mu), intent(inout) :: sim
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
      sll_real64, dimension(:, :), pointer :: f
      sll_real64, dimension(:, :), pointer :: Jf
      sll_real64, dimension(:, :), pointer :: f_old
      sll_real64, dimension(:, :), pointer :: Jf_old
      sll_real64, dimension(:, :), pointer :: phi
      sll_real64, dimension(:, :), pointer :: A1 !advection fields
      sll_real64, dimension(:, :), pointer :: A2
      sll_int32 :: ierr
      sll_int32 :: nb_step
      sll_int32 :: step
      sll_real64 :: dt
      sll_int32 :: thdiag_id = 99
      sll_int32             :: IO_stat
      sll_int32 :: iplot

      Nc_x1 = sim%mesh_2d%num_cells1
      Nc_x2 = sim%mesh_2d%num_cells2
      delta_x1 = sim%mesh_2d%delta_eta1
      delta_x2 = sim%mesh_2d%delta_eta2
      x1_min = sim%mesh_2d%eta1_min
      x2_min = sim%mesh_2d%eta2_min
      nb_step = sim%num_iterations
      dt = sim%dt

      !allocation
      SLL_ALLOCATE(f(Nc_x1 + 1, Nc_x2 + 1), ierr)
      SLL_ALLOCATE(Jf(Nc_x1 + 1, Nc_x2 + 1), ierr)
      SLL_ALLOCATE(f_old(Nc_x1 + 1, Nc_x2 + 1), ierr)
      SLL_ALLOCATE(Jf_old(Nc_x1 + 1, Nc_x2 + 1), ierr)
      SLL_ALLOCATE(phi(Nc_x1 + 1, Nc_x2 + 1), ierr)
      SLL_ALLOCATE(A1(Nc_x1 + 1, Nc_x2 + 1), ierr)
      SLL_ALLOCATE(A2(Nc_x1 + 1, Nc_x2 + 1), ierr)

      !initialisation of distribution function
      do i2 = 1, Nc_x2 + 1
         x2 = x2_min + real(i2 - 1, f64)*delta_x2
         do i1 = 1, Nc_x1 + 1
            x1 = x1_min + real(i1 - 1, f64)*delta_x1
            f(i1, i2) = sim%init_func(x1, x2, sim%params)
         end do
      end do

      Jf = f

      call sim%gyroaverage%compute_gyroaverage( &
         sqrt(2*sim%mu), &
         Jf(1:Nc_x1 + 1, 1:Nc_x2 + 1))

      call sim%poisson%compute_phi_from_rho(phi, Jf)

      call sim%gyroaverage%compute_gyroaverage( &
         sqrt(2*sim%mu), &
         phi(1:Nc_x1 + 1, 1:Nc_x2 + 1))

      call compute_field_from_phi_2d_polar(phi, sim%mesh_2d, A1, A2, sim%phi_interp2d)

      open (unit=thdiag_id, file='thdiag.dat', IOStat=IO_stat)
      if (IO_stat /= 0) then
         print *, '#run_gc2d_polar(sim) failed to open file thdiag.dat'
         STOP
      end if

      iplot = 0

      do step = 1, nb_step + 1
         f_old = f

         Jf_old = f_old

         call sim%gyroaverage%compute_gyroaverage( &
            sqrt(2*sim%mu), &
            Jf_old(1:Nc_x1 + 1, 1:Nc_x2 + 1))

         call sim%poisson%compute_phi_from_rho(phi, Jf_old)

         call sim%gyroaverage%compute_gyroaverage( &
            sqrt(2*sim%mu), &
            phi(1:Nc_x1 + 1, 1:Nc_x2 + 1))

         call compute_field_from_phi_2d_polar(phi, sim%mesh_2d, A1, A2, sim%phi_interp2d)

         if (modulo(step - 1, sim%freq_diag_time) == 0) then
            call time_history_diagnostic_gc_polar( &
               thdiag_id, &
               step - 1, &
               dt, &
               sim%mesh_2d, &
               f, &
               phi, &
               A1, &
               A2)
         end if

#ifndef NOHDF5
         if (modulo(step - 1, sim%freq_diag) == 0) then
            print *, "#step= ", step
            call plot_f_polar(iplot, f, sim%mesh_2d)
            iplot = iplot + 1
         end if
#endif

         select case (sim%time_loop_case)
         case (SLL_EULER)
            call sim%advect_2d%advect_2d(A1, A2, sim%dt, f_old, f)
         case (SLL_PREDICTOR_CORRECTOR)
            call sim%advect_2d%advect_2d(A1, A2, 0.5_f64*sim%dt, f_old, f)
            Jf = f
            call sim%gyroaverage%compute_gyroaverage( &
               sqrt(2*sim%mu), &
               Jf(1:Nc_x1 + 1, 1:Nc_x2 + 1))
            call sim%poisson%compute_phi_from_rho(phi, Jf)
            call sim%gyroaverage%compute_gyroaverage( &
               sqrt(2*sim%mu), &
               phi(1:Nc_x1 + 1, 1:Nc_x2 + 1))
            call compute_field_from_phi_2d_polar(phi, sim%mesh_2d, A1, A2, sim%phi_interp2d)
            call sim%advect_2d%advect_2d(A1, A2, sim%dt, f_old, f)
         case default
            print *, '#bad time_loop_case', sim%time_loop_case
            print *, '#not implemented'
            print *, '#in run_gc2d_polar'
            print *, '#available options are:'
            print *, '#SLL_EULER=', SLL_EULER
            print *, '#SLL_PREDICTOR_CORRECTOR=', SLL_PREDICTOR_CORRECTOR

         end select

      end do

      close (thdiag_id)

      !print *,'#not implemented for the moment!'
   end subroutine run_gc2d_polar_one_mu

   subroutine compute_field_from_phi_2d_polar(phi, mesh_2d, A1, A2, interp2d)
      sll_real64, dimension(:, :), intent(in) :: phi
      sll_real64, dimension(:, :), intent(out) :: A1
      sll_real64, dimension(:, :), intent(out) :: A2
      type(sll_t_cartesian_mesh_2d), pointer :: mesh_2d
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

      do i2 = 1, Nc_x2 + 1
         x2 = x2_min + real(i2 - 1, f64)*delta_x2
         do i1 = 1, Nc_x1 + 1
            x1 = x1_min + real(i1 - 1, f64)*delta_x1
            A1(i1, i2) = interp2d%interpolate_from_interpolant_derivative_eta2(x1, x2)/x1
            A2(i1, i2) = -interp2d%interpolate_from_interpolant_derivative_eta1(x1, x2)/x1
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
      type(sll_t_cartesian_mesh_2d), pointer :: mesh_2d
      sll_real64, dimension(:, :), intent(in) :: f
      sll_real64, dimension(:, :), intent(in) :: phi
      sll_real64, dimension(:, :), intent(in) :: A1
      sll_real64, dimension(:, :), intent(in) :: A2
      sll_real64 :: time_mode(8)
      !sll_real64 :: mode_slope(8)
      sll_real64 :: w
      sll_real64 :: l1
      sll_real64 :: l2
      sll_real64 :: e
      sll_real64, dimension(:), allocatable :: int_r
      sll_real64, dimension(:), allocatable :: data
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
      type(sll_t_fft)        :: pfwd

      Nc_x1 = mesh_2d%num_cells1
      Nc_x2 = mesh_2d%num_cells2

      x1_min = mesh_2d%eta1_min
      x1_max = mesh_2d%eta1_max

      delta_x1 = mesh_2d%delta_eta1
      delta_x2 = mesh_2d%delta_eta2

      SLL_ALLOCATE(int_r(Nc_x2), ierr)
      SLL_ALLOCATE(data(Nc_x1 + 1), ierr)
      call sll_s_fft_init_r2r_1d(pfwd, Nc_x2, int_r, int_r, sll_p_fft_forward, normalized=.TRUE.)

      w = 0.0_f64
      l1 = 0.0_f64
      l2 = 0.0_f64
      e = 0.0_f64
      int_r = 0.0_f64

      do i2 = 1, Nc_x2
         do i1 = 1, Nc_x1 + 1
            x1 = x1_min + real(i1 - 1, f64)*delta_x1
            data(i1) = x1*f(i1, i2)
         end do
         w = w + sll_f_compute_integral_trapezoid_1d(data, Nc_x1 + 1, delta_x1)

         do i1 = 1, Nc_x1 + 1
            x1 = x1_min + real(i1 - 1, f64)*delta_x1
            data(i1) = x1*abs(f(i1, i2))
         end do
         l1 = l1 + sll_f_compute_integral_trapezoid_1d(data, Nc_x1 + 1, delta_x1)

         do i1 = 1, Nc_x1 + 1
            x1 = x1_min + real(i1 - 1, f64)*delta_x1
            data(i1) = x1*(f(i1, i2))**2
         end do
         l2 = l2 + sll_f_compute_integral_trapezoid_1d(data, Nc_x1 + 1, delta_x1)

         do i1 = 1, Nc_x1 + 1
            x1 = x1_min + real(i1 - 1, f64)*delta_x1
            data(i1) = x1*((x1*A2(i1, i2))**2 + A1(i1, i2)**2)
         end do
         e = e + sll_f_compute_integral_trapezoid_1d(data, Nc_x1 + 1, delta_x1)

         do i1 = 1, Nc_x1 + 1
            x1 = x1_min + real(i1 - 1, f64)*delta_x1
            data(i1) = x1*phi(i1, i2)
         end do
         int_r(i2) = sll_f_compute_integral_trapezoid_1d(data, Nc_x1 + 1, delta_x1)
      end do

      w = w*delta_x2
      l1 = l1*delta_x2
      l2 = sqrt(l2*delta_x2)
      e = 0.5_f64*e*delta_x2
      call sll_s_fft_exec_r2r_1d(pfwd, int_r, int_r)
      do i1 = 1, 8
         !mode_slope(i1) = time_mode(i1)
         time_mode(i1) = abs(sll_f_fft_get_mode_r2c_1d(pfwd, int_r, i1 - 1))
         !mode_slope(i1) = &
         !  (log(0*time_mode(i1)+1.e-40_f64)-log(0*mode_slope(i1)+1.e-40_f64))/(dt+1.e-40_f64)
      end do

      write (file_id, *) &
         dt*real(step, f64), &
         w, &
         l1, &
         l2, &
         e, &
         maxval(abs(phi(1:Nc_x1 + 1, 1:Nc_x2 + 1))), &
         time_mode(1:8)!,mode_slope

      call sll_s_fft_free(pfwd)

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
   subroutine plot_f_polar(iplot, f, mesh_2d)
      sll_int32, intent(in) :: iplot
      sll_real64, dimension(:, :), intent(in) :: f
      type(sll_t_cartesian_mesh_2d), pointer :: mesh_2d

      sll_int32 :: file_id
      type(sll_t_hdf5_ser_handle) :: hfile_id
      sll_int32 :: error
      sll_real64, dimension(:, :), allocatable :: x1
      sll_real64, dimension(:, :), allocatable :: x2
      sll_int32 :: i, j
      character(len=4)      :: cplot
      sll_int32             :: nnodes_x1, nnodes_x2
      sll_real64 :: r
      sll_real64 :: theta
      sll_real64 :: rmin
      sll_real64 :: rmax
      sll_real64 :: dr
      sll_real64 :: dtheta

      nnodes_x1 = mesh_2d%num_cells1 + 1
      nnodes_x2 = mesh_2d%num_cells2 + 1
      rmin = mesh_2d%eta1_min
      rmax = mesh_2d%eta1_max
      dr = mesh_2d%delta_eta1
      dtheta = mesh_2d%delta_eta2

      !print *,'#maxf=',iplot,maxval(f),minval(f)

      if (iplot == 1) then

         SLL_ALLOCATE(x1(nnodes_x1, nnodes_x2), error)
         SLL_ALLOCATE(x2(nnodes_x1, nnodes_x2), error)
         do j = 1, nnodes_x2
            do i = 1, nnodes_x1
               r = rmin + real(i - 1, f64)*dr
               theta = real(j - 1, f64)*dtheta
               x1(i, j) = r*cos(theta)
               x2(i, j) = r*sin(theta)
            end do
         end do
         call sll_s_hdf5_ser_file_create("polar_mesh-x1.h5", hfile_id, error)
         call sll_o_hdf5_ser_write_array(hfile_id, x1, "/x1", error)
         call sll_s_hdf5_ser_file_close(hfile_id, error)
         call sll_s_hdf5_ser_file_create("polar_mesh-x2.h5", hfile_id, error)
         call sll_o_hdf5_ser_write_array(hfile_id, x2, "/x2", error)
         call sll_s_hdf5_ser_file_close(hfile_id, error)
         deallocate (x1)
         deallocate (x2)

      end if

      call sll_s_int2string(iplot, cplot)
      call sll_o_xdmf_open("f"//cplot//".xmf", "polar_mesh", &
                           nnodes_x1, nnodes_x2, file_id, error)
      call sll_o_xdmf_write_array("f"//cplot, f, "values", &
                                  error, file_id, "Node")
      call sll_s_xdmf_close(file_id, error)
   end subroutine plot_f_polar

#endif

end module sll_m_sim_bsl_gc_2d0v_polar_one_mu
