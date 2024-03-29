!------------------------------------------------------------------------------
! SELALIB
!------------------------------------------------------------------------------
! MODULE: simulation_2d_analytic_field_cartesian
!
! DESCRIPTION:
!> @author  Michel Mehrenberger (mehrenbe@math.unistra.fr)
!> @brief   Swirling deformation flow test - uniform Cartesian grid
!> @details 2D advection with analytical flow field, time dependent, which
!>          causes a periodic swirling deformation (see [REF]). The Cartesian
!>          grid is obtained as a tensor product of two uniform 1D grids.
!------------------------------------------------------------------------------
module sll_m_sim_bsl_ad_2d0v_cart

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

   use sll_m_advection_1d_base, only: &
      sll_c_advector_1d

   use sll_m_advection_1d_bsl, only: &
      sll_f_new_advector_1d_bsl

   use sll_m_advection_1d_csl, only: &
      sll_f_new_csl_1d_advector

   use sll_m_advection_1d_psm, only: &
      sll_f_new_psm_1d_advector

   use sll_m_advection_2d_base, only: &
      sll_c_advector_2d

   use sll_m_advection_2d_bsl, only: &
      sll_f_new_advector_2d_bsl

   use sll_m_advection_2d_tensor_product, only: &
      sll_f_new_advector_2d_tensor_product

   use sll_m_boundary_condition_descriptors, only: &
      sll_p_periodic

   use sll_m_cartesian_meshes, only: &
      sll_o_get_node_positions, &
      sll_f_new_cartesian_mesh_1d, &
      sll_t_cartesian_mesh_1d, &
      sll_t_cartesian_mesh_2d, &
      operator(*)

   use sll_m_characteristics_1d_base, only: &
      sll_c_characteristics_1d_base

   use sll_m_characteristics_1d_explicit_euler, only: &
      sll_f_new_charac_1d_explicit_euler

   use sll_m_characteristics_1d_explicit_euler_conservative, only: &
      sll_f_new_explicit_euler_conservative_1d_charac

   use sll_m_characteristics_1d_trapezoid, only: &
      sll_f_new_trapezoid_1d_charac

   use sll_m_characteristics_1d_trapezoid_conservative, only: &
      sll_f_new_trapezoid_conservative_1d_charac

   use sll_m_characteristics_2d_base, only: &
      sll_c_characteristics_2d_base

   use sll_m_characteristics_2d_explicit_euler, only: &
      sll_f_new_explicit_euler_2d_charac

   use sll_m_characteristics_2d_verlet, only: &
      sll_f_new_verlet_2d_charac

   use sll_m_common_array_initializers, only: &
      sll_f_cos_bell_initializer_2d, &
      sll_f_gaussian_initializer_2d, &
      sll_f_khp1_2d, &
      sll_i_scalar_initializer_1d, &
      sll_i_scalar_initializer_2d, &
      sll_f_sdf_a1_initializer_2d, &
      sll_f_sdf_a2_initializer_2d, &
      sll_f_sdf_time_initializer_1d

   use sll_m_constants, only: &
      sll_p_pi

   use sll_m_cubic_spline_interpolator_1d, only: &
      sll_f_new_cubic_spline_interpolator_1d

   use sll_m_cubic_spline_interpolator_2d, only: &
      sll_t_cubic_spline_interpolator_2d

   use sll_m_gnuplot, only: &
      sll_o_gnuplot_1d

   use sll_m_hdf5_io_serial, only: &
      sll_t_hdf5_ser_handle, &
      sll_s_hdf5_ser_file_create, &
      sll_s_hdf5_ser_file_close, &
      sll_o_hdf5_ser_write_array

   use sll_m_interpolators_1d_base, only: &
      sll_c_interpolator_1d

   use sll_m_interpolators_2d_base, only: &
      sll_c_interpolator_2d

   use sll_m_sim_base, only: &
      sll_c_simulation_base_class

   use sll_m_utilities, only: &
      sll_s_int2string

   use sll_m_xdmf, only: &
      sll_s_xdmf_close, &
      sll_o_xdmf_open, &
      sll_o_xdmf_write_array

   implicit none

   public :: &
      sll_f_new_analytic_field_2d_cartesian

   private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   sll_int32, parameter :: SLL_EULER = 0
   sll_int32, parameter :: SLL_PREDICTOR_CORRECTOR = 1

   !============================================================================
   !> New simulation class (derived from abstract simulation)
   !============================================================================
   type, extends(sll_c_simulation_base_class) :: &
      sll_simulation_2d_analytic_field_cartesian

      !>@name Geometry
      type(sll_t_cartesian_mesh_2d), pointer :: mesh_2d

      !>@name Initial Conditions
      procedure(sll_i_scalar_initializer_2d), nopass, pointer :: init_func
      sll_real64, dimension(:), pointer :: params

      !>@name Advector
      class(sll_c_advector_2d), pointer    :: advect_2d
      procedure(sll_i_scalar_initializer_2d), nopass, pointer :: A1_func
      procedure(sll_i_scalar_initializer_2d), nopass, pointer :: A2_func
      sll_real64, dimension(:), pointer :: A_func_params
      procedure(sll_i_scalar_initializer_1d), nopass, pointer :: A_time_func
      sll_real64, dimension(:), pointer :: A_time_func_params

      !>@name Time Iterations
      sll_real64 :: dt
      sll_int32  :: num_iterations
      sll_int32  :: freq_diag
      sll_int32  :: freq_diag_time

      !>@name Time Loop
      sll_int32 :: time_loop_case
      !>@}

   contains
      procedure, pass(sim) :: run => run_af2d_cartesian   !< Run simulation
      procedure, pass(sim) :: init_from_file => init_fake !< Initialize from file

   end type sll_simulation_2d_analytic_field_cartesian
   !============================================================================

contains

   !----------------------------------------------------------------------------
   !> @brief     Factory function: creates simulation and returns pointer to it
   !> @details   (WIP)
   !> @param[in] filename   Name of the (optional) input file
   !> @return    Pointer to simulation object of proper derived type
   !----------------------------------------------------------------------------
   function sll_f_new_analytic_field_2d_cartesian(filename) result(sim)
      character(len=*), intent(in), optional :: filename
      type(sll_simulation_2d_analytic_field_cartesian), pointer :: sim
      sll_int32 :: ierr

      SLL_ALLOCATE(sim, ierr)

      call initialize_analytic_field_2d_cartesian(sim, filename)

   end function sll_f_new_analytic_field_2d_cartesian

   !----------------------------------------------------------------------------
   !> @brief     Initialize simulation according to defaults or input file
   !> @details   (WIP)
   !> @param[in,out]  sim   Simulation object to be initialized
   !> @param[in] filename   Name of the (optional) input file
   !----------------------------------------------------------------------------
   subroutine initialize_analytic_field_2d_cartesian(sim, filename)
      class(sll_simulation_2d_analytic_field_cartesian), &
         intent(inout)        :: sim
      character(len=*), intent(in), optional :: filename

      sll_int32             :: IO_stat
      sll_int32, parameter  :: input_file = 99

      !geometry
      character(len=256) :: mesh_case_x1
      sll_int32 :: num_cells_x1
      sll_real64 :: x1_min
      sll_real64 :: x1_max
      character(len=256) :: mesh_case_x2
      sll_int32 :: num_cells_x2
      sll_real64 :: x2_min
      sll_real64 :: x2_max

      !initial_function
      character(len=256) :: initial_function_case
      sll_real64 :: kmode_x1
      sll_real64 :: kmode_x2
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
      character(len=256) :: charac2d_case
      character(len=256) :: A_interp_case
      character(len=256) :: advection_field_case
      character(len=256) :: advect1d_x1_case
      character(len=256) :: advect1d_x2_case
      character(len=256) :: charac1d_x1_case
      character(len=256) :: charac1d_x2_case
      character(len=256) :: f_interp1d_x1_case
      character(len=256) :: f_interp1d_x2_case

      sll_real64 :: time_period

      !local variables
      sll_int32 :: Nc_x1
      sll_int32 :: Nc_x2
      type(sll_t_cartesian_mesh_1d), pointer :: mesh_x1
      type(sll_t_cartesian_mesh_1d), pointer :: mesh_x2
      class(sll_c_characteristics_2d_base), pointer :: charac2d
      class(sll_c_characteristics_1d_base), pointer :: charac1d_x1
      class(sll_c_characteristics_1d_base), pointer :: charac1d_x2
      class(sll_c_interpolator_2d), pointer :: f_interp2d
      class(sll_c_interpolator_2d), pointer   :: A1_interp2d
      class(sll_c_interpolator_2d), pointer   :: A2_interp2d
      type(sll_t_cubic_spline_interpolator_2d), target :: f_cs2d
      type(sll_t_cubic_spline_interpolator_2d), target :: A1_cs2d
      type(sll_t_cubic_spline_interpolator_2d), target :: A2_cs2d
      class(sll_c_interpolator_1d), pointer   :: A1_interp1d_x1
      class(sll_c_interpolator_1d), pointer   :: A2_interp1d_x1
      class(sll_c_interpolator_1d), pointer   :: A1_interp1d_x2
      class(sll_c_interpolator_1d), pointer   :: A2_interp1d_x2
      class(sll_c_interpolator_1d), pointer :: f_interp1d_x1
      class(sll_c_interpolator_1d), pointer :: f_interp1d_x2
      class(sll_c_advector_1d), pointer    :: advect_1d_x1
      class(sll_c_advector_1d), pointer    :: advect_1d_x2
      sll_int32 :: ierr
      sll_real64 :: x1_min_bis
      sll_real64 :: x1_max_bis
      sll_real64 :: x2_min_bis
      sll_real64 :: x2_max_bis
      sll_int32 :: Nc_x1_bis
      sll_int32 :: Nc_x2_bis

      namelist /geometry/ &
         mesh_case_x1, &
         num_cells_x1, &
         x1_min, &
         x1_max, &
         mesh_case_x2, &
         num_cells_x2, &
         x2_min, &
         x2_max

      namelist /initial_function/ &
         initial_function_case, &
         kmode_x1, &
         kmode_x2, &
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
         charac2d_case, &
         A_interp_case, &
         advection_field_case, &
         charac1d_x1_case, &
         charac1d_x2_case, &
         advect1d_x1_case, &
         advect1d_x2_case, &
         time_period

    !! set default parameters

      !geometry
      mesh_case_x1 = "SLL_CARTESIAN_MESH"
      num_cells_x1 = 32
      x1_min = 0.0_f64
      x1_max = 2._f64*sll_p_pi
      mesh_case_x2 = "SLL_CARTESIAN_MESH"
      num_cells_x2 = 32
      x2_min = 0.0_f64
      x2_max = 2._f64*sll_p_pi

      !initial function
      initial_function_case = "SLL_KHP1"
      kmode_x1 = 0.5_f64
      kmode_x2 = 1._f64
      eps = 0.015_f64
      sigma_1 = 0.70710678118654752440_f64
      sigma_2 = 0.70710678118654752440_f64
      !initial_function_case="SLL_COS_BELL"
      xc_1 = 1.0_f64
      xc_2 = -0.20_f64

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
      !charac2d_case = "SLL_EULER"
      charac2d_case = "SLL_VERLET"
      A_interp_case = "SLL_CUBIC_SPLINES"
      advection_field_case = "SLL_SWIRLING_DEFORMATION_FLOW"

      advect1d_x1_case = "SLL_BSL"
      advect1d_x2_case = "SLL_BSL"
      charac1d_x1_case = "SLL_EULER"
      !charac1d_x1_case = "SLL_TRAPEZOID"
      charac1d_x2_case = "SLL_EULER"
      !charac1d_x2_case = "SLL_TRAPEZOID"
      f_interp1d_x1_case = "SLL_CUBIC_SPLINES"
      f_interp1d_x2_case = "SLL_CUBIC_SPLINES"
      time_period = 1.5_f64

    !! If available, load data from input file (parameters organized
    !! into namelists and read in blocks)
      if (present(filename)) then
         open (unit=input_file, &
               file=trim(filename)//'.nml', &
               IOStat=IO_stat)
         if (IO_stat /= 0) then
            print *, '#initialize_guiding_center_2d_cartesian()'// &
               ' failed to open file ', trim(filename)//'.nml'
            stop
         end if
         print *, '#initialization with filename:'
         print *, '#', trim(filename)//'.nml'
         read (input_file, geometry)
         read (input_file, initial_function)
         read (input_file, time_iterations)
         read (input_file, advector)
         close (input_file)
      else
         print *, '#initialization with default parameters'
      end if

      Nc_x1 = num_cells_x1
      Nc_x2 = num_cells_x2

      sim%dt = dt
      sim%num_iterations = number_iterations
      sim%freq_diag = freq_diag
      sim%freq_diag_time = freq_diag_time

      select case (mesh_case_x1)
      case ("SLL_CARTESIAN_MESH")
         mesh_x1 => sll_f_new_cartesian_mesh_1d(num_cells_x1, eta_min=x1_min, eta_max=x1_max)
      case default
         print *, '#mesh_case_x1', mesh_case_x1, ' not implemented'
         stop
      end select
      select case (mesh_case_x2)
      case ("SLL_CARTESIAN_MESH")
         mesh_x2 => sll_f_new_cartesian_mesh_1d(num_cells_x2, eta_min=x2_min, eta_max=x2_max)
      case default
         print *, '#mesh_case_x2', mesh_case_x2, ' not implemented'
         stop
      end select
      !Tensor product
      sim%mesh_2d => mesh_x1*mesh_x2

      select case (advect1d_x1_case) ! Advection method along x1 coordinate
      case ("SLL_BSL")
         x1_min_bis = x1_min
         x1_max_bis = x1_max
         Nc_x1_bis = Nc_x1
      case ("SLL_PSM")
         x1_min_bis = x1_min
         x1_max_bis = x1_max
         Nc_x1_bis = Nc_x1
      case ("SLL_CSL")
         x1_min_bis = x1_min - 0.5_f64*mesh_x1%delta_eta
         x1_max_bis = x1_max - 0.5_f64*mesh_x1%delta_eta
         Nc_x1_bis = Nc_x1
      case default
         print *, '#bad value of advect1d_x1_case'
         stop
      end select

      select case (advect1d_x2_case) ! Advection method along x2 coordinate
      case ("SLL_BSL")
         x2_min_bis = x2_min
         x2_max_bis = x2_max
         Nc_x2_bis = Nc_x2
      case ("SLL_PSM")
         x2_min_bis = x2_min
         x2_max_bis = x2_max
         Nc_x2_bis = Nc_x2
      case ("SLL_CSL")
         x2_min_bis = x2_min - 0.5_f64*mesh_x2%delta_eta
         x2_max_bis = x2_max - 0.5_f64*mesh_x2%delta_eta
         Nc_x2_bis = Nc_x2
      case default
         print *, '#bad value of advect1d_x2_case'
         stop
      end select

      select case (f_interp2d_case) ! 2D interpolation method for f(x1,x2)
      case ("SLL_CUBIC_SPLINES")
         call f_cs2d%init( &
            Nc_x1 + 1, &
            Nc_x2 + 1, &
            x1_min, &
            x1_max, &
            x2_min, &
            x2_max, &
            sll_p_periodic, &
            sll_p_periodic)
         f_interp2d => f_cs2d
      case default
         print *, '#bad f_interp2d_case', f_interp2d_case
         print *, '#not implemented'
         print *, '#in initialize_guiding_center_2d_polar'
         stop
      end select

      select case (A_interp_case) ! Interp. method for velocity field A
      case ("SLL_CUBIC_SPLINES")
         call A1_cs2d%init( &
            Nc_x1 + 1, &
            Nc_x2 + 1, &
            x1_min, &
            x1_max, &
            x2_min, &
            x2_max, &
            sll_p_periodic, &
            sll_p_periodic)

         A1_interp2d => A1_cs2d

         call A2_cs2d%init( &
            Nc_x1 + 1, &
            Nc_x2 + 1, &
            x1_min, &
            x1_max, &
            x2_min, &
            x2_max, &
            sll_p_periodic, &
            sll_p_periodic)

         A2_interp2d => A2_cs2d

         A1_interp1d_x1 => sll_f_new_cubic_spline_interpolator_1d( &
                           Nc_x1 + 1, &
                           x1_min, &
                           x1_max, &
                           sll_p_periodic)
         A2_interp1d_x1 => sll_f_new_cubic_spline_interpolator_1d( &
                           Nc_x1 + 1, &
                           x1_min, &
                           x1_max, &
                           sll_p_periodic)
         A1_interp1d_x2 => sll_f_new_cubic_spline_interpolator_1d( &
                           Nc_x2 + 1, &
                           x2_min, &
                           x2_max, &
                           sll_p_periodic)
         A2_interp1d_x2 => sll_f_new_cubic_spline_interpolator_1d( &
                           Nc_x2 + 1, &
                           x2_min, &
                           x2_max, &
                           sll_p_periodic)
      case default
         print *, '#bad A_interp_case', A_interp_case
         print *, '#not implemented'
         print *, '#in initialize_guiding_center_2d_polar'
         stop
      end select

      select case (f_interp1d_x1_case) ! 1D interp. along x1 for f(x1,x2)
      case ("SLL_CUBIC_SPLINES")
         f_interp1d_x1 => sll_f_new_cubic_spline_interpolator_1d( &
                          Nc_x1_bis + 1, &
                          x1_min_bis, &
                          x1_max_bis, &
                          sll_p_periodic)
      case default
         print *, '#bad f_interp1d_x1_case', f_interp1d_x1_case
         print *, '#not implemented'
         print *, '#in initialize_guiding_center_2d_polar'
         stop
      end select

      select case (f_interp1d_x2_case) ! 1D interp. along x2 for f(x1,x2)
      case ("SLL_CUBIC_SPLINES")
         f_interp1d_x2 => sll_f_new_cubic_spline_interpolator_1d( &
                          Nc_x2_bis + 1, &
                          x2_min_bis, &
                          x2_max_bis, &
                          sll_p_periodic)
      case default
         print *, '#bad f_interp1d_x2_case', f_interp1d_x2_case
         print *, '#not implemented'
         print *, '#in initialize_guiding_center_2d_polar'
         stop
      end select

      select case (charac1d_x1_case)
      case ("SLL_EULER")
         charac1d_x1 => sll_f_new_charac_1d_explicit_euler( &
                        Nc_x1_bis + 1, &
                        eta_min=x1_min_bis, &
                        eta_max=x1_max_bis, &
                        bc_type=sll_p_periodic)
      case ("SLL_TRAPEZOID")
         charac1d_x1 => &
            sll_f_new_trapezoid_1d_charac( &
            Nc_x1_bis + 1, &
            A1_interp1d_x1, &
            bc_type=sll_p_periodic, &
            eta_min=x1_min_bis, &
            eta_max=x1_max_bis)
      case ("SLL_EULER_CONSERVATIVE")
         charac1d_x1 => sll_f_new_explicit_euler_conservative_1d_charac( &
                        Nc_x1_bis + 1, &
                        eta_min=x1_min_bis, &
                        eta_max=x1_max_bis, &
                        bc_type=sll_p_periodic)
      case ("SLL_TRAPEZOID_CONSERVATIVE")
         charac1d_x1 => &
            sll_f_new_trapezoid_conservative_1d_charac( &
            Nc_x1_bis + 1, &
            A1_interp1d_x1, &
            bc_type=sll_p_periodic, &
            eta_min=x1_min_bis, &
            eta_max=x1_max_bis)
      case default
         print *, '#bad charac1d_x1_case', charac1d_x1_case
         print *, '#not implemented'
         print *, '#in initialize_guiding_center_2d_polar'
         stop
      end select

      select case (charac1d_x2_case)
      case ("SLL_EULER")
         charac1d_x2 => sll_f_new_charac_1d_explicit_euler( &
                        Nc_x2_bis + 1, &
                        eta_min=x2_min_bis, &
                        eta_max=x2_max_bis, &
                        bc_type=sll_p_periodic)
      case ("SLL_TRAPEZOID")
         charac1d_x2 => &
            sll_f_new_trapezoid_1d_charac( &
            Nc_x2_bis + 1, &
            A2_interp1d_x2, &
            bc_type=sll_p_periodic, &
            eta_min=x2_min_bis, &
            eta_max=x2_max_bis)
      case ("SLL_EULER_CONSERVATIVE")
         charac1d_x2 => sll_f_new_explicit_euler_conservative_1d_charac( &
                        Nc_x2_bis + 1, &
                        eta_min=x2_min_bis, &
                        eta_max=x2_max_bis, &
                        bc_type=sll_p_periodic)
      case ("SLL_TRAPEZOID_CONSERVATIVE")
         charac1d_x2 => &
            sll_f_new_trapezoid_conservative_1d_charac( &
            Nc_x2_bis + 1, &
            A2_interp1d_x2, &
            bc_type=sll_p_periodic, &
            eta_min=x2_min_bis, &
            eta_max=x2_max_bis)
      case default
         print *, '#bad charac1d_x2_case', charac1d_x2_case
         print *, '#not implemented'
         print *, '#in initialize_guiding_center_2d_polar'
         stop
      end select

      select case (advect1d_x1_case)
      case ("SLL_BSL")
         advect_1d_x1 => sll_f_new_advector_1d_bsl( &
                         f_interp1d_x1, &
                         charac1d_x1, &
                         Nc_x1_bis + 1, &
                         eta_min=x1_min_bis, &
                         eta_max=x1_max_bis)
      case ("SLL_CSL")
         advect_1d_x1 => sll_f_new_csl_1d_advector( &
                         f_interp1d_x1, &
                         charac1d_x1, &
                         Nc_x1_bis + 1, &
                         eta_min=x1_min_bis, &
                         eta_max=x1_max_bis, &
                         bc_type=sll_p_periodic)
      case ("SLL_PSM")
         advect_1d_x1 => sll_f_new_psm_1d_advector( &
                         Nc_x1 + 1, &
                         eta_min=x1_min, &
                         eta_max=x1_max)
      case default
         print *, '#bad advect_case', advect1d_x1_case
         print *, '#not implemented'
         print *, '#in initialize_guiding_center_2d_polar'
         stop
      end select

      select case (advect1d_x2_case)
      case ("SLL_BSL")
         advect_1d_x2 => sll_f_new_advector_1d_bsl( &
                         f_interp1d_x2, &
                         charac1d_x2, &
                         Nc_x2_bis + 1, &
                         eta_min=x2_min_bis, &
                         eta_max=x2_max_bis)
      case ("SLL_CSL")
         advect_1d_x2 => sll_f_new_csl_1d_advector( &
                         f_interp1d_x2, &
                         charac1d_x2, &
                         Nc_x2_bis + 1, &
                         eta_min=x2_min_bis, &
                         eta_max=x2_max_bis, &
                         bc_type=sll_p_periodic)
      case ("SLL_PSM")
         advect_1d_x2 => sll_f_new_psm_1d_advector( &
                         Nc_x2 + 1, &
                         eta_min=x2_min, &
                         eta_max=x2_max)
      case default
         print *, '#bad advect_case', advect1d_x2_case
         print *, '#not implemented'
         print *, '#in initialize_guiding_center_2d_polar'
         stop
      end select

      select case (charac2d_case)
      case ("SLL_EULER")
         charac2d => sll_f_new_explicit_euler_2d_charac( &
           Nc_x1+1, &
           Nc_x2+1, &
           eta1_min=x1_min, &
           eta1_max=x1_max, &
           eta2_min=x2_min, &
           eta2_max=x2_max, &
           bc_type_1=sll_p_periodic, &!&sll_p_set_to_limit, &
           bc_type_2=sll_p_periodic)
      case ("SLL_VERLET")
         charac2d => sll_f_new_verlet_2d_charac( &
           Nc_x1+1, &
           Nc_x2+1, &
           A1_interp2d, &
           A2_interp2d, &
           A1_interp1d_x1, &
           A2_interp1d_x1, &
           bc_type_1=sll_p_periodic, &!&sll_p_set_to_limit, &
           bc_type_2=sll_p_periodic, &
           eta1_min=x1_min, &
           eta1_max=x1_max, &
           eta2_min=x2_min, &
           eta2_max=x2_max )
      case default
         print *, '#bad charac2d_case', charac2d_case
         print *, '#not implemented'
         print *, '#in initialize_guiding_center_2d_polar'
         stop
      end select

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
      case ("SLL_TENSOR_PRODUCT")
         sim%advect_2d => sll_f_new_advector_2d_tensor_product( &
                          advect_1d_x1, &
                          advect_1d_x2, &
                          Nc_x1 + 1, &
                          Nc_x2 + 1)
      case default
         print *, '#bad advect_case', advect2d_case
         print *, '#not implemented'
         print *, '#in initialize_guiding_center_2d_polar'
         stop
      end select

      select case (advection_field_case)
      case ("SLL_SWIRLING_DEFORMATION_FLOW")
         sim%A1_func => sll_f_sdf_a1_initializer_2d
         sim%A2_func => sll_f_sdf_a2_initializer_2d
         SLL_ALLOCATE(sim%A_func_params(2), ierr)
         sim%A_time_func => sll_f_sdf_time_initializer_1d
         SLL_ALLOCATE(sim%A_time_func_params(1), ierr)
         sim%A_time_func_params(1) = time_period
      case default
         print *, '#bad advect_case', advect2d_case
         print *, '#not implemented'
         print *, '#in initialize_guiding_center_2d_polar'
         stop
      end select

      select case (initial_function_case) ! Initial conditions
      case ("SLL_KHP1")
         sim%init_func => sll_f_khp1_2d
         SLL_ALLOCATE(sim%params(3), ierr)
         sim%params(1) = eps
         sim%params(2) = kmode_x1
         sim%params(3) = kmode_x2
      case ("SLL_GAUSSIAN")
         sim%init_func => sll_f_gaussian_initializer_2d
         SLL_ALLOCATE(sim%params(4), ierr)
         sim%params(1) = xc_1
         sim%params(2) = xc_2
         sim%params(3) = sigma_1
         sim%params(4) = sigma_2
      case ("SLL_COS_BELL")
         sim%init_func => sll_f_cos_bell_initializer_2d
         SLL_ALLOCATE(sim%params(2), ierr)
         sim%params(1) = xc_1
         sim%params(2) = xc_2
      case default
         print *, '#bad initial_function_case', initial_function_case
         print *, '#not implemented'
         print *, '#in initialize_analytic_field_2d_cartesian'
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

   end subroutine initialize_analytic_field_2d_cartesian

   !----------------------------------------------------------------------------
   !> @brief         Fake initialization of simulation (do nothing)
   !> @details       (WIP)
   !> @param[in,out] sim       Simulation object
   !> @param[in]     filename  Name of input file
   !----------------------------------------------------------------------------
   subroutine init_fake(sim, filename)
      class(sll_simulation_2d_analytic_field_cartesian), &
         intent(inout) :: sim
      character(len=*), intent(in)    :: filename

      print *, '# Do not use the routine init_vp4d_fake'
      print *, '#use instead initialize_vlasov_par_poisson_seq_cart'
      print *, '#filename=', filename
      print *, sim%dt
      stop

   end subroutine init_fake

   !----------------------------------------------------------------------------
   !> @brief         Run simulation
   !> @details       (WIP)
   !> @param[in,out] sim  Simulation object
   !----------------------------------------------------------------------------
   subroutine run_af2d_cartesian(sim)
      ! Argument declaration
      class(sll_simulation_2d_analytic_field_cartesian), intent(inout) :: sim

      ! Declarations: arrays allocated within this subroutine
      sll_real64, dimension(:, :), allocatable :: f, f_old, f_init
      sll_real64, dimension(:, :), allocatable :: A1, A1_init !advection fields
      sll_real64, dimension(:, :), allocatable :: A2, A2_init !advection fields
      sll_real64, dimension(:), allocatable :: f_visu_buf1d
      !
      ! Declarations: objects allocated by Selalib library functions
      sll_real64, dimension(:), allocatable :: node_positions_x1
      sll_real64, dimension(:), allocatable :: node_positions_x2
      type(sll_t_cartesian_mesh_1d), pointer :: mesh_x1
      type(sll_t_cartesian_mesh_1d), pointer :: mesh_x2
      !
      ! Declarations: mesh-related scalars
      sll_int32  :: Nc_x1, Nc_x2
      sll_real64 :: delta_x1, delta_x2
      sll_real64 :: x1_min, x1_max, x2_min, x2_max
      sll_real64 :: x1, x2
      sll_int32  :: i1, i2, i
      !
      ! Declarations: time-related scalars
      sll_int32  :: step, nb_step
      sll_real64 :: dt, time_factor
      sll_int32  :: iplot
      !
      ! Declarations: error codes
      sll_int32 :: ierr      ! allocation
      sll_int32 :: IO_stat   ! input/output operations
      !
      ! Declarations: other scalars
      sll_int32 :: diag_id = 88

      ! Preprocessing
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

      ! Create 1D cartesian meshes along x1 and x2 directions, then extract nodes
      mesh_x1 => sll_f_new_cartesian_mesh_1d(Nc_x1, eta_min=x1_min, eta_max=x1_max)
      mesh_x2 => sll_f_new_cartesian_mesh_1d(Nc_x2, eta_min=x2_min, eta_max=x2_max)
      call sll_o_get_node_positions(mesh_x1, node_positions_x1)
      call sll_o_get_node_positions(mesh_x2, node_positions_x2)

      ! Allocations
      SLL_ALLOCATE(f(Nc_x1 + 1, Nc_x2 + 1), ierr)
      SLL_ALLOCATE(f_old(Nc_x1 + 1, Nc_x2 + 1), ierr)
      SLL_ALLOCATE(f_init(Nc_x1 + 1, Nc_x2 + 1), ierr)
      SLL_ALLOCATE(A1(Nc_x1 + 1, Nc_x2 + 1), ierr)
      SLL_ALLOCATE(A2(Nc_x1 + 1, Nc_x2 + 1), ierr)
      SLL_ALLOCATE(A1_init(Nc_x1 + 1, Nc_x2 + 1), ierr)
      SLL_ALLOCATE(A2_init(Nc_x1 + 1, Nc_x2 + 1), ierr)
      SLL_ALLOCATE(f_visu_buf1d(Nc_x2 + 1), ierr)

      ! Initialisation of distribution function
      do i2 = 1, Nc_x2 + 1
         x2 = x2_min + real(i2 - 1, f64)*delta_x2
         do i1 = 1, Nc_x1 + 1
            x1 = x1_min + real(i1 - 1, f64)*delta_x1
            f(i1, i2) = sim%init_func(x1, x2, sim%params)
            f_init(i1, i2) = sim%init_func(x1, x2, sim%params)
            A1_init(i1, i2) = sim%A1_func(x1, x2, sim%A_func_params)
            A2_init(i1, i2) = sim%A2_func(x1, x2, sim%A_func_params)
         end do
      end do

      ! Open diagnostics file (never used)
      open (unit=diag_id, file='thdiag.dat', IOStat=IO_stat)
      if (IO_stat /= 0) then
         print *, '#run_af2d_cartesian (sim) failed to open file thdiag.dat'
         STOP
      end if

      ! Initialize plot index
      iplot = 0

      !++++++++++++++++++++++++++   BEGIN TIME LOOP   +++++++++++++++++++++++++++
      do step = 0, nb_step
         f_old = f

#ifndef NOHDF5
         if (modulo(step, sim%freq_diag) == 0) then
            print *, "#step= ", step
            call plot_f_cartesian(iplot, f, sim%mesh_2d)
            do i = 1, Nc_x2 + 1
               f_visu_buf1d(i) = sum(f(1:Nc_x1, i))*delta_x1
            end do
            call sll_o_gnuplot_1d( &
               f_visu_buf1d(1:Nc_x2 + 1), &
               node_positions_x2(1:Nc_x2 + 1), &
               'intfdx', &
               iplot)
            call sll_o_gnuplot_1d( &
               f_visu_buf1d(1:Nc_x2 + 1), &
               node_positions_x2(1:Nc_x2 + 1), &
               'intfdx')
            iplot = iplot + 1
         end if
#endif

         select case (sim%time_loop_case)
         case (SLL_EULER)
            time_factor = sim%A_time_func( &
                          real(step, f64)*sim%dt, &
                          sim%A_time_func_params)
            A1 = time_factor*A1_init
            A2 = time_factor*A2_init
            call sim%advect_2d%advect_2d(A1, A2, sim%dt, f_old, f)
         case (SLL_PREDICTOR_CORRECTOR)
            time_factor = sim%A_time_func( &
                          (real(step, f64) + 0.5_f64)*sim%dt, &
                          sim%A_time_func_params)
            A1 = time_factor*A1_init
            A2 = time_factor*A2_init
            call sim%advect_2d%advect_2d(A1, A2, sim%dt, f_old, f)
         case default
            print *, '#bad time_loop_case', sim%time_loop_case
            print *, '#not implemented'
            print *, '#in run_af2d_cartesian'
            print *, '#available options are:'
            print *, '#SLL_EULER=', SLL_EULER
            print *, '#SLL_PREDICTOR_CORRECTOR=', SLL_PREDICTOR_CORRECTOR
         end select

      end do
      !+++++++++++++++++++++++++++   END TIME LOOP   ++++++++++++++++++++++++++++

      close (diag_id)

      print *, maxval(abs(f - f_init))
      print *, '#run_af2d_cartesian PASSED'
   end subroutine run_af2d_cartesian

!*********************
#ifndef NOHDF5
!*********************

   !----------------------------------------------------------------------------
   !> @brief     Save the mesh structure
   !> @details   (WIP)
   !> @param[in] iplot    Plot index
   !> @param[in] f        Distribution function (2D array)
   !> @param[in] mesh_2d  2D cartesian mesh object
   !----------------------------------------------------------------------------
   subroutine plot_f_cartesian(iplot, f, mesh_2d)

      ! Function arguments
      sll_int32, intent(in) :: iplot
      sll_real64, dimension(:, :), intent(in) :: f
      type(sll_t_cartesian_mesh_2d), intent(in) :: mesh_2d

      ! Local variables declarations
      sll_int32 :: file_id
      type(sll_t_hdf5_ser_handle) :: hfile_id
      sll_int32 :: error
      sll_real64, dimension(:, :), allocatable :: x1
      sll_real64, dimension(:, :), allocatable :: x2
      sll_int32 :: i, j
      character(len=4)      :: cplot
      sll_int32             :: nnodes_x1, nnodes_x2
      sll_real64 ::  x1_min, x2_min
      sll_real64 ::  x1_max, x2_max
      sll_real64 :: dx1
      sll_real64 :: dx2

      nnodes_x1 = mesh_2d%num_cells1 + 1
      nnodes_x2 = mesh_2d%num_cells2 + 1
      x1_min = mesh_2d%eta1_min
      x1_max = mesh_2d%eta1_max
      x2_min = mesh_2d%eta2_min
      x2_max = mesh_2d%eta2_max
      dx1 = mesh_2d%delta_eta1
      dx2 = mesh_2d%delta_eta2

      !print *,'#maxf=',iplot,maxval(f),minval(f)

      if (iplot == 1) then

         SLL_ALLOCATE(x1(nnodes_x1, nnodes_x2), error)
         SLL_ALLOCATE(x2(nnodes_x1, nnodes_x2), error)
         do j = 1, nnodes_x2
            do i = 1, nnodes_x1
               x1(i, j) = x1_min + real(i - 1, f64)*dx1
               x2(i, j) = x2_min + real(j - 1, f64)*dx2
            end do
         end do
         call sll_s_hdf5_ser_file_create("cartesian_mesh-x1.h5", hfile_id, error)
         call sll_o_hdf5_ser_write_array(hfile_id, x1, "/x1", error)
         call sll_s_hdf5_ser_file_close(hfile_id, error)
         call sll_s_hdf5_ser_file_create("cartesian_mesh-x2.h5", hfile_id, error)
         call sll_o_hdf5_ser_write_array(hfile_id, x2, "/x2", error)
         call sll_s_hdf5_ser_file_close(hfile_id, error)
         deallocate (x1)
         deallocate (x2)

      end if

      call sll_s_int2string(iplot, cplot)
      call sll_o_xdmf_open("f"//cplot//".xmf", "cartesian_mesh", &
                           nnodes_x1, nnodes_x2, file_id, error)
      call sll_o_xdmf_write_array("f"//cplot, f, "values", &
                                  error, file_id, "Node")
      call sll_s_xdmf_close(file_id, error)
   end subroutine plot_f_cartesian

!*********************
#endif
!*********************

end module sll_m_sim_bsl_ad_2d0v_cart
