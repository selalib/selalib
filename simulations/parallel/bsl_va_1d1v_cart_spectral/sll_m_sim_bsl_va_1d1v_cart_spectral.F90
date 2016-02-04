!**************************************************************
!  Copyright INRIA
!  Authors : 
!     CALVI project team
!  
!  This code SeLaLib (for Semi-Lagrangian-Library) 
!  is a parallel library for simulating the plasma turbulence 
!  in a tokamak.
!  
!  This software is governed by the CeCILL-B license 
!  under French law and abiding by the rules of distribution 
!  of free software.  You can  use, modify and redistribute 
!  the software under the terms of the CeCILL-B license as 
!  circulated by CEA, CNRS and INRIA at the following URL
!  "http://www.cecill.info". 
!**************************************************************
!
! Vlasov-Poisson simulation in 1Dx1D
! Landau and KEEN waves test cases
!
! contact: Michel Mehrenberger (mehrenbe@math.unistra.fr)
!
! current investigations:
!   High order splitting in time
!   KEEN waves with uniform and non uniform grid in velocity

module sll_m_sim_bsl_va_1d1v_cart_spectral

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_errors.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_advection_1d_ampere, only: &
    sll_t_ampere_1d_advector_ptr, &
    sll_f_new_ampere_1d_advector

  use sll_m_advection_1d_base, only: &
    sll_t_advection_1d_base_ptr

  use sll_m_advection_1d_non_uniform_cubic_splines, only: &
    sll_f_new_non_uniform_cubic_splines_1d_advector

  use sll_m_advection_1d_periodic, only: &
    sll_f_new_periodic_1d_advector

  use sll_m_advection_1d_spectral, only: &
    sll_f_new_spectral_1d_advector

  use sll_m_ascii_io, only: &
    sll_s_ascii_file_close, &
    sll_s_ascii_file_create

  use sll_m_binary_io, only: &
    sll_s_binary_file_close, &
    sll_s_binary_file_create, &
    sll_s_binary_read_array_0d, &
    sll_s_binary_read_array_2d, &
    sll_s_binary_write_array_0d, &
    sll_s_binary_write_array_1d, &
    sll_s_binary_write_array_2d

  use sll_m_buffer_loader_utilities, only: &
    sll_s_compute_displacements_array_2d, &
    sll_s_load_buffer_2d, &
    sll_f_receive_counts_array_2d

  use sll_m_cartesian_meshes, only: &
    sll_o_get_node_positions, &
    sll_f_new_cartesian_mesh_1d, &
    sll_t_cartesian_mesh_1d, &
    sll_t_cartesian_mesh_2d, &
    operator(*)

  use sll_m_collective, only: &
    sll_o_collective_allreduce, &
    sll_s_collective_gatherv_real64, &
    sll_f_get_collective_rank, &
    sll_f_get_collective_size, &
    sll_v_world_collective

  use sll_m_common_array_initializers, only: &
    sll_f_beam_initializer_2d, &
    sll_f_bump_on_tail_initializer_2d, &
    sll_f_landau_initializer_2d, &
    sll_i_scalar_initializer_2d, &
    sll_f_two_stream_instability_initializer_2d

  use sll_m_constants, only: &
    sll_p_pi

  use sll_m_fft, only: &
    sll_s_fft_exec_c2r_1d, &
    sll_s_fft_exec_r2c_1d, &
    sll_s_fft_exec_r2r_1d, &
    sll_p_fft_forward, &
    sll_f_fft_get_mode_r2c_1d, &
    sll_s_fft_init_r2r_1d, &
    sll_t_fft

  use sll_m_gnuplot, only: &
    sll_o_gnuplot_1d

  use sll_m_parallel_array_initializer, only: &
    sll_o_2d_parallel_array_initializer_cartesian

  use sll_m_periodic_interp, only: &
    sll_p_lagrange, &
    sll_p_spline, &
    sll_p_trigo

  use sll_m_poisson_1d_base, only: &
    sll_c_poisson_1d_base

  use sll_m_poisson_1d_periodic, only: &
    sll_f_new_poisson_1d_periodic

  use sll_m_poisson_1d_polar, only: &
    sll_f_new_poisson_1d_polar

  use sll_m_primitives, only: &
    sll_s_function_to_primitive, &
    sll_s_primitive_to_function

  use sll_m_remapper, only: &
    sll_o_apply_remap_2d, &
    sll_o_compute_local_sizes, &
    sll_o_initialize_layout_with_distributed_array, &
    sll_t_layout_2d, &
    sll_o_local_to_global, &
    sll_f_new_layout_2d, &
    sll_o_new_remap_plan, &
    sll_t_remap_plan_2d_real64, &
    sll_o_view_lims

  use sll_m_sim_base, only: &
    sll_c_simulation_base_class

  use sll_m_time_splitting_coeff, only: &
    sll_f_new_time_splitting_coeff, &
    sll_p_lie_tv, &
    sll_p_lie_vt, &
    sll_p_order6_vtv, &
    sll_p_order6vp_tvt, &
    sll_p_order6vp_vtv, &
    sll_p_order6vpnew1_vtv, &
    sll_p_order6vpnew2_vtv, &
    sll_p_order6vpnew_tvt, &
    sll_p_strang_tvt, &
    sll_p_strang_vtv, &
    sll_p_triple_jump_tvt, &
    sll_p_triple_jump_vtv, &
    sll_t_splitting_coeff

  use sll_m_utilities, only: &
    sll_s_compute_bloc, &
    sll_s_compute_mesh_from_bloc, &
    sll_s_int2string, &
    sll_s_pfenvelope, &
    sll_s_new_file_id

  use sll_m_xdmf, only: &
    sll_s_plot_f_cartesian

  use sll_mpi, only: &
    mpi_sum

#ifdef _OPENMP
  use omp_lib, only: &
    omp_get_num_threads, &
    omp_get_thread_num

#endif
  implicit none

  public :: &
    sll_f_new_va2d_par_cart, &
    sll_t_simulation_2d_vlasov_ampere_cart

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  integer, parameter :: SLL_ADVECTIVE    = 0
  integer, parameter :: SLL_CONSERVATIVE = 1

  type, extends(sll_c_simulation_base_class) :: &
       sll_t_simulation_2d_vlasov_ampere_cart

   sll_int32                            :: num_threads
   type(sll_t_cartesian_mesh_2d), pointer :: mesh2d
   sll_int32                            :: num_dof_x2

   sll_real64, dimension(:),   pointer  :: x1_array
   sll_real64, dimension(:),   pointer  :: x2_array
   sll_real64, dimension(:,:), pointer  :: x2_array_omp
   sll_real64, dimension(:),   pointer  :: integration_weight
   sll_int32,  dimension(:),   pointer  :: every_x1
   sll_int32,  dimension(:),   pointer  :: every_x2
   sll_int32                            :: num_bloc_x1
   sll_int32                            :: num_bloc_x2
   sll_int32,  dimension(:),   pointer  :: bloc_index_x1
   sll_int32,  dimension(:),   pointer  :: bloc_index_x2
   
   sll_real64                 :: kx
   sll_real64                 :: eps
   character(len=256)         :: restart_file
   logical                    :: time_init_from_restart_file

   !MM initial function
   procedure(sll_i_scalar_initializer_2d), nopass, pointer :: init_func
   !ES equilibrium function
   procedure(sll_i_scalar_initializer_2d), nopass, pointer :: equil_func
   sll_real64, dimension(:), pointer                     :: equil_func_params
   sll_real64, dimension(:), pointer                     :: params

   sll_real64 :: nrj0
   sll_real64 :: dt
   sll_int32  :: num_iterations 
   sll_int32  :: freq_diag
   sll_int32  :: freq_diag_time
   sll_int32  :: freq_diag_restart
   sll_int32  :: nb_mode
   sll_real64 :: time_init

   type(sll_t_splitting_coeff), pointer :: split
   character(len=256)      :: thdiag_filename
  
   logical    :: driven 
   sll_real64 :: t0
   sll_real64 :: twL
   sll_real64 :: twR
   sll_real64 :: tflat
   sll_real64 :: tL
   sll_real64 :: tR
   sll_real64 :: Edrmax
   sll_real64 :: omegadr
   logical    :: turn_drive_off

   type(sll_t_advection_1d_base_ptr), dimension(:), pointer :: advect_x1 
   type(sll_t_advection_1d_base_ptr), dimension(:), pointer :: advect_x2
   type(sll_t_ampere_1d_advector_ptr),    dimension(:), pointer :: advect_ampere_x1

   sll_int32  :: advection_form_x2
   sll_real64 :: factor_x1
   sll_real64 :: factor_x2_rho
   sll_real64 :: factor_x2_1

   class(sll_c_poisson_1d_base), pointer :: poisson 
   logical :: ampere = .false.

   sll_real64 :: L

   sll_real64, dimension(:), allocatable :: node_positions_x2
           
   contains !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     procedure, pass(sim) :: run => run_va2d_cartesian
     procedure, pass(sim) :: init_from_file => init_va2d_fake !init_va2d_par_cart

  end type sll_t_simulation_2d_vlasov_ampere_cart

  interface sll_o_delete
     module procedure delete_va2d_par_cart
  end interface sll_o_delete


  sll_int32                          :: istep
  sll_int32                          :: iplot
  sll_int32                          :: adr_id
  sll_int32                          :: edr_id
  sll_int32                          :: t_id
  sll_int32                          :: deltaf_id
  sll_int32                          :: rhotot_id
  sll_int32                          :: efield_id
  sll_int32                          :: thdiag_id
  sll_real64                         :: adr
  sll_real64                         :: time

  sll_int32, dimension(:),    allocatable :: collective_displs
  sll_int32, dimension(:),    allocatable :: collective_recvcnts

  sll_real64, dimension(:),   allocatable :: buf_fft
  sll_real64, dimension(:),   allocatable :: x2_array_unit
  sll_real64, dimension(:),   allocatable :: f_x1_buf1d
  sll_real64, dimension(:),   allocatable :: f_visu_buf1d
  sll_real64, dimension(:,:), allocatable :: f_visu 
  sll_real64, dimension(:,:), allocatable :: f1d_omp_in
  sll_real64, dimension(:,:), allocatable :: f1d_omp_out

  type(sll_t_fft)             :: pfwd
  logical                     :: MPI_MASTER

contains

  function sll_f_new_va2d_par_cart( filename, num_run ) result(sim)    

    type(sll_t_simulation_2d_vlasov_ampere_cart), pointer :: sim    
    character(len=*), intent(in), optional               :: filename
    sll_int32, intent(in), optional                      :: num_run

    sll_int32                                            :: ierr

    SLL_ALLOCATE(sim, ierr)
    call init_va2d_par_cart( sim, filename, num_run )
       
  end function sll_f_new_va2d_par_cart

  subroutine change_initial_function_va2d_par_cart( sim,       &
                                                    init_func, &
                                                    params,    &
                                                    num_params)

    class(sll_t_simulation_2d_vlasov_ampere_cart)  :: sim
    procedure(sll_i_scalar_initializer_2d), pointer :: init_func
    sll_real64, dimension(:), pointer             :: params
    sll_int32, intent(in)                         :: num_params

    character(len=*), parameter :: this_sub_name = &
      'change_initial_function_va2d_par_cart'

    sll_int32 :: ierr
    sll_int32 :: i
    
    sim%init_func => init_func
    if (associated(sim%params)) SLL_DEALLOCATE(sim%params,ierr)

    if (num_params<1) then
      SLL_ERROR( this_sub_name, '#num_params should be >=1' )
    endif

    SLL_ALLOCATE(sim%params(num_params),ierr)

    if (size(params)<num_params) then
      SLL_ERROR( this_sub_name, '#size of params is not good' )
    endif

    do i=1,num_params
      sim%params(i) = params(i)
    end do
    
  end subroutine change_initial_function_va2d_par_cart

  !ES enable externally defined equilibrium function  
  subroutine change_equilibrium_function_va2d_par_cart( sim,        &
                                                        equil_func, &
                                                        params,     &
                                                        num_params)

    class(sll_t_simulation_2d_vlasov_ampere_cart)  :: sim
    procedure(sll_i_scalar_initializer_2d), pointer :: equil_func
    sll_real64, dimension(:), pointer             :: params
    sll_int32, intent(in)                         :: num_params

    character(len=*), parameter :: this_sub_name = &
      'change_equilibrium_function_va2d_par_cart'

    sll_int32 :: ierr
    sll_int32 :: i
    
    sim%equil_func => equil_func
    if (associated(sim%equil_func_params)) then
      SLL_DEALLOCATE(sim%equil_func_params,ierr)
    end if

    if (num_params<1) then
      SLL_ERROR( this_sub_name, '#num_params should be >=1' )
    endif

    SLL_ALLOCATE(sim%equil_func_params(num_params),ierr)

    if (size(params)<num_params) then
      SLL_ERROR( this_sub_name, '#size of params is not good' )
    endif

    do i=1,num_params
      sim%equil_func_params(i) = params(i)
    enddo
    
  end subroutine change_equilibrium_function_va2d_par_cart

  subroutine init_va2d_par_cart( sim, filename, num_run )

    class(sll_t_simulation_2d_vlasov_ampere_cart) :: sim
    character(len=*), optional                   :: filename
    sll_int32, intent(in), optional :: num_run

    character(len=*), parameter :: this_sub_name = 'init_va2d_par_cart'
    character(len=291)          :: err_msg

    character(len=256) :: mesh_case_x1
    sll_int32          :: num_cells_x1
    sll_real64         :: x1_min
    sll_real64         :: x1_max
    sll_int32          :: nbox_x1

    character(len=256) :: mesh_case_x2
    sll_int32          :: num_cells_x2
    sll_real64         :: x2_min
    sll_real64         :: x2_max
    sll_real64         :: x2_fine_min
    sll_real64         :: x2_fine_max
    sll_real64         :: density_x2_min_to_x2_fine_min
    sll_real64         :: density_x2_fine_min_to_x2_fine_max
    sll_real64         :: density_x2_fine_max_to_x2_max
    sll_int32          :: every_x2_min_to_x2_fine_min
    sll_int32          :: every_x2_fine_min_to_x2_fine_max
    sll_int32          :: every_x2_fine_max_to_x2_max
    sll_int32          :: every_x1
    sll_int32          :: every_x2
    
    
    !initial_function
    character(len=256) :: initial_function_case
    sll_real64         :: kmode
    sll_real64         :: eps
    sll_real64         :: alpha_gaussian
    character(len=256) :: restart_file
    logical            :: time_init_from_restart_file
    
    !time_iterations
    sll_real64         :: dt
    sll_int32          :: number_iterations
    sll_int32          :: freq_diag
    sll_int32          :: freq_diag_restart
    sll_int32          :: freq_diag_time
    sll_int32          :: nb_mode
    sll_real64         :: time_init
    character(len=256) :: split_case

    !advector
    character(len=256) :: advector_x1
    sll_int32          :: order_x1
    character(len=256) :: advector_x2
    sll_int32          :: order_x2
    character(len=256) :: advection_form_x2
    character(len=256) :: integration_case
    sll_real64         :: factor_x1
    sll_real64         :: factor_x2_rho
    sll_real64         :: factor_x2_1

    character(len=256) :: poisson_solver
    character(len=256) :: ampere_solver
    character(len=256) :: drive_type
    sll_real64         :: keen_t0
    sll_real64         :: keen_tL
    sll_real64         :: keen_tR
    sll_real64         :: keen_twL
    sll_real64         :: keen_twR
    sll_real64         :: keen_tflat
    logical            :: keen_turn_drive_off
    sll_real64         :: keen_Edrmax
    sll_real64         :: keen_omegadr
    
    sll_int32                            :: IO_stat
    sll_int32                            :: input_file
    type(sll_t_cartesian_mesh_1d), pointer :: mesh_x1
    type(sll_t_cartesian_mesh_1d), pointer :: mesh_x2
    sll_int32                            :: ierr
    sll_int32, parameter                 :: param_out = 37
    sll_int32, parameter                 :: param_out_drive = 40
    sll_real64                           :: bloc_coord(2)
    sll_int32                            :: bloc_index(3)
    sll_int32                            :: i
    sll_int32                            :: num_threads
    sll_int32                            :: tid
    character(len=256)      :: str_num_run
    character(len=256)      :: filename_loc


  
    namelist /geometry/                   &
      mesh_case_x1,                       &
      num_cells_x1,                       &
      x1_min,                             &
      x1_max,                             &
      nbox_x1,                            &
      mesh_case_x2,                       &
      num_cells_x2,                       &
      x2_min,                             &
      x2_max,                             &
      x2_fine_min,                        &
      x2_fine_max,                        &
      density_x2_min_to_x2_fine_min,      &
      density_x2_fine_min_to_x2_fine_max, &
      density_x2_fine_max_to_x2_max,      &
      every_x2_min_to_x2_fine_min,        &
      every_x2_fine_min_to_x2_fine_max,   &
      every_x2_fine_max_to_x2_max,        &
      every_x1,                           &
      every_x2

    namelist /initial_function/           &
      initial_function_case,              &
      kmode,                              &
      eps,                                &
      alpha_gaussian,                     &
      restart_file,                       &
      time_init_from_restart_file

    namelist /time_iterations/            &
      dt,                                 &
      number_iterations,                  &
      freq_diag,                          &
      freq_diag_time,                     &
      freq_diag_restart,                  &
      nb_mode,                            &
      time_init,                          &
      split_case

    namelist /advector/                   &
      advector_x1,                        &
      order_x1,                           &
      advector_x2,                        &
      order_x2,                           &
      advection_form_x2,                  &
      factor_x1,                          &
      factor_x2_rho,                      &
      factor_x2_1,                        &
      integration_case

    namelist /poisson/                    &
      ampere_solver,                      &
      poisson_solver


    namelist / drive /                    &
      drive_type,                         &
      keen_t0,                            &
      keen_twL,                           &
      keen_twR,                           &
      keen_tflat,                         &
      keen_tL,                            &
      keen_tR,                            &
      keen_turn_drive_off,                &
      keen_Edrmax,                        &
      keen_omegadr

    num_threads = 1

    if (sll_f_get_collective_rank(sll_v_world_collective)==0) then
      MPI_MASTER = .true.
    else
      MPI_MASTER = .false.
    end if
    
#ifdef _OPENMP
    !$OMP PARALLEL SHARED(num_threads)
    if (omp_get_thread_num()==0) then
      num_threads =  omp_get_num_threads()      
    endif
    !$OMP END PARALLEL
#endif

    sim%num_threads = num_threads
    print *,'#num_threads=',num_threads

    !! set default parameters
    
    !geometry
    mesh_case_x1 = "SLL_LANDAU_MESH"
    num_cells_x1 = 32
    x1_min       = 0.0_f64
    nbox_x1      = 1
    mesh_case_x2 = "SLL_CARTESIAN_MESH"
    num_cells_x2 = 64
    x2_min       = -6._f64
    x2_max       =  6._f64
    
    mesh_case_x2 = "SLL_TWO_GRID_MESH"
    num_cells_x2 = 64
    x2_min       = -6._f64
    x2_max       = 6._f64
    x2_fine_min  = 0.36_f64
    x2_fine_max  = 2.28_f64

    density_x2_min_to_x2_fine_min      = 1.0_f64
    density_x2_fine_min_to_x2_fine_max = 1.0_f64
    density_x2_fine_max_to_x2_max      = 1.0_f64
    every_x2_min_to_x2_fine_min        = 1
    every_x2_fine_min_to_x2_fine_max   = 1
    every_x2_fine_max_to_x2_max        = 1
    every_x1                           = 1
    every_x2                           = 1

        
    !initial_function
    initial_function_case       = "SLL_LANDAU"   !"SLL_BEAM"
    kmode                       = 0.5_f64
    eps                         = 0.001_f64
    alpha_gaussian              = 0.2_f64
    restart_file                = "no_restart_file"
    time_init_from_restart_file = .false.
    
    !time_iterations
    dt                          = 0.1_f64
    number_iterations           = 600
    freq_diag                   = 100
    freq_diag_time              = 1
    freq_diag_restart           = 5000
    nb_mode                     = 5
    time_init                   = 0._f64
    split_case                  = "SLL_STRANG_VTV" 

    !split_case                 = "SLL_STRANG_TVT" 
    !split_case                 = "SLL_ORDER6VPNEW1_VTV" 
    !split_case                 = "SLL_ORDER6VPNEW2_VTV" 
    !split_case                 = "SLL_ORDER6_VTV"
    !split_case                 = "SLL_LIE_TV"

    !advector
    advector_x1                 = "SLL_LAGRANGE"
    order_x1                    = 4
    advector_x2                 = "SLL_LAGRANGE"
    order_x2                    = 4
    advection_form_x2           = "SLL_ADVECTIVE"
    factor_x1                   = 1._f64
    factor_x2_rho               = 1._f64
    factor_x2_1                 = 1._f64

    integration_case = "SLL_TRAPEZOID" !"SLL_RECTANGLE"
    poisson_solver   = "SLL_FFT"
    drive_type       = "SLL_NO_DRIVE"  !"SLL_KEEN_DRIVE"  

    !keen_t0             = 0.
    !keen_tL             = 69.
    !keen_tR             = 307.
    !keen_twL            = 20.
    !keen_twR            = 20.
    !keen_tflat          = 100.
    !keen_turn_drive_off = .true.
    !keen_Edrmax         = 0.2
    !keen_omegadr        = 0.37	

    if(present(num_run))then
      !call sll_s_int2string(num_run, str_num_run)
      write(str_num_run, *) num_run
      str_num_run = adjustl(str_num_run) 
      sim%thdiag_filename = "thdiag_"//trim(str_num_run)//".dat"
    else      
      sim%thdiag_filename = "thdiag.dat"
    endif

    if (present(filename)) then

      filename_loc = filename
      filename_loc = adjustl(filename_loc)
      
      if(present(num_run)) then
        filename_loc = trim(filename)//"_"//trim(str_num_run)
        !filename_loc = adjustl(filename_loc)
        !print *,'filename_loc=',filename_loc
      endif

      call sll_s_new_file_id(input_file, ierr)

      open(unit = input_file, file=trim(filename_loc)//'.nml',IOStat=IO_stat)
      if ( IO_stat /= 0 ) then
        err_msg = 'failed to open file '//trim(filename_loc)//'.nml'
        SLL_ERROR( this_sub_name, err_msg )
      end if

      if (MPI_MASTER) then
        print *,'#initialization with filename:'
        print *,'#',trim(filename_loc)//'.nml'
      endif

      read(input_file, geometry) 
      read(input_file, initial_function)
      read(input_file, time_iterations)
      read(input_file, advector)
      read(input_file, poisson)
      read(input_file, drive)
      close(input_file)

    else

      if (sll_f_get_collective_rank(sll_v_world_collective)==0) then
        print *,'#initialization with default parameters'
      endif      

    endif

    select case (mesh_case_x1) !geometry
      case ("SLL_LANDAU_MESH")
        x1_max = real(nbox_x1,f64) * 2._f64 * sll_p_pi / kmode
        mesh_x1 => sll_f_new_cartesian_mesh_1d(num_cells_x1,eta_min=x1_min, eta_max=x1_max)
        call sll_o_get_node_positions( mesh_x1, sim%x1_array )
      case ("SLL_CARTESIAN_MESH")
        mesh_x1 => sll_f_new_cartesian_mesh_1d(num_cells_x1,eta_min=x1_min, eta_max=x1_max)  
        call sll_o_get_node_positions( mesh_x1, sim%x1_array )
      case default
        err_msg = '#mesh_case_x1 '//mesh_case_x1//' not implemented'
        SLL_ERROR( this_sub_name, err_msg )
    end select

    sim%num_bloc_x1 = 1
    SLL_ALLOCATE(sim%every_x1(sim%num_bloc_x1),ierr)
    SLL_ALLOCATE(sim%bloc_index_x1(sim%num_bloc_x1),ierr)

    sim%every_x1(1) = every_x1
    sim%bloc_index_x1(sim%num_bloc_x1) = num_cells_x1+1
        
    select case (mesh_case_x2)

      case ("SLL_CARTESIAN_MESH")
        mesh_x2 => sll_f_new_cartesian_mesh_1d(num_cells_x2,eta_min=x2_min, eta_max=x2_max)
        call sll_o_get_node_positions( mesh_x2, sim%x2_array )
        SLL_ALLOCATE(sim%x2_array_omp(num_cells_x2+1,0:sim%num_threads-1),ierr)
        do i=0,sim%num_threads-1
          sim%x2_array_omp(:,i) = sim%x2_array(:)
        enddo
        sim%num_bloc_x2 = 1
        SLL_ALLOCATE(sim%every_x2(sim%num_bloc_x2),ierr)
        SLL_ALLOCATE(sim%bloc_index_x2(sim%num_bloc_x2),ierr)
        sim%every_x2(1) = every_x2
        sim%bloc_index_x2(sim%num_bloc_x2) = num_cells_x2+1

      case ("SLL_TWO_GRID_MESH")

        bloc_coord(1) = (x2_fine_min-x2_min)/(x2_max-x2_min)
        bloc_coord(2) = (x2_fine_max-x2_min)/(x2_max-x2_min)
        bloc_index(1) = floor(density_x2_min_to_x2_fine_min)
        bloc_index(2) = floor(density_x2_fine_min_to_x2_fine_max)
        bloc_index(3) = floor(density_x2_fine_max_to_x2_max)
                
        call sll_s_compute_bloc(bloc_coord,bloc_index,num_cells_x2)
        SLL_ALLOCATE(sim%x2_array(num_cells_x2+1),ierr)
        call sll_s_compute_mesh_from_bloc(bloc_coord,bloc_index,sim%x2_array)
        sim%x2_array = x2_min+sim%x2_array*(x2_max-x2_min)
        SLL_ALLOCATE(sim%x2_array_omp(num_cells_x2+1,0:sim%num_threads-1),ierr)
        do i=0,sim%num_threads-1
          sim%x2_array_omp(:,i) = sim%x2_array(:)
        enddo
        mesh_x2 => sll_f_new_cartesian_mesh_1d(num_cells_x2,eta_min=x2_min, eta_max=x2_max)
        sim%num_bloc_x2 = 3
        SLL_ALLOCATE(sim%every_x2(sim%num_bloc_x2),ierr)
        SLL_ALLOCATE(sim%bloc_index_x2(sim%num_bloc_x2),ierr)

        sim%bloc_index_x2(1:3) = bloc_index(1:3)
        sim%every_x2(1) = every_x2_min_to_x2_fine_min
        sim%every_x2(2) = every_x2_fine_min_to_x2_fine_max
        sim%every_x2(3) = every_x2_fine_max_to_x2_max

      case default
        err_msg = '#mesh_case_x2 '//mesh_case_x2//' not implemented'
        SLL_ERROR( this_sub_name, err_msg )

    end select

    sim%mesh2d => mesh_x1 * mesh_x2 ! tensor product
    
    !initial function
    sim%nrj0 = 0._f64
    sim%kx   = kmode
    sim%eps  = eps

    select case (initial_function_case)

      case ("SLL_LANDAU")
        sim%init_func => sll_f_landau_initializer_2d
        SLL_ALLOCATE(sim%params(2),ierr)
        sim%params(1) = kmode
        sim%params(2) = eps
        sim%nrj0      = 0._f64  !compute the right value
        !(0.5_f64*eps*sll_p_pi)**2/(kmode_x1*kmode_x2) &
          !*(1._f64/kmode_x1**2+1._f64/kmode_x2**2)
        !for the moment
        sim%kx = kmode
        sim%eps = eps

      case ("SLL_BUMP_ON_TAIL")
        sim%init_func => sll_f_bump_on_tail_initializer_2d
        SLL_ALLOCATE(sim%params(2),ierr)
        sim%params(1) = kmode
        sim%params(2) = eps
        sim%nrj0 = 0._f64  !compute the right value
        !(0.5_f64*eps*sll_p_pi)**2/(kmode_x1*kmode_x2) &
          !*(1._f64/kmode_x1**2+1._f64/kmode_x2**2)
        !for the moment
        sim%kx = kmode
        sim%eps = eps

      case ("SLL_TWO_STREAM_INSTABILITY")
        sim%init_func => sll_f_two_stream_instability_initializer_2d
        SLL_ALLOCATE(sim%params(2),ierr)
        sim%params(1) = kmode
        sim%params(2) = eps
        sim%nrj0 = 0._f64  !compute the right value
        !(0.5_f64*eps*sll_p_pi)**2/(kmode_x1*kmode_x2) &
          !*(1._f64/kmode_x1**2+1._f64/kmode_x2**2)
        !for the moment
        sim%kx = kmode
        sim%eps = eps

      case ("SLL_BEAM")  
        sim%init_func => sll_f_beam_initializer_2d
        SLL_ALLOCATE(sim%params(1),ierr)
        sim%params(1) = alpha_gaussian             

      case default
        err_msg = '#init_func_case not implemented'
        SLL_ERROR( this_sub_name, err_msg )

    end select

    !ES Equilibrium is initialised to unperturbed normalised Maxwellian by default
    sim%equil_func => sll_f_landau_initializer_2d
    SLL_ALLOCATE(sim%equil_func_params(2),ierr)
    sim%equil_func_params(1) = sim%kx
    sim%equil_func_params(2) = 0._f64

    ! restart file
    sim%time_init_from_restart_file = time_init_from_restart_file
    sim%restart_file = restart_file
    
    !time iterations
    sim%dt                = dt
    sim%num_iterations    = number_iterations
    sim%freq_diag         = freq_diag
    sim%freq_diag_restart = freq_diag_restart
    sim%freq_diag_time    = freq_diag_time
    sim%nb_mode           = nb_mode
    sim%time_init         = time_init
    
    if (sim%nb_mode<0) then
      write( err_msg,* ) '#bad value of nb_mode=', nb_mode, '; #should be >=0'
      SLL_ERROR( this_sub_name, err_msg )
    endif
    
    select case (split_case)    
      case ("SLL_LIE_TV")
        sim%split => sll_f_new_time_splitting_coeff(sll_p_lie_tv)
      case ("SLL_LIE_VT") 
        sim%split => sll_f_new_time_splitting_coeff(sll_p_lie_vt)
      case ("SLL_STRANG_TVT") 
        sim%split => sll_f_new_time_splitting_coeff(sll_p_strang_tvt)
      case ("SLL_STRANG_VTV") 
        sim%split => sll_f_new_time_splitting_coeff(sll_p_strang_vtv)
      case ("SLL_TRIPLE_JUMP_TVT") 
        sim%split => sll_f_new_time_splitting_coeff(sll_p_triple_jump_tvt)
      case ("SLL_TRIPLE_JUMP_VTV") 
        sim%split => sll_f_new_time_splitting_coeff(sll_p_triple_jump_vtv)
      case ("SLL_ORDER6_VTV") 
        sim%split => sll_f_new_time_splitting_coeff(sll_p_order6_vtv)
      case ("SLL_ORDER6VP_TVT") 
        sim%split => sll_f_new_time_splitting_coeff(sll_p_order6vp_tvt,dt=dt)
      case ("SLL_ORDER6VP_VTV") 
        sim%split => sll_f_new_time_splitting_coeff(sll_p_order6vp_vtv,dt=dt)
      case ("SLL_ORDER6VPNEW_TVT") 
        sim%split => sll_f_new_time_splitting_coeff(sll_p_order6vpnew_tvt,dt=dt)
      case ("SLL_ORDER6VPNEW1_VTV") 
        sim%split => sll_f_new_time_splitting_coeff(sll_p_order6vpnew1_vtv,dt=dt)
      case ("SLL_ORDER6VPNEW2_VTV") 
        sim%split => sll_f_new_time_splitting_coeff(sll_p_order6vpnew2_vtv,dt=dt)
      case default
        err_msg = '#split_case not defined'
        SLL_ERROR( this_sub_name, err_msg )
    end select

    !advector
    SLL_ALLOCATE(sim%advect_x1(num_threads),ierr)
    SLL_ALLOCATE(sim%advect_x2(num_threads),ierr)

    !$OMP PARALLEL DEFAULT(SHARED) &
    !$OMP PRIVATE(tid)

#ifdef _OPENMP
    tid = omp_get_thread_num()+1
#else
    tid = 1
#endif

    select case (advector_x1)

      case ("SLL_SPLINES") ! arbitrary order periodic splines

        sim%advect_x1(tid)%ptr => sll_f_new_periodic_1d_advector( &
          num_cells_x1,                                     &
          x1_min,                                           &
          x1_max,                                           &
          sll_p_spline,                                           &  
          order_x1) 

      case("SLL_LAGRANGE") ! arbitrary order sll_p_lagrange periodic interpolation

        sim%advect_x1(tid)%ptr => sll_f_new_periodic_1d_advector( &
          num_cells_x1,                                     &
          x1_min,                                           &
          x1_max,                                           &
          sll_p_lagrange,                                         & 
          order_x1)

      case("SLL_TRIGO") ! sll_p_trigo periodic advection

        sim%advect_x1(tid)%ptr => sll_f_new_periodic_1d_advector( &
          num_cells_x1,                                     &
          x1_min,                                           &
          x1_max,                                           &
          sll_p_trigo,                                       &
          order_x1)

      case("SLL_SPECTRAL") ! spectral periodic advection

        sim%advect_x1(tid)%ptr => sll_f_new_spectral_1d_advector( &
          num_cells_x1,                                     &
          x1_min,                                           &
          x1_max)


      case default

        err_msg = '#advector in x1 '//advector_x1//' not implemented'
        SLL_ERROR( this_sub_name, err_msg )

    end select    

    select case (advector_x2)

      case ("SLL_SPLINES") ! arbitrary order periodic splines

        sim%advect_x2(tid)%ptr => sll_f_new_periodic_1d_advector( &
          num_cells_x2,                                     &
          x2_min,                                           &
          x2_max,                                           &
          sll_p_spline,                                           & 
          order_x2) 

      case("SLL_LAGRANGE") ! arbitrary order sll_p_lagrange periodic interpolation

        sim%advect_x2(tid)%ptr => sll_f_new_periodic_1d_advector( &
          num_cells_x2,                                     &
          x2_min,                                           &
          x2_max,                                           &
          sll_p_lagrange,                                         & 
          order_x2)

      case("SLL_NON_UNIFORM_CUBIC_SPLINES") ! arbitrary order sll_p_lagrange 
                                            ! periodic interpolation

        sim%advect_x2(tid)%ptr =>                    &
          sll_f_new_non_uniform_cubic_splines_1d_advector( &
            num_cells_x2,                            &
            x2_min,                                  &
            x2_max,                                  &
            order_x2,                                &
            sim%x2_array)

      case default

        err_msg = '#advector in x2 '//advector_x2//' not implemented'
        SLL_ERROR( this_sub_name, err_msg )

    end select

    !$OMP END PARALLEL

    select case (advection_form_x2)
      case ("SLL_ADVECTIVE")
        sim%advection_form_x2 = SLL_ADVECTIVE
        sim%num_dof_x2        = num_cells_x2+1
      case ("SLL_CONSERVATIVE")
        sim%advection_form_x2 = SLL_CONSERVATIVE
        sim%num_dof_x2        = num_cells_x2
      case default
        err_msg = '#advection_form_x2 '//advection_form_x2//' not implemented'
        SLL_ERROR( this_sub_name, err_msg )
    end select  
    
    sim%factor_x1     = factor_x1
    sim%factor_x2_rho = factor_x2_rho
    sim%factor_x2_1   = factor_x2_1
    
    SLL_ALLOCATE(sim%integration_weight(sim%num_dof_x2),ierr)

    select case (integration_case)
      case ("SLL_RECTANGLE")
        do i=1,num_cells_x2
          sim%integration_weight(i) = sim%x2_array(i+1)-sim%x2_array(i)
        enddo
        sim%integration_weight(num_cells_x2+1) = 0._f64   
      case ("SLL_TRAPEZOID")
        sim%integration_weight(1)=0.5_f64*(sim%x2_array(2)-sim%x2_array(1))
        do i=2,num_cells_x2
          sim%integration_weight(i) = 0.5_f64*(sim%x2_array(i+1)-sim%x2_array(i-1))
        enddo  
        sim%integration_weight(num_cells_x2+1) = &
          0.5_f64*(sim%x2_array(num_cells_x2+1)-sim%x2_array(num_cells_x2))
      case ("SLL_CONSERVATIVE")
        do i=1,num_cells_x2
          sim%integration_weight(i)=sim%x2_array(i+1)-sim%x2_array(i)
        enddo  
      case default
        err_msg = '#integration_case not implemented'
        SLL_ERROR( this_sub_name, err_msg )
    end select  
    
    select case (poisson_solver)
      case ("SLL_FFT")
        sim%poisson => sll_f_new_poisson_1d_periodic( &
          x1_min, &
          x1_max, &
          num_cells_x1)
      case ("SLL_POLAR")
        sim%poisson => sll_f_new_poisson_1d_polar( &
          x1_min, &
          x1_max, &
          num_cells_x1)
      case default
        err_msg = '#poisson_solver '//poisson_solver//' not implemented'
        SLL_ERROR( this_sub_name, err_msg )
    end select

    select case (ampere_solver)
      case ("SLL_VLASOV_AMPERE")
        print*,'########################'
        print*,'# Vlasov-Ampere scheme #'
        print*,'########################'
        sim%ampere = .true.
        SLL_ALLOCATE(sim%advect_ampere_x1(num_threads),ierr)
        tid = 1
        !$OMP PARALLEL DEFAULT(SHARED) &
        !$OMP PRIVATE(tid)
        !$ tid = omp_get_thread_num()+1
        sim%advect_ampere_x1(tid)%ptr => sll_f_new_ampere_1d_advector( &
          num_cells_x1, &
          x1_min,       &
          x1_max )
        !$OMP END PARALLEL
      case default
        continue
    end select
    
    select case (drive_type)

      case ("SLL_NO_DRIVE")
        sim%driven         = .false.
        sim%t0             = 0.0_f64
        sim%twL            = 0.0_f64
        sim%twR            = 0.0_f64
        sim%tflat          = 0.0_f64
        sim%tL             = 0.0_f64
        sim%tR             = 0.0_f64
        sim%turn_drive_off = .false.
        sim%Edrmax         = 0.0_f64
        sim%omegadr        = 0.0_f64

      case("SLL_KEEN_DRIVE")
        sim%driven         = .true.
        sim%t0             = keen_t0
        sim%twL            = keen_twL
        sim%twR            = keen_twR
        sim%tflat          = keen_tflat
        sim%tL             = keen_tL
        sim%tR             = keen_tR
        sim%turn_drive_off = keen_turn_drive_off
        sim%Edrmax         = keen_Edrmax
        sim%omegadr        = keen_omegadr        

      case default
        err_msg = '#drive_type '//drive_type//' not implemented'
        SLL_ERROR( this_sub_name, err_msg )

    end select

    if (sll_f_get_collective_rank(sll_v_world_collective)==0) then
            
        open(unit = param_out, file = 'parameters.dat')
      
        write(param_out, *) real(number_iterations,f64)*dt !Tmax
        write(param_out, *) dt                             !dt
        write(param_out, *) number_iterations              !Nt
        write(param_out, *) kmode                          !k0
        write(param_out, *) sim%omegadr                    !omega0
        write(param_out, *) x1_max-x1_min                  !L
        write(param_out, *) nbox_x1                        !mbox
        write(param_out, *) sim%Edrmax                     !Edrmax
        write(param_out, *) freq_diag_time                 !nsave
        write(param_out, *) freq_diag                      !nsavef1
        write(param_out, *) freq_diag                      !nsavef2
        write(param_out, *) sim%tR+sim%tflat               !Tsetup
        write(param_out, *) x2_max                         !vxmax
        write(param_out, *) x2_min                         !vxmin
        write(param_out, *) num_cells_x1                   !N
        write(param_out, *) num_cells_x2                   !Nv
        ! KEEN drive
        write(param_out, *) sim%tL                         !tL
        write(param_out, *) sim%tR                         !tR
        write(param_out, *) sim%twL                        !twL
        write(param_out, *) sim%twR                        !twR
        write(param_out, *) sim%tflat                      !tflat
            
        close(param_out)
      
    endif

    sim%L = x1_max - x1_min

  end subroutine init_va2d_par_cart

  subroutine init_va2d_fake(sim, filename)

    class(sll_t_simulation_2d_vlasov_ampere_cart), intent(inout) :: sim
    character(len=*), intent(in)                               :: filename
  
    character(len=*), parameter :: this_sub_name = 'init_va2d_fake'
    character(len=288)          :: err_msg

    print *,sim%dt
    print *,filename
    err_msg = '# Do not use the routine init_va2d_fake\n'// & 
              '# Use instead init_va2d_par_cart'
    SLL_ERROR( this_sub_name, err_msg )
  
  end subroutine init_va2d_fake

  subroutine run_va2d_cartesian(sim)

    class(sll_t_simulation_2d_vlasov_ampere_cart), intent(inout) :: sim
    character(len=*), parameter :: this_sub_name = 'run_va2d_cartesian'
    character(len=*), parameter :: this_prog_name = 'vlasov_ampere_2d'
    character(len=128)          :: err_msg
    
    sll_int32           :: ierr
    sll_int32           :: i
    
    procedure(sll_i_scalar_initializer_2d), pointer :: init_func
    
    sll_real64, dimension(:,:), pointer :: f_x1
    sll_real64, dimension(:,:), pointer :: f_x2
    sll_real64, dimension(:,:), pointer :: f_x1_equil
    
    sll_real64, dimension(:),   pointer :: rho
    sll_real64, dimension(:),   pointer :: efield
    sll_real64, dimension(:),   pointer :: e_app
    sll_real64, dimension(:),   pointer :: current
    
    
    type(sll_t_layout_2d),            pointer :: layout_x1
    type(sll_t_layout_2d),            pointer :: layout_x2
    type(sll_t_remap_plan_2d_real64), pointer :: remap_plan_x1_x2
    type(sll_t_remap_plan_2d_real64), pointer :: remap_plan_x2_x1
    sll_real64, dimension(:),   pointer :: f1d
    sll_int32                           :: np_x1
    sll_int32                           :: np_x2
    sll_int32                           :: nproc_x1
    sll_int32                           :: nproc_x2
    sll_int32                           :: global_indices(2)
    sll_int32                           :: local_size_x1
    sll_int32                           :: local_size_x2
    
    sll_real64, dimension(:), allocatable :: x2_array_middle
    
    
    sll_int32  :: nb_mode 
    sll_int32  :: split_istep
    sll_int32  :: num_dof_x2 
    sll_real64 :: time_init
    sll_real64 :: t_step
     
    logical :: split_T
    
    ! for parallelization (output of distribution function in one single file)
    sll_int32                              :: collective_size
    character(len=4)                       :: cproc
    character(len=4)                       :: cplot
    sll_int32                              :: iproc
    
    
    iplot = 1
    
    nb_mode          = sim%nb_mode
    time_init        = sim%time_init
    np_x1            = sim%mesh2d%num_cells1+1
    np_x2            = sim%mesh2d%num_cells2+1
    num_dof_x2       = sim%num_dof_x2

    collective_size = sll_f_get_collective_size(sll_v_world_collective)

    if (MPI_MASTER) then
    
      print *,'#collective_size=', collective_size
      SLL_ALLOCATE(f_visu(np_x1,num_dof_x2),ierr)
      SLL_ALLOCATE(f_visu_buf1d(np_x1*num_dof_x2),ierr)
    
    else
    
      SLL_ALLOCATE(f_visu(1:1,1:1),ierr)          
      SLL_ALLOCATE(f_visu_buf1d(1:1),ierr)
    
    endif
    
    SLL_ALLOCATE(collective_displs(collective_size),ierr)
    SLL_ALLOCATE(collective_recvcnts(collective_size),ierr)
    
    SLL_ALLOCATE(buf_fft(np_x1-1),ierr)
    call sll_s_fft_init_r2r_1d(pfwd,np_x1-1,buf_fft,buf_fft,sll_p_fft_forward,normalized = .TRUE.)
    
    layout_x1       => sll_f_new_layout_2d( sll_v_world_collective )
    layout_x2       => sll_f_new_layout_2d( sll_v_world_collective )    
    nproc_x1 = sll_f_get_collective_size( sll_v_world_collective )
    nproc_x2 = 1
    call sll_o_initialize_layout_with_distributed_array( &
      np_x1, num_dof_x2, nproc_x1, nproc_x2, layout_x2 )
    call sll_o_initialize_layout_with_distributed_array( &
      np_x1, num_dof_x2, nproc_x2, nproc_x1, layout_x1 )
    
    call sll_o_view_lims( layout_x1 )
    call sll_o_view_lims( layout_x2 )
    
    call sll_o_compute_local_sizes( layout_x2, local_size_x1, local_size_x2 )
    SLL_ALLOCATE(f_x2(local_size_x1,local_size_x2),ierr)
    
    call sll_o_compute_local_sizes( layout_x1, local_size_x1, local_size_x2 )
    global_indices(1:2) = sll_o_local_to_global( layout_x1, (/1, 1/) )
    SLL_ALLOCATE(f_x1(local_size_x1,local_size_x2),ierr)    
    SLL_ALLOCATE(f_x1_equil(local_size_x1,local_size_x2),ierr)    
    SLL_ALLOCATE(f_x1_buf1d(local_size_x1*local_size_x2),ierr)    
    
    remap_plan_x1_x2 => sll_o_new_remap_plan(layout_x1, layout_x2, f_x1)
    remap_plan_x2_x1 => sll_o_new_remap_plan(layout_x2, layout_x1, f_x2)
    
    SLL_ALLOCATE(rho(np_x1),ierr)
    SLL_ALLOCATE(current(np_x1),ierr)
    SLL_ALLOCATE(efield(np_x1),ierr)
    SLL_ALLOCATE(e_app(np_x1),ierr)
    SLL_ALLOCATE(f1d(max(np_x1,np_x2)),ierr)
    SLL_ALLOCATE(f1d_omp_in(max(np_x1,np_x2),sim%num_threads),ierr)
    SLL_ALLOCATE(f1d_omp_out(max(np_x1,np_x2),sim%num_threads),ierr)
    SLL_ALLOCATE(x2_array_unit(np_x2),ierr)
    SLL_ALLOCATE(x2_array_middle(np_x2),ierr)
    SLL_ALLOCATE(sim%node_positions_x2(num_dof_x2),ierr)
    
    x2_array_unit(1:np_x2) = &
      (sim%x2_array(1:np_x2)-sim%x2_array(1))/(sim%x2_array(np_x2)-sim%x2_array(1))
    do i = 1, np_x2-1
       x2_array_middle(i) = 0.5_f64*(sim%x2_array(i)+sim%x2_array(i+1))
    end do
    x2_array_middle(np_x2) = x2_array_middle(1)+sim%x2_array(np_x2)-sim%x2_array(1)
    
    select case (sim%advection_form_x2)
      case (SLL_ADVECTIVE)
        sim%node_positions_x2(1:num_dof_x2) = sim%x2_array(1:num_dof_x2)
      case (SLL_CONSERVATIVE)
        sim%node_positions_x2(1:num_dof_x2) = x2_array_middle(1:num_dof_x2)
      case default
        write(err_msg,*) '#sim%advection_form_x2=', sim%advection_form_x2, &
                         ' not implemented'
        SLL_ERROR( this_prog_name, err_msg )
    end select  
        
    call sll_o_2d_parallel_array_initializer_cartesian( &
       layout_x1,                                     &
       sim%x1_array,                                  &
       sim%node_positions_x2,                         &
       f_x1,                                          &
       sim%init_func,                                 &
       sim%params)
    
    iproc = sll_f_get_collective_rank(sll_v_world_collective)
    call sll_s_int2string(iproc, cproc)
    call sll_s_int2string(iplot, cplot)    
    
    call check_restart(sim, f_x1)
    
    !ES initialise f_x1_equil that is used to define deltaf: 
    !ES deltaf =  f_x1 - f_x1_equil
    call sll_o_2d_parallel_array_initializer_cartesian( &
       layout_x1,                                     &
       sim%x1_array,                                  &
       sim%node_positions_x2,                         &
       f_x1_equil,                                    &
       sim%equil_func,                                &
       sim%equil_func_params)
    
    call sll_s_compute_displacements_array_2d( layout_x1,       &
                                         collective_size, &
                                         collective_displs )
    
    collective_recvcnts = sll_f_receive_counts_array_2d( layout_x1, &
                                                   collective_size )
    
    !ES replace f_x1 by f_x1_equil and write f_x1_equil into 
    !ES f0.bdat instead of f_x1  
    
    call sll_s_load_buffer_2d( layout_x1, f_x1_equil, f_x1_buf1d )
    
    call sll_s_collective_gatherv_real64( sll_v_world_collective,        &
                                        f_x1_buf1d,                  &
                                        local_size_x1*local_size_x2, &
                                        collective_recvcnts,         &
                                        collective_displs,           &
                                        0,                           &
                                        f_visu_buf1d )
        
    f_visu = reshape(f_visu_buf1d, shape(f_visu))
    
    call compute_rho(sim, layout_x1, f_x1, rho)
    
    call sim%poisson%compute_E_from_rho( efield, rho )
        
    ! Ponderomotive force at initial time. We use a sine wave
    ! with parameters k_dr and omega_dr.
    istep = 0
    
    if (sim%driven) then
      call set_e_app(sim, sim%time_init, e_app)
    else
      e_app = 0._f64
    end if
    
    if (MPI_MASTER) call write_init_files(sim, efield, rho, e_app)
    
    call sll_s_load_buffer_2d( layout_x1, f_x1-f_x1_equil, f_x1_buf1d )
    
    call sll_s_collective_gatherv_real64( sll_v_world_collective,        &
                                        f_x1_buf1d,                  &
                                        local_size_x1*local_size_x2, &
                                        collective_recvcnts,         &
                                        collective_displs,           &
                                        0,                           &
                                        f_visu_buf1d )
    
    f_visu = reshape(f_visu_buf1d, shape(f_visu))
    
    if (MPI_MASTER) then
      call sll_s_binary_write_array_2d(deltaf_id,f_visu(1:np_x1-1,1:np_x2-1),ierr)
      print *,'#step=',0,sim%time_init,'iplot=',iplot
    endif
    
    iplot = iplot+1  
    
    do istep = 1, sim%num_iterations
    
       split_T = sim%split%split_begin_T
       t_step = real(istep-1,f64)
    
       do split_istep=1,sim%split%nb_split_step
    
         if (split_T) then
    
           if (sim%ampere) then
             call advection_ampere_x(sim,       &
                                     layout_x1, &
                                     efield,    &
                                     f_x1,      &
                                     sim%split%split_step(split_istep)*sim%dt)
           else
             call advection_poisson_x( sim,            &
               layout_x1,                              &
               f_x1,                                   &
               efield,                                 &
               rho,                                    &
               sim%split%split_step(split_istep)*sim%dt)
           end if
           t_step = t_step+sim%split%split_step(split_istep)
    
         else
    
           if (sim%driven) call set_e_app(sim, sim%time_init+(istep-1)*sim%dt, e_app)
           call sll_o_apply_remap_2d( remap_plan_x1_x2, f_x1, f_x2 )
           call advection_v(sim, layout_x2, f_x2, efield, e_app, &
                               sim%split%split_step(split_istep)*sim%dt)
           call sll_o_apply_remap_2d( remap_plan_x2_x1, f_x2, f_x1 )
    
         endif
         split_T = .not.(split_T)
       enddo
    
      if (mod(istep,sim%freq_diag_time)==0) then
    
        call diagnostics(sim, layout_x1, f_x1, rho, efield, e_app)
        
        if (mod(istep,sim%freq_diag_restart)==0) then          
          call save_for_restart(layout_x1, f_x1)
        endif 
    
        if (mod(istep,sim%freq_diag)==0) then          
    
          if (MPI_MASTER) then        
            print *,'#step=',istep,sim%time_init+real(istep,f64)*sim%dt,'iplot=',iplot
          endif
    
          !call gnuplot_write(sim, layout_x1, f_x1-f_x1_equil, 'deltaf', 'intdeltafdx')
          !call gnuplot_write(sim, layout_x1, f_x1,                 'f',      'intfdx')
    
          iplot = iplot+1  
                    
        endif
    
     end if
    
    enddo
    
    if (MPI_MASTER) then
      call sll_s_ascii_file_close(thdiag_id,ierr) 
      call sll_s_binary_file_close(deltaf_id,ierr) 
      call sll_s_binary_file_close(efield_id,ierr)
      call sll_s_binary_file_close(rhotot_id,ierr)
      call sll_s_binary_file_close(t_id,ierr)
      if (sim%driven) then
        call sll_s_binary_file_close(Edr_id,ierr)
        call sll_s_binary_file_close(adr_id,ierr)
      endif   
    endif
    
  end subroutine run_va2d_cartesian
    
  subroutine advection_poisson_x(sim, &
     layout_x1,                       &
     f_x1,                            &
     efield,                          &
     rho,                             &
     delta_t)

    class(sll_t_simulation_2d_vlasov_ampere_cart), intent(inout) :: sim
    type(sll_t_layout_2d), intent(in),    pointer        :: layout_x1
    sll_real64,      intent(inout), dimension(:,:) :: f_x1
    sll_real64,      intent(out),   dimension(:)   :: efield
    sll_real64,      intent(out),   dimension(:)   :: rho
    sll_real64 :: delta_t
    sll_int32  :: ig_omp
    sll_int32  :: i_omp
    sll_int32  :: global_indices(2)
    sll_int32  :: local_size_x1
    sll_int32  :: local_size_x2
    sll_int32  :: np_x1
    sll_int32  :: tid
    sll_real64 :: alpha_omp
    
    np_x1 = sim%mesh2d%num_cells1+1
    call sll_o_compute_local_sizes( layout_x1, local_size_x1, local_size_x2 )
    global_indices = sll_o_local_to_global( layout_x1, (/1, 1/) )
    
    tid=1          
    !$OMP PARALLEL DEFAULT(SHARED) &
    !$OMP PRIVATE(i_omp,ig_omp,alpha_omp,tid) 
    !advection in x
    !$ tid = omp_get_thread_num()+1
    !$OMP DO
    do i_omp = 1, local_size_x2
    
      ig_omp    = i_omp+global_indices(2)-1
      alpha_omp = sim%factor_x1*sim%node_positions_x2(ig_omp)
      f1d_omp_in(1:np_x1,tid) = f_x1(1:np_x1,i_omp)
      
      call sim%advect_x1(tid)%ptr%advect_1d_constant(  &
        alpha_omp,                                     &
        delta_t,                                       &
        f1d_omp_in(1:np_x1,tid),                       &
        f1d_omp_out(1:np_x1,tid))
    
      f_x1(1:np_x1,i_omp)=f1d_omp_out(1:np_x1,tid)
    
    end do
    !$OMP END DO          
    !$OMP END PARALLEL
    
    call compute_rho(sim, layout_x1, f_x1, rho)
    call sim%poisson%compute_E_from_rho( efield, rho )
    
  end subroutine advection_poisson_x
    
  subroutine advection_ampere_x(sim, layout_x1, efield, f_x1, delta_t)
    class(sll_t_simulation_2d_vlasov_ampere_cart), intent(inout) :: sim
    type(sll_t_layout_2d) , pointer :: layout_x1
    sll_real64 :: efield(:)
    sll_real64 :: f_x1(:,:)
    sll_int32  :: local_size_x1
    sll_int32  :: local_size_x2
    sll_int32  :: global_indices(2)
    
    sll_real64 :: delta_t
    sll_int32  :: nc_x1
    sll_int32  :: np_x1
    sll_comp64 :: s
    sll_int32  :: tid, ig_omp, i, i_omp
    sll_real64 :: alpha_omp
    
    call sll_o_compute_local_sizes( layout_x1, local_size_x1, local_size_x2 )
    global_indices = sll_o_local_to_global( layout_x1, (/1, 1/) )
    
    np_x1 = sim%mesh2d%num_cells1+1
    nc_x1 = np_x1-1
    
    tid=1          
    
    
    !$OMP PARALLEL DEFAULT(SHARED) &
    !$OMP PRIVATE(i_omp,ig_omp,alpha_omp,tid) 
    !advection in x
    !$ tid = omp_get_thread_num()+1
    
    sim%advect_ampere_x1(tid)%ptr%r1 = cmplx(0.0,0.0,kind=f64)
    !$OMP DO 
    do i_omp = 1, local_size_x2
    
      ig_omp    = i_omp+global_indices(2)-1
      alpha_omp = sim%factor_x1*sim%node_positions_x2(ig_omp)*delta_t
      f1d_omp_in(1:np_x1,tid) = f_x1(1:np_x1,i_omp)
      
      sim%advect_ampere_x1(tid)%ptr%d_dx = f1d_omp_in(1:nc_x1,tid)
    
      call sll_s_fft_exec_r2c_1d(sim%advect_ampere_x1(tid)%ptr%fwx,  &
           sim%advect_ampere_x1(tid)%ptr%d_dx, &
           sim%advect_ampere_x1(tid)%ptr%fk)
      do i = 2, nc_x1/2+1
        sim%advect_ampere_x1(tid)%ptr%fk(i) = &
           sim%advect_ampere_x1(tid)%ptr%fk(i) & 
           * cmplx(cos(sim%advect_ampere_x1(tid)%ptr%kx(i)*alpha_omp), &
                  -sin(sim%advect_ampere_x1(tid)%ptr%kx(i)*alpha_omp),kind=f64)
      end do
    
      sim%advect_ampere_x1(tid)%ptr%r1(2:nc_x1/2+1) = &
           sim%advect_ampere_x1(tid)%ptr%r1(2:nc_x1/2+1) &
         + sim%advect_ampere_x1(tid)%ptr%fk(2:nc_x1/2+1) * sim%integration_weight(ig_omp)
    
      call sll_s_fft_exec_c2r_1d(sim%advect_ampere_x1(tid)%ptr%bwx, &
           sim%advect_ampere_x1(tid)%ptr%fk,  &
           sim%advect_ampere_x1(tid)%ptr%d_dx)
    
      f1d_omp_out(1:nc_x1, tid) = sim%advect_ampere_x1(tid)%ptr%d_dx/nc_x1
      f1d_omp_out(np_x1, tid)   = f1d_omp_out(1, tid) 
    
      f_x1(1:np_x1,i_omp)=f1d_omp_out(1:np_x1,tid)
    
    end do
    !$OMP END DO          
    
    !$OMP END PARALLEL
    
    
    sim%advect_ampere_x1(tid)%ptr%d_dx = efield(1:nc_x1)
    call sll_s_fft_exec_r2c_1d(sim%advect_ampere_x1(1)%ptr%fwx,  &
         sim%advect_ampere_x1(1)%ptr%d_dx, &
         sim%advect_ampere_x1(1)%ptr%ek)
    
    do i = 2, nc_x1/2+1
      s = cmplx(0.0,0.0,kind=f64)
      do tid = 1, sim%num_threads
        s = s + sim%advect_ampere_x1(tid)%ptr%r1(i)
      end do
      sim%advect_ampere_x1(1)%ptr%r1(i) = s
    end do
    
    do i = 2, nc_x1/2+1
      sim%advect_ampere_x1(1)%ptr%ek(i) =  &
         - sim%advect_ampere_x1(1)%ptr%r1(i) * sim%L / cmplx(0.,2.*sll_p_pi*(i-1),kind=f64)
    end do
    
    call sll_s_fft_exec_c2r_1d(sim%advect_ampere_x1(1)%ptr%bwx, &
         sim%advect_ampere_x1(1)%ptr%ek,  &
         efield)
    
    efield(1:nc_x1) = efield(1:nc_x1) / nc_x1
    efield(np_x1) = efield(1)
    
    
    !call compute_rho()
    !call sim%poisson%compute_E_from_rho( efield, rho )
    
  end subroutine advection_ampere_x
    
  subroutine compute_rho(sim, layout_x1, f_x1, rho)
    class(sll_t_simulation_2d_vlasov_ampere_cart), intent(inout) :: sim
    type(sll_t_layout_2d), pointer :: layout_x1
    sll_real64      :: rho(:)
    sll_real64      :: f_x1(:,:)
    sll_int32       :: local_size_x1, local_size_x2
    sll_int32       :: global_indices(2)
    sll_int32       :: np_x1
    sll_int32       :: i
    sll_int32       :: ig
    sll_int32       :: ierr

    sll_real64, dimension(:), allocatable :: rho_loc

    np_x1 = sim%mesh2d%num_cells1+1
    SLL_ALLOCATE(rho_loc(np_x1),ierr)
    
    call sll_o_compute_local_sizes( layout_x1, local_size_x1, local_size_x2 )
    global_indices = sll_o_local_to_global( layout_x1, (/1, 1/) )
    !computation of electric field
    rho_loc = 0._f64
    ig = global_indices(2)-1
    do i=1,np_x1
      rho_loc(i)=rho_loc(i)                              &
        +sum(f_x1(i,1:local_size_x2)                     &
        *sim%integration_weight(1+ig:local_size_x2+ig))
    end do
        
    call sll_o_collective_allreduce( sll_v_world_collective, &
                                   rho_loc,              &
                                   np_x1,                &
                                   MPI_SUM,              &
                                   rho )
    
    rho = sim%factor_x2_1*1._f64-sim%factor_x2_rho*rho
    
  end subroutine compute_rho
    
! subroutine compute_current(sim, layout_x1, f_x1, current)
!   
!   class(sll_t_simulation_2d_vlasov_ampere_cart), intent(inout) :: sim
!   type(sll_t_layout_2d), pointer :: layout_x1
!   sll_real64      :: current(:)
!   sll_real64      :: f_x1(:,:)
!   sll_int32       :: np_x1
!   sll_int32       :: i
!   sll_int32       :: ig
!   sll_int32       :: ierr
!   sll_int32  :: global_indices(2)
!   sll_int32  :: local_size_x1
!   sll_int32  :: local_size_x2
!   sll_real64 :: v
!   sll_int32  :: j
!   sll_int32  :: gj
!   sll_real64, dimension(:), allocatable :: rho_loc
!   
!   np_x1 = sim%mesh2d%num_cells1+1
!   SLL_ALLOCATE(rho_loc(np_x1),ierr)
!   call sll_o_compute_local_sizes( layout_x1, local_size_x1, local_size_x2 )
!   global_indices = sll_o_local_to_global( layout_x1, (/1, 1/) )
!   
!   do i = 1,local_size_x1
!     rho_loc(i) = 0._f64
!     do j = 1,local_size_x2
!       global_indices = sll_o_local_to_global( layout_x1, (/i, j/) )
!       gj = global_indices(2)
!       v  = sim%node_positions_x2(gj)
!       rho_loc(i) = rho_loc(i)+f_x1(i,j)*sim%integration_weight(gj)*v
!     end do
!   end do
!   
!   call sll_o_collective_allreduce( sll_v_world_collective, &
!                                  rho_loc,              &
!                                  np_x1,                &
!                                  MPI_SUM,              &
!                                  current )
!   
! end subroutine compute_current
!    
  subroutine advection_v(sim, layout_x2, f_x2, efield, e_app, delta_t)

    class(sll_t_simulation_2d_vlasov_ampere_cart), intent(inout) :: sim
    type(sll_t_layout_2d), pointer :: layout_x2
    sll_real64 :: f_x2(:,:)
    sll_real64 :: efield(:)
    sll_real64 :: e_app(:)
    sll_real64 :: delta_t, alpha_omp
    sll_int32  :: i_omp, ig_omp
    sll_int32  :: local_size_x1, local_size_x2
    sll_int32  :: global_indices(2)
    sll_int32  :: tid, np_x2
    sll_real64 :: mean_omp
    
    np_x2 = sim%mesh2d%num_cells2+1
    call sll_o_compute_local_sizes( layout_x2, local_size_x1, local_size_x2 )
    global_indices = sll_o_local_to_global( layout_x2, (/1, 1/) )
    tid = 1
    !$OMP PARALLEL DEFAULT(SHARED) &
    !$OMP PRIVATE(i_omp,ig_omp,alpha_omp,tid,mean_omp) 
    !advection in v
    !$ tid = omp_get_thread_num()+1
    !advection in v
    !$OMP DO
    do i_omp = 1,local_size_x1
    
      ig_omp=i_omp+global_indices(1)-1
    
      alpha_omp = -(efield(ig_omp)+e_app(ig_omp))
    
      f1d_omp_in(1:sim%num_dof_x2,tid) = f_x2(i_omp,1:sim%num_dof_x2)
    
      if (sim%advection_form_x2==SLL_CONSERVATIVE) then
    
        call sll_s_function_to_primitive(f1d_omp_in(:,tid),    &
                                   x2_array_unit,        &
                                   np_x2-1,mean_omp)
      endif
    
      call sim%advect_x2(tid)%ptr%advect_1d_constant(    &
        alpha_omp,                                       &
        delta_t,                                         &
        f1d_omp_in(1:sim%num_dof_x2,tid),                &
        f1d_omp_out(1:sim%num_dof_x2,tid))
    
      if (sim%advection_form_x2==SLL_CONSERVATIVE) then
    
        call sll_s_primitive_to_function(f1d_omp_out(:,tid),   &
                                   x2_array_unit,        &
                                   np_x2-1,              &
                                   mean_omp)
      endif
      f_x2(i_omp,1:sim%num_dof_x2) = f1d_omp_out(1:sim%num_dof_x2,tid)
    end do
    !$OMP END DO          
    !$OMP END PARALLEL
    
  end subroutine advection_v
    
    
  subroutine gnuplot_write( sim, layout_x1, f, fname, intfname)
    
    class(sll_t_simulation_2d_vlasov_ampere_cart), intent(inout) :: sim
    type(sll_t_layout_2d), pointer   :: layout_x1
    sll_real64, dimension(:,:) :: f
    character(len=*)           :: fname
    character(len=*)           :: intfname
    
    sll_int32 :: np_x1, np_x2
    sll_int32 :: local_size_x1, local_size_x2
    sll_int32 :: i, ierr
    
    np_x1 = sim%mesh2d%num_cells1+1
    np_x2 = sim%mesh2d%num_cells2+1
    call sll_s_load_buffer_2d( layout_x1, f, f_x1_buf1d )
    call sll_o_compute_local_sizes( layout_x1, local_size_x1, local_size_x2 )
    
    call sll_s_collective_gatherv_real64( &
      sll_v_world_collective,             &
      f_x1_buf1d,                       &
      local_size_x1*local_size_x2,      &
      collective_recvcnts,              &
      collective_displs,                &
      0,                                &
      f_visu_buf1d )
    
    f_visu = reshape(f_visu_buf1d, shape(f_visu))
    
    if (MPI_MASTER) then
    
      do i=1,sim%num_dof_x2
        f_visu_buf1d(i) = sum(f_visu(1:np_x1-1,i))*sim%mesh2d%delta_eta1
      enddo
    
      call sll_o_gnuplot_1d(                       &
        f_visu_buf1d(1:sim%num_dof_x2),          &
        sim%node_positions_x2(1:sim%num_dof_x2), &
        intfname,                                &
        iplot )
    
      call sll_s_binary_write_array_2d(deltaf_id,                      &
                                     f_visu(1:np_x1-1,1:np_x2-1),    &
                                     ierr)  
    
      call sll_s_plot_f_cartesian(     &
            iplot,                   &
            f_visu,                  &
            sim%x1_array,            &
            np_x1,                   &
            sim%node_positions_x2,   &
            sim%num_dof_x2,          &
            fname,time)                    
    
    
    endif
    
  end subroutine gnuplot_write
    
  subroutine diagnostics(sim, layout_x1, f_x1, rho, efield, e_app)
    
    class(sll_t_simulation_2d_vlasov_ampere_cart), intent(in) :: sim

    type(sll_t_layout_2d), pointer :: layout_x1
    sll_real64, intent(in)   :: f_x1(:,:)
    sll_real64, intent(in)   :: rho(:)
    sll_real64, intent(in)   :: efield(:)
    sll_real64, intent(in)   :: e_app(:)

    sll_int32  :: local_size_x1, local_size_x2
    sll_int32  :: global_indices(2)
    sll_real64 :: tmp_loc(5), tmp(5)
    sll_real64 :: mass, l1norm, l2norm, momentum
    sll_real64 :: potential_energy, kinetic_energy
    sll_int32  :: i, ig, k, np_x1, ierr

    sll_real64, dimension(:), allocatable :: f_hat_x2_loc
    sll_real64, dimension(:), allocatable :: f_hat_x2
    sll_comp64, dimension(:), allocatable :: rho_mode

    SLL_ALLOCATE(rho_mode(0:sim%nb_mode),ierr) ; rho_mode = cmplx(0.,0.,f64)
    SLL_CLEAR_ALLOCATE(f_hat_x2(1:sim%nb_mode+1),ierr)
    SLL_CLEAR_ALLOCATE(f_hat_x2_loc(1:sim%nb_mode+1),ierr)

    np_x1            = sim%mesh2d%num_cells1+1
    call sll_o_compute_local_sizes( layout_x1, local_size_x1, local_size_x2 )
    global_indices = sll_o_local_to_global( layout_x1, (/1, 1/) )
    
    time             = sim%time_init+real(istep,f64)*sim%dt
    mass             = 0._f64
    momentum         = 0._f64
    l1norm           = 0._f64
    l2norm           = 0._f64
    kinetic_energy   = 0._f64
    potential_energy = 0._f64
    tmp_loc          = 0._f64
    ig               = global_indices(2)-1               
    
    
    do i = 1, np_x1-1        
      tmp_loc(1)= tmp_loc(1)+sum(f_x1(i,1:local_size_x2) &
        *sim%integration_weight(1+ig:local_size_x2+ig))
      tmp_loc(2)= tmp_loc(2)+sum(abs(f_x1(i,1:local_size_x2)) &
        *sim%integration_weight(1+ig:local_size_x2+ig))
      tmp_loc(3)= tmp_loc(3)+sum((f_x1(i,1:local_size_x2))**2 &
        *sim%integration_weight(1+ig:local_size_x2+ig))
      tmp_loc(4)= tmp_loc(4) +sum(f_x1(i,1:local_size_x2) &
        *sim%x2_array(global_indices(2)-1+1:global_indices(2)-1+local_size_x2) &
        *sim%integration_weight(1+ig:local_size_x2+ig))          
      tmp_loc(5)= tmp_loc(5)+sum(f_x1(i,1:local_size_x2) &
        *sim%x2_array(global_indices(2)-1+1:global_indices(2)-1+local_size_x2)**2 &
        *sim%integration_weight(1+ig:local_size_x2+ig) )          
    end do
    
    call sll_o_collective_allreduce( sll_v_world_collective, &
                                   tmp_loc,              &
                                   5,                    &
                                   MPI_SUM,              &
                                   tmp )
    
    mass           = tmp(1)  * sim%mesh2d%delta_eta1 !* sim%mesh2d%delta_eta2
    l1norm         = tmp(2)  * sim%mesh2d%delta_eta1 !* sim%mesh2d%delta_eta2
    l2norm         = tmp(3)  * sim%mesh2d%delta_eta1 !* sim%mesh2d%delta_eta2
    momentum       = tmp(4)  * sim%mesh2d%delta_eta1 !* sim%mesh2d%delta_eta2
    kinetic_energy = 0.5_f64 * tmp(5) * sim%mesh2d%delta_eta1 !*sim%mesh2d%delta_eta2
    potential_energy = 0._f64
    do i=1, np_x1-1
      potential_energy = potential_energy+(efield(i)+e_app(i))**2
    enddo
    potential_energy = 0.5_f64*potential_energy* sim%mesh2d%delta_eta1
    
    f_hat_x2_loc(1:sim%nb_mode+1) = 0._f64
    do i=1,local_size_x2
      buf_fft = f_x1(1:np_x1-1,i)
      call sll_s_fft_exec_r2r_1d(pfwd,buf_fft,buf_fft)
      do k=0,sim%nb_mode
        f_hat_x2_loc(k+1) = f_hat_x2_loc(k+1) &
          +abs(sll_f_fft_get_mode_r2c_1d(pfwd,buf_fft,k))**2 &
          *sim%integration_weight(ig+i)
      enddo
    enddo
    
    call sll_o_collective_allreduce( sll_v_world_collective, &
                                   f_hat_x2_loc,         &
                                   sim%nb_mode+1,        &
                                   MPI_SUM,              &
                                   f_hat_x2 )
    
    if (MPI_MASTER) then                  
    
      buf_fft = rho(1:np_x1-1)
      call sll_s_fft_exec_r2r_1d(pfwd,buf_fft,buf_fft)
    
      do k=0,sim%nb_mode
        rho_mode(k)=sll_f_fft_get_mode_r2c_1d(pfwd,buf_fft,k)
      enddo  
    
      write(thdiag_id,'(f12.5,7g20.12)',advance='no') &
        time,                                         &
        mass,                                         &
        l1norm,                                       &
        momentum,                                     &
        l2norm,                                       &
        kinetic_energy,                               &
        potential_energy,                             &
        kinetic_energy+potential_energy
    
      do k=0,sim%nb_mode
        write(thdiag_id,'(g20.12)',advance='no') abs(rho_mode(k))
      enddo
    
      do k=0,sim%nb_mode-1
        write(thdiag_id,'(g20.12)',advance='no') f_hat_x2(k+1)
      enddo
    
      write(thdiag_id,'(g20.12)') f_hat_x2(sim%nb_mode+1)
    
      call sll_s_binary_write_array_1d(efield_id,efield(1:np_x1-1),ierr)
      call sll_s_binary_write_array_1d(rhotot_id,rho(1:np_x1-1),ierr)
      call sll_s_binary_write_array_0d(t_id,time,ierr)
      if (sim%driven) then
        call sll_s_binary_write_array_1d(Edr_id,e_app(1:np_x1-1),ierr)
        call sll_s_binary_write_array_0d(adr_id,adr,ierr)
      endif   

    endif
    
  end subroutine diagnostics
    
  subroutine write_init_files(sim, efield, rho, e_app)
    class(sll_t_simulation_2d_vlasov_ampere_cart), intent(inout) :: sim
    sll_real64 :: efield(:), rho(:), e_app(:)
    sll_int32 :: file_id
    sll_int32 :: np_x1
    sll_int32 :: np_x2
    sll_int32 :: ierr
    
    np_x1 = sim%mesh2d%num_cells1+1
    np_x2 = sim%mesh2d%num_cells2+1

    call sll_s_binary_file_create('f0.bdat', file_id, ierr)
    call sll_s_binary_write_array_2d(file_id,f_visu(1:np_x1-1,1:np_x2-1),ierr)
    call sll_s_binary_file_close(file_id,ierr)
    
    call sll_s_plot_f_cartesian( iplot,             &
                               f_visu,            &
                               sim%x1_array,      &
                               np_x1,             &
                               sim%node_positions_x2, &
                               sim%num_dof_x2,    &
                               'f',               &
                               sim%time_init )        
    
    print *,'#maxf',maxval(f_visu), minval(f_visu) 
    
    call sll_s_binary_file_create("x.bdat", file_id, ierr)
    call sll_s_binary_write_array_1d(file_id,sim%x1_array(1:np_x1-1),ierr)
    call sll_s_binary_file_close(file_id,ierr)                    
    call sll_s_binary_file_create("v.bdat", file_id, ierr)
    call sll_s_binary_write_array_1d(file_id,sim%node_positions_x2(1:np_x2-1),ierr)
    call sll_s_binary_file_close(file_id,ierr)
    
    call sll_s_ascii_file_create(sim%thdiag_filename, thdiag_id, ierr)
    
    call sll_s_binary_file_create('deltaf.bdat', deltaf_id, ierr)
    call sll_s_binary_file_create('rhotot.bdat', rhotot_id, ierr)
    call sll_s_binary_file_create('efield.bdat', efield_id, ierr)
    call sll_s_binary_file_create('t.bdat', t_id, ierr)
    call sll_s_binary_write_array_1d(efield_id,efield(1:np_x1-1),ierr)
    call sll_s_binary_write_array_1d(rhotot_id,rho(1:np_x1-1),ierr)
    call sll_s_binary_write_array_0d(t_id,real(istep,f64)*sim%dt,ierr)
    
    if (sim%driven) then
      call sll_s_binary_file_create('adr.bdat', adr_id, ierr)
      call sll_s_binary_file_create('Edr.bdat', Edr_id, ierr)
      call sll_s_binary_write_array_1d(Edr_id,e_app(1:np_x1-1),ierr)
      call sll_s_binary_write_array_0d(adr_id,adr,ierr)
    endif                    
    
  end subroutine write_init_files
    
  subroutine set_e_app(sim,t,e_app)
    class(sll_t_simulation_2d_vlasov_ampere_cart), intent(in) :: sim
    sll_real64 :: t
    sll_int32  :: i
    sll_int32  :: np_x1
    sll_real64, dimension(:) :: e_app
    
    np_x1 = sim%mesh2d%num_cells1+1
    call sll_s_pfenvelope(adr,                    &
                    t,                      &
                    sim%tflat,              &
                    sim%tL,                 &
                    sim%tR,                 &
                    sim%twL,                &
                    sim%twR,                &
                    sim%t0,                 &
                    sim%turn_drive_off)
    
    do i = 1, np_x1
      e_app(i) = sim%Edrmax*adr*sim%kx                          &
               * sin(sim%kx*real(i-1,f64)*sim%mesh2d%delta_eta1 &
               - sim%omegadr*t)
    enddo
    
  end subroutine set_e_app
    
  subroutine save_for_restart(layout_x1, f_x1)
    type(sll_t_layout_2d), pointer :: layout_x1
    sll_real64               :: f_x1(:,:)
    character(len=4) :: cplot
    sll_int32        :: iproc
    character(len=4) :: cproc
    sll_int32        :: restart_id
    sll_int32        :: ierr
    sll_int32        :: local_size_x1, local_size_x2

    call sll_o_compute_local_sizes( layout_x1, local_size_x1, local_size_x2 )

    iproc = sll_f_get_collective_rank(sll_v_world_collective)
    call sll_s_int2string(iproc, cproc)
    call sll_s_int2string(iplot, cplot) 

    call sll_s_binary_file_create('f_plot_'//cplot//'_proc_'//cproc//'.rst', &
                                restart_id, ierr )
    call sll_s_binary_write_array_0d(restart_id,time,ierr)
    call sll_s_binary_write_array_2d(restart_id, &
                                   f_x1(1:local_size_x1,1:local_size_x2),ierr)
    call sll_s_binary_file_close(restart_id,ierr)    
    
  end subroutine save_for_restart
    
  subroutine check_restart(sim, f_x1)

    class(sll_t_simulation_2d_vlasov_ampere_cart), intent(inout) :: sim
    sll_real64                  :: f_x1(:,:)
    character(len=*), parameter :: this_sub_name = 'check_restart'
    character(len=288)          :: err_msg
    sll_int32                   :: iproc
    character(len=4)            :: cproc
    logical                     :: file_exists
    sll_int32                   :: ierr
    sll_real64                  :: time_init
    sll_int32                   :: restart_id
    
    iproc = sll_f_get_collective_rank(sll_v_world_collective)
    call sll_s_int2string(iproc, cproc)

    if (trim(sim%restart_file) /= "no_restart_file" ) then
    
    
      INQUIRE(FILE=trim(sim%restart_file)//'_proc_'//cproc//'.rst', EXIST=file_exists)
    
      if (.not. file_exists) then
        err_msg = '#file '//trim(sim%restart_file)//'_proc_'//cproc//'.rst &
                  & does not exist'
        SLL_ERROR( this_sub_name, err_msg )
      endif
    
      open(unit=restart_id, &
           file=trim(sim%restart_file)//'_proc_'//cproc//'.rst', ACCESS="STREAM", &
           form='unformatted', IOStat=ierr)      
    
      if ( ierr .ne. 0 ) then
        err_msg = 'ERROR while opening file &
                   &'//trim(sim%restart_file)//'_proc_'//cproc//'.rst &
                   & Called from run_va2d_cartesian().'
        SLL_ERROR( this_sub_name, err_msg )
      end if
    
      print *,'#read restart file '//trim(sim%restart_file)//'_proc_'//cproc//'.rst'      
      call sll_s_binary_read_array_0d(restart_id,time_init,ierr)
      call sll_s_binary_read_array_2d(restart_id,f_x1,ierr)
      call sll_s_binary_file_close(restart_id,ierr)
    
    endif      
    
    if (sim%time_init_from_restart_file .eqv. .true.) then
      sim%time_init = time_init  
    endif
    
  end subroutine check_restart


  subroutine delete_va2d_par_cart( sim )

    class(sll_t_simulation_2d_vlasov_ampere_cart) :: sim
    sll_int32 :: ierr
    
    if(associated(sim%x1_array)) then
      SLL_DEALLOCATE(sim%x1_array,ierr)
      nullify(sim%x1_array)
    endif
    if(associated(sim%x2_array)) then
      SLL_DEALLOCATE(sim%x2_array,ierr)
      nullify(sim%x2_array)
    endif
        
  end subroutine delete_va2d_par_cart

  elemental function f_equilibrium(v)

    sll_real64, intent(in) :: v
    sll_real64             :: f_equilibrium

    f_equilibrium = 1.0_f64/sqrt(2*sll_p_pi)*exp(-0.5_f64*v*v)

  end function f_equilibrium

end module sll_m_sim_bsl_va_1d1v_cart_spectral
