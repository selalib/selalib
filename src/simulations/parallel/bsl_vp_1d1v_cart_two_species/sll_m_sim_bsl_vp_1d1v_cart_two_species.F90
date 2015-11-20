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

module sll_m_sim_bsl_vp_1d1v_cart_two_species

#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_errors.h"
use sll_m_collective
use sll_m_remapper
use sll_m_buffer_loader_utilities
use sll_m_constants
use sll_m_cartesian_meshes  
use sll_m_gnuplot
use sll_m_coordinate_transformation_2d_base
use sll_m_coordinate_transformations_2d
use sll_m_common_coordinate_transformations
use sll_m_common_array_initializers
use sll_m_parallel_array_initializer
use sll_m_advection_1d_periodic
use sll_m_advection_1d_non_uniform_cubic_splines
use sll_m_fft
use sll_m_sim_base
use sll_m_time_splitting_coeff
use sll_m_poisson_1d_periodic  
use sll_m_poisson_1d_periodic_solver
use sll_m_poisson_1d_polar_solver
use sll_m_advection_1d_ampere
use sll_m_primitives
use sll_m_xdmf
!#ifdef _OPENMP
!use omp_lib
!#endif


  implicit none

  integer, parameter :: SLL_ADVECTIVE = 0
  integer, parameter :: SLL_CONSERVATIVE = 1

  type, extends(sll_simulation_base_class) :: &
       sll_simulation_2d_vlasov_poisson_cart_two_species

   
   sll_int32 :: num_threads
   !geometry
   type(sll_cartesian_mesh_2d), pointer :: mesh2d_sp1
   type(sll_cartesian_mesh_2d), pointer :: mesh2d_sp2
   sll_int32 :: num_dof_x2_sp1
   sll_int32 :: num_dof_x2_sp2
   sll_real64, dimension(:), pointer :: x1_array
   sll_real64, dimension(:), pointer :: x2_array_sp1
   sll_real64, dimension(:), pointer :: x2_array_sp2
   sll_real64, dimension(:,:), pointer :: x2_array_omp_sp1
   sll_real64, dimension(:,:), pointer :: x2_array_omp_sp2
   sll_real64, dimension(:), pointer :: integration_weight_sp1
   sll_real64, dimension(:), pointer :: integration_weight_sp2
   sll_int32, dimension(:), pointer :: every_x1
   sll_int32, dimension(:), pointer :: every_x2
   sll_int32 :: num_bloc_x1
   sll_int32 :: num_bloc_x2_sp1
   sll_int32 :: num_bloc_x2_sp2
   sll_int32, dimension(:), pointer :: bloc_index_x1
   sll_int32, dimension(:), pointer :: bloc_index_x2_sp1
   sll_int32, dimension(:), pointer :: bloc_index_x2_sp2
      
   !initial function
   sll_real64  :: kx_sp1
   sll_real64  :: eps_sp1
   sll_real64  :: kx_sp2
   sll_real64  :: eps_sp2
   character(len=256) :: restart_file
   logical :: time_init_from_restart_file
   
   !initial function
   procedure(sll_scalar_initializer_2d), nopass, pointer :: init_func_sp1
   sll_real64, dimension(:), pointer :: params_sp1
   sll_real64 :: nrj0_sp1
   procedure(sll_scalar_initializer_2d), nopass, pointer :: init_func_sp2
   sll_real64, dimension(:), pointer :: params_sp2
   sll_real64 :: nrj0_sp2
   
   !time_iterations
   sll_real64 :: dt
   sll_int32 :: num_iterations
   sll_int32 :: freq_diag
   sll_int32 :: freq_diag_time
   sll_int32 :: freq_diag_restart
   sll_int32 :: nb_mode
   sll_real64 :: time_init
   !sll_int32  :: split_case
   type(splitting_coeff), pointer :: split
   character(len=256)      :: thdiag_filename
 
  
   !parameters for drive
   logical :: driven
   sll_real64 :: t0
   sll_real64 :: twL
   sll_real64 :: twR
   !sll_real64 :: tstart
   sll_real64 :: tflat
   sll_real64 :: tL
   sll_real64 :: tR
   sll_real64 :: Edrmax
   sll_real64 :: omegadr
   logical :: turn_drive_off

   !advector
   type(sll_advection_1d_base_ptr), dimension(:), pointer :: advect_x1_sp1
   type(sll_advection_1d_base_ptr), dimension(:), pointer :: advect_x1_sp2
   type(sll_advection_1d_base_ptr), dimension(:), pointer :: advect_x2_sp1
   type(sll_advection_1d_base_ptr), dimension(:), pointer :: advect_x2_sp2
   sll_int32 :: advection_form_x2_sp1
   sll_int32 :: advection_form_x2_sp2
   sll_real64 :: factor_x1
   sll_real64 :: factor_x2_rho
   sll_real64 :: factor_x2_1

   !poisson solver
   class(sll_poisson_1d_base), pointer   :: poisson
   sll_real64 :: mass_ratio
   sll_real64 :: alpha_sp1
   sll_real64 :: alpha_sp2
   !sll_real64, dimension(:), pointer :: mixt_bc
   
   contains
     procedure, pass(sim) :: run => run_vp2d_cartesian_two_species
     procedure, pass(sim) :: init_from_file => init_vp2d_fake !init_vp2d_par_cart_two_species
  end type sll_simulation_2d_vlasov_poisson_cart_two_species

  interface delete
     module procedure delete_vp2d_par_cart_two_species
  end interface delete

contains

  function new_vp2d_par_cart_two_species( &
    filename, &
    num_run) &
    result(sim)    
    type(sll_simulation_2d_vlasov_poisson_cart_two_species), pointer :: sim    
    character(len=*), intent(in), optional                                :: filename
    sll_int32, intent(in), optional                      :: num_run
    sll_int32 :: ierr   

    SLL_ALLOCATE(sim,ierr)            
    call init_vp2d_par_cart_two_species( &
      sim, &
      filename, &
      num_run)
       
  end function new_vp2d_par_cart_two_species
  
  
  subroutine change_initial_function_vp2d_par_cart_two_species( &
    sim, &
    init_func_sp1, &
    params_sp1, &
    num_params_sp1, &
    init_func_sp2, &
    params_sp2, &
    num_params_sp2)
    class(sll_simulation_2d_vlasov_poisson_cart_two_species), intent(inout) :: sim
    procedure(sll_scalar_initializer_2d), pointer :: init_func_sp1
    sll_real64, dimension(:), pointer :: params_sp1
    sll_int32, intent(in) :: num_params_sp1
    procedure(sll_scalar_initializer_2d), pointer :: init_func_sp2
    sll_real64, dimension(:), pointer :: params_sp2
    sll_int32, intent(in) :: num_params_sp2

    !local variables
    sll_int32 :: ierr
    sll_int32 :: i
    
    sim%init_func_sp1 => init_func_sp1

    if(associated(sim%params_sp1))then
      SLL_DEALLOCATE(sim%params_sp1,ierr)
    endif
    if(num_params_sp1<1)then
      print *,'#num_params should be >=1 in change_initial_function_vp2d_par_cart'
      stop
    endif
    SLL_ALLOCATE(sim%params_sp1(num_params_sp1),ierr)
    if(size(params_sp1)<num_params_sp1)then
      print *,'#size of params is not good in change_initial_function_vp2d_par_cart'
      stop
    endif
    do i=1,num_params_sp1
      sim%params_sp1(i) = params_sp1(i)
    enddo

    sim%init_func_sp2 => init_func_sp2
    if(associated(sim%params_sp2))then
      SLL_DEALLOCATE(sim%params_sp2,ierr)
    endif
    if(num_params_sp2<1)then
      print *,'#num_params should be >=1 in change_initial_function_vp2d_par_cart'
      stop
    endif
    SLL_ALLOCATE(sim%params_sp2(num_params_sp2),ierr)
    if(size(params_sp2)<num_params_sp2)then
      print *,'#size of params is not good in change_initial_function_vp2d_par_cart'
      stop
    endif
    do i=1,num_params_sp2
      sim%params_sp2(i) = params_sp2(i)
    enddo

  end subroutine change_initial_function_vp2d_par_cart_two_species

  subroutine init_vp2d_par_cart_two_species( sim, filename, num_run )
    class(sll_simulation_2d_vlasov_poisson_cart_two_species), intent(inout) :: sim
    character(len=*), intent(in), optional :: filename
    sll_int32, intent(in), optional :: num_run

    character(len=*), parameter :: this_sub_name = 'init_vp2d_par_cart_two_species'
    character(len=128)          :: err_msg


    !geometry
    character(len=256) :: mesh_case_x1
    sll_int32 :: num_cells_x1
    sll_real64 :: x1_min
    sll_real64 :: x1_max
    sll_int32 :: nbox_x1
    character(len=256) :: mesh_case_x2_sp1
    character(len=256) :: mesh_case_x2_sp2
    sll_int32 :: num_cells_x2_sp1
    sll_int32 :: num_cells_x2_sp2
    sll_real64 :: x2_min_sp1
    sll_real64 :: x2_max_sp1
    sll_real64 :: x2_min_sp2
    sll_real64 :: x2_max_sp2

    !physical_params
    sll_real64 :: mass_ratio
    sll_real64 :: vd
    sll_real64 :: alpha_sp1
    sll_real64 :: alpha_sp2
    
    !initial_function
    character(len=256) :: initial_function_case_sp1
    sll_real64 :: kmode_sp1
    sll_real64 :: eps_sp1
    sll_real64 :: sigma_sp1
    sll_real64 :: v0_sp1
    sll_real64 :: factor1_sp1
    sll_real64 :: alpha_gaussian_sp1
    character(len=256) :: initial_function_case_sp2
    sll_real64 :: kmode_sp2
    sll_real64 :: eps_sp2
    sll_real64 :: sigma_sp2
    sll_real64 :: v0_sp2
    sll_real64 :: factor1_sp2
    sll_real64 :: alpha_gaussian_sp2
    character(len=256) :: restart_file
    logical :: time_init_from_restart_file
    
    !time_iterations
    sll_real64 :: dt
    sll_int32 :: number_iterations
    sll_int32 :: freq_diag
    sll_int32 :: freq_diag_restart
    sll_int32 :: freq_diag_time
    sll_int32 :: nb_mode
    sll_real64 :: time_init
    character(len=256) :: split_case

    !advector
    character(len=256) :: advector_x1_sp1
    sll_int32 :: order_x1_sp1
    character(len=256) :: advector_x1_sp2
    sll_int32 :: order_x1_sp2
    character(len=256) :: advector_x2_sp1
    sll_int32 :: order_x2_sp1
    character(len=256) :: advector_x2_sp2
    sll_int32 :: order_x2_sp2
    character(len=256) :: advection_form_x2_sp1
    character(len=256) :: advection_form_x2_sp2
    character(len=256) :: integration_case
    sll_real64 :: factor_x1
    sll_real64 :: factor_x2_rho
    sll_real64 :: factor_x2_1

    !poisson
    character(len=256) :: poisson_solver
    
    !drive
    character(len=256) :: drive_type
    sll_real64 :: keen_t0
    sll_real64 :: keen_tL
    sll_real64 :: keen_tR
    sll_real64 :: keen_twL
    sll_real64 :: keen_twR
    sll_real64 :: keen_tflat
    logical :: keen_turn_drive_off
    sll_real64 :: keen_Edrmax
    sll_real64 :: keen_omegadr
    
    
    !local variables
    intrinsic :: trim
    sll_int32             :: IO_stat
    sll_int32                            :: input_file
    type(sll_cartesian_mesh_1d), pointer :: mesh_x1
    type(sll_cartesian_mesh_1d), pointer :: mesh_x2_sp1
    type(sll_cartesian_mesh_1d), pointer :: mesh_x2_sp2
    sll_int32 :: ierr
    sll_int32, parameter  :: param_out = 37, param_out_drive = 40
    !sll_real64 :: bloc_coord(2)
    !sll_int32 :: bloc_index(3)
    sll_int32 :: i
    sll_int32 :: num_threads
    sll_int32 :: tid
    character(len=256)      :: str_num_run
    character(len=256)      :: filename_loc
  
    ! namelists for data input
    namelist /geometry/ &
      mesh_case_x1, &
      num_cells_x1, &
      x1_min, &
      x1_max, &
      nbox_x1, &
      mesh_case_x2_sp1, &
      mesh_case_x2_sp2, &
      num_cells_x2_sp1, &
      num_cells_x2_sp2, &
      x2_min_sp1, &
      x2_max_sp1, &
      x2_min_sp2, &
      x2_max_sp2

    namelist /physical_params/ &
         mass_ratio, &
         vd, &
         alpha_sp1, &
         alpha_sp2
    
    namelist /initial_function/ &
      initial_function_case_sp1, &
      kmode_sp1, &
      eps_sp1, &
      sigma_sp1, &
      v0_sp1, &
      factor1_sp1, &
      alpha_gaussian_sp1, &
      initial_function_case_sp2, &
      kmode_sp2, &
      eps_sp2, &
      sigma_sp2, &
      v0_sp2, &
      factor1_sp2, &
      alpha_gaussian_sp2, &
      restart_file, &
      time_init_from_restart_file

    namelist /time_iterations/ &
      dt, &
      number_iterations, &
      freq_diag, &
      freq_diag_time, &
      freq_diag_restart, &
      nb_mode, &
      time_init, &
      split_case

    namelist /advector/ &
      advector_x1_sp1, &
      order_x1_sp1, &
      advector_x1_sp2, &
      order_x1_sp2, &
      advector_x2_sp1, &
      order_x2_sp1, &
      advector_x2_sp2, &
      order_x2_sp2, &
      advection_form_x2_sp1, &
      advection_form_x2_sp2, &
      factor_x1, &
      factor_x2_rho, &
      factor_x2_1, &
      integration_case

    namelist /poisson/ &
      poisson_solver

    namelist /drive/ &
      drive_type, &
      keen_t0, &
      keen_twL, &
      keen_twR, &
      !keen_tstart, &
      keen_tflat, &
      keen_tL, &
      keen_tR, &
      keen_turn_drive_off, &
      keen_Edrmax, &
      keen_omegadr

    num_threads = 1

!#ifdef _OPENMP
!!$OMP PARALLEL SHARED(num_threads)
!    if(omp_get_thread_num()==0)then
!      num_threads =  omp_get_num_threads()      
!    endif
!!$OMP END PARALLEL
!#endif
   sim%num_threads = num_threads
   print *,'#num_threads=',num_threads

    !! set default parameters

    !geometry
    mesh_case_x1 = "SLL_LANDAU_MESH"
    num_cells_x1 = 32
    x1_min = 0.0_f64
    nbox_x1 = 1
    mesh_case_x2_sp1 = "SLL_CARTESIAN_MESH"
    num_cells_x2_sp1 = 64
    x2_min_sp1 = -6._f64
    x2_max_sp1 = 6._f64
    mesh_case_x2_sp2 = "SLL_CARTESIAN_MESH"
    num_cells_x2_sp2 = 64
    x2_min_sp2 = -0.5_f64
    x2_max_sp2 = 0.5_f64

    !physical_params
    mass_ratio = 0.0005_f64
    vd = 0._f64
    alpha_sp1 = 1._f64
    alpha_sp2 = 1._f64
        
    !initial_function
    initial_function_case_sp1 = "SLL_LANDAU"
    kmode_sp1 = 0.5_f64
    eps_sp1 = 0.001_f64
    sigma_sp1 = 1._f64
    v0_sp1 = 0._f64
    factor1_sp1 = 1._f64/sqrt(2._f64*sll_pi)
    !initial_function_case = "SLL_BEAM"
    alpha_gaussian_sp1 = 0.2_f64
    !initial_function_case = "SLL_PLASMA_SHEATH"
    !sigma_sp1 = 5.52d-4
    !v0_sp1 = 0._f64
    
    initial_function_case_sp2 = "SLL_LANDAU"
    kmode_sp2 = 0.5_f64
    eps_sp2 = 0.001_f64
    sigma_sp2 = 1._f64
    v0_sp2 = 0._f64
    factor1_sp2 = 1._f64/sqrt(2._f64*sll_pi)
    !initial_function_case = "SLL_BEAM"
    alpha_gaussian_sp2 = 0.2_f64
    !initial_function_case = "SLL_PLASMA_SHEATH"
    !sigma_sp2 = 5.52d-4
    !v0_sp2 = 0._f64
    
    restart_file = "no_restart_file"
    time_init_from_restart_file = .false.
    
    !time_iterations
    dt = 0.1_f64
    number_iterations = 600
    freq_diag = 100
    freq_diag_time = 1
    freq_diag_restart = 5000
    nb_mode = 5
    time_init = 0._f64
    split_case = "SLL_STRANG_VTV" 
    !split_case = "SLL_STRANG_TVT" 
    !split_case = "SLL_ORDER6VPnew1_VTV" 
    !split_case = "SLL_ORDER6VPnew2_VTV" 
    !split_case = "SLL_ORDER6_VTV"
    !split_case = "SLL_LIE_TV"

    !advector
    advector_x1_sp1 = "SLL_LAGRANGE"
    order_x1_sp1 = 4
    advector_x1_sp2 = "SLL_LAGRANGE"
    order_x1_sp2 = 4
    advector_x2_sp1 = "SLL_LAGRANGE"
    order_x2_sp1 = 4
    advector_x2_sp2 = "SLL_LAGRANGE"
    order_x2_sp2 = 4
    advection_form_x2_sp1 = "SLL_ADVECTIVE"
    advection_form_x2_sp2 = "SLL_ADVECTIVE"
    factor_x1 = 1._f64
    factor_x2_rho = 1._f64
    factor_x2_1 = 1._f64

    !integration_case = "SLL_RECTANGLE"
    integration_case = "SLL_TRAPEZOID"
    
    !poisson
    poisson_solver = "SLL_FFT"
    
    !drive
    drive_type = "SLL_NO_DRIVE"
    !drive_type = "SLL_KEEN_DRIVE"  
      !keen_t0 = 0.
      !keen_tL = 69.
      !keen_tR = 307.
      !keen_twL = 20.
      !keen_twR = 20.
      !keen_tflat = 100.
      !keen_turn_drive_off = .true.
      !keen_Edrmax = 0.2
      !keen_omegadr = 0.37	

    if(present(num_run))then
      !call int2string(num_run, str_num_run)
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

      call sll_new_file_id(input_file, ierr)

      open(unit = input_file, file=trim(filename_loc)//'.nml',IOStat=IO_stat)
      if ( IO_stat /= 0 ) then
        err_msg = 'failed to open file '//trim(filename_loc)//'.nml'
        SLL_ERROR( this_sub_name, err_msg )
      end if


      if(sll_get_collective_rank(sll_world_collective)==0)then
        print *,'#initialization with filename:'
        print *,'#',trim(filename_loc)//'.nml'
      endif
      read(input_file, geometry)
      read(input_file, physical_params)
      read(input_file, initial_function)
      read(input_file, time_iterations)
      read(input_file, advector)
      read(input_file, poisson)
      read(input_file, drive)
      close(input_file)
    else
      if(sll_get_collective_rank(sll_world_collective)==0)then
        print *,'#initialization with default parameters'
      endif      
    endif


    !geometry
    
    select case (mesh_case_x1)
      case ("SLL_LANDAU_MESH")
        x1_max = real(nbox_x1,f64) * 2._f64 * sll_pi / kmode_sp1
        mesh_x1 => new_cartesian_mesh_1d(num_cells_x1,eta_min=x1_min, eta_max=x1_max)
        call get_node_positions( mesh_x1, sim%x1_array )
      case ("SLL_CARTESIAN_MESH")
        mesh_x1 => new_cartesian_mesh_1d(num_cells_x1,eta_min=x1_min, eta_max=x1_max)  
        call get_node_positions( mesh_x1, sim%x1_array )
      case default
        print*,'#mesh_case_x1', mesh_case_x1, ' not implemented'
        print*,'#in init_vp2d_par_cart'
        stop 
    end select
    
    select case (mesh_case_x2_sp1)
      case ("SLL_CARTESIAN_MESH")
        mesh_x2_sp1 => new_cartesian_mesh_1d(num_cells_x2_sp1,eta_min=x2_min_sp1, eta_max=x2_max_sp1)
        call get_node_positions( mesh_x2_sp1, sim%x2_array_sp1 )
        SLL_ALLOCATE(sim%x2_array_omp_sp1(num_cells_x2_sp1+1,0:sim%num_threads-1),ierr)
        do i=0,sim%num_threads-1
          sim%x2_array_omp_sp1(:,i) = sim%x2_array_sp1(:)
        enddo
      case default
        print*,'#mesh_case_x2', mesh_case_x2_sp1, ' not implemented'
        print*,'#in init_vp2d_par_cart'
        stop 
    end select
    sim%mesh2d_sp1 => tensor_product_1d_1d( mesh_x1, mesh_x2_sp1)

    select case (mesh_case_x2_sp2)
      case ("SLL_CARTESIAN_MESH")
        mesh_x2_sp2 => new_cartesian_mesh_1d(num_cells_x2_sp2,eta_min=x2_min_sp2, eta_max=x2_max_sp2)
        call get_node_positions( mesh_x2_sp2, sim%x2_array_sp2 )
        SLL_ALLOCATE(sim%x2_array_omp_sp2(num_cells_x2_sp2+1,0:sim%num_threads-1),ierr)
        do i=0,sim%num_threads-1
          sim%x2_array_omp_sp2(:,i) = sim%x2_array_sp2(:)
        enddo
      case default
        print*,'#mesh_case_x2', mesh_case_x2_sp2, ' not implemented'
        print*,'#in init_vp2d_par_cart'
        stop 
    end select
    sim%mesh2d_sp2 => tensor_product_1d_1d( mesh_x1, mesh_x2_sp2)



    
    !initial function
    sim%nrj0_sp1 = 0._f64
    sim%kx_sp1 = kmode_sp1
    sim%eps_sp1 = eps_sp1
    select case (initial_function_case_sp1)
    case ("SLL_LANDAU")
       sim%init_func_sp1 => sll_landau_initializer_2d
       SLL_ALLOCATE(sim%params_sp1(4),ierr)
       sim%params_sp1(1) = kmode_sp1
       sim%params_sp1(2) = eps_sp1
       sim%params_sp1(3) = v0_sp1
       sim%params_sp1(4) = sigma_sp1
       sim%nrj0_sp1 = 0._f64  !compute the right value
       !(0.5_f64*eps*sll_pi)**2/(kmode_x1*kmode_x2) &
       !*(1._f64/kmode_x1**2+1._f64/kmode_x2**2)
       !for the moment
       sim%kx_sp1 = kmode_sp1
       sim%eps_sp1 = eps_sp1
    case ("SLL_BUMP_ON_TAIL")
       sim%init_func_sp1 => sll_bump_on_tail_initializer_2d
       SLL_ALLOCATE(sim%params_sp1(2),ierr)
       sim%params_sp1(1) = kmode_sp1
       sim%params_sp1(2) = eps_sp1
       sim%nrj0_sp1 = 0._f64  !compute the right value
       !(0.5_f64*eps*sll_pi)**2/(kmode_x1*kmode_x2) &
       !*(1._f64/kmode_x1**2+1._f64/kmode_x2**2)
       !for the moment
       sim%kx_sp1 = kmode_sp1
       sim%eps_sp1 = eps_sp1
    case ("SLL_TWO_STREAM_INSTABILITY")
       sim%init_func_sp1 => sll_two_stream_instability_initializer_2d
       SLL_ALLOCATE(sim%params_sp1(4),ierr)
       sim%params_sp1(1) = kmode_sp1
       sim%params_sp1(2) = eps_sp1
       sim%params_sp1(3) = sigma_sp1
       sim%params_sp1(4) = factor1_sp1
       sim%nrj0_sp1 = 0._f64  !compute the right value
       !(0.5_f64*eps*sll_pi)**2/(kmode_x1*kmode_x2) &
       !*(1._f64/kmode_x1**2+1._f64/kmode_x2**2)
       !for the moment
       sim%kx_sp1 = kmode_sp1
       sim%eps_sp1 = eps_sp1
    case ("SLL_BEAM")  
       sim%init_func_sp1 => sll_beam_initializer_2d
       SLL_ALLOCATE(sim%params_sp1(1),ierr)
       sim%params_sp1(1) = alpha_gaussian_sp1
!    case ("SLL_PLASMA_SHEATH")
!       sim%init_func_sp1 => sll_plasma_sheath_initializer_2d
!       SLL_ALLOCATE(sim%params_sp1(2),ierr)
!       sim%params_sp1(1) = theta_sp1
!       sim%params_sp1(2) = v0_sp1
    case default
        print *,'#init_func_case not implemented'
        print *,'#in init_vp2d_par_cart'  
        stop
    end select

    !initial function
    sim%nrj0_sp2 = 0._f64
    sim%kx_sp2 = kmode_sp2
    sim%eps_sp2 = eps_sp2
    select case (initial_function_case_sp2)
    case ("SLL_LANDAU")
       sim%init_func_sp2 => sll_landau_initializer_2d
       SLL_ALLOCATE(sim%params_sp2(4),ierr)
       sim%params_sp2(1) = kmode_sp2
       sim%params_sp2(2) = eps_sp2
       sim%params_sp2(3) = v0_sp2
       sim%params_sp2(4) = sigma_sp2
       sim%nrj0_sp2 = 0._f64  !compute the right value
       !(0.5_f64*eps*sll_pi)**2/(kmode_x1*kmode_x2) &
       !*(1._f64/kmode_x1**2+1._f64/kmode_x2**2)
       !for the moment
       sim%kx_sp2 = kmode_sp2
       sim%eps_sp2 = eps_sp2
    case ("SLL_BUMP_ON_TAIL")
       sim%init_func_sp2 => sll_bump_on_tail_initializer_2d
       SLL_ALLOCATE(sim%params_sp2(2),ierr)
       sim%params_sp2(1) = kmode_sp2
       sim%params_sp2(2) = eps_sp2
       sim%nrj0_sp2 = 0._f64  !compute the right value
       !(0.5_f64*eps*sll_pi)**2/(kmode_x1*kmode_x2) &
       !*(1._f64/kmode_x1**2+1._f64/kmode_x2**2)
       !for the moment
       sim%kx_sp2 = kmode_sp2
       sim%eps_sp2 = eps_sp2
    case ("SLL_TWO_STREAM_INSTABILITY")
       sim%init_func_sp2 => sll_two_stream_instability_initializer_2d
       SLL_ALLOCATE(sim%params_sp2(4),ierr)
       sim%params_sp2(1) = kmode_sp2
       sim%params_sp2(2) = eps_sp2
       sim%params_sp2(3) = sigma_sp2
       sim%params_sp2(4) = factor1_sp2
       sim%nrj0_sp2 = 0._f64  !compute the right value
       !(0.5_f64*eps*sll_pi)**2/(kmode_x1*kmode_x2) &
       !*(1._f64/kmode_x1**2+1._f64/kmode_x2**2)
       !for the moment
       sim%kx_sp2 = kmode_sp2
       sim%eps_sp2 = eps_sp2
    case ("SLL_BEAM")  
       sim%init_func_sp2 => sll_beam_initializer_2d
       SLL_ALLOCATE(sim%params_sp2(1),ierr)
       sim%params_sp2(1) = alpha_gaussian_sp2
!    case ("SLL_PLASMA_SHEATH")
!       sim%init_func_sp2 => sll_plasma_sheath_initializer_2d
!       SLL_ALLOCATE(sim%params_sp2(2),ierr)
!       sim%params_sp2(1) = theta_sp2
!       sim%params_sp2(2) = v0_sp2
!    case ("SLL_STATIONARY_IONS")
!       sim%init_func_sp2 => sll_stationary_ions_initializer_2d
!       SLL_ALLOCATE(sim%params_sp2(3),ierr)
!       sim%params_sp2(1) = theta_sp2
!       sim%params_sp2(2) = mass_ratio
!       sim%params_sp2(3) = 0
    case default
       print *,'#init_func_case not implemented'
       print *,'#in init_vp2d_par_cart'  
       stop
    end select

    sim%time_init_from_restart_file = time_init_from_restart_file
    sim%restart_file = restart_file
    
    !time iterations
    sim%dt=dt
    sim%num_iterations=number_iterations
    sim%freq_diag=freq_diag
    sim%freq_diag_restart=freq_diag_restart
    sim%freq_diag_time=freq_diag_time
    sim%nb_mode = nb_mode
    sim%time_init = time_init
    
    if(sim%nb_mode<0)then
      print *,'#bad value of nb_mode=',nb_mode      
      print *,'#should be >=0'
      print *,'#in init_vp2d_par_cart'
      stop      
    endif
    
    select case (split_case)    
      case ("SLL_LIE_TV")
        sim%split => new_time_splitting_coeff(SLL_LIE_TV)
      case ("SLL_LIE_VT") 
        sim%split => new_time_splitting_coeff(SLL_LIE_VT)
      case ("SLL_STRANG_TVT") 
        sim%split => new_time_splitting_coeff(SLL_STRANG_TVT)
      case ("SLL_STRANG_VTV") 
        sim%split => new_time_splitting_coeff(SLL_STRANG_VTV)
      case ("SLL_TRIPLE_JUMP_TVT") 
        sim%split => new_time_splitting_coeff(SLL_TRIPLE_JUMP_TVT)
      case ("SLL_TRIPLE_JUMP_VTV") 
        sim%split => new_time_splitting_coeff(SLL_TRIPLE_JUMP_VTV)
      case ("SLL_ORDER6_VTV") 
        sim%split => new_time_splitting_coeff(SLL_ORDER6_VTV)
      case ("SLL_ORDER6VP_TVT") 
        sim%split => new_time_splitting_coeff(SLL_ORDER6VP_TVT,dt=dt)
      case ("SLL_ORDER6VP_VTV") 
        sim%split => new_time_splitting_coeff(SLL_ORDER6VP_VTV,dt=dt)
      case ("SLL_ORDER6VPnew_TVT") 
        sim%split => new_time_splitting_coeff(SLL_ORDER6VPnew_TVT,dt=dt)
      case ("SLL_ORDER6VPnew1_VTV") 
        sim%split => new_time_splitting_coeff(SLL_ORDER6VPnew1_VTV,dt=dt)
      case ("SLL_ORDER6VPnew2_VTV") 
        sim%split => new_time_splitting_coeff(SLL_ORDER6VPnew2_VTV,dt=dt)
      case default
        print *,'#split_case not defined'
        print *,'#in initialize_vlasov_par_poisson_seq_cart'
        stop       
    end select

    !advector
    SLL_ALLOCATE(sim%advect_x1_sp1(num_threads),ierr)
    SLL_ALLOCATE(sim%advect_x2_sp1(num_threads),ierr)
    SLL_ALLOCATE(sim%advect_x1_sp2(num_threads),ierr)
    SLL_ALLOCATE(sim%advect_x2_sp2(num_threads),ierr)

    tid = 1
!#ifdef _OPENMP
!!$OMP PARALLEL DEFAULT(SHARED) &
!!$OMP PRIVATE(tid)
!    tid = omp_get_thread_num()+1
!#endif
    select case (advector_x1_sp1)
      case ("SLL_SPLINES") ! arbitrary order periodic splines
        sim%advect_x1_sp1(tid)%ptr => new_periodic_1d_advector( &
          num_cells_x1, &
          x1_min, &
          x1_max, &
          SPLINE, & 
          order_x1_sp1) 

     case("SLL_LAGRANGE") ! arbitrary order Lagrange periodic interpolation
        sim%advect_x1_sp1(tid)%ptr => new_periodic_1d_advector( &
          num_cells_x1, &
          x1_min, &
          x1_max, &
          LAGRANGE, & 
          order_x1_sp1)

!     case("SLL_BSL_0INFLOW") ! BSL advection with no inflow BC
!        sim%advect_x1_sp1(tid)%ptr => new_BSL_1d_advector_0inflow(&
!             num_cells_x1+1, &
!             x1_min, &
!             x1_max, &
!             SLL_HERMITE, &
!             SLL_SET_TO_LIMIT)
!
!     case("SLL_SOURCE_COLLECTOR") !BSL advection with collector at xmin and source at xmax
!        sim%advect_x1_sp1(tid)%ptr => new_BSL_1d_advector_sc(&
!             num_cells_x1+1, &
!             x1_min, &
!             x1_max, &
!             SLL_HERMITE, &
!             SLL_SET_TO_LIMIT, &
!             vd, &
!             1._f64)
        
     case default
        print*,'#advector in x1 for species 1 ', advector_x1_sp1, ' not implemented'
        stop 
     end select
     
     select case (advector_x1_sp2)
      case ("SLL_SPLINES") ! arbitrary order periodic splines
        sim%advect_x1_sp2(tid)%ptr => new_periodic_1d_advector( &
          num_cells_x1, &
          x1_min, &
          x1_max, &
          SPLINE, & 
          order_x1_sp2) 

     case("SLL_LAGRANGE") ! arbitrary order Lagrange periodic interpolation
        sim%advect_x1_sp2(tid)%ptr => new_periodic_1d_advector( &
          num_cells_x1, &
          x1_min, &
          x1_max, &
          LAGRANGE, & 
          order_x1_sp2)

!     case("SLL_BSL_0INFLOW") ! BSL advection with no inflow BC
!        sim%advect_x1_sp2(tid)%ptr => new_BSL_1d_advector_0inflow(&
!             num_cells_x1+1, &
!             x1_min, &
!             x1_max, &
!             SLL_HERMITE, &
!             SLL_SET_TO_LIMIT)
!
!     case("SLL_SOURCE_COLLECTOR") !BSL advection with collector at xmin and source at xmax
!        sim%advect_x1_sp2(tid)%ptr => new_BSL_1d_advector_sc(&
!             num_cells_x1+1, &
!             x1_min, &
!             x1_max, &
!             SLL_HERMITE, &
!             SLL_SET_TO_LIMIT, &
!             vd, &
!             mass_ratio)
!
      case default
        print*,'#advector in x1 for species 2 ', advector_x1_sp2, ' not implemented'
        stop 
    end select    

    select case (advector_x2_sp2)
      case ("SLL_SPLINES") ! arbitrary order periodic splines
        sim%advect_x2_sp2(tid)%ptr => new_periodic_1d_advector( &
          num_cells_x2_sp2, &
          x2_min_sp2, &
          x2_max_sp2, &
          SPLINE, & 
          order_x2_sp2) 
      case("SLL_LAGRANGE") ! arbitrary order Lagrange periodic interpolation
        sim%advect_x2_sp2(tid)%ptr => new_periodic_1d_advector( &
          num_cells_x2_sp2, &
          x2_min_sp2, &
          x2_max_sp2, &
          LAGRANGE, & 
          order_x2_sp2)
      case("SLL_NON_UNIFORM_CUBIC_SPLINES") ! arbitrary order Lagrange periodic interpolation
        sim%advect_x2_sp2(tid)%ptr => new_non_uniform_cubic_splines_1d_advector( &
          num_cells_x2_sp2, &
          x2_min_sp2, &
          x2_max_sp2, &
          order_x2_sp2, &
          sim%x2_array_sp2)           
          !sim%x2_array_omp(:,tid))           

     case default
        print*,'#advector in x2 for species 2', advector_x2_sp2, ' not implemented'
        stop 
     end select
     
     select case (advector_x2_sp1)
      case ("SLL_SPLINES") ! arbitrary order periodic splines
        sim%advect_x2_sp1(tid)%ptr => new_periodic_1d_advector( &
          num_cells_x2_sp1, &
          x2_min_sp1, &
          x2_max_sp1, &
          SPLINE, & 
          order_x2_sp1) 
      case("SLL_LAGRANGE") ! arbitrary order Lagrange periodic interpolation
        sim%advect_x2_sp1(tid)%ptr => new_periodic_1d_advector( &
          num_cells_x2_sp1, &
          x2_min_sp1, &
          x2_max_sp1, &
          LAGRANGE, & 
          order_x2_sp1)
      case("SLL_NON_UNIFORM_CUBIC_SPLINES") ! arbitrary order Lagrange periodic interpolation
        sim%advect_x2_sp1(tid)%ptr => new_non_uniform_cubic_splines_1d_advector( &
          num_cells_x2_sp1, &
          x2_min_sp1, &
          x2_max_sp1, &
          order_x2_sp1, &
          sim%x2_array_sp1)           
          !sim%x2_array_omp(:,tid))           

     case default
        print*,'#advector in x2 for species 1', advector_x2_sp1, ' not implemented'
        stop 
    end select

!#ifdef _OPENMP
!!$OMP END PARALLEL
!#endif
    select case (advection_form_x2_sp1)
      case ("SLL_ADVECTIVE")
        sim%advection_form_x2_sp1 = SLL_ADVECTIVE
        sim%num_dof_x2_sp1 = num_cells_x2_sp1+1
      case ("SLL_CONSERVATIVE")
        sim%advection_form_x2_sp1 = SLL_CONSERVATIVE
        sim%num_dof_x2_sp1 = num_cells_x2_sp1
      case default
        print*,'#advection_form_x2_sp1', advection_form_x2_sp1, ' not implemented'
        print *,'#in init_vp2d_par_cart'
        stop 
    end select  

    select case (advection_form_x2_sp2)
      case ("SLL_ADVECTIVE")
        sim%advection_form_x2_sp2 = SLL_ADVECTIVE
        sim%num_dof_x2_sp2 = num_cells_x2_sp2+1
      case ("SLL_CONSERVATIVE")
        sim%advection_form_x2_sp2 = SLL_CONSERVATIVE
        sim%num_dof_x2_sp2 = num_cells_x2_sp2
      case default
        print*,'#advection_form_x2', advection_form_x2_sp2, ' not implemented'
        print *,'#in init_vp2d_par_cart'
        stop 
    end select  

    sim%factor_x1 = factor_x1
    sim%factor_x2_rho = factor_x2_rho
    sim%factor_x2_1 = factor_x2_1
    
    
    
    SLL_ALLOCATE(sim%integration_weight_sp1(sim%num_dof_x2_sp1),ierr)
    SLL_ALLOCATE(sim%integration_weight_sp2(sim%num_dof_x2_sp2),ierr)
    select case (integration_case)
      case ("SLL_RECTANGLE")
        do i=1,num_cells_x2_sp1
          sim%integration_weight_sp1(i) = sim%x2_array_sp1(i+1)-sim%x2_array_sp1(i)
        enddo
        sim%integration_weight_sp1(num_cells_x2_sp1+1) = 0._f64   

        do i=1,num_cells_x2_sp2
          sim%integration_weight_sp2(i) = sim%x2_array_sp2(i+1)-sim%x2_array_sp2(i)
        enddo
        sim%integration_weight_sp2(num_cells_x2_sp2+1) = 0._f64
        
      case ("SLL_TRAPEZOID")
        sim%integration_weight_sp1(1)=0.5_f64*(sim%x2_array_sp1(2)-sim%x2_array_sp1(1))
        do i=2,num_cells_x2_sp1
          sim%integration_weight_sp1(i) = 0.5_f64*(sim%x2_array_sp1(i+1)-sim%x2_array_sp1(i-1))
        enddo  
        sim%integration_weight_sp1(num_cells_x2_sp1+1) = &
          0.5_f64*(sim%x2_array_sp1(num_cells_x2_sp1+1)-sim%x2_array_sp1(num_cells_x2_sp1))

        sim%integration_weight_sp2(1)=0.5_f64*(sim%x2_array_sp2(2)-sim%x2_array_sp2(1))
        do i=2,num_cells_x2_sp2
          sim%integration_weight_sp2(i) = 0.5_f64*(sim%x2_array_sp2(i+1)-sim%x2_array_sp2(i-1))
        enddo  
        sim%integration_weight_sp2(num_cells_x2_sp2+1) = &
          0.5_f64*(sim%x2_array_sp2(num_cells_x2_sp2+1)-sim%x2_array_sp2(num_cells_x2_sp2))

      case ("SLL_CONSERVATIVE")
        do i=1,num_cells_x2_sp1
          sim%integration_weight_sp1(i)=sim%x2_array_sp1(i+1)-sim%x2_array_sp1(i)
        enddo  

        do i=1,num_cells_x2_sp2
          sim%integration_weight_sp2(i)=sim%x2_array_sp2(i+10)-sim%x2_array_sp2(i)
        enddo  

      case default
        print *,'#integration_case not implemented'
        print *,'#in init_vp2d_par_cart'  
        stop      
    end select  

    sim%mass_ratio = mass_ratio
    sim%alpha_sp1 = alpha_sp1    
    sim%alpha_sp2 = alpha_sp2
        
    !poisson
    !SLL_ALLOCATE(sim%mixt_bc(2),ierr)
    select case (poisson_solver)
    case ("SLL_FFT")
       sim%poisson => new_poisson_1d_periodic_solver( &
            x1_min, &
            x1_max, &
            num_cells_x1)
    case ("SLL_POLAR")
       sim%poisson => new_poisson_1d_polar_solver( &
            x1_min, &
            x1_max, &
            num_cells_x1)
!    case ("SLL_DIRICHLET")
!       sim%poisson => new_poisson_1d_dirichlet_solver( &
!            x1_min, &
!            x1_max, &
!            num_cells_x1)
!    case ("SLL_MIXT")
!       sim%poisson => new_poisson_1d_mixt_solver( &
!            x1_min, &
!            x1_max, &
!            num_cells_x1, &
!            sim%mixt_bc)
!#ifdef PASTIX
!    case ("SLL_FML_DIRICHLET")
!       sim%poisson => new_poisson_1d_fml_solver( &
!            x1_min, &
!            x1_max, &
!            num_cells_x1, &
!            'dirichlet')
!    case ("SLL_FML_MIXT")
!       sim%poisson => new_poisson_1d_fml_solver( &
!            x1_min, &
!            x1_max, &
!            num_cells_x1, &
!            'mixt', &
!            sim%mixt_bc)
!#endif
    case default
       print*,'#poisson_solver', poisson_solver, ' not implemented'
       print *,'#in init_vp2d_par_cart'
       stop 
    end select
    
    

    !drive
    select case (drive_type)
      case ("SLL_NO_DRIVE")
        sim%driven = .false.
      case("SLL_KEEN_DRIVE")
        sim%driven = .true.
        sim%t0=keen_t0
        sim%twL=keen_twL
        sim%twR=keen_twR
        !sim%tstart=keen_tstart
        sim%tflat=keen_tflat
        sim%tL=keen_tL
        sim%tR=keen_tR
        sim%turn_drive_off=keen_turn_drive_off
        sim%Edrmax=keen_Edrmax
        sim%omegadr=keen_omegadr        
      case default
        print*,'#drive_type', drive_type, ' not implemented'
        stop 
    end select




    


    
     
!!!--- ### ---!!!    
    if(sll_get_collective_rank(sll_world_collective)==0)then
            
      !to be compatible with VPpostprocessing_drive_KEEN.f
      if(sim%driven) then     
      open(unit = param_out, file = 'parameters.dat')
      
      
        write(param_out, *) real(number_iterations,f64)*dt !Tmax
        write(param_out, *) dt                  !dt
        write(param_out, *) number_iterations              !Nt
        write(param_out, *) kmode_sp1               !k0
        write(param_out, *) sim%omegadr             !omega0
        write(param_out, *) x1_max-x1_min           !L
        write(param_out, *) nbox_x1                !mbox
        write(param_out, *) sim%Edrmax              !Edrmax
        write(param_out, *) freq_diag_time       !nsave
        write(param_out, *) freq_diag            !nsavef1
        write(param_out, *) freq_diag            !nsavef2
        write(param_out, *) sim%tR+sim%tflat            !Tsetup
        write(param_out, *) x2_max_sp1                !vxmax (sp1)
        write(param_out, *) x2_min_sp1                !vxmin (sp1)
        write(param_out, *) x2_max_sp2                !vxmax (sp2)
        write(param_out, *) x2_min_sp2                !vxmin (sp2)
        write(param_out, *) num_cells_x1                 !N
        write(param_out, *) num_cells_x2_sp1                 !Nv
        write(param_out, *) num_cells_x2_sp2                 !Nv
        ! added ES
        write(param_out, *) sim%tL                  !tL
        write(param_out, *) sim%tR                  !tR
        write(param_out, *) sim%twL                 !twL
        write(param_out, *) sim%twR                 !twR
        write(param_out, *) sim%tflat               !tflat
      
      
      close(param_out)
      endif  
      

    endif
  end subroutine init_vp2d_par_cart_two_species


  subroutine init_vp2d_fake(sim, filename)
    class(sll_simulation_2d_vlasov_poisson_cart_two_species), intent(inout) :: sim
    character(len=*), intent(in)                                :: filename
  
    print *,'# Do not use the routine init_vp2d_fake'
    print *,'#use instead init_vp2d_par_cart'
    print *,sim%dt
    print *,filename
    stop
  
  end subroutine init_vp2d_fake



  subroutine run_vp2d_cartesian_two_species(sim)
    class(sll_simulation_2d_vlasov_poisson_cart_two_species), intent(inout) :: sim
    sll_real64,dimension(:,:),pointer :: f_x1_sp1,f_x2_sp1,f_x1_init_sp1
    sll_real64,dimension(:,:),pointer :: f_x1_sp2,f_x2_sp2,f_x1_init_sp2
    sll_real64,dimension(:),pointer :: rho_sp1,rho_sp2,efield,e_app,rho_loc_sp1,rho_loc_sp2
    !sll_real64 :: qel(1),qel_loc(1)
    !sll_real64, dimension(:), allocatable :: rho_split
    !sll_real64, dimension(:), allocatable :: rho_full
    
    sll_int32 :: rhotote_id,rhototi_id
    sll_int32 :: efield_id     
    sll_int32 :: th_diag_id
    sll_int32 :: restart_id
    type(layout_2D), pointer :: layout_x1_sp1
    type(layout_2D), pointer :: layout_x1_sp2
    type(layout_2D), pointer :: layout_x2_sp1
    type(layout_2D), pointer :: layout_x2_sp2
    type(remap_plan_2D_real64), pointer :: remap_plan_x1_x2_sp1
    type(remap_plan_2D_real64), pointer :: remap_plan_x2_x1_sp1
    type(remap_plan_2D_real64), pointer :: remap_plan_x1_x2_sp2
    type(remap_plan_2D_real64), pointer :: remap_plan_x2_x1_sp2
    sll_real64, dimension(:), pointer     :: f1d_sp1
    sll_real64, dimension(:,:), pointer     :: f1d_omp_sp1
    sll_real64, dimension(:,:), pointer     :: f1d_omp_in_sp1
    sll_real64, dimension(:,:), pointer     :: f1d_omp_out_sp1
    sll_real64, dimension(:), pointer     :: f1d_sp2
    sll_real64, dimension(:,:), pointer     :: f1d_omp_sp2
    sll_real64, dimension(:,:), pointer     :: f1d_omp_in_sp2
    sll_real64, dimension(:,:), pointer     :: f1d_omp_out_sp2
    sll_int32 :: np_x1,np_x2_sp1,np_x2_sp2
    sll_int32 :: nproc_x1,nproc_x2
    sll_int32 :: global_indices_sp1(2)
    sll_int32 :: global_indices_sp2(2)
    sll_int32 :: ierr
    sll_int32 :: local_size_x1_sp1,local_size_x1_sp2,local_size_x2_sp1,local_size_x2_sp2
    !sll_real64 :: adr
    sll_real64 :: tmp_loc(10),tmp(10)
    sll_int32 :: i,istep,ig,k
    
    sll_real64  ::   time
    sll_real64  ::   mass_sp1, momentum_sp1, l1norm_sp1, l2norm_sp1
    sll_real64  ::   mass_sp2, momentum_sp2, l1norm_sp2, l2norm_sp2
    sll_real64  ::   kinetic_energy_sp1,kinetic_energy_sp2,potential_energy

    !sll_real64, dimension(:), allocatable :: x2_array
    sll_real64, dimension(:), allocatable :: x2_array_unit_sp1
    sll_real64, dimension(:), allocatable :: x2_array_unit_sp2
    sll_real64, dimension(:), allocatable :: x2_array_middle_sp1
    sll_real64, dimension(:), allocatable :: x2_array_middle_sp2
    !sll_real64, dimension(:), allocatable :: x1_array
    sll_real64, dimension(:), allocatable :: node_positions_x2_sp1
    sll_real64, dimension(:), allocatable :: node_positions_x2_sp2
    !character(len=4)           :: fin   
    sll_int32                  :: file_id
    
    type(sll_fft_plan), pointer         :: pfwd
    sll_real64, dimension(:), allocatable :: buf_fft
    sll_comp64,dimension(:),allocatable :: rho_mode
    
    sll_int32 :: nb_mode 
    sll_real64 :: t_step
    sll_real64 :: time_init
    sll_int32 :: split_istep
    !sll_int32 :: split_x
    !sll_int32 :: split_x_init
    sll_int32 :: num_dof_x2_sp1
    sll_int32 :: num_dof_x2_sp2
     
    logical :: split_T
    !sll_int32 ::conservative_case
    
    
    ! for parallelization (output of distribution function in one single file)
    sll_int32, dimension(:), allocatable :: collective_displs_sp1
    sll_int32, dimension(:), allocatable :: collective_displs_sp2
    sll_int32, dimension(:), allocatable :: collective_recvcnts_sp1
    sll_int32, dimension(:), allocatable :: collective_recvcnts_sp2
    sll_int32 :: collective_size
    sll_real64,dimension(:,:),pointer :: f_visu_sp1 
    sll_real64,dimension(:,:),pointer :: f_visu_sp2
    sll_real64,dimension(:),pointer :: f_visu_buf1d_sp1
    sll_real64,dimension(:),pointer :: f_visu_buf1d_sp2
    sll_real64,dimension(:),pointer :: f_x1_buf1d_sp1
    sll_real64,dimension(:),pointer :: f_x1_buf1d_sp2
    sll_real64,dimension(:),pointer :: f_hat_x2_sp1_loc
    sll_real64,dimension(:),pointer :: f_hat_x2_sp1
    sll_real64,dimension(:),pointer :: f_hat_x2_sp2_loc
    sll_real64,dimension(:),pointer :: f_hat_x2_sp2

    sll_int32 :: iplot
    character(len=4) :: cproc
    character(len=4) :: cplot
    sll_int32 :: iproc
    logical :: file_exists
    sll_int32 :: tid
    sll_int32 :: i_omp
    sll_int32 :: ig_omp
    sll_real64 :: alpha_omp
    sll_real64 :: mean_omp
    !for temporary poisson
    !sll_int32 :: N_buf_poisson
    !sll_real64, dimension(:), allocatable :: buf_poisson
    
    

    iplot = 1

    nb_mode = sim%nb_mode
    time_init = sim%time_init
    np_x1 = sim%mesh2d_sp1%num_cells1+1
    np_x2_sp1 = sim%mesh2d_sp1%num_cells2+1
    num_dof_x2_sp1 = sim%num_dof_x2_sp1

    np_x2_sp2 = sim%mesh2d_sp2%num_cells2+1
    num_dof_x2_sp2 = sim%num_dof_x2_sp2

    !qel(1) = 0._f64
    
    if(sll_get_collective_rank(sll_world_collective)==0)then
       print *,'#collective_size=',sll_get_collective_size(sll_world_collective)
       SLL_ALLOCATE(f_visu_sp1(np_x1,num_dof_x2_sp1),ierr)
       SLL_ALLOCATE(f_visu_sp2(np_x1,num_dof_x2_sp2),ierr)
       SLL_ALLOCATE(f_visu_buf1d_sp1(np_x1*num_dof_x2_sp1),ierr)
       SLL_ALLOCATE(f_visu_buf1d_sp2(np_x1*num_dof_x2_sp2),ierr)
    else
       SLL_ALLOCATE(f_visu_sp1(1:1,1:1),ierr)
       SLL_ALLOCATE(f_visu_sp2(1:1,1:1),ierr)
       SLL_ALLOCATE(f_visu_buf1d_sp1(1:1),ierr)
       SLL_ALLOCATE(f_visu_buf1d_sp2(1:1),ierr)
    endif

 
    collective_size = sll_get_collective_size(sll_world_collective)
    SLL_ALLOCATE(collective_displs_sp1(collective_size),ierr)
    SLL_ALLOCATE(collective_displs_sp2(collective_size),ierr)
    SLL_ALLOCATE(collective_recvcnts_sp1(collective_size),ierr)
    SLL_ALLOCATE(collective_recvcnts_sp2(collective_size),ierr)

        
    !if(sll_get_collective_rank(sll_world_collective)==0)then
      SLL_ALLOCATE(buf_fft(np_x1-1),ierr)
      pfwd => fft_new_plan(np_x1-1,buf_fft,buf_fft,FFT_FORWARD,FFT_NORMALIZE)
      SLL_ALLOCATE(rho_mode(0:nb_mode),ierr)      
    !endif
    ! allocate and initialize the layouts...

      layout_x1_sp1       => new_layout_2D( sll_world_collective )
      layout_x1_sp2       => new_layout_2D( sll_world_collective )
      layout_x2_sp1       => new_layout_2D( sll_world_collective )
      layout_x2_sp2       => new_layout_2D( sll_world_collective )    
      nproc_x1 = sll_get_collective_size( sll_world_collective )
      nproc_x2 = 1
      call initialize_layout_with_distributed_array( &
           np_x1, num_dof_x2_sp1, nproc_x1, nproc_x2, layout_x2_sp1 )
      call initialize_layout_with_distributed_array( &
           np_x1, num_dof_x2_sp2, nproc_x1, nproc_x2, layout_x2_sp2 )
      call initialize_layout_with_distributed_array( &
           np_x1, num_dof_x2_sp1, nproc_x2, nproc_x1, layout_x1_sp1 )
      call initialize_layout_with_distributed_array( &
           np_x1, num_dof_x2_sp2, nproc_x2, nproc_x1, layout_x1_sp2 )
      !call sll_view_lims_2D( layout_x1 )

    
      !allocation of distribution functions f_x1 and f_x2
      call compute_local_sizes( layout_x2_sp1, local_size_x1_sp1, local_size_x2_sp1 )
      call compute_local_sizes( layout_x2_sp2, local_size_x1_sp2, local_size_x2_sp2 )
      
      SLL_ALLOCATE(f_x2_sp1(local_size_x1_sp1,local_size_x2_sp1),ierr)
      SLL_ALLOCATE(f_x2_sp2(local_size_x1_sp2,local_size_x2_sp2),ierr)

      call compute_local_sizes( layout_x1_sp1, local_size_x1_sp1, local_size_x2_sp1 )
      call compute_local_sizes( layout_x1_sp2, local_size_x1_sp2, local_size_x2_sp2 )

      global_indices_sp1(1:2) = local_to_global( layout_x1_sp1, (/1, 1/) )
      SLL_ALLOCATE(f_x1_sp1(local_size_x1_sp1,local_size_x2_sp1),ierr)    
      SLL_ALLOCATE(f_x1_init_sp1(local_size_x1_sp1,local_size_x2_sp1),ierr)    
      SLL_ALLOCATE(f_x1_buf1d_sp1(local_size_x1_sp1*local_size_x2_sp1),ierr)    

      global_indices_sp2(1:2) = local_to_global( layout_x1_sp2, (/1, 1/) )
      SLL_ALLOCATE(f_x1_sp2(local_size_x1_sp2,local_size_x2_sp2),ierr)    
      SLL_ALLOCATE(f_x1_init_sp2(local_size_x1_sp2,local_size_x2_sp2),ierr)    
      SLL_ALLOCATE(f_x1_buf1d_sp2(local_size_x1_sp2*local_size_x2_sp2),ierr)    



    !definition of remap
    remap_plan_x1_x2_sp1 => NEW_REMAP_PLAN(layout_x1_sp1, layout_x2_sp1, f_x1_sp1)
    remap_plan_x2_x1_sp1 => NEW_REMAP_PLAN(layout_x2_sp1, layout_x1_sp1, f_x2_sp1)
    remap_plan_x1_x2_sp2 => NEW_REMAP_PLAN(layout_x1_sp2, layout_x2_sp2, f_x1_sp2)
    remap_plan_x2_x1_sp2 => NEW_REMAP_PLAN(layout_x2_sp2, layout_x1_sp2, f_x2_sp2)
    
    !      print *,'hello',size(f_x1,1),size(f_x1,2),size(f_x2,1),size(f_x2,2)
    !      call apply_remap_2D( remap_plan_x2_x1, f_x2, f_x1 )
    !      print *,'hello2',size(f_x1,1),size(f_x1,2),size(f_x2,1),size(f_x2,2)


    
    !allocation of 1d arrays
    SLL_ALLOCATE(rho_sp1(np_x1),ierr)
    SLL_ALLOCATE(rho_sp2(np_x1),ierr)
    SLL_ALLOCATE(rho_loc_sp1(np_x1),ierr)
    SLL_ALLOCATE(rho_loc_sp2(np_x1),ierr)
    SLL_ALLOCATE(efield(np_x1),ierr)
    SLL_ALLOCATE(e_app(np_x1),ierr)
    e_app = 0._f64
    
    SLL_ALLOCATE(f1d_sp1(max(np_x1,np_x2_sp1)),ierr)
    SLL_ALLOCATE(f1d_omp_sp1(max(np_x1,np_x2_sp1),sim%num_threads),ierr)
    SLL_ALLOCATE(f1d_omp_in_sp1(max(np_x1,np_x2_sp1),sim%num_threads),ierr)
    SLL_ALLOCATE(f1d_omp_out_sp1(max(np_x1,np_x2_sp1),sim%num_threads),ierr)
    SLL_ALLOCATE(x2_array_unit_sp1(np_x2_sp1),ierr)
    SLL_ALLOCATE(x2_array_middle_sp1(np_x2_sp1),ierr)
    SLL_ALLOCATE(node_positions_x2_sp1(num_dof_x2_sp1),ierr)
    SLL_ALLOCATE(f_hat_x2_sp1_loc(nb_mode+1),ierr)
    SLL_ALLOCATE(f_hat_x2_sp1(nb_mode+1),ierr)

    SLL_ALLOCATE(f1d_sp2(max(np_x1,np_x2_sp2)),ierr)
    SLL_ALLOCATE(f1d_omp_sp2(max(np_x1,np_x2_sp2),sim%num_threads),ierr)
    SLL_ALLOCATE(f1d_omp_in_sp2(max(np_x1,np_x2_sp2),sim%num_threads),ierr)
    SLL_ALLOCATE(f1d_omp_out_sp2(max(np_x1,np_x2_sp2),sim%num_threads),ierr)
    SLL_ALLOCATE(x2_array_unit_sp2(np_x2_sp2),ierr)
    SLL_ALLOCATE(x2_array_middle_sp2(np_x2_sp2),ierr)
    SLL_ALLOCATE(node_positions_x2_sp2(num_dof_x2_sp2),ierr)
    SLL_ALLOCATE(f_hat_x2_sp2_loc(nb_mode+1),ierr)
    SLL_ALLOCATE(f_hat_x2_sp2(nb_mode+1),ierr)


    !temporary poisson
    !N_buf_poisson=2*(np_x1-1)+15  
    !allocate(buf_poisson(N_buf_poisson))
    !call poisson1dper_init(buf_poisson,np_x1-1)




    x2_array_unit_sp1(1:np_x2_sp1) = &
      (sim%x2_array_sp1(1:np_x2_sp1)-sim%x2_array_sp1(1))/(sim%x2_array_sp1(np_x2_sp1)-sim%x2_array_sp1(1))
    do i = 1, np_x2_sp1-1
       x2_array_middle_sp1(i) = 0.5_f64*(sim%x2_array_sp1(i)+sim%x2_array_sp1(i+1))
    end do
    x2_array_middle_sp1(np_x2_sp1) = x2_array_middle_sp1(1)+sim%x2_array_sp1(np_x2_sp1)-sim%x2_array_sp1(1)
    
    select case (sim%advection_form_x2_sp1)
      case (SLL_ADVECTIVE)
        node_positions_x2_sp1(1:num_dof_x2_sp1) = sim%x2_array_sp1(1:num_dof_x2_sp1)
      case (SLL_CONSERVATIVE)
        node_positions_x2_sp1(1:num_dof_x2_sp1) = x2_array_middle_sp1(1:num_dof_x2_sp1)
      case default
        print *,'#sim%advection_form_x2_sp1=',sim%advection_form_x2_sp1
        print *,'#not implemented'
        print *,'#in run_vp2d_cartesian_two_species'
        stop
    end select  

    x2_array_unit_sp2(1:np_x2_sp2) = &
         (sim%x2_array_sp2(1:np_x2_sp2)-sim%x2_array_sp2(1))/(sim%x2_array_sp2(np_x2_sp2)-sim%x2_array_sp2(1))
    do i = 1, np_x2_sp2-1
       x2_array_middle_sp2(i) = 0.5_f64*(sim%x2_array_sp2(i)+sim%x2_array_sp2(i+1))
    end do
    x2_array_middle_sp2(np_x2_sp2) = x2_array_middle_sp2(1)+sim%x2_array_sp2(np_x2_sp2)-sim%x2_array_sp2(1)

    select case (sim%advection_form_x2_sp2)
      case (SLL_ADVECTIVE)
        node_positions_x2_sp2(1:num_dof_x2_sp2) = sim%x2_array_sp2(1:num_dof_x2_sp2)
      case (SLL_CONSERVATIVE)
        node_positions_x2_sp2(1:num_dof_x2_sp2) = x2_array_middle_sp2(1:num_dof_x2_sp2)
      case default
        print *,'#sim%advection_form_x2_sp2=',sim%advection_form_x2_sp2
        print *,'#not implemented'
        print *,'#in run_vp2d_cartesian_two_species'
        stop
    end select  
        

    
    !initialize distribution function     
    call sll_2d_parallel_array_initializer_cartesian( &
       layout_x1_sp1, &
       sim%x1_array, &
       node_positions_x2_sp1, &
       f_x1_sp1, &
       sim%init_func_sp1, &
       sim%params_sp1)


    call sll_2d_parallel_array_initializer_cartesian( &
       layout_x1_sp2, &
       sim%x1_array, &
       node_positions_x2_sp2, &
       f_x1_sp2, &
       sim%init_func_sp2, &
       sim%params_sp2)



    iproc = sll_get_collective_rank(sll_world_collective)
    call int2string(iproc, cproc)
    call int2string(iplot,cplot)    

    if(sim%restart_file/="no_restart_file")then
      INQUIRE(FILE=trim(sim%restart_file)//'_proc_'//cproc//'.rst', EXIST=file_exists)
      if(.not.(file_exists))then
        print *,'#file ',trim(sim%restart_file)//'_proc_'//cproc//'.rst'
        print *,'does not exist'
        stop
      endif
      open(unit=restart_id, &
        file=trim(sim%restart_file)//'_proc_'//cproc//'.rst', ACCESS="STREAM", &
        form='unformatted', IOStat=ierr)      
      if( ierr .ne. 0 ) then
        print *, 'ERROR while opening file ', &
          trim(sim%restart_file)//'_proc_'//cproc//'.rst', &
           '. Called from run_vp2d_cartesian().'
       stop
      end if
      print *,'#read restart file '//trim(sim%restart_file)//'_proc_'//cproc//'.rst'      
      call sll_binary_read_array_0d(restart_id,time_init,ierr)
      call sll_binary_read_array_2d(restart_id,f_x1_sp1(1:local_size_x1_sp1,1:local_size_x2_sp1),ierr)
      call sll_binary_file_close(restart_id,ierr)
    endif      
    
    if(sim%time_init_from_restart_file .eqv. .true.) then
      sim%time_init = time_init  
    endif
    time_init = sim%time_init

    !call sll_binary_file_create('f_plot_'//cplot//'_proc_'//cproc//'.rst', restart_id, ierr )
    !call sll_binary_write_array_0d(restart_id,time_init,ierr)
    !call sll_binary_write_array_2d(restart_id,f_x1(1:local_size_x1,1:local_size_x2),ierr)
    !call sll_binary_file_close(restart_id,ierr)    


    
    call sll_2d_parallel_array_initializer_cartesian( &
       layout_x1_sp1, &
       sim%x1_array, &
       node_positions_x2_sp1, &
       f_x1_init_sp1, &
       sll_landau_initializer_2d, &
       (/ sim%params_sp1(1),0._f64,sim%params_sp1(3),sim%params_sp1(4) /))
       !(/ sim%kx_sp1,0._f64 /))

    call sll_2d_parallel_array_initializer_cartesian( &
       layout_x1_sp2, &
       sim%x1_array, &
       node_positions_x2_sp2, &
       f_x1_init_sp2, &
       sll_landau_initializer_2d, &
       (/ sim%params_sp2(1),0._f64,sim%params_sp2(3),sim%params_sp2(4) /))
       !(/ sim%kx_sp2,0._f64 /))
    
    
    call compute_displacements_array_2d( &
      layout_x1_sp1, &
      collective_size, &
      collective_displs_sp1 )
    collective_recvcnts_sp1 = receive_counts_array_2d( &
      layout_x1_sp1, &
      collective_size )

    call compute_displacements_array_2d( &
      layout_x1_sp2, &
      collective_size, &
      collective_displs_sp2 )
    collective_recvcnts_sp2 = receive_counts_array_2d( &
      layout_x1_sp2, &
      collective_size )

    call load_buffer_2d( layout_x1_sp1, f_x1_sp1, f_x1_buf1d_sp1 )

    call load_buffer_2d( layout_x1_sp2, f_x1_sp2, f_x1_buf1d_sp2 )

    call sll_collective_gatherv_real64( &
      sll_world_collective, &
      f_x1_buf1d_sp1, &
      local_size_x1_sp1*local_size_x2_sp1, &
      collective_recvcnts_sp1, &
      collective_displs_sp1, &
      0, &
      f_visu_buf1d_sp1 )


    call sll_collective_gatherv_real64( &
      sll_world_collective, &
      f_x1_buf1d_sp2, &
      local_size_x1_sp2*local_size_x2_sp2, &
      collective_recvcnts_sp2, &
      collective_displs_sp2, &
      0, &
      f_visu_buf1d_sp2 )

    


    f_visu_sp1 = reshape(f_visu_buf1d_sp1, shape(f_visu_sp1))
    if(sll_get_collective_rank(sll_world_collective)==0)then
      call sll_binary_file_create('f0.bdat', file_id, ierr)
      call sll_binary_write_array_2d(file_id,f_visu_sp1(1:np_x1-1,1:np_x2_sp1-1),ierr)
      call sll_binary_file_close(file_id,ierr)
#ifndef NOHDF5
      call plot_f_cartesian( &
        iplot, &
        f_visu_sp1, &
        sim%x1_array, &
        np_x1, &
        node_positions_x2_sp1, &
        sim%num_dof_x2_sp1, &
        'fe', 'e', &
        time_init )        

#endif



      print *,'#maxfe = maxf_sp1',maxval(f_visu_sp1), minval(f_visu_sp1) 

    endif


    f_visu_sp2 = reshape(f_visu_buf1d_sp2, shape(f_visu_sp2))
    if(sll_get_collective_rank(sll_world_collective)==0)then
      call sll_binary_file_create('f0.bdat', file_id, ierr)
      call sll_binary_write_array_2d(file_id,f_visu_sp2(1:np_x1-1,1:np_x2_sp2-1),ierr)
      call sll_binary_file_close(file_id,ierr)
#ifndef NOHDF5
      call plot_f_cartesian( &
        iplot, &
        f_visu_sp2, &
        sim%x1_array, &
        np_x1, &
        node_positions_x2_sp2, &
        sim%num_dof_x2_sp2, &
        'fi', 'i', &
        time_init )
#endif



      print *,'#maxfi = maxf_sp2',maxval(f_visu_sp2), minval(f_visu_sp2) 

    endif
    
    



    
    


    !computation of electric field
    rho_loc_sp1 = 0._f64
    rho_loc_sp2 = 0._f64
    ig = global_indices_sp1(2)-1
    do i=1,np_x1
      rho_loc_sp1(i)=rho_loc_sp1(i)&
        +sum(f_x1_sp1(i,1:local_size_x2_sp1)*sim%integration_weight_sp1(1+ig:local_size_x2_sp1+ig))
   end do
   ig = global_indices_sp2(2)-1
   do i = 1,np_x1
      rho_loc_sp2(i)=rho_loc_sp2(i)&
        +sum(f_x1_sp2(i,1:local_size_x2_sp2)*sim%integration_weight_sp2(1+ig:local_size_x2_sp2+ig))
   end do
   

    call sll_collective_allreduce( &
      sll_world_collective, &
      rho_loc_sp1, &
      np_x1, &
      MPI_SUM, &
      rho_sp1 )

    call sll_collective_allreduce( &
      sll_world_collective, &
      rho_loc_sp2, &
      np_x1, &
      MPI_SUM, &
      rho_sp2 )


    !sim%mixt_bc = (/ 0._f64, 0._f64 /)
    
    call sim%poisson%compute_E_from_rho( efield, rho_sp2-rho_sp1 )
    
    istep = 0
    if(sim%driven)then
      call compute_e_app(sim,e_app,time_init+real(istep,f64)*sim%dt)
    endif
    
    ! write initial fields
    if(sll_get_collective_rank(sll_world_collective)==0)then
      call sll_ascii_file_create('thdiag.dat', th_diag_id, ierr)
      call sll_ascii_file_create('rhotote.dat', rhotote_id, ierr)
      call sll_ascii_file_create('rhototi.dat', rhototi_id, ierr)
      call sll_ascii_file_create('efield.dat', efield_id, ierr)

      write(rhototi_id,*) sim%x1_array
      write(rhotote_id,*) sim%x1_array
      write(efield_id,*) sim%x1_array
    endif
    
    
    !write also initial deltaf function
    call load_buffer_2d( layout_x1_sp1, f_x1_sp1-f_x1_init_sp1, f_x1_buf1d_sp1 )
    call sll_collective_gatherv_real64( &
      sll_world_collective, &
      f_x1_buf1d_sp1, &
      local_size_x1_sp1*local_size_x2_sp1, &
      collective_recvcnts_sp1, &
      collective_displs_sp1, &
      0, &
      f_visu_buf1d_sp1 )
    f_visu_sp1 = reshape(f_visu_buf1d_sp1, shape(f_visu_sp1))


    call load_buffer_2d( layout_x1_sp2, f_x1_sp2-f_x1_init_sp2, f_x1_buf1d_sp2 )
    call sll_collective_gatherv_real64( &
      sll_world_collective, &
      f_x1_buf1d_sp2, &
      local_size_x1_sp2*local_size_x2_sp2, &
      collective_recvcnts_sp2, &
      collective_displs_sp2, &
      0, &
      f_visu_buf1d_sp2 )
    f_visu_sp2 = reshape(f_visu_buf1d_sp2, shape(f_visu_sp2))




    if(sll_get_collective_rank(sll_world_collective)==0) then        
      print *,'#step=',0,time_init+real(0,f64)*sim%dt,'iplot=',iplot
    endif
    iplot = iplot+1  
    !sim%num_iterations = 0
    do istep = 1, sim%num_iterations
      if (mod(istep,sim%freq_diag)==0) then
        if(sll_get_collective_rank(sll_world_collective)==0) then        
          print *,'#step=',istep,time_init+real(istep,f64)*sim%dt,'iplot=',iplot
        endif
      endif  

      split_T = sim%split%split_begin_T
      t_step = real(istep-1,f64)
      do split_istep=1,sim%split%nb_split_step
        if(split_T) then
          !! T ADVECTION 
          tid=1          
!#ifdef _OPENMP
!!$OMP PARALLEL DEFAULT(SHARED) &
!!$OMP PRIVATE(i_omp,ig_omp,alpha_omp,tid) 
!          !advection in x
!          tid = omp_get_thread_num()+1
!!$OMP DO
!#endif

          do i_omp = 1, local_size_x2_sp1
             ig_omp = i_omp+global_indices_sp1(2)-1
             alpha_omp = sim%factor_x1*node_positions_x2_sp1(ig_omp) * sim%split%split_step(split_istep) 
             f1d_omp_in_sp1(1:np_x1,tid) = f_x1_sp1(1:np_x1,i_omp)
             call sim%advect_x1_sp1(tid)%ptr%advect_1d_constant(&
                  alpha_omp, &
                  sim%dt, &
                  f1d_omp_in_sp1(1:np_x1,tid), &
                  f1d_omp_out_sp1(1:np_x1,tid))
             
             f_x1_sp1(1:np_x1,i_omp)=f1d_omp_out_sp1(1:np_x1,tid)
          end do
          
          do i_omp = 1,local_size_x2_sp2
             ig_omp = i_omp+global_indices_sp2(2)-1   
             alpha_omp = sim%factor_x1*node_positions_x2_sp2(ig_omp) * sim%split%split_step(split_istep) 
             f1d_omp_in_sp2(1:np_x1,tid) = f_x1_sp2(1:np_x1,i_omp)
             call sim%advect_x1_sp2(tid)%ptr%advect_1d_constant(&
                  alpha_omp, &
                  sim%dt, &
                  f1d_omp_in_sp2(1:np_x1,tid), &
                  f1d_omp_out_sp2(1:np_x1,tid))
             
             f_x1_sp2(1:np_x1,i_omp)=f1d_omp_out_sp2(1:np_x1,tid)
             
          end do
!#ifdef _OPENMP
!!$OMP END DO          
!!$OMP END PARALLEL
!#endif
          t_step = t_step+sim%split%split_step(split_istep)
          !computation of electric field
          rho_loc_sp1 = 0._f64
          rho_loc_sp2 = 0._f64
          !qel_loc(1) = 0._f64
          ig = global_indices_sp1(2)-1
          do i=1,np_x1
             rho_loc_sp1(i)=rho_loc_sp1(i)&
                  +sum(f_x1_sp1(i,1:local_size_x2_sp1)*sim%integration_weight_sp1(1+ig:local_size_x2_sp1+ig))
          end do

          !qel_loc(1) = qel_loc(1) + sum(f_x1_sp1(1,1:local_size_x2_sp1)*(-abs(node_positions_x2_sp1(ig+1:ig+local_size_x2_sp1))+node_positions_x2_sp1(ig+1:ig+local_size_x2_sp1))*0.5_f64*sim%integration_weight_sp1(1+ig:local_size_x2_sp1+ig))

          ig = global_indices_sp2(2)-1
          do i = 1,np_x1
             rho_loc_sp2(i)=rho_loc_sp2(i)&
                  +sum(f_x1_sp2(i,1:local_size_x2_sp2)*sim%integration_weight_sp2(1+ig:local_size_x2_sp2+ig))
          end do

          !qel_loc(1) = qel_loc(1) - sum(f_x1_sp2(1,1:local_size_x2_sp2)*(-abs(node_positions_x2_sp2(ig+1:ig+local_size_x2_sp2))+node_positions_x2_sp2(ig+1:ig+local_size_x2_sp2))*0.5_f64*sim%integration_weight_sp2(1+ig:local_size_x2_sp2+ig))
          
          call sll_collective_allreduce( &
            sll_world_collective, &
            rho_loc_sp1, &
            np_x1, &
            MPI_SUM, &
            rho_sp1 )

          call sll_collective_allreduce( &
            sll_world_collective, &
            rho_loc_sp2, &
            np_x1, &
            MPI_SUM, &
            rho_sp2 )

!          call sll_collective_allreduce( &
!               sll_world_collective, &
!               qel_loc, &
!               1, &
!               MPI_SUM, &
!               qel )

          !sim%mixt_bc(1) = sim%mixt_bc(1) - 0.5_f64*qel(1)*sim%dt

          call sim%poisson%compute_E_from_rho( efield, rho_sp2-rho_sp1 )
          
          t_step = t_step+sim%split%split_step(split_istep)
          
          if(sim%driven)then
            call compute_e_app(sim,e_app,time_init+t_step*sim%dt)
          endif

          
          
       else

          !! V ADVECTION 
          !transposition
           call apply_remap_2D( remap_plan_x1_x2_sp1, f_x1_sp1, f_x2_sp1 )
           call apply_remap_2D( remap_plan_x1_x2_sp2, f_x1_sp2, f_x2_sp2 )
           call compute_local_sizes( layout_x2_sp1, local_size_x1_sp1, local_size_x2_sp1 )
           call compute_local_sizes( layout_x2_sp2, local_size_x1_sp2, local_size_x2_sp2 )
           global_indices_sp1(1:2) = local_to_global( layout_x2_sp1, (/1, 1/) )
           global_indices_sp2(1:2) = local_to_global( layout_x2_sp2, (/1, 1/) )
           tid = 1
           
!#ifdef _OPENMP
!!$OMP PARALLEL DEFAULT(SHARED) &
!!$OMP PRIVATE(i_omp,ig_omp,alpha_omp,tid,mean_omp,f1d) 
!          !advection in v
!          tid = omp_get_thread_num()+1
!!$OMP DO
!#endif
          !advection in v
          do i_omp = 1,local_size_x1_sp1
             ig_omp=i_omp+global_indices_sp1(1)-1
             alpha_omp = -(efield(ig_omp)+sim%alpha_sp1*e_app(ig_omp)) * sim%split%split_step(split_istep)
             f1d_omp_in_sp1(1:num_dof_x2_sp1,tid) = f_x2_sp1(i_omp,1:num_dof_x2_sp1)
             if(sim%advection_form_x2_sp1==SLL_CONSERVATIVE)then
                call function_to_primitive(f1d_omp_in_sp1(:,tid),x2_array_unit_sp1,np_x2_sp1-1,mean_omp)
             endif
             call sim%advect_x2_sp1(tid)%ptr%advect_1d_constant(&
                  alpha_omp, &
                  sim%dt, &
                  f1d_omp_in_sp1(1:num_dof_x2_sp1,tid), &
                  f1d_omp_out_sp1(1:num_dof_x2_sp1,tid))
             if(sim%advection_form_x2_sp1==SLL_CONSERVATIVE)then
                call primitive_to_function(f1d_omp_out_sp1(:,tid),x2_array_unit_sp1,np_x2_sp1-1,mean_omp)
             endif
             f_x2_sp1(i_omp,1:num_dof_x2_sp1) = f1d_omp_out_sp1(1:num_dof_x2_sp1,tid)
          end do

          do i_omp = 1,local_size_x1_sp2
             ig_omp=i_omp+global_indices_sp2(1)-1
             alpha_omp = sim%mass_ratio*(efield(ig_omp)-sim%alpha_sp2*e_app(ig_omp)) * sim%split%split_step(split_istep)
             f1d_omp_in_sp2(1:num_dof_x2_sp2,tid) = f_x2_sp2(i_omp,1:num_dof_x2_sp2)
             if(sim%advection_form_x2_sp2==SLL_CONSERVATIVE)then
                call function_to_primitive(f1d_omp_in_sp2(:,tid),x2_array_unit_sp2,np_x2_sp2-1,mean_omp)
             endif
             
             call sim%advect_x2_sp2(tid)%ptr%advect_1d_constant(&
                  alpha_omp, &
                  sim%dt, &
                  f1d_omp_in_sp2(1:num_dof_x2_sp2,tid), &
                  f1d_omp_out_sp2(1:num_dof_x2_sp2,tid))
             
             if(sim%advection_form_x2_sp2==SLL_CONSERVATIVE)then
                call primitive_to_function(f1d_omp_out_sp2(:,tid),x2_array_unit_sp2,np_x2_sp2-1,mean_omp)
             endif
             f_x2_sp2(i_omp,1:num_dof_x2_sp2) = f1d_omp_out_sp2(1:num_dof_x2_sp2,tid)
             
          end do
!#ifdef _OPENMP
!!$OMP END DO          
!!$OMP END PARALLEL
!#endif
         !transposition
         call apply_remap_2D( remap_plan_x2_x1_sp1, f_x2_sp1, f_x1_sp1 )
         call apply_remap_2D( remap_plan_x2_x1_sp2, f_x2_sp2, f_x1_sp2 )
         call compute_local_sizes( layout_x1_sp1, local_size_x1_sp1, local_size_x2_sp1 )
         call compute_local_sizes( layout_x1_sp2, local_size_x1_sp2, local_size_x2_sp2 )
         global_indices_sp1(1:2) = local_to_global( layout_x1_sp1, (/1, 1/) )
         global_indices_sp2(1:2) = local_to_global( layout_x1_sp2, (/1, 1/) )


       endif
       !split_x= 1-split_x
       split_T = .not.(split_T)
      enddo
      

     !!DIAGNOSTICS
     if (mod(istep,sim%freq_diag_time)==0) then
        time = time_init+real(istep,f64)*sim%dt
        mass_sp1 = 0._f64
        momentum_sp1 = 0._f64
        l1norm_sp1 = 0._f64
        l2norm_sp1 = 0._f64
        kinetic_energy_sp1 = 0._f64

        mass_sp2 = 0._f64
        momentum_sp2 = 0._f64
        l1norm_sp2 = 0._f64
        l2norm_sp2 = 0._f64
        kinetic_energy_sp2 = 0._f64

        potential_energy = 0._f64
        tmp_loc = 0._f64
        ig = global_indices_sp1(2)-1               
        do i = 1, np_x1-1        
           tmp_loc(1)= tmp_loc(1)+sum(f_x1_sp1(i,1:local_size_x2_sp1) &
                *sim%integration_weight_sp1(1+ig:local_size_x2_sp1+ig))
           tmp_loc(2)= tmp_loc(2)+sum(abs(f_x1_sp1(i,1:local_size_x2_sp1)) &
                *sim%integration_weight_sp1(1+ig:local_size_x2_sp1+ig))
           tmp_loc(3)= tmp_loc(3)+sum((f_x1_sp1(i,1:local_size_x2_sp1))**2 &
                *sim%integration_weight_sp1(1+ig:local_size_x2_sp1+ig))
           tmp_loc(4)= tmp_loc(4) +sum(f_x1_sp1(i,1:local_size_x2_sp1) &
                *sim%x2_array_sp1(global_indices_sp1(2)-1+1:global_indices_sp1(2)-1+local_size_x2_sp1) &
                *sim%integration_weight_sp1(1+ig:local_size_x2_sp1+ig))          
           tmp_loc(5)= tmp_loc(5)+sum(f_x1_sp1(i,1:local_size_x2_sp1) &
                *sim%x2_array_sp1(global_indices_sp1(2)-1+1:global_indices_sp1(2)-1+local_size_x2_sp1)**2 &
                *sim%integration_weight_sp1(1+ig:local_size_x2_sp1+ig) )          
        end do
        
        ig = global_indices_sp2(2)-1
        do i = 1, np_x1-1        
           tmp_loc(6)= tmp_loc(6)+sum(f_x1_sp2(i,1:local_size_x2_sp2) &
                *sim%integration_weight_sp2(1+ig:local_size_x2_sp2+ig))
           tmp_loc(7)= tmp_loc(7)+sum(abs(f_x1_sp2(i,1:local_size_x2_sp2)) &
                *sim%integration_weight_sp2(1+ig:local_size_x2_sp2+ig))
           tmp_loc(8)= tmp_loc(8)+sum((f_x1_sp2(i,1:local_size_x2_sp2))**2 &
                *sim%integration_weight_sp2(1+ig:local_size_x2_sp2+ig))
           tmp_loc(9)= tmp_loc(9) +sum(f_x1_sp2(i,1:local_size_x2_sp2) &
                *sim%x2_array_sp2(global_indices_sp2(2)-1+1:global_indices_sp2(2)-1+local_size_x2_sp2) &
                *sim%integration_weight_sp2(1+ig:local_size_x2_sp2+ig))          
           tmp_loc(10)= tmp_loc(10)+sum(f_x1_sp2(i,1:local_size_x2_sp2) &
                *sim%x2_array_sp2(global_indices_sp2(2)-1+1:global_indices_sp2(2)-1+local_size_x2_sp2)**2 &
                *sim%integration_weight_sp2(1+ig:local_size_x2_sp2+ig) )          
        end do

        call sll_collective_allreduce( &
          sll_world_collective, &
          tmp_loc, &
          10, &
          MPI_SUM, &
          tmp )


        mass_sp1 = tmp(1) * sim%mesh2d_sp1%delta_eta1 !* sim%mesh2d%delta_eta2
        l1norm_sp1 = tmp(2)  * sim%mesh2d_sp1%delta_eta1 !* sim%mesh2d%delta_eta2
        l2norm_sp1 = tmp(3)  * sim%mesh2d_sp1%delta_eta1 !* sim%mesh2d%delta_eta2
        momentum_sp1 = tmp(4) * sim%mesh2d_sp1%delta_eta1 !* sim%mesh2d%delta_eta2
        kinetic_energy_sp1 = 0.5_f64 *tmp(5) * sim%mesh2d_sp1%delta_eta1 !* sim%mesh2d%delta_eta2

        mass_sp2 = tmp(6) * sim%mesh2d_sp2%delta_eta1 !* sim%mesh2d%delta_eta2
        l1norm_sp2 = tmp(7)  * sim%mesh2d_sp2%delta_eta1 !* sim%mesh2d%delta_eta2
        l2norm_sp2 = tmp(8)  * sim%mesh2d_sp2%delta_eta1 !* sim%mesh2d%delta_eta2
        momentum_sp2 = tmp(9) * sim%mesh2d_sp2%delta_eta1 !* sim%mesh2d%delta_eta2
        kinetic_energy_sp2 = 0.5_f64 *tmp(10) * sim%mesh2d_sp2%delta_eta1 !* sim%mesh2d%delta_eta2

        potential_energy = 0._f64
        do i=1, np_x1-1
          potential_energy = potential_energy+(efield(i)+e_app(i))**2
        enddo
        potential_energy = 0.5_f64*potential_energy* sim%mesh2d_sp1%delta_eta1
        
        f_hat_x2_sp1_loc(1:nb_mode+1) = 0._f64
        do i=1,local_size_x2_sp1
          buf_fft = f_x1_sp1(1:np_x1-1,i)
          call fft_apply_plan(pfwd,buf_fft,buf_fft)
          do k=0,nb_mode
            f_hat_x2_sp1_loc(k+1) = f_hat_x2_sp1_loc(k+1) &
              +abs(fft_get_mode(pfwd,buf_fft,k))**2 &
              *sim%integration_weight_sp1(ig+i)
          enddo
        enddo
        call sll_collective_allreduce( &
          sll_world_collective, &
          f_hat_x2_sp1_loc, &
          nb_mode+1, &
          MPI_SUM, &
          f_hat_x2_sp1 )

        f_hat_x2_sp2_loc(1:nb_mode+1) = 0._f64
        do i=1,local_size_x2_sp2
          buf_fft = f_x1_sp2(1:np_x1-1,i)
          call fft_apply_plan(pfwd,buf_fft,buf_fft)
          do k=0,nb_mode
            f_hat_x2_sp2_loc(k+1) = f_hat_x2_sp2_loc(k+1) &
              +abs(fft_get_mode(pfwd,buf_fft,k))**2 &
              *sim%integration_weight_sp2(ig+i)
          enddo
        enddo
        call sll_collective_allreduce( &
          sll_world_collective, &
          f_hat_x2_sp2_loc, &
          nb_mode+1, &
          MPI_SUM, &
          f_hat_x2_sp2 )

        
        
        if (mod(istep,sim%freq_diag_restart)==0) then          
          call int2string(iplot,cplot) 
          call sll_binary_file_create('f_plot_'//cplot//'_proc_'//cproc//'.rst', restart_id, ierr )
          call sll_binary_write_array_0d(restart_id,time,ierr)
          call sll_binary_write_array_2d(restart_id,f_x1_sp1(1:local_size_x1_sp1,1:local_size_x2_sp1),ierr)
          call sll_binary_write_array_2d(restart_id,f_x1_sp2(1:local_size_x1_sp2,1:local_size_x2_sp2),ierr)
          call sll_binary_file_close(restart_id,ierr)    
        endif 

        if(sll_get_collective_rank(sll_world_collective)==0)then                  
          buf_fft = rho_sp1(1:np_x1-1)-rho_sp2(1:np_x1-1)
          call fft_apply_plan(pfwd,buf_fft,buf_fft)
          do k=0,nb_mode
            rho_mode(k)=fft_get_mode(pfwd,buf_fft,k)
          enddo  
          write(th_diag_id,'(f12.5,12g20.12)',advance='no') &
            time, &
            mass_sp1, &
            l1norm_sp1, &
            momentum_sp1, &
            l2norm_sp1, &
            kinetic_energy_sp1, &
            mass_sp2, &
            l1norm_sp2, &
            momentum_sp2, &
            l2norm_sp2, &
            kinetic_energy_sp2, &
            potential_energy, &
            kinetic_energy_sp1 + kinetic_energy_sp2 + potential_energy
          do k=0,nb_mode
            write(th_diag_id,'(g20.12)',advance='no') &
              abs(rho_mode(k))
          enddo
          do k=0,nb_mode-1
            write(th_diag_id,'(2g20.12)',advance='no') &
                 f_hat_x2_sp1(k+1), &
                 f_hat_x2_sp2(k+1)
          enddo
          write(th_diag_id,'(2g20.12)') &
               f_hat_x2_sp1(nb_mode+1), &
               f_hat_x2_sp2(nb_mode+1)

          write(efield_id,*) efield
          write(rhotote_id,*) rho_sp1
          write(rhototi_id,*) rho_sp2
        endif
          
        if (mod(istep,sim%freq_diag)==0) then          
          !we substract f0
          !we gather in one file

          !sp1
          call load_buffer_2d( layout_x1_sp1, f_x1_sp1-f_x1_init_sp1, f_x1_buf1d_sp1 )
          call sll_collective_gatherv_real64( &
            sll_world_collective, &
            f_x1_buf1d_sp1, &
            local_size_x1_sp1*local_size_x2_sp1, &
            collective_recvcnts_sp1, &
            collective_displs_sp1, &
            0, &
            f_visu_buf1d_sp1 )
          f_visu_sp1 = reshape(f_visu_buf1d_sp1, shape(f_visu_sp1))
          if(sll_get_collective_rank(sll_world_collective)==0) then
            do i=1,num_dof_x2_sp1
              f_visu_buf1d_sp1(i) = sum(f_visu_sp1(1:np_x1-1,i))*sim%mesh2d_sp1%delta_eta1
            enddo
            call sll_gnuplot_1d( &
              f_visu_buf1d_sp1(1:num_dof_x2_sp1), &
              node_positions_x2_sp1(1:num_dof_x2_sp1), &
              'intdeltafedx', &
              iplot )

#ifndef NOHDF5
            call plot_f_cartesian( &
              iplot, &
              f_visu_sp1, &
              sim%x1_array, &
              np_x1, &
              node_positions_x2_sp1, &
              sim%num_dof_x2_sp1, &
              'deltafe', 'e', &
              time)                    
#endif
          
          endif
          !we store f for visu
          call load_buffer_2d( layout_x1_sp1, f_x1_sp1, f_x1_buf1d_sp1 )
          call sll_collective_gatherv_real64( &
            sll_world_collective, &
            f_x1_buf1d_sp1, &
            local_size_x1_sp1*local_size_x2_sp1, &
            collective_recvcnts_sp1, &
            collective_displs_sp1, &
            0, &
            f_visu_buf1d_sp1 )
          f_visu_sp1 = reshape(f_visu_buf1d_sp1, shape(f_visu_sp1))
          if(sll_get_collective_rank(sll_world_collective)==0) then
            do i=1,num_dof_x2_sp1
              f_visu_buf1d_sp1(i) = sum(f_visu_sp1(1:np_x1-1,i))*sim%mesh2d_sp1%delta_eta1
            enddo
            call sll_gnuplot_1d( &
              f_visu_buf1d_sp1(1:num_dof_x2_sp1), &
              node_positions_x2_sp1(1:num_dof_x2_sp1), &
              'intfedx', &
              iplot )
#ifndef NOHDF5
        call plot_f_cartesian( &
          iplot, &
          f_visu_sp1, &
          sim%x1_array, &
          np_x1, &
          node_positions_x2_sp1, &
          sim%num_dof_x2_sp1, &
          'fe', 'e', &
          time)                    
#endif
          endif

          
          !sp2
          call load_buffer_2d( layout_x1_sp2, f_x1_sp2-f_x1_init_sp2, f_x1_buf1d_sp2 )
          call sll_collective_gatherv_real64( &
            sll_world_collective, &
            f_x1_buf1d_sp2, &
            local_size_x1_sp2*local_size_x2_sp2, &
            collective_recvcnts_sp2, &
            collective_displs_sp2, &
            0, &
            f_visu_buf1d_sp2 )
          f_visu_sp2 = reshape(f_visu_buf1d_sp2, shape(f_visu_sp2))
          if(sll_get_collective_rank(sll_world_collective)==0) then
            do i=1,num_dof_x2_sp2
              f_visu_buf1d_sp2(i) = sum(f_visu_sp2(1:np_x1-1,i))*sim%mesh2d_sp2%delta_eta1
            enddo
            call sll_gnuplot_1d( &
              f_visu_buf1d_sp2(1:num_dof_x2_sp2), &
              node_positions_x2_sp2(1:num_dof_x2_sp2), &
              'intdeltafedx', &
              iplot )

#ifndef NOHDF5
            call plot_f_cartesian( &
              iplot, &
              f_visu_sp2, &
              sim%x1_array, &
              np_x1, &
              node_positions_x2_sp2, &
              sim%num_dof_x2_sp2, &
              'deltafi', 'i', &
              time)                    
#endif
          
          endif
          !we store f for visu
          call load_buffer_2d( layout_x1_sp2, f_x1_sp2, f_x1_buf1d_sp2 )
          call sll_collective_gatherv_real64( &
            sll_world_collective, &
            f_x1_buf1d_sp2, &
            local_size_x1_sp2*local_size_x2_sp2, &
            collective_recvcnts_sp2, &
            collective_displs_sp2, &
            0, &
            f_visu_buf1d_sp2 )
          f_visu_sp2 = reshape(f_visu_buf1d_sp2, shape(f_visu_sp2))
          if(sll_get_collective_rank(sll_world_collective)==0) then
            do i=1,num_dof_x2_sp2
              f_visu_buf1d_sp2(i) = sum(f_visu_sp2(1:np_x1-1,i))*sim%mesh2d_sp2%delta_eta1
            enddo
            call sll_gnuplot_1d( &
              f_visu_buf1d_sp2(1:num_dof_x2_sp2), &
              node_positions_x2_sp2(1:num_dof_x2_sp2), &
              'intfedx', &
              iplot )
#ifndef NOHDF5
        call plot_f_cartesian( &
          iplot, &
          f_visu_sp2, &
          sim%x1_array, &
          np_x1, &
          node_positions_x2_sp2, &
          sim%num_dof_x2_sp2, &
          'fi','i', &
          time)                    
#endif
          endif
     
          iplot = iplot+1  
                    
        endif
          

     end if

    
    enddo
    
    ! close files
    if(sll_get_collective_rank(sll_world_collective)==0)then
      call sll_ascii_file_close(th_diag_id,ierr) 
      call sll_ascii_file_close(efield_id,ierr)
      call sll_ascii_file_close(rhotote_id,ierr)
      call sll_ascii_file_close(rhototi_id,ierr)
    endif


    print*, 176.00010668708197, 820.34117552361215 
    print*, sum(f_x1_sp1), sum(f_x1_sp2)
    
    
  end subroutine run_vp2d_cartesian_two_species

  subroutine delete_vp2d_par_cart_two_species( sim )
    class(sll_simulation_2d_vlasov_poisson_cart_two_species) :: sim
    !sll_int32 :: ierr
    
    print *,'#delete_vp2d_par_cart not implemented'
    print *,sim%dt
    
    
  end subroutine delete_vp2d_par_cart_two_species








  elemental function f_equilibrium(v)
    sll_real64, intent(in) :: v
    sll_real64 :: f_equilibrium

    f_equilibrium = 1.0_f64/sqrt(2*sll_pi)*exp(-0.5_f64*v*v)
  end function f_equilibrium

  subroutine plot_f_cartesian( &
    iplot, &
    f, &
    node_positions_x1, &
    nnodes_x1, &
    node_positions_x2, &
    nnodes_x2, &
    array_name, &
    spec_name, &
    time)
    use sll_m_xdmf
    use sll_m_hdf5_io_serial
    sll_int32 :: file_id
    sll_int32 :: error
    sll_real64, dimension(:), intent(in) :: node_positions_x1
    sll_real64, dimension(:), intent(in) :: node_positions_x2    
    character(len=*), intent(in) :: array_name !< field name
    character(len=*), intent(in) :: spec_name !< name of the species
    sll_real64, dimension(:,:), allocatable :: x1
    sll_real64, dimension(:,:), allocatable :: x2
    sll_int32, intent(in) :: nnodes_x1
    sll_int32, intent(in) :: nnodes_x2
    sll_int32 :: i, j
    sll_int32, intent(in) :: iplot
    character(len=4)      :: cplot
    sll_real64, dimension(:,:), intent(in) :: f
    sll_real64 :: time
    
    if (iplot == 1) then

      SLL_ALLOCATE(x1(nnodes_x1,nnodes_x2), error)
      SLL_ALLOCATE(x2(nnodes_x1,nnodes_x2), error)
      do j = 1,nnodes_x2
        do i = 1,nnodes_x1
          x1(i,j) = node_positions_x1(i) !x1_min+real(i-1,f32)*dx1
          x2(i,j) = node_positions_x2(j) !x2_min+real(j-1,f32)*dx2
        end do
     end do
      call sll_hdf5_file_create("cartesian_mesh_"//trim(spec_name)//"-x1.h5",file_id,error)
      call sll_hdf5_write_array(file_id,x1,"/x1",error)
      call sll_hdf5_file_close(file_id, error)
      call sll_hdf5_file_create("cartesian_mesh_"//trim(spec_name)//"-x2.h5",file_id,error)
      call sll_hdf5_write_array(file_id,x2,"/x2",error)
      call sll_hdf5_file_close(file_id, error)
      deallocate(x1)
      deallocate(x2)

    end if

    call int2string(iplot,cplot)
    call sll_xdmf_open(trim(array_name)//cplot//".xmf","cartesian_mesh_"//trim(spec_name), &
      nnodes_x1,nnodes_x2,file_id,error)
    write(file_id,"(a,f8.3,a)") "<Time Value='",time,"'/>"
    call sll_xdmf_write_array(trim(array_name)//cplot,f,"values", &
      error,file_id,"Node")
    call sll_xdmf_close(file_id,error)
  end subroutine plot_f_cartesian

  subroutine compute_e_app(sim,e_app,t)
    class(sll_simulation_2d_vlasov_poisson_cart_two_species) :: sim
    sll_real64, dimension(:), intent(inout) :: e_app
    sll_real64, intent(in) :: t
    sll_real64 :: adr
    sll_int32 :: i
    sll_int32 :: np_x1
    !t = time_init+istep*sim%dt
    e_app = 0._f64
    np_x1 = sim%mesh2d_sp1%num_cells1+1      
    call PFenvelope(adr, &
      t, &
      sim%tflat, &
      sim%tL, &
      sim%tR, &
      sim%twL, &
      sim%twR, &
      sim%t0, &
      sim%turn_drive_off)

    do i = 1, np_x1
     e_app(i) = sim%Edrmax*adr*sim%kx_sp1 &
       * sin(sim%kx_sp1*real(i-1,f64)*sim%mesh2d_sp1%delta_eta1 &
       - sim%omegadr*t)
    enddo

  end subroutine compute_e_app


end module sll_m_sim_bsl_vp_1d1v_cart_two_species
