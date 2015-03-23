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

module sll_simulation_2d_vlasov_ampere_cartesian

#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_field_2d.h"
#include "sll_utilities.h"
#include "sll_poisson_solvers.h"

use sll_collective
use sll_remapper
use sll_buffer_loader_utilities_module
use sll_constants
use sll_cartesian_meshes  
use sll_gnuplot_parallel
use sll_coordinate_transformation_2d_base_module
use sll_module_coordinate_transformations_2d
use sll_common_coordinate_transformations
use sll_common_array_initializers_module
use sll_parallel_array_initializer_module
use sll_module_advection_1d_periodic
use sll_module_advection_1d_non_uniform_cubic_splines
use sll_poisson_1d_periodic  
use sll_fft
use sll_simulation_base
use sll_time_splitting_coeff_module
use sll_module_poisson_1d_periodic_solver
use sll_module_poisson_1d_polar_solver
use sll_module_advection_1d_spectral
use sll_module_advection_1d_ampere

#ifdef _OPENMP
use omp_lib
#endif

  implicit none

  integer, parameter :: SLL_ADVECTIVE    = 0
  integer, parameter :: SLL_CONSERVATIVE = 1

  type, extends(sll_simulation_base_class) :: &
       sll_simulation_2d_vlasov_ampere_cart

   sll_int32                            :: num_threads
   type(sll_cartesian_mesh_2d), pointer :: mesh2d
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
   procedure(sll_scalar_initializer_2d), nopass, pointer :: init_func
   !ES equilibrium function
   procedure(sll_scalar_initializer_2d), nopass, pointer :: equil_func
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

   type(splitting_coeff), pointer :: split
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

   type(sll_advection_1d_base_ptr), dimension(:), pointer :: advect_x1 
   type(sll_advection_1d_base_ptr), dimension(:), pointer :: advect_x2
   type(ampere_1d_advector_ptr),    dimension(:), pointer :: advect_ampere_x1

   sll_int32  :: advection_form_x2
   sll_real64 :: factor_x1
   sll_real64 :: factor_x2_rho
   sll_real64 :: factor_x2_1

   class(sll_poisson_1d_base), pointer :: poisson 
   logical :: ampere = .false.

   sll_real64 :: L
           
   contains !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     procedure, pass(sim) :: run => run_va2d_cartesian
     procedure, pass(sim) :: init_from_file => init_va2d_fake !init_va2d_par_cart

  end type sll_simulation_2d_vlasov_ampere_cart

  interface delete
     module procedure delete_va2d_par_cart
  end interface delete

contains

  function new_va2d_par_cart( filename, num_run ) result(sim)    

    type(sll_simulation_2d_vlasov_ampere_cart), pointer :: sim    
    character(len=*), intent(in), optional               :: filename
    sll_int32, intent(in), optional                      :: num_run

    sll_int32                                            :: ierr

    SLL_ALLOCATE(sim, ierr)
    call init_va2d_par_cart( sim, filename, num_run )
       
  end function new_va2d_par_cart

  subroutine change_initial_function_va2d_par_cart( sim,       &
                                                    init_func, &
                                                    params,    &
                                                    num_params)

    class(sll_simulation_2d_vlasov_ampere_cart)  :: sim
    procedure(sll_scalar_initializer_2d), pointer :: init_func
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

    class(sll_simulation_2d_vlasov_ampere_cart)  :: sim
    procedure(sll_scalar_initializer_2d), pointer :: equil_func
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

    class(sll_simulation_2d_vlasov_ampere_cart) :: sim
    character(len=*), optional                   :: filename
    sll_int32, intent(in), optional :: num_run

    character(len=*), parameter :: this_sub_name = 'init_va2d_par_cart'
    character(len=128)          :: err_msg

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
    type(sll_cartesian_mesh_1d), pointer :: mesh_x1
    type(sll_cartesian_mesh_1d), pointer :: mesh_x2
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

    density_x2_min_to_x2_fine_min      = 1
    density_x2_fine_min_to_x2_fine_max = 1
    density_x2_fine_max_to_x2_max      = 1
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
    !split_case                 = "SLL_ORDER6VPnew1_VTV" 
    !split_case                 = "SLL_ORDER6VPnew2_VTV" 
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
      !call int2string(num_run, str_num_run)
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


      call sll_new_file_id(input_file, ierr)

      open(unit = input_file, file=trim(filename_loc)//'.nml',IOStat=IO_stat)
      if ( IO_stat /= 0 ) then
        err_msg = 'failed to open file '//trim(filename_loc)//'.nml'
        SLL_ERROR( this_sub_name, err_msg )
      end if

      if (sll_get_collective_rank(sll_world_collective)==0) then
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

      if (sll_get_collective_rank(sll_world_collective)==0) then
        print *,'#initialization with default parameters'
      endif      

    endif

    select case (mesh_case_x1) !geometry
      case ("SLL_LANDAU_MESH")
        x1_max = real(nbox_x1,f64) * 2._f64 * sll_pi / kmode
        mesh_x1 => new_cartesian_mesh_1d(num_cells_x1,eta_min=x1_min, eta_max=x1_max)
        call get_node_positions( mesh_x1, sim%x1_array )
      case ("SLL_CARTESIAN_MESH")
        mesh_x1 => new_cartesian_mesh_1d(num_cells_x1,eta_min=x1_min, eta_max=x1_max)  
        call get_node_positions( mesh_x1, sim%x1_array )
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
        mesh_x2 => new_cartesian_mesh_1d(num_cells_x2,eta_min=x2_min, eta_max=x2_max)
        call get_node_positions( mesh_x2, sim%x2_array )
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
                
        call compute_bloc(bloc_coord,bloc_index,num_cells_x2)
        SLL_ALLOCATE(sim%x2_array(num_cells_x2+1),ierr)
        call compute_mesh_from_bloc(bloc_coord,bloc_index,sim%x2_array)
        sim%x2_array = x2_min+sim%x2_array*(x2_max-x2_min)
        SLL_ALLOCATE(sim%x2_array_omp(num_cells_x2+1,0:sim%num_threads-1),ierr)
        do i=0,sim%num_threads-1
          sim%x2_array_omp(:,i) = sim%x2_array(:)
        enddo
        mesh_x2 => new_cartesian_mesh_1d(num_cells_x2,eta_min=x2_min, eta_max=x2_max)
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
        sim%init_func => sll_landau_initializer_2d
        SLL_ALLOCATE(sim%params(2),ierr)
        sim%params(1) = kmode
        sim%params(2) = eps
        sim%nrj0      = 0._f64  !compute the right value
        !(0.5_f64*eps*sll_pi)**2/(kmode_x1*kmode_x2) &
          !*(1._f64/kmode_x1**2+1._f64/kmode_x2**2)
        !for the moment
        sim%kx = kmode
        sim%eps = eps

      case ("SLL_BUMP_ON_TAIL")
        sim%init_func => sll_bump_on_tail_initializer_2d
        SLL_ALLOCATE(sim%params(2),ierr)
        sim%params(1) = kmode
        sim%params(2) = eps
        sim%nrj0 = 0._f64  !compute the right value
        !(0.5_f64*eps*sll_pi)**2/(kmode_x1*kmode_x2) &
          !*(1._f64/kmode_x1**2+1._f64/kmode_x2**2)
        !for the moment
        sim%kx = kmode
        sim%eps = eps

      case ("SLL_TWO_STREAM_INSTABILITY")
        sim%init_func => sll_two_stream_instability_initializer_2d
        SLL_ALLOCATE(sim%params(2),ierr)
        sim%params(1) = kmode
        sim%params(2) = eps
        sim%nrj0 = 0._f64  !compute the right value
        !(0.5_f64*eps*sll_pi)**2/(kmode_x1*kmode_x2) &
          !*(1._f64/kmode_x1**2+1._f64/kmode_x2**2)
        !for the moment
        sim%kx = kmode
        sim%eps = eps

      case ("SLL_BEAM")  
        sim%init_func => sll_beam_initializer_2d
        SLL_ALLOCATE(sim%params(1),ierr)
        sim%params(1) = alpha_gaussian             

      case default
        err_msg = '#init_func_case not implemented'
        SLL_ERROR( this_sub_name, err_msg )

    end select

    !ES Equilibrium is initialised to unperturbed normalised Maxwellian by default
    sim%equil_func => sll_landau_initializer_2d
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

        sim%advect_x1(tid)%ptr => new_periodic_1d_advector( &
          num_cells_x1,                                     &
          x1_min,                                           &
          x1_max,                                           &
          SPLINE,                                           &  
          order_x1) 

      case("SLL_LAGRANGE") ! arbitrary order Lagrange periodic interpolation

        sim%advect_x1(tid)%ptr => new_periodic_1d_advector( &
          num_cells_x1,                                     &
          x1_min,                                           &
          x1_max,                                           &
          LAGRANGE,                                         & 
          order_x1)

      case("SLL_TRIGO") ! trigo periodic advection

        sim%advect_x1(tid)%ptr => new_periodic_1d_advector( &
          num_cells_x1,                                     &
          x1_min,                                           &
          x1_max,                                           &
          TRIGO,                                       &
          order_x1)

      case("SLL_SPECTRAL") ! spectral periodic advection

        sim%advect_x1(tid)%ptr => new_spectral_1d_advector( &
          num_cells_x1,                                     &
          x1_min,                                           &
          x1_max)


      case default

        err_msg = '#advector in x1 '//advector_x1//' not implemented'
        SLL_ERROR( this_sub_name, err_msg )

    end select    

    select case (advector_x2)

      case ("SLL_SPLINES") ! arbitrary order periodic splines

        sim%advect_x2(tid)%ptr => new_periodic_1d_advector( &
          num_cells_x2,                                     &
          x2_min,                                           &
          x2_max,                                           &
          SPLINE,                                           & 
          order_x2) 

      case("SLL_LAGRANGE") ! arbitrary order Lagrange periodic interpolation

        sim%advect_x2(tid)%ptr => new_periodic_1d_advector( &
          num_cells_x2,                                     &
          x2_min,                                           &
          x2_max,                                           &
          LAGRANGE,                                         & 
          order_x2)

      case("SLL_NON_UNIFORM_CUBIC_SPLINES") ! arbitrary order Lagrange 
                                            ! periodic interpolation

        sim%advect_x2(tid)%ptr =>                    &
          new_non_uniform_cubic_splines_1d_advector( &
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
        sim%poisson => new_poisson_1d_periodic_solver( &
          x1_min, &
          x1_max, &
          num_cells_x1)
      case ("SLL_POLAR")
        sim%poisson => new_poisson_1d_polar_solver( &
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
        sim%advect_ampere_x1(tid)%ptr => new_ampere_1d_advector( &
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
        sim%t0             = 0.
        sim%twL            = 0.
        sim%twR            = 0.
        sim%tflat          = 0.
        sim%tL             = 0.
        sim%tR             = 0.
        sim%turn_drive_off = .false.
        sim%Edrmax         = 0.
        sim%omegadr        = 0.    

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

    if (sll_get_collective_rank(sll_world_collective)==0) then
            
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

    class(sll_simulation_2d_vlasov_ampere_cart), intent(inout) :: sim
    character(len=*), intent(in)                               :: filename
  
    character(len=*), parameter :: this_sub_name = 'init_va2d_fake'
    character(len=128)          :: err_msg

    print *,sim%dt
    print *,filename
    err_msg = '# Do not use the routine init_va2d_fake\n'// & 
              '# Use instead init_va2d_par_cart'
    SLL_ERROR( this_sub_name, err_msg )
  
  end subroutine init_va2d_fake

  subroutine run_va2d_cartesian(sim)

    class(sll_simulation_2d_vlasov_ampere_cart), intent(inout) :: sim

    character(len=*), parameter :: this_sub_name = 'run_va2d_cartesian'

    print *,sim%dt
    SLL_WARNING( this_sub_name, '# Do not use the routine run_va2d_fake' )

  end subroutine run_va2d_cartesian

  subroutine delete_va2d_par_cart( sim )

    class(sll_simulation_2d_vlasov_ampere_cart) :: sim
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

    f_equilibrium = 1.0_f64/sqrt(2*sll_pi)*exp(-0.5_f64*v*v)

  end function f_equilibrium

end module sll_simulation_2d_vlasov_ampere_cartesian
