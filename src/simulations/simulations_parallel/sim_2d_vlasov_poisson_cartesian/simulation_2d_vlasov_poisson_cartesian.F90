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

module sll_simulation_2d_vlasov_poisson_cartesian

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
use sll_module_advection_1d_spectral
use sll_poisson_1d_periodic  
use sll_fft
use sll_simulation_base
use sll_time_splitting_coeff_module
use sll_module_poisson_1d_periodic_solver
use sll_module_poisson_1d_polar_solver
use sll_module_ampere_1d_pstd

#ifdef _OPENMP
use omp_lib
#endif

  implicit none

  integer, parameter :: SLL_ADVECTIVE    = 0
  integer, parameter :: SLL_CONSERVATIVE = 1

  type, extends(sll_simulation_base_class) :: &
       sll_simulation_2d_vlasov_poisson_cart

   sll_int32                            :: num_threads
   type(sll_cartesian_mesh_2d), pointer :: mesh2d
   sll_int32                            :: num_dof_x2

   sll_real64, dimension(:),   pointer  :: x1_array
   sll_real64, dimension(:),   pointer  :: x2_array
   sll_real64, dimension(:,:), pointer  :: x2_array_omp
   sll_real64, dimension(:),   pointer  :: x1_array_light
   sll_real64, dimension(:),   pointer  :: x2_array_light
   sll_real64, dimension(:),   pointer  :: integration_weight
   sll_int32,  dimension(:),   pointer  :: every_x1
   sll_int32,  dimension(:),   pointer  :: every_x2
   sll_int32                            :: num_bloc_x1
   sll_int32                            :: num_bloc_x2
   sll_int32,  dimension(:),   pointer  :: bloc_index_x1
   sll_int32,  dimension(:),   pointer  :: bloc_index_x2
   sll_int32                            :: light_size_x1
   sll_int32                            :: light_size_x2
   sll_int32                            :: num_dof_x2_light
   
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

   sll_int32  :: advection_form_x2
   sll_real64 :: factor_x1
   sll_real64 :: factor_x2_rho
   sll_real64 :: factor_x2_1

   class(sll_poisson_1d_base), pointer :: poisson 
   class(sll_ampere_1d_pstd),  pointer :: ampere
           
   contains !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     procedure, pass(sim) :: run => run_vp2d_cartesian
     procedure, pass(sim) :: init_from_file => init_vp2d_fake !init_vp2d_par_cart

  end type sll_simulation_2d_vlasov_poisson_cart

  interface delete
     module procedure delete_vp2d_par_cart
  end interface delete

contains

  function new_vp2d_par_cart( filename ) result(sim)    

    type(sll_simulation_2d_vlasov_poisson_cart), pointer :: sim    
    character(len=*), intent(in), optional               :: filename

    sll_int32                                            :: ierr

    SLL_ALLOCATE(sim, ierr)
    call init_vp2d_par_cart( sim, filename)
       
  end function new_vp2d_par_cart

  subroutine change_initial_function_vp2d_par_cart( sim,       &
                                                    init_func, &
                                                    params,    &
                                                    num_params)

    class(sll_simulation_2d_vlasov_poisson_cart)  :: sim
    procedure(sll_scalar_initializer_2d), pointer :: init_func
    sll_real64, dimension(:), pointer             :: params
    sll_int32, intent(in)                         :: num_params

    sll_int32 :: ierr
    sll_int32 :: i
    
    sim%init_func => init_func
    if (associated(sim%params)) SLL_DEALLOCATE(sim%params,ierr)

    if (num_params<1) then
      SLL_ERROR('#num_params should be >=1')
    endif

    SLL_ALLOCATE(sim%params(num_params),ierr)

    if (size(params)<num_params) then
      SLL_ERROR('#size of params is not good')
    endif

    do i=1,num_params
      sim%params(i) = params(i)
    end do
    
  end subroutine change_initial_function_vp2d_par_cart

  !ES enable externally defined equilibrium function  
  subroutine change_equilibrium_function_vp2d_par_cart( sim,        &
                                                        equil_func, &
                                                        params,     &
                                                        num_params)

    class(sll_simulation_2d_vlasov_poisson_cart)  :: sim
    procedure(sll_scalar_initializer_2d), pointer :: equil_func
    sll_real64, dimension(:), pointer             :: params
    sll_int32, intent(in)                         :: num_params

    sll_int32 :: ierr
    sll_int32 :: i
    
    sim%equil_func => equil_func
    if (associated(sim%equil_func_params)) then
      SLL_DEALLOCATE(sim%equil_func_params,ierr)
    end if

    if (num_params<1) then
      SLL_ERROR('#num_params should be >=1')
    endif

    SLL_ALLOCATE(sim%equil_func_params(num_params),ierr)

    if (size(params)<num_params) then
      SLL_ERROR('#size of params is not good')
    endif

    do i=1,num_params
      sim%equil_func_params(i) = params(i)
    enddo
    
  end subroutine change_equilibrium_function_vp2d_par_cart

  subroutine init_vp2d_par_cart( sim, filename )

    class(sll_simulation_2d_vlasov_poisson_cart) :: sim
    character(len=*), optional                   :: filename

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

    if (present(filename)) then

      call sll_new_file_id(input_file, ierr)

      open(unit = input_file, file=trim(filename)//'.nml',IOStat=IO_stat)
      if ( IO_stat /= 0 ) then
        SLL_ERROR( 'failed to open file '//trim(filename)//'.nml')
      end if

      if (sll_get_collective_rank(sll_world_collective)==0) then
        print *,'#initialization with filename:'
        print *,'#',trim(filename)//'.nml'
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
        SLL_ERROR('#mesh_case_x1 '//mesh_case_x1//' not implemented')
    end select

    sim%num_bloc_x1 = 1
    SLL_ALLOCATE(sim%every_x1(sim%num_bloc_x1),ierr)
    SLL_ALLOCATE(sim%bloc_index_x1(sim%num_bloc_x1),ierr)

    sim%every_x1(1) = every_x1
    sim%bloc_index_x1(sim%num_bloc_x1) = num_cells_x1+1
    sim%light_size_x1 = compute_light_size(num_cells_x1+1,every_x1)   
    SLL_ALLOCATE(sim%x1_array_light(sim%light_size_x1),ierr)

    call compute_light_mesh( num_cells_x1+1,     &
                             sim%x1_array,       &
                             sim%x1_array_light, &
                             sim%light_size_x1,  &
                             every_x1 )
        
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
        sim%light_size_x2 = compute_light_size(num_cells_x2+1,every_x2)   
        !sim%light_size_x1 = compute_light_size_bloc(sim%bloc_index_x1,sim%every_x1)   
        SLL_ALLOCATE(sim%x2_array_light(sim%light_size_x2),ierr)

        call compute_light_mesh( num_cells_x2+1,     &
                                 sim%x2_array,       & 
                                 sim%x2_array_light, &
                                 sim%light_size_x2,  &
                                 every_x2 )

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
        sim%light_size_x2 = compute_light_size_3_bloc( bloc_index, sim%every_x2 )
        SLL_ALLOCATE(sim%x2_array_light(sim%light_size_x2),ierr)

        call compute_light_mesh_3_bloc( bloc_index,         &
                                        sim%x2_array,       &
                                        sim%x2_array_light, &
                                        sim%light_size_x2,  &
                                        sim%every_x2 )
      case default

        SLL_ERROR('#mesh_case_x2 '//mesh_case_x2//' not implemented')

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

        SLL_ERROR('#init_func_case not implemented')

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
      print *,'#bad value of nb_mode=',nb_mode      
      SLL_ERROR('#should be >=0')
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
        SLL_ERROR('#split_case not defined')
    end select

    !advector
    SLL_ALLOCATE(sim%advect_x1(num_threads),ierr)
    SLL_ALLOCATE(sim%advect_x2(num_threads),ierr)
    tid = 1

    !$OMP PARALLEL DEFAULT(SHARED) &
    !$OMP PRIVATE(tid)

#ifdef _OPENMP
    tid = omp_get_thread_num()+1
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

      case("SLL_SPECTRAL") ! spectral advection BC must be periodic

        sim%advect_x1(tid)%ptr => new_spectral_1d_advector( &
          num_cells_x1,                                     &
          x1_min,                                           &
          x1_max,                                           &
          SLL_PERIODIC)

      case default

        SLL_ERROR('#advector in x1 '//advector_x1//' not implemented')
        stop 

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

        SLL_ERROR('#advector in x2 '//advector_x2//' not implemented')

    end select

    !$OMP END PARALLEL

    select case (advection_form_x2)
      case ("SLL_ADVECTIVE")
        sim%advection_form_x2 = SLL_ADVECTIVE
        sim%num_dof_x2        = num_cells_x2+1
        sim%num_dof_x2_light  = sim%light_size_x2
      case ("SLL_CONSERVATIVE")
        sim%advection_form_x2 = SLL_CONSERVATIVE
        sim%num_dof_x2        = num_cells_x2
        sim%num_dof_x2_light  = sim%light_size_x2-1
      case default
        SLL_ERROR('#advection_form_x2 '//advection_form_x2//' not implemented')
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
        SLL_ERROR('#integration_case not implemented')
    end select  
    
    !poisson
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
      case ("SLL_AMPERE_PSTD")
        sim%ampere => new_ampere_1d_pstd( &
          x1_min, &
          x1_max, &
          num_cells_x1)
      case default
        SLL_ERROR('#poisson_solver '//poisson_solver//' not implemented')
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
        SLL_ERROR('#drive_type '//drive_type//' not implemented')

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

  end subroutine init_vp2d_par_cart

  subroutine init_vp2d_fake(sim, filename)

    class(sll_simulation_2d_vlasov_poisson_cart), intent(inout) :: sim
    character(len=*), intent(in)                                :: filename
  
    print *,sim%dt
    print *,filename
    SLL_WARNING('# Do not use the routine init_vp2d_fake')
    SLL_ERROR('#use instead init_vp2d_par_cart')
  
  end subroutine init_vp2d_fake

  subroutine run_vp2d_cartesian(sim)

    class(sll_simulation_2d_vlasov_poisson_cart), intent(inout) :: sim

    sll_real64, dimension(:,:), pointer :: f_x1
    sll_real64, dimension(:,:), pointer :: f_x2
    sll_real64, dimension(:,:), pointer :: f_x1_equil

    sll_real64, dimension(:),   pointer :: rho
    sll_real64, dimension(:),   pointer :: efield
    sll_real64, dimension(:),   pointer :: e_app
    sll_real64, dimension(:),   pointer :: rho_loc
    sll_real64, dimension(:),   pointer :: current
    
    sll_int32 :: rhotot_id
    sll_int32 :: efield_id     
    sll_int32 :: adr_id
    sll_int32 :: Edr_id
    sll_int32 :: deltaf_id
    sll_int32 :: t_id
    sll_int32 :: th_diag_id
    sll_int32 :: restart_id

    type(layout_2D),            pointer :: layout_x1
    type(layout_2D),            pointer :: layout_x2
    type(remap_plan_2D_real64), pointer :: remap_plan_x1_x2
    type(remap_plan_2D_real64), pointer :: remap_plan_x2_x1
    sll_real64, dimension(:),   pointer :: f1d
    sll_real64, dimension(:,:), pointer :: f1d_omp
    sll_real64, dimension(:,:), pointer :: f1d_omp_in
    sll_real64, dimension(:,:), pointer :: f1d_omp_out
    sll_int32                           :: np_x1
    sll_int32                           :: np_x2
    sll_int32                           :: np_x1_light
    sll_int32                           :: nproc_x1
    sll_int32                           :: nproc_x2
    sll_int32                           :: global_indices(2)
    sll_int32                           :: ierr
    sll_int32                           :: local_size_x1
    sll_int32                           :: local_size_x2
    type(poisson_1d_periodic)           :: poisson_1d
    sll_real64                          :: adr
    sll_real64                          :: tmp_loc(5)
    sll_real64                          :: tmp(5)
    sll_int32                           :: i
    sll_int32                           :: j
    sll_int32                           :: istep
    sll_int32                           :: ig
    sll_int32                           :: k
    
    sll_real64  ::   time, mass, momentum, l1norm, l2norm
    sll_real64  ::   kinetic_energy,potential_energy

    sll_real64, dimension(:), allocatable :: x2_array_unit
    sll_real64, dimension(:), allocatable :: x2_array_middle
    sll_real64, dimension(:), allocatable :: node_positions_x2
    sll_real64, dimension(:), allocatable :: node_positions_x2_light

    sll_int32                             :: file_id
    
    type(sll_fft_plan), pointer           :: pfwd
    sll_real64, dimension(:), allocatable :: buf_fft
    sll_comp64,dimension(:),allocatable   :: rho_mode
    
    sll_int32  :: nb_mode 
    sll_real64 :: t_step
    sll_real64 :: time_init
    sll_int32  :: split_istep
    sll_int32  :: num_dof_x2 
    sll_int32  :: num_dof_x2_light
     
    logical :: split_T
    
    ! for parallelization (output of distribution function in one single file)
    sll_int32                              :: collective_size
    sll_int32, dimension(:),   allocatable :: collective_displs
    sll_int32, dimension(:),   allocatable :: collective_recvcnts
    sll_real64,dimension(:,:), pointer     :: f_visu 
    sll_real64,dimension(:,:), pointer     :: f_visu_light 
    sll_real64,dimension(:),   pointer     :: f_visu_buf1d
    sll_real64,dimension(:),   pointer     :: f_x1_buf1d
    sll_real64,dimension(:),   pointer     :: f_hat_x2_loc
    sll_real64,dimension(:),   pointer     :: f_hat_x2
    sll_int32                              :: iplot
    character(len=4)                       :: cproc
    character(len=4)                       :: cplot
    sll_int32                              :: iproc
    logical                                :: file_exists
    sll_int32                              :: tid
    sll_int32                              :: i_omp
    sll_int32                              :: ig_omp
    sll_real64                             :: alpha_omp
    sll_real64                             :: mean_omp

    logical                                :: MPI_MASTER

    if (sll_get_collective_rank(sll_world_collective)==0) then
      MPI_MASTER = .true.
    else
      MPI_MASTER = .false.
    end if

    iplot = 1

    nb_mode          = sim%nb_mode
    time_init        = sim%time_init
    np_x1            = sim%mesh2d%num_cells1+1
    np_x2            = sim%mesh2d%num_cells2+1
    num_dof_x2       = sim%num_dof_x2
    num_dof_x2_light = sim%num_dof_x2_light
    np_x1_light      = sim%light_size_x1
    
    if (MPI_MASTER) then

      print *,'#collective_size=',sll_get_collective_size(sll_world_collective)
      SLL_ALLOCATE(f_visu(np_x1,num_dof_x2),ierr)
      SLL_ALLOCATE(f_visu_buf1d(np_x1*num_dof_x2),ierr)
      SLL_ALLOCATE(f_visu_light(np_x1_light,num_dof_x2_light),ierr)

    else

      SLL_ALLOCATE(f_visu(1:1,1:1),ierr)          
      SLL_ALLOCATE(f_visu_buf1d(1:1),ierr)
      SLL_ALLOCATE(f_visu_light(1:1,1:1),ierr)          

    endif
    
    collective_size = sll_get_collective_size(sll_world_collective)
    SLL_ALLOCATE(collective_displs(collective_size),ierr)
    SLL_ALLOCATE(collective_recvcnts(collective_size),ierr)

    SLL_ALLOCATE(buf_fft(np_x1-1),ierr)
    pfwd => fft_new_plan(np_x1-1,buf_fft,buf_fft,FFT_FORWARD,FFT_NORMALIZE)
    SLL_ALLOCATE(rho_mode(0:nb_mode),ierr)      

    layout_x1       => new_layout_2D( sll_world_collective )
    layout_x2       => new_layout_2D( sll_world_collective )    
    nproc_x1 = sll_get_collective_size( sll_world_collective )
    nproc_x2 = 1
    call initialize_layout_with_distributed_array( &
      np_x1, num_dof_x2, nproc_x1, nproc_x2, layout_x2 )
    call initialize_layout_with_distributed_array( &
      np_x1, num_dof_x2, nproc_x2, nproc_x1, layout_x1 )

    call sll_view_lims( layout_x1 )
    call sll_view_lims( layout_x2 )
    
    call compute_local_sizes( layout_x2, local_size_x1, local_size_x2 )
    SLL_ALLOCATE(f_x2(local_size_x1,local_size_x2),ierr)

    call compute_local_sizes( layout_x1, local_size_x1, local_size_x2 )
    global_indices(1:2) = local_to_global( layout_x1, (/1, 1/) )
    SLL_ALLOCATE(f_x1(local_size_x1,local_size_x2),ierr)    
    SLL_ALLOCATE(f_x1_equil(local_size_x1,local_size_x2),ierr)    
    SLL_ALLOCATE(f_x1_buf1d(local_size_x1*local_size_x2),ierr)    

    remap_plan_x1_x2 => NEW_REMAP_PLAN(layout_x1, layout_x2, f_x1)
    remap_plan_x2_x1 => NEW_REMAP_PLAN(layout_x2, layout_x1, f_x2)
    
    SLL_ALLOCATE(rho(np_x1),ierr)
    SLL_ALLOCATE(rho_loc(np_x1),ierr)
    SLL_ALLOCATE(current(np_x1),ierr)
    SLL_ALLOCATE(efield(np_x1),ierr)
    SLL_ALLOCATE(e_app(np_x1),ierr)
    SLL_ALLOCATE(f1d(max(np_x1,np_x2)),ierr)
    SLL_ALLOCATE(f1d_omp(max(np_x1,np_x2),sim%num_threads),ierr)
    SLL_ALLOCATE(f1d_omp_in(max(np_x1,np_x2),sim%num_threads),ierr)
    SLL_ALLOCATE(f1d_omp_out(max(np_x1,np_x2),sim%num_threads),ierr)
    SLL_ALLOCATE(x2_array_unit(np_x2),ierr)
    SLL_ALLOCATE(x2_array_middle(np_x2),ierr)
    SLL_ALLOCATE(node_positions_x2(num_dof_x2),ierr)
    SLL_ALLOCATE(node_positions_x2_light(num_dof_x2),ierr)
    SLL_ALLOCATE(f_hat_x2_loc(nb_mode+1),ierr)
    SLL_ALLOCATE(f_hat_x2(nb_mode+1),ierr)

    x2_array_unit(1:np_x2) = &
      (sim%x2_array(1:np_x2)-sim%x2_array(1))/(sim%x2_array(np_x2)-sim%x2_array(1))
    do i = 1, np_x2-1
       x2_array_middle(i) = 0.5_f64*(sim%x2_array(i)+sim%x2_array(i+1))
    end do
    x2_array_middle(np_x2) = x2_array_middle(1)+sim%x2_array(np_x2)-sim%x2_array(1)
    
    select case (sim%advection_form_x2)
      case (SLL_ADVECTIVE)
        node_positions_x2(1:num_dof_x2) = sim%x2_array(1:num_dof_x2)
      case (SLL_CONSERVATIVE)
        node_positions_x2(1:num_dof_x2) = x2_array_middle(1:num_dof_x2)
      case default
        print *,'#sim%advection_form_x2=',sim%advection_form_x2
        print *,'#not implemented'
        print *,'#in run_vp2d_cartesian'
        stop
    end select  
        

    if (MPI_MASTER) then

      call sll_binary_file_create("x.bdat", file_id, ierr)
      call sll_binary_write_array_1d(file_id,sim%x1_array(1:np_x1-1),ierr)
      call sll_binary_file_close(file_id,ierr)                    
      call sll_binary_file_create("v.bdat", file_id, ierr)
      call sll_binary_write_array_1d(file_id,node_positions_x2(1:np_x2-1),ierr)
      call sll_binary_file_close(file_id,ierr)                                             
    endif

    call sll_2d_parallel_array_initializer_cartesian( &
       layout_x1,                                     &
       sim%x1_array,                                  &
       node_positions_x2,                             &
       f_x1,                                          &
       sim%init_func,                                 &
       sim%params)

    iproc = sll_get_collective_rank(sll_world_collective)
    call int2string(iproc, cproc)
    call int2string(iplot, cplot)    

    if (trim(sim%restart_file) /= "no_restart_file" ) then

      INQUIRE(FILE=trim(sim%restart_file)//'_proc_'//cproc//'.rst', EXIST=file_exists)

      if (.not. file_exists) then
        SLL_ERROR('#file '//trim(sim%restart_file)//'_proc_'//cproc//'.rst &
                  & does not exist')
      endif

      open(unit=restart_id, &
           file=trim(sim%restart_file)//'_proc_'//cproc//'.rst', ACCESS="STREAM", &
           form='unformatted', IOStat=ierr)      

      if ( ierr .ne. 0 ) then
        SLL_ERROR('ERROR while opening file &
                   &'//trim(sim%restart_file)//'_proc_'//cproc//'.rst &
                   & Called from run_vp2d_cartesian().')
      end if

      print *,'#read restart file '//trim(sim%restart_file)//'_proc_'//cproc//'.rst'      
      call sll_binary_read_array_0d(restart_id,time_init,ierr)
      call sll_binary_read_array_2d(restart_id,f_x1,ierr)
      call sll_binary_file_close(restart_id,ierr)

    endif      
    
    if (sim%time_init_from_restart_file .eqv. .true.) then
      sim%time_init = time_init  
    endif
    time_init = sim%time_init

    !ES initialise f_x1_equil that is used to define deltaf: 
    !ES deltaf =  f_x1 - f_x1_equil
    call sll_2d_parallel_array_initializer_cartesian( &
       layout_x1,                                     &
       sim%x1_array,                                  &
       node_positions_x2,                             &
       f_x1_equil,                                    &
       sim%equil_func,                                &
       sim%equil_func_params)

    call compute_displacements_array_2d( layout_x1,       &
                                         collective_size, &
                                         collective_displs )

    collective_recvcnts = receive_counts_array_2d( layout_x1, &
                                                   collective_size )

    !ES replace f_x1 by f_x1_equil and write f_x1_equil into 
    !ES f0.bdat instead of f_x1  

    call load_buffer_2d( layout_x1, f_x1_equil, f_x1_buf1d )

    call sll_collective_gatherv_real64( sll_world_collective,        &
                                        f_x1_buf1d,                  &
                                        local_size_x1*local_size_x2, &
                                        collective_recvcnts,         &
                                        collective_displs,           &
                                        0,                           &
                                        f_visu_buf1d )
        
    f_visu = reshape(f_visu_buf1d, shape(f_visu))

    if (MPI_MASTER) then

      call sll_binary_file_create('f0.bdat', file_id, ierr)
      call sll_binary_write_array_2d(file_id,f_visu(1:np_x1-1,1:np_x2-1),ierr)
      call sll_binary_file_close(file_id,ierr)

      call plot_f_cartesian( iplot,             &
                             f_visu,            &
                             sim%x1_array,      &
                             np_x1,             &
                             node_positions_x2, &
                             sim%num_dof_x2,    &
                             'f',               &
                             time_init )        

      print *,'#maxf',maxval(f_visu), minval(f_visu) 

    endif

    call initialize(poisson_1d,sim%mesh2d%eta1_min,sim%mesh2d%eta1_max,np_x1-1,ierr)

    !computation of electric field
    rho_loc = 0._f64
    ig = global_indices(2)-1
    do i=1,np_x1
      rho_loc(i) = rho_loc(i)&
        +sum(f_x1(i,:)*sim%integration_weight(1+ig:local_size_x2+ig))
    end do

    call sll_collective_allreduce( sll_world_collective, &
                                   rho_loc,              &
                                   np_x1,                &
                                   MPI_SUM,              &
                                   rho )

    rho = sim%factor_x2_1*1._f64-sim%factor_x2_rho*rho        
    
    call sim%poisson%compute_E_from_rho( efield, rho )
        
    ! Ponderomotive force at initial time. We use a sine wave
    ! with parameters k_dr and omega_dr.
    istep = 0
    e_app = 0._f64

    if (sim%driven) then

      call PFenvelope(adr,                    &
                      time_init+istep*sim%dt, &
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
                 - sim%omegadr*(time_init+real(istep,f64)*sim%dt))
      enddo

    endif

    if (MPI_MASTER) then

      call sll_ascii_file_create('thdiag.dat', th_diag_id, ierr)

      call sll_binary_file_create('deltaf.bdat', deltaf_id, ierr)
      call sll_binary_file_create('rhotot.bdat', rhotot_id, ierr)
      call sll_binary_file_create('efield.bdat', efield_id, ierr)
      call sll_binary_file_create('t.bdat', t_id, ierr)
      call sll_binary_write_array_1d(efield_id,efield(1:np_x1-1),ierr)
      call sll_binary_write_array_1d(rhotot_id,rho(1:np_x1-1),ierr)
      call sll_binary_write_array_0d(t_id,real(istep,f64)*sim%dt,ierr)

      if (sim%driven) then
        call sll_binary_file_create('adr.bdat', adr_id, ierr)
        call sll_binary_file_create('Edr.bdat', Edr_id, ierr)
        call sll_binary_write_array_1d(Edr_id,e_app(1:np_x1-1),ierr)
        call sll_binary_write_array_0d(adr_id,adr,ierr)
      endif                    

    endif
    
    !write also initial deltaf function
    call load_buffer_2d( layout_x1, f_x1-f_x1_equil, f_x1_buf1d )

    call sll_collective_gatherv_real64( sll_world_collective,        &
                                        f_x1_buf1d,                  &
                                        local_size_x1*local_size_x2, &
                                        collective_recvcnts,         &
                                        collective_displs,           &
                                        0,                           &
                                        f_visu_buf1d )

    f_visu = reshape(f_visu_buf1d, shape(f_visu))

    if (MPI_MASTER) then
      call sll_binary_write_array_2d(deltaf_id,f_visu(1:np_x1-1,1:np_x2-1),ierr)
      print *,'#step=',0,time_init+real(0,f64)*sim%dt,'iplot=',iplot
    endif

    iplot = iplot+1  

    do istep = 1, sim%num_iterations

      if (mod(istep,sim%freq_diag)==0) then
        if (MPI_MASTER) then        
          print *,'#step=',istep,time_init+real(istep,f64)*sim%dt,'iplot=',iplot
        endif
      endif  

      split_T = sim%split%split_begin_T
      t_step = real(istep-1,f64)

      do split_istep=1,sim%split%nb_split_step

        if (split_T) then
          !! T ADVECTION 
          tid=1          
          !$OMP PARALLEL DEFAULT(SHARED) &
          !$OMP PRIVATE(i_omp,ig_omp,alpha_omp,tid) 
          !advection in x
#ifdef _OPENMP
          tid = omp_get_thread_num()+1
#endif
          !$OMP DO
          do i_omp = 1, local_size_x2

            ig_omp    = i_omp+global_indices(2)-1
            alpha_omp = sim%factor_x1*node_positions_x2(ig_omp) &
                      * sim%split%split_step(split_istep) 
            f1d_omp_in(1:np_x1,tid) = f_x1(1:np_x1,i_omp)
            
            call sim%advect_x1(tid)%ptr%advect_1d_constant(  &
              alpha_omp,                                     &
              sim%dt,                                        &
              f1d_omp_in(1:np_x1,tid),                       &
              f1d_omp_out(1:np_x1,tid))

            f_x1(1:np_x1,i_omp)=f1d_omp_out(1:np_x1,tid)

          end do
          !$OMP END DO          
          !$OMP END PARALLEL

          t_step = t_step+sim%split%split_step(split_istep)

          if (associated(sim%ampere)) then
            !Compute current
            rho_loc = 0._f64
            do i=1,np_x1
              do j= 1, local_size_x2
                global_indices(1:2) = local_to_global( layout_x1, (/i, j/) )
                rho_loc(i)=rho_loc(i)                                &
                +f_x1(i,j)*sim%integration_weight(global_indices(2)) &
                *node_positions_x2(global_indices(2))
              end do
            end do
                
            call sll_collective_allreduce( sll_world_collective, &
                                           rho_loc,              &
                                           np_x1,                &
                                           MPI_SUM,              &
                                           current )

            call sim%ampere%compute_e_from_j( sim%dt, current, efield )

          else

            !computation of electric field
            rho_loc = 0._f64
            ig = global_indices(2)-1
            do i=1,np_x1
              rho_loc(i)=rho_loc(i)                              &
                +sum(f_x1(i,1:local_size_x2)                     &
                *sim%integration_weight(1+ig:local_size_x2+ig))
            end do
                
            call sll_collective_allreduce( sll_world_collective, &
                                           rho_loc,              &
                                           np_x1,                &
                                           MPI_SUM,              &
                                           rho )
  
            rho = sim%factor_x2_1*1._f64-sim%factor_x2_rho*rho
  
            call sim%poisson%compute_E_from_rho( efield, rho )

          end if
          
          if (sim%driven) then

            call PFenvelope(adr,                     &
                            time_init+t_step*sim%dt, &
                            sim%tflat,               &
                            sim%tL,                  &
                            sim%tR,                  &
                            sim%twL,                 &
                            sim%twR,                 &
                            sim%t0,                  &
                            sim%turn_drive_off)

            do i = 1, np_x1
              e_app(i) = sim%Edrmax*adr*sim%kx      &
                       * sin(sim%kx*real(i-1,f64)   &
                       * sim%mesh2d%delta_eta1      &
                       - sim%omegadr*(time_init+t_step*sim%dt))
            enddo

          endif

        else

          !! V ADVECTION 
          !transposition
          call apply_remap_2D( remap_plan_x1_x2, f_x1, f_x2 )
          call compute_local_sizes( layout_x2, local_size_x1, local_size_x2 )
          global_indices(1:2) = local_to_global( layout_x2, (/1, 1/) )
          tid = 1

          !$OMP PARALLEL DEFAULT(SHARED) &
          !$OMP PRIVATE(i_omp,ig_omp,alpha_omp,tid,mean_omp,f1d) 
          !advection in v
#ifdef _OPENMP
          tid = omp_get_thread_num()+1
#endif
          !advection in v
          !$OMP DO
          do i_omp = 1,local_size_x1

            ig_omp=i_omp+global_indices(1)-1

            alpha_omp = -(efield(ig_omp)+e_app(ig_omp)) &
                      * sim%split%split_step(split_istep)

            f1d_omp_in(1:num_dof_x2,tid) = f_x2(i_omp,1:num_dof_x2)

            if (sim%advection_form_x2==SLL_CONSERVATIVE) then

              call function_to_primitive(f1d_omp_in(:,tid),    &
                                         x2_array_unit,        &
                                         np_x2-1,mean_omp)
            endif

            call sim%advect_x2(tid)%ptr%advect_1d_constant(    &
              alpha_omp,                                       &
              sim%dt,                                          &
              f1d_omp_in(1:num_dof_x2,tid),                    &
              f1d_omp_out(1:num_dof_x2,tid))

            if (sim%advection_form_x2==SLL_CONSERVATIVE) then

              call primitive_to_function(f1d_omp_out(:,tid),   &
                                         x2_array_unit,        &
                                         np_x2-1,              &
                                         mean_omp)
            endif
            f_x2(i_omp,1:num_dof_x2) = f1d_omp_out(1:num_dof_x2,tid)
          end do
          !$OMP END DO          
          !$OMP END PARALLEL

          !transposition
          call apply_remap_2D( remap_plan_x2_x1, f_x2, f_x1 )
          call compute_local_sizes( layout_x1, local_size_x1, local_size_x2 )
          global_indices(1:2) = local_to_global( layout_x1, (/1, 1/) )

        endif
        !split_x= 1-split_x
        split_T = .not.(split_T)
      enddo
      
     !!DIAGNOSTICS
     if (mod(istep,sim%freq_diag_time)==0) then

        time             = time_init+real(istep,f64)*sim%dt
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
        
        call sll_collective_allreduce( sll_world_collective, &
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
        
        f_hat_x2_loc(1:nb_mode+1) = 0._f64
        do i=1,local_size_x2
          buf_fft = f_x1(1:np_x1-1,i)
          call fft_apply_plan(pfwd,buf_fft,buf_fft)
          do k=0,nb_mode
            f_hat_x2_loc(k+1) = f_hat_x2_loc(k+1) &
              +abs(fft_get_mode(pfwd,buf_fft,k))**2 &
              *sim%integration_weight(ig+i)
          enddo
        enddo

        call sll_collective_allreduce( sll_world_collective, &
                                       f_hat_x2_loc,         &
                                       nb_mode+1,            &
                                       MPI_SUM,              &
                                       f_hat_x2 )
        
        if (mod(istep,sim%freq_diag_restart)==0) then          
          call int2string(iplot,cplot) 
          call sll_binary_file_create('f_plot_'//cplot//'_proc_'//cproc//'.rst', &
                                      restart_id, ierr )
          call sll_binary_write_array_0d(restart_id,time,ierr)
          call sll_binary_write_array_2d(restart_id, &
                                         f_x1(1:local_size_x1,1:local_size_x2),ierr)
          call sll_binary_file_close(restart_id,ierr)    
        endif 

        if (sll_get_collective_rank(sll_world_collective)==0) then                  

          buf_fft = rho(1:np_x1-1)
          call fft_apply_plan(pfwd,buf_fft,buf_fft)

          do k=0,nb_mode
            rho_mode(k)=fft_get_mode(pfwd,buf_fft,k)
          enddo  

          write(th_diag_id,'(f12.5,7g20.12)',advance='no') &
            time,                                          &
            mass,                                          &
            l1norm,                                        &
            momentum,                                      &
            l2norm,                                        &
            kinetic_energy,                                &
            potential_energy,                              &
            kinetic_energy + potential_energy

          do k=0,nb_mode
            write(th_diag_id,'(g20.12)',advance='no') abs(rho_mode(k))
          enddo

          do k=0,nb_mode-1
            write(th_diag_id,'(g20.12)',advance='no') f_hat_x2(k+1)
          enddo

          write(th_diag_id,'(g20.12)') f_hat_x2(nb_mode+1)

          call sll_binary_write_array_1d(efield_id,efield(1:np_x1-1),ierr)
          call sll_binary_write_array_1d(rhotot_id,rho(1:np_x1-1),ierr)
          call sll_binary_write_array_0d(t_id,time,ierr)
          if (sim%driven) then
            call sll_binary_write_array_1d(Edr_id,e_app(1:np_x1-1),ierr)
            call sll_binary_write_array_0d(adr_id,adr,ierr)
          endif   
        endif
          
        if (mod(istep,sim%freq_diag)==0) then          
          !we substract f0
          !we gather in one file
          call load_buffer_2d( layout_x1, f_x1-f_x1_equil, f_x1_buf1d )

          call sll_collective_gatherv_real64( &
            sll_world_collective,             &
            f_x1_buf1d,                       &
            local_size_x1*local_size_x2,      &
            collective_recvcnts,              &
            collective_displs,                &
            0,                                &
            f_visu_buf1d )

          f_visu = reshape(f_visu_buf1d, shape(f_visu))

          if (MPI_MASTER) then
            do i=1,num_dof_x2
              f_visu_buf1d(i) = sum(f_visu(1:np_x1-1,i))*sim%mesh2d%delta_eta1
            enddo

            call sll_gnuplot_write_1d(         &
              f_visu_buf1d(1:num_dof_x2),      &
              node_positions_x2(1:num_dof_x2), &
              'intdeltafdx',                   &
              iplot )

            call sll_gnuplot_write_1d(         &
              f_visu_buf1d(1:num_dof_x2),      &
              node_positions_x2(1:num_dof_x2), &
              'intdeltafdx')                        

            call sll_binary_write_array_2d(deltaf_id,           &
                                           f_visu(1:np_x1-1,1:np_x2-1),ierr)  

            call plot_f_cartesian( &
              iplot,               &
              f_visu,              &
              sim%x1_array,        &
              np_x1,               &
              node_positions_x2,   &
              sim%num_dof_x2,      &
              'deltaf',time)                    
          
          
          endif
          !we store f for visu
          call load_buffer_2d( layout_x1, f_x1, f_x1_buf1d )
          call sll_collective_gatherv_real64( &
            sll_world_collective,             &
            f_x1_buf1d,                       &
            local_size_x1*local_size_x2,      &
            collective_recvcnts,              &
            collective_displs,                &
            0,                                &
            f_visu_buf1d )

          f_visu = reshape(f_visu_buf1d, shape(f_visu))

          if (MPI_MASTER) then

            do i=1,num_dof_x2
              f_visu_buf1d(i) = sum(f_visu(1:np_x1-1,i))*sim%mesh2d%delta_eta1
            enddo

            call sll_gnuplot_write_1d(         &
              f_visu_buf1d(1:num_dof_x2),      &
              node_positions_x2(1:num_dof_x2), &
              'intfdx',                        &
              iplot )

            call sll_gnuplot_write_1d(         &
              f_visu_buf1d(1:num_dof_x2),      &
              node_positions_x2(1:num_dof_x2), &
              'intfdx')

            call plot_f_cartesian( iplot,             &
                                   f_visu,            &
                                   sim%x1_array,      &
                                   np_x1,             &
                                   node_positions_x2, &
                                   sim%num_dof_x2,    &
                                   'f', time)                    
          endif

          iplot = iplot+1  
                    
        endif

     end if

    enddo
    
    if (MPI_MASTER) then
      call sll_ascii_file_close(th_diag_id,ierr) 
      call sll_binary_file_close(deltaf_id,ierr) 
      call sll_binary_file_close(efield_id,ierr)
      call sll_binary_file_close(rhotot_id,ierr)
      call sll_binary_file_close(t_id,ierr)
      if (sim%driven) then
        call sll_binary_file_close(Edr_id,ierr)
        call sll_binary_file_close(adr_id,ierr)
      endif   
    endif

  end subroutine run_vp2d_cartesian

  subroutine delete_vp2d_par_cart( sim )

    class(sll_simulation_2d_vlasov_poisson_cart) :: sim
    
    print *,sim%dt
    SLL_WARNING('#delete_vp2d_par_cart not implemented')
    
  end subroutine delete_vp2d_par_cart

  subroutine function_to_primitive(f,node_positions,N,M)

    sll_real64,dimension(:),intent(inout) :: f
    sll_real64,dimension(:),intent(in)    :: node_positions
    sll_int32,              intent(in)    :: N
    sll_real64,             intent(out)   :: M

    sll_int32  :: i
    sll_real64 :: dx,tmp,tmp2

    dx = 1._f64/real(N,f64)
        
    !from f compute the mean
    M=0._f64
    do i=1,N
      M=M+f(i)*(node_positions(i+1)-node_positions(i))
    enddo
    
    f(1)=(f(1)-M)*(node_positions(2)-node_positions(1))
    tmp=f(1)
    f(1)=0._f64
    do i=2,N!+1
      f(i)=(f(i)-M)*(node_positions(i+1)-node_positions(i))
      tmp2=f(i)
      f(i)=f(i-1)+tmp
      tmp=tmp2
    enddo    
    f(N+1)=f(N)+tmp

  end subroutine function_to_primitive

  subroutine primitive_to_function(f,node_positions,N,M)
    sll_real64, dimension(:), intent(inout) :: f
    sll_real64, dimension(:), intent(in)    :: node_positions
    sll_int32,                intent(in)    :: N
    sll_real64,               intent(in)    :: M

    sll_int32                               :: i
    sll_real64                              :: tmp
    sll_real64                              :: dx 

    dx = 1._f64/real(N,f64)
    
    tmp=f(1)
    do i=1,N-1
      f(i)=f(i+1)-f(i)+M*(node_positions(i+1)-node_positions(i))
    enddo
    f(N)=tmp-f(N)+M*(node_positions(1)+1._f64-node_positions(N))

    !from mean compute f
    do i=1,N
      f(i)=f(i)/(node_positions(i+1)-node_positions(i))
    enddo

    f(N+1) = f(1)

  end subroutine primitive_to_function

  ! Input: 
  !  a=bloc_coord(1) b=bloc_coord(2)
  !  (a,b) \subset (0,1) is the refine zone
  !  bloc_index(1) = density of points in (0,a) 
  !  bloc_index(2) = density of points in (a,b) 
  !  bloc_index(3) = density of points in (b,1)
  !
  ! Output:
  !  0<=i1<i1+N_fine<=N and x(i1)=a, x(i1+N_fine)=b (approx), x(0)=0, x(N)=1
  !  bloc_coord(1) = x(i1)
  !  bloc_coord(2) = x(i1+N_fine)
  !  bloc_index(1) = i1 
  !  bloc_index(2) = N_fine 
  !  bloc_index(3) = N-i1-N_fine

  subroutine compute_bloc(bloc_coord,bloc_index,N)

    sll_real64, intent(inout)  :: bloc_coord(2)
    sll_int32,  intent(inout)  :: bloc_index(3)
    sll_int32,  intent(in)     :: N

    sll_real64 :: a,b
    sll_int32  :: i1,i2,N_coarse,N_local,N_fine
    
    a=bloc_coord(1)
    b=bloc_coord(2)
    
    !case of uniform mesh with refined zone
    !we have a coarse mesh with N_coarse
    !N=i1+N_local*(i2-i1)+N_coarse-i2
    !N_fine=N_local*(i2-i1)
    !x(i1)=i1/N_coarse x(i1+N_fine)=i2/N_coarse

    if ((bloc_index(1)==1).and.(bloc_index(3)==1)) then      

      N_local = bloc_index(2)
      N_coarse = floor(real(N,f64)/(1._f64+(b-a)*(real(N_local,f64)-1._f64)))
      if (N_local/=1) then
        i2 = (N-N_coarse)/(N_local-1)
      else
        i2 = floor((b-a)*N_coarse)  
      endif   
      N_coarse      = N-i2*(N_local-1)
      i1            = floor(a*N_coarse)
      i2            = i2+i1
      bloc_index(1) = i1
      N_fine        = N_local*(i2-i1)
      bloc_index(2) = N_fine
      bloc_index(3) = N-i1-N_fine
      bloc_coord(1) = real(i1,f64)/real(N_coarse,f64)
      bloc_coord(2) = real(i2,f64)/real(N_coarse,f64)
           
      print *,'#uniform fine mesh would be:',N_coarse*N_local
      print *,'#N_coarse=',N_coarse
      print *,'#saving:',real(N,f64)/real(N_coarse*N_local,f64)
      print *,'#new x(i1),x(i1+N_fine)=',bloc_coord(1),bloc_coord(2)
      print *,'#error for x(i1),x(i1+N_fine)=',bloc_coord(1)-a,bloc_coord(2)-b
      print *,'#i1,i1+N_fine,N_fine,N=',i1,i1+N_fine,N_fine,N

    else

      SLL_ERROR('case in compute_bloc not implemented yet')

    endif
    
    
  end subroutine compute_bloc
  
  
  !Input:   
  !  x1=bloc_coord(1),x2=bloc_coord(2)
  !  with 0<i1<i2<N i1=bloc_index(1), i2=i1+bloc_index(2)
  !  N=bloc_index(1)+bloc_index(2)+bloc_index(3)
  !
  !Output:  
  !  node_positions(1:N+1)
  !  with constraints node_positions(i1+1)=x1,node_positions(i2+1)=x2
  !  node_positions(1)=0, node_positions(N+1)=1
  subroutine compute_mesh_from_bloc(bloc_coord,bloc_index,node_positions)

    sll_int32,               intent(in)  :: bloc_index(3)
    sll_real64,              intent(in)  :: bloc_coord(2)
    sll_real64,dimension(:), intent(out) :: node_positions

    sll_int32  :: i, i1, i2, N
    sll_real64 :: dx
    
    N                     = bloc_index(1)+bloc_index(2)+bloc_index(3)
    i1                    = bloc_index(1)
    i2                    = i1+bloc_index(2)
    node_positions(1:N+1) = -1._f64
    node_positions(1)     = 0._f64
    node_positions(i1+1)  = bloc_coord(1)
    node_positions(i2+1)  = bloc_coord(2)
    node_positions(N+1)   = 1._f64
    
    !piecewise linear mapping (maybe enhanced like in complete mesh)
    dx=bloc_coord(1)/real(bloc_index(1),f64)
    do i=2,bloc_index(1)
       node_positions(i) = (real(i,f64)-1._f64)*dx
    enddo
    dx=(bloc_coord(2)-bloc_coord(1))/real(bloc_index(2),f64)
    do i=2,bloc_index(2)
      node_positions(i+i1)=bloc_coord(1)+(real(i,f64)-1._f64)*dx
    enddo
    dx=(1._f64-bloc_coord(2))/real(bloc_index(3),f64)
    do i=2,bloc_index(3)
      node_positions(i+i2)=bloc_coord(2)+(real(i,f64)-1._f64)*dx
    enddo
          
  end subroutine compute_mesh_from_bloc

  function compute_light_size( bloc_index, every ) result(light_size)

    sll_int32, intent(in) :: bloc_index
    sll_int32, intent(in) :: every

    sll_int32 :: light_size
    sll_int32 :: i
    sll_int32 :: s
        
    s = 1
    do i=1,bloc_index-1,every
      s = s+1
    enddo
    light_size = s      

  end function compute_light_size

  function compute_light_size_3_bloc( bloc_index, every ) result(light_size)

    sll_int32, intent(in) :: bloc_index(3)
    sll_int32, intent(in) :: every(3)

    sll_int32 :: light_size
    sll_int32 :: i1
    sll_int32 :: i2
    sll_int32 :: i
    sll_int32 :: s
        
    i1 = bloc_index(1)
    i2 = i1 + bloc_index(2)
    s = 0
    do i=1,bloc_index(1)-1,every(1)
      s = s+1
    enddo
    s = s+1 !for bloc_index(1)
    do i=2,bloc_index(2)-1,every(2)
      s = s+1
    enddo
    s = s+1 !for bloc_index(2)
    do i=2,bloc_index(3)-1,every(3)
      s = s+1
    enddo
    s = s+1 !for bloc_index(3)
    light_size = s      

  end function compute_light_size_3_bloc

  function compute_light_size_bloc( bloc_index, num_bloc, every ) result(light_size)

    sll_int32, dimension(:), intent(in) :: bloc_index
    sll_int32,               intent(in) :: num_bloc 
    sll_int32, dimension(:), intent(in) :: every

    sll_int32 :: light_size
    sll_int32 :: i
    sll_int32 :: j
    sll_int32 :: s
    
    if (size(bloc_index)<num_bloc) then
      print *,'#bad size for every',size(bloc_index),num_bloc
      SLL_ERROR( '#in function compute_light_size_bloc')
    endif

    if (size(every)<num_bloc) then
      print *,'#bad size for every',size(every),num_bloc
      SLL_ERROR( '#in function compute_light_size_bloc')
    endif

    s = 0
    do j=1,num_bloc
      do i=1,bloc_index(j)-1,every(j)
        s = s+1
      enddo    
      s = s+1 !for bloc_index(1)
    enddo  
    light_size = s      

  end function compute_light_size_bloc

  subroutine compute_light_mesh_3_bloc( bloc_index,           &
                                        node_positions,       &
                                        light_node_positions, &
                                        light_size,           &
                                        every )

    sll_int32,                intent(in)  :: bloc_index(3)
    sll_real64, dimension(:), intent(in)  :: node_positions
    sll_real64, dimension(:), intent(out) :: light_node_positions
    sll_int32,                intent(in)  :: light_size
    sll_int32,                intent(in)  :: every(3)

    sll_int32                             :: i,i1,i2,s
    
    if (size(light_node_positions)<light_size) then
      print *,'#bad value of light_size',light_size
      print *,'#or allocation of light_node_positions',size(light_node_positions)
      SLL_ERROR('#in compute_light_mesh')
    endif
    i1=bloc_index(1)
    i2=i1+bloc_index(2)

    s = 0
    do i=1,bloc_index(1)-1,every(1)
      s = s+1
      light_node_positions(s) = node_positions(i)
    enddo
    i=bloc_index(1)
    s = s+1
    light_node_positions(s) = node_positions(i)
    do i=2,bloc_index(2)-1,every(2)
      s = s+1
      light_node_positions(s) = node_positions(i+i1)
    enddo
    i=bloc_index(2)
    s = s+1
    light_node_positions(s) = node_positions(i+i1)
    do i=2,bloc_index(3)-1,every(3)
      s = s+1
      light_node_positions(s) = node_positions(i+i2)
    enddo
    i=bloc_index(3)
    s = s+1 !for bloc_index(3)
    light_node_positions(s) = node_positions(i+i2)

    if (s /= light_size) then
      print *,'#bad value of light_size',light_size,s
      print *,'#in compute_light_mesh'
      stop
    endif
    
          
  end subroutine compute_light_mesh_3_bloc


  subroutine compute_light_mesh_bloc( bloc_index,           &
                                      num_bloc,             &
                                      node_positions,       & 
                                      light_node_positions, &
                                      light_size,           &
                                      every )

    sll_int32,  dimension(:), intent(in)  :: bloc_index
    sll_int32                             :: num_bloc
    sll_real64, dimension(:), intent(in)  :: node_positions
    sll_real64, dimension(:), intent(out) :: light_node_positions
    sll_int32,                intent(in)  :: light_size
    sll_int32,  dimension(:), intent(in)  :: every

    sll_int32 :: i, i1, s, j

    if (size(light_node_positions)<light_size) then
      print *,'#bad value of light_size',light_size
      print *,'#or allocation of light_node_positions',size(light_node_positions)
      SLL_ERROR('#in compute_light_mesh')
    endif
    
    s = 0    
    do i=1,bloc_index(1)-1,every(1)
      s = s+1
      light_node_positions(s) = node_positions(i)
    enddo
    i=bloc_index(1)
    s = s+1
    light_node_positions(s) = node_positions(i)
    i1=0
    do j=2,num_bloc
      i1 = i1+bloc_index(j-1)
      do i=2,bloc_index(j)-1,every(j)
        s = s+1
        light_node_positions(s) = node_positions(i+i1)
      enddo
      i=bloc_index(j)
      s = s+1
      light_node_positions(s) = node_positions(i+i1)
    enddo  

    if (s /= light_size) then
      print *,'#bad value of light_size',light_size,s
      print *,'#in compute_light_mesh'
      stop
    endif
    
  end subroutine compute_light_mesh_bloc

  subroutine compute_light_mesh( bloc_index,           &
                                 node_positions,       &
                                 light_node_positions, &
                                 light_size,           &
                                 every                 )

    sll_int32,                intent(in)  :: bloc_index
    sll_real64, dimension(:), intent(in)  :: node_positions
    sll_real64, dimension(:), intent(out) :: light_node_positions
    sll_int32, intent(in)                 :: light_size
    sll_int32, intent(in)                 :: every
    sll_int32::i,s
    
    if (size(light_node_positions)<light_size) then
      print *,'#bad value of light_size',light_size
      print *,'#or allocation of light_node_positions',size(light_node_positions)
      print *,'#in compute_light_mesh'
      stop
    endif
    s = 0
    do i=1,bloc_index-1,every
      s = s+1
      light_node_positions(s) = node_positions(i)
    enddo
    i=bloc_index
    s = s+1
    light_node_positions(s) = node_positions(i)

    if (s /= light_size) then
      print *,'#bad value of light_size',light_size,s
      print *,'#in compute_light_mesh'
      stop
    endif
          
  end subroutine compute_light_mesh
  
  subroutine compute_full_to_light_2d( bloc_index_x1, &
                                       every_x1,      &
                                       num_bloc_x1,   &
                                       bloc_index_x2, &
                                       every_x2,      &
                                       num_bloc_x2,   &
                                       f_full,        &
                                       np_x1,         &
                                       np_x2,         &
                                       f_light,       &
                                       light_size_x1, &
                                       light_size_x2  )

    sll_int32,  dimension(:),   intent(in)  :: bloc_index_x1
    sll_int32,  dimension(:),   intent(in)  :: every_x1
    sll_int32,                  intent(in)  :: num_bloc_x1
    sll_int32,  dimension(:),   intent(in)  :: bloc_index_x2
    sll_int32,  dimension(:),   intent(in)  :: every_x2
    sll_int32,                  intent(in)  :: num_bloc_x2
    sll_real64, dimension(:,:), intent(in)  :: f_full
    sll_int32,                  intent(in)  :: np_x1
    sll_int32,                  intent(in)  :: np_x2
    sll_real64, dimension(:,:), intent(out) :: f_light
    sll_int32,                  intent(in)  :: light_size_x1
    sll_int32,                  intent(in)  :: light_size_x2

    if ((size(f_full,1)<np_x1).or.(size(f_full,2)<np_x2)) then
      print *,'#bad size for f_full 1',size(f_full,1)<np_x1
      print *,'f_full 2',size(f_full,1)<np_x1
      SLL_ERROR('#in compute_full_to_light_2d')
    endif

    if ((size(f_light,1)<light_size_x1).or.(size(f_light,2)<light_size_x2)) then
      print *,'#bad size for f_light 1',size(f_light,1)<light_size_x1
      print *,'f_light 2',size(f_light,2)<light_size_x2
      SLL_ERROR('#in compute_full_to_light_2d')
    endif

    if (size(bloc_index_x1)<num_bloc_x1) then
      print *,'#bad size for bloc_index_x1',bloc_index_x1,num_bloc_x1
      SLL_ERROR('#in compute_full_to_light_2d')
    endif

    if (size(bloc_index_x2)<num_bloc_x2) then
      print *,'#bad size for bloc_index_x2',bloc_index_x2,num_bloc_x2
      SLL_ERROR('#in compute_full_to_light_2d')
    endif

    if (size(every_x1)<num_bloc_x1) then
      print *,'#bad size for every_x1',every_x1,num_bloc_x1
      SLL_ERROR('#in compute_full_to_light_2d')
    endif

    if (size(every_x2)<num_bloc_x2) then
      print *,'#bad size for every_x2',every_x2,num_bloc_x2
      SLL_ERROR('#in compute_full_to_light_2d')
    endif
    
  end subroutine compute_full_to_light_2d
  
  elemental function f_equilibrium(v)

    sll_real64, intent(in) :: v
    sll_real64             :: f_equilibrium

    f_equilibrium = 1.0_f64/sqrt(2*sll_pi)*exp(-0.5_f64*v*v)

  end function f_equilibrium

  ! DESCRIPTION
  ! -----------
  ! S: the wave form at a given point in time. This wave form is 
  !    not scaled (its maximum value is 1).
  ! t: the time at which the envelope is being evaluated
  ! tflat, tL, tR, twL, twR, tstart, t0: the parameters defining the
  !    envelope, defined in the main portion of this program.
  ! turn_drive_off: 1 if the drive should be turned off after a time
  !    tflat, and 0 otherwise
  subroutine PFenvelope(S,               &
                        t,               &
                        tflat,           &
                        tL,              &
                        tR,              &
                        twL,             &
                        twR,             &
                        t0,              &
                        turn_drive_off)

    sll_real64, intent(in)  :: t
    sll_real64, intent(in)  :: tflat
    sll_real64, intent(in)  :: tL
    sll_real64, intent(in)  :: tR
    sll_real64, intent(in)  :: twL
    sll_real64, intent(in)  :: twR
    sll_real64, intent(in)  :: t0
    sll_real64, intent(out) :: S
    logical,    intent(in)  :: turn_drive_off

    sll_real64 :: epsilon

    ! The envelope function is defined such that it is zero at t0,
    ! rises to 1 smoothly, stay constant for tflat, and returns
    ! smoothly to zero.
    if (turn_drive_off) then
       epsilon = 0.5*(tanh((t0-tL)/twL) - tanh((t0-tR)/twR))
       S = 0.5*(tanh((t-tL)/twL) - tanh((t-tR)/twR)) - epsilon
       S = S / (1-epsilon)
    else
       epsilon = 0.5*(tanh((t0-tL)/twL) + 1)
       S = 0.5*(tanh((t-tL)/twL) + 1) - epsilon
       S = S / (1-epsilon)
    endif
    if (S<0) then
       S = 0.
    endif
    S = S + 0.*tflat ! for use of unused
    return

  end subroutine PFenvelope

  !---------------------------------------------------
  ! Save the mesh structure
  !---------------------------------------------------
  subroutine plot_f_cartesian( iplot,             &
                               f,                 &
                               node_positions_x1, &
                               nnodes_x1,         &
                               node_positions_x2, &
                               nnodes_x2,         &
                               array_name,        &
                               time)    

    use sll_xdmf
    use sll_hdf5_io_serial

    sll_int32,                  intent(in) :: iplot
    sll_real64, dimension(:,:), intent(in) :: f
    sll_real64, dimension(:),   intent(in) :: node_positions_x1
    sll_int32,                  intent(in) :: nnodes_x1
    sll_real64, dimension(:),   intent(in) :: node_positions_x2    
    sll_int32,                  intent(in) :: nnodes_x2
    character(len=*),           intent(in) :: array_name 
    sll_real64                             :: time

    sll_int32                               :: file_id
    sll_int32                               :: error
    sll_real64, dimension(:,:), allocatable :: x1
    sll_real64, dimension(:,:), allocatable :: x2
    sll_int32                               :: i
    sll_int32                               :: j
    character(len=4)                        :: cplot
    
    if (iplot == 1) then

      SLL_ALLOCATE(x1(nnodes_x1,nnodes_x2), error)
      SLL_ALLOCATE(x2(nnodes_x1,nnodes_x2), error)
      do j = 1,nnodes_x2
        do i = 1,nnodes_x1
          x1(i,j) = node_positions_x1(i) !x1_min+real(i-1,f32)*dx1
          x2(i,j) = node_positions_x2(j) !x2_min+real(j-1,f32)*dx2
        end do
      end do

#ifndef NOHDF5
      call sll_hdf5_file_create("cartesian_mesh-x1.h5",file_id,error)
      call sll_hdf5_write_array(file_id,x1,"/x1",error)
      call sll_hdf5_file_close(file_id, error)
      call sll_hdf5_file_create("cartesian_mesh-x2.h5",file_id,error)
      call sll_hdf5_write_array(file_id,x2,"/x2",error)
      call sll_hdf5_file_close(file_id, error)
#endif

      deallocate(x1)
      deallocate(x2)

    end if

    call int2string(iplot,cplot)
    call sll_xdmf_open(trim(array_name)//cplot//".xmf","cartesian_mesh", &
                       nnodes_x1,nnodes_x2,file_id,error)

    write(file_id,"(a,f8.3,a)") "<Time Value='",time,"'/>"

    call sll_xdmf_write_array(trim(array_name)//cplot,f,"values", &
                              error,file_id,"Node")
    call sll_xdmf_close(file_id,error)

  end subroutine plot_f_cartesian

end module sll_simulation_2d_vlasov_poisson_cartesian
