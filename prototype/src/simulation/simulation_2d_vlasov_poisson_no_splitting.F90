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
!

module sll_simulation_2d_vlasov_poisson_no_splitting

#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_field_2d.h"
#include "sll_utilities.h"
#include "sll_poisson_solvers.h"
  use sll_constants
  use sll_logical_meshes  
  use sll_gnuplot_parallel
  use sll_coordinate_transformation_2d_base_module
  use sll_module_coordinate_transformations_2d
  use sll_common_coordinate_transformations
  use sll_common_array_initializers_module
  use sll_module_advection_2d_BSL
  use sll_module_characteristics_2d_explicit_euler
  use sll_module_characteristics_2d_verlet
  use sll_cubic_spline_interpolator_2d
  use sll_cubic_spline_interpolator_1d
  use sll_poisson_1d_periodic  
  use sll_fft
  use sll_simulation_base
  use sll_module_poisson_1d_periodic_solver
  use sll_module_poisson_1d_polar_solver




  implicit none

  sll_int32, parameter :: SLL_EULER = 0 
  sll_int32, parameter :: SLL_PREDICTOR_CORRECTOR = 1 
  sll_int32, parameter :: SLL_LEAP_FROG = 2 


  type, extends(sll_simulation_base_class) :: &
       sll_simulation_2d_vlasov_poisson_no_split

   !geometry
   type(sll_logical_mesh_2d), pointer :: mesh_2d
   sll_real64, dimension(:), pointer :: x1_array
   sll_real64, dimension(:), pointer :: x2_array
   sll_real64, dimension(:), pointer :: integration_weight
      
   !initial function
   sll_real64  :: kx
   sll_real64  :: eps

   !initial function
   procedure(sll_scalar_initializer_2d), nopass, pointer :: init_func
   sll_real64, dimension(:), pointer :: params
   sll_real64 :: nrj0
   
   !time_iterations
   sll_real64 :: dt
   sll_int32 :: num_iterations
   sll_int32 :: freq_diag
   sll_int32 :: freq_diag_time
   sll_int32 :: time_loop_case
   sll_int32 :: freq_leap_frog
  
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
   class(sll_advection_2d_base), pointer    :: advect_2d
   !interpolator for derivatives
   class(sll_interpolator_2d_base), pointer   :: phi_interp2d
   sll_real64 :: factor_x1
   sll_real64 :: factor_x2_rho
   sll_real64 :: factor_x2_1

   !poisson solver
   class(sll_poisson_1d_base), pointer   :: poisson
           
   contains
     procedure, pass(sim) :: run => run_vp2d_no_split
     procedure, pass(sim) :: init_from_file => init_vp2d_fake !init_vp2d_par_cart
  end type sll_simulation_2d_vlasov_poisson_no_split

  interface delete
     module procedure delete_vp2d_no_split
  end interface delete

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

  function new_vp2d_no_split( &
    filename ) &
    result(sim)    
    type(sll_simulation_2d_vlasov_poisson_no_split), pointer :: sim    
    character(len=*), intent(in), optional :: filename
    sll_int32 :: ierr   

    SLL_ALLOCATE(sim,ierr)            
    call init_vp2d_no_split( &
      sim, &
      filename)
       
  end function new_vp2d_no_split



  subroutine init_vp2d_no_split( sim, filename )
    class(sll_simulation_2d_vlasov_poisson_no_split), intent(inout) :: sim
    character(len=*), intent(in), optional :: filename

    !geometry
    character(len=256) :: mesh_case_x1
    sll_int32 :: num_cells_x1
    sll_real64 :: x1_min
    sll_real64 :: x1_max
    sll_int32 :: nbox_x1
    character(len=256) :: mesh_case_x2
    sll_int32 :: num_cells_x2
    sll_real64 :: x2_min
    sll_real64 :: x2_max
    
    !initial_function
    character(len=256) :: initial_function_case
    sll_real64 :: kmode
    sll_real64 :: eps
    sll_real64 :: alpha_gaussian
    
    !time_iterations
    sll_real64 :: dt
    sll_int32 :: number_iterations
    sll_int32 :: freq_diag
    sll_int32 :: freq_diag_time
    character(len=256) :: time_loop_case
    sll_int32 :: freq_leap_frog

    !advector
    character(len=256) :: advect2d_case 
    character(len=256) :: f_interp2d_case
    character(len=256) :: phi_interp2d_case
    character(len=256) ::  charac2d_case
    character(len=256) ::  A_interp_case
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
    sll_int32, parameter  :: input_file = 99
    type(sll_logical_mesh_1d), pointer :: mesh_x1
    type(sll_logical_mesh_1d), pointer :: mesh_x2
    sll_int32 :: Nc_x1
    sll_int32 :: Nc_x2
    class(sll_interpolator_2d_base), pointer :: f_interp2d
    class(sll_interpolator_2d_base), pointer :: phi_interp2d
    class(sll_characteristics_2d_base), pointer :: charac2d
    class(sll_interpolator_2d_base), pointer   :: A1_interp2d
    class(sll_interpolator_2d_base), pointer   :: A2_interp2d
    class(sll_interpolator_1d_base), pointer   :: A1_interp1d_x1
    class(sll_interpolator_1d_base), pointer   :: A2_interp1d_x1
    sll_int32 :: ierr
    sll_int32, parameter  :: param_out = 37, param_out_drive = 40
    sll_int32 :: i

    ! namelists for data input
    namelist /geometry/ &
      mesh_case_x1, &
      num_cells_x1, &
      x1_min, &
      x1_max, &
      nbox_x1, &
      mesh_case_x2, &
      num_cells_x2, &
      x2_min, &
      x2_max


    namelist /initial_function/ &
      initial_function_case, &
      kmode, &
      eps, &
      alpha_gaussian

    namelist /time_iterations/ &
      dt, &
      number_iterations, &
      freq_diag, &
      freq_diag_time, &
      time_loop_case, &
      freq_leap_frog

    namelist /advector/ &
      advect2d_case, &   
      f_interp2d_case, &
      phi_interp2d_case, &
      charac2d_case, &
      A_interp_case, &
      factor_x1, &
      factor_x2_rho, &
      factor_x2_1, &
      integration_case

    namelist /poisson/ &
      poisson_solver


    namelist / drive / &
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


    
    
    !! set default parameters
    
    !geometry
    mesh_case_x1 = "SLL_LANDAU_MESH"
    num_cells_x1 = 32
    x1_min = 0.0_f64
    nbox_x1 = 1
    mesh_case_x2 = "SLL_LOGICAL_MESH"
    num_cells_x2 = 64
    x2_min = -6._f64
    x2_max = 6._f64
    
        
    !initial_function
    initial_function_case = "SLL_LANDAU"
    kmode = 0.5_f64
    eps = 0.001_f64
    !initial_function_case = "SLL_BEAM"
    alpha_gaussian = 0.2_f64
    
    !time_iterations
    dt = 0.1_f64
    number_iterations = 600
    freq_diag = 100
    freq_diag_time = 1
    !time_loop_case = "SLL_EULER"
    time_loop_case = "SLL_PREDICTOR_CORRECTOR" 
    freq_leap_frog = 1e9 

    !advector
    advect2d_case = "SLL_BSL"    
    f_interp2d_case = "SLL_CUBIC_SPLINES"
    phi_interp2d_case = "SLL_CUBIC_SPLINES"
    !charac2d_case = "SLL_EULER"
    charac2d_case = "SLL_VERLET"
    A_interp_case = "SLL_CUBIC_SPLINES"
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

    if(present(filename))then
      open(unit = input_file, file=trim(filename)//'.nml',IOStat=IO_stat)
      if( IO_stat /= 0 ) then
        print *, '#init_vp2d_par_cart() failed to open file ', trim(filename)//'.nml'
        stop
      end if
      !if(sll_get_collective_rank(sll_world_collective)==0)then
        print *,'#initialization with filename:'
        print *,'#',trim(filename)//'.nml'
      !endif
      read(input_file, geometry) 
      read(input_file, initial_function)
      read(input_file, time_iterations)
      read(input_file, advector)
      read(input_file, poisson)
      read(input_file, drive)
      close(input_file)
    else
      !if(sll_get_collective_rank(sll_world_collective)==0)then
        print *,'#initialization with default parameters'
      !endif      
    endif

    Nc_x1 = num_cells_x1
    Nc_x2 = num_cells_x2


    !geometry
    
    select case (mesh_case_x1)
      case ("SLL_LANDAU_MESH")
        x1_max = real(nbox_x1,f64) * 2._f64 * sll_pi / kmode
        mesh_x1 => new_logical_mesh_1d(num_cells_x1,eta_min=x1_min, eta_max=x1_max)
        call initialize_eta1_node_1d( mesh_x1, sim%x1_array )
      case ("SLL_LOGICAL_MESH")
        mesh_x1 => new_logical_mesh_1d(num_cells_x1,eta_min=x1_min, eta_max=x1_max)  
        call initialize_eta1_node_1d( mesh_x1, sim%x1_array )
      case default
        print*,'#mesh_case_x1', mesh_case_x1, ' not implemented'
        print*,'#in init_vp2d_par_cart'
        stop 
    end select
    select case (mesh_case_x2)
      case ("SLL_LOGICAL_MESH")
        mesh_x2 => new_logical_mesh_1d(num_cells_x2,eta_min=x2_min, eta_max=x2_max)
        call initialize_eta1_node_1d( mesh_x2, sim%x2_array )
      case default
        print*,'#mesh_case_x2', mesh_case_x2, ' not implemented'
        print*,'#in init_vp2d_par_cart'
        stop 
    end select
    !sim%mesh2d => tensor_product_1d_1d( mesh_x1, mesh_x2)
    
    
    !initial function
    sim%nrj0 = 0._f64
    sim%kx = kmode
    sim%eps = eps
    select case (initial_function_case)
      case ("SLL_LANDAU")
        sim%init_func => sll_landau_initializer_2d
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
        print *,'#init_func_case not implemented'
        print *,'#in init_vp2d_par_cart'  
        stop
    end select


    !time iterations
    sim%dt=dt
    sim%num_iterations=number_iterations
    sim%freq_diag=freq_diag
    sim%freq_diag_time=freq_diag_time
    !time_loop
    select case(time_loop_case)
      case ("SLL_EULER")
        sim%time_loop_case = SLL_EULER
      case ("SLL_PREDICTOR_CORRECTOR")
        sim%time_loop_case = SLL_PREDICTOR_CORRECTOR
      case ("SLL_LEAP_FROG")
        sim%time_loop_case = SLL_LEAP_FROG
        sim%freq_leap_frog = freq_leap_frog
      case default
        print *,'#bad time_loop_case',time_loop_case
        print *,'#not implemented'
        print *,'#in initialize_guiding_center_2d_polar'
        stop
    end select

    !advector 
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
    
    sim%factor_x1 = factor_x1
    sim%factor_x2_rho = factor_x2_rho
    sim%factor_x2_1 = factor_x2_1
    
    
    
    SLL_ALLOCATE(sim%integration_weight(num_cells_x2+1),ierr)
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
        print *,'#integration_case not implemented'
        print *,'#in init_vp2d_par_cart'  
        stop      
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





    


    
     
    
            
      !to be compatible with VPpostprocessing_drive_KEEN.f
      if(sim%driven) then     
      open(unit = param_out, file = 'parameters.dat')
      
      
        write(param_out, *) real(number_iterations,f64)*dt !Tmax
        write(param_out, *) dt                  !dt
        write(param_out, *) number_iterations              !Nt
        write(param_out, *) kmode               !k0
        write(param_out, *) sim%omegadr             !omega0
        write(param_out, *) x1_max-x1_min           !L
        write(param_out, *) nbox_x1                !mbox
        write(param_out, *) sim%Edrmax              !Edrmax
        write(param_out, *) freq_diag_time       !nsave
        write(param_out, *) freq_diag            !nsavef1
        write(param_out, *) freq_diag            !nsavef2
        write(param_out, *) sim%tR+sim%tflat            !Tsetup
        write(param_out, *) x2_max                !vxmax
        write(param_out, *) x2_min                !vxmin
        write(param_out, *) num_cells_x1                 !N
        write(param_out, *) num_cells_x2                 !Nv
        ! added ES
        write(param_out, *) sim%tL                  !tL
        write(param_out, *) sim%tR                  !tR
        write(param_out, *) sim%twL                 !twL
        write(param_out, *) sim%twR                 !twR
        write(param_out, *) sim%tflat               !tflat
      
      
      close(param_out)
      endif  
      

  end subroutine init_vp2d_no_split


  subroutine init_vp2d_fake(sim, filename)
    class(sll_simulation_2d_vlasov_poisson_no_split), intent(inout) :: sim
    character(len=*), intent(in)                                :: filename
  
    print *,'# Do not use the routine init_vp2d_fake'
    print *,'#use instead init_vp2d_par_cart'
    stop
  
  end subroutine init_vp2d_fake



  subroutine run_vp2d_no_split(sim)
    class(sll_simulation_2d_vlasov_poisson_no_split), intent(inout) :: sim


    sll_real64,dimension(:,:),pointer :: f_x1,f_x2,f_x1_init
    
    sll_real64,dimension(:,:),pointer :: f_store
    
    sll_real64,dimension(:),pointer :: rho,efield,e_app,rho_loc
    sll_real64, dimension(:), allocatable :: rho_split
    sll_real64, dimension(:), allocatable :: rho_full
    
    sll_int32 :: rhotot_id
    sll_int32 :: efield_id     
    sll_int32 :: adr_id
    sll_int32 :: Edr_id
    sll_int32 :: deltaf_id
    sll_int32 :: t_id
    sll_int32 :: th_diag_id
    sll_real64, dimension(:), pointer     :: f1d
    sll_int32 :: np_x1,np_x2
    sll_int32 :: nproc_x1,nproc_x2
    sll_int32 :: global_indices(2)
    sll_int32 :: ierr
    sll_int32 :: local_size_x1,local_size_x2
    type(poisson_1d_periodic)  :: poisson_1d
    sll_real64 :: adr
    sll_real64::alpha
    sll_real64 ::tmp_loc(5),tmp(5)
    sll_int32  ::i,istep,ig,k
    
    sll_real64  ::   time, mass, momentum, l1norm, l2norm
    sll_real64  ::   kinetic_energy,potential_energy

    sll_real64, dimension(:), allocatable :: x2_array
    sll_real64, dimension(:), allocatable :: x2_array_unit
    sll_real64, dimension(:), allocatable :: x2_array_middle
    sll_real64, dimension(:), allocatable :: x1_array
    character(len=4)           :: fin   
    sll_int32                  :: file_id
    
    type(sll_fft_plan), pointer         :: pfwd
    sll_real64, dimension(:), allocatable :: buf_fft
    sll_comp64,dimension(:),allocatable :: rho_mode

    sll_real64,dimension(:,:), pointer :: A1 !advection fields
    sll_real64,dimension(:,:), pointer :: A2


    
    sll_int32 :: nb_mode = 5
    sll_int32 :: i1
    sll_int32 :: i2
    sll_real64 :: t_step
    sll_real64 :: x1
    sll_real64 :: x2
     
    
    logical :: split_T
    sll_int32 ::conservative_case
    
    
    ! for parallelization (output of distribution function in one single file)
    sll_real64,dimension(:,:),pointer :: f_visu 
    sll_int32 :: iplot
    
    !for temporary poisson
    !sll_int32 :: N_buf_poisson
    !sll_real64, dimension(:), allocatable :: buf_poisson



    np_x1 = sim%mesh_2d%num_cells1+1
    np_x2 = sim%mesh_2d%num_cells2+1

    SLL_ALLOCATE(f_visu(np_x1,np_x2),ierr)
 

        
      SLL_ALLOCATE(buf_fft(np_x1-1),ierr)
      pfwd => fft_new_plan(np_x1-1,buf_fft,buf_fft,FFT_FORWARD,FFT_NORMALIZE)
      SLL_ALLOCATE(rho_mode(0:nb_mode),ierr)      


    
    SLL_ALLOCATE(f_x2(np_x1,np_x2),ierr)

    SLL_ALLOCATE(f_x1(np_x1,np_x2),ierr)    
    SLL_ALLOCATE(f_store(np_x1,np_x2),ierr)    
    SLL_ALLOCATE(f_x1_init(np_x1,np_x2),ierr)    



    
    !allocation of 1d arrays
    SLL_ALLOCATE(rho(np_x1),ierr)
    SLL_ALLOCATE(rho_loc(np_x1),ierr)
    SLL_ALLOCATE(efield(np_x1),ierr)
    SLL_ALLOCATE(e_app(np_x1),ierr)
    SLL_ALLOCATE(f1d(max(np_x1,np_x2)),ierr)
    !SLL_ALLOCATE(x2_array(np_x2),ierr)
    !SLL_ALLOCATE(x1_array(np_x1),ierr)
    SLL_ALLOCATE(x2_array_unit(np_x2),ierr)
    SLL_ALLOCATE(x2_array_middle(np_x2),ierr)

    SLL_ALLOCATE(A1(np_x1,np_x2),ierr)
    SLL_ALLOCATE(A2(np_x1,np_x2),ierr)




    x2_array_unit(1:np_x2) = &
      (sim%x2_array(1:np_x2)-sim%x2_array(1))/(sim%x2_array(np_x2)-sim%x2_array(1))
    do i = 1, np_x2-1
       x2_array_middle(i) = 0.5_f64*(sim%x2_array(i)+sim%x2_array(i+1))
    end do
    x2_array_middle(np_x2) = x2_array_middle(1)+sim%x2_array(np_x2)-sim%x2_array(1)
    
        

    if(sim%driven)then
      call sll_binary_file_create("x.bdat", file_id, ierr)
      call sll_binary_write_array_1d(file_id,sim%x1_array(1:np_x1-1),ierr)
      call sll_binary_file_close(file_id,ierr)                    
      call sll_binary_file_create("v.bdat", file_id, ierr)
      !call sll_binary_write_array_1d(file_id,x2_array(1:np_x2-1),ierr)
      call sll_binary_write_array_1d(file_id,sim%x2_array(1:np_x2-1),ierr)
      call sll_binary_file_close(file_id,ierr)                                             
    endif


    
    !initialize distribution function     

    !initialisation of distribution function
    do i2=1,np_x2
      x2=sim%x2_array(i2)
      do i1=1,np_x1
        x1=sim%x1_array(i1)
        f_x1(i1,i2) =  sim%init_func(x1,x2,sim%params)
      end do
    end do


    do i2=1,np_x2
      x2=sim%x2_array(i2)
      do i1=1,np_x1
        x1=sim%x1_array(i1)
        f_x1_init(i1,i2) =  sll_landau_initializer_2d(x1,x2,(/ sim%kx,0._f64 /))
      end do
    end do




      call sll_binary_file_create('f0.bdat', file_id, ierr)
      call sll_binary_write_array_2d(file_id,f_x1(1:np_x1-1,1:np_x2-1),ierr)
      call sll_binary_file_close(file_id,ierr)
#ifndef NOHDF5
      iplot = 1
      call plot_f_cartesian(iplot,f_x1,sim%mesh_2d)
      iplot = iplot+1  
#endif
      print *,'#maxf',maxval(f_x1), minval(f_x1) 
    !computation of electric field
    rho = 0._f64
    do i=1,np_x1
      rho(i)=rho(i)+sum(f_x1(i,1:np_x2)*sim%integration_weight(1:np_x2))
    end do
    rho = sim%factor_x2_1*1._f64-sim%factor_x2_rho*rho        
    
    call sim%poisson%compute_E_from_rho( efield, rho )
        
    ! Ponderomotive force at initial time. We use a sine wave
    ! with parameters k_dr and omega_dr.
    istep = 0
    e_app = 0._f64
    if (sim%driven) then
      call PFenvelope(adr, istep*sim%dt, sim%tflat, sim%tL, sim%tR, sim%twL, sim%twR, &
          sim%t0, sim%turn_drive_off)
      do i = 1, np_x1
        e_app(i) = sim%Edrmax*adr*sim%kx&
          *sin(sim%kx*real(i-1,f64)*sim%mesh_2d%delta_eta1&
          -sim%omegadr*real(istep,f64)*sim%dt)
      enddo
    endif

    ! write initial fields
      call sll_ascii_file_create('thdiag.dat', th_diag_id, ierr)
      call sll_binary_file_create('deltaf.bdat', deltaf_id, ierr)
      if(sim%driven)then
        call sll_binary_file_create('rhotot.bdat', rhotot_id, ierr)
        call sll_binary_file_create('efield.bdat', efield_id, ierr)
        call sll_binary_file_create('adr.bdat', adr_id, ierr)
        call sll_binary_file_create('Edr.bdat', Edr_id, ierr)
        call sll_binary_file_create('t.bdat', t_id, ierr)
        call sll_binary_write_array_1d(efield_id,efield(1:np_x1-1),ierr)
        call sll_binary_write_array_1d(rhotot_id,rho(1:np_x1-1),ierr)
        call sll_binary_write_array_1d(Edr_id,e_app(1:np_x1-1),ierr)
        call sll_binary_write_array_0d(adr_id,adr,ierr)
        call sll_binary_write_array_0d(t_id,real(istep,f64)*sim%dt,ierr)
      endif                    
    
    
    call sll_binary_write_array_2d( &
        deltaf_id, &
        f_x1(1:np_x1-1,1:np_x2-1)-f_x1_init(1:np_x1-1,1:np_x2-1), &
        ierr )




    

    do istep = 1, sim%num_iterations+1
      if (mod(istep-1,sim%freq_diag)==0) then
        print *,'#step=',istep-1,real(istep-1,f64)*sim%dt
      endif  

      f_x2 = f_x1


         !computation of electric field
          rho = 0._f64
          do i=1,np_x1
            rho(i)=rho(i)+sum(f_x1(i,1:np_x2)*sim%integration_weight(1:np_x2))
          end do
          rho = sim%factor_x2_1*1._f64-sim%factor_x2_rho*rho

          call sim%poisson%compute_E_from_rho( efield, rho )
          
          
          t_step = real(istep-1,f64)
          if (sim%driven) then
            call PFenvelope(adr, t_step*sim%dt, sim%tflat, sim%tL, sim%tR, sim%twL, sim%twR, &
              sim%t0, sim%turn_drive_off)
            do i = 1, np_x1
              e_app(i) = sim%Edrmax*adr*sim%kx&
              *sin(sim%kx*real(i-1,f64)*sim%mesh_2d%delta_eta1&
              -sim%omegadr*t_step*sim%dt)
            enddo
          endif

      
      
      !f_old = f_x2
      !f = f_x1

          do i2=1,np_x2
            do i1=1,np_x1
              A1(i1,i2) = sim%factor_x1*sim%x2_array(i2)
              A2(i1,i2) = -efield(i1)-e_app(i1)
            enddo
          enddo
      
      !print *,'#max A1',istep,maxval(A1(:,:)),minval(A1(:,:))
      !print *,'#max A2',istep,maxval(A2(:,:)),minval(A2(:,:))
      

      select case (sim%time_loop_case)
        case (SLL_EULER)
          call sim%advect_2d%advect_2d(A1, A2, sim%dt, f_x2, f_x1)
        case (SLL_PREDICTOR_CORRECTOR)
          call sim%advect_2d%advect_2d(A1, A2, 0.5_f64*sim%dt, f_x2, f_x1)


         !computation of electric field
          rho = 0._f64
          do i=1,np_x1
            rho(i)=rho(i)+sum(f_x1(i,1:np_x2)*sim%integration_weight(1:np_x2))
          end do
          rho = sim%factor_x2_1*1._f64-sim%factor_x2_rho*rho

          call sim%poisson%compute_E_from_rho( efield, rho )
          
          
          t_step = real(istep-1,f64) + 0.5_f64
          if (sim%driven) then
            call PFenvelope(adr, t_step*sim%dt, sim%tflat, sim%tL, sim%tR, sim%twL, sim%twR, &
              sim%t0, sim%turn_drive_off)
            do i = 1, np_x1
              e_app(i) = sim%Edrmax*adr*sim%kx&
              *sin(sim%kx*real(i-1,f64)*sim%mesh_2d%delta_eta1&
              -sim%omegadr*t_step*sim%dt)
            enddo
          endif

          do i2=1,np_x2
            do i1=1,np_x1
              A1(i1,i2) = sim%factor_x1*sim%x2_array(i2)
              A2(i1,i2) = -efield(i1)-e_app(i1)
            enddo
          enddo




          call sim%advect_2d%advect_2d(A1, A2, sim%dt, f_x2, f_x1)
        case (SLL_LEAP_FROG)
          if(istep==1)then
            call sim%advect_2d%advect_2d(A1, A2, sim%dt, f_x2, f_x1)
            f_store = f_x2
          else
            call sim%advect_2d%advect_2d(A1, A2, 2*sim%dt, f_store, f_x1)            
            if(mod(istep,sim%freq_leap_frog)==0) then
              f_store = 0.5_f64*(f_store+f_x1)
            else
              f_store = f_x2   
            endif                        
          endif              
        case default  
          print *,'#bad time_loop_case',sim%time_loop_case
          print *,'#not implemented'
          print *,'#in run_vp2d_no_split'
          print *,'#available options are:'
          print *,'#SLL_EULER=',SLL_EULER
          print *,'#SLL_PREDICTOR_CORRECTOR=',SLL_PREDICTOR_CORRECTOR
          
      end select






        
      

     !!DIAGNOSTICS
     if (mod(istep,sim%freq_diag_time)==0) then
        time = real(istep,f64)*sim%dt
        mass = 0._f64
        momentum = 0._f64
        l1norm = 0._f64
        l2norm = 0._f64
        kinetic_energy = 0._f64
        potential_energy = 0._f64
        tmp = 0._f64        
        do i = 1, np_x1-1        
          tmp(1)= tmp(1)+sum(f_x1(i,1:np_x2))
          tmp(2)= tmp(2)+sum(abs(f_x1(i,1:np_x2)))
          tmp(3)= tmp(3)+sum((f_x1(i,1:np_x2))**2)
          tmp(4)= tmp(4) +sum(f_x1(i,1:np_x2) &
            *sim%x2_array(1:np_x2))          
          tmp(5)= tmp(5)+sum(f_x1(i,1:np_x2) &
            *sim%x2_array(1:np_x2)**2)          
        end do
        
        mass = tmp(1) * sim%mesh_2d%delta_eta1 * sim%mesh_2d%delta_eta2
        l1norm = tmp(2)  * sim%mesh_2d%delta_eta1 * sim%mesh_2d%delta_eta2
        l2norm = tmp(3)  * sim%mesh_2d%delta_eta1 * sim%mesh_2d%delta_eta2
        momentum = tmp(4) * sim%mesh_2d%delta_eta1 * sim%mesh_2d%delta_eta2
        kinetic_energy = 0.5_f64 *tmp(5) * sim%mesh_2d%delta_eta1 * sim%mesh_2d%delta_eta2
        potential_energy = 0._f64
        do i=1, np_x1-1
          potential_energy = potential_energy+(efield(i)+e_app(i))**2
        enddo
        potential_energy = 0.5_f64*potential_energy* sim%mesh_2d%delta_eta1
          buf_fft = rho(1:np_x1-1)
          call fft_apply_plan(pfwd,buf_fft,buf_fft)
          do k=0,nb_mode
            rho_mode(k)=fft_get_mode(pfwd,buf_fft,k)
          enddo  
          write(th_diag_id,'(f12.5,13g20.12)') &
            time, &
            mass, &
            l1norm, &
            momentum, &
            l2norm, &
            kinetic_energy, &
            potential_energy, &
            kinetic_energy + potential_energy, &
            abs(rho_mode(0)), &
            abs(rho_mode(1)), &
            abs(rho_mode(2)), &
            abs(rho_mode(3)), &
            abs(rho_mode(4)), &
            abs(rho_mode(5))
          if(sim%driven)then
            call sll_binary_write_array_1d(efield_id,efield(1:np_x1-1),ierr)
            call sll_binary_write_array_1d(rhotot_id,rho(1:np_x1-1),ierr)
            call sll_binary_write_array_1d(Edr_id,e_app(1:np_x1-1),ierr)
            call sll_binary_write_array_0d(adr_id,adr,ierr)
            call sll_binary_write_array_0d(t_id,real(istep,f64)*sim%dt,ierr)
          endif   
          
        if (mod(istep,sim%freq_diag)==0) then          
          !we substract f0
          call sll_binary_write_array_2d( &
            deltaf_id, &
            f_x1(1:np_x1-1,1:np_x2-1)-f_x1_init(1:np_x1-1,1:np_x2-1), &
            ierr)  
          !we store f for visu
#ifndef NOHDF5
            call plot_f_cartesian( &
              iplot, &
              f_x1, &
              sim%mesh_2d)
            iplot = iplot+1  
#endif
                    
        endif
          

     end if


    
    
    enddo
    
    ! close files
      call sll_ascii_file_close(th_diag_id,ierr) 
      call sll_binary_file_close(deltaf_id,ierr) 
      if(sim%driven)then
        call sll_binary_file_close(efield_id,ierr)
        call sll_binary_file_close(rhotot_id,ierr)
        call sll_binary_file_close(Edr_id,ierr)
        call sll_binary_file_close(adr_id,ierr)
        call sll_binary_file_close(t_id,ierr)
      endif   


    
    
    
  end subroutine run_vp2d_no_split

  subroutine delete_vp2d_no_split( sim )
    class(sll_simulation_2d_vlasov_poisson_no_split) :: sim
    sll_int32 :: ierr
  end subroutine delete_vp2d_no_split










  elemental function f_equilibrium(v)
    sll_real64, intent(in) :: v
    sll_real64 :: f_equilibrium

    f_equilibrium = 1.0_f64/sqrt(2*sll_pi)*exp(-0.5_f64*v*v)
  end function f_equilibrium

  subroutine PFenvelope(S, t, tflat, tL, tR, twL, twR, t0, &
       turn_drive_off)

    ! DESCRIPTION
    ! -----------
    ! S: the wave form at a given point in time. This wave form is 
    !    not scaled (its maximum value is 1).
    ! t: the time at which the envelope is being evaluated
    ! tflat, tL, tR, twL, twR, tstart, t0: the parameters defining the
    !    envelope, defined in the main portion of this program.
    ! turn_drive_off: 1 if the drive should be turned off after a time
    !    tflat, and 0 otherwise

    sll_real64, intent(in) :: t, tflat, tL, tR, twL, twR, t0
    sll_real64, intent(out) :: S
    logical, intent(in) :: turn_drive_off
    ! local variables
    sll_int32 :: i 
    sll_real64 :: epsilon

    ! The envelope function is defined such that it is zero at t0,
    ! rises to 1 smoothly, stay constant for tflat, and returns
    ! smoothly to zero.
    if(turn_drive_off) then
       epsilon = 0.5*(tanh((t0-tL)/twL) - tanh((t0-tR)/twR))
       S = 0.5*(tanh((t-tL)/twL) - tanh((t-tR)/twR)) - epsilon
       S = S / (1-epsilon)
    else
       epsilon = 0.5*(tanh((t0-tL)/twL) + 1)
       S = 0.5*(tanh((t-tL)/twL) + 1) - epsilon
       S = S / (1-epsilon)
    endif
    if(S<0) then
       S = 0.
    endif
    return
  end subroutine PFenvelope


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



end module sll_simulation_2d_vlasov_poisson_no_splitting
