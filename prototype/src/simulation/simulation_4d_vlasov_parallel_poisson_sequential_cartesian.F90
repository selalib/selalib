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
!> @author
!> Michel Mehrenberger (mehrenbe@math.unistra.fr)
!> Nicolas Crouseilles
!> Edwin Chacon Golcher (for checking code writing/interfaces/reusability)
!> @brief 
!> Simulation class to solve 4D vlasov poisson equation
!> parallelization follows simulation_4d_qns_general.F90
!> differs from simulation_4d_vlasov_poisson_cartesian.F90
!> as Poisson solvers is not intended to be parallel for the moment
!> should include new advectors and input files with explicit names
!> implementation and test of high order splitting is today's application goal  


module sll_simulation_4d_vlasov_parallel_poisson_sequential_cartesian

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
  !for mesh
  use sll_logical_meshes
  
  use sll_gnuplot_parallel
  use sll_coordinate_transformation_2d_base_module
  use sll_module_coordinate_transformations_2d
  use sll_common_coordinate_transformations
  use sll_common_array_initializers_module
  use sll_parallel_array_initializer_module
  
  use sll_module_advection_1d_periodic

  use sll_poisson_2d_periodic
  
  use sll_fft

  use sll_simulation_base
  implicit none


  type, extends(sll_simulation_base_class) :: &
       sll_simulation_4d_vlasov_par_poisson_seq_cart

   !geometry
   type(sll_logical_mesh_1d), pointer :: mesh_x1
   type(sll_logical_mesh_1d), pointer :: mesh_x2
   type(sll_logical_mesh_1d), pointer :: mesh_x3
   type(sll_logical_mesh_1d), pointer :: mesh_x4


   !initial function
   procedure(sll_scalar_initializer_4d), nopass, pointer :: init_func
   sll_real64, dimension(:), pointer :: params
      
   !advector (should replace interpolator)
   class(sll_advection_1d_base), pointer    :: advect_x1
   class(sll_advection_1d_base), pointer    :: advect_x2
   class(sll_advection_1d_base), pointer    :: advect_x3
   class(sll_advection_1d_base), pointer    :: advect_x4
   
   !poisson solver
   !class(sll_poisson_2d_base), pointer   :: poisson
   type(poisson_2d_periodic), pointer   :: poisson

   
   !time_iterations
   sll_real64 :: dt
   sll_int32  :: num_iterations
   sll_int32  :: freq_diag
   sll_int32  :: freq_diag_time
  
   
   !parameters for splitting
   sll_real64, dimension(:), pointer :: split_step
   sll_int32 :: nb_split_step
   logical :: split_begin_T

   
              
   contains
     procedure, pass(sim) :: run => run_vp4d_cartesian
     procedure, pass(sim) :: init_from_file => init_vp4d_fake
     
  end type sll_simulation_4d_vlasov_par_poisson_seq_cart

contains
  
  
  function new_vlasov_par_poisson_seq_cart( &
    filename, &
    mesh_x1, &
    mesh_x2, &
    mesh_x3, &
    mesh_x4, &
    init_func, &
    params, &
    advect_x1, &
    advect_x2, &
    advect_x3, &
    advect_x4, &
    dt, &
    num_iterations, &
    freq_diag, &
    freq_diag_time, &
    split_step, &
    nb_split_step, &
    split_begin_T) &
    result(sim)    
    type(sll_simulation_4d_vlasov_par_poisson_seq_cart), pointer :: sim    
    character(len=*), intent(in), optional                                :: filename
    type(sll_logical_mesh_1d), pointer, optional :: mesh_x1
    type(sll_logical_mesh_1d), pointer, optional :: mesh_x2
    type(sll_logical_mesh_1d), pointer, optional :: mesh_x3
    type(sll_logical_mesh_1d), pointer, optional :: mesh_x4
    procedure(sll_scalar_initializer_4d), optional :: init_func
    sll_real64, dimension(:), pointer, optional :: params
    class(sll_advection_1d_base), pointer, optional    :: advect_x1
    class(sll_advection_1d_base), pointer, optional    :: advect_x2
    class(sll_advection_1d_base), pointer, optional    :: advect_x3
    class(sll_advection_1d_base), pointer, optional    :: advect_x4
    sll_real64, intent(in), optional :: dt
    sll_int32, intent(in), optional  :: num_iterations
    sll_int32, intent(in), optional  :: freq_diag
    sll_int32, intent(in), optional  :: freq_diag_time
    sll_real64, dimension(:), pointer, optional :: split_step
    sll_int32, intent(in), optional :: nb_split_step
    logical, intent(in), optional :: split_begin_T



    
    sll_int32 :: ierr   
    SLL_ALLOCATE(sim,ierr)
    
    sim%mesh_x1 => null()
    sim%mesh_x2 => null()
    sim%mesh_x3 => null()
    sim%mesh_x4 => null()
    sim%init_func => null()
        
    call initialize_vlasov_par_poisson_seq_cart( &
      sim, &
      filename, &
      mesh_x1, &
      mesh_x2, &
      mesh_x3, &
      mesh_x4, &
      init_func, &
      params, &
      advect_x1, &
      advect_x2, &
      advect_x3, &
      advect_x4, &
      dt, &
      num_iterations, &
      freq_diag, &
      freq_diag_time, &
      split_step, &
      nb_split_step, &
      split_begin_T)

     
  
  end function new_vlasov_par_poisson_seq_cart


  subroutine initialize_vlasov_par_poisson_seq_cart( &
    sim, &
    filename, &
    mesh_x1, &
    mesh_x2, &
    mesh_x3, &
    mesh_x4, &
    init_func, &
    params, &
    advect_x1, &
    advect_x2, &
    advect_x3, &
    advect_x4, &
    user_dt, &
    user_num_iterations, &
    user_freq_diag, &
    user_freq_diag_time, &
    split_step, &
    nb_split_step, &
    split_begin_T)
    
    class(sll_simulation_4d_vlasov_par_poisson_seq_cart), intent(inout) :: sim
    character(len=*), intent(in), optional                                :: filename
    type(sll_logical_mesh_1d), pointer, optional :: mesh_x1
    type(sll_logical_mesh_1d), pointer, optional :: mesh_x2
    type(sll_logical_mesh_1d), pointer, optional :: mesh_x3
    type(sll_logical_mesh_1d), pointer, optional :: mesh_x4
    procedure(sll_scalar_initializer_4d), optional :: init_func
    sll_real64, dimension(:), pointer, optional :: params
    class(sll_advection_1d_base), pointer, optional    :: advect_x1
    class(sll_advection_1d_base), pointer, optional    :: advect_x2
    class(sll_advection_1d_base), pointer, optional    :: advect_x3
    class(sll_advection_1d_base), pointer, optional    :: advect_x4
    sll_real64, intent(in), optional :: user_dt
    sll_int32, intent(in), optional  :: user_num_iterations
    sll_int32, intent(in), optional  :: user_freq_diag
    sll_int32, intent(in), optional  :: user_freq_diag_time
    sll_real64, dimension(:), pointer, optional :: split_step
    sll_int32, intent(in), optional :: nb_split_step
    logical, intent(in), optional :: split_begin_T
    intrinsic :: trim
    sll_int32             :: IO_stat
    sll_int32, parameter  :: input_file = 99
    character(len=256)      :: mesh_case_x1 
    character(len=256)      :: mesh_case_x2 
    character(len=256)      :: mesh_case_x3 
    character(len=256)      :: mesh_case_x4 
    sll_real64            :: x1_min
    sll_real64            :: x1_max
    sll_real64            :: x2_min
    sll_real64            :: x2_max
    sll_real64            :: x3_min
    sll_real64            :: x3_max
    sll_real64            :: x4_min
    sll_real64            :: x4_max
    sll_int32             :: num_cells_x1
    sll_int32             :: num_cells_x2
    sll_int32             :: num_cells_x3
    sll_int32             :: num_cells_x4    
    sll_int32             :: nbox_x1
    sll_int32             :: nbox_x2
    character(len=256)      :: advector_x1 
    character(len=256)      :: advector_x2 
    character(len=256)      :: advector_x3 
    character(len=256)      :: advector_x4 
    sll_int32             :: order_x1
    sll_int32             :: order_x2
    sll_int32             :: order_x3
    sll_int32             :: order_x4
    sll_real64            :: dt
    sll_int32             :: nbiter
    sll_int32             :: freqdiag
    sll_int32             :: freqdiag_time
    character(len=256)      :: split_case
    character(len=256)      :: init_func_case
    sll_real64            :: kmode_x1
    sll_real64            :: kmode_x2
    sll_real64            :: eps
    sll_int32             :: ierr 
      
    ! namelists for data input
    namelist / geometry /   &
      num_cells_x1, &
      num_cells_x2, &
      num_cells_x3, &
      num_cells_x4, &
      x1_min, &
      x1_max, &
      x2_min, &
      x2_max, &
      x3_min, &
      x3_max, &
      x4_min, &
      x4_max, &
      nbox_x1, &
      nbox_x2
    namelist / advector /   &
      advector_x1, &
      order_x1, &
      advector_x2, &
      order_x2, &
      advector_x3, &
      order_x3, &
      advector_x4, &
      order_x4

    namelist / time_iterations / &
      dt, &
      nbiter, &
      freqdiag, &
      freqdiag_time, &
      split_case

    namelist / initial_function / &
      init_func_case, &
      kmode_x1, &
      kmode_x2, &
      eps
 
 
 
   
    if(present(mesh_x1))then
      sim%mesh_x1 => mesh_x1
    endif
    if(present(mesh_x2))then
      sim%mesh_x2 => mesh_x2
    else
      sim%mesh_x2 => null()
    endif
    if(present(mesh_x3))then
      sim%mesh_x3 => mesh_x3
    else
      sim%mesh_x3 => null()
    endif
    if(present(mesh_x4))then
      sim%mesh_x4 => mesh_x4
    else
      sim%mesh_x4 => null()
    endif
    
    if(present(init_func))then
      sim%init_func => init_func
    else
      sim%init_func => null()  
    endif

    if(present(params))then
      sim%params => params
    else
      sim%params => null()  
    endif

    if(present(advect_x1))then
      sim%advect_x1 => advect_x1
    else
      sim%advect_x1 => advect_x1  
    endif

    if(present(advect_x2))then
      sim%advect_x2 => advect_x2
    else
      sim%advect_x2 => advect_x2  
    endif

    if(present(advect_x3))then
      sim%advect_x3 => advect_x3
    else
      sim%advect_x3 => advect_x3  
    endif

    if(present(advect_x4))then
      sim%advect_x4 => advect_x4
    else
      sim%advect_x4 => advect_x4  
    endif
    
    if(present(user_dt))then
      sim%dt = user_dt
    else
      sim%dt = 0._f64  
    endif

    if(present(user_num_iterations))then
      sim%num_iterations = user_num_iterations
    else
      sim%num_iterations = 0._f64  
    endif

    if(present(user_freq_diag))then
      sim%freq_diag = user_freq_diag
    else
      sim%freq_diag_time = user_freq_diag_time  
    endif

    if(present(split_step))then
      sim%split_step => split_step
    else
      sim%split_step => null()  
    endif
    
    if(present(nb_split_step))then
      sim%nb_split_step = nb_split_step
    else
      sim%nb_split_step = 0 
    endif

    if(present(split_begin_T))then
      sim%split_begin_T = split_begin_T
    else
      sim%split_begin_T = .true. 
    endif
        
 
    if(present(filename))then
 
    open(unit = input_file, file=trim(filename)//'.nml',IOStat=IO_stat)
    if( IO_stat /= 0 ) then
       print *, '#initialize_vlasov_par_poisson_seq_cart() failed to open file ', &
         trim(filename)//'.nml'
       STOP
    end if
    read(input_file, geometry) 
    read(input_file, advector)
    read(input_file, time_iterations)
    read(input_file, initial_function)
    close(input_file)
    
    select case (mesh_case_x1)
      case ("SLL_LANDAU")
        x1_max = real(nbox_x1,f64) * 2._f64 * sll_pi / kmode_x1
        sim%mesh_x1 => new_logical_mesh_1d(num_cells_x1,eta_min=x1_min, eta_max=x1_max)
      case ("LOGICAL_MESH")
        sim%mesh_x1 => new_logical_mesh_1d(num_cells_x1,eta_min=x1_min, eta_max=x1_max)  
      case default
         print*,'#mesh_case_x1', mesh_case_x1, ' not implemented'
         print*,'#in initialize_vlasov_par_poisson_seq_cart'
         stop 
    end select

    select case (mesh_case_x2)
      case ("SLL_LANDAU")
        x2_max = real(nbox_x2,f64) * 2._f64 * sll_pi / kmode_x1
        sim%mesh_x2 => new_logical_mesh_1d(num_cells_x2,eta_min=x2_min, eta_max=x2_max)
      case ("LOGICAL_MESH")
        sim%mesh_x2 => new_logical_mesh_1d(num_cells_x2,eta_min=x2_min, eta_max=x2_max)
      case default
         print*,'#mesh_case_x2', mesh_case_x2, ' not implemented'
         print*,'#in initialize_vlasov_par_poisson_seq_cart'
         stop 
    end select

    select case (mesh_case_x3)
      case ("LOGICAL_MESH")
        sim%mesh_x3 => new_logical_mesh_1d(num_cells_x3,eta_min=x3_min, eta_max=x3_max)
      case default
         print*,'#mesh_case_x3', mesh_case_x3, ' not implemented'
         print*,'#in initialize_vlasov_par_poisson_seq_cart'
         stop 
    end select

    select case (mesh_case_x4)
      case ("LOGICAL_MESH")
        sim%mesh_x4 => new_logical_mesh_1d(num_cells_x4,eta_min=x4_min, eta_max=x4_max)
      case default
         print*,'#mesh_case_x4', mesh_case_x4, ' not implemented'
         print*,'#in initialize_vlasov_par_poisson_seq_cart'
         stop 
    end select
    


    select case (advector_x1)
      case ("SLL_SPLINES") ! arbitrary order periodic splines
        sim%advect_x1 => new_periodic_1d_advector( &
          num_cells_x1, &
          x1_min, &
          x1_max, &
          SPLINE, & 
          order_x1) 
      case("SLL_LAGRANGE") ! arbitrary order Lagrange periodic interpolation
        sim%advect_x1 => new_periodic_1d_advector( &
          num_cells_x1, &
          x1_min, &
          x1_max, &
          LAGRANGE, & 
          order_x1)
       case default
         print*,'#advector in x1', advector_x1, ' not implemented'
         stop 
    end select



    select case (advector_x2)
      case ("SLL_SPLINES") ! arbitrary order periodic splines
        sim%advect_x2 => new_periodic_1d_advector( &
          num_cells_x2, &
          x2_min, &
          x2_max, &
          SPLINE, & 
          order_x2) 
      case("SLL_LAGRANGE") ! arbitrary order Lagrange periodic interpolation
        sim%advect_x2 => new_periodic_1d_advector( &
          num_cells_x2, &
          x2_min, &
          x2_max, &
          LAGRANGE, & 
          order_x2) 
       case default
         print*,'#advector in x2', advector_x2, ' not implemented'
         stop 
    end select

    select case (advector_x3)
      case ("SLL_SPLINES") ! arbitrary order periodic splines
        sim%advect_x3 => new_periodic_1d_advector( &
          num_cells_x3, &
          x3_min, &
          x3_max, &
          SPLINE, & 
          order_x3) 
      case("SLL_LAGRANGE") ! arbitrary order Lagrange periodic interpolation
        sim%advect_x3 => new_periodic_1d_advector( &
          num_cells_x3, &
          x3_min, &
          x3_max, &
          LAGRANGE, & 
          order_x3) 
       case default
         print*,'#advector in x3', advector_x3, ' not implemented'
         stop 
    end select

    select case (advector_x4)
      case ("SLL_SPLINES") ! arbitrary order periodic splines
        sim%advect_x4 => new_periodic_1d_advector( &
          num_cells_x4, &
          x4_min, &
          x4_max, &
          SPLINE, & 
          order_x4) 
      case("SLL_LAGRANGE") ! arbitrary order Lagrange periodic interpolation
        sim%advect_x4 => new_periodic_1d_advector( &
          num_cells_x4, &
          x4_min, &
          x4_max, &
          LAGRANGE, & 
          order_x4) 
       case default
         print*,'#advector in x4', advector_x4, ' not implemented'
         stop 
    end select

    sim%dt=dt
    sim%num_iterations=nbiter
    sim%freq_diag=freqdiag
    sim%freq_diag_time=freqdiag_time



    select case (split_case)    
      case ("SLL_STRANG_TVT") 
        sim%nb_split_step = 3
        SLL_ALLOCATE(sim%split_step(sim%nb_split_step),ierr)
        sim%split_begin_T = .true.
        sim%split_step(1) = 0.5_f64
        sim%split_step(2) = 1._f64
        sim%split_step(3) = sim%split_step(1)
      case ("SLL_STRANG_VTV") 
        sim%nb_split_step = 3
        SLL_ALLOCATE(sim%split_step(sim%nb_split_step),ierr)
        sim%split_begin_T = .false.
        sim%split_step(1) = 0.5_f64
        sim%split_step(2) = 1._f64
        sim%split_step(3) = sim%split_step(1)
      case ("SLL_TRIPLE_JUMP_TVT") 
        sim%nb_split_step = 7
        SLL_ALLOCATE(sim%split_step(sim%nb_split_step),ierr)
        sim%split_begin_T = .true.
        sim%split_step(1) = 0.675603595979829_f64
        sim%split_step(2) = 1.351207191959658_f64
        sim%split_step(3) = -0.17560359597982855_f64
        sim%split_step(4) = -1.702414383919315_f64
        sim%split_step(5) = sim%split_step(3)
        sim%split_step(6) = sim%split_step(2)
        sim%split_step(7) = sim%split_step(1)
      case ("SLL_TRIPLE_JUMP_VTV") 
        sim%nb_split_step = 7
        SLL_ALLOCATE(sim%split_step(sim%nb_split_step),ierr)
        sim%split_begin_T = .false.
        sim%split_step(1) = 0.675603595979829_f64
        sim%split_step(2) = 1.351207191959658_f64
        sim%split_step(3) = -0.17560359597982855_f64
        sim%split_step(4) = -1.702414383919315_f64
        sim%split_step(5) = sim%split_step(3)
        sim%split_step(6) = sim%split_step(2)
        sim%split_step(7) = sim%split_step(1)
      case ("SLL_ORDER6_VTV") 
        sim%nb_split_step = 23
        SLL_ALLOCATE(sim%split_step(sim%nb_split_step),ierr)
        sim%split_begin_T = .false.
        sim%split_step(1) = 0.0414649985182624_f64
        sim%split_step(2) = 0.123229775946271_f64
        sim%split_step(3) = 0.198128671918067_f64
        sim%split_step(4) = 0.290553797799558_f64
        sim%split_step(5) = -0.0400061921041533_f64
        sim%split_step(6) = -0.127049212625417_f64
        sim%split_step(7) = 0.0752539843015807_f64          
        sim%split_step(8) = -0.246331761062075_f64
        sim%split_step(9) = -0.0115113874206879_f64
        sim%split_step(10) = 0.357208872795928_f64
        sim%split_step(11) = 0.23666992478693111_f64
        sim%split_step(12) = 0.20477705429147008_f64
        sim%split_step(13) = sim%split_step(11)
        sim%split_step(14) = sim%split_step(10)
        sim%split_step(15) = sim%split_step(9)          
        sim%split_step(16) = sim%split_step(8)
        sim%split_step(17) = sim%split_step(7)
        sim%split_step(18) = sim%split_step(6)
        sim%split_step(19) = sim%split_step(5)
        sim%split_step(20) = sim%split_step(4)
        sim%split_step(21) = sim%split_step(3)
        sim%split_step(22) = sim%split_step(2)
        sim%split_step(23) = sim%split_step(1)          
      case ("SLL_ORDER6VP_TVT") 
        sim%nb_split_step = 9
        SLL_ALLOCATE(sim%split_step(sim%nb_split_step),ierr)
        sim%split_begin_T = .true.
        sim%split_step(1) = 0.1095115577513980413559540_f64
        sim%split_step(2) = 0.268722208204814693684441_f64&
          -2._f64*sim%dt**2*0.000805681667096178271312_f64&
          +4._f64*sim%dt**4*0.000017695766224036466792_f64
        sim%split_step(3) = 0.4451715080955340951457244_f64
        sim%split_step(4) = 0.2312777917951853063155588_f64&
          -2._f64*sim%dt**2*0.003955911930042478239763_f64&
          +4._f64*sim%dt**4*0.000052384078562246674986_f64
        sim%split_step(5) = -0.1093661316938642730033570_f64
        sim%split_step(6) = sim%split_step(4)
        sim%split_step(7) = sim%split_step(3)
        sim%split_step(8) = sim%split_step(2)
        sim%split_step(9) = sim%split_step(1)
      case ("SLL_ORDER6VP_VTV") 
        sim%nb_split_step = 9
        SLL_ALLOCATE(sim%split_step(sim%nb_split_step),ierr)
        sim%split_begin_T = .false.
        sim%split_step(1) = 0.359950808794143627485664_f64&
          -2._f64*sim%dt**2*(-0.01359558332625151635_f64)&
          +4._f64*sim%dt**4*(-8.562814848565929e-6_f64)
        sim%split_step(2) = 1.079852426382430882456991_f64
        sim%split_step(3) = -0.1437147273026540434771131_f64&
          -2._f64*sim%dt**2*(-0.00385637757897273261_f64)&
          +4._f64*sim%dt**4*(0.0004883788785819335822_f64)
        sim%split_step(4) = -0.579852426382430882456991_f64
        sim%split_step(5) = 0.567527837017020831982899_f64&
          -2._f64*sim%dt**2*(-0.03227361602037480885_f64)&
          +4._f64*sim%dt**4*0.002005141087312622342_f64
        sim%split_step(6) = sim%split_step(4)
        sim%split_step(7) = sim%split_step(3)
        sim%split_step(8) = sim%split_step(2)
        sim%split_step(9) = sim%split_step(1)
      case ("SLL_ORDER6VPnew_TVT") 
        sim%nb_split_step = 9
        SLL_ALLOCATE(sim%split_step(sim%nb_split_step),ierr)
        sim%split_begin_T = .true.
        sim%split_step(1) = 0.1095115577513980413559540_f64
        sim%split_step(2) = 0.268722208204814693684441_f64&
          -2._f64*sim%dt**2*0.000805681667096178271312_f64&
          +4._f64*sim%dt**4*(8.643923349886021963e-6_f64)&
          -8._f64*sim%dt**6*(1.4231479258353431522e-6_f64)
        sim%split_step(3) = 0.4451715080955340951457244_f64
        sim%split_step(4) = 0.2312777917951853063155588_f64&
          -2._f64*sim%dt**2*0.003955911930042478239763_f64&
          +4._f64*sim%dt**4*(0.000061435921436397119815_f64)
        sim%split_step(5) = -0.1093661316938642730033570_f64
        sim%split_step(6) = sim%split_step(4)
        sim%split_step(7) = sim%split_step(3)
        sim%split_step(8) = sim%split_step(2)
        sim%split_step(9) = sim%split_step(1)
      case ("SLL_ORDER6VPnew1_VTV") 
        sim%nb_split_step = 11
        SLL_ALLOCATE(sim%split_step(sim%nb_split_step),ierr)
        sim%split_begin_T = .false.
        sim%split_step(1) = 0.0490864609761162454914412_f64&
          -2._f64*sim%dt**2*(0.0000697287150553050840999_f64)
        sim%split_step(2) = 0.1687359505634374224481957_f64
        sim%split_step(3) = 0.2641776098889767002001462_f64&
          -2._f64*sim%dt**2*(0.000625704827430047189169_f64)&
          +4._f64*sim%dt**4*(-2.91660045768984781644e-6_f64)
        sim%split_step(4) = 0.377851589220928303880766_f64
        sim%split_step(5) = 0.1867359291349070543084126_f64&
          -2._f64*sim%dt**2*(0.00221308512404532556163_f64)&
          +4._f64*sim%dt**4*(0.0000304848026170003878868_f64)&
          -8._f64*sim%dt**6*(4.98554938787506812159e-7_f64)
        sim%split_step(6) = -0.0931750795687314526579244_f64
        sim%split_step(7) = sim%split_step(5)
        sim%split_step(8) = sim%split_step(4)
        sim%split_step(9) = sim%split_step(3)
        sim%split_step(10) = sim%split_step(2)
        sim%split_step(11) = sim%split_step(1)
      case ("SLL_ORDER6VPnew2_VTV") 
        sim%nb_split_step = 11
        SLL_ALLOCATE(sim%split_step(sim%nb_split_step),ierr)
        sim%split_begin_T = .false.
        sim%split_step(1) = 0.083335463273305120964507_f64&
          -2._f64*sim%dt**2*(-0.00015280483587048489661_f64)&
          +4._f64*sim%dt**4*( -0.0017675734111895638156_f64)&
          -8._f64*sim%dt**6*( 0.00021214072262165668039_f64)
        sim%split_step(2) = 0.72431592569108212422250_f64
        sim%split_step(3) = 0.827694857845135145869413_f64&
          -2._f64*sim%dt**2*(-0.010726848627286273332_f64)&
          +4._f64*sim%dt**4*(0.012324362982853212700_f64)
        sim%split_step(4) = -0.4493507217041624582458844_f64
        sim%split_step(5) = -0.4110303211184402668339201_f64&
          -2._f64*sim%dt**2*(0.014962337009932798678_f64)
        !  +4._f64*sim%dt**4*(-0.20213028317750837901_f64)
        sim%split_step(6) = 0.4500695920261606680467717_f64
        sim%split_step(7) = sim%split_step(5)
        sim%split_step(8) = sim%split_step(4)
        sim%split_step(9) = sim%split_step(3)
        sim%split_step(10) = sim%split_step(2)
        sim%split_step(11) = sim%split_step(1)
      case default
        print *,'#sim%split_case not defined'
        print *,'#in initialize_vlasov_par_poisson_seq_cart'
        stop       
    end select


    sim%poisson => new(&
      x1_min, &
      x1_max, &
      num_cells_x1, &
      x2_min, &
      x2_max, &
      num_cells_x2, &
      ierr)
    select case (init_func_case)
      case ("SLL_LANDAU")
        sim%init_func => sll_landau_mode_initializer_4d
        SLL_ALLOCATE(sim%params(3),ierr)
        sim%params(1) = kmode_x1
        sim%params(2) = kmode_x2
        sim%params(3) = eps
      case default
        print *,'#init_func_case not implemented'
        print *,'#in initialize_vlasov_par_poisson_seq_cart'  
        stop
    end select
    
    endif    
  end subroutine initialize_vlasov_par_poisson_seq_cart
  
  subroutine init_vp4d_fake(sim, filename)
    class(sll_simulation_4d_vlasov_par_poisson_seq_cart), intent(inout) :: sim
    character(len=*), intent(in)                                :: filename
  
    print *,'# Do not use the routine init_vp4d_fake'
    print *,'#use instead initialize_vlasov_par_poisson_seq_cart'
    stop
  
  end subroutine init_vp4d_fake
  
  subroutine run_vp4d_cartesian(sim)
    class(sll_simulation_4d_vlasov_par_poisson_seq_cart), intent(inout) :: sim
    
    if(sll_get_collective_rank(sll_world_collective)==0) then
      print *,'#not implemented for the moment!'
    endif
  end subroutine run_vp4d_cartesian    
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  

end module sll_simulation_4d_vlasov_parallel_poisson_sequential_cartesian
