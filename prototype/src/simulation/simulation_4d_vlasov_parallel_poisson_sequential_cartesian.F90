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
! Vlasov-Poisson simulation in 2Dx2D
! Poisson is sequential and Vlasov parallel
!
! contact: Michel Mehrenberger (mehrenbe@math.unistra.fr)
!
! current investigations:
!   High order splitting in time


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
  use sll_reduction_module
  

  use sll_simulation_base
  use sll_time_splitting_coeff_module
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
   sll_real64 :: nrj0   

   !time_iterations
   sll_real64 :: dt
   sll_int32  :: num_iterations
   sll_int32  :: freq_diag
   sll_int32  :: freq_diag_time
   type(splitting_coeff), pointer :: split
   
   !advector
   class(sll_advection_1d_base), pointer    :: advect_x1
   class(sll_advection_1d_base), pointer    :: advect_x2
   class(sll_advection_1d_base), pointer    :: advect_x3
   class(sll_advection_1d_base), pointer    :: advect_x4
   
   !poisson solver
   !class(sll_poisson_2d_base), pointer   :: poisson
   type(poisson_2d_periodic), pointer   :: poisson
                 
   contains
     procedure, pass(sim) :: run => run_vp4d_cartesian
     procedure, pass(sim) :: init_from_file => init_vp4d_fake
     
  end type sll_simulation_4d_vlasov_par_poisson_seq_cart

contains
  
  
  function new_vlasov_par_poisson_seq_cart( &
    filename ) &
    result(sim)    
    type(sll_simulation_4d_vlasov_par_poisson_seq_cart), pointer :: sim    
    character(len=*), intent(in), optional                                :: filename
    sll_int32 :: ierr   

    SLL_ALLOCATE(sim,ierr)            
    call initialize_vlasov_par_poisson_seq_cart( &
      sim, &
      filename)
       
  end function new_vlasov_par_poisson_seq_cart


  subroutine initialize_vlasov_par_poisson_seq_cart( &
    sim, &
    filename)
        
    class(sll_simulation_4d_vlasov_par_poisson_seq_cart), intent(inout) :: sim
    character(len=*), intent(in), optional                                :: filename
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
    sll_int32             :: number_iterations
    sll_int32             :: freq_diag
    sll_int32             :: freq_diag_time
    character(len=256)      :: split_case
    character(len=256)      :: initial_function_case
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
      nbox_x2, &
      mesh_case_x1, &
      mesh_case_x2, &
      mesh_case_x3, &
      mesh_case_x4

    namelist / initial_function / &
      initial_function_case, &
      kmode_x1, &
      kmode_x2, &
      eps

    namelist / time_iterations / &
      dt, &
      number_iterations, &
      freq_diag, &
      freq_diag_time, &
      split_case

      
    namelist / advector /   &
      advector_x1, &
      order_x1, &
      advector_x2, &
      order_x2, &
      advector_x3, &
      order_x3, &
      advector_x4, &
      order_x4


 
    !!set default parameters

    !geometry
    mesh_case_x1="SLL_LANDAU_MESH"
    num_cells_x1 = 16
    x1_min = 0.0_f64
    nbox_x1 = 1
    mesh_case_x2="SLL_LANDAU_MESH"
    num_cells_x2 = 16
    x2_min = 0.0_f64
    nbox_x2 = 1
    mesh_case_x3="SLL_LOGICAL_MESH"   
    num_cells_x3 = 16
    x3_min = -6._f64
    x3_max = 6._f64
    mesh_case_x4="SLL_LOGICAL_MESH"   
    num_cells_x4 = 16
    x4_min = -6._f64
    x4_max = 6._f64

    !initial_function
    initial_function_case="SLL_LANDAU"
    kmode_x1 = 0.5_f64
    kmode_x2 = 0.5_f64
    eps = 1e-3_f64


    !time_iterations
    dt = 2._f64
    number_iterations = 5
    freq_diag = 100
    freq_diag_time = 1
    !split_case = "SLL_STRANG_VTV" 
    !split_case = "SLL_STRANG_TVT" 
    split_case = "SLL_ORDER6VPnew1_VTV" 
    !split_case = "SLL_ORDER6VPnew2_VTV" 
    !split_case = "SLL_ORDER6_VTV"
    !split_case = "SLL_LIE_TV"



    !advector
    advector_x1 = "SLL_LAGRANGE"
    order_x1 = 4
    advector_x2 = "SLL_LAGRANGE"
    order_x2 = 4
    advector_x3 = "SLL_LAGRANGE"
    order_x3 = 4
    advector_x4 = "SLL_LAGRANGE"
    order_x4 = 4
 
    if(present(filename))then
 
      open(unit = input_file, file=trim(filename)//'.nml',IOStat=IO_stat)
      if( IO_stat /= 0 ) then
        print *, '#initialize_vlasov_par_poisson_seq_cart() failed to open file ', &
          trim(filename)//'.nml'
        stop
      end if
      if(sll_get_collective_rank(sll_world_collective)==0)then
        print *,'#initialization with filename:'
        print *,'#',trim(filename)//'.nml'
      endif      
      read(input_file, geometry) 
      read(input_file, initial_function)
      read(input_file, time_iterations)
      read(input_file, advector)
      close(input_file)
    else
      if(sll_get_collective_rank(sll_world_collective)==0)then
        print *,'#initialization with default parameters'
      endif      
    endif   
    !geometry
    select case (mesh_case_x1)
      case ("SLL_LANDAU_MESH")
        x1_max = real(nbox_x1,f64) * 2._f64 * sll_pi / kmode_x1
        sim%mesh_x1 => new_logical_mesh_1d(num_cells_x1,eta_min=x1_min, eta_max=x1_max)
      case ("SLL_LOGICAL_MESH")
        sim%mesh_x1 => new_logical_mesh_1d(num_cells_x1,eta_min=x1_min, eta_max=x1_max)  
      case default
        print*,'#mesh_case_x1', mesh_case_x1, ' not implemented'
        print*,'#in initialize_vlasov_par_poisson_seq_cart'
        stop 
    end select
    select case (mesh_case_x2)
      case ("SLL_LANDAU_MESH")
        x2_max = real(nbox_x2,f64) * 2._f64 * sll_pi / kmode_x1
        sim%mesh_x2 => new_logical_mesh_1d(num_cells_x2,eta_min=x2_min, eta_max=x2_max)
      case ("SLL_LOGICAL_MESH")
        sim%mesh_x2 => new_logical_mesh_1d(num_cells_x2,eta_min=x2_min, eta_max=x2_max)
      case default
        print*,'#mesh_case_x2', mesh_case_x2, ' not implemented'
        print*,'#in initialize_vlasov_par_poisson_seq_cart'
        stop 
    end select
    select case (mesh_case_x3)
      case ("SLL_LOGICAL_MESH")
        sim%mesh_x3 => new_logical_mesh_1d(num_cells_x3,eta_min=x3_min, eta_max=x3_max)
      case default
        print*,'#mesh_case_x3', mesh_case_x3, ' not implemented'
        print*,'#in initialize_vlasov_par_poisson_seq_cart'
        stop 
    end select
    select case (mesh_case_x4)
      case ("SLL_LOGICAL_MESH")
        sim%mesh_x4 => new_logical_mesh_1d(num_cells_x4,eta_min=x4_min, eta_max=x4_max)
      case default
        print*,'#mesh_case_x4', mesh_case_x4, ' not implemented'
        print*,'#in initialize_vlasov_par_poisson_seq_cart'
        stop 
    end select

    !initial function
    select case (initial_function_case)
      case ("SLL_LANDAU")
        sim%init_func => sll_landau_mode_initializer_4d
        SLL_ALLOCATE(sim%params(3),ierr)
        sim%params(1) = kmode_x1
        sim%params(2) = kmode_x2
        sim%params(3) = eps
        sim%nrj0 = (0.5_f64*eps*sll_pi)**2/(kmode_x1*kmode_x2) &
          *(1._f64/kmode_x1**2+1._f64/kmode_x2**2) 
      case default
        print *,'#init_func_case not implemented'
        print *,'#in initialize_vlasov_par_poisson_seq_cart'  
        stop
    end select


    !time iterations
    sim%dt=dt
    sim%num_iterations=number_iterations
    sim%freq_diag=freq_diag
    sim%freq_diag_time=freq_diag_time
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
      
    !poisson: for the moment no choice
    sim%poisson => new(&
      x1_min, &
      x1_max, &
      num_cells_x1, &
      x2_min, &
      x2_max, &
      num_cells_x2, &
      ierr)
        
  end subroutine initialize_vlasov_par_poisson_seq_cart
  
  subroutine init_vp4d_fake(sim, filename)
    class(sll_simulation_4d_vlasov_par_poisson_seq_cart), intent(inout) :: sim
    character(len=*), intent(in)                                :: filename
  
    print *,'# Do not use the routine init_vp4d_fake'
    print *,'#use instead initialize_vlasov_par_poisson_seq_cart'
    print *,sim%dt
    print *,filename
    stop
  
  end subroutine init_vp4d_fake
  
  subroutine run_vp4d_cartesian(sim)
    class(sll_simulation_4d_vlasov_par_poisson_seq_cart), intent(inout) :: sim
    sll_int32 :: world_size
    sll_int32 :: my_rank
    sll_int32 :: ierr
    sll_int32, dimension(:), allocatable :: recv_sz
    sll_int32, dimension(:), allocatable :: disps
    sll_int32 :: send_size
    sll_int32 :: nproc_x1
    sll_int32 :: nproc_x2
    sll_int32 :: nproc_x3
    sll_int32 :: nproc_x4
    sll_real64, dimension(:,:,:,:), allocatable :: f_seq_x3x4
    sll_real64, dimension(:,:,:,:), allocatable :: f_seq_x1x2
    sll_real64, dimension(:,:), allocatable :: rho_split
    sll_real64, dimension(:,:), allocatable :: rho_full
    sll_real64, dimension(:,:), allocatable :: intfdx_split
    sll_real64, dimension(:,:), allocatable :: intfdx_full
    sll_real64, dimension(:,:), allocatable :: E_x1
    sll_real64, dimension(:,:), allocatable :: E_x2
    sll_real64, dimension(:), allocatable :: send_buf_x1x2
    sll_real64, dimension(:), allocatable :: recv_buf_x1x2
    sll_real64, dimension(:), allocatable :: send_buf_x3x4
    sll_real64, dimension(:), allocatable :: recv_buf_x3x4
    sll_int32 :: loc_sz_x1
    sll_int32 :: loc_sz_x2
    sll_int32 :: loc_sz_x3
    sll_int32 :: loc_sz_x4
    sll_int32 :: nc_x1
    sll_int32 :: nc_x2
    sll_int32 :: nc_x3
    sll_int32 :: nc_x4
    sll_int32 :: itemp
    type(layout_4D), pointer :: sequential_x1x2
    type(layout_4D), pointer :: sequential_x3x4
    type(layout_2D), pointer :: layout2d_par_x1x2
    type(layout_2D), pointer :: layout2d_par_x3x4
    sll_int32 :: power2
    sll_real64 :: delta1
    sll_real64 :: delta2
    sll_real64 :: delta3
    sll_real64 :: delta4
    sll_real64 :: nrj 
    sll_real64 :: ekin
    type(remap_plan_4D_real64), pointer :: seqx1x2_to_seqx3x4
    type(remap_plan_4D_real64), pointer :: seqx3x4_to_seqx1x2
    sll_int32 :: istep
    logical :: split_T
    sll_int32 :: split_istep
    sll_real64, dimension(:), allocatable     :: f1d
    sll_int32 :: nc_max
    sll_real64 :: alpha
    sll_int32 :: i1
    sll_int32 :: i2
    sll_int32 :: i3
    sll_int32 :: i4
    sll_real64 :: time
    sll_int32 :: th_diag_id
    sll_int32 :: intfdx_id
    sll_int32 :: rho_id
    sll_int32 :: E_x1_id
    sll_int32 :: E_x2_id
    sll_int32 :: global_indices(4)
    sll_int32 :: ig
    sll_real64 :: tmp
    sll_real64 :: ekin0
    sll_real64 :: nrj0
    
    world_size = sll_get_collective_size(sll_world_collective)
    my_rank    = sll_get_collective_rank(sll_world_collective)
    
    
    nc_x1 = sim%mesh_x1%num_cells
    nc_x2 = sim%mesh_x2%num_cells
    nc_x3 = sim%mesh_x3%num_cells
    nc_x4 = sim%mesh_x4%num_cells
    
    nc_max = maxval((/nc_x1,nc_x2,nc_x3,nc_x4/))
    
    SLL_ALLOCATE(f1d(nc_max+1),ierr)
    
    
    delta1 = sim%mesh_x1%delta_eta
    delta2 = sim%mesh_x2%delta_eta
    delta3 = sim%mesh_x3%delta_eta
    delta4 = sim%mesh_x4%delta_eta
    
   
    
    if(sll_get_collective_rank(sll_world_collective)==0) then
      !print *,'world_size=',world_size
      !print *,'#not implemented for the moment!'
    endif
    SLL_ALLOCATE(recv_sz(world_size),ierr)
    SLL_ALLOCATE(disps(world_size),ierr)


    ! allocate the layouts...
    sequential_x1x2  => new_layout_4D( sll_world_collective )
    sequential_x3x4  => new_layout_4D( sll_world_collective )
    layout2d_par_x1x2 => new_layout_2D( sll_world_collective )
    layout2d_par_x3x4 => new_layout_2D( sll_world_collective )
    
    power2 = int(log(real(world_size))/log(2.0))
    ! special case N = 1, so power2 = 0
    if(power2 == 0) then
       nproc_x1 = 1
       nproc_x2 = 1
       nproc_x3 = 1
       nproc_x4 = 1
    end if    
    if(is_even(power2)) then
       nproc_x1 = 2**(power2/2)
       nproc_x2 = 2**(power2/2)
       nproc_x3 = 1
       nproc_x4 = 1
    else 
       nproc_x1 = 2**((power2-1)/2)
       nproc_x2 = 2**((power2+1)/2)
       nproc_x3 = 1
       nproc_x4 = 1
    end if

    call initialize_layout_with_distributed_2D_array( &
         nc_x1+1, & 
         nc_x2+1, & 
         nproc_x1, &
         nproc_x2, &
         layout2d_par_x1x2 )
    call compute_local_sizes( layout2d_par_x1x2, loc_sz_x1, loc_sz_x2)
    SLL_ALLOCATE(rho_split(loc_sz_x1,loc_sz_x2),ierr)
    SLL_ALLOCATE(send_buf_x1x2(loc_sz_x1*loc_sz_x2), ierr)



    
    
    
    SLL_ALLOCATE(rho_full(nc_x1+1,nc_x2+1),ierr)
    SLL_ALLOCATE(intfdx_full(nc_x3+1,nc_x4+1),ierr)
    SLL_ALLOCATE(E_x1(nc_x1+1,nc_x2+1),ierr)
    SLL_ALLOCATE(E_x2(nc_x1+1,nc_x2+1),ierr)
    SLL_ALLOCATE(recv_buf_x1x2((nc_x1+1)*(nc_x2+1)),ierr)
    SLL_ALLOCATE(recv_buf_x3x4((nc_x3+1)*(nc_x4+1)),ierr)



    call initialize_layout_with_distributed_4D_array( &
         nc_x1+1, &
         nc_x2+1, &
         nc_x3+1, &
         nc_x4+1, &
         nproc_x1, &
         nproc_x2, &
         nproc_x3, &
         nproc_x4, &
         sequential_x3x4 )
    
    call compute_local_sizes_4d( sequential_x3x4, &
         loc_sz_x1, &
         loc_sz_x2, &
         loc_sz_x3, &
         loc_sz_x4 )
    
    SLL_ALLOCATE(f_seq_x3x4(loc_sz_x1,loc_sz_x2,loc_sz_x3,loc_sz_x4),ierr)

    itemp = nproc_x3
    nproc_x3 = nproc_x1
    nproc_x1 = itemp
    itemp = nproc_x4
    nproc_x4 = nproc_x2 
    nproc_x2 = itemp
    
    call initialize_layout_with_distributed_4D_array( &
         nc_x1+1, & 
         nc_x2+1, & 
         nc_x3+1, &
         nc_x4+1, &
         nproc_x1, &
         nproc_x2, &
         nproc_x3, &
         nproc_x4, &
         sequential_x1x2 )
    call compute_local_sizes_4d( sequential_x1x2, &
         loc_sz_x1, &
         loc_sz_x2, &
         loc_sz_x3, &
         loc_sz_x4 )
    SLL_ALLOCATE(f_seq_x1x2(loc_sz_x1,loc_sz_x2,loc_sz_x3,loc_sz_x4),ierr)

    call initialize_layout_with_distributed_2D_array( &
         nc_x3+1, & 
         nc_x4+1, & 
         nproc_x3, &
         nproc_x4, &
         layout2d_par_x3x4 )
    call compute_local_sizes( layout2d_par_x3x4, loc_sz_x3, loc_sz_x4)
    SLL_ALLOCATE(intfdx_split(loc_sz_x3,loc_sz_x4),ierr)
    SLL_ALLOCATE(send_buf_x3x4(loc_sz_x3*loc_sz_x4), ierr)




    call sll_4d_parallel_array_initializer_cartesian( &
         sequential_x3x4, &
         sim%mesh_x1, &
         sim%mesh_x2, &
         sim%mesh_x3, &
         sim%mesh_x4, &
         f_seq_x3x4, &
         sim%init_func, &
         sim%params)
    
    call compute_local_sizes( layout2d_par_x1x2, loc_sz_x1, loc_sz_x2)    
    call compute_reduction_4d_to_2d_direction34(&
      f_seq_x3x4, &
      rho_split, &
      loc_sz_x1, &
      loc_sz_x2, &
      nc_x3+1, &
      nc_x4+1, &
      delta3, &    
      delta4) 
    call load_buffer_2d( layout2d_par_x1x2, rho_split, send_buf_x1x2 )
    recv_sz(:) = receive_counts_array_2d( layout2d_par_x1x2, world_size )    
    send_size = size(send_buf_x1x2)
    call compute_displacements_array_2d( &
         layout2d_par_x1x2, &
         world_size, &
         disps )
    call sll_collective_allgatherv_real64( &
         sll_world_collective, &
         send_buf_x1x2, &
         send_size, &
         recv_sz, &
         disps, &
         recv_buf_x1x2 )
    call unload_buffer_2d(layout2d_par_x1x2, recv_buf_x1x2, rho_full)

    call solve(sim%poisson,E_x1,E_x2,rho_full,nrj)


    seqx3x4_to_seqx1x2 => &
         NEW_REMAP_PLAN(sequential_x3x4,sequential_x1x2,f_seq_x3x4)
    
    seqx1x2_to_seqx3x4 => &
         NEW_REMAP_PLAN(sequential_x1x2,sequential_x3x4,f_seq_x1x2)


    ekin = (sim%mesh_x1%eta_max-sim%mesh_x1%eta_min) &
      *(sim%mesh_x2%eta_max-sim%mesh_x2%eta_min) ! analytic expression
    intfdx_full = 0._f64 ! to compute correctly
    
    nrj0=sim%nrj0
    ekin0=ekin
    
     
    if(sll_get_collective_rank(sll_world_collective)==0) then
      call sll_ascii_file_create('thdiag.dat', th_diag_id, ierr)
      call sll_binary_file_create('rho.bdat', rho_id, ierr)
      call sll_binary_file_create('E_x1.bdat', E_x1_id, ierr)
      call sll_binary_file_create('E_x2.bdat', E_x2_id, ierr)
      call sll_binary_file_create('intfdx.bdat', rho_id, ierr)
      call sll_binary_write_array_2d(rho_id,rho_full(1:nc_x1,1:nc_x2),ierr)  
      call sll_binary_write_array_2d(E_x1_id,E_x1(1:nc_x1,1:nc_x2),ierr)  
      call sll_binary_write_array_2d(E_x2_id,E_x2(1:nc_x1,1:nc_x2),ierr)  
      call sll_binary_write_array_2d(intfdx_id,intfdx_full(1:nc_x3+1,1:nc_x4+1),ierr)  
      write(th_diag_id,'(f12.5,5g20.12)') time, nrj,ekin,nrj0,ekin0!0.5*nrj+ekin
    endif

    
    !print *,'#nrj=',nrj
    call apply_remap_4D( seqx3x4_to_seqx1x2, f_seq_x3x4, f_seq_x1x2 )
    call compute_local_sizes_4d( sequential_x1x2, &
         loc_sz_x1, &
         loc_sz_x2, &
         loc_sz_x3, &
         loc_sz_x4 )
    
    do istep = 1, sim%num_iterations
      if (mod(istep-1,sim%freq_diag)==0) then
        if(sll_get_collective_rank(sll_world_collective)==0) then        
          print *,'#step=',istep-1,real(istep-1,f64)*sim%dt
        endif
      endif  
      
      split_T = sim%split%split_begin_T      
      do split_istep=1,sim%split%nb_split_step
        if(split_T)then
          !T advection
          global_indices(1:4) = local_to_global_4D( sequential_x1x2, (/1, 1, 1, 1/) )
          do i4=1,loc_sz_x4
            do i3=1,loc_sz_x3
              !advection in x1
              ig=global_indices(3)
              alpha = (sim%mesh_x3%eta_min + real(i3+ig-2,f64) * delta3) &
                * sim%split%split_step(split_istep)
              do i2=1,nc_x2+1
                f1d(1:nc_x1+1)=f_seq_x1x2(1:nc_x1+1,i2,i3,i4) 
                call sim%advect_x1%advect_1d_constant(&
                  alpha, &
                  sim%dt, &
                  f1d(1:nc_x1+1), &
                  f1d(1:nc_x1+1))
                f_seq_x1x2(1:nc_x1+1,i2,i3,i4)=f1d(1:nc_x1+1)
              enddo
              !advection in x2
              ig=global_indices(4)
              alpha = (sim%mesh_x4%eta_min + real(i4+ig-2,f64) * delta4) &
                * sim%split%split_step(split_istep)
              do i1=1,nc_x1+1
                f1d(1:nc_x2+1)=f_seq_x1x2(i1,1:nc_x2+1,i3,i4) 
                call sim%advect_x2%advect_1d_constant(&
                  alpha, &
                  sim%dt, &
                  f1d(1:nc_x2+1), &
                  f1d(1:nc_x2+1))
                f_seq_x1x2(i1,1:nc_x2+1,i3,i4)=f1d(1:nc_x2+1)
              enddo                            
            enddo 
          enddo
        else
          !V advection
          call apply_remap_4D( seqx1x2_to_seqx3x4, f_seq_x1x2, f_seq_x3x4 )
          call compute_local_sizes_4d( sequential_x3x4, &
           loc_sz_x1, &
           loc_sz_x2, &
           loc_sz_x3, &
           loc_sz_x4 )
          
          !begin compute poisson
          call compute_reduction_4d_to_2d_direction34(&
            f_seq_x3x4, &
            rho_split, &
            loc_sz_x1, &
            loc_sz_x2, &
            nc_x3+1, &
            nc_x4+1, &
            delta3, &    
            delta4) 
          call load_buffer_2d( layout2d_par_x1x2, rho_split, send_buf_x1x2 )
          recv_sz(:) = receive_counts_array_2d( layout2d_par_x1x2, world_size )    
          send_size = size(send_buf_x1x2)
          call compute_displacements_array_2d( &
            layout2d_par_x1x2, &
            world_size, &
            disps )
          call sll_collective_allgatherv_real64( &
            sll_world_collective, &
            send_buf_x1x2, &
            send_size, &
            recv_sz, &
            disps, &
            recv_buf_x1x2 )
          call unload_buffer_2d(layout2d_par_x1x2, recv_buf_x1x2, rho_full)
          call solve(sim%poisson,E_x1,E_x2,rho_full,nrj)
          !end compute poisson 
          
          
          
          
          
          
          global_indices(1:4) = local_to_global_4D( sequential_x3x4, (/1, 1, 1, 1/) ) 
          do i2=1,loc_sz_x2
            do i1=1,loc_sz_x1

              !advection in x3
              alpha = E_x1(i1-1+global_indices(1),i2-1+global_indices(2)) &
                * sim%split%split_step(split_istep)
              do i4=1,nc_x4+1
                f1d(1:nc_x3+1)=f_seq_x3x4(i1,i2,1:nc_x3+1,i4) 
                call sim%advect_x3%advect_1d_constant(&
                  alpha, &
                  sim%dt, &
                  f1d(1:nc_x3+1), &
                  f1d(1:nc_x3+1))
                f_seq_x3x4(i1,i2,1:nc_x3+1,i4)=f1d(1:nc_x3+1)
              enddo
              !advection in x4
              alpha = E_x2(i1-1+global_indices(1),i2-1+global_indices(2)) &
                * sim%split%split_step(split_istep)
              do i3=1,nc_x3+1
                f1d(1:nc_x4+1)=f_seq_x3x4(i1,i2,i3,1:nc_x4+1) 
                call sim%advect_x4%advect_1d_constant(&
                  alpha, &
                  sim%dt, &
                  f1d(1:nc_x4+1), &
                  f1d(1:nc_x4+1))
                f_seq_x3x4(i1,i2,i3,1:nc_x4+1)=f1d(1:nc_x4+1)
              enddo                            






              
            enddo 
          enddo

          call apply_remap_4D( seqx3x4_to_seqx1x2, f_seq_x3x4, f_seq_x1x2 )
          call compute_local_sizes_4d( sequential_x1x2, &
           loc_sz_x1, &
           loc_sz_x2, &
           loc_sz_x3, &
           loc_sz_x4 )
        endif
        split_T = .not.(split_T)  
      enddo
      if (mod(istep,sim%freq_diag_time)==0) then
        time = real(istep,f64)*sim%dt
        
        
        
        call compute_reduction_4d_to_2d_direction12(&
            f_seq_x1x2, &
            intfdx_split, &
            loc_sz_x1, &
            loc_sz_x2, &
            loc_sz_x3, &
            loc_sz_x4, &
            delta1, &    
            delta2) 
        call load_buffer_2d( layout2d_par_x3x4, intfdx_split, send_buf_x3x4 )
        recv_sz(:) = receive_counts_array_2d( layout2d_par_x3x4, world_size )    
        send_size = size(send_buf_x3x4)
        call compute_displacements_array_2d( &
            layout2d_par_x3x4, &
            world_size, &
            disps )
        call sll_collective_allgatherv_real64( &
            sll_world_collective, &
            send_buf_x3x4, &
            send_size, &
            recv_sz, &
            disps, &
            recv_buf_x3x4 )
        call unload_buffer_2d(layout2d_par_x3x4, recv_buf_x3x4, intfdx_full)
        if (mod(istep,sim%freq_diag)==0) then
          if(sll_get_collective_rank(sll_world_collective)==0) then
            call sll_binary_write_array_2d(intfdx_id,intfdx_full(1:nc_x3+1,1:nc_x4+1),ierr)
          endif
        endif


        do i4=1,nc_x4+1
          do i3=1,nc_x3+1
            tmp = (sim%mesh_x4%eta_min + real(i4-1,f64) * delta4)**2
            tmp = tmp+(sim%mesh_x3%eta_min + real(i3-1,f64) * delta3)**2
            intfdx_full(i3,i4) = intfdx_full(i3,i4) * 0.5_f64*tmp 
          enddo
        enddo
        
        call compute_reduction_2d_to_0d(&
            intfdx_full, &
            ekin, &
            nc_x3+1, &
            nc_x4+1, &
            delta3, &    
            delta4) 
        
        
        
        
        
        if(sll_get_collective_rank(sll_world_collective)==0) then
          write(th_diag_id,'(f12.5,5g20.12)') time, nrj,ekin,nrj0,ekin0!,ekin+0.5_f64*nrj
        endif
      endif
      if (mod(istep,sim%freq_diag)==0) then
        if(sll_get_collective_rank(sll_world_collective)==0) then
          call sll_binary_write_array_2d(rho_id,rho_full(1:nc_x1,1:nc_x2),ierr)  
          call sll_binary_write_array_2d(E_x1_id,E_x1(1:nc_x1,1:nc_x2),ierr)  
          call sll_binary_write_array_2d(E_x2_id,E_x2(1:nc_x1,1:nc_x2),ierr)  
        endif



      endif   
    enddo
    
    
    if(sll_get_collective_rank(sll_world_collective)==0)then
      call sll_ascii_file_close(th_diag_id,ierr) 
    endif
    
  end subroutine run_vp4d_cartesian    
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  

end module sll_simulation_4d_vlasov_parallel_poisson_sequential_cartesian
