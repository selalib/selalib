module sll_simulation_2d_analytic_field_cartesian_module

! for swirling deformation flow test

!contact: Michel Mehrenberger (mehrenbe@math.unistra.fr)


#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_field_2d.h"
#include "sll_utilities.h"
#include "sll_poisson_solvers.h"
  use sll_constants
  use sll_logical_meshes  
  use sll_module_advection_1d_periodic
  use sll_module_advection_2d_BSL
  use sll_module_characteristics_2d_explicit_euler
  use sll_module_characteristics_2d_verlet
  use sll_reduction_module
  use sll_simulation_base
  use sll_cubic_spline_interpolator_2d
  use sll_cubic_spline_interpolator_1d
  use sll_coordinate_transformation_2d_base_module
  use sll_module_coordinate_transformations_2d
  use sll_common_coordinate_transformations
  use sll_common_array_initializers_module
  !use sll_mudpack_curvilinear
  use sll_module_poisson_2d_mudpack_solver
  use sll_module_poisson_2d_mudpack_curvilinear_solver_old
  use sll_module_poisson_2d_elliptic_solver
  use sll_module_scalar_field_2d_base
  use sll_module_scalar_field_2d_alternative
  use sll_timer
  use sll_fft
  use sll_module_poisson_2d_periodic_solver

  implicit none

  
  sll_int32, parameter :: SLL_EULER = 0 
  sll_int32, parameter :: SLL_PREDICTOR_CORRECTOR = 1 
  sll_int32, parameter :: SLL_LEAP_FROG = 2 
  sll_int32, parameter :: SLL_PHI_FROM_RHO = 0
  sll_int32, parameter :: SLL_E_FROM_RHO = 1


  type, extends(sll_simulation_base_class) :: &
    sll_simulation_2d_analytic_field_cartesian

   !geometry
   type(sll_logical_mesh_2d), pointer :: mesh_2d


   !initial function
   procedure(sll_scalar_initializer_2d), nopass, pointer :: init_func
   sll_real64, dimension(:), pointer :: params
      
   !advector
   class(sll_advection_2d_base), pointer    :: advect_2d
   procedure(sll_scalar_initializer_2d), nopass, pointer :: A1_func
   procedure(sll_scalar_initializer_2d), nopass, pointer :: A2_func
   sll_real64, dimension(:), pointer :: A_func_params
   procedure(sll_scalar_initializer_1d), nopass, pointer :: A_time_func
   sll_real64, dimension(:), pointer :: A_time_func_params
   
   !interpolator for derivatives
   class(sll_interpolator_2d_base), pointer   :: phi_interp2d

   !time_iterations
   sll_real64 :: dt
   sll_int32  :: num_iterations
   sll_int32  :: freq_diag
   sll_int32  :: freq_diag_time

   !time_loop
   sll_int32 :: time_loop_case
   
       
  contains
    procedure, pass(sim) :: run => run_af2d_cartesian
    procedure, pass(sim) :: init_from_file => init_fake
     
  end type sll_simulation_2d_analytic_field_cartesian


  abstract interface
    function sll_scalar_initializer_2d( x1, x2, params )
      use sll_working_precision
      sll_real64                                     :: sll_scalar_initializer_2d
      sll_real64, intent(in)                         :: x1
      sll_real64, intent(in)                         :: x2
      sll_real64, dimension(:), intent(in), optional :: params
    end function sll_scalar_initializer_2d
  end interface
  abstract interface
    function sll_scalar_initializer_1d( x1,  params )
      use sll_working_precision
      sll_real64                                     :: sll_scalar_initializer_1d
      sll_real64, intent(in)                         :: x1
      sll_real64, dimension(:), intent(in), optional :: params
    end function sll_scalar_initializer_1d
  end interface



contains

  function new_analytic_field_2d_cartesian(filename) result(sim)
    type(sll_simulation_2d_analytic_field_cartesian), pointer :: sim    
    character(len=*), intent(in), optional :: filename
    sll_int32 :: ierr
    
    SLL_ALLOCATE(sim,ierr)
    
    call initialize_analytic_field_2d_cartesian(sim,filename)
    
  
  
  end function new_analytic_field_2d_cartesian
  
  subroutine initialize_analytic_field_2d_cartesian(sim, filename)
    class(sll_simulation_2d_analytic_field_cartesian), intent(inout) :: sim
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
    character(len=256) :: phi_interp2d_case
    character(len=256) ::  charac2d_case
    character(len=256) ::  A_interp_case
    character(len=256) ::  advection_field_case
    sll_real64 :: time_period

 

    !local variables
    sll_int32 :: Nc_x1
    sll_int32 :: Nc_x2
    type(sll_logical_mesh_1d), pointer :: mesh_x1
    type(sll_logical_mesh_1d), pointer :: mesh_x2
    class(sll_interpolator_2d_base), pointer :: f_interp2d
    class(sll_interpolator_2d_base), pointer :: phi_interp2d
    class(sll_characteristics_2d_base), pointer :: charac2d
    class(sll_interpolator_2d_base), pointer   :: A1_interp2d
    class(sll_interpolator_2d_base), pointer   :: A2_interp2d
    class(sll_interpolator_1d_base), pointer   :: A1_interp1d_x1
    class(sll_interpolator_1d_base), pointer   :: A2_interp1d_x1
    sll_int32 :: ierr

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
      phi_interp2d_case, &
      charac2d_case, &
      A_interp_case, &
      advection_field_case, &
      time_period


    !! set default parameters
    
    !geometry
    mesh_case_x1="SLL_LOGICAL_MESH"
    num_cells_x1 = 32
    x1_min = 0.0_f64
    x1_max = 2._f64*sll_pi
    mesh_case_x2="SLL_LOGICAL_MESH"
    num_cells_x2 = 32
    x2_min = 0.0_f64
    x2_max = 2._f64*sll_pi
    
    !initial function
    initial_function_case="SLL_KHP1"
    kmode_x1 = 0.5_f64
    kmode_x2 = 1._f64
    eps = 0.015_f64
    
    
    !time_iterations
    dt = 0.1_f64
    number_iterations  = 600
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
    advection_field_case = "SLL_SWIRLING_DEFORMATION_FLOW"
    

    if(present(filename))then
      open(unit = input_file, file=trim(filename)//'.nml',IOStat=IO_stat)
        if( IO_stat /= 0 ) then
          print *, '#initialize_guiding_center_2d_cartesian() failed to open file ', &
          trim(filename)//'.nml'
          STOP
        end if
      print *,'#initialization with filename:'
      print *,'#',trim(filename)//'.nml'
      read(input_file, geometry) 
      read(input_file, initial_function)
      read(input_file, time_iterations)
      read(input_file, advector)
      close(input_file)
    else
      print *,'#initialization with default parameters'    
    endif

    Nc_x1 = num_cells_x1
    Nc_x2 = num_cells_x2
     
    
    sim%dt = dt
    sim%num_iterations = number_iterations
    sim%freq_diag = freq_diag
    sim%freq_diag_time = freq_diag_time

    select case (mesh_case_x1)
      case ("SLL_LOGICAL_MESH")
        mesh_x1 => new_logical_mesh_1d(num_cells_x1,eta_min=x1_min, eta_max=x1_max)  
      case default
        print*,'#mesh_case_x1', mesh_case_x1, ' not implemented'
        stop 
    end select
    select case (mesh_case_x2)
      case ("SLL_LOGICAL_MESH")
        mesh_x2 => new_logical_mesh_1d(num_cells_x2,eta_min=x2_min, eta_max=x2_max)
      case default
        print*,'#mesh_case_x2', mesh_case_x2, ' not implemented'
        stop 
    end select
    sim%mesh_2d => tensor_product_1d_1d( mesh_x1, mesh_x2)



!    sim%mesh_2d => new_logical_mesh_2d( &
!      num_cells_x1, &
!      num_cells_x2, &
!      eta1_min = x1_min, &
!      eta1_max = x1_max, &
!      eta2_min = x2_min, &
!      eta2_max = x2_max)      
      
      
      
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

    select case(advection_field_case)
      case ("SLL_SWIRLING_DEFORMATION_FLOW")
        sim%A1_func => sll_SDF_A1_initializer_2d 
        sim%A2_func => sll_SDF_A2_initializer_2d 
        SLL_ALLOCATE(sim%A_func_params(2),ierr)
        sim%A_time_func => sll_SDF_time_initializer_1d 
        SLL_ALLOCATE(sim%A_time_func_params(1),ierr)
        sim%A_time_func_params(1) = time_period
      case default
        print *,'#bad advect_case',advect2d_case
        print *,'#not implemented'
        print *,'#in initialize_guiding_center_2d_polar'
        stop
    end select

    
    
    select case(initial_function_case)
      case ("SLL_KHP1")
        sim%init_func => sll_KHP1_2d
        SLL_ALLOCATE(sim%params(2),ierr)
        sim%params(1) = eps
        sim%params(2) = kmode_x1
      case ("SLL_GAUSSIAN")
        sim%init_func => sll_gaussian_initializer_2d
        SLL_ALLOCATE(sim%params(4),ierr)
        sim%params(1) = xc_1
        sim%params(2) = xc_2
        sim%params(3) = sigma_1
        sim%params(4) = sigma_2
      case ("SLL_COS_BELL")
        sim%init_func => sll_cos_bell_initializer_2d
        SLL_ALLOCATE(sim%params(2),ierr)
        sim%params(1) = xc_1
        sim%params(2) = xc_2
      case default
        print *,'#bad initial_function_case',initial_function_case
        print *,'#not implemented'
        print *,'#in initialize_analytic_field_2d_cartesian'
        stop
    end select
    
    
    !time_loop
    select case(time_loop_case)
      case ("SLL_EULER")
        sim%time_loop_case = SLL_EULER
      case ("SLL_PREDICTOR_CORRECTOR")
        sim%time_loop_case = SLL_PREDICTOR_CORRECTOR
      case default
        print *,'#bad time_loop_case',time_loop_case
        print *,'#not implemented'
        print *,'#in initialize_guiding_center_2d_polar'
        stop
    end select

   
  end subroutine initialize_analytic_field_2d_cartesian
  


  subroutine init_fake(sim, filename)
    class(sll_simulation_2d_analytic_field_cartesian), intent(inout) :: sim
    character(len=*), intent(in)                                :: filename
    
    print *,'# Do not use the routine init_vp4d_fake'
    print *,'#use instead initialize_vlasov_par_poisson_seq_cart'
    print *,'#filename=',filename
    print *,sim%dt
    stop
  
  end subroutine init_fake
  
  subroutine run_af2d_cartesian(sim)
    class(sll_simulation_2d_analytic_field_cartesian), intent(inout) :: sim
    sll_int32 :: Nc_x1
    sll_int32 :: Nc_x2
    sll_real64 :: delta_x1
    sll_real64 :: delta_x2
    sll_real64 :: x1_min,x1_max
    sll_real64 :: x2_min,x2_max   
    sll_real64 :: x1
    sll_real64 :: x2
    sll_int32 :: i1 
    sll_int32 :: i2
    sll_real64,dimension(:,:), pointer :: f
    sll_real64,dimension(:,:), pointer :: f_old
    sll_real64,dimension(:,:), pointer :: f_init
    sll_real64,dimension(:,:), pointer :: phi
    sll_real64,dimension(:,:), pointer :: A1 !advection fields
    sll_real64,dimension(:,:), pointer :: A2
    sll_real64,dimension(:,:), pointer :: A1_init !advection fields
    sll_real64,dimension(:,:), pointer :: A2_init
    sll_int32 :: ierr
    sll_int32 :: nb_step
    sll_int32 :: step
    sll_real64 :: dt
    sll_int32 :: diag_id = 88 
    sll_int32             :: IO_stat
    sll_int32 :: iplot
    sll_real64 :: time_factor
    
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
    
    
    !allocation
    SLL_ALLOCATE(f(Nc_x1+1,Nc_x2+1),ierr)
    SLL_ALLOCATE(f_old(Nc_x1+1,Nc_x2+1),ierr)
    SLL_ALLOCATE(f_init(Nc_x1+1,Nc_x2+1),ierr)
    SLL_ALLOCATE(phi(Nc_x1+1,Nc_x2+1),ierr)
    SLL_ALLOCATE(A1(Nc_x1+1,Nc_x2+1),ierr)
    SLL_ALLOCATE(A2(Nc_x1+1,Nc_x2+1),ierr)
    SLL_ALLOCATE(A1_init(Nc_x1+1,Nc_x2+1),ierr)
    SLL_ALLOCATE(A2_init(Nc_x1+1,Nc_x2+1),ierr)

    

    
    !initialisation of distribution function
    do i2=1,Nc_x2+1
      x2=x2_min+real(i2-1,f64)*delta_x2
      do i1=1,Nc_x1+1
        x1=x1_min+real(i1-1,f64)*delta_x1
        f(i1,i2) =  sim%init_func(x1,x2,sim%params)
        f_init(i1,i2) =  sim%init_func(x1,x2,sim%params)
        A1_init(i1,i2) =  sim%A1_func(x1,x2,sim%A_func_params)
        A2_init(i1,i2) =  sim%A2_func(x1,x2,sim%A_func_params)
      end do
    end do
        

    
    
    open(unit = diag_id, file='thdiag.dat',IOStat=IO_stat)
    if( IO_stat /= 0 ) then
       print *, '#run_af2d_cartesian (sim) failed to open file thdiag.dat'
       STOP
    end if
    
    iplot = 0

    do step=0,nb_step-1
      f_old = f

     
      
      if(modulo(step-1,sim%freq_diag_time)==0)then
      endif            
     
#ifndef NOHDF5
      if(modulo(step,sim%freq_diag)==0)then
        print*,"#step= ", step
        call plot_f_cartesian(iplot,f,sim%mesh_2d)
        iplot = iplot+1  
      endif            
#endif
      
      select case (sim%time_loop_case)
        case (SLL_EULER)
          time_factor = sim%A_time_func( &
            real(step,f64)*sim%dt, &
            sim%A_time_func_params )
          A1 = time_factor*A1_init
          A2 = time_factor*A2_init          
          call sim%advect_2d%advect_2d(A1, A2, sim%dt, f_old, f)
        case (SLL_PREDICTOR_CORRECTOR)
          time_factor = sim%A_time_func( &
            (real(step,f64)+0.5_f64)*sim%dt, &
            sim%A_time_func_params )
          A1 = time_factor*A1_init
          A2 = time_factor*A2_init          
          call sim%advect_2d%advect_2d(A1, A2, sim%dt, f_old, f)
        case default  
          print *,'#bad time_loop_case',sim%time_loop_case
          print *,'#not implemented'
          print *,'#in run_af2d_cartesian'
          print *,'#available options are:'
          print *,'#SLL_EULER=',SLL_EULER
          print *,'#SLL_PREDICTOR_CORRECTOR=',SLL_PREDICTOR_CORRECTOR
          
      end select
         
    enddo
    
    close(diag_id)
    
    print *,maxval(abs(f-f_init))
     
    print *,'#run_af2d_cartesian PASSED'
  end subroutine run_af2d_cartesian  
  

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
    !sll_real64 :: r
    !sll_real64 :: theta
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


end module sll_simulation_2d_analytic_field_cartesian_module
