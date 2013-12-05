module sll_simulation_2d_guiding_center_curvilinear_mudpack_module

!the aim is to create guiding center cartesian in simulation class
!related to
!simulation_guiding_center_2D_generalized_coords_mudpack.F90
!but here geometry and test is specifically cartesian


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
  use sll_mudpack_curvilinear
  !use sll_parallel_array_initializer_module

  implicit none
  
  
  sll_int32, parameter :: SLL_EULER = 0 
  sll_int32, parameter :: SLL_PREDICTOR_CORRECTOR = 1 


  type, extends(sll_simulation_base_class) :: &
    sll_simulation_2d_guiding_center_curvilinear_mudpack

   !geometry
   type(sll_logical_mesh_2d), pointer :: mesh_2d
  
   !transformation 
    class(sll_coordinate_transformation_2d_base), pointer :: transformation
  
   !initial function
   procedure(sll_scalar_initializer_2d), nopass, pointer :: init_func
   sll_real64, dimension(:), pointer :: params
      
   !advector
   class(sll_advection_2d_base), pointer    :: advect_2d
   
   !interpolator for derivatives
   class(sll_interpolator_2d_base), pointer   :: phi_interp2d
   !coef
   sll_real64, dimension(:,:), pointer :: b11
   sll_real64, dimension(:,:), pointer :: b12
   sll_real64, dimension(:,:), pointer :: b21
   sll_real64, dimension(:,:), pointer :: b22
   sll_real64, dimension(:,:), pointer :: c
   !poisson solver
   !class(sll_poisson_2d_base), pointer   :: poisson
   !type(poisson_2d_periodic), pointer   :: poisson
   !type(sll_plan_poisson_polar), pointer :: poisson 
    type(mudpack_2d) :: poisson 
   !time_iterations
   sll_real64 :: dt
   sll_int32  :: num_iterations
   sll_int32  :: freq_diag
   sll_int32  :: freq_diag_time

   !time_loop
   sll_int32 :: time_loop_case
   
       
  contains
    procedure, pass(sim) :: run => run_gc2d_curvilinear_mudpack
    procedure, pass(sim) :: init_from_file => init_fake
     
  end type sll_simulation_2d_guiding_center_curvilinear_mudpack


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

  function new_guiding_center_2d_curvilinear_mudpack() result(sim)
    type(sll_simulation_2d_guiding_center_curvilinear_mudpack), pointer :: sim    
    sll_int32 :: ierr
    
    SLL_ALLOCATE(sim,ierr)
    
    call initialize_guiding_center_2d_curvilinear_mudpack(sim)
    
  
  
  end function new_guiding_center_2d_curvilinear_mudpack
  
  subroutine initialize_guiding_center_2d_curvilinear_mudpack(sim)
    class(sll_simulation_2d_guiding_center_curvilinear_mudpack), intent(inout) :: sim
    sll_int32 :: Nc_x1
    sll_int32 :: Nc_x2
    sll_real64 :: x1_min
    sll_real64 :: x1_max
    sll_real64 :: x2_min
    sll_real64 :: x2_max
    sll_real64 :: r_minus
    sll_real64 :: r_plus
    sll_real64 :: eps
    sll_real64 :: k_mode
    sll_int32  :: nb_step
    sll_int32  :: i,j
    sll_real64 :: dt
    sll_int32 :: visu_step
    class(sll_interpolator_2d_base), pointer :: f_interp2d
    class(sll_interpolator_2d_base), pointer :: phi_interp2d
    class(sll_characteristics_2d_base), pointer :: charac2d
    class(sll_interpolator_2d_base), pointer   :: A1_interp2d
    class(sll_interpolator_2d_base), pointer   :: A2_interp2d
    class(sll_interpolator_1d_base), pointer   :: A1_interp1d_x1
    class(sll_interpolator_1d_base), pointer   :: A2_interp1d_x1
    character(len=256)      :: advect2d_case 
    character(len=256)      :: charac2d_case
    character(len=256)      :: f_interp2d_case 
    character(len=256)      :: phi_interp2d_case 
    character(len=256)      :: A_interp_case 
    character(len=256)      :: initial_function_case 
    character(len=256)      :: time_loop_case 
    character(len=256)      :: poisson_case
    sll_int32 :: ierr
    !character(len=256)      :: interp1d_x2_case 
    sll_real64, dimension(4) :: params_mesh
    !here we do all the initialization
    !in future, we will use namelist file

    
    Nc_x1 = 128
    Nc_x2 = 128
    k_mode = 0.5_f64
    eps = 0.015_f64
    x1_min = 0._f64
    x1_max = 2._f64*sll_pi/k_mode
    x2_min = 0._f64
    x2_max = 2._f64*sll_pi
    nb_step = 600
    
    dt = 0.1_f64
    visu_step = 100
    f_interp2d_case = "SLL_CUBIC_SPLINES"
    phi_interp2d_case = "SLL_CUBIC_SPLINES"
    A_interp_case = "SLL_CUBIC_SPLINES"
    charac2d_case = "SLL_VERLET"
    !charac2d_case = "SLL_EULER"
    advect2d_case = "SLL_BSL"
    initial_function_case = "SLL_KHP1" 
    !time_loop_case = "SLL_EULER"
    time_loop_case = "SLL_PREDICTOR_CORRECTOR" 
    poisson_case = "MUDPACK"   
    
    sim%dt = dt
    sim%num_iterations = nb_step
    sim%freq_diag = visu_step
    sim%freq_diag_time = 1
    
    
    SLL_ALLOCATE(sim%b11(Nc_x1+1,Nc_x2+1),ierr)
    SLL_ALLOCATE(sim%b12(Nc_x1+1,Nc_x2+1),ierr)
    SLL_ALLOCATE(sim%b21(Nc_x1+1,Nc_x2+1),ierr)
    SLL_ALLOCATE(sim%b22(Nc_x1+1,Nc_x2+1),ierr)
    SLL_ALLOCATE(sim%c(Nc_x1+1,Nc_x2+1),ierr)
    do j=1,Nc_x2+1
     do i=1,Nc_x1+1
        sim%b11(i,j)= 1._f64
        sim%b22(i,j)= 1._f64
        sim%b12(i,j)= 0._f64
        sim%b21(i,j)= 0._f64
        sim%c(i,j) = 0._f64
     enddo
    enddo    
    
    !  In collela  mesh params_mesh =( alpha1, alpha2, L1, L2 ) such that :
    !  x1= eta1 + alpha1*sin(2*pi*eta1/L1)*sin(2*pi*eta2/L2)
    params_mesh = (/ 0.1_f64, 0.1_f64, 1.0_f64, 1.0_f64/)
    
    sim%mesh_2d => new_logical_mesh_2d( &
      Nc_x1, &
      Nc_x2, &
      eta1_min = x1_min, &
      eta1_max = x1_max, &
      eta2_min = x2_min, &
      eta2_max = x2_max) 
           
    sim%transformation => new_coordinate_transformation_2d_analytic( &
       "analytic_identity_transformation", &
       sim%mesh_2d, &
       identity_x1, &
       identity_x2, &
       identity_jac11, &
       identity_jac12, &
       identity_jac21, &
       identity_jac22, &
       params_mesh   )  
       
    ! transformation => new_coordinate_transformation_2d_analytic( &
!       "analytic_polar_transformation", &
!       sim%mesh_2d, &
!       polar_x1, &
!       polar_x2, &
!       polar_jac11, &
!       polar_jac12, &
!       polar_jac21, &
!       polar_jac22, &
!       params_mesh  )     

! transformation => new_coordinate_transformation_2d_analytic( &
!       "analytic_collela_transformation", &
!       sim%mesh_2d, &
!       sinprod_x1, &
!       sinprod_x2, &
!       sinprod_jac11, &
!       sinprod_jac12, &
!       sinprod_jac21, &
!       sinprod_jac22, &
!       params_mesh  )  
      
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
    
    
    select case(initial_function_case)
      case ("SLL_KHP1")
        sim%init_func => sll_KHP1_2d
        SLL_ALLOCATE(sim%params(2),ierr)
        sim%params(1) = eps
        sim%params(2) = k_mode
      case default
        print *,'#bad initial_function_case',initial_function_case
        print *,'#not implemented'
        print *,'#in initialize_guiding_center_2d_polar'
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
    
    
    
    
    !poisson solver
     !poisson solver
    select case(poisson_case)    
      case ("MUDPACK")     
        call initialize_mudpack_curvilinear(sim%poisson,&
         sim%transformation, &
         sim%b11,&
         sim%b12,&
         sim%b21,&
         sim%b22,&
         sim%c,& 
         x1_min,&
         x1_max,&
         Nc_x1,&
         x2_min,&
         x2_max,&
         Nc_x2,&
         SLL_PERIODIC,& 
         SLL_PERIODIC,& 
         SLL_PERIODIC,& 
         SLL_PERIODIC)
      case default
        print *,'#bad poisson_case',poisson_case
        print *,'#not implemented'
        print *,'#in initialize_mudpack_curvilinear'
        stop
    end select

   
  end subroutine initialize_guiding_center_2d_curvilinear_mudpack
  


  subroutine init_fake(sim, filename)
    class(sll_simulation_2d_guiding_center_curvilinear_mudpack), intent(inout) :: sim
    character(len=*), intent(in)                                :: filename
  
    print *,'# Do not use the routine init_vp4d_fake'
    print *,'#use instead initialize_vlasov_par_poisson_seq_cart'
    stop
  
  end subroutine init_fake
  
  subroutine run_gc2d_curvilinear_mudpack(sim)
    class(sll_simulation_2d_guiding_center_curvilinear_mudpack), intent(inout) :: sim
    sll_int32 :: Nc_x1
    sll_int32 :: Nc_x2
    sll_real64 :: delta_eta1
    sll_real64 :: delta_eta2
    sll_real64 :: eta1_min,eta1_max
    sll_real64 :: eta2_min,eta2_max 
    sll_real64 :: eta1
    sll_real64 :: eta2  
    sll_real64 :: x1
    sll_real64 :: x2
    sll_int32 :: i1 
    sll_int32 :: i2
    sll_real64,dimension(:,:), pointer :: f
    sll_real64,dimension(:,:), pointer :: f_old
    sll_real64,dimension(:,:), pointer :: phi
    sll_real64,dimension(:,:), pointer :: A1 !advection fields
    sll_real64,dimension(:,:), pointer :: A2
    sll_int32 :: ierr
    sll_int32 :: nb_step
    sll_int32 :: step
    sll_real64 :: dt
    sll_int32 :: diag_id = 88 
    sll_int32             :: IO_stat
    sll_int32 :: iplot
    
    Nc_x1 = sim%mesh_2d%num_cells1
    Nc_x2 = sim%mesh_2d%num_cells2
    delta_eta1 = sim%mesh_2d%delta_eta1
    delta_eta2 = sim%mesh_2d%delta_eta2
    eta1_min = sim%mesh_2d%eta1_min
    eta2_min = sim%mesh_2d%eta2_min
    eta1_max = sim%mesh_2d%eta1_max
    eta2_max = sim%mesh_2d%eta2_max
    nb_step = sim%num_iterations
    dt = sim%dt
    
    
    !allocation
    SLL_ALLOCATE(f(Nc_x1+1,Nc_x2+1),ierr)
    SLL_ALLOCATE(f_old(Nc_x1+1,Nc_x2+1),ierr)
    SLL_ALLOCATE(phi(Nc_x1+1,Nc_x2+1),ierr)
    SLL_ALLOCATE(A1(Nc_x1+1,Nc_x2+1),ierr)
    SLL_ALLOCATE(A2(Nc_x1+1,Nc_x2+1),ierr)

    

    
    !initialisation of distribution function    
     do i2=1,Nc_x2+1
        eta2=eta2_min+real(i2-1,f64)*delta_eta2
        do i1=1,Nc_x1+1
          eta1=eta1_min+real(i1-1,f64)*delta_eta1
          x1 = sim%transformation%x1(eta1,eta2)
          x2 = sim%transformation%x2(eta1,eta2)
          f(i1,i2) =  sim%init_func(x1,x2,sim%params) 
        end do
     end do
        
    !solve poisson
    !call poisson_solve_cartesian(sim%poisson,f,phi)
    call solve_mudpack_curvilinear(sim%poisson,phi,-f)
    call compute_field_from_phi_2d_curvilinear_mudpack(phi,sim%mesh_2d,sim%transformation,A1,A2,sim%phi_interp2d)
    print*,"PASSED"
    
    
    !print *,A1
    !print *,A2
    
    open(unit = diag_id, file='diag_curvilinear_mudpack.dat',IOStat=IO_stat)
    if( IO_stat /= 0 ) then
       print *, '#run_gc2d_cartesian (sim) failed to open file diag_curvilinear_mudpack.dat'
       STOP
    end if
    
    iplot = 0

    do step=1,nb_step+1
      print*,"step= ", step
      f_old = f
      
      !call poisson_solve_cartesian(sim%poisson,f_old,phi)
      call solve_mudpack_curvilinear(sim%poisson, phi, -f_old)
      
      call compute_field_from_phi_2d_curvilinear_mudpack(phi,sim%mesh_2d,sim%transformation,A1,A2,sim%phi_interp2d)      
      
      if(modulo(step-1,sim%freq_diag_time)==0)then
        call time_history_diagnostic_gc_cartesian( &
          diag_id, &    
          step-1, &
          dt, &
          sim%mesh_2d, &
          f, &
          phi, &
          A1, &
          A2)
      endif            
     
      if(modulo(step-1,sim%freq_diag)==0)then
        call plot_f_curvilinear(iplot,f,sim%mesh_2d,sim%transformation)
        iplot = iplot+1  
      endif            
      
      select case (sim%time_loop_case)
        case (SLL_EULER)
          call sim%advect_2d%advect_2d(A1, A2, sim%dt, f_old, f)
        case (SLL_PREDICTOR_CORRECTOR)
          call sim%advect_2d%advect_2d(A1, A2, 0.5_f64*sim%dt, f_old, f)
          !call poisson_solve_cartesian(sim%poisson,f,phi)
          call solve_mudpack_curvilinear(sim%poisson, phi, -f)
          call compute_field_from_phi_2d_curvilinear_mudpack(phi,sim%mesh_2d,sim%transformation,A1,A2,sim%phi_interp2d)      
          f_old = f
          call sim%advect_2d%advect_2d(A1, A2, 0.5_f64*sim%dt, f_old, f)
        case default  
          print *,'#bad time_loop_case',sim%time_loop_case
          print *,'#not implemented'
          print *,'#in run_gc2d_cartesian'
          print *,'#available options are:'
          print *,'#SLL_EULER=',SLL_EULER
          print *,'#SLL_PREDICTOR_CORRECTOR=',SLL_PREDICTOR_CORRECTOR
          
      end select
         
    enddo
    
    close(diag_id)

    !print *,'#not implemented for the moment!'
  end subroutine run_gc2d_curvilinear_mudpack 
  
  
  subroutine compute_field_from_phi_2d_curvilinear_mudpack(phi,mesh_2d,transformation,A1,A2,interp2d)
    sll_real64, dimension(:,:), intent(in) :: phi
    sll_real64, dimension(:,:), intent(out) :: A1
    sll_real64, dimension(:,:), intent(out) :: A2
    type(sll_logical_mesh_2d), pointer :: mesh_2d
    class(sll_coordinate_transformation_2d_base), pointer :: transformation
    class(sll_interpolator_2d_base), pointer   :: interp2d
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

    do i2=1,Nc_x2+1
      x2=x2_min+real(i2-1,f64)*delta_x2
      do i1=1,Nc_x1+1
        x1=x1_min+real(i1-1,f64)*delta_x1
        A1(i1,i2)=interp2d%interpolate_derivative_eta2(x1,x2)/transformation%jacobian(x1,x2)
        A2(i1,i2)=-interp2d%interpolate_derivative_eta1(x1,x2)/transformation%jacobian(x1,x2)
      end do
    end do
    
    
    
  end subroutine compute_field_from_phi_2d_curvilinear_mudpack
  
  subroutine time_history_diagnostic_gc_cartesian( &
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
    type(sll_logical_mesh_2d), pointer :: mesh_2d
    sll_real64, dimension(:,:), intent(in) :: f
    sll_real64, dimension(:,:), intent(in) :: phi
    sll_real64, dimension(:,:), intent(in) :: A1
    sll_real64, dimension(:,:), intent(in) :: A2 
    sll_real64 :: mass
    sll_real64 :: linf
    sll_real64 :: l1
    sll_real64 :: l2
    sll_real64 :: e
    sll_real64, dimension(:),allocatable :: data
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

    
    Nc_x1 = mesh_2d%num_cells1
    Nc_x2 = mesh_2d%num_cells2
    
    
    x1_min = mesh_2d%eta1_min
    x1_max = mesh_2d%eta1_max

    delta_x1 = mesh_2d%delta_eta1
    delta_x2 = mesh_2d%delta_eta2

    SLL_ALLOCATE(data(Nc_x1+1),ierr)
 
    linf  = 0.0_f64
    l1    = 0.0_f64
    l2    = 0.0_f64
    mass  = 0.0_f64
     e     = 0.0_f64
    
    do i2 = 1, Nc_x2+1
      do i1=1,Nc_x1+1
        data(i1) = f(i1,i2)
      enddo
      mass = mass + compute_integral_trapezoid_1d(data, Nc_x1+1, delta_x1)

      do i1=1,Nc_x1+1
        data(i1) = abs(f(i1,i2))
      enddo
      l1 = l1 + compute_integral_trapezoid_1d(data, Nc_x1+1, delta_x1)

      do i1=1,Nc_x1+1
        data(i1) = (f(i1,i2))**2
      enddo
      l2 = l2 + compute_integral_trapezoid_1d(data, Nc_x1+1, delta_x1)

      do i1=1,Nc_x1+1
        data(i1) = A2(i1,i2)**2+A1(i1,i2)**2
      enddo
      e = e + compute_integral_trapezoid_1d(data, Nc_x1+1, delta_x1)

      do i1=1,Nc_x1+1
       linf = max(linf,abs(f(i1,i2)))
      enddo
         
    enddo     

    mass = mass*delta_x2
    l1 = l1*delta_x2
    l2 = sqrt(l2*delta_x2)
    e  = e*delta_x2
    
    write(file_id,*) dt*real(step,f64),linf,l1,l2,mass,e


    
    
  end subroutine time_history_diagnostic_gc_cartesian


#ifndef NOHDF5
!*********************
!*********************

  !---------------------------------------------------
  ! Save the mesh structure
  !---------------------------------------------------
  subroutine plot_f_curvilinear(iplot,f,mesh_2d,transf)
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
    class(sll_coordinate_transformation_2d_base), pointer :: transf
    sll_real64, dimension(:,:), intent(in) :: f
    sll_real64 :: eta1
    sll_real64 :: eta2
    sll_real64 ::  eta1_min, eta2_min
    sll_real64 ::  eta1_max, eta2_max  
    sll_real64 :: deta1
    sll_real64 :: deta2
    
    
    nnodes_x1 = mesh_2d%num_cells1+1
    nnodes_x2 = mesh_2d%num_cells2+1
    eta1_min = mesh_2d%eta1_min
    eta1_max = mesh_2d%eta1_max
    eta2_min = mesh_2d%eta2_min
    eta2_max = mesh_2d%eta2_max
    deta1 = mesh_2d%delta_eta1
    deta2 = mesh_2d%delta_eta2
    
    !print *,'#maxf=',iplot,maxval(f),minval(f)
    

    
    if (iplot == 1) then

      SLL_ALLOCATE(x1(nnodes_x1,nnodes_x2), error)
      SLL_ALLOCATE(x2(nnodes_x1,nnodes_x2), error)
      do j = 1,nnodes_x2
        do i = 1,nnodes_x1
          eta1 = eta1_min+real(i-1,f32)*deta1
          eta2 = eta2_min+real(j-1,f32)*deta2
          x1(i,j) = transf%x1(eta1,eta2)
          x2(i,j) = transf%x2(eta1,eta2)
        end do
      end do
      call sll_hdf5_file_create("curvilinear_mesh-x1.h5",file_id,error)
      call sll_hdf5_write_array(file_id,x1,"/x1",error)
      call sll_hdf5_file_close(file_id, error)
      call sll_hdf5_file_create("curvilinear_mesh-x2.h5",file_id,error)
      call sll_hdf5_write_array(file_id,x2,"/x2",error)
      call sll_hdf5_file_close(file_id, error)
      deallocate(x1)
      deallocate(x2)

    end if

    call int2string(iplot,cplot)
    call sll_xdmf_open("f"//cplot//".xmf","curvilinear_mesh", &
      nnodes_x1,nnodes_x2,file_id,error)
    call sll_xdmf_write_array("f"//cplot,f,"values", &
      error,file_id,"Node")
    call sll_xdmf_close(file_id,error)
  end subroutine plot_f_curvilinear

#endif


end module sll_simulation_2d_guiding_center_curvilinear_mudpack_module
