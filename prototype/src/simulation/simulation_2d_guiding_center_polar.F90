module sll_simulation_2d_guiding_center_polar_module

!the aim is to translate CG_polar in simulation class
!related to
!simulation_guiding_center_2D_generalized_coords.F90
!but here geometry and test is specifically polar


#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_field_2d.h"
#include "sll_utilities.h"
#include "sll_poisson_solvers.h"
  use sll_logical_meshes  
  use sll_module_advection_1d_periodic
  use sll_module_advection_2d_BSL
  use sll_module_characteristics_2d_explicit_euler
  use sll_module_characteristics_2d_verlet
  !use sll_poisson_2d_periodic  
  use sll_fft
  use sll_reduction_module
  use sll_simulation_base
  use sll_cubic_spline_interpolator_2d
  use sll_cubic_spline_interpolator_1d
  use sll_coordinate_transformation_2d_base_module
  use sll_module_coordinate_transformations_2d
  use sll_common_coordinate_transformations
  use sll_common_array_initializers_module
  use sll_module_poisson_2d_polar_solver
  use sll_module_poisson_2d_elliptic_solver
  use sll_module_scalar_field_2d_base
  use sll_module_scalar_field_2d_alternative


  
  !use sll_parallel_array_initializer_module

  implicit none

!#define OLD_POISSON  
!#define NEW_POISSON  
  
  sll_int32, parameter :: SLL_EULER = 0 
  sll_int32, parameter :: SLL_PREDICTOR_CORRECTOR = 1 


  type, extends(sll_simulation_base_class) :: &
    sll_simulation_2d_guiding_center_polar

   !geometry
   type(sll_logical_mesh_2d), pointer :: mesh_2d


   !initial function
   procedure(sll_scalar_initializer_2d), nopass, pointer :: init_func
   sll_real64, dimension(:), pointer :: params
      
   !advector
   class(sll_advection_2d_base), pointer    :: advect_2d
   
   !interpolator for derivatives
   class(sll_interpolator_2d_base), pointer   :: phi_interp2d

   
   !poisson solver

!#ifdef NEW_POISSON
   class(sll_poisson_2d_base), pointer   :: poisson
!#endif

!#ifdef OLD_POISSON
!   !type(poisson_2d_periodic), pointer   :: poisson
!   type(sll_plan_poisson_polar), pointer :: poisson 
!#endif   
   !time_iterations
   sll_real64 :: dt
   sll_int32  :: num_iterations
   sll_int32  :: freq_diag
   sll_int32  :: freq_diag_time

   !time_loop
   sll_int32 :: time_loop_case
   
       
  contains
    procedure, pass(sim) :: run => run_gc2d_polar
    procedure, pass(sim) :: init_from_file => init_fake
     
  end type sll_simulation_2d_guiding_center_polar


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

  function new_guiding_center_2d_polar() result(sim)
    type(sll_simulation_2d_guiding_center_polar), pointer :: sim    
    sll_int32 :: ierr
    
    SLL_ALLOCATE(sim,ierr)
    
    call initialize_guiding_center_2d_polar(sim)
    
  
  
  end function new_guiding_center_2d_polar
  
  subroutine initialize_guiding_center_2d_polar(sim)
    class(sll_simulation_2d_guiding_center_polar), intent(inout) :: sim
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
    sll_real64 :: dt
    sll_int32 :: visu_step
    class(sll_interpolator_2d_base), pointer :: f_interp2d
    class(sll_interpolator_2d_base), pointer :: phi_interp2d
    class(sll_characteristics_2d_base), pointer :: charac2d
    class(sll_interpolator_2d_base), pointer   :: A1_interp2d
    class(sll_interpolator_2d_base), pointer   :: A2_interp2d
    class(sll_interpolator_1d_base), pointer   :: A1_interp1d_x1
    class(sll_interpolator_1d_base), pointer   :: A2_interp1d_x1
    sll_real64, dimension(:,:), pointer :: b11
    sll_real64, dimension(:,:), pointer :: b12
    sll_real64, dimension(:,:), pointer :: b21
    sll_real64, dimension(:,:), pointer :: b22
    sll_real64, dimension(:,:), pointer :: c
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
    
    !here we do all the initialization
    !in future, we will use namelist file

    x1_min = 1._f64
    x1_max = 10._f64
    x2_min = 0._f64
    x2_max = 2._f64*sll_pi
    Nc_x1 = 128
    Nc_x2 = 128
    r_minus = 4._f64
    r_plus = 5._f64
    k_mode = 3
    eps = 1.e-6_f64
    nb_step = 600
    
    dt = 0.1_f64
    !dt = 0.05_f64
    visu_step = 20
    f_interp2d_case = "SLL_CUBIC_SPLINES"
    phi_interp2d_case = "SLL_CUBIC_SPLINES"
    A_interp_case = "SLL_CUBIC_SPLINES"
    charac2d_case = "SLL_VERLET"
    !charac2d_case = "SLL_EULER"
    advect2d_case = "SLL_BSL"
    initial_function_case = "SLL_DIOCOTRON" 
    !time_loop_case = "SLL_EULER"
    time_loop_case = "SLL_PREDICTOR_CORRECTOR"    
    poisson_case = "POLAR_FFT"
    
    
    sim%dt = dt
    sim%num_iterations = nb_step
    sim%freq_diag = visu_step
    sim%freq_diag_time = 1

    sim%mesh_2d => new_logical_mesh_2d( &
      Nc_x1, &
      Nc_x2, &
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
          SLL_HERMITE, &
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
          SLL_HERMITE, &
          SLL_PERIODIC)
        A2_interp2d => new_cubic_spline_2d_interpolator( &
          Nc_x1+1, &
          Nc_x2+1, &
          x1_min, &
          x1_max, &
          x2_min, &
          x2_max, &
          SLL_HERMITE, &
          SLL_PERIODIC)  
        A1_interp1d_x1 => new_cubic_spline_1d_interpolator( &
          Nc_x1+1, &
          x1_min, &
          x1_max, &
          SLL_HERMITE)
        A2_interp1d_x1 => new_cubic_spline_1d_interpolator( &
          Nc_x1+1, &
          x1_min, &
          x1_max, &
          SLL_HERMITE)
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
          SLL_HERMITE, &
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
          bc_type_1=SLL_SET_TO_LIMIT, &
          bc_type_2=SLL_PERIODIC)    
      case ("SLL_VERLET")      
        charac2d => new_verlet_2d_charac(&
          Nc_x1+1, &
          Nc_x2+1, &
          A1_interp2d, &
          A2_interp2d, &
          A1_interp1d_x1, &
          A2_interp1d_x1, &
          bc_type_1=SLL_SET_TO_LIMIT, &
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
      case ("SLL_DIOCOTRON")
        sim%init_func => sll_diocotron_initializer_2d
        SLL_ALLOCATE(sim%params(4),ierr)
        sim%params(1) = r_minus
        sim%params(2) = r_plus
        sim%params(3) = eps
        sim%params(4) = k_mode
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
    
    
!#ifdef OLD_POISSON
!        sim%poisson => new_plan_poisson_polar( &
!      sim%mesh_2d%delta_eta1, &
!      x1_min, &
!      Nc_x1, &
!      Nc_x2, &
!      (/ SLL_NEUMANN_MODE_0,SLL_DIRICHLET/))
!#endif
    
!#ifdef NEW_POISSON    
    !poisson solver
    select case(poisson_case)    
      case ("POLAR_FFT")     
        sim%poisson =>new_poisson_2d_polar_solver( &
          x1_min, &
          x1_max, &
          Nc_x1, &
          Nc_x2, &
          (/SLL_NEUMANN_MODE_0, SLL_DIRICHLET/))
      case default
        print *,'#bad poisson_case',poisson_case
        print *,'#not implemented'
        print *,'#in initialize_guiding_center_2d_polar'
        stop
    end select
!#endif
    
  
  end subroutine initialize_guiding_center_2d_polar
  


  subroutine init_fake(sim, filename)
    class(sll_simulation_2d_guiding_center_polar), intent(inout) :: sim
    character(len=*), intent(in)                                :: filename
  
    print *,'# Do not use the routine init_vp4d_fake'
    print *,'#use instead initialize_vlasov_par_poisson_seq_cart'
    stop
  
  end subroutine init_fake
  
  subroutine run_gc2d_polar(sim)
    class(sll_simulation_2d_guiding_center_polar), intent(inout) :: sim
    sll_int32 :: Nc_x1
    sll_int32 :: Nc_x2
    sll_real64 :: delta_x1
    sll_real64 :: delta_x2
    sll_real64 :: x1_min
    sll_real64 :: x2_min    
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
    sll_int32 :: thdiag_id = 99 
    sll_int32             :: IO_stat
    sll_int32 :: iplot
    
    Nc_x1 = sim%mesh_2d%num_cells1
    Nc_x2 = sim%mesh_2d%num_cells2
    delta_x1 = sim%mesh_2d%delta_eta1
    delta_x2 = sim%mesh_2d%delta_eta2
    x1_min = sim%mesh_2d%eta1_min
    x2_min = sim%mesh_2d%eta2_min
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
      x2=x2_min+real(i2-1,f64)*delta_x2
      do i1=1,Nc_x1+1
        x1=x1_min+real(i1-1,f64)*delta_x1
        f(i1,i2) =  sim%init_func(x1,x2,sim%params)
      end do
    end do
        
    !solve poisson
!#ifdef OLD_POISSON
!    call poisson_solve_polar(sim%poisson,f,phi)
!#endif
!#ifdef NEW_POISSON
    call sim%poisson%compute_phi_from_rho( phi, f )
!#endif
    call compute_field_from_phi_2d_polar(phi,sim%mesh_2d,A1,A2,sim%phi_interp2d)
    
    
    !print *,A1
    !print *,A2
    
    open(unit = thdiag_id, file='thdiag.dat',IOStat=IO_stat)
    if( IO_stat /= 0 ) then
       print *, '#run_gc2d_polar(sim) failed to open file thdiag.dat'
       STOP
    end if
    
    iplot = 0

    do step=1,nb_step+1
      f_old = f
!#ifdef OLD_POISSON      
!      call poisson_solve_polar(sim%poisson,f_old,phi)
!#endif
!#ifdef NEW_POISSON
      call sim%poisson%compute_phi_from_rho( phi, f_old )
!#endif      
      call compute_field_from_phi_2d_polar(phi,sim%mesh_2d,A1,A2,sim%phi_interp2d)      
      
      if(modulo(step-1,sim%freq_diag_time)==0)then
        call time_history_diagnostic_gc_polar( &
          thdiag_id, &    
          step-1, &
          dt, &
          sim%mesh_2d, &
          f, &
          phi, &
          A1, &
          A2)
      endif            
      
      if(modulo(step-1,sim%freq_diag)==0)then
        call plot_f_polar(iplot,f,sim%mesh_2d)
        iplot = iplot+1  
      endif            
      
      select case (sim%time_loop_case)
        case (SLL_EULER)
          call sim%advect_2d%advect_2d(A1, A2, sim%dt, f_old, f)
        case (SLL_PREDICTOR_CORRECTOR)
          call sim%advect_2d%advect_2d(A1, A2, 0.5_f64*sim%dt, f_old, f)
!#ifdef NEW_POISSON
          call sim%poisson%compute_phi_from_rho( phi, f )
!#endif
!#ifdef OLD_POISSON
!          call poisson_solve_polar(sim%poisson,f,phi)
!#endif
          call compute_field_from_phi_2d_polar(phi,sim%mesh_2d,A1,A2,sim%phi_interp2d)      
          f_old = f
          call sim%advect_2d%advect_2d(A1, A2, 0.5_f64*sim%dt, f_old, f)
        case default  
          print *,'#bad time_loop_case',sim%time_loop_case
          print *,'#not implemented'
          print *,'#in run_gc2d_polar'
          print *,'#available options are:'
          print *,'#SLL_EULER=',SLL_EULER
          print *,'#SLL_PREDICTOR_CORRECTOR=',SLL_PREDICTOR_CORRECTOR
          
      end select
         
    enddo
    
    close(thdiag_id)

    !print *,'#not implemented for the moment!'
  end subroutine run_gc2d_polar    
  
  
  subroutine compute_field_from_phi_2d_polar(phi,mesh_2d,A1,A2,interp2d)
    sll_real64, dimension(:,:), intent(in) :: phi
    sll_real64, dimension(:,:), intent(out) :: A1
    sll_real64, dimension(:,:), intent(out) :: A2
    type(sll_logical_mesh_2d), pointer :: mesh_2d
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
        A1(i1,i2)=interp2d%interpolate_derivative_eta2(x1,x2)/x1
        A2(i1,i2)=-interp2d%interpolate_derivative_eta1(x1,x2)/x1
      end do
    end do
    
    
    
  end subroutine compute_field_from_phi_2d_polar
  
  subroutine time_history_diagnostic_gc_polar( &
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
    sll_real64 :: time_mode(8) 
    !sll_real64 :: mode_slope(8) 
    sll_real64 :: w
    sll_real64 :: l1
    sll_real64 :: l2
    sll_real64 :: e
    sll_real64, dimension(:),allocatable :: int_r
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
    type(sll_fft_plan), pointer         :: pfwd
    
    Nc_x1 = mesh_2d%num_cells1
    Nc_x2 = mesh_2d%num_cells2
    
    
    x1_min = mesh_2d%eta1_min
    x1_max = mesh_2d%eta1_max

    delta_x1 = mesh_2d%delta_eta1
    delta_x2 = mesh_2d%delta_eta2

    
    SLL_ALLOCATE(int_r(Nc_x2),ierr)
    SLL_ALLOCATE(data(Nc_x1+1),ierr)
    pfwd => fft_new_plan(Nc_x2,int_r,int_r,FFT_FORWARD,FFT_NORMALIZE)
 
    w     = 0.0_f64
    l1    = 0.0_f64
    l2    = 0.0_f64
    e     = 0.0_f64
    int_r = 0.0_f64
    
    
    do i2 = 1, Nc_x2
      do i1=1,Nc_x1+1
        x1 = x1_min+real(i1-1,f64)*delta_x1
        data(i1) = x1*f(i1,i2)
      enddo
      w = w + compute_integral_trapezoid_1d(data, Nc_x1+1, delta_x1)

      do i1=1,Nc_x1+1
        x1 = x1_min+real(i1-1,f64)*delta_x1
        data(i1) = x1*abs(f(i1,i2))
      enddo
      l1 = l1 + compute_integral_trapezoid_1d(data, Nc_x1+1, delta_x1)

      do i1=1,Nc_x1+1
        x1 = x1_min+real(i1-1,f64)*delta_x1
        data(i1) = x1*(f(i1,i2))**2
      enddo
      l2 = l2 + compute_integral_trapezoid_1d(data, Nc_x1+1, delta_x1)

      do i1=1,Nc_x1+1
        x1 = x1_min+real(i1-1,f64)*delta_x1
        data(i1) = x1*((x1*A2(i1,i2))**2+A1(i1,i2)**2)
      enddo
      e = e + compute_integral_trapezoid_1d(data, Nc_x1+1, delta_x1)

      do i1=1,Nc_x1+1
        x1 = x1_min+real(i1-1,f64)*delta_x1
        data(i1) = x1*phi(i1,i2)
      enddo
      int_r(i2) = compute_integral_trapezoid_1d(data, Nc_x1+1, delta_x1)      
    enddo     

    w = w*delta_x2
    l1 = l1*delta_x2
    l2 = sqrt(l2*delta_x2)
    e  = 0.5_f64*e*delta_x2
    call fft_apply_plan(pfwd,int_r,int_r)
    do i1=1,8
      !mode_slope(i1) = time_mode(i1)
      time_mode(i1) = abs(fft_get_mode(pfwd,int_r,i1-1))**2
      !mode_slope(i1) = &
      !  (log(0*time_mode(i1)+1.e-40_f64)-log(0*mode_slope(i1)+1.e-40_f64))/(dt+1.e-40_f64)
    enddo
    
    write(file_id,*) dt*real(step,f64),w,l1,l2,e,time_mode(1:8)!,mode_slope



    call fft_delete_plan(pfwd)
    
!    call fft_apply_plan(plan_sl%poisson%pfwd,int_r,int_r)
!
!
    
    
    
  end subroutine time_history_diagnostic_gc_polar


#ifndef NOHDF5
!*********************
!*********************

  !---------------------------------------------------
  ! Save the mesh structure
  !---------------------------------------------------
  subroutine plot_f_polar(iplot,f,mesh_2d)
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
    sll_real64 :: rmin
    sll_real64 :: rmax
    sll_real64 :: dr
    sll_real64 :: dtheta
    
    
    nnodes_x1 = mesh_2d%num_cells1+1
    nnodes_x2 = mesh_2d%num_cells2+1
    rmin = mesh_2d%eta1_min
    rmax = mesh_2d%eta1_max
    dr = mesh_2d%delta_eta1
    dtheta = mesh_2d%delta_eta2
    
    !print *,'#maxf=',iplot,maxval(f),minval(f)
    

    
    if (iplot == 1) then

      SLL_ALLOCATE(x1(nnodes_x1,nnodes_x2), error)
      SLL_ALLOCATE(x2(nnodes_x1,nnodes_x2), error)
      do j = 1,nnodes_x2
        do i = 1,nnodes_x1
          r       = rmin+real(i-1,f32)*dr
          theta   = real(j-1,f32)*dtheta
          x1(i,j) = r*cos(theta)
          x2(i,j) = r*sin(theta)
        end do
      end do
      call sll_hdf5_file_create("polar_mesh-x1.h5",file_id,error)
      call sll_hdf5_write_array(file_id,x1,"/x1",error)
      call sll_hdf5_file_close(file_id, error)
      call sll_hdf5_file_create("polar_mesh-x2.h5",file_id,error)
      call sll_hdf5_write_array(file_id,x2,"/x2",error)
      call sll_hdf5_file_close(file_id, error)
      deallocate(x1)
      deallocate(x2)

    end if

    call int2string(iplot,cplot)
    call sll_xdmf_open("f"//cplot//".xmf","polar_mesh", &
      nnodes_x1,nnodes_x2,file_id,error)
    call sll_xdmf_write_array("f"//cplot,f,"values", &
      error,file_id,"Node")
    call sll_xdmf_close(file_id,error)
  end subroutine plot_f_polar

#endif



end module sll_simulation_2d_guiding_center_polar_module
