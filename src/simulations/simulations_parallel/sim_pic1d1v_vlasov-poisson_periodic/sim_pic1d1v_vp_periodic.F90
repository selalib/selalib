module sll_module_simulation_pic1d1v_vp_periodic

#include "sll_working_precision.h"
#include "sll_errors.h"
#include "sll_memory.h"
#include "sll_assert.h"

  use sll_simulation_base, only: sll_simulation_base_class
  
  use sll_collective , only :       sll_world_collective , sll_collective_barrier ,&
     sll_boot_collective
  
  use sll_pic_1d_field_solver
  
  
      use pic_1d_particle_loading, only : sll_normal_rnd ,& 
         sll_initialize_intrinsic_mpi_random ,&   
         load_particle_species ,&
         sll_pic1d_ensure_boundary_conditions_species  ,&
         sll_pic1d_ensure_boundary_conditions ,&
         sll_pic1d_ensure_periodicity  ,& 
         sll_pic_1d_landaudamp_PDFxv ,&
         sll_local_maxwellian ,&
         control_variate_xv ,& 
         SLL_PIC1D_TESTCASE_IONBEAM, SLL_PIC1D_TESTCASE_LANDAU, pic1d_testcase  ,&
         num_species, landau_mode ,landau_alpha, enable_deltaf, enable_deltaf ,&
         SLL_PIC1D_TESTCASE_IONBEAM_ELECTRONS ,SLL_PIC1D_TESTCASE_QUIET ,&
         control_variate_v ,SLL_PIC1D_TESTCASE_BUMPONTAIL, set_loading_parameters
  
  implicit none
   
!==============================================================================
  
  ! Enumeration: particle pusher
  sll_int32, parameter :: SLL_PIC1D_PPUSHER_NONE       = 0
  sll_int32, parameter :: SLL_PIC1D_PPUSHER_EULER      = 1    !Explicit EULER
  sll_int32, parameter :: SLL_PIC1D_PPUSHER_VERLET     = 2
  sll_int32, parameter :: SLL_PIC1D_PPUSHER_RK2        = 3 !Runge Kutta 2
  sll_int32, parameter :: SLL_PIC1D_PPUSHER_RK4        = 4   !Runge Kutta 4
  sll_int32, parameter :: SLL_PIC1D_PPUSHER_HEUN       = 5
  sll_int32, parameter :: SLL_PIC1D_PPUSHER_SHIFT      = 6 !Shift by velocity
  sll_int32, parameter :: SLL_PIC1D_PPUSHER_LEAPFROG   = 7
  sll_int32, parameter :: SLL_PIC1D_PPUSHER_LEAPFROG_V = 8 !Variational Leapfrog
  sll_int32, parameter :: SLL_PIC1D_PPUSHER_RK3        = 9     ! Runge Kutta 3
  sll_int32, parameter :: SLL_PIC1D_PPUSHER_MERSON     = 10 !Merson 3-5 with integrated Err estimate

  !For parallelization MPI Rank and collective size
  sll_int32 :: coll_rank,coll_size
  sll_int32 :: ierr,i
   
  type, extends( sll_simulation_base_class ) :: &
      sll_simulation_pic1d1v_vp_periodic

    CHARACTER(LEN=256) :: path

    ! TODO: use new enumerations
    character(len=32) :: testcase    , ppusher    , psolver    , scenario
    sll_int32         :: testcase_int, ppusher_int, psolver_int, scenario_int

    sll_int32 ::  tsteps
    sll_int32 ::  sdeg
    sll_int32 ::  nmark
    sll_int32 ::  nstreams 
    sll_int32  ::  femp
    sll_int32 :: mesh_cells
    
    sll_real64 :: lalpha
    sll_real64 :: lmode
    sll_real64 :: tstepw

    sll_real64 :: interval_a
    sll_real64 :: interval_b

    logical :: pi_unit
    logical :: deltaf
    logical :: gnuplot_inline_output_user
    
    class( pic_1d_field_solver ), pointer :: fsolver
    sll_real64, allocatable               :: knots(:) 
    sll_real64, allocatable               :: electricpotential_interp(:) 
    sll_int32                             :: pushed_species
    sll_real64, allocatable               :: particle(:,:)
    procedure( sll_pic_1d_electric_field_external ), pointer :: Eex
     
  contains
  
    procedure :: run            => run_fake
    procedure :: init_from_file => init_fake
    procedure :: new_pic        => new_sll_pic_1d
    
  end type sll_simulation_pic1d1v_vp_periodic

  abstract interface
      function sll_pic_1d_electric_field_external( sim, x, t ) result(E)
          use sll_working_precision
          import sll_simulation_pic1d1v_vp_periodic
          class( sll_simulation_pic1d1v_vp_periodic ), intent( in ) :: sim  !< Simulation obj.
          sll_real64                                 , intent( in ) :: x(:) !< Position
          sll_real64                                 , intent( in ) :: t    !< Time
          sll_real64                                                :: E( size(x) )
      end function
  end interface
    
  interface sll_delete
     module procedure delete_pid1d_vp_periodic
  end interface sll_delete

!==============================================================================
contains
!==============================================================================




  subroutine new_sll_pic_1d( sim, mesh_dx )
            
    class( sll_simulation_pic1d1v_vp_periodic ), intent(inout) :: sim
    sll_real64                                 , intent(inout) :: mesh_dx
   
      !  SLL_ASSERT( is_power_of_two( int( sim%mesh_cells,i64)))

        !particle_pusher=trim(particle_pusher_user)

        !############## PARALLEL ##################


        !really global variables
        coll_rank = sll_get_collective_rank( sll_world_collective )
        coll_size = sll_get_collective_size( sll_world_collective )
        call sll_collective_barrier(sll_world_collective)


        !The Mesh is needed on every Node, it is global
        !space set and initiation of the knots and electrocpotential
        
        mesh_dx = (sim%interval_b - sim%interval_a)/sim%mesh_cells
        SLL_CLEAR_ALLOCATE(sim%knots(1:sim%mesh_cells+1),ierr)
        SLL_CLEAR_ALLOCATE(sim%electricpotential_interp(1:sim%mesh_cells),ierr)


        !Knots coordinates
        sim%knots(1)=sim%interval_a
        do i=2,size(sim%knots)
            sim%knots(i)=  sim%knots(1) + (i-1)*mesh_dx*1.0_f64
        enddo
        sim%knots(size(sim%knots))=sim%interval_b
        SLL_ASSERT(sim%knots(size(sim%knots))==sim%interval_b)


        if (coll_rank==0) then
            print *, "Size of MPI-Collective: ", coll_size   !number of cores involved in calculus (mpirun argument)
        endif
!parallel configuration
        call sll_collective_barrier(sll_world_collective)
        call sll_initialize_intrinsic_mpi_random(sll_world_collective)
        call sll_collective_barrier(sll_world_collective)

        !Set up Quasineutral solver
        
        !fsolver is of class pic_1d_field_solver
        
        select case( sim%scenario )
        case("landau")
            pic1d_testcase = SLL_PIC1D_TESTCASE_LANDAU
        case("ionbeam")
            pic1d_testcase = SLL_PIC1D_TESTCASE_IONBEAM
            sim%interval_a=0
            sim%interval_b=200  !0.0022_f64 !20mm
        case("quiet")
            pic1d_testcase = SLL_PIC1D_TESTCASE_QUIET
        case("bump")
            pic1d_testcase = SLL_PIC1D_TESTCASE_BUMPONTAIL
            landau_alpha=0.001_f64
            landau_mode=0.5_f64
        end select
         
        select case(pic1d_testcase)
            case(SLL_PIC1D_TESTCASE_IONBEAM)
                sim%fsolver => new_pic_1d_field_solver( sim%interval_a, sim%interval_b,& 
                    sim%sdeg, sim%mesh_cells, sim%psolver_int, sll_world_collective, SLL_DIRICHLET)
            case default
                sim%fsolver=>new_pic_1d_field_solver(sim%interval_a, sim%interval_b, sim%sdeg, &
                    sim%mesh_cells, sim%psolver_int, sll_world_collective,SLL_PERIODIC)
        end select

        !Set pointer to external electric field
        sim%Eex => pic1d_Eex_zero

        !Check Marker distribution on Cores
        if ( coll_rank==0 .and. mod(sim%nmark, coll_size)/=0) then
            print *, "Number of Markers per core: ", sim%nmark/(coll_size*1.0_f64)
            print *, "Choose appropriate number of markers!"
            stop
        endif
   
    end subroutine new_sll_pic_1d
    
    
    function pic1d_Eex_zero( sim, x, t ) result( E )
      class( sll_simulation_pic1d1v_vp_periodic ), intent( in ) :: sim  !< Simulation obj.
      sll_real64                                 , intent( in ) :: x(:) !< Position
      sll_real64                                 , intent( in ) :: t    !< Time
      sll_real64                                                :: E( size(x) )
      E = 0.0_f64
    end function

    
!    
!        function new_pic_1d_field_solver(eta_min, eta_max, &
!            spline_degree, num_cells, poisson_solver_type, collective ,boundary_type ) &
!            result(qn_solver)
!        sll_int32, intent(in):: spline_degree
!        sll_int32, intent(in)::  poisson_solver_type
!        sll_int32, intent(in)::  num_cells
!        sll_real64, intent(in) :: eta_min, eta_max
!        sll_int32, intent(in) :: boundary_type
!        type(sll_collective_t), pointer , intent(in):: collective

!        class(pic_1d_field_solver), pointer :: qn_solver

!        SLL_ALLOCATE(qn_solver,ierr)

!        call   pic_1d_field_solver_initialize( qn_solver, eta_min, eta_max, &
!            spline_degree, num_cells, poisson_solver_type, collective,boundary_type  )
!    endfunction
    
    
    
    
    
    
    



  subroutine run_fake( sim )
    class( sll_simulation_pic1d1v_vp_periodic ), intent( inout ) :: sim
    character( len=64 ), parameter :: this_sub_name = "run_fake"
    SLL_WARNING( this_sub_name, "'run' method not implemented" )   
  end subroutine run_fake





  subroutine init_fake( sim, filename )
    class( sll_simulation_pic1d1v_vp_periodic ), intent( inout ) :: sim
    character( len=* )                         , intent( in    ) :: filename
    character( len=64 ), parameter :: this_sub_name = "init_fake"
!    SLL_WARNING( this_sub_name, "'init_from_file' method not implemented" )
     
    CHARACTER(LEN=256) :: path

    CHARACTER(LEN=32) ::  testcase
    CHARACTER(LEN=32) ::  ppusher
    CHARACTER(LEN=32) ::  psolver
    CHARACTER(LEN=32) :: scenario

    sll_int32 ::  tsteps
    sll_int32 ::  sdeg
    sll_int32 ::  nmark
    sll_int32 ::  nstreams 
    sll_int32  ::  femp

    sll_real64 :: lalpha
    sll_real64 :: lmode
    sll_real64 :: tstepw

    sll_real64 :: interval_a
    sll_real64 :: interval_b

    logical :: pi_unit
    logical :: deltaf
    logical :: gnuplot_inline_output_user
    sll_int32, parameter  :: input_file =99
    sll_int32             :: IO_stat
    
    ! General parameters
   
   
    namelist /params/ nmark, tstepw, ppusher, scenario, psolver,gnuplot_inline_output_user,deltaf,nstreams
    
    ! Numerical parameters

    namelist /numerical_params/sdeg,tstepw,femp
    
    ! Landau parameters
   

    namelist /landau_params/ lalpha,lmode,pi_unit,interval_a,interval_b
 
    
   call getarg(1,filename)
    open(unit=input_file, file=trim(filename), IOStat=IO_stat)
    if( IO_stat /= 0 ) then
        print *, 'init_file() failed to open file ', filename
        STOP
    end if
    read(input_file,landau_params)
    read(input_file,numerical_params)
    read(input_file,params)
    close(input_file) 
  
  
      sim%lalpha=lalpha
      sim%lmode=lmode
      sim%pi_unit=pi_unit
      sim%interval_a=interval_a
      sim%interval_b=interval_b
            
      sim%nmark=nmark
      sim%tstepw  =tstepw
      sim%ppusher=ppusher
      sim%scenario=scenario
      sim%psolver=psolver
      sim%gnuplot_inline_output_user=gnuplot_inline_output_user
      
      sim%deltaf=deltaf
      sim%nstreams=nstreams
    
      sim%sdeg=sdeg
      sim%tstepw=tstepw
      sim%femp=femp
  
      
  end subroutine init_fake




!  interface initialize
!     module procedure initialize_pid1d_vlasov_poisson_periodic
!  end interface initialize
!contains
! subroutine initialize_pid1d_vlasov_poisson_periodi( sim,lalpha,lmode,pi_unit,interval_a,interval_b,nmark,sdeg,tstepw,femp,tsteps,ppusher,scenario,psolver,gnuplot_inline_output_user,deltaf,nstreams)
!    type(ssl_simulation_pic_1d_sheath), intent(inout)     :: sim


!end subroutine initialize_pid1d_vlasov_poisson_periodic


  subroutine delete_pid1d_vp_periodic( sim )
    class( sll_simulation_pic1d1v_vp_periodic ), intent( inout ) :: sim
    character( len=64 ), parameter :: this_sub_name = &
      "delete_pid1d_vp_periodic"
    SLL_WARNING( this_sub_name, "'delete' method not implemented" )   
  end subroutine delete_pid1d_vp_periodic

  !===========================================================================
  ! Matching enumerations (clumsy!).  TODO: use new enumeration classes
  
  function match_enumeration( enum_name, enum_string ) result( enum_int )
    character( len=32 ), intent( in ) :: enum_name
    character( len=32 ), intent( in ) :: enum_string
    sll_int32                         :: enum_int
    
    select case( enum_name )
      !-----------------------------------------------------------------------
      case( "psolver" )
        select case( enum_string )
          case("fem");     enum_int = SLL_SOLVER_FEM
          case("fd");      enum_int = SLL_SOLVER_FD
          case("fourier"); enum_int = SLL_SOLVER_FOURIER
          case("spec");    enum_int = SLL_SOLVER_SPECTRAL
          case default;    enum_int = SLL_SOLVER_FEM
        end select
      !-----------------------------------------------------------------------
      case( "ppusher" )
        select case( enum_string )
          case( "rk4" );   enum_int = SLL_PIC1D_PPUSHER_RK4
          case("verlet");  enum_int = SLL_PIC1D_PPUSHER_VERLET
          case("euler");   enum_int = SLL_PIC1D_PPUSHER_EULER
          case("lfrog_v"); enum_int = SLL_PIC1D_PPUSHER_LEAPFROG_V
          case("lfrog");   enum_int = SLL_PIC1D_PPUSHER_LEAPFROG
          case("rk2");     enum_int = SLL_PIC1D_PPUSHER_RK2
          case("rk3");     enum_int = SLL_PIC1D_PPUSHER_RK2
          case("merson");  enum_int = SLL_PIC1D_PPUSHER_MERSON
          case("heun");    enum_int = SLL_PIC1D_PPUSHER_HEUN
          case("none");    enum_int = SLL_PIC1D_PPUSHER_NONE
          case("shift ");  enum_int = SLL_PIC1D_PPUSHER_SHIFT
          case default;    enum_int = SLL_PIC1D_PPUSHER_NONE
        end select
      !-----------------------------------------------------------------------
      case( "testcase" )
      !-----------------------------------------------------------------------
      case( "scenario" )
      !-----------------------------------------------------------------------
      case default
    end select
    

    
  end function match_enumeration

end module sll_module_simulation_pic1d1v_vp_periodic












