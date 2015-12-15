module sll_m_sim_pic_vp_1d1v_cart

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_errors.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_boundary_condition_descriptors, only: &
    sll_p_dirichlet, &
    sll_p_periodic

  use sll_m_collective, only: &
    sll_s_boot_collective, &
    sll_s_collective_barrier, &
    sll_s_collective_bcast_real64, &
    sll_o_collective_globalsum, &
    sll_s_collective_globalsum_array_comp64, &
    sll_s_collective_globalsum_array_real64, &
    sll_t_collective_t, &
    sll_f_get_collective_rank, &
    sll_f_get_collective_size, &
    sll_v_world_collective

  use sll_m_constants, only: &
    sll_p_pi

  use sll_m_particle_1d_description, only: &
    sll_t_particle_1d_group

  use sll_m_pic_1d_distribution, only: &
    sll_t_pic1d_eulerian_distribution

  use sll_m_pic_1d_field_solver, only: &
    sll_f_new_pic_1d_field_solver, &
    sll_t_pic_1d_field_solver

  use sll_m_pic_1d_particle_loading, only: &
    sll_f_control_variate_xv, &
    sll_v_enable_deltaf, &
    sll_s_load_particle_species, &
    sll_v_num_species, &
    sll_s_set_loading_parameters, &
    sll_s_initialize_intrinsic_mpi_random, &
    sll_s_pic1d_ensure_boundary_conditions, &
    sll_s_pic1d_ensure_boundary_conditions_species, &
    sll_f_pic1d_ensure_periodicity, &
    sll_p_pic1d_testcase_bumpontail, &
    sll_p_pic1d_testcase_ionbeam, &
    sll_p_pic1d_testcase_ionbeam_electrons, &
    sll_p_pic1d_testcase_landau, &
    sll_p_pic1d_testcase_quiet

  use sll_m_pic_postprocessing, only: &
    sll_s_det_landau_damping

  use sll_m_pic_visu, only: &
    sll_s_distribution_xdmf, &
    sll_s_electricpotential_gnuplot_inline, &
    sll_s_energies_electrostatic_gnuplot_inline, &
    sll_s_particles_center_gnuplot_inline

  use sll_m_sim_base, only: &
    sll_c_simulation_base_class

  use sll_m_timer, only: &
    sll_s_set_time_mark, &
    sll_f_time_elapsed_between, &
    sll_t_time_mark

  use sll_m_utilities, only: &
    sll_s_int2string, &
    sll_s_new_file_id

  implicit none

  public :: &
    sll_o_delete, &
    sll_t_simulation_pic1d1v_vp_periodic

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   
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
  sll_int32, parameter :: SLL_SOLVER_FEM = 1
  sll_int32, parameter :: SLL_SOLVER_FD = 2
  sll_int32, parameter :: SLL_SOLVER_SPECTRAL = 3
  sll_int32, parameter :: SLL_SOLVER_FOURIER = 4
  sll_int32 :: n_x,n_v
  !For parallelization MPI Rank and collective size
  sll_int32 :: coll_rank,coll_size
  sll_int32 :: ierr,i,j
   
  type, extends( sll_c_simulation_base_class ) :: &
      sll_t_simulation_pic1d1v_vp_periodic

    character(len=256) :: root_path

    
    character(len=32) ::  ppusher    , psolver    , scenario
    sll_int32         ::  ppusher_int, psolver_int, scenario_int

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

    logical :: deltaf
    logical :: gnuplot_inline_output
    sll_real64 ::  particle_qm !<mass to electron mass ratio, with the sign of the charge
    class( sll_t_pic_1d_field_solver ), pointer :: fsolver
    sll_real64, allocatable               :: knots(:) 
    sll_real64, allocatable               :: electricpotential_interp(:) 
    sll_int32                             :: pushed_species
    sll_real64, allocatable               :: particle(:,:)
    type(sll_t_particle_1d_group)           :: species(10)
    sll_real64 :: kineticenergy_offset, impulse_offset
    sll_real64 :: initial_fem_inhom , fem_inhom
    sll_int64  :: fastest_particle_idx
    sll_int32  :: frequency_of_printing ! number of timesteps needed to generate one output_file 
    sll_real64, allocatable:: eval_solution(:)
    sll_real64,  allocatable:: particleposition(:)
    sll_real64,  allocatable :: particlespeed(:)
    sll_real64,  allocatable:: steadyparticleposition(:) !Steady non-moving objects aka Ions
    sll_real64, dimension(:), allocatable :: fieldenergy, kineticenergy, impulse, &
    thermal_velocity_estimate, particleweight_mean ,particleweight_var, inhom_var ,&
    push_error_mean, push_error_var, total_momentum
     
     
     !error 
     sll_real64 ::  num_err_noise_max, num_err_noise_l2, num_err_seminorm
    
    
    procedure( sll_pic_1d_electric_field_external ), pointer :: Eex   
  contains
  
    procedure :: run                           => run_pic_1d
    procedure :: init_from_file                => init_from_file
    procedure :: new_pic                       => new_sll_pic_1d
    procedure :: error_estimation              => sll_pic_1d_initial_error_estimates  
    procedure :: pic1d_write_phasespace        => sll_pic1d_write_phasespace
    !ode integrators
    
    procedure :: pic_1d_rungekutta4_step_array => sll_pic_1d_rungekutta4_step_array
    procedure :: pic_1d_Verlet_scheme          => sll_pic_1d_Verlet_scheme 
    procedure :: pic_1d_leap_frog              => sll_pic_1d_leap_frog
    procedure :: pic_1d_explicit_euler         => sll_pic_1d_explicit_euler
    procedure :: pic_1d_variational_leap_frog  => sll_pic_1d_variational_leap_frog
    procedure :: pic_1d_rungekutta2            => sll_pic_1d_rungekutta2
    procedure :: pic_1d_rungekutta3            => sll_pic_1d_rungekutta3
    procedure :: pic_1d_heun                   => sll_pic_1d_heun                     
    procedure :: pic_1d_merson4                => sll_pic_1d_merson4 
    procedure :: pic_1d_solve_qn               => sll_pic_1d_solve_qn
    procedure :: pic1d_adjustweights           => sll_pic1d_adjustweights
    procedure :: pic1d_write_result            => sll_pic1d_write_result
    procedure :: pic_1d_calc_push_error        => sll_pic_1d_calc_push_error
    procedure :: get_total_momentum            => sll_total_momentum
  end type sll_t_simulation_pic1d1v_vp_periodic

  abstract interface
      function sll_pic_1d_electric_field_external( sim, x, t ) result(E)
          use sll_m_working_precision
          import sll_t_simulation_pic1d1v_vp_periodic
          class( sll_t_simulation_pic1d1v_vp_periodic ), intent( in ) :: sim  !< Simulation obj.
          sll_real64                                 , intent( in ) :: x(:) !< Position
          sll_real64                                 , intent( in ) :: t    !< Time
          sll_real64                                                :: E( size(x) )
      end function
  end interface
    
  interface sll_o_delete
     module procedure delete_pid1d_vp_periodic
  end interface sll_o_delete

!==============================================================================
contains
!==============================================================================




  subroutine new_sll_pic_1d( sim )
            
    class( sll_t_simulation_pic1d1v_vp_periodic ), intent(inout) :: sim
    sll_real64   :: mesh_dx
   
      !  SLL_ASSERT( sll_f_is_power_of_two( int( sim%mesh_cells,i64)))

        !particle_pusher=trim(particle_pusher_user)

        !############## PARALLEL ##################


        !really global variables
        coll_rank = sll_f_get_collective_rank( sll_v_world_collective )
        coll_size = sll_f_get_collective_size( sll_v_world_collective )

        !number of cores involved in calculus (mpirun argument)
        if (coll_rank==0) print *, "Size of MPI-Collective: ", coll_size

        call sll_s_collective_barrier(sll_v_world_collective)
  
        !The Mesh is needed on every Node, it is global
        !space set and initiation of the knots and electrocpotential
        sim%mesh_cells=2**sim%femp
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

        !parallel configuration
        call sll_s_collective_barrier(sll_v_world_collective)
        call sll_s_initialize_intrinsic_mpi_random(sll_v_world_collective)
        call sll_s_collective_barrier(sll_v_world_collective)

        !Set up Quasineutral solver
        !fsolver is of class sll_t_pic_1d_field_solver
        select case(sim%scenario_int)
            case(sll_p_pic1d_testcase_ionbeam)
                sim%fsolver => sll_f_new_pic_1d_field_solver( sim%interval_a, sim%interval_b,& 
                    sim%sdeg, sim%mesh_cells, sim%psolver_int, sll_v_world_collective, sll_p_dirichlet)
            case default
                sim%fsolver=>sll_f_new_pic_1d_field_solver(sim%interval_a, sim%interval_b, sim%sdeg, &
                    sim%mesh_cells, sim%psolver_int, sll_v_world_collective,sll_p_periodic)
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
      class( sll_t_simulation_pic1d1v_vp_periodic ), intent( in ) :: sim  !< Simulation obj.
      sll_real64                                 , intent( in ) :: x(:) !< Position
      sll_real64                                 , intent( in ) :: t    !< Time
      sll_real64                                                :: E( size(x) )
      E = 0.0_f64
    end function

    
!    
!        function sll_f_new_pic_1d_field_solver(eta_min, eta_max, &
!            spline_degree, num_cells, poisson_solver_type, collective ,boundary_type ) &
!            result(qn_solver)
!        sll_int32, intent(in):: spline_degree
!        sll_int32, intent(in)::  poisson_solver_type
!        sll_int32, intent(in)::  num_cells
!        sll_real64, intent(in) :: eta_min, eta_max
!        sll_int32, intent(in) :: boundary_type
!        type(sll_t_collective_t), pointer , intent(in):: collective

!        class(sll_t_pic_1d_field_solver), pointer :: qn_solver

!        SLL_ALLOCATE(qn_solver,ierr)

!        call   pic_1d_field_solver_initialize( qn_solver, eta_min, eta_max, &
!            spline_degree, num_cells, poisson_solver_type, collective,boundary_type  )
!    endfunction
    



  subroutine init_from_file( sim, filename )
    class( sll_t_simulation_pic1d1v_vp_periodic ), intent( inout ) :: sim
    character( len=* )                         , intent( in )    :: filename
    character( len=64 ), parameter :: this_sub_name = "init_fake"
!    SLL_WARNING( this_sub_name, "'init_from_file' method not implemented" )
    character( len=256 )                                         :: arg_filename
     
    !character(len=256) :: path

    !character(len=32) ::  testcase
    character(len=32) ::  ppusher
    character(len=32) ::  psolver
    character(len=32) :: scenario

    sll_int32 ::  tsteps
    sll_int32 ::  sdeg
    sll_int32 ::  nmark
    sll_int32 ::  nstreams 
    sll_int32 ::  femp
    sll_int32 ::  frequency_of_printing
    
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
    namelist /params/ scenario, gnuplot_inline_output_user,  nstreams,frequency_of_printing,n_x,n_v
    
    ! Numerical parameters
    namelist /numerical_params/ nmark, ppusher, psolver, tstepw, tsteps, sdeg,&
      femp,deltaf
    
    ! Landau parameters
    namelist /landau_params/ lalpha, lmode, pi_unit, interval_a, interval_b
 
    ! Read namelists from input file
    call get_command_argument(1,arg_filename)
    open(unit=input_file, file=trim(arg_filename), IOStat=IO_stat)
    if( IO_stat /= 0 ) then
        print *, 'init_file() failed to open file ', arg_filename
        STOP
    end if
    read(input_file,landau_params)
    read(input_file,numerical_params)
    read(input_file,params)
    close(input_file) 
  
    ! Define global variables in "sll_m_pic_1d_particle_loading" module
    ! (needed for loading initial distribution function)
    call sll_s_set_loading_parameters( lalpha, lmode, nstreams )

    ! Store various parameters in PIC simulation object
    sim%lalpha = lalpha
    sim%lmode  = lmode
    if ( pi_unit ) then
      sim%interval_a = interval_a*sll_p_pi
      sim%interval_b = interval_b*sll_p_pi
    else
      sim%interval_a = interval_a
      sim%interval_b = interval_b
    endif
    sim%nmark    = nmark
    sim%tstepw   = tstepw
    sim%ppusher  = ppusher
    sim%scenario = scenario
    sim%psolver  = psolver
    sim%gnuplot_inline_output=gnuplot_inline_output_user
    
    sim%deltaf=deltaf
    sim%nstreams=nstreams
    sll_v_enable_deltaf=deltaf
    sim%sdeg=sdeg
    sim%tsteps=tsteps
    sim%femp=femp
    sim%frequency_of_printing=frequency_of_printing
    sim%ppusher_int= match_enumeration('ppusher' ,sim%ppusher) 
    sim%psolver_int= match_enumeration('psolver' ,sim%psolver)
    sim%scenario_int= match_enumeration('scenario' ,sim%scenario)
    
  end subroutine init_from_file

!  interface initialize
!     module procedure initialize_pid1d_vlasov_poisson_periodic
!  end interface initialize
!contains
! subroutine initialize_pid1d_vlasov_poisson_periodi( sim,lalpha,lmode,pi_unit,interval_a,interval_b,nmark,sdeg,tstepw,femp,tsteps,ppusher,scenario,psolver,gnuplot_inline_output_user,deltaf,nstreams)
!    type(ssl_simulation_pic_1d_sheath), intent(inout)     :: sim


!end subroutine initialize_pid1d_vlasov_poisson_periodic

  !============================================================================

  subroutine delete_pid1d_vp_periodic( sim )
    type( sll_t_simulation_pic1d1v_vp_periodic ), intent( inout ) :: sim
    character( len=64 ), parameter :: this_sub_name = &
      "delete_pid1d_vp_periodic"
    SLL_WARNING( this_sub_name, "'delete' method not implemented" )   
  end subroutine delete_pid1d_vp_periodic

  !============================================================================

  subroutine run_pic_1d( sim )
     class( sll_t_simulation_pic1d1v_vp_periodic ), intent( inout ) :: sim
     type( sll_t_pic1d_eulerian_distribution )  :: phase_space_distribution
     integer              :: timestep
     sll_int32            :: idx,jdx,kdx,nx,nv
     sll_real64           :: time
     type(sll_t_time_mark)  :: tstart, tstop
     logical              :: gnuplot_now
     sll_real64,   dimension(:), allocatable       ::   total_mass     
    
        sim%root_path=""
        sim%nmark=sim%nmark/coll_size
        print*, "#Core ", coll_rank, " handles particles", coll_rank*sim%nmark +1, "-", (coll_rank+1)*sim%nmark
        call sll_s_collective_barrier(sll_v_world_collective)
         
        if (coll_rank==0) print*, "#Total Number of particles: ", sim%nmark*coll_size
            
            !Scalar values to be recorded  
            !outputs : fieldenery kineticenergy  
            
            ! allocating output pointers 
            !global outputs
            SLL_CLEAR_ALLOCATE(sim%fieldenergy(1:sim%tsteps+1), ierr)
            SLL_CLEAR_ALLOCATE(sim%kineticenergy(1:+sim%tsteps+1), ierr)
            SLL_CLEAR_ALLOCATE(sim%impulse(1:sim%tsteps+1), ierr)
            SLL_CLEAR_ALLOCATE(sim%thermal_velocity_estimate(1:sim%tsteps+1), ierr)
            SLL_CLEAR_ALLOCATE(sim%particleweight_mean(1:sim%tsteps+1), ierr)
            SLL_CLEAR_ALLOCATE(sim%particleweight_var(1:sim%tsteps+1), ierr)
            SLL_CLEAR_ALLOCATE(sim%push_error_mean(1:sim%tsteps+1), ierr)
            SLL_CLEAR_ALLOCATE(sim%push_error_var(1:sim%tsteps+1), ierr)
            SLL_CLEAR_ALLOCATE(sim%inhom_var(1:sim%tsteps+1), ierr)
            SLL_CLEAR_ALLOCATE(sim%total_momentum(1:sim%tsteps+1), ierr)
            !particle caracterized  outputs
            SLL_CLEAR_ALLOCATE( sim%steadyparticleposition(1:sim%nmark),ierr)
            SLL_CLEAR_ALLOCATE(sim%particleposition(1:sim%nmark),ierr)
            SLL_CLEAR_ALLOCATE(sim%particlespeed(1:sim%nmark),ierr)
            SLL_CLEAR_ALLOCATE(sim%eval_solution(1:size(sim%knots)),ierr)
            SLL_CLEAR_ALLOCATE(total_mass(1:sim%tsteps+1), ierr)
            
            if (coll_rank==0) write(*,*) "#PIC1D: Loading particles..."
            
            call sim%fsolver%set_num_sample(sim%nmark*coll_size)
            call sll_s_load_particle_species (sim%nmark, sim%interval_a, sim%interval_b ,&
                                        sim%species, sim%scenario_int )
            sim%particle_qm = sim%species(1)%qm
            
            ! Wait here for all processors to be done with the instructions above
            call sll_s_collective_barrier( sll_v_world_collective )
            
            if (coll_rank==0) write(*,*) "#PIC1D: Particles loaded..."
               
            !Determine index of fastest particle
            sim%fastest_particle_idx = int(maxloc( sim%particlespeed, dim=1 ),8)

            !Do precalculations, especially for the FEM solver  
            !call sll_bspline_fem_solver_1d_initialize(knots, spline_degree, steadyparticleposition, -(sll_p_charge/sll_p_epsilon_0) )
           
            if (coll_rank==0) print *, "#Number of FEM-Cells ", sim%mesh_cells

            
            ! Add ion background density to Poisson solver, if needed
            ! TODO: ask Jakob about exact behavior here...
            
            
            select case( sim%scenario_int )
              case( sll_p_pic1d_testcase_quiet )
                call sim%fsolver%set_ions_constant(0.0_f64)
                sim%particle_qm = -1.0_f64
              case( sll_p_pic1d_testcase_landau )
                call sim%fsolver%set_ions_constant(0.0_f64)
              case( sll_p_pic1d_testcase_ionbeam )
                call sim%fsolver%set_ions_constant(0.0_f64)
              case( sll_p_pic1d_testcase_ionbeam_electrons )
                call sim%fsolver%set_ions_constant_particles(sim%steadyparticleposition)
              case default
                sim%particle_qm = -1.0_f64 !By default we simulate Electrons
                call sim%fsolver%set_ions_constant(0.0_f64)
            end select

            ! Control variate method
            if (sim%deltaf .eqv. .TRUE.) then
                  !The Control Variate, as it is constant over time gos into the constant inhomogenity
                  !In the constant electron case this should be normalized to zero
                  !But since the fem solver ignores the constant ion background anyway
                  !this is going to be ok.
                  !Also this should be done analytically and not like that as an Monte Carlo estimate
                  !Here suppose the control variate is the initial state by default
                  call  sim%fsolver%set_bg_particles(sim%species(1)%particle%dx, &
                        sign(1.0_f64, sim%species(1)%qm)*sim%species(1)%particle%weight_const)

                sim%kineticenergy_offset=sll_pic1d_calc_kineticenergy_offset(sim%species(1:sll_v_num_species))
                sim%impulse_offset=sll_pic1d_calc_impulse_offset(sim%species(1:sll_v_num_species))
            else
                !Here should be a quick deactivate since every code is primarily deltaf
                sim%kineticenergy_offset=0.0_f64
                sim%impulse_offset=0.0_f64
            endif

           !Start with PIC Method
           !Do the initial solve of the field
           call sll_s_pic1d_ensure_boundary_conditions_species(sim%species(1:sll_v_num_species) , sim%scenario_int  )


           !Initial Solve
           call sim%fsolver%reset_particles()
           call sim%fsolver%set_species(sim%species(1:sll_v_num_species))
           call sim%fsolver%solve()
           print *, "#Initial field solve"
           print *, "#Test initial field, and loading condition for landau damping"
           call sim%fsolver%evalE(sim%knots(1:sim%mesh_cells)+ (sll_p_pi/6)*(sim%knots(2)-sim%knots(1)), sim%eval_solution(1:sim%mesh_cells))
           
           call  sim%error_estimation()        
              ! modification veryfing 
  !Information on initial field

          print *, "#Initial Field average: " , sum( sim%fsolver%getPotential())
          print *, "#Initial Field variance: " , &
          sum((sum(sim%fsolver%getPotential())/real(sim%mesh_cells,f64) &
          -  sim%fsolver%getPotential())**2)/real(sim%mesh_cells,f64)


        !---------------------------Initial Plots----------------------------------------
        if (coll_rank==0) then
           !  open(unit=20, file="initial_field.txt")
        open(unit=20, file=trim(ADJUSTL(sim%root_path))//"initial_field.txt")
          write (20,*) sim%eval_solution(1:sim%mesh_cells)
          close(20)
          print *, "#Field Energy                  " , &
                  "Kinetic Energy                  ",    "Total energy          ","Total mass  ",  " Total momentum "
        endif
        
        !Write initial phasespace
        do jdx=0, coll_size
            call sll_s_collective_barrier(sll_v_world_collective)
            if( coll_rank==jdx) then
                if (jdx==0) then
                    open(unit=20, file=trim(ADJUSTL(sim%root_path))//"initial_phasespace.dat")
                else
                    open(unit=20, file=trim(sim%root_path)//"initial_phasespace.dat",position="append")
                endif
                do kdx=1, sll_v_num_species
                    do idx=1, size(sim%species(kdx)%particle)
                        write (20,*) sim%species(kdx)%particle(idx)%dx, sim%species(kdx)%particle(idx)%vx, sim%species(kdx)%particle(idx)%weight
                    enddo
                enddo
                close(20)
            endif
        enddo
        
        ! TODO: find out meaning of this
        call sll_s_set_time_mark( tstart )
        
!---------------------------### Main Loop
        time = 0.0_f64
!        nx=floor(sqrt(sim%nmark/10.0_f64))+1
!        nv=nx
        nx=n_x
        nv=n_v 
        call phase_space_distribution%initialize(sim%interval_a,sim%interval_b,nx ,& 
        minval(sim%species(1)%particle%vx), maxval(sim%species(1)%particle%vx),nv)  
        
        do timestep=1, sim%tsteps
        
            ! Compute time instant
            time = sim%tstepw*(timestep-1) ! WARNING: bug was here
            
            ! Determine if gnuplot should create plots
            gnuplot_now = (sim%nmark*coll_size<=10000 .AND. (sim%gnuplot_inline_output .eqv. .true.) .AND. (coll_rank==0) .AND. mod(timestep-1,sim%tsteps/100)==0)
            
            ! Write currrent phase-space to file
           
!            call sim%pic1d_write_phasespace( timestep )

            
            if (mod(timestep-1,sim%frequency_of_printing)==0) then
                            
                 call phase_space_distribution%compute_f(sim%species(1)%particle%dx, sim%species(1)%particle%vx, &
                sim%species(1)%particle%weight)
                 call phase_space_distribution%print_f('distribution_',timestep)
                 call phase_space_distribution%compute_moments()
                 call phase_space_distribution%print_moments('moment_output',timestep,sim%root_path)
           end if  
            ! Gnuplot real-time plots
            if (gnuplot_now) then
                call sll_s_particles_center_gnuplot_inline(pic_1d_allparticlepos(sim%species(1:sll_v_num_species),sim%nmark), &
                    pic_1d_allparticlev(sim%species(1:sll_v_num_species),sim%nmark), &
                    sim%interval_a, sim%interval_b,&
                    -1.1_f64*maxval(pic_1d_allparticlev(sim%species(1:sll_v_num_species),sim%nmark)),&
                    1.1_f64*maxval( pic_1d_allparticlev(sim%species(1:sll_v_num_species),sim%nmark)),&
                    time)
            end if
            
            ! Compute mean value and variance 
            sim%particleweight_mean(timestep)       = sum(sim%species(1)%particle%weight) &
               / real(size(sim%species(1)%particle),f64) ! WARNING: bug here
            
            sim%particleweight_var(timestep)        = sum(sim%species(1)%particle%weight**2) &
               / real(size(sim%species(1)%particle),f64)          &
                                                    - sim%particleweight_mean(timestep)**2

            ! Other diagnostic quantities: kin. energy, el. energy, momentum, etc..
            sim%kineticenergy(timestep)             = sll_pic1d_calc_kineticenergy( sim%species(1:sll_v_num_species),sim%deltaf )
            sim%fieldenergy(timestep)               = sll_pic1d_calc_fieldenergy(sim%species(1:sll_v_num_species),sim)
            sim%impulse(timestep)                   = sll_pic1d_calc_impulse(sim%species(1:sll_v_num_species),sim%deltaf)
            sim%total_momentum(timestep)            = sim%get_total_momentum()
            sim%thermal_velocity_estimate(timestep) = sll_pic1d_calc_thermal_velocity(sim%species(1)%particle%vx,sim%species(1)%particle%weight)
            sim%inhom_var(timestep)                 = sim%fsolver%calc_variance_rhs()
            total_mass(timestep)                    = calculate_total_mass(sim)
            ! Gnuplot real-time plots
            if ( gnuplot_now ) then
                call sll_s_energies_electrostatic_gnuplot_inline(sim%kineticenergy(1:timestep), sim%fieldenergy(1:timestep), sim%impulse(1:timestep),sim%tstepw)
                call sim%fsolver%evalPhi(sim%knots(1:sim%mesh_cells), sim%electricpotential_interp)
                call sll_s_electricpotential_gnuplot_inline( sim%electricpotential_interp  &
                    !+0.5_f64*knots(1:mesh_cells)&
                    ,sim%knots(1:sim%mesh_cells) )
            endif
            
            ! Stop if Error is catastrophic
            if ( coll_rank==0) then
                if ( (sim%kineticenergy(timestep)+sim%fieldenergy(timestep) -(sim%kineticenergy(1)+sim%fieldenergy(1)))&
                        /(sim%kineticenergy(1)+sim%fieldenergy(1)) > 10000.0_f64) then
                    !                    print *, "Its Over, Giving Up!"
                    !                    stop
                endif
            endif
             
             
            ! Print electrostatic and kinetic energy to terminal
            if ((sim%gnuplot_inline_output .eqv. .FALSE.) .AND. coll_rank==0) then
                print *, timestep, sim%fieldenergy(timestep),sim%kineticenergy(timestep),sim%kineticenergy(timestep)+sim%fieldenergy(timestep),&
                total_mass(timestep),sim%total_momentum(timestep)
            endif
            
            ! TODO: discover meaning of this line!
            sim%pushed_species=1
            
            ! Evolve solution from t to t+dt (i.e. push particles)
            select case (sim%ppusher_int)
            
                case(SLL_PIC1D_PPUSHER_RK4)         ;call sim%pic_1d_rungekutta4_step_array(time)                  
                case(SLL_PIC1D_PPUSHER_VERLET)      ;call sim%pic_1d_Verlet_scheme(time)
                case(SLL_PIC1D_PPUSHER_EULER)       ;call sim%pic_1d_explicit_euler(time)
                case(SLL_PIC1D_PPUSHER_LEAPFROG_V)  ;call sim%pic_1d_variational_leap_frog()    
                case(SLL_PIC1D_PPUSHER_LEAPFROG)    ;call sim%pic_1d_leap_frog()
                case(SLL_PIC1D_PPUSHER_RK2)         ;call sim%pic_1d_rungekutta2(time)
                case(SLL_PIC1D_PPUSHER_RK3)         ;call sim%pic_1d_rungekutta3(time)
                case(SLL_PIC1D_PPUSHER_HEUN)        ;call sim%pic_1d_heun(time)
                case(SLL_PIC1D_PPUSHER_MERSON)      ;call sim%pic_1d_merson4(time)
                case(SLL_PIC1D_PPUSHER_NONE)
                    print *, "No Particle pusher choosen"
                    !                case(SLL_PIC1D_PPUSHER_SHIFT)
                    !                    particleposition=particleposition+timestepwidth*(interval_b-interval_a)
                    !                    !particleposition=sll_f_pic1d_ensure_periodicity( particleposition, &
                        !                     !   interval_a, interval_b)
                    !                     call sll_s_pic1d_ensure_boundary_conditions(particleposition, particlespeed)
                    !
                    !                    !call sll_bspline_fem_solver_1d_solve(particleposition)
                    !                    call sll_pic_1d_solve_qn(particleposition)
            end select
     
            ! Estimate total simulation runtime
            if (coll_rank==0 .AND. timestep==1) then
                call sll_s_set_time_mark(tstop)
                print *, "Remaining Time: " , &
       (sll_f_time_elapsed_between(tstart,tstop)/real(timestep,f64))*real(sim%tsteps-timestep,i64)
            endif
            phase_space_distribution%vmin = minval(sim%species(1)%particle%vx)
            phase_space_distribution%vmax = maxval(sim%species(1)%particle%vx)
   !         print *,"min and max values of velocity:",phase_space_distribution%vmin,phase_space_distribution%vmax

        enddo

        ! Print final value of various diagnostic quantities (in dedicated arrays)
        sim%kineticenergy(timestep)   = sll_pic1d_calc_kineticenergy( sim%species(1:sll_v_num_species),sim%deltaf )
        sim%fieldenergy(timestep)     = sll_pic1d_calc_fieldenergy(sim%species(1:sll_v_num_species),sim)
        sim%impulse(timestep)         = sll_pic1d_calc_impulse(sim%species(1:sll_v_num_species),sim%deltaf)
        sim%thermal_velocity_estimate(timestep) = sll_pic1d_calc_thermal_velocity(sim%species(1)%particle%vx,sim%species(1)%particle%weight)
        sim%particleweight_mean(timestep) = sum(sim%species(1)%particle%weight)&
          /real(size(sim%species(1)%particle),f64)
        sim%particleweight_var(timestep)  = sum(sim%species(1)%particle%weight**2)&
         /real(size(sim%species(1)%particle),f64)-sim%particleweight_mean(timestep)**2

        close(25)

        ! Save Results to file
        call sim%pic1d_write_result("pic1dresult")

        ! Print total simulation time
        if (coll_rank==0) then
            call sll_s_set_time_mark(tstop)
            print *, "Overall Time: " , sll_f_time_elapsed_between(tstart,tstop)
            print *, "Particles Pushed/s: " , real(sim%nmark*coll_size*sim%tsteps,f64)/sll_f_time_elapsed_between(tstart,tstop)
        endif
        
        ! Compute damping (exponential factor) of electrostatic energy (this makes sense only for Landau damping)
        if (coll_rank==0) then
         call sll_s_det_landau_damping( (/ ( timestep*sim%tstepw, timestep = 0, sim%tsteps) /), sim%fieldenergy )
        endif
        
        ! Deallocate field solver
        !call sll_bspline_fem_solver_1d_destroy
        call sim%fsolver%delete()

  end subroutine run_pic_1d

  
  !===========================================================================
  ! Matching enumerations (clumsy!).  TODO: use new enumeration classes
  
  function match_enumeration( enum_name, enum_string ) result( enum_int )
    character( len=* ), intent( in ) :: enum_name
    character( len=* ), intent( in ) :: enum_string
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
      case( "scenario" )
        select case( enum_string )
          case("landau");  enum_int = sll_p_pic1d_testcase_landau
          case("ionbeam"); enum_int = sll_p_pic1d_testcase_ionbeam
          case("quiet");   enum_int = sll_p_pic1d_testcase_quiet
          case("bump");    enum_int = sll_p_pic1d_testcase_bumpontail
        end select

    end select
    

    
  end function match_enumeration







function sll_pic1d_calc_kineticenergy_offset(p_species ) &
            result(energy)
        type( sll_t_particle_1d_group), dimension(:), intent(in) :: p_species
        sll_int32 :: idx
        sll_int32 :: num_sp
        sll_real64 :: energy
        num_sp=size(p_species)
        energy=0.0_f64
        do idx=1,num_sp
            energy=energy + &
            sll_pic1d_calc_kineticenergy_weighted(p_species(idx)%particle%vx, &
                            p_species(idx)%particle%weight_const,1.0_f64/abs(p_species(idx)%qm))
        enddo

        call sll_o_collective_globalsum(sll_v_world_collective, energy, 0)
    end function



    function sll_pic1d_calc_impulse(p_species,bool) &
            result(impulse)
        type( sll_t_particle_1d_group), dimension(:), intent(in) :: p_species
        LOGICAL, intent(in) :: bool
        sll_real64:: impulse
        sll_int32 :: idx, num_sp
        num_sp=size(p_species)

        impulse=0.0_f64
        do idx=1,num_sp
            impulse=impulse + sll_pic1d_calc_impulse_weighted&
               (p_species(idx)%particle%weight, p_species(idx)%particle%vx, &
                            1.0_f64/abs(p_species(idx)%qm))
        enddo
        call sll_o_collective_globalsum(sll_v_world_collective, impulse, 0)
        if(bool .eqv. .TRUE.) then
        impulse=impulse+sll_pic1d_calc_impulse_offset(p_species) !Warning: bug here...
        endif
          
    end function

 function sll_pic1d_calc_impulse_weighted(particle_v, particle_weight, particlemass) &
            result(impulse)
                sll_real64, dimension(:),intent(in) :: particle_v,particle_weight
        sll_real64, intent(in) ::particlemass
        sll_real64 :: impulse
        SLL_ASSERT(size(particle_v)==size(particle_weight))

        impulse= 0.5_f64*particlemass*&
                dot_product(particle_weight, particle_v)
    end function


  function sll_pic1d_calc_impulse_offset(p_species) &
            result(impulse)
        type( sll_t_particle_1d_group), dimension(:), intent(in) :: p_species
        sll_real64:: impulse
        sll_int32 :: idx, num_sp
        num_sp=size(p_species)

        impulse=0.0_f64
        do idx=1,num_sp
            impulse=impulse + sll_pic1d_calc_impulse_weighted&
               (p_species(idx)%particle%weight_const, p_species(idx)%particle%vx, &
                            1.0_f64/abs(p_species(idx)%qm))
        enddo

        call sll_o_collective_globalsum(sll_v_world_collective, impulse, 0)
    end function



    function pic_1d_allparticlepos(p_species,npartic) result( ppos)
        type( sll_t_particle_1d_group), dimension(:), intent(in) :: p_species
        sll_int32, intent(inout) :: npartic
        sll_int32 :: idx, num, npart
        sll_int32 :: off
        sll_real64, dimension(npartic) :: ppos
        num=size(p_species)

        off=0

        do idx=1,num
            npart=size(p_species(idx)%particle)
            ppos(1+off:npart+off)=p_species(idx)%particle%dx
            off=npart
        enddo
    end function


    function sll_pic1d_calc_kineticenergy_weighted( particle_v, particle_weight, particlemass ) &
            result(energy)
        sll_real64, dimension(:),intent(in) :: particle_v,particle_weight
        sll_real64, intent(in) ::particlemass
        sll_real64 :: energy
        SLL_ASSERT(size(particle_v)==size(particle_weight))

        energy= 0.5_f64*particlemass*&
                dot_product(particle_weight, particle_v**2)
    end function


    function pic_1d_allparticlev(p_species,npartic) result( pv)
        type( sll_t_particle_1d_group), dimension(:), intent(in) :: p_species
        sll_int32, intent(inout) :: npartic
        sll_int32 :: idx, num, npart
        sll_int32 :: off
        sll_real64, dimension(npartic) :: pv
        num=size(p_species)

        off=0

        do idx=1,num
            npart=size(p_species(idx)%particle)
            pv(1+off:npart+off)=p_species(idx)%particle%vx
            off=npart
        enddo
    end function
    
    function calculate_total_mass(sim) result (total_massk)
       !integer   ::  i,k 
       class(sll_t_simulation_pic1d1v_vp_periodic), intent(inout) :: sim
       sll_real64                                              :: total_massk

     
       
        total_massk=sum((sim%species(1)%particle%weight))
       
          
       
  end function     
    
    
       
   function sll_total_momentum(sim) result (momentum)   
       class(sll_t_simulation_pic1d1v_vp_periodic), intent(inout) :: sim
       sll_real64   :: momentum
       !integer      :: K
       
       momentum=dot_product(sim%species(1)%particle%vx,sim%species(1)%particle%weight) 
       
  end function     
       
       
      
       
       
       
       
    
        function sll_pic1d_calc_kineticenergy(p_species, bool) &
            result(energy)
        type( sll_t_particle_1d_group), dimension(:), intent(in) :: p_species
        LOGICAL, intent(in) :: bool
        sll_int32 :: idx
        sll_int32 :: num_sp
        sll_real64 :: energy
        num_sp=size(p_species)
        energy=0.0_f64
       do idx=1,num_sp
       
            energy=energy + &
            sll_pic1d_calc_kineticenergy_weighted(p_species(idx)%particle%vx, &
                            p_species(idx)%particle%weight,1.0_f64/abs(p_species(idx)%qm))
        enddo

        call sll_o_collective_globalsum(sll_v_world_collective, energy, 0)
        !energy=energy+kineticenergy_offset
        if(bool .eqv. .TRUE.) then 
        energy=energy+sll_pic1d_calc_kineticenergy_offset(p_species)
        endif
    end function




!warning : can this a method in the class ??

  function sll_pic1d_calc_fieldenergy(p_species,sim) &
            result(energy)
        type( sll_t_particle_1d_group), dimension(:), intent(in) :: p_species
        class(sll_t_simulation_pic1d1v_vp_periodic), intent(inout) :: sim
        sll_real64 :: energy
        !Only for constant case
        sll_int32 :: idx
        sll_int32 :: num_sp
        num_sp=size(p_species)

        energy=0.0_f64
        do idx=1,num_sp
            energy=energy+ &
                dot_product(sim%Eex( p_species(idx)%particle%dx, 0.0_f64), &
                p_species(idx)%particle%dx)*(-1.0_f64)*p_species(idx)%qm/abs(p_species(idx)%qm)
        enddo

        call sll_o_collective_globalsum(sll_v_world_collective, energy, 0)

        !if (coll_rank==0)  energy=energy+0.5_f64*bspline_fem_solver_1d_H1seminorm_solution()/(interval_b-interval_a)
        if (coll_rank==0)  energy=energy+sim%fsolver%fieldenergy()


    end function




  function sll_pic1d_calc_thermal_velocity(particlespeed, particleweight) &
            result(vth)
        sll_real64, DIMENSION(:), intent(in) :: particlespeed
        sll_real64, DIMENSION(:), intent(in):: particleweight
        sll_real64 :: vth
        sll_real64 :: mean
        mean=dot_product(particlespeed,particleweight)
        call sll_o_collective_globalsum(sll_v_world_collective, mean, 0)


        vth=dot_product((particlespeed-mean)**2,particleweight)
        call sll_o_collective_globalsum(sll_v_world_collective, vth, 0)
    end function






subroutine sll_pic_1d_initial_error_estimates(sim)

            class( sll_t_simulation_pic1d1v_vp_periodic ), intent( inout ) :: sim
!           sll_real64 ::  num_err_noise_max, num_err_noise_l2, num_err_seminorm
!           sll_real64 :: num_err_mean
            sll_real64, dimension(size(sim%knots)-1) :: analytical_solution
        
        
        selectcase(sim%scenario_int)
               case(sll_p_pic1d_testcase_landau)
                !The first field solve is a Monte Carlo estimator
                   if (coll_rank==0) then
                    !analytical_solution=landau_alpha/((interval_b-interval_a)*landau_mode)*&
                            if (sim%lalpha/=0) then
                        analytical_solution=sim%lalpha/(sim%lmode)*&
                            sin(sim%lmode*(sim%knots(1:sim%mesh_cells)+(sll_p_pi/6)*(sim%knots(2)-sim%knots(1))))
                        sim%num_err_noise_max= maxval(abs((sim%eval_solution(1:sim%mesh_cells)-analytical_solution)))/&
                            maxval(abs(analytical_solution))
                        sim%num_err_noise_l2=sum(((sim%eval_solution(1:sim%mesh_cells)-analytical_solution))**2)/&
                            sum(analytical_solution**2)
                    else
                        sim%num_err_noise_max= maxval(abs((sim%eval_solution(1:sim%mesh_cells))))
                    endif


                    print *, "Error in MC estimate of E-field (Max):", sim%num_err_noise_max, "(L2):" ,sim%num_err_noise_l2
                endif

                !num_err_seminorm=landau_alpha**2/(2*(interval_b-interval_a)*landau_mode**2)
                if (sim%lalpha/=0) then
                    sim%num_err_seminorm=sim%lalpha**2/(2*sim%lmode**2)/2.0_f64 !*(interval_b-interval_a)/2.0_f64
                    sim%num_err_seminorm=(sll_pic1d_calc_fieldenergy(sim%species(1:sll_v_num_species),sim)-sim%num_err_seminorm)/sim%num_err_seminorm
                else
                    sim%num_err_seminorm=sll_pic1d_calc_fieldenergy(sim%species(1:sll_v_num_species),sim)
                endif
                if (coll_rank==0) print *, "Error in MC estimate of E-Energy: (rel.)", sim%num_err_seminorm

        endselect

    end subroutine







    subroutine sll_pic1d_write_phasespace(sim, timestep)
        sll_int32, intent(in) :: timestep
        class( sll_t_simulation_pic1d1v_vp_periodic ), intent( inout ) :: sim

        sll_int32:: grid=100
        if (coll_rank==0) then
            !            if (phasespace_file_id==-1) then
            !                        call sll_s_new_file_id(phasespace_file_id, ierr)
            !                    open(phasespace_file_id, file = trim(root_path)//'phasespace.dat')
            !            endif
            !
            !            write(phasespace_file_id,*)  particleposition, particlespeed
            ! if (phasespace_file_id==-1) then
            SLL_ASSERT(    minval(sim%particleposition) >=sim%interval_a)
            SLL_ASSERT(    maxval(sim%particleposition) <=sim%interval_b)
            grid=floor(sqrt(sim%nmark/10.0_f64))+1
            call sll_s_distribution_xdmf(trim(sim%root_path)//"phasespace", sim%species(1)%particle%dx, sim%species(1)%particle%vx, &
                sim%species(1)%particle%weight,  &
                sim%interval_a-1.0_f64, sim%interval_b+1.0_f64, grid, &
                minval(sim%species(1)%particle%vx), maxval(sim%species(1)%particle%vx), grid, timestep)
            !          phasespace_file_id=0
            !endif
        endif
    end subroutine








subroutine sll_pic_1d_Verlet_scheme(sim, t)
        class( sll_t_simulation_pic1d1v_vp_periodic ), intent( inout ) :: sim
        sll_int32 :: jdx, num_part
        sll_real64, intent(in):: t
        type(sll_t_particle_1d_group), dimension(sll_v_num_species) :: species_1, species_05

        sll_real64 , dimension(:),allocatable :: DPhidx

        !----ALLOCATE STAGES-----------------------------------------------
        call sll_pic_1d_copy_species(sim%species(1:sll_v_num_species) , species_1)
        call sll_pic_1d_copy_species(sim%species(1:sll_v_num_species) , species_05)


        !--------------------Stage 1-------------------------------------------------
        do jdx=1,size(sim%species(1:sll_v_num_species))
            num_part=size(sim%species(jdx)%particle)
            SLL_ALLOCATE(DPhidx(num_part),ierr)

            call sim%fsolver%evalE(sim%species(jdx)%particle%dx, DPhidx)

            species_05(jdx)%particle%vx=sim%species(jdx)%particle%vx+ &
                0.5_f64*sim%tstepw*(DPhidx+sim%Eex(sim%species(jdx)%particle%dx,t))*(-sim%species(jdx)%qm)

            species_05(jdx)%particle%dx=sim%species(jdx)%particle%dx  + sim%tstepw*species_05(jdx)%particle%vx

            SLL_DEALLOCATE_ARRAY(DPhidx,ierr)
        enddo

        call sll_pic1d_adjustweights_advection_species(sim%fsolver,sim%species(1:sll_v_num_species), species_05)

        !call sll_s_pic1d_ensure_boundary_conditions_species(species_05, sim%scenario_int)
        call sim%fsolver%set_species(species_05)
        call sim%fsolver%solve()

        do jdx=1,size(sim%species(1:sll_v_num_species))
            num_part=size(sim%species(jdx)%particle)
            SLL_ALLOCATE(DPhidx(num_part),ierr)

            species_1(jdx)%particle%dx=species_05(jdx)%particle%dx

            call sim%fsolver%evalE(species_1(jdx)%particle%dx, DPhidx)

            species_1(jdx)%particle%vx=species_05(jdx)%particle%vx + &
                0.5_f64*sim%tstepw*(DPhidx+sim%Eex( &
                species_1(jdx)%particle%dx,t+sim%tstepw))*(-species_1(jdx)%qm)

            SLL_DEALLOCATE_ARRAY(DPhidx,ierr)
        enddo

        call sll_s_pic1d_ensure_boundary_conditions_species(species_1, sim%scenario_int)


        call sll_pic1d_adjustweights_advection_species(sim%fsolver,sim%species(1:sll_v_num_species), species_1)
        call sll_pic_1d_copy_species(species_1, sim%species(1:sll_v_num_species))


    end subroutine




    subroutine sll_pic_1d_copy_species(original_species, species )
        type(sll_t_particle_1d_group), dimension(:), intent(in) :: original_species
        type(sll_t_particle_1d_group), dimension(:), intent(inout) :: species
        sll_int32 :: sll_v_num_species, numpart, jdx

        sll_v_num_species=size(original_species)
        SLL_ASSERT(sll_v_num_species==size(species))

        do jdx=1,sll_v_num_species
            numpart= size(original_species(jdx)%particle)
            if (allocated(species(jdx)%particle)) then
                if (size(species(jdx)%particle)==numpart) then
                else
                    SLL_DEALLOCATE_ARRAY(species(jdx)%particle,ierr)
                    SLL_ALLOCATE(species(jdx)%particle(1:numpart),ierr)
                endif
            else
                SLL_ALLOCATE(species(jdx)%particle(1:numpart),ierr)
            endif
            species(jdx)%particle=original_species(jdx)%particle
            species(jdx)%qm=original_species(jdx)%qm
        enddo
    end subroutine

  subroutine sll_pic1d_adjustweights_advection_species(fsolver,species_old, species_new)
        type(sll_t_particle_1d_group), dimension(:), intent(in) :: species_old
        type(sll_t_particle_1d_group), dimension(:), intent(inout) :: species_new
        class( sll_t_pic_1d_field_solver ), pointer,intent(inout) :: fsolver
        sll_real64, dimension(:), allocatable:: ratio
        sll_int32 :: sll_v_num_species, numpart, jdx

        sll_v_num_species=size(species_old)
        SLL_ASSERT(sll_v_num_species==size(species_new))

        !Adjust weights
        if (sll_v_enable_deltaf .eqv. .TRUE.) then
            do jdx=1, sll_v_num_species
                numpart=size(species_old(jdx)%particle)

                SLL_ALLOCATE(ratio(1:numpart),ierr)
                SLL_ASSERT(size(species_new(jdx)%particle)==numpart)
                ratio=1.0_f64
                ratio=sll_f_control_variate_xv(fsolver%BC(species_new(jdx)%particle%dx), &
                                        species_new(jdx)%particle%vx) &
                        /sll_f_control_variate_xv(fsolver%BC(species_old(jdx)%particle%dx), &
                                        species_old(jdx)%particle%vx)
                !                ratio=control_variate_x(species_new(jdx)%particle%dx)&
                    !                                    /control_variate_x(species_old(jdx)%particle%dx)

                species_new(jdx)%particle%weight=species_new(jdx)%particle%weight_const  - &
                    ( species_new(jdx)%particle%weight_const  - species_new(jdx)%particle%weight )*ratio

                !species_new(jdx)%particle%weight_const= species_old(jdx)%particle%weight_const

                ! species_old(jdx)%particle%weight=species_new(jdx)%particle%weight

                SLL_DEALLOCATE_ARRAY(ratio,ierr)
            enddo
        endif

    end subroutine
    
    
    
    
    
    
    
    
    
    
    
    
    
    
     subroutine sll_pic_1d_rungekutta4_step_array(sim, t)
       
        sll_real64 :: h 
        sll_real64, intent(in):: t
        class( sll_t_simulation_pic1d1v_vp_periodic ), intent( inout ) :: sim
        sll_real64, dimension(size(sim%species(1)%particle%dx)) :: x_0
        sll_real64, dimension(size(sim%species(1)%particle%vx)) :: v_0
        sll_real64, dimension(size(x_0)) :: k_x1, k_x2,&
            k_x3, k_x4, k_v1, &
            k_v2, k_v3, k_v4, stage_DPhidx,x_1,v_1

        !x_0=sll_f_pic1d_ensure_periodicity(x_0,  interval_a, interval_b)
        
        h=sim%tstepw 
        x_0=sim%species(1)%particle%dx
        v_0=sim%species(1)%particle%vx
        call sll_s_pic1d_ensure_boundary_conditions(x_0, v_0,  sim%scenario_int)
      
        
               !--------------------Stage 1-------------------------------------------------
        call sim%fsolver%evalE(x_0, stage_DPhidx)
        k_x1= h*v_0
        k_v1= h*(stage_DPhidx+sim%Eex(x_0,t))*(-sim%particle_qm)
        !--------------------Stage 2-------------------------------------------------
        x_1=x_0 + 0.5_f64*k_x1
        v_1=v_0 + 0.5_f64*k_v1
        call sll_s_pic1d_ensure_boundary_conditions(x_1, v_1,  sim%scenario_int)
        call sim%pic1d_adjustweights(x_0,x_1,v_0,v_1)
        call sim%pic_1d_solve_qn(x_1)
        call sim%fsolver%evalE( x_1, stage_DPhidx)
        k_x2= h*(v_1)
        k_v2= h*( stage_DPhidx + sim%Eex(x_1, t+0.5_f64*h))*(-sim%particle_qm)

        !--------------------Stage 3-------------------------------------------------
        x_1=x_0 + 0.5_f64*k_x2
        v_1=v_0 + 0.5_f64*k_v2
        call sll_s_pic1d_ensure_boundary_conditions(x_1, v_1,  sim%scenario_int)
        call sim%pic1d_adjustweights(x_0 + 0.5_f64*k_x1,x_1,v_0 + 0.5_f64*k_v1,v_1)
        call sim%pic_1d_solve_qn(x_1)
        call sim%fsolver%evalE(x_1, stage_DPhidx)
        k_x3= h*(v_1)
        k_v3= h*(stage_DPhidx+ sim%Eex(x_1, t+0.5_f64*h)) *(-sim%particle_qm)
        !--------------------Stage 4-------------------------------------------------
        x_1=x_0 + k_x3
        v_1=v_0 + k_v3
        call sll_s_pic1d_ensure_boundary_conditions(x_1, v_1,  sim%scenario_int)
        call sim%pic1d_adjustweights(x_0 + 0.5_f64*k_x2,x_1,v_0 + 0.5_f64*k_v2,v_1)
        call sim%pic_1d_solve_qn(x_0 +  k_x3)
        call sim%fsolver%evalE(x_1, stage_DPhidx)
        k_x4= h*(v_1)
        k_v4= h*(stage_DPhidx+ sim%Eex(x_1, t+h)) *(-sim%particle_qm)
        !Perform step---------------------------------------------------------------
        x_1= x_0 + (  k_x1 +  2.0_f64 *k_x2 +  2.0_f64*k_x3 + k_x4 )/6.0_f64
        v_1= v_0 + (  k_v1 +  2.0_f64 *k_v2 +  2.0_f64*k_v3 + k_v4)/6.0_f64
        call sll_s_pic1d_ensure_boundary_conditions(x_1, v_1,  sim%scenario_int)

        call sim%pic1d_adjustweights(x_0 + 0.5_f64*k_x3,x_1,v_0 + 0.5_f64*k_v3,v_1)

        !call sll_pic1d_adjustweights(x_0,x_1,v_0,v_1)
        sim%species(1)%particle%dx=x_1
        sim%species(1)%particle%vx=v_1

        call sim%pic_1d_solve_qn(sim%species(1)%particle%dx)

    end subroutine



    !<Merson scheme with built in error estimate for each particle
    subroutine sll_pic_1d_merson4(sim, t)
        sll_real64  :: h
        sll_real64, intent(in):: t
        class( sll_t_simulation_pic1d1v_vp_periodic ), intent( inout ) :: sim
        sll_real64, dimension(size(sim%species(sim%pushed_species)%particle%dx))                :: x_0
        sll_real64, dimension(size(sim%species(sim%pushed_species)%particle%vx))                :: v_0
        sll_real64, dimension(size(sim%species(sim%pushed_species)%particle%weight))                :: weight
        
        sll_real64, dimension(size(x_0)) :: k_x1, k_x2,&
            k_x3, k_x4, k_v1, &
            k_v2, k_v3, k_v4, k_v5, k_x5, stage_DPhidx,x_1,v_1, err_x, err_v
           
        h=sim%tstepw
        x_0=sim%species(sim%pushed_species)%particle%dx
        v_0=sim%species(sim%pushed_species)%particle%vx
        weight=sim%species(sim%pushed_species)%particle%weight
        
        !x_0=sll_f_pic1d_ensure_periodicity(x_0,  interval_a, interval_b)
        call sll_s_pic1d_ensure_boundary_conditions(x_0, v_0,  sim%scenario_int)

        !--------------------Stage 1-------------------------------------------------
        call sim%fsolver%evalE(x_0, stage_DPhidx)
        k_x1= h*v_0
        k_v1= h*(stage_DPhidx+sim%Eex(x_0,t))*(-sim%particle_qm)
        !--------------------Stage 2-------------------------------------------------
        call sim%pic_1d_solve_qn(x_0 + (1.0_f64/3.0_f64)*k_x1)
        call sim%fsolver%evalE(x_0 + (1.0_f64/3.0_f64) *k_x1, stage_DPhidx)

        k_x2= h*(v_0  + (1.0_f64/3.0_f64)*k_v1)
        k_v2= h*( stage_DPhidx + sim%Eex(x_0 + (1.0_f64/3.0_f64)*k_x1, t+(1.0_f64/3.0_f64)*h))*(-sim%particle_qm)
        !--------------------Stage 3-------------------------------------------------
        k_x3= h*(v_0+ (1.0_f64/3.0_f64) * k_v2)
        call sim%pic_1d_solve_qn(x_0 +  (1.0_f64/6.0_f64)*(k_x1 + k_x2 ))

        call sim%fsolver%evalE(x_0 + (1.0_f64/6.0_f64)*(k_x1 + k_x2 ), stage_DPhidx)
        k_v3= h*(stage_DPhidx+ sim%Eex(x_0 + (1.0_f64/6.0_f64)*(k_x1 + k_x2 ), &
            t+(1.0_f64/3.0_f64)*h)) *(-sim%particle_qm)
        !--------------------Stage 4-------------------------------------------------
        k_x4= h*(v_0+ 0.5_f64*k_v3)
        call sim%pic_1d_solve_qn(x_0 +  (k_x1 + 3*k_x3)/8.0_f64  )
        call sim%fsolver%evalE(x_0 +  (k_x1 + 3*k_x3)/8.0_f64  , stage_DPhidx)
        k_v4= h*(stage_DPhidx+ sim%Eex(x_0 +  (k_x1 + 3*k_x3)/8.0_f64  , t+0.5_f64*h)) *(-sim%particle_qm)
        !--------------------Stage 5-------------------------------------------------
        k_x5= h*(v_0+ k_v4)
        call sim%pic_1d_solve_qn(x_0 +  0.5_f64*k_x1  -1.5_f64*k_x3 + 2.0_f64*k_x1)
        call sim%fsolver%evalE(x_0 +  0.5_f64*k_x1  -1.5_f64*k_x3 + 2.0_f64*k_x1, stage_DPhidx)
        k_v4= h*(stage_DPhidx+ sim%Eex(x_0 +  0.5_f64*k_x1  -1.5_f64*k_x3 + 2.0_f64*k_x1, t+h)) *(-sim%particle_qm)

        !Perform step of Order 4-----------------------------------------------------------
        x_1= x_0 + (  k_x1 +  0.0_f64 *k_x2 +  0.0_f64*k_x3 + 4.0_f64*k_x4 + k_x5 )/6.0_f64
        v_1= v_0 + (  k_v1 +  0.0_f64 *k_v2 +  0.0_f64*k_v3 + 4.0_f64*k_v4 + k_v5 )/6.0_f64

        !Perform error step of Order 3 or 5 in linear case
        err_x= x_0 + (  k_x1 +  3*k_x3 + 4*k_x4 + 2*k_x5 )/10.0_f64
        err_v= x_0 + (  k_v1 +  3*k_v3 + 4*k_v4 + 2*k_v5 )/10.0_f64
        err_x=x_1 - err_x
        err_v=v_1 - err_v

        call sim%pic_1d_calc_push_error(err_x,weight,t)

        call sll_s_pic1d_ensure_boundary_conditions(x_1, v_1,  sim%scenario_int)

        call sim%pic1d_adjustweights(x_0,x_1,v_0,v_1)
        x_0=x_1
        v_0=v_1

        call sim%pic_1d_solve_qn(x_0)
        sim%species(sim%pushed_species)%particle%dx=x_0
        sim%species(sim%pushed_species)%particle%vx=v_0

    end subroutine

    subroutine sll_pic_1d_rungekutta2(sim,  t)
        class( sll_t_simulation_pic1d1v_vp_periodic ), intent( inout ) :: sim
        sll_real64   :: h
        sll_real64, intent(in):: t
        
        sll_real64, dimension(size(sim%species(sim%pushed_species)%particle%dx))                :: x_0
        sll_real64, dimension(size(sim%species(sim%pushed_species)%particle%vx))                :: v_0
        sll_real64, dimension(size(x_0)) :: k_x1, k_x2,&
            k_v1, k_v2, stage_DPhidx,x_1,v_1
        
         h=sim%tstepw
         x_0=sim%species(sim%pushed_species)%particle%dx
         v_0=sim%species(sim%pushed_species)%particle%vx
        !x_0=sll_f_pic1d_ensure_periodicity(x_0,  interval_a, interval_b)
        call sll_s_pic1d_ensure_boundary_conditions(x_0, v_0,  sim%scenario_int)

        !--------------------Stage 1-------------------------------------------------
        !call sll_bspline_fem_solver_1d_solve(x_0)
        call sim%fsolver%evalE(x_0, stage_DPhidx)
        k_x1= h*v_0
        k_v1= h*(stage_DPhidx+sim%Eex(x_0,t)) *(-sim%particle_qm)
        !--------------------Stage 2-------------------------------------------------
        call sim%pic_1d_solve_qn(x_0 + 0.5_f64 *k_x1)
        call sim%fsolver%evalE(x_0 + 0.5_f64 *k_x1, stage_DPhidx)
        k_x2= h*(v_0  + 0.5_f64  *k_v1)
        k_v2= h*( stage_DPhidx + sim%Eex(x_0 + 0.5_f64 *k_x1,t + 0.5_f64*h))&
            *(-sim%particle_qm)
        !Perform step---------------------------------------------------------------
        x_1= x_0 + k_x2
        v_1= v_0 + k_v2

        !        x_0=sll_f_pic1d_ensure_periodicity(x_0,  interval_a, interval_b)
        call sll_s_pic1d_ensure_boundary_conditions(x_1, v_1,  sim%scenario_int)

        call sim%pic1d_adjustweights(x_0,x_1,v_0,v_1)
        x_0=x_1
        v_0=v_1

        call sim%pic_1d_solve_qn(x_0)
        sim%species(sim%pushed_species)%particle%dx=x_0
        sim%species(sim%pushed_species)%particle%vx=v_0
    end subroutine


    subroutine sll_pic_1d_rungekutta3(sim, t)
        class( sll_t_simulation_pic1d1v_vp_periodic ), intent( inout ) :: sim
        sll_real64            :: h
        sll_real64, intent(in):: t
        sll_real64, dimension(size(sim%species(sim%pushed_species)%particle%dx))                :: x_0
        sll_real64, dimension(size(sim%species(sim%pushed_species)%particle%vx))                :: v_0
        sll_real64, dimension(size(x_0)) :: k_x1, k_x2,&
            k_x3, k_v3, k_v1,k_v2, stage_DPhidx

        
        
        h=sim%tstepw
        x_0=sim%species(sim%pushed_species)%particle%dx
        v_0=sim%species(sim%pushed_species)%particle%vx
        call sll_s_pic1d_ensure_boundary_conditions(x_0, v_0,  sim%scenario_int)

        !--------------------Stage 1-------------------------------------------------
        call sim%fsolver%evalE(x_0, stage_DPhidx)
        k_x1= h*v_0
        k_v1= h*(stage_DPhidx+sim%Eex(x_0,t)) *(-sim%particle_qm)
        !--------------------Stage 2-------------------------------------------------
        k_x2= h*(v_0  + 0.5_f64  *k_v1)
        call sim%pic_1d_solve_qn(x_0 + 0.5_f64 *k_x1)
        call sim%fsolver%evalE(x_0 + 0.5_f64 *k_x1, stage_DPhidx)
        k_v2= h*( stage_DPhidx + sim%Eex(x_0 + 0.5_f64 *k_x1,t + 0.5_f64*h))&
            *(-sim%particle_qm)
        !--------------------Stage 3-------------------------------------------------
        k_x3= h*(v_0  - k_v1 + 2.0_f64*k_v2)
        call sim%pic_1d_solve_qn(x_0 - k_x1 + 2.0_f64*k_x2)
        call sim%fsolver%evalE(x_0 - k_x1 + 2.0_f64*k_x2, stage_DPhidx)
        k_v3= h*( stage_DPhidx + sim%Eex(x_0 - k_x1 + 2.0_f64*k_x2,t + h))&
            *(-sim%particle_qm)

        !Perform step---------------------------------------------------------------
        x_0= x_0 + (k_x1 + 4.0_f64 * k_x2 + k_x3)/6.0_f64
        v_0= v_0 + (k_v1 + 4.0_f64 * k_v2 + k_v3)/6.0_f64

        call sll_s_pic1d_ensure_boundary_conditions(x_0, v_0,  sim%scenario_int)

        call sim%pic_1d_solve_qn(x_0)
        sim%species(sim%pushed_species)%particle%vx=v_0
        sim%species(sim%pushed_species)%particle%dx=x_0
    end subroutine


    subroutine sll_pic_1d_heun(sim, t)
        class( sll_t_simulation_pic1d1v_vp_periodic ), intent( inout ) :: sim
        sll_real64            :: h
        sll_real64, intent(in):: t
        sll_real64, dimension(size(sim%species(sim%pushed_species)%particle%dx))                :: x_0
        sll_real64, dimension(size(sim%species(sim%pushed_species)%particle%vx))                :: v_0
        sll_real64, dimension(size(x_0)) :: k_x1, k_x2 , &
            k_v1 , k_v2 ,stage_DPhidx,x_1,v_1

        !x_0=sll_f_pic1d_ensure_periodicity(x_0,  interval_a, interval_b)
        
         h=sim%tstepw
         x_0=sim%species(sim%pushed_species)%particle%dx
         v_0=sim%species(sim%pushed_species)%particle%vx
         
        call sll_s_pic1d_ensure_boundary_conditions(x_0, v_0,  sim%scenario_int)
        !--------------------Stage 1-------------------------------------------------
        !call sll_bspline_fem_solver_1d_solve(x_0)
        call sim%fsolver%evalE(x_0, stage_DPhidx)
        k_x1= h*v_0
        k_v1= h*(stage_DPhidx+sim%Eex(x_0,t))*(-sim%particle_qm)
        !call sll_pic1d_adjustweights(x_0,x_0+ k_x1,v_0,v_0 + k_v1)
        !--------------------Stage 2-------------------------------------------------
        call sim%pic_1d_solve_qn(x_0 + k_x1)
        call sim%fsolver%evalE(x_0 + k_x1, stage_DPhidx)
        k_x2= h*(v_0  + k_v1)
        k_v2= h*( stage_DPhidx+sim%Eex(x_0 +k_x1,t+h))*(-sim%particle_qm)

        !Perform step---------------------------------------------------------------
        x_1= x_0 + 0.5_f64 *(k_x1+k_x2   )
        v_1= v_0 + 0.5_f64 *(k_v1+k_v2  )
        !x_0=sll_f_pic1d_ensure_periodicity(x_0,  interval_a, interval_b)
        call sll_s_pic1d_ensure_boundary_conditions(x_1, v_1,  sim%scenario_int)

        !call sll_pic1d_adjustweights(x_0+ k_x1,x_1,v_0 + k_v1,v_1)
        call sim%pic1d_adjustweights(x_0,x_1,v_0,v_1)
        x_0=x_1
        x_0=v_1

        call sim%pic_1d_solve_qn(x_0)
        sim%species(sim%pushed_species)%particle%dx=x_0
        sim%species(sim%pushed_species)%particle%vx=v_0
    end subroutine

    subroutine sll_pic_1d_explicit_euler(sim,  t)
    
 
        !type(sll_t_particle_1d_group), dimension(:), intent(inout) :: species_0
        sll_int32 :: jdx, num_part
        sll_real64            :: h
        sll_real64, intent(in):: t
        class( sll_t_simulation_pic1d1v_vp_periodic ), intent( inout ) :: sim
        
        type(sll_t_particle_1d_group), dimension(sll_v_num_species) :: species_1

        sll_real64 , dimension(:),allocatable :: stage_DPhidx
        
         
         h=sim%tstepw

        !----ALLOCATE STAGES-----------------------------------------------
        call sll_pic_1d_copy_species(sim%species(1:sll_v_num_species), species_1)

        !--------------------Stage 1-------------------------------------------------
        do jdx=1,sll_v_num_species
            num_part=size(sim%species(jdx)%particle)
            SLL_ALLOCATE(stage_DPhidx(num_part),ierr)

            call sim%fsolver%evalE(sim%species(jdx)%particle%dx, stage_DPhidx)

            species_1(jdx)%particle%dx=sim%species(jdx)%particle%dx + h*sim%species(jdx)%particle%vx
            species_1(jdx)%particle%vx=sim%species(jdx)%particle%vx + h*(stage_DPhidx + &
                sim%Eex( sim%species(jdx)%particle%dx,t) )*(-sim%species(jdx)%qm)

            SLL_DEALLOCATE_ARRAY(stage_DPhidx,ierr)

        enddo

        call sll_s_pic1d_ensure_boundary_conditions_species(species_1, sim%scenario_int)

        call sll_pic1d_adjustweights_advection_species(sim%fsolver,sim%species(1:sll_v_num_species), species_1)

        call sll_pic_1d_copy_species(species_1, sim%species(1:sll_v_num_species))

        call sim%fsolver%set_species(sim%species(1:sll_v_num_species))
        call sim%fsolver%solve()

    end subroutine
    
    
        subroutine sll_pic_1d_variational_leap_frog(sim)
        
        class( sll_t_simulation_pic1d1v_vp_periodic ), intent( inout ) :: sim
        sll_real64                              :: h
        sll_real64, dimension(size(sim%species(sim%pushed_species)%particle%dx))  :: x_0
        sll_real64, dimension(size(sim%species(sim%pushed_species)%particle%dx))  :: v_0
        sll_real64, dimension(size(x_0)) ::  DPhidx_0, DPhidx_1, x_1, v_1
!       sll_real64, dimension(size(x_0)) ::  DPhidx_1
       
        h=sim%tstepw
        x_0=sim%species(sim%pushed_species)%particle%dx
        v_0=sim%species(sim%pushed_species)%particle%vx
        call sim%fsolver%evalE(x_0, DPhidx_0)
        x_1=x_0+ h*v_0 - ((h**2)/2.0_f64) *DPhidx_0*(-sim%particle_qm)
        x_1=sll_f_pic1d_ensure_periodicity(x_1,  sim%interval_a, sim%interval_b)
        call sll_s_pic1d_ensure_boundary_conditions(x_1, v_0,  sim%scenario_int)

        call sim%pic_1d_solve_qn(x_1)
        call sim%fsolver%evalE(x_1, DPhidx_1)

        v_1=v_0 + 0.5_f64*h* (DPhidx_0 + DPhidx_1)*(-sim%particle_qm)

        !Push
        sim%species(sim%pushed_species)%particle%dx=x_1
        sim%species(sim%pushed_species)%particle%vx=v_1
    end subroutine

    subroutine sll_pic_1d_leap_frog(sim)
    
         class( sll_t_simulation_pic1d1v_vp_periodic ), intent( inout ) :: sim
        sll_real64                              :: h
        sll_real64, dimension(size(sim%species(sim%pushed_species)%particle%dx))                :: x_0
        sll_real64, dimension(size(sim%species(sim%pushed_species)%particle%vx))                :: v_0
        sll_real64, dimension(size(x_0)) ::  DPhidx_0, x_1, v_1
       
        h=sim%tstepw
        x_0 =sim%species(sim%pushed_species)%particle%dx
        v_0 =sim%species(sim%pushed_species)%particle%vx
        x_0=sll_f_pic1d_ensure_periodicity(x_0,  sim%interval_a, sim%interval_b)
        call sim%fsolver%evalE(x_0, DPhidx_0)
        !DPhidx_0=0

        v_1= v_0 - 0.5_f64 *DPhidx_0*h*(-sim%particle_qm)
        v_1=v_1+  DPhidx_0*h
        !v_1=v_0+ 0.5_f64* DPhidx_0*h
        x_1=x_0+ v_1*h

        !Check the CFL timestep criterion
        !if  ((maxval(v_1) *h)/2>1) print *, "CFL!!"

        !        x_1=x_0 + (2.0/(128.0))
        !        v_1=v_0
        !Push
        x_0=x_1
        v_0=v_1

        !x_0=sll_f_pic1d_ensure_periodicity(x_0,  interval_a, interval_b)
        call sll_s_pic1d_ensure_boundary_conditions(x_0, v_0,  sim%scenario_int)

        call sim%pic_1d_solve_qn(x_0)
        sim%species(sim%pushed_species)%particle%dx=x_0
        sim%species(sim%pushed_species)%particle%vx=v_0
    end subroutine
    
    
    
    
    subroutine sll_pic_1d_solve_qn(sim, particleposition_selected)
        
        class( sll_t_simulation_pic1d1v_vp_periodic ), intent( inout ) :: sim
        sll_real64, dimension(:),intent(in) :: particleposition_selected
        sll_int32 ::jdx
        !Add all the other species as constant
        call sim%fsolver%reset_particles()

        do jdx=1,sll_v_num_species
            if (jdx/=sim%pushed_species) then
                call sim%fsolver%add_species(sim%species(jdx))
            endif
        enddo

        selectcase(sim%scenario_int)
            !            case(sll_p_pic1d_testcase_quiet)
            case(sll_p_pic1d_testcase_landau)
                !All particles, are only electrons
                !load constant ion background
                SLL_ASSERT(size(particleposition_selected)==size(sim%species(sim%pushed_species)%particle))

                call sim%fsolver%set_electrons_only_weighted(particleposition_selected, sim%species(sim%pushed_species)%particle%weight)
                call sim%fsolver%solve()
            case(sll_p_pic1d_testcase_ionbeam)
                SLL_ASSERT(size(particleposition_selected)==size(sim%species(sim%pushed_species)%particle))

                call sim%fsolver%add_particles_weighted(particleposition_selected, &
                    sign(sim%species(sim%pushed_species)%particle%weight,sim%species(sim%pushed_species)%qm) )
                call sim%fsolver%solve()
            case default
                SLL_ASSERT(sll_v_num_species==1)
                !All particles, are only electrons
                SLL_ASSERT(size(particleposition_selected)==size(sim%species(sim%pushed_species)%particle))
                call sim%fsolver%set_electrons_only_weighted(particleposition_selected, sim%species(sim%pushed_species)%particle%weight)
                call sim%fsolver%solve()

                !load constant ion background
                !                SLL_ASSERT(size(allparticleposition)==size(particleweight))
                !                call fsolver%set_electrons_only_weighted(allparticleposition, particleweight)
                !                call fsolver%solve()
        endselect

    end subroutine    
    
    
    
    
    
    
    
    
    subroutine sll_pic1d_adjustweights(sim,xold, xnew, vold, vnew)
        class( sll_t_simulation_pic1d1v_vp_periodic ), intent( inout ) :: sim
        sll_real64, dimension(:),intent(in) ::vold, xold
        sll_real64, dimension(:),intent(in) ::vnew, xnew
        sll_real64, dimension(size(sim%species(sim%pushed_species)%particle)):: ratio
        sll_int32 :: N

        N=size(xold)
        SLL_ASSERT(N==size(xnew))
        SLL_ASSERT(N==size(vold))
        SLL_ASSERT(N==size(vnew))

        SLL_ASSERT(N==size( sim%species(sim%pushed_species)%particle))
        !Adjust weights
        if (sim%deltaf .eqv. .TRUE.) then
            ratio=sll_f_control_variate_xv(sim%fsolver%BC(xnew), vnew )/&
                                sll_f_control_variate_xv(sim%fsolver%BC(xold), vold )
            !ratio=ratio/(N*coll_size)
            !ratio=ratio
            sim%species(sim%pushed_species)%particle%weight=sim%species(sim%pushed_species)%particle%weight_const &
                -(sim%species(sim%pushed_species)%particle%weight_const-sim%species(sim%pushed_species)%particle%weight)*ratio
            
        endif

    end subroutine
    
    
    
    
    
    
    
    ! write result in file subroutine
    
    
        subroutine sll_pic1d_write_result(sim,filename)
        
        class( sll_t_simulation_pic1d1v_vp_periodic ), intent( inout ) :: sim
        character(len=*), intent(in) :: filename
!        sll_real64, dimension(:), intent(in) :: kinetic_energy, electrostatic_energy, impulse,&
!            particleweight_mean,particleweight_var, perror_mean  , perror_var
        integer :: k,file_id,file_id_err
        sll_real64,  dimension(size(sim%kineticenergy)) :: total_energy


        if (coll_rank==0) then

            SLL_ASSERT(size(sim%kineticenergy)==sim%tsteps+1)
            SLL_ASSERT(size(sim%fieldenergy)==sim%tsteps+1)
            SLL_ASSERT(size(sim%impulse)==sim%tsteps+1)

            total_energy=sim%kineticenergy +sim%fieldenergy

            

            call sll_s_new_file_id(file_id, ierr)

            !Write Data File
            !            open(file_id, file = plot_name//"_"//fin//'.dat' )
            open(file_id, file = trim(sim%root_path)//filename//'.dat')

            write (file_id, *)  "#Full 1d1v Electrostatic PIC"
            write (file_id,*)  "#Time steps:", sim%tsteps
            write (file_id,*)  "#Time stepwidth:", sim%tstepw
            write (file_id,*)  "#Marko particles:", coll_size*sim%nmark
            write (file_id,*)  "#Particle Pusher:", sim%ppusher
            write (file_id,*)  "#Finite Elements: 2^(", log(real(sim%mesh_cells,i64))/log(2.0_f64),")"
            write (file_id,*)  "#Size of MPI Collective: ", coll_size

            write (file_id,*)  "time  ", "kineticenergy  ", "electrostaticenergy  ","total energy   ", "impulse  ", &
                "vthermal  ", "weightmean   ", "weightvar  ", "perrormean  ", "perrorvar   ","total momentum           "
            do k=1,sim%tsteps+1
                write (file_id,*)  sim%tstepw*(k-1), sim%kineticenergy(k), sim%fieldenergy(k) ,total_energy(k), sim%impulse(k), &
                    sim%thermal_velocity_estimate(k), sim%particleweight_mean(k),sim%particleweight_var(k),&
                    sim%push_error_mean(k), sim%push_error_var(k), sim%total_momentum(k)
            enddo
            close(file_id)

            call sll_s_new_file_id(file_id_err, ierr)

            open(file_id_err, file = trim(sim%root_path)//filename//'-errors.dat')

            if (sim%impulse(1)/=0 .AND. total_energy(1)/=0) then
                !Write the relative errors
                do k=1,sim%tsteps+1
                    write (file_id_err,*)  sim%tstepw*(k-1),&
                        abs((total_energy(k)-total_energy(1))/total_energy(1)), &
                        abs((  sim%impulse(k)-sim%impulse(1)))
                enddo
            endif
            close(file_id_err)


            !Write Gnuplot file
            open(file_id, file = trim(sim%root_path)//filename//'.gnu')

            !First Plot
            write (file_id, *) "set term x11 1"

            write(file_id,"(A17,G10.3,A4,G10.3,A1)")"set title 'alpha=",sim%lalpha,", k=",sim%lmode,"'"
            write(file_id, "(A14,G10.5,A1)") "set xrange [0:",sim%tstepw*sim%tsteps,"]"
            write (file_id,*)  "set logscale y"
            write (file_id,*)  "unset logscale x"
            write (file_id,*)  "l_alpha=",sim%lalpha
            write (file_id,"(A21,G4.3,A18,G10.3,A1)")  "set xlabel 'Time t in", sim%tsteps, &
                " steps, stepwidth=", sim%tstepw, "'"
            write (file_id,*)  "set ylabel 'Energy'"


            !Plot expected damping curve if known
            if (sim%lalpha==0.001_f64 .AND. sim%lmode==0.4_f64) then

                write (file_id,*) "E(x)=abs(0.002*0.424666*exp(-0.0661*x)*cos(1.2850*x-0.3357725))**2"
                write (file_id,*)  "plot '"//filename//".dat' using 1:3 with lines,\"
                write (file_id,*) "E(x) with lines linestyle 2;"
            elseif (sim%lalpha/=0) then
                write (file_id,*) "E(x)=0.1*exp(-0.1533*2*x)"
                write (file_id,*)  "plot '"//filename//".dat' using 1:3 with lines,\"
                write (file_id,*) "E(x) with lines linestyle 2;"

            else

                write (file_id,*)  "plot '"//filename//".dat' using 1:3 with lines"


            endif
            write (file_id,*) "set autoscale x; unset logscale;unset xlabel;unset ylabel;unset xrange;unset yrange;"


            !Plot Kernel density estimates for initial distribution
            !!!Enter v_thermal and two stream here
            write (file_id, *) "set term x11 2"
            write (file_id, *)  "set multiplot layout 2,1 rowsfirst title 'Distribution Kernel Density Estimates'"
            write (file_id, *) "set autoscale x; set autoscale y"
            write (file_id,*)  "vel_pdf(x)=1/sqrt(2*pi)*exp(-0.5*(x**2))"
            write (file_id,*)  "set title 'Velocity Distribution'"
            write (file_id,*)  "plot '"//trim(sim%root_path)//"initial_phasespace.dat' using 2:3 smooth kdensity, vel_pdf(x)"
            write (file_id,*)  "set title 'Spacial Distribution'"
            write (file_id,*)  "plot '"//trim(sim%root_path)//"initial_phasespace.dat' using 1:3 smooth kdensity"
            write (file_id,*) "unset multiplot"


            write (file_id, *) "set term x11 3"
            write (file_id, *) "set autoscale x; set autoscale y"
            write (file_id,*)  "set title 'Thermal velocity estimate'"
            write (file_id,*)  "set xlabel 'time'"
            write (file_id,*)  "set ylabel 'v_th'"
            write (file_id,*)  "plot '"//filename//".dat' using 1:5 with lines"

            write (file_id, *) "set term x11 4"
            write (file_id, *) "set autoscale x; set autoscale y; unset logscale y"
            write (file_id,*)  "set title 'Absolute Energies'"
            write (file_id,*)  "set xlabel 'time'"
            write (file_id,*)  "set ylabel 'energy'"
            write (file_id,*)  "plot '"//filename//".dat' using 1:($2+$3) with lines\"
            write (file_id,*)  "title 'total (kinetic+electrostatic)',\"
            write (file_id,*)  "'"//filename//".dat' using 1:2 with lines \"
            write (file_id,*)  "title 'kinetic'"

            write (file_id, *) "set term x11 5"
            write (file_id, *) "set autoscale x; set logscale y"
            write (file_id,*)  "set title 'Relative Energy Error'"
            write (file_id,*)  "set xlabel 'time'"
            write (file_id,*)  "set ylabel 'rel. error'"
            write (file_id,*)  "plot '"//filename//"-errors.dat' using 1:2 with lines"

            write (file_id, *) "set term x11 6"
            write (file_id, *) "set autoscale x; set logscale y"
            write (file_id,*)  "set title 'Total Impulse Error'"
            write (file_id,*)  "set xlabel 'time'"
            write (file_id,*)  "set ylabel 'total. error'"
            write (file_id,*)  "plot '"//filename//"-errors.dat' using 1:3 with lines"

        if (sim%deltaf .eqv. .TRUE.) then
            write (file_id, *) "set term x11 7"
            write (file_id, *)  "set multiplot layout 2,1 rowsfirst title 'Time Development of Particle Weights'"
            write (file_id, *) "set autoscale x; set autoscale y"
            write (file_id,*)  "set title 'Mean'"
            write (file_id,*)  "plot '"//filename//".dat' using 1:6 with lines"
            write (file_id,*)  "set title 'Variance'"
            write (file_id,*)  "plot '"//filename//".dat' using 1:7 with lines"
            write (file_id,*) "unset multiplot"
        endif

            write (file_id, *) "set term x11 8"
            write (file_id, *)  "set multiplot layout 2,1 rowsfirst title 'Time Development of local Pusher Error in x'"
            write (file_id, *) "set autoscale x; set autoscale y"
            write (file_id,*)  "set title 'Mean'"
            write (file_id,*)  "plot '"//filename//".dat' using 1:8 with lines"
            write (file_id,*)  "set title 'Variance'"
            write (file_id,*)  "plot '"//filename//".dat' using 1:9 with lines"
            write (file_id,*) "unset multiplot"
        endif
        
    end subroutine

    
    subroutine sll_pic_1d_calc_push_error(sim,perror, pweight,t)
        class( sll_t_simulation_pic1d1v_vp_periodic ), intent( inout ) :: sim
        sll_real64, dimension(:), intent(in) :: perror
        sll_real64, dimension(:), intent(in) :: pweight
        sll_real64, intent(in)              :: t
        integer                              :: timestep
        timestep=1+NINT(t/sim%tstepw)
        sim%push_error_mean(timestep)=dot_product(perror, pweight)
        sim%push_error_var(timestep)=dot_product(perror**2, pweight)&
            - sim%push_error_mean(timestep)

        sim%push_error_mean=maxval(perror)
    end subroutine

    

end module sll_m_sim_pic_vp_1d1v_cart












