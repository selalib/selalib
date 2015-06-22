module sll_module_simulation_pic1d1v_vp_periodic

#include "sll_working_precision.h"
#include "sll_errors.h"
#include "sll_memory.h"
#include "sll_assert.h"

  use sll_simulation_base, only: sll_simulation_base_class
  
  use sll_collective , only :       sll_world_collective , sll_collective_barrier ,&
     sll_boot_collective
  use sll_timer , only :   sll_set_time_mark, sll_time_mark,  sll_time_elapsed_between
  
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
  sll_int32 :: ierr,i,j
   
  type, extends( sll_simulation_base_class ) :: &
      sll_simulation_pic1d1v_vp_periodic

    CHARACTER(LEN=256) :: path

    ! TODO: use new enumerations
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

    logical :: pi_unit
    logical :: deltaf
    logical :: gnuplot_inline_output
    sll_real64 ::  particle_qm !<mass to electron mass ratio, with the sign of the charge
    class( pic_1d_field_solver ), pointer :: fsolver
    sll_real64, allocatable               :: knots(:) 
    sll_real64, allocatable               :: electricpotential_interp(:) 
    sll_int32                             :: pushed_species
    sll_real64, allocatable               :: particle(:,:)
    type(sll_particle_1d_group),dimension(10) :: species
    sll_real64 :: kineticenergy_offset, impulse_offset
    sll_real64 :: initial_fem_inhom , fem_inhom
    sll_int64    ::  fastest_particle_idx
    sll_real64, allocatable:: eval_solution(:)
    sll_real64,  allocatable:: particleposition(:)
    sll_real64,  allocatable :: particlespeed(:)
    sll_real64,  allocatable:: steadyparticleposition(:) !Steady non-moving objects aka Ions
    sll_real64, dimension(:), allocatable :: fieldenergy, kineticenergy, impulse, &
    thermal_velocity_estimate, particleweight_mean ,particleweight_var, inhom_var ,&
    push_error_mean, push_error_var
    
    
    procedure( sll_pic_1d_electric_field_external ), pointer :: Eex   
  contains
  
    procedure :: run            => run_pic_1d
    procedure :: init_from_file => init_from_file
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
            print *, "Size of MPI-Collective: ", coll_size               !number of cores involved in calculus (mpirun argument)
        endif
!parallel configuration
        call sll_collective_barrier(sll_world_collective)
        call sll_initialize_intrinsic_mpi_random(sll_world_collective)
        call sll_collective_barrier(sll_world_collective)

        !Set up Quasineutral solver
        
        !fsolver is of class pic_1d_field_solver
        

         
        select case(sim%scenario_int)
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





  subroutine init_from_file( sim, filename )
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
      sim%gnuplot_inline_output=gnuplot_inline_output_user
      
      sim%deltaf=deltaf
      sim%nstreams=nstreams
    
      sim%sdeg=sdeg
      sim%tstepw=tstepw
      sim%femp=femp
  
      sim%ppusher_int= match_enumeration("ppusher",sim%ppusher) 
      sim%psolver_int= match_enumeration("psolver",sim%psolver)
      sim%scenario_int= match_enumeration("scenario",sim%scenario)
      
  end subroutine init_from_file




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

  !============================================================================









  subroutine run_pic_1d(sim)
     class( sll_simulation_pic1d1v_vp_periodic ), intent( inout ) :: sim
     sll_int32 :: idx,jdx,kdx
     sll_real64 :: time,timestep
     type(sll_time_mark)  :: tstart, tstop


        sim%nmark=nmark/coll_size
        print*, "#Core ", coll_rank, " handles particles", coll_rank*nparticles +1, "-", (coll_rank+1)*nparticles
        call sll_collective_barrier(sll_world_collective)
        if (coll_rank==0) then
            print*, "#Total Number of particles: ", sim%nmark*coll_size
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
            
            !particle caracterized  outputs
            SLL_CLEAR_ALLOCATE( sim%steadyparticleposition(1:sim%nmark),ierr)
            SLL_CLEAR_ALLOCATE(sim%particleposition(1:sim%nmark),ierr)
            SLL_CLEAR_ALLOCATE(sim%particlespeed(1:sim%nmark),ierr)
            SLL_CLEAR_ALLOCATE(sim%eval_solution(1:size(sim%knots)),ierr)   
    
            if (coll_rank==0) then
                write(*,*) "#PIC1D: Loading particles..."
                call sim%fsolver%set_num_sample(nparticles*coll_size)
                call load_particle_species (sim%nmark, sim%interval_a, sim%interval_b ,&
                 sim%species)
                sim%particle_qm=sim%species(1)%qm
                call sll_collective_barrier(sll_world_collective)
                if (coll_rank==0) then 
                   write(*,*) "#PIC1D: Particles loaded..."
               
                  !Determine index of fastest particle
                   fastest_particle_idx=MAXLOC(particlespeed, dim=1)   


               !Do precalculations, especially for the FEM solver  
               !call sll_bspline_fem_solver_1d_initialize(knots, spline_degree, steadyparticleposition, -(sll_e_charge/sll_epsilon_0) )
                     if (coll_rank==0) then
                            print *, "#Number of FEM-Cells ",sim%mesh_cells



                            select case(sim%scenario_int)
                              case(SLL_PIC1D_TESTCASE_QUIET)
                                call sim%fsolver%set_ions_constant(0.0_f64)
                                sim%particle_qm=-1.0_f64
                              case(SLL_PIC1D_TESTCASE_LANDAU)
                                call sim%fsolver%set_ions_constant(0.0_f64)
                              case(SLL_PIC1D_TESTCASE_IONBEAM)
                                call sim%fsolver%set_ions_constant(0.0_f64)
                              case(SLL_PIC1D_TESTCASE_IONBEAM_ELECTRONS)
                                call sim%fsolver%set_ions_constant_particles(steadyparticleposition)
                              case default
                                sim%particle_qm=-1.0_f64 !By default we simulate Electrons
                                call sim%fsolver%set_ions_constant(0.0_f64)
                            end select



                           if (sim%deltaf .eqv. .TRUE.) then
                                  !The Control Variate, as it is constant over time gos into the constant inhomogenity
                                  !In the constant electron case this should be normalized to zero
                                  !But since the fem solver ignores the constant ion background anyway
                                  !this is going to be ok.
                                  !Also this should be done analytically and not like that as an Monte Carlo estimate
                                  !Here suppose the control variate is the initial state by default
                                  call  sim%fsolver%set_bg_particles(sim%species(1)%particle%dx, &
                                  sign(1.0_f64, sim%species(1)%qm)*sim%species(1)%particle%weight_const)

                                sim%kineticenergy_offset=sll_pic1d_calc_kineticenergy_offset(sim%species(1:num_species))
                                sim%impulse_offset=sll_pic1d_calc_impulse_offset(sim%species(1:num_species))
                          else
                                !Here should be a quick deactivate since every code is primarily deltaf
                                sim%kineticenergy_offset=0.0_f64
                                sim%impulse_offset=0.0_f64
                         endif


!Start with PIC Method
        !Do the initial solve of the field
        !particleposition=sll_pic1d_ensure_periodicity(particleposition,  interval_a, interval_b)
        !call sll_pic1d_ensure_boundary_conditions(species(1)%particle%dx, species(1)%particle%vx)


           call sll_pic1d_ensure_boundary_conditions_species(sim%species(1:num_species)  )


        !Initial Solve
           call sim%fsolver%reset_particles()
           call sim%fsolver%set_species(sim%species(1:num_species))
           call sim%fsolver%solve()
           print *, "#Initial field solve"
           print *, "#Test initial field, and loading condition for landau damping"
           call sim%fsolver%evalE(sim%knots(1:sim%mesh_cells)+ (sll_pi/6)*(sim%knots(2)-sim%knots(1)), sim%eval_solution(1:sim%mesh_cells))
           
           call  sll_pic_1d_initial_error_estimates()
  !Information on initial field

          print *, "#Initial Field average: " , sum( sim%fsolver%getPotential())
          print *, "#Initial Field variance: " , sum((sum(  sim%fsolver%getPotential())/sim%mesh_cells &
            -  sim%fsolver%getPotential())**2)/sim%mesh_cells


        !---------------------------Initial Plots------------------------------------------------------------


          open(unit=20, file=trim(root_path)//"initial_field.txt")
          write (20,*) sim%eval_solution(1:sim%mesh_cells)
          close(20)
          print *, "#Field Energy                  " , &
                "Kinetic Energy                  ", "Inhom. Variance       "
        endif
!Write initial phasespace
        do jdx=0, coll_size
            call sll_collective_barrier(sll_world_collective)
            if( coll_rank==jdx) then
                if (jdx==0) then
                    open(unit=20, file=trim(root_path)//"initial_phasespace.dat")
                else
                    open(unit=20, file=trim(root_path)//"initial_phasespace.dat",position="append")
                endif
                do kdx=1, num_species
                    do idx=1, size(sim%species(kdx)%particle)
                        write (20,*) sim%species(kdx)%particle(idx)%dx, sim%species(kdx)%particle(idx)%vx, sim%species(kdx)%particle(idx)%weight
                    enddo
                enddo
                close(20)
            endif
        enddo
        call sll_set_time_mark(tstart)
!---------------------------### Main Loop
        time=0.0_f64
        do timestep=1, tsteps
            time=tstepw*tstepw*(timestep-1)    
            call sll_pic1d_write_phasespace(timestep,sim%species(1:num_species))
            if (sim%nmark*coll_size<=10000 .AND. (sim%gnuplot_inline_output .eqv. .true.) .AND. (coll_rank==0) .AND. mod(timestep-1,sim%tsteps/100)==0) then
                call particles_center_gnuplot_inline(pic_1d_allparticlepos(sim%species(1:num_species)), &
                    pic_1d_allparticlev(sim%species(1:num_species)), &
                    sim%interval_a, sim%interval_b,&
                    -1.1_f64*maxval(pic_1d_allparticlev(sim%species(1:num_species))),&
                    1.1_f64*maxval( pic_1d_allparticlev(sim%species(1:num_species))),&
                    time)

            end if
        
            sim%particleweight_mean(timestep)=sum(sim%species(1)%particle%weight)/size(sim%species(1)%particle%weight)
            
            sim%particleweight_var(timestep)=sum(sim%species(1)%particle%weight**2)/size(sim%species(1)%particle)-sim%particleweight_mean(timestep)**2


             sim%kineticenergy(timestep)=sll_pic1d_calc_kineticenergy( sim%species(1:num_species) )
             sim%fieldenergy(timestep)=sll_pic1d_calc_fieldenergy(sim%species(1:num_species))
!            fieldenergy(1)=landau_alpha**2/(2*landau_mode**2)/2.0_f64
!            kineticenergy(1)=0.5_f64

             sim%impulse(timestep)=sll_pic1d_calc_impulse(sim%species(1:num_species))

             sim%thermal_velocity_estimate(timestep)=sll_pic1d_calc_thermal_velocity(sim%species(1)%particle%vx,sim%species(1)%particle%weight)
             sim%inhom_var(timestep)=sim%fsolver%calc_variance_rhs()
               
        
            if ( (sim%gnuplot_inline_output.eqv. .true.) .AND. coll_rank==0 .AND. mod(timestep-1,tsteps/100)==0  ) then
                call energies_electrostatic_gnuplot_inline(sim%kineticenergy(1:timestep), sim%fieldenergy(1:timestep), sim%impulse(1:timestep),sim%tstepw)
            else

            endif
            if ( (sim%gnuplot_inline_output.eqv. .true.) .AND. coll_rank==0 .AND. mod(timestep-1,tsteps/100)==0 ) then
                !call sll_bspline_fem_solver_1d_eval_solution(knots(1:mesh_cells), electricpotential_interp)
                call sim%fsolver%evalPhi(sim%knots(1:sim%mesh_cells), sim%electricpotential_interp)
                !Scale potential to acutal values
                !                    electricpotential_interp=electricpotential_interp*&
                    !                    (sll_mass_e/sll_e_charge)*vthermal_e**2
                call electricpotential_gnuplot_inline( sim%electricpotential_interp  &
                    !+0.5_f64*knots(1:mesh_cells)&
                    ,sim%knots(1:sim%mesh_cells) )

            endif
            if ( coll_rank==0) then
                !Stop if Error ist catastrophal
                if ( (sim%kineticenergy(timestep)+sim%fieldenergy(timestep) -(sim%kineticenergy(1)+sim%fieldenergy(1)))&
                        /(sim%kineticenergy(1)+sim%fieldenergy(1)) > 10000.0_f64) then
                    !                    print *, "Its Over, Giving Up!"
                    !                    stop
                endif
            endif


            if ((sim%gnuplot_inline_output .eqv. .FALSE.) .AND. coll_rank==0) then
                print *, timestep, sim%fieldenergy(timestep),sim%kineticenergy(timestep),sim%inhom_var(timestep)

            endif
            sim%pushed_species=1
            selectcase (ppusher_int)
                case(SLL_PIC1D_PPUSHER_RK4)
                    call  sll_pic_1d_rungekutta4_step_array( sim%species(1)%particle%dx, &
                        sim%species(1)%particle%vx ,sim%tstepw, time)
                case(SLL_PIC1D_PPUSHER_VERLET)
                    call sll_pic_1d_Verlet_scheme(  sim%species(1:num_species) ,sim%tstepw, time)
                case(SLL_PIC1D_PPUSHER_EULER)
                    call sll_pic_1d_explicit_euler( sim%species(1:num_species),sim%tstepw, time)
                case(SLL_PIC1D_PPUSHER_LEAPFROG_V)
                    call  sll_pic_1d_variational_leap_frog( sim%species(sim%pushed_species)%particle%dx, &
                        sim%species(sim%pushed_species)%particle%vx ,sim%tstepw)
                case(SLL_PIC1D_PPUSHER_LEAPFROG)
                    call  sll_pic_1d_leap_frog( sim%species(sim%pushed_species)%particle%dx, &
                        sim%species(sim%pushed_species)%particle%vx ,sim%tstepw)
                case(SLL_PIC1D_PPUSHER_RK2)
                    call sll_pic_1d_rungekutta2( sim%species(sim%pushed_species)%particle%dx, &
                        sim%species(sim%pushed_species)%particle%vx ,sim%tstepw, time)
                case(SLL_PIC1D_PPUSHER_RK3)
                    call sll_pic_1d_rungekutta3( sim%species(sim%pushed_species)%particle%dx, &
                        sim%species(sim%pushed_species)%particle%vx ,sim%tstepw, time)
                case(SLL_PIC1D_PPUSHER_HEUN)
                    call sll_pic_1d_heun( sim%species(sim%pushed_species)%particle%dx, &
                        sim%species(sim%pushed_species)%particle%vx ,sim%tstepw, time)
                case(SLL_PIC1D_PPUSHER_MERSON)
                    call sll_pic_1d_merson4( sim%species(sim%pushed_species)%particle%dx, &
                        sim%species(sim%pushed_species)%particle%vx ,sim%tstepw, time, &
                        sim%species(sim%pushed_species)%particle%weight      )


                case(SLL_PIC1D_PPUSHER_NONE)
                    print *, "No Particle pusher choosen"
                    !                case(SLL_PIC1D_PPUSHER_SHIFT)
                    !                    particleposition=particleposition+timestepwidth*(interval_b-interval_a)
                    !                    !particleposition=sll_pic1d_ensure_periodicity( particleposition, &
                        !                     !   interval_a, interval_b)
                    !                     call sll_pic1d_ensure_boundary_conditions(particleposition, particlespeed)
                    !
                    !                    !call sll_bspline_fem_solver_1d_solve(particleposition)
                    !                    call sll_pic_1d_solve_qn(particleposition)

            end select
     
            if (coll_rank==0 .AND. timestep==1) then
                call sll_set_time_mark(tstop)
                !Calculate remaining Time
                print *, "Remaining Time: " , (sll_time_elapsed_between(tstart,tstop)/timestep)*real(sim%tsteps-timestep,i64)
            endif

        enddo
        




sim%kineticenergy(timestep)=sll_pic1d_calc_kineticenergy( sim%species(1:num_species) )
        sim%fieldenergy(timestep)=sll_pic1d_calc_fieldenergy(sim%species(1:num_species))
        sim%impulse(timestep)=sll_pic1d_calc_impulse(sim%species(1:num_species))
        sim%thermal_velocity_estimate(timestep)=sll_pic1d_calc_thermal_velocity(sim%species(1)%particle%vx,sim%species(1)%particle%weight)
        sim%particleweight_mean(timestep)=sum(sim%species(1)%particle%weight)/size(sim%species(1)%particle)
        sim%particleweight_var(timestep)=sum(sim%species(1)%particle%weight**2)/size(sim%species(1)%particle)-sim%particleweight_mean(timestep)**2

        close(25)

        !Save Results to file
        call sll_pic1d_write_result("pic1dresult", sim%kineticenergy, sim%fieldenergy, sim%impulse, &
            sim%particleweight_mean, sim%particleweight_var, sim%push_error_mean, sim%push_error_var )


        if (coll_rank==0) then
            call sll_set_time_mark(tstop)
            !Calculate remaining Time
            print *, "Overall Time: " , sll_time_elapsed_between(tstart,tstop)
            print *, "Particles Pushed/s: " , (sim%nmark*coll_size*sim%tsteps)/sll_time_elapsed_between(tstart,tstop)
        endif
        if (coll_rank==0)then
         call det_landau_damping((/ ( timestep*timestepwidth, timestep = 0, sim%tsteps) /),sim%fieldenergy)

        !call sll_bspline_fem_solver_1d_destroy
        call sim%fsolver%delete()












  end subroutine run_pic_1d
  
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
      case( "scenario" )
        select case( enum_string )
          case("landau");  enum_int = SLL_PIC1D_TESTCASE_LANDAU
          case("ionbeam"); enum_int = SLL_PIC1D_TESTCASE_IONBEAM
          case("quiet");   enum_int = SLL_PIC1D_TESTCASE_QUIET
          case("bump");    enum_int = SLL_PIC1D_TESTCASE_BUMPONTAIL
        end select

    end select
    

    
  end function match_enumeration

end module sll_module_simulation_pic1d1v_vp_periodic












