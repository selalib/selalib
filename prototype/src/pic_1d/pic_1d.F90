module sll_pic_1d_Class
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_utilities.h"

    use sll_constants
    !use sll_logical_meshes
    !use sll_poisson_1d_periodic
    !use sll_poisson_solvers
    use sll_collective !Parallel operations

    !  use hdf5, only: HID_T,HSIZE_T,HSSIZE_T
    !  ! #ifdef HDF5_PARALLEL
    !          use sll_hdf5_io_parallel
    !        use sll_xdmf_parallel
    ! #endif

    use sll_arbitrary_degree_splines
    use gauss_legendre_integration
    !use sll_bspline_fem_solver_1d
    use sll_pic_1d_quasi_neutral_solver
    use sll_visu_pic !Visualization with gnuplot
    use pic_1d_particle_loading !should be integrated here
    use sll_timer
    use pic_postprocessing

    implicit none

    integer :: ierr
    integer, private :: i, j


    !    type particle
    !        sll_real64, DIMENSION(:), allocatable:: particleposition
    !        sll_real64, DIMENSION(:), allocatable :: particlespeed
    !        sll_real64, DIMENSION(:), allocatable :: particleweight          ! cf. y_k in delta f and c_k in full f pic
    !
    !    endtype


    type(sll_particle_1d_group),dimension(10) :: species
    !1d1v data for electrostatic PIC
    sll_real64, DIMENSION(:), allocatable:: particleposition
    sll_real64, DIMENSION(:), allocatable :: particlespeed
    !sll_real64, DIMENSION(:), allocatable :: particleweight          ! cf. y_k in delta f and c_k in full f pic
    !sll_real64, DIMENSION(:), allocatable :: particleweight_constant !cf. c_k
    sll_real64 ::  particle_qm !<mass to electron mass ratio, with the sign of the charge
    sll_real64, DIMENSION(:), allocatable:: steadyparticleposition !Steady non-moving objects aka Ions

    sll_int32 :: particle_pusher, poisson_solver

    integer, private :: timestep

    !Data to be collected throughout a run
    sll_real64, dimension(:), allocatable :: fieldenergy, kineticenergy, impulse, &
        thermal_velocity_estimate, particleweight_mean ,particleweight_var, inhom_var ,&
        push_error_mean, push_error_var

    !Initial offets
    sll_real64 :: kineticenergy_offset, impulse_offset

    !sll_real64 :: initial_fieldenergy , initial_kineticenergy , initial_total_energy
    sll_real64 :: initial_fem_inhom , fem_inhom
    sll_int64    ::  fastest_particle_idx
    sll_real64 , DIMENSION(:), allocatable:: eval_solution


    sll_int32 :: timesteps = 100
    sll_real64 :: timestepwidth=0.1_f64
    sll_int32 :: nparticles = 1000  !GLOBAL, marker particles
    sll_int32 :: particles_system= 1 !Number of the particles in the system
    LOGICAL :: gnuplot_inline_output=.FALSE. !write gnuplot output during the simulation
    !sll_real64,parameter :: plasma_frequency = 1.0_f64 != sqrt( sll_e_charge**2/sll_e_mass/sll_epsilon_0)
    !sll_real64 ,parameter:: thermal_velocity = 1.0_f64

    !Mesh parameters
    sll_int32 :: mesh_cells = 2**8
    sll_int32                   :: spline_degree=3

    sll_real :: interval_a=0, interval_b=4.0_f64*sll_pi
    sll_real64 :: mesh_delta_t
    sll_real64, dimension(:), allocatable, private   :: knots
    sll_int32                  :: num_pts

    sll_int32 :: pushed_species

    class(pic_1d_quasi_neutral_solver), pointer :: qnsolver  !<Quasi neutral solver

    integer, private ::  phasespace_file_id=-1

    CHARACTER(LEN=255) :: root_path

    !type(sll_logical_mesh_1d), pointer :: mesh1d
    !    type(poisson_1d_periodic), pointer :: poisson1dsolver


    sll_real64, dimension(:), allocatable ::electricpotential_interp
    !For parallelization
    sll_int32, private :: coll_rank, coll_size

    ! mesh1d=> new_logical_mesh_1d(mesh_cells, interval_a,interval_b)
    !call sll_display(mesh1d)

    !Definitions for different particle pushers, enumerator
    sll_int32, parameter :: SLL_PIC1D_PPUSHER_NONE       = 0
    sll_int32, parameter :: SLL_PIC1D_PPUSHER_EULER   = 1    !Explicit EULER
    sll_int32, parameter :: SLL_PIC1D_PPUSHER_VERLET   = 2
    sll_int32, parameter :: SLL_PIC1D_PPUSHER_RK2       = 3 !Runge Kutta 2
    sll_int32, parameter :: SLL_PIC1D_PPUSHER_RK4       = 4   !Runge Kutta 4
    sll_int32, parameter :: SLL_PIC1D_PPUSHER_HEUN       = 5
    sll_int32, parameter :: SLL_PIC1D_PPUSHER_SHIFT       = 6 !Shift by velocity
    sll_int32, parameter :: SLL_PIC1D_PPUSHER_LEAPFROG       = 7
    sll_int32, parameter :: SLL_PIC1D_PPUSHER_LEAPFROG_V       = 8 !Variational Leapfrog
    sll_int32, parameter :: SLL_PIC1D_PPUSHER_RK3        =9     ! Runge Kutta 3
    sll_int32, parameter ::  SLL_PIC1D_PPUSHER_MERSON =10 !Merson 3-5 with integrated Err estimate
    sll_int32, parameter :: SLL_PIC1D_FULLF       = 1
    sll_int32, parameter :: SLL_PIC1D_DELTAF       = 2


    !Interface for one dimensional function for right hand side
    abstract interface
        function sll_pic_1d_electric_field_external(x,t) result(E)
            use sll_working_precision
            sll_real64, dimension(:),intent(in) :: x !<Position
            sll_real64, intent(in)  :: t  !<Time
            sll_real64, dimension(size(x)) :: E
        endfunction
    endinterface

    procedure(sll_pic_1d_electric_field_external), pointer:: Eex !<External (time dependent) electric field


contains

    !<
    subroutine new_sll_pic_1d(mesh_cells_user, spline_degree_user, numberofparticles_user,&
            timesteps_user, timestepwidth_user,particle_pusher_user, psolver_user )
        implicit none
        sll_int32, intent(in) :: mesh_cells_user
        sll_int32 , intent(in):: spline_degree_user
        sll_int32 , intent(in):: numberofparticles_user
        sll_int32, intent(in) :: timesteps_user
        sll_real64 , intent(in):: timestepwidth_user
        sll_int32 :: idx
        !character(len=255), intent(in):: particle_pusher_user
        sll_int32, intent(in):: particle_pusher_user , psolver_user
        SLL_ASSERT( is_power_of_two( int( mesh_cells,i64)))
        mesh_cells=mesh_cells_user
        spline_degree=spline_degree_user
        nparticles=numberofparticles_user
        timesteps=timesteps_user
        timestepwidth=timestepwidth_user
        !particle_pusher=trim(particle_pusher_user)
        particle_pusher=particle_pusher_user
        poisson_solver=psolver_user
        !############## PARALLEL ##################


        !really global variables
        coll_rank = sll_get_collective_rank( sll_world_collective )
        coll_size = sll_get_collective_size( sll_world_collective )
        call sll_collective_barrier(sll_world_collective)


        !The Mesh is needed on every Node, it is global
        mesh_delta_t = (interval_b - interval_a)/mesh_cells
        SLL_CLEAR_ALLOCATE(knots(1:mesh_cells+1),ierr)
        SLL_CLEAR_ALLOCATE(electricpotential_interp(1:mesh_cells),ierr)



        knots(1)=interval_a
        do i=2,size(knots)
            knots(i)=  knots(1) + (i-1)*mesh_delta_t*1.0_f64
        enddo
        knots(size(knots))=interval_b
        SLL_ASSERT(knots(size(knots))==interval_b)


        if (coll_rank==0) then
            print *, "Size of MPI-Collective: ", coll_size
        endif

        call sll_collective_barrier(sll_world_collective)


        !Set up Quasineutral solver
        selectcase(pic1d_testcase)
            case(SLL_PIC1D_TESTCASE_IONBEAM)
                qnsolver=>new_pic_1d_quasi_neutral_solver(interval_a, interval_b, spline_degree, &
                    mesh_cells, poisson_solver, sll_world_collective, SLL_DIRICHLET)
            case default
                qnsolver=>new_pic_1d_quasi_neutral_solver(interval_a, interval_b, spline_degree, &
                    mesh_cells, poisson_solver, sll_world_collective,SLL_PERIODIC)
        endselect

        !electric_field_external=>sll_pic_1d_electric_field
        Eex=>pic_1d_electric_field_external


        !Check Marker distribution on Cores
        if ( coll_rank==0 .and. mod(nparticles, coll_size)/=0) then
            print *, "Number of Markers per core: ", nparticles/(coll_size*1.0_f64)
            print *, "Choose appropriate number of markers!"
            stop
        endif

    endsubroutine new_sll_pic_1d

    !<Destructor
    subroutine destroy_sll_pic_1d()
        call qnsolver%delete()
        call sll_halt_collective()

    endsubroutine

    !TODO fix for delta f
    subroutine sll_pic1d_write_phasespace(timestep, p_species)
        sll_int32, intent(in) :: timestep
        type( sll_particle_1d_group), dimension(:), intent(in) :: p_species

        sll_int32:: grid=100
        if (coll_rank==0) then
            !            if (phasespace_file_id==-1) then
            !                        call sll_new_file_id(phasespace_file_id, ierr)
            !                    open(phasespace_file_id, file = trim(root_path)//'phasespace.dat')
            !            endif
            !
            !            write(phasespace_file_id,*)  particleposition, particlespeed
            ! if (phasespace_file_id==-1) then
            SLL_ASSERT(    minval(particleposition) >=interval_a)
            SLL_ASSERT(    maxval(particleposition) <=interval_b)
            grid=floor(sqrt(nparticles/10.0_f64))+1
            call distribution_xdmf(trim(root_path)//"phasespace", p_species(1)%particle%dx, p_species(1)%particle%vx, &
                p_species(1)%particle%weight,  &
                interval_a-1.0_f64, interval_b+1.0_f64, grid, &
                minval(p_species(1)%particle%vx), maxval(p_species(1)%particle%vx), grid, timestep)
            !          phasespace_file_id=0
            !endif
        endif
    endsubroutine



    !<Calculates Error estimates for different Testcases using preknown analytical
    !<solutions. Call after the first field solve.
    subroutine sll_pic_1d_initial_error_estimates()
        sll_real64 :: num_err_mean, num_err_noise_max, num_err_noise_l2, num_err_seminorm
        sll_real64, dimension(size(knots)-1) :: analytical_solution
        sll_real64 :: flag=1.0_f64
        selectcase(pic1d_testcase)
            case(SLL_PIC1D_TESTCASE_LANDAU)
                !The first field solve is a Monte Carlo estimator
                if (coll_rank==0) then
                    !analytical_solution=landau_alpha/((interval_b-interval_a)*landau_mode)*&
                            if (landau_alpha/=0) then
                        analytical_solution=landau_alpha/(landau_mode)*&
                            sin(landau_mode*(knots(1:mesh_cells)+(sll_pi/6)*(knots(2)-knots(1))))
                        num_err_noise_max= maxval(abs((eval_solution(1:mesh_cells)-analytical_solution)))/&
                            maxval(abs(analytical_solution))
                        num_err_noise_l2=sum(((eval_solution(1:mesh_cells)-analytical_solution))**2)/&
                            sum(analytical_solution**2)
                    else
                        num_err_noise_max= maxval(abs((eval_solution(1:mesh_cells))))
                    endif


                    print *, "Error in MC estimate of E-field (Max):", num_err_noise_max, "(L2):" ,num_err_noise_l2
                endif

                !num_err_seminorm=landau_alpha**2/(2*(interval_b-interval_a)*landau_mode**2)
                if (landau_alpha/=0) then
                    num_err_seminorm=landau_alpha**2/(2*landau_mode**2)/2.0_f64 !*(interval_b-interval_a)/2.0_f64
                    num_err_seminorm=(sll_pic1d_calc_fieldenergy(species(1:num_species))-num_err_seminorm)/num_err_seminorm
                else
                    num_err_seminorm=sll_pic1d_calc_fieldenergy(species(1:num_species))
                endif
                if (coll_rank==0) print *, "Error in MC estimate of E-Energy: (rel.)", num_err_seminorm
        endselect

    endsubroutine

    subroutine sll_pic_1d_run(gnuplot_inline_output_user)
        logical, intent(in), optional ::gnuplot_inline_output_user
        sll_real64, dimension(size(knots)-1) :: analytical_solution
        sll_int32 :: idx,jdx,kdx
        sll_real64 :: time
        !Timers
        type(sll_time_mark)  :: tstart, tstop

        gnuplot_inline_output=gnuplot_inline_output_user
        !For Parallelization we take the number of particles nparticles and distribute these
        !particles on the different cores
        !this also means we only plot the particles on one core
        !energy calculation and everything has to be paralellized
        nparticles=nparticles/coll_size
        print*, "#Core ", coll_rank, " handles particles", coll_rank*nparticles +1, "-", (coll_rank+1)*nparticles
        call sll_collective_barrier(sll_world_collective)
        if (coll_rank==0) print*, "#Total Number of particles: ", nparticles*coll_size


        !Scalar values to be recorded
        SLL_CLEAR_ALLOCATE(fieldenergy(1:timesteps+1), ierr)
        SLL_CLEAR_ALLOCATE(kineticenergy(1:timesteps+1), ierr)
        SLL_CLEAR_ALLOCATE(impulse(1:timesteps+1), ierr)
        SLL_CLEAR_ALLOCATE(thermal_velocity_estimate(1:timesteps+1), ierr)
        SLL_CLEAR_ALLOCATE(particleweight_mean(1:timesteps+1), ierr)
        SLL_CLEAR_ALLOCATE(particleweight_var(1:timesteps+1), ierr)
        SLL_CLEAR_ALLOCATE(push_error_mean(1:timesteps+1), ierr)
        SLL_CLEAR_ALLOCATE(push_error_var(1:timesteps+1), ierr)
        SLL_CLEAR_ALLOCATE(inhom_var(1:timesteps+1), ierr)


        SLL_CLEAR_ALLOCATE( steadyparticleposition(1:nparticles),ierr)
        SLL_CLEAR_ALLOCATE(particleposition(1:nparticles),ierr)
        SLL_CLEAR_ALLOCATE(particlespeed(1:nparticles),ierr)

        !SLL_CLEAR_ALLOCATE(particleweight(1:nparticles),ierr)
        !SLL_CLEAR_ALLOCATE(particleweight_constant(1:nparticles),ierr)

        SLL_CLEAR_ALLOCATE(eval_solution(1:size(knots)),ierr)


        if (coll_rank==0) write(*,*) "#PIC1D: Loading particles..."


        call load_particle_species (nparticles, interval_a, interval_b, species)

        !        call load_particles (nparticles, interval_a, interval_b,&
            !            steadyparticleposition, particleposition, &
            !            particlespeed, particleweight, particleweight_constant, particle_qm)
        !        particleposition=species(1)%particle%dx
        !        particlespeed=species(1)%particle%vx
        !        particleweight=species(1)%particle%weight
        !        particleweight_constant =species(1)%particle%weight_const
        particle_qm=species(1)%qm

        call sll_collective_barrier(sll_world_collective)
        if (coll_rank==0) write(*,*) "#PIC1D: Particles loaded..."

        !Determine index of fastest particle
        fastest_particle_idx=MAXLOC(particlespeed, dim=1)


        !Do precalculations, especially for the FEM solver
        !call sll_bspline_fem_solver_1d_initialize(knots, spline_degree, steadyparticleposition, -(sll_e_charge/sll_epsilon_0) )
        if (coll_rank==0) print *, "#Number of FEM-Cells ", mesh_cells


        selectcase(pic1d_testcase)
            case(SLL_PIC1D_TESTCASE_QUIET)
                call qnsolver%set_ions_constant(0.0_f64)
                particle_qm=-1.0_f64
                !call sll_bspline_fem_solver_1d_set_inhomogenity_constant(0.0_f64 )
            case(SLL_PIC1D_TESTCASE_LANDAU)
                !load constant ion background
                call qnsolver%set_ions_constant(0.0_f64)
                !call sll_bspline_fem_solver_1d_set_inhomogenity_constant(0.0_f64 )
            case(SLL_PIC1D_TESTCASE_IONBEAM)
                call qnsolver%set_ions_constant(0.0_f64)
                !call qnsolver%set_ions_constant_particles(steadyparticleposition)
            case(SLL_PIC1D_TESTCASE_IONBEAM_ELECTRONS)
                call qnsolver%set_ions_constant_particles(steadyparticleposition)
            case default
                particle_qm=-1.0_f64 !By default we simulate Electrons
                call qnsolver%set_ions_constant(0.0_f64)
        endselect


        !Delta F Method, initialize the control variate
        if (enable_deltaf .eqv. .TRUE.) then
            !The Control Variate, as it is constant over time gos into the constant inhomogenity
            !In the constant electron case this should be normalized to zero
            !But since the fem solver ignores the constant ion background anyway
            !this is going to be ok.
            !Also this should be done analytically and not like that as an Monte Carlo estimate
            !Here suppose the control variate is the initial state by default
            call  qnsolver%set_bg_particles(species(1)%particle%dx, &
                sign(1.0_f64, species(1)%qm)*species(1)%particle%weight_const)

            kineticenergy_offset=sll_pic1d_calc_kineticenergy_offset(species(1:num_species))
            impulse_offset=sll_pic1d_calc_impulse_offset(species(1:num_species))
        else
            !Here should be a quick deactivate since every code is primarily deltaf
            kineticenergy_offset=0.0_f64
            impulse_offset=0.0_f64
        endif




        !Start with PIC Method
        !Do the initial solve of the field
        !particleposition=sll_pic1d_ensure_periodicity(particleposition,  interval_a, interval_b)
        !call sll_pic1d_ensure_boundary_conditions(species(1)%particle%dx, species(1)%particle%vx)
        call sll_pic1d_ensure_boundary_conditions_species(species(1:num_species)  )


        !Initial Solve
        call qnsolver%reset_particles()
        call qnsolver%set_species(species(1:num_species))
        call qnsolver%solve()

        print *, "#Initial field solve"


        print *, "#Test initial field, and loading condition for landau damping"
        !call sll_bspline_fem_solver_1d_eval_solution_derivative(knots(1:mesh_cells)+ (sll_pi/6)*(knots(2)-knots(1)), eval_solution(1:mesh_cells))
        call qnsolver%evalE(knots(1:mesh_cells)+ (sll_pi/6)*(knots(2)-knots(1)), eval_solution(1:mesh_cells))
        call  sll_pic_1d_initial_error_estimates()


        !Information on initial field

        print *, "#Initial Field average: " , sum( qnsolver%getPotential())
        print *, "#Initial Field variance: " , sum((sum(  qnsolver%getPotential())/mesh_cells &
            -  qnsolver%getPotential())**2)/mesh_cells



        !---------------------------Initial Plots------------------------------------------------------------
        if (coll_rank==0) then



            open(unit=20, file=trim(root_path)//"initial_field.txt")
            write (20,*) eval_solution(1:mesh_cells)
            close(20)
            !
            !            !Particle density
            !            open(unit=20, file=trim(root_path)//"initial_electrondensity.txt")
            !
            !            write (20,*) bspline_particle_density(particleposition)
            !            close(20)
            !
            !            !Particle density
            !            open(unit=20, file=trim(root_path)//"initial_iondensity.txt")
            !            write (20,*) bspline_particle_density(steadyparticleposition)
            !            close(20)

            !            open(unit=20, file=trim(root_path)//"initial_phasespace.dat")
            !            do idx=1, size(particleposition)
            !                write (20,*) particleposition(idx), particlespeed(idx), particleweight(idx)
            !            enddo
            !            close(20)

            !            open(unit=20, file=trim(root_path)//"initial_electronposition.txt")
            !            do idx=1, size(particleposition)
            !                write (20,*) particleposition(idx)
            !            enddo
            !            close(20)
            !
            !            open(unit=20, file=trim(root_path)//"initial_electronspeed.txt")
            !            do idx=1, size(particleposition)
            !                write (20,*) particlespeed(idx)
            !            enddo
            !            close(20)
            !        print *, sum(abs(particleposition - sll_pic1d_ensure_periodicity( particleposition + 1.0D-3, &
                !            interval_a, interval_b)))

            !                call xv_particles_center_gnuplot('/tmp/initial_phasedensity', particleposition, particlespeed, &
                !                            interval_a, interval_b,&
                !                            -1.1_f64*maxval(particlespeed), 1.1_f64*maxval(particlespeed),&
                !                            ierr, 1.0_f64*timestep )


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
                    do idx=1, size(species(kdx)%particle)
                        write (20,*) species(kdx)%particle(idx)%dx, species(kdx)%particle(idx)%vx, species(kdx)%particle(idx)%weight
                    enddo
                enddo
                close(20)
            endif
        enddo

        !----------------------------------------------------------------------------------------------------------

        !open file

        call sll_set_time_mark(tstart)
        !Open HDF5 file for parallel write


        !---------------------------### Main Loop ###----------------------------------------------------
        time=0.0_f64
        do timestep=1, timesteps
            time=timestepwidth*timestepwidth*(timestep-1)

            call sll_pic1d_write_phasespace(timestep,species(1:num_species))

            !print *, sum(abs((sll_bspline_fem_solver_1d_get_inhomogenity())))
            !            call xv_particles_center_gnuplot('/tmp/phasedensity', particleposition, particlespeed, &
                !                interval_a, interval_b,&
                !                -1.1_f64*maxval(particlespeed), 1.1_f64*maxval(particlespeed),&
                !                ierr, 1.0_f64*timestep )


            if (nparticles*coll_size<=10000 .AND. (gnuplot_inline_output .eqv. .true.) .AND. (coll_rank==0) .AND. mod(timestep-1,timesteps/100)==0) then
                call particles_center_gnuplot_inline(pic_1d_allparticlepos(species(1:num_species)), &
                    pic_1d_allparticlev(species(1:num_species)), &
                    interval_a, interval_b,&
                    -1.1_f64*maxval(pic_1d_allparticlev(species(1:num_species))),&
                    1.1_f64*maxval( pic_1d_allparticlev(species(1:num_species))),&
                    time)

            end if


            particleweight_mean(timestep)=sum(species(1)%particle%weight)/size(species(1)%particle%weight)
            particleweight_var(timestep)=sum(species(1)%particle%weight**2)/size(species(1)%particle)-particleweight_mean(timestep)**2


            kineticenergy(timestep)=sll_pic1d_calc_kineticenergy( species(1:num_species) )
            fieldenergy(timestep)=sll_pic1d_calc_fieldenergy(species(1:num_species))
            impulse(timestep)=sll_pic1d_calc_impulse(species(1:num_species))

            thermal_velocity_estimate(timestep)=sll_pic1d_calc_thermal_velocity(species(1)%particle%vx,species(1)%particle%weight)
            inhom_var(timestep)=qnsolver%calc_variance_rhs()


            if ( (gnuplot_inline_output.eqv. .true.) .AND. coll_rank==0 .AND. mod(timestep-1,timesteps/100)==0  ) then
                call energies_electrostatic_gnuplot_inline(kineticenergy(1:timestep), fieldenergy(1:timestep), impulse(1:timestep),timestepwidth)
            else

            endif


            if ( (gnuplot_inline_output.eqv. .true.) .AND. coll_rank==0 .AND. mod(timestep-1,timesteps/100)==0 ) then
                !call sll_bspline_fem_solver_1d_eval_solution(knots(1:mesh_cells), electricpotential_interp)
                call qnsolver%evalPhi(knots(1:mesh_cells), electricpotential_interp)
                !Scale potential to acutal values
                !                    electricpotential_interp=electricpotential_interp*&
                    !                    (sll_mass_e/sll_e_charge)*vthermal_e**2
                call electricpotential_gnuplot_inline( electricpotential_interp  &
                    !+0.5_f64*knots(1:mesh_cells)&
                    ,knots(1:mesh_cells) )

            endif

            if ( coll_rank==0) then
                !Stop if Error ist catastrophal
                if ( (kineticenergy(timestep)+fieldenergy(timestep) -(kineticenergy(1)+fieldenergy(1)))&
                        /(kineticenergy(1)+fieldenergy(1)) > 10000.0_f64) then
                    !                    print *, "Its Over, Giving Up!"
                    !                    stop
                endif
            endif


            if ((gnuplot_inline_output .eqv. .FALSE.) .AND. coll_rank==0) then
                print *, timestep, fieldenergy(timestep),kineticenergy(timestep),inhom_var(timestep)

            endif


            !Push particles by species
            ! particle_mask=(/I,(I=1,nparticles)/)
            !do pushed_species=1,num_species

            !   particle_qm=species(pushed_species)%qm

            pushed_species=1
            selectcase (particle_pusher)
                case(SLL_PIC1D_PPUSHER_RK4)
                    call  sll_pic_1d_rungekutta4_step_array( species(1)%particle%dx, &
                        species(1)%particle%vx ,timestepwidth, time)
                case(SLL_PIC1D_PPUSHER_VERLET)
                    call sll_pic_1d_Verlet_scheme(  species(1:num_species) ,timestepwidth, time)
                case(SLL_PIC1D_PPUSHER_EULER)
                    call sll_pic_1d_explicit_euler( species(1:num_species),timestepwidth, time)
                case(SLL_PIC1D_PPUSHER_LEAPFROG_V)
                    call  sll_pic_1d_variational_leap_frog( species(pushed_species)%particle%dx, &
                        species(pushed_species)%particle%vx ,timestepwidth)
                case(SLL_PIC1D_PPUSHER_LEAPFROG)
                    call  sll_pic_1d_leap_frog( species(pushed_species)%particle%dx, &
                        species(pushed_species)%particle%vx ,timestepwidth)
                case(SLL_PIC1D_PPUSHER_RK2)
                    call sll_pic_1d_rungekutta2( species(pushed_species)%particle%dx, &
                        species(pushed_species)%particle%vx ,timestepwidth, time)
                case(SLL_PIC1D_PPUSHER_RK3)
                    call sll_pic_1d_rungekutta3( species(pushed_species)%particle%dx, &
                        species(pushed_species)%particle%vx ,timestepwidth, time)
                case(SLL_PIC1D_PPUSHER_HEUN)
                    call sll_pic_1d_heun( species(pushed_species)%particle%dx, &
                        species(pushed_species)%particle%vx ,timestepwidth, time)
                case(SLL_PIC1D_PPUSHER_MERSON)
                    call sll_pic_1d_merson4( species(pushed_species)%particle%dx, &
                        species(pushed_species)%particle%vx ,timestepwidth, time, &
                        species(pushed_species)%particle%weight      )


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

            ! enddo
            !Collisions
            !call sll_pic1d_collision_step()

            !Inject and remove particles
            !call sll_pic1d_injection_removing()



            if (coll_rank==0 .AND. timestep==1) then
                call sll_set_time_mark(tstop)
                !Calculate remaining Time
                print *, "Remaining Time: " , (sll_time_elapsed_between(tstart,tstop)/timestep)*real(timesteps-timestep,i64)
            endif

        enddo

        !Last update of relevant data
        kineticenergy(timestep)=sll_pic1d_calc_kineticenergy( species(1:num_species) )
        fieldenergy(timestep)=sll_pic1d_calc_fieldenergy(species(1:num_species))
        impulse(timestep)=sll_pic1d_calc_impulse(species(1:num_species))
        thermal_velocity_estimate(timestep)=sll_pic1d_calc_thermal_velocity(species(1)%particle%vx,species(1)%particle%weight)
        particleweight_mean(timestep)=sum(species(1)%particle%weight)/size(species(1)%particle)
        particleweight_var(timestep)=sum(species(1)%particle%weight**2)/size(species(1)%particle)-particleweight_mean(timestep)**2

        close(25)

        !Save Results to file
        call sll_pic1d_write_result("pic1dresult", kineticenergy, fieldenergy, impulse, &
            particleweight_mean, particleweight_var, push_error_mean, push_error_var )


        if (coll_rank==0) then
            call sll_set_time_mark(tstop)
            !Calculate remaining Time
            print *, "Overall Time: " , sll_time_elapsed_between(tstart,tstop)
            print *, "Particles Pushed/s: " , (nparticles*coll_size*timesteps)/sll_time_elapsed_between(tstart,tstop)
        endif
        if (coll_rank==0) call det_landau_damping((/ ( timestep*timestepwidth, timestep = 0, timesteps) /),fieldenergy)

        !call sll_bspline_fem_solver_1d_destroy
        call qnsolver%delete()

    endsubroutine

    !    subroutine sll_pic1d_collision_step()
    !      !  sll_real64, dimension(size(particleweight)) :: collision_weights
    !        sll_real64 :: ueq, Teq
    !        !       sll_real64:: feq,v
    !        sll_int32 :: idx
    !        !feq(v)=1.0_f64/sqrt(sll_kx)*exp(-0.5_f64*(v-ueq)**2/Teq)
    !        ueq=dot_product(particlespeed, particleweight)
    !        Teq=dot_product(particlespeed**2, particleweight)-ueq**2
    !
    !        !                collision_weights=particleweight&
        !            !                         /sqrt(sll_kx)*exp(-0.5_f64*(particlespeed-ueq)**2/Teq)/(interval_length)
    !        !                         initial_dist_xv(particleposition, particlespeed)
    !        !
    !
    !
    !        !        do idx=1,size(particleweight)
    !        !                collision_weights(idx)=particleweight(idx)*&
        !            !                         feq(particlespeed(idx))/&
        !            !                         initial_dist_xv((/particleposition(idx)/), (/particlespeed(idx)/))
    !        !        enddo
    !
    !
    !
    !
    !    endsubroutine


    subroutine sll_pic1d_injection_removing()
        sll_real64 :: weight
        sll_real64 :: injectrate
        sll_int32 :: numinject, idx ,jdx

        injectrate=1D-1
        weight=2.0_f64/(nparticles*coll_size*1.0_f64)
        selectcase(pic1d_testcase)
            case(SLL_PIC1D_TESTCASE_IONBEAM)
                numinject=floor(injectrate*timestepwidth*nparticles*1.0_f64)
                jdx=0
                do idx=1,size(species(3)%particle)
                    if (species(3)%particle(idx)%weight==0.0_f64) then
                        jdx=jdx+1
                        species(3)%particle(idx)%weight=weight
                    endif
                    if (jdx==numinject) then
                        !print *, "injected", numinject
                        exit
                    endif
                enddo
        endselect

    endsubroutine

    !    !<Wrapper for the field solver, here particle destruction and injection can be made
    !    !<during the push
    !    subroutine sll_pic_1d_solve_field(p_species)
    !        type( sll_particle_1d_group), dimension(:), intent(in) :: p_species
    !
    !
    !
    !
    !
    !    endsubroutine



    !<Only relevant for solving the QN in the pusher
    subroutine sll_pic_1d_solve_qn(particleposition_selected)
        sll_real64, dimension(:),intent(in) :: particleposition_selected
        !sll_real64, dimension(nparticles) :: allparticleposition
        sll_real64, dimension(nparticles) :: allparticleposition
        !allparticleposition=particleposition
        !allparticleposition( particle_mask)=particleposition_selected
        sll_int32 ::jdx
        !Add all the other species as constant
        call qnsolver%reset_particles()

        do jdx=1,num_species
            if (jdx/=pushed_species) then
                call qnsolver%add_species(species(jdx))
            endif
        enddo

        selectcase(pic1d_testcase)
            !            case(SLL_PIC1D_TESTCASE_QUIET)
            case(SLL_PIC1D_TESTCASE_LANDAU)
                !All particles, are only electrons
                !load constant ion background
                SLL_ASSERT(size(particleposition_selected)==size(species(pushed_species)%particle))

                call qnsolver%set_electrons_only_weighted(particleposition_selected, species(pushed_species)%particle%weight)
                call qnsolver%solve()
            case(SLL_PIC1D_TESTCASE_IONBEAM)
                SLL_ASSERT(size(particleposition_selected)==size(species(pushed_species)%particle))

                call qnsolver%add_particles_weighted(particleposition_selected, &
                    sign(species(pushed_species)%particle%weight,species(pushed_species)%qm) )
                call qnsolver%solve()
            case default
                SLL_ASSERT(num_species==1)
                !All particles, are only electrons
                SLL_ASSERT(size(particleposition_selected)==size(species(pushed_species)%particle))
                call qnsolver%set_electrons_only_weighted(particleposition_selected, species(pushed_species)%particle%weight)
                call qnsolver%solve()

                !load constant ion background
                !                SLL_ASSERT(size(allparticleposition)==size(particleweight))
                !                call qnsolver%set_electrons_only_weighted(allparticleposition, particleweight)
                !                call qnsolver%solve()
        endselect

    endsubroutine

    !Dummy function for solving the potential Phi
    !Input x \in [a, b] is the position in the electric field
    !Result is $\frac{DPhi}{dx}(x)$
    !This evaluates with the solution of the weak Poisson equation
    !the solution vector at the
    !    function sll_pic_1d_DPhidx(x, bspline, knots_mesh, fem_Phi_solution ) &
        !            result(DPhidx)
    !        sll_real64, dimension(:), intent(in) :: x
    !        type(arbitrary_degree_spline_1d), pointer, intent(in) :: bspline
    !        sll_real64, dimension(:), intent(in)     :: knots_mesh
    !        sll_real64, dimension(:), intent(in) :: fem_Phi_solution
    !        sll_real :: DPhidx(size(x))
    !        sll_real, dimension(:) :: x_trunc(size(x))
    !        !Account for periodic geometry
    !        x_trunc=sll_pic1d_ensure_periodicity(x,knots_mesh(1), knots_mesh(size(knots_mesh)))
    !
    !        DPhidx=bspline_basis_to_realvals(bspline, knots_mesh, fem_Phi_solution, x_trunc, b_spline_derivatives_at_x)
    !    endfunction



    !    function sll_pic_1d_solve_DPhidx(bspline, knots_mesh, fem_Phi_solution, particleposition, fem_inhomogenity_steady ) &
        !        sll_real64, dimension(:), intent(in) :: eval_points
    !        type(arbitrary_degree_spline_1d), pointer, intent(in) :: bspline
    !        sll_real64, dimension(:), intent(in)     :: knots_mesh
    !        sll_real64, dimension(:), intent(out) :: fem_Phi_solution
    !        sll_real64, dimension(:), intent(inout) particleposition
    !        sll_real :: DPhidx(size(x))
    !        sll_real64, dimension(:), intent(in) :: fem_inhomogenity_steady
    !        sll_real64, dimension(:), intent(inout) :: fem_inhomogenity
    !
    !
    !        !Account for periodic geometry
    !        x_trunc=sll_pic1d_ensure_periodicity(particleposition,knots_mesh(1), knots_mesh(size(knots_mesh)))
    !
    !        fem_inhomogenity=(sll_e_charge/sll_epsilon_0)   &
        !                *( fem_inhomogenity_steady + &
        !                interpolate_particles_bsplines( arbitrarydegree1D, knots, particleposition, b_splines_at_x))
    !
    !
    !            fem_Phi_solution=sll_bspline_fem_solve_poisson_1d(arbitrarydegree1D, knots,fem_inhomogenity)
    !
    !    endfunction




    subroutine sll_pic_1d_rungekutta4_step_array(x_0, v_0, h, t)
        sll_real64, intent(in):: h
        sll_real64, intent(in):: t
        sll_real64, dimension(:) ,intent(inout) :: x_0
        sll_real64, dimension(:) ,intent(inout) :: v_0
        sll_real64, dimension(size(x_0)) :: k_x1, k_x2,&
            k_x3, k_x4, k_v1, &
            k_v2, k_v3, k_v4, stage_DPhidx,x_1,v_1

        !x_0=sll_pic1d_ensure_periodicity(x_0,  interval_a, interval_b)
        call sll_pic1d_ensure_boundary_conditions(x_0, v_0)

        !--------------------Stage 1-------------------------------------------------
        !call sll_bspline_fem_solver_1d_solve(x_0)
        call qnsolver%evalE(x_0, stage_DPhidx)
        k_x1= h*v_0
        k_v1= h*(stage_DPhidx+Eex(x_0,t))*(-particle_qm)
        !--------------------Stage 2-------------------------------------------------
        k_x2= h*(v_0  + 0.5_f64  *k_v1)
        call sll_pic_1d_solve_qn(x_0 + 0.5_f64*k_x2)
        !call sll_pic_1d_solve_qn(x_0 + 0.5_f64*k_x1)
        call qnsolver%evalE(x_0 + 0.5_f64 *k_x1, stage_DPhidx)
        k_v2= h*( stage_DPhidx + Eex(x_0 + 0.5_f64*k_x1, t+0.5_f64*h))*(-particle_qm)
        !--------------------Stage 3-------------------------------------------------
        k_x3= h*(v_0+ 0.5_f64 * k_v2)
        call sll_pic_1d_solve_qn(x_0 + 0.5_f64 *k_x3)
        !call sll_bspline_fem_solver_1d_solve(sll_pic1d_ensure_periodicity(x_0 + (1*k_x1 + 2.0_f64*k_x3)/6.0_f64,  interval_a, interval_b))
        !call sll_pic_1d_solve_qn(x_0 + 0.5_f64 *k_x2)
        call qnsolver%evalE(x_0 + 0.5_f64 *k_x2, stage_DPhidx)
        k_v3= h*(stage_DPhidx+ Eex(x_0 + 0.5_f64*k_x2, t+0.5_f64*h)) *(-particle_qm)
        !--------------------Stage 4-------------------------------------------------
        k_x4= h*(v_0+ k_v3)
        !call sll_bspline_fem_solver_1d_solve(sll_pic1d_ensure_periodicity(x_0 + (k_x1 + 2.0_f64*k_x2 + 2.0_f64*k_x3 +k_x4  )/6.0_f64,  interval_a, interval_b))
        call sll_pic_1d_solve_qn(x_0 +  k_x4)
        !call sll_pic_1d_solve_qn(x_0 +  k_x3)
        call qnsolver%evalE(x_0 +  k_x3, stage_DPhidx)
        k_v4= h*(stage_DPhidx+ Eex(x_0 + k_x3, t+h)) *(-particle_qm)
        !Perform step---------------------------------------------------------------
        x_1= x_0 + (  k_x1 +  2.0_f64 *k_x2 +  2.0_f64*k_x3 + k_x4 )/6.0_f64
        v_1= v_0 + (  k_v1 +  2.0_f64 *k_v2 +  2.0_f64*k_v3 + k_v4)/6.0_f64




        !x_0=sll_pic1d_ensure_periodicity(x_0,  interval_a, interval_b)
        call sll_pic1d_ensure_boundary_conditions(x_1, v_1)

        call sll_pic1d_adjustweights(x_0,x_1,v_0,v_1)
        x_0=x_1
        v_0=v_1

        call sll_pic_1d_solve_qn(x_0)

    endsubroutine



    !<Merson scheme with built in error estimate for each particle
    subroutine sll_pic_1d_merson4(x_0, v_0, h, t, weight)
        sll_real64, intent(in):: h
        sll_real64, intent(in):: t
        sll_real64, dimension(:) ,intent(inout) :: x_0
        sll_real64, dimension(:) ,intent(inout) :: v_0
        sll_real64, dimension(:) ,intent(in) :: weight

        sll_real64, dimension(size(x_0)) :: k_x1, k_x2,&
            k_x3, k_x4, k_v1, &
            k_v2, k_v3, k_v4, k_v5, k_x5, stage_DPhidx,x_1,v_1, err_x, err_v

        !x_0=sll_pic1d_ensure_periodicity(x_0,  interval_a, interval_b)
        call sll_pic1d_ensure_boundary_conditions(x_0, v_0)

        !--------------------Stage 1-------------------------------------------------
        !call sll_bspline_fem_solver_1d_solve(x_0)
        call qnsolver%evalE(x_0, stage_DPhidx)
        k_x1= h*v_0
        k_v1= h*(stage_DPhidx+Eex(x_0,t))*(-particle_qm)
        !--------------------Stage 2-------------------------------------------------
        k_x2= h*(v_0  + (1.0_f64/3.0_f64)*k_v1)
        call sll_pic_1d_solve_qn(x_0 + (1.0_f64/3.0_f64)*k_x1)
        call qnsolver%evalE(x_0 + (1.0_f64/3.0_f64) *k_x1, stage_DPhidx)
        k_v2= h*( stage_DPhidx + Eex(x_0 + (1.0_f64/3.0_f64)*k_x1, t+(1.0_f64/3.0_f64)*h))*(-particle_qm)
        !--------------------Stage 3-------------------------------------------------
        k_x3= h*(v_0+ (1.0_f64/3.0_f64) * k_v2)
        call sll_pic_1d_solve_qn(x_0 +  (1.0_f64/6.0_f64)*(k_x1 + k_x2 ))

        call qnsolver%evalE(x_0 + (1.0_f64/6.0_f64)*(k_x1 + k_x2 ), stage_DPhidx)
        k_v3= h*(stage_DPhidx+ Eex(x_0 + (1.0_f64/6.0_f64)*(k_x1 + k_x2 ), &
            t+(1.0_f64/3.0_f64)*h)) *(-particle_qm)
        !--------------------Stage 4-------------------------------------------------
        k_x4= h*(v_0+ 0.5_f64*k_v3)
        call sll_pic_1d_solve_qn(x_0 +  (k_x1 + 3*k_x3)/8.0_f64  )
        call qnsolver%evalE(x_0 +  (k_x1 + 3*k_x3)/8.0_f64  , stage_DPhidx)
        k_v4= h*(stage_DPhidx+ Eex(x_0 +  (k_x1 + 3*k_x3)/8.0_f64  , t+0.5_f64*h)) *(-particle_qm)
        !--------------------Stage 5-------------------------------------------------
        k_x5= h*(v_0+ k_v4)
        call sll_pic_1d_solve_qn(x_0 +  0.5_f64*k_x1  -1.5_f64*k_x3 + 2.0_f64*k_x1)
        call qnsolver%evalE(x_0 +  0.5_f64*k_x1  -1.5_f64*k_x3 + 2.0_f64*k_x1, stage_DPhidx)
        k_v4= h*(stage_DPhidx+ Eex(x_0 +  0.5_f64*k_x1  -1.5_f64*k_x3 + 2.0_f64*k_x1, t+h)) *(-particle_qm)

        !Perform step of Order 4-----------------------------------------------------------
        x_1= x_0 + (  k_x1 +  0.0_f64 *k_x2 +  0.0_f64*k_x3 + 4.0_f64*k_x4 + k_x5 )/6.0_f64
        v_1= v_0 + (  k_v1 +  0.0_f64 *k_v2 +  0.0_f64*k_v3 + 4.0_f64*k_v4 + k_v5 )/6.0_f64

        !Perform error step of Order 3 or 5 in linear case
        err_x= x_0 + (  k_x1 +  3*k_x3 + 4*k_x4 + 2*k_x5 )/10.0_f64
        err_v= x_0 + (  k_v1 +  3*k_v3 + 4*k_v4 + 2*k_v5 )/10.0_f64
        err_x=x_1 - err_x
        err_v=v_1 - err_v

        call sll_pic_1d_calc_push_error(err_x,weight)

        call sll_pic1d_ensure_boundary_conditions(x_1, v_1)

        call sll_pic1d_adjustweights(x_0,x_1,v_0,v_1)
        x_0=x_1
        v_0=v_1

        call sll_pic_1d_solve_qn(x_0)

    endsubroutine

    subroutine sll_pic_1d_rungekutta2(x_0, v_0, h, t)
        sll_real64, intent(in):: h
        sll_real64, intent(in):: t
        sll_real64, dimension(:) ,intent(inout) :: x_0
        sll_real64, dimension(:) ,intent(inout) :: v_0
        sll_real64, dimension(size(x_0)) :: k_x1, k_x2,&
            k_v1, k_v2, stage_DPhidx,x_1,v_1

        !x_0=sll_pic1d_ensure_periodicity(x_0,  interval_a, interval_b)
        call sll_pic1d_ensure_boundary_conditions(x_0, v_0)

        !--------------------Stage 1-------------------------------------------------
        !call sll_bspline_fem_solver_1d_solve(x_0)
        call qnsolver%evalE(x_0, stage_DPhidx)
        k_x1= h*v_0
        k_v1= h*(stage_DPhidx+Eex(x_0,t)) *(-particle_qm)
        !--------------------Stage 2-------------------------------------------------
        k_x2= h*(v_0  + 0.5_f64  *k_v1)
        call sll_pic_1d_solve_qn(x_0 + 0.5_f64 *k_x1)
        call qnsolver%evalE(x_0 + 0.5_f64 *k_x1, stage_DPhidx)
        k_v2= h*( stage_DPhidx + Eex(x_0 + 0.5_f64 *k_x1,t + 0.5_f64*h))&
            *(-particle_qm)
        !Perform step---------------------------------------------------------------
        x_1= x_0 + k_x2
        v_1= v_0 + k_v2

        !        x_0=sll_pic1d_ensure_periodicity(x_0,  interval_a, interval_b)
        call sll_pic1d_ensure_boundary_conditions(x_1, v_1)

        call sll_pic1d_adjustweights(x_0,x_1,v_0,v_1)
        x_0=x_1
        v_0=v_1

        call sll_pic_1d_solve_qn(x_0)
    endsubroutine


    subroutine sll_pic_1d_rungekutta3(x_0, v_0, h, t)
        sll_real64, intent(in):: h
        sll_real64, intent(in):: t
        sll_real64, dimension(:) ,intent(inout) :: x_0
        sll_real64, dimension(:) ,intent(inout) :: v_0
        sll_real64, dimension(size(x_0)) :: k_x1, k_x2,&
            k_x3, k_v3, k_v1,k_v2, stage_DPhidx

        !x_0=sll_pic1d_ensure_periodicity(x_0,  interval_a, interval_b)
        call sll_pic1d_ensure_boundary_conditions(x_0, v_0)

        !--------------------Stage 1-------------------------------------------------
        call qnsolver%evalE(x_0, stage_DPhidx)
        k_x1= h*v_0
        k_v1= h*(stage_DPhidx+Eex(x_0,t)) *(-particle_qm)
        !--------------------Stage 2-------------------------------------------------
        k_x2= h*(v_0  + 0.5_f64  *k_v1)
        call sll_pic_1d_solve_qn(x_0 + 0.5_f64 *k_x1)
        call qnsolver%evalE(x_0 + 0.5_f64 *k_x1, stage_DPhidx)
        k_v2= h*( stage_DPhidx + Eex(x_0 + 0.5_f64 *k_x1,t + 0.5_f64*h))&
            *(-particle_qm)
        !--------------------Stage 3-------------------------------------------------
        k_x3= h*(v_0  - k_v1 + 2.0_f64*k_v2)
        call sll_pic_1d_solve_qn(x_0 - k_x1 + 2.0_f64*k_x2)
        call qnsolver%evalE(x_0 - k_x1 + 2.0_f64*k_x2, stage_DPhidx)
        k_v3= h*( stage_DPhidx + Eex(x_0 - k_x1 + 2.0_f64*k_x2,t + h))&
            *(-particle_qm)

        !Perform step---------------------------------------------------------------
        x_0= x_0 + (k_x1 + 4.0_f64 * k_x2 + k_x3)/6.0_f64
        v_0= v_0 + (k_v1 + 4.0_f64 * k_v2 + k_v3)/6.0_f64

        call sll_pic1d_ensure_boundary_conditions(x_0, v_0)

        call sll_pic_1d_solve_qn(x_0)
    endsubroutine


    subroutine sll_pic_1d_heun(x_0, v_0, h, t)
        sll_real64, intent(in):: h
        sll_real64, intent(in):: t
        sll_real64, dimension(:) ,intent(inout) :: x_0
        sll_real64, dimension(:) ,intent(inout) :: v_0
        sll_real64, dimension(:) :: k_x1(size(x_0)), k_x2(size(x_0)), &
            k_v1(size(x_0)), k_v2(size(x_0)),stage_DPhidx(size(x_0))

        !x_0=sll_pic1d_ensure_periodicity(x_0,  interval_a, interval_b)
        call sll_pic1d_ensure_boundary_conditions(x_0, v_0)
        !--------------------Stage 1-------------------------------------------------
        !call sll_bspline_fem_solver_1d_solve(x_0)
        call qnsolver%evalE(x_0, stage_DPhidx)
        k_x1= h*v_0
        k_v1= h*(stage_DPhidx+Eex(x_0,t))*(-particle_qm)
        !--------------------Stage 2-------------------------------------------------
        call sll_pic_1d_solve_qn(x_0 + k_x1)
        call qnsolver%evalE(x_0 + k_x1, stage_DPhidx)
        k_x2= h*(v_0  + k_v1)
        k_v2= h*( stage_DPhidx+Eex(x_0 +k_x1,t+h))*(-particle_qm)

        !Perform step---------------------------------------------------------------
        x_0= x_0 + 0.5_f64 *(k_x1+k_x2   )
        v_0= v_0 + 0.5_f64 *(k_v1+k_v2  )

        !x_0=sll_pic1d_ensure_periodicity(x_0,  interval_a, interval_b)
        call sll_pic1d_ensure_boundary_conditions(x_0, v_0)
        call sll_pic_1d_solve_qn(x_0)
    endsubroutine

    subroutine sll_pic_1d_explicit_euler(species_0, h, t)
        type(sll_particle_1d_group), dimension(:), intent(inout) :: species_0
        sll_int32 :: jdx, num_part
        sll_real64, intent(in):: h
        sll_real64, intent(in):: t
        type(sll_particle_1d_group), dimension(size(species_0)) :: species_1

        sll_real64 , dimension(:),allocatable :: stage_DPhidx

        !----ALLOCATE STAGES-----------------------------------------------
        call sll_pic_1d_copy_species(species_0, species_1)

        !--------------------Stage 1-------------------------------------------------
        do jdx=1,size(species_0)
            num_part=size(species_0(jdx)%particle)
            SLL_ALLOCATE(stage_DPhidx(num_part),ierr)

            call qnsolver%evalE(species_0(jdx)%particle%dx, stage_DPhidx)

            species_1(jdx)%particle%dx=species_0(jdx)%particle%dx + h*species_0(jdx)%particle%vx
            species_1(jdx)%particle%vx=species_0(jdx)%particle%vx + h*(stage_DPhidx + &
                Eex( species_0(jdx)%particle%dx,t) )*(-species_0(jdx)%qm)

            SLL_DEALLOCATE_ARRAY(stage_DPhidx,ierr)

        enddo

        call sll_pic1d_ensure_boundary_conditions_species(species_1)

        call sll_pic1d_adjustweights_advection_species(species_0, species_1)

        call sll_pic_1d_copy_species(species_1, species_0)

        call qnsolver%set_species(species_0)
        call qnsolver%solve()

    endsubroutine


    subroutine sll_pic_1d_copy_species(original_species, species )
        type(sll_particle_1d_group), dimension(:), intent(in) :: original_species
        type(sll_particle_1d_group), dimension(:), intent(inout) :: species
        sll_int32 :: num_species, numpart, jdx

        num_species=size(original_species)
        SLL_ASSERT(num_species==size(species))

        do jdx=1,num_species
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
    endsubroutine


    !    subroutine sll_pic_1d_explicit_euler(x_0, v_0, h, t, particle_qm)
    !        sll_real64, intent(in):: h
    !        sll_real64, intent(in):: t
    !
    !        sll_real64, dimension(:) ,intent(inout) :: x_0
    !        sll_real64, dimension(:) ,intent(inout) :: v_0
    !        sll_real64.  dimension(:) ,intent(in), pointer :: particle_qm
    !        sll_real64, dimension(size(x_0)) :: stage_DPhidx,x_1,v_1
    !
    !        !x_0=sll_pic1d_ensure_periodicity(x_0,  interval_a, interval_b)
    !
    !        !--------------------Stage 1-------------------------------------------------
    !        !call sll_bspline_fem_solver_1d_solve(x_0)
    !        call qnsolver%evalE(x_0, stage_DPhidx)
    !        x_1=x_0 + h*v_0
    !        v_1=v_0 + h*(stage_DPhidx + Eex(x_0,t) )*(-particle_qm)
    !
    !        call sll_pic1d_ensure_boundary_conditions(x_0, v_0)
    !
    !        call sll_pic1d_adjustweights(x_0,x_1,v_0,v_1)
    !        x_0=x_1
    !        v_0=v_1
    !
    !        call sll_pic_1d_solve_qn(x_0)
    !    endsubroutine

    !<From Erics lecture notes on Vaslov Equations
    subroutine sll_pic_1d_Verlet_scheme(species_0, h, t)
        type(sll_particle_1d_group), dimension(:), intent(inout) :: species_0
        sll_int32 :: jdx, num_part
        sll_real64, intent(in):: h
        sll_real64, intent(in):: t
        type(sll_particle_1d_group), dimension(size(species_0)) :: species_1, species_05

        sll_real64 , dimension(:),allocatable :: DPhidx

        !----ALLOCATE STAGES-----------------------------------------------
        call sll_pic_1d_copy_species(species_0, species_1)
        call sll_pic_1d_copy_species(species_0, species_05)


        !--------------------Stage 1-------------------------------------------------
        do jdx=1,size(species_0)
            num_part=size(species_0(jdx)%particle)
            SLL_ALLOCATE(DPhidx(num_part),ierr)

            call qnsolver%evalE(species_0(jdx)%particle%dx, DPhidx)

            species_05(jdx)%particle%vx=species_0(jdx)%particle%vx+ &
                0.5_f64*h*(DPhidx+Eex(species_0(jdx)%particle%dx,t))*(-species_0(jdx)%qm)

            species_05(jdx)%particle%dx=species_0(jdx)%particle%dx  + h*species_05(jdx)%particle%vx

            SLL_DEALLOCATE_ARRAY(DPhidx,ierr)
        enddo

        !call sll_pic1d_adjustweights_advection_species(species_0, species_05)

        call sll_pic1d_ensure_boundary_conditions_species(species_05)
        call qnsolver%set_species(species_05)
        call qnsolver%solve()

        do jdx=1,size(species_0)
            num_part=size(species_0(jdx)%particle)
            SLL_ALLOCATE(DPhidx(num_part),ierr)

            call qnsolver%evalE(species_05(jdx)%particle%dx, DPhidx)

            species_1(jdx)%particle%dx=species_05(jdx)%particle%dx
            species_1(jdx)%particle%vx=species_05(jdx)%particle%vx + &
                0.5_f64*h*(DPhidx+Eex( &
                species_0(jdx)%particle%dx,t+h))*(-species_1(jdx)%qm)

            SLL_DEALLOCATE_ARRAY(DPhidx,ierr)
        enddo

        call sll_pic1d_ensure_boundary_conditions_species(species_1)


        call sll_pic1d_adjustweights_advection_species(species_0, species_1)
        call sll_pic_1d_copy_species(species_1, species_0)
        !sll_real64, intent(in):: h
        !        sll_real64, intent(in):: t
        !        sll_real64, dimension(:) ,intent(inout) :: x_0
        !        sll_real64, dimension(:) ,intent(inout) :: v_0
        !        sll_real64, dimension(size(x_0)) ::  DPhidx, x_1, v_1, v_05
        !
        !        !x_0=sll_pic1d_ensure_periodicity(x_0,  interval_a, interval_b)
        !
        !        call qnsolver%evalE(x_0, DPhidx)
        !        v_05=v_0+ 0.5_f64*h*(DPhidx+Eex(x_0,t))*(-particle_qm)
        !        x_1=x_0 + h*v_05
        !
        !        call sll_pic1d_ensure_boundary_conditions(x_1, v_05)
        !        call sll_pic_1d_solve_qn(x_1)
        !
        !        call qnsolver%evalE(x_1, DPhidx)
        !        v_1=v_05 + 0.5_f64*h*(DPhidx+Eex(x_1,t+h))*(-particle_qm)
        !        !Commit Push
        !        call sll_pic1d_adjustweights(x_0,x_1,v_0,v_1)
        !        x_0=x_1
        !        v_0=v_1

    endsubroutine

    subroutine sll_pic_1d_variational_leap_frog(x_0, v_0, h)
        sll_real64, intent(in):: h
        sll_real64, dimension(:) ,intent(inout) :: x_0
        sll_real64, dimension(:) ,intent(inout) :: v_0
        sll_real64, dimension(size(x_0)) ::  DPhidx_0, DPhidx_1, x_1, v_1


        call qnsolver%evalE(x_0, DPhidx_0)
        x_1=x_0+ h*v_0 - ((h**2)/2.0_f64) *DPhidx_0*(-particle_qm)
        x_1=sll_pic1d_ensure_periodicity(x_1,  interval_a, interval_b)
        call sll_pic1d_ensure_boundary_conditions(x_1, v_0)

        call sll_pic_1d_solve_qn(x_1)
        call qnsolver%evalE(x_1, DPhidx_1)

        v_1=v_0 + 0.5_f64*h* (DPhidx_0 + DPhidx_1)*(-particle_qm)

        !Push
        x_0=x_1
        v_0=v_1
    endsubroutine

    subroutine sll_pic_1d_leap_frog(x_0, v_0, h)
        sll_real64, intent(in):: h
        sll_real64, dimension(:) ,intent(inout) :: x_0
        sll_real64, dimension(:) ,intent(inout) :: v_0
        sll_real64, dimension(size(x_0)) ::  DPhidx_0, DPhidx_1, x_1, v_1

        x_0=sll_pic1d_ensure_periodicity(x_0,  interval_a, interval_b)
        call qnsolver%evalE(x_0, DPhidx_0)
        !DPhidx_0=0

        v_1= v_0 - 0.5_f64 *DPhidx_0*h*(-particle_qm)
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

        !x_0=sll_pic1d_ensure_periodicity(x_0,  interval_a, interval_b)
        call sll_pic1d_ensure_boundary_conditions(x_0, v_0)

        call sll_pic_1d_solve_qn(x_0)
    endsubroutine


    subroutine sll_pic_1d_calc_push_error(perror, pweight)
        sll_real64, dimension(:), intent(in) :: perror
        sll_real64, dimension(:), intent(in) :: pweight

        push_error_mean(timestep)=dot_product(perror, pweight)
        push_error_var(timestep)=dot_product(perror**2, pweight)&
            - push_error_mean(timestep)

        push_error_mean=maxval(perror)
    endsubroutine
    !    !<Nonrelativistic kinetic energy 0.5*m*v^2
    !    function sll_pic1d_calc_kineticenergy(relmass,  particlespeed , particleweight ) &
        !            result(energy)
    !        sll_real64, DIMENSION(:), intent(in) :: particlespeed
    !        sll_real64, DIMENSION(:), intent(in):: particleweight
    !        sll_real64, DIMENSION(:), intent(in) ::relmass
    !
    !        sll_real64 :: energy
    !        SLL_ASSERT(size(particlespeed)==size(particleweight))
    !
    !        energy=0.5_f64*dot_product(particleweight*relmass,particlespeed**2)
    !        call sll_collective_globalsum(sll_world_collective, energy, 0)
    !    endfunction


    !<Nonrelativistic kinetic energy 0.5*m*v^2
    function sll_pic1d_calc_kineticenergy(p_species ) &
            result(energy)
        type( sll_particle_1d_group), dimension(:), intent(in) :: p_species
        sll_int32 :: idx
        sll_int32 :: num_sp
        sll_real64 :: energy
        num_sp=size(p_species)
        energy=0
       do idx=1,num_sp
            energy=energy + &
            sll_pic1d_calc_kineticenergy_weighted(p_species(idx)%particle%vx, &
                            p_species(idx)%particle%weight,1.0_f64/abs(p_species(idx)%qm))
        enddo

        call sll_collective_globalsum(sll_world_collective, energy, 0)
        energy=energy+kineticenergy_offset
    endfunction

    function sll_pic1d_calc_kineticenergy_offset(p_species ) &
            result(energy)
        type( sll_particle_1d_group), dimension(:), intent(in) :: p_species
        sll_int32 :: idx
        sll_int32 :: num_sp
        sll_real64 :: energy
        num_sp=size(p_species)
        energy=0
        do idx=1,num_sp
            energy=energy + &
            sll_pic1d_calc_kineticenergy_weighted(p_species(idx)%particle%vx, &
                            p_species(idx)%particle%weight_const,1.0_f64/abs(p_species(idx)%qm))
        enddo

        call sll_collective_globalsum(sll_world_collective, energy, 0)
    endfunction


    !<Calculate the kinetic energy for markers, in a consistent way
    function sll_pic1d_calc_kineticenergy_weighted( particle_v, particle_weight, particlemass ) &
            result(energy)
        sll_real64, dimension(:),intent(in) :: particle_v,particle_weight
        sll_real64, intent(in) ::particlemass
        sll_real64 :: energy
        SLL_ASSERT(size(particle_v)==size(particle_weight))

        energy= 0.5_f64*particlemass*&
                dot_product(particle_weight, particle_v**2)
    endfunction

    !<overall impulse 0.5*m*v
    function sll_pic1d_calc_impulse(p_species) &
            result(impulse)
        type( sll_particle_1d_group), dimension(:), intent(in) :: p_species
        sll_real64:: impulse
        sll_int32 :: idx, num_sp
        num_sp=size(p_species)

        impulse=0
        do idx=1,num_sp
            impulse=impulse + sll_pic1d_calc_impulse_weighted&
               (p_species(idx)%particle%weight, p_species(idx)%particle%vx, &
                            1.0_f64/abs(p_species(idx)%qm))
        enddo
        call sll_collective_globalsum(sll_world_collective, impulse, 0)

        impulse=impulse+impulse_offset
    endfunction

    function sll_pic1d_calc_impulse_offset(p_species) &
            result(impulse)
        type( sll_particle_1d_group), dimension(:), intent(in) :: p_species
        sll_real64:: impulse
        sll_int32 :: idx, num_sp
        num_sp=size(p_species)

        impulse=0
        do idx=1,num_sp
            impulse=impulse + sll_pic1d_calc_impulse_weighted&
               (p_species(idx)%particle%weight_const, p_species(idx)%particle%vx, &
                            1.0_f64/abs(p_species(idx)%qm))
        enddo

        call sll_collective_globalsum(sll_world_collective, impulse, 0)
    endfunction



    function sll_pic1d_calc_impulse_weighted(particle_v, particle_weight, particlemass) &
            result(impulse)
                sll_real64, dimension(:),intent(in) :: particle_v,particle_weight
        sll_real64, intent(in) ::particlemass
        sll_real64 :: impulse
        SLL_ASSERT(size(particle_v)==size(particle_weight))

        impulse= 0.5_f64*particlemass*&
                dot_product(particle_weight, particle_v)
    endfunction


    !<Thermal velocity, mostly defined as the variance of the particle velocity
    function sll_pic1d_calc_thermal_velocity(particlespeed, particleweight) &
            result(vth)
        sll_real64, DIMENSION(:), intent(in) :: particlespeed
        sll_real64, DIMENSION(:), intent(in):: particleweight
        sll_real64 :: vth
        sll_real64 :: mean
        mean=dot_product(particlespeed,particleweight)
        call sll_collective_globalsum(sll_world_collective, mean, 0)


        vth=dot_product((particlespeed-mean)**2,particleweight)
        call sll_collective_globalsum(sll_world_collective, vth, 0)
    endfunction


    !<Nonrelativistic field energy, namely $|\nabla \Phi$|^2$
    !<which is nothing else than the H1 seminorm of the self generated field
    !<additionally we add the potential energy by the external field E
    function sll_pic1d_calc_fieldenergy(p_species) &
            result(energy)
        type( sll_particle_1d_group), dimension(:), intent(in) :: p_species
        sll_real64 :: energy
        !Only for constant case
        sll_int32 :: idx
        sll_int32 :: num_sp
        num_sp=size(p_species)

        energy=0
        do idx=1,num_sp
            energy=energy+ &
                dot_product(Eex( p_species(idx)%particle%dx, 0.0_f64), &
                p_species(idx)%particle%dx)*(-1.0_f64)*p_species(idx)%qm/abs(p_species(idx)%qm)
        enddo

        call sll_collective_globalsum(sll_world_collective, energy, 0)

        !if (coll_rank==0)  energy=energy+0.5_f64*bspline_fem_solver_1d_H1seminorm_solution()/(interval_b-interval_a)
        if (coll_rank==0)  energy=energy+qnsolver%fieldenergy()


    endfunction




    !Write all results of the pic simulation to file
    subroutine sll_pic1d_write_result(filename, kinetic_energy, electrostatic_energy, impulse, &
            particleweight_mean, particleweight_var, perror_mean  , perror_var)
        character(len=*), intent(in) :: filename
        sll_real64, dimension(:), intent(in) :: kinetic_energy, electrostatic_energy, impulse,&
            particleweight_mean,particleweight_var, perror_mean  , perror_var
        integer :: k,file_id,file_id_err
        sll_real64,  dimension(size(kinetic_energy)) :: total_energy


        if (coll_rank==0) then

            SLL_ASSERT(size(kineticenergy)==timesteps+1)
            SLL_ASSERT(size(fieldenergy)==timesteps+1)
            SLL_ASSERT(size(impulse)==timesteps+1)

            total_energy=kinetic_energy +electrostatic_energy



            call sll_new_file_id(file_id, ierr)

            !Write Data File
            !            open(file_id, file = plot_name//"_"//fin//'.dat' )
            open(file_id, file = trim(root_path)//filename//'.dat')

            write (file_id, *)  "#Full 1d1v Electrostatic PIC"
            write (file_id,*)  "#Time steps:", timesteps
            write (file_id,*)  "#Time stepwidth:", timestepwidth
            write (file_id,*)  "#Marko particles:", coll_size*nparticles
            write (file_id,*)  "#Particle Pusher:", particle_pusher
            write (file_id,*)  "#Finite Elements: 2^(", log(real(mesh_cells,i64))/log(2.0_f64),")"
            write (file_id,*)  "#Size of MPI Collective: ", coll_size

            write (file_id,*)  "time  ", "kineticenergy  ", "electrostaticenergy  ", "impulse  ", &
                "vthermal  ", "weightmean   ", "weightvar  ", "perrormean  ", "perrorvar   "
            do k=1,timesteps+1
                write (file_id,*)  timestepwidth*(k-1), kineticenergy(k), fieldenergy(k), impulse(k), &
                    thermal_velocity_estimate(k), particleweight_mean(k),particleweight_var(k),&
                    perror_mean(k), perror_var(k)
            enddo
            close(file_id)

            call sll_new_file_id(file_id_err, ierr)

            open(file_id_err, file = trim(root_path)//filename//'-errors.dat')

            if (impulse(1)/=0 .AND. total_energy(1)/=0) then
                !Write the relative errors
                do k=1,timesteps+1
                    write (file_id_err,*)  timestepwidth*(k-1),&
                        abs((total_energy(k)-total_energy(1))/total_energy(1)), &
                        abs((  impulse(k)-impulse(1)))
                enddo
            endif
            close(file_id_err)


            !Write Gnuplot file
            open(file_id, file = trim(root_path)//filename//'.gnu')

            !First Plot
            write (file_id, *) "set term x11 1"

            write(file_id,"(A17,G10.3,A4,G10.3,A1)")"set title 'alpha=",landau_alpha,", k=",landau_mode,"'"
            write(file_id, "(A14,G10.5,A1)") "set xrange [0:",timestepwidth*timesteps,"]"
            write (file_id,*)  "set logscale y"
            write (file_id,*)  "unset logscale x"
            write (file_id,*)  "l_alpha=",landau_alpha
            write (file_id,"(A21,G4.3,A18,G10.3,A1)")  "set xlabel 'Time t in", timesteps, &
                " steps, stepwidth=", timestepwidth, "'"
            write (file_id,*)  "set ylabel 'Energy'"


            !Plot expected damping curve if known
            if (landau_alpha==0.001_f64 .AND. landau_mode==0.4_f64) then

                write (file_id,*) "E(x)=abs(0.002*0.424666*exp(-0.0661*x)*cos(1.2850*x-0.3357725))**2"
                write (file_id,*)  "plot '"//filename//".dat' using 1:3 with lines,\"
                write (file_id,*) "E(x) with lines linestyle 2;"
            elseif (landau_alpha/=0) then
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
            write (file_id,*)  "plot '"//trim(root_path)//"initial_phasespace.dat' using 2:3 smooth kdensity, vel_pdf(x)"
            write (file_id,*)  "set title 'Spacial Distribution'"
            write (file_id,*)  "plot '"//trim(root_path)//"initial_phasespace.dat' using 1:3 smooth kdensity"
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

        if (enable_deltaf .eqv. .TRUE.) then
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
    endsubroutine

    !    !This function should not have any side-effects
    !    function sll_pic1d_ensure_periodicity( particle_position, &
        !                    interval_a, interval_b) result( particle_position_out)
    !        sll_real, intent(in)::interval_a, interval_b
    !        sll_real, intent(in) :: particle_position
    !        sll_real :: interval_length
    !        sll_real :: particle_position_out
    !        SLL_ASSERT(interval_a < interval_b)
    !
    !        interval_length=interval_b - interval_a
    !        particle_position_out=modulo(particle_position- interval_a, interval_length)&
        !                            +interval_a
    !
    !
    !        SLL_ASSERT(particle_position>=interval_a)
    !        SLL_ASSERT(particle_position<=interval_b)
    !    endfunction

    !> \brief Check whether timestepwidth is accurate and scale down, when particles are to fast
    subroutine sll_pic1d_check_timescale( particle_position, &
            interval_a, interval_b, timestepwidth)
        sll_real64, intent(in)::interval_a, interval_b
        sll_real64, dimension(:) ,intent(in) :: particle_position
        sll_real64, intent(inout) :: timestepwidth
        sll_real64 :: max_offset

        !Calculate maximum relative offset
        max_offset = maxval(max( abs(particle_position- interval_a),  abs(particle_position- interval_b)))&
            /(interval_b-interval_a)

        if (1.0_f64/max_offset> 2_f64) then
            print *, "Maximal offset is:",max_offset
            timestepwidth=timestepwidth*0.5_f64
            print *, "Scaling down timestepwidth to:", timestepwidth
        endif
    endsubroutine




    !> \brief Adjust the weights after advection step
    subroutine sll_pic1d_adjustweights_advection_species(species_old, species_new)
        type(sll_particle_1d_group), dimension(:), intent(in) :: species_old
        type(sll_particle_1d_group), dimension(:), intent(inout) :: species_new
        sll_real64, dimension(:), allocatable:: ratio
        sll_int32 :: num_species, numpart, jdx

        num_species=size(species_old)
        SLL_ASSERT(num_species==size(species_new))

        !Adjust weights
        if (enable_deltaf .eqv. .TRUE.) then
            do jdx=1, num_species
                numpart=size(species_old(jdx)%particle)

                SLL_ALLOCATE(ratio(1:numpart),ierr)
                SLL_ASSERT(size(species_new(jdx)%particle)==numpart)
                ratio=1.0_f64
                ratio=control_variate_xv(species_new(jdx)%particle%dx, species_new(jdx)%particle%vx)&
                    /control_variate_xv(species_old(jdx)%particle%dx, species_old(jdx)%particle%vx)
                !                ratio=control_variate_x(species_new(jdx)%particle%dx)&
                    !                                    /control_variate_x(species_old(jdx)%particle%dx)

                species_new(jdx)%particle%weight=species_new(jdx)%particle%weight_const  - &
                    ( species_new(jdx)%particle%weight_const  - species_new(jdx)%particle%weight )*ratio

                !species_new(jdx)%particle%weight_const= species_old(jdx)%particle%weight_const

                ! species_old(jdx)%particle%weight=species_new(jdx)%particle%weight

                SLL_DEALLOCATE_ARRAY(ratio,ierr)
            enddo
        endif

    endsubroutine




    subroutine sll_pic1d_adjustweights(xold, xnew, vold, vnew)
        sll_real64, dimension(:),intent(in) ::vold, xold
        sll_real64, dimension(:),intent(out) ::vnew, xnew
        sll_real64, dimension(size(species(pushed_species)%particle)):: ratio
        sll_int32 :: N

        N=size(xold)
        SLL_ASSERT(N==size(xnew))
        SLL_ASSERT(N==size(vold))
        SLL_ASSERT(N==size(vnew))

        SLL_ASSERT(N==size( species(pushed_species)%particle))
        !Adjust weights
        if (enable_deltaf .eqv. .TRUE.) then
            ratio=control_variate_xv(xnew, vnew )/control_variate_xv(xold, vold )
            !ratio=ratio/(N*coll_size)
            !ratio=ratio
            species(pushed_species)%particle%weight=species(pushed_species)%particle%weight_const &
                -(species(pushed_species)%particle%weight_const-species(pushed_species)%particle%weight)*ratio
            !particleweight = sll_pic_1d_landaudamp_PDFxv(particleposition, vnew)/ sll_pic_1d_landaudamp_PDFxv(particleposition, vold)
            !particleweight = particleweight_constant - (particleweight_constant)
        endif

    endsubroutine


    !    !<Dummy function for Electric field
    !    function sll_pic_1d_electric_field(x) result(E)
    !        sll_real64, dimension(:), intent(in) :: x
    !        !sll_real64, intent(in) :: t
    !        sll_real64, dimension(size(x)) :: E
    !        E=0
    !        if (pic1d_testcase == SLL_PIC1D_TESTCASE_IONBEAM)         E=1.0_f64
    !    endfunction



    !<Dummy function for Electric field
    function pic_1d_electric_field_external(x,t) result(E)
        sll_real64, dimension(:), intent(in) :: x
        sll_real64, intent(in) :: t
        sll_real64, dimension(size(x)) :: E
        E=0.0_f64
        if (pic1d_testcase == SLL_PIC1D_TESTCASE_IONBEAM)   E=0

    endfunction



    function pic_1d_allparticlepos(p_species) result( ppos)
        type( sll_particle_1d_group), dimension(:), intent(in) :: p_species
        sll_int32 :: idx, num, npart
        sll_int32 :: off
        sll_real64, dimension(nparticles) :: ppos
        num=size(p_species)

        off=0

        do idx=1,num
            npart=size(p_species(idx)%particle)
            ppos(1+off:npart+off)=p_species(idx)%particle%dx
            off=npart
        enddo
    endfunction

    function pic_1d_allparticlev(p_species) result( pv)
        type( sll_particle_1d_group), dimension(:), intent(in) :: p_species
        sll_int32 :: idx, num, npart
        sll_int32 :: off
        sll_real64, dimension(nparticles) :: pv
        num=size(p_species)

        off=0

        do idx=1,num
            npart=size(p_species(idx)%particle)
            pv(1+off:npart+off)=p_species(idx)%particle%vx
            off=npart
        enddo
    endfunction

    function pic_1d_allparticleqm(p_species) result( qm)
        type( sll_particle_1d_group), dimension(:), intent(in) :: p_species
        sll_int32 :: idx, num, npart
        sll_int32 :: off
        sll_real64, dimension(nparticles) :: qm
        num=size(p_species)
        off=0
        do idx=1,num
            npart=size(p_species(idx)%particle)
            qm(1+off:npart+off)=p_species(idx)%qm
            off=npart
        enddo
    endfunction

endmodule  sll_pic_1d_Class


